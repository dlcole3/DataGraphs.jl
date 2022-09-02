function get_EC(g::DataGraph)
    nodes = g.nodes
    edges = g.edges

    EC = length(nodes) - length(edges)

    return EC
end

function matrix_to_graph(matrix, weight_name::String="weight")

    dim1, dim2 = size(matrix)
    g = DataGraph()

    fadjlist     = [Vector{Int}() for i in 1:(dim1 * dim2)]
    edges        = Vector{Tuple{Int, Int}}()
    nodes        = Vector{Any}()
    nodes_index  = Dict{Any, Int}()
    edges_index  = Dict{Tuple{Int, Int}, Int}()
    node_data    = NamedArray(fill(NaN, (dim1*dim2,1)))
    setnames!(node_data, [weight_name], 2)

    for j in 1:dim2
        for i in 1:dim1
            push!(nodes, (i, j))
            nodes_index[(i, j)] = length(nodes)
            node_data[length(nodes), weight_name] = matrix[i, j]
        end
    end

    if dim1 > 1 && dim2 > 1
        for j in 1:dim2
            column_offset = dim2 * (j - 1)
            for i in 1:dim1
                if j != dim2
                    edge = (i + column_offset, i + column_offset + dim1)
                    push!(edges, edge)
                    push!(fadjlist[i + column_offset], i + column_offset + dim1)
                    push!(fadjlist[i + column_offset + dim1], i + column_offset)
                    edges_index[edge] = length(edges)
                end

                if i != dim1
                    edge = (i + column_offset, i + column_offset + 1)
                    push!(fadjlist[i + column_offset], i + column_offset + 1)
                    push!(fadjlist[i + column_offset + 1], i + column_offset)
                    push!(edges, edge)
                    edges_index[edge] = length(edges)
                end

            end
        end
    end

    g.ne              = length(edges)
    g.fadjlist        = fadjlist
    g.nodes           = nodes
    g.nodes_index     = nodes_index
    g.edges           = edges
    g.edges_index     = edges_index
    g.node_data       = node_data
    g.node_attributes = [weight_name]

    return g
end

function sym_matrix_to_graph(matrix, weight_name::String="weight", tol = 1e-9)

    dim1, dim2 = size(matrix)
    if dim1 != dim2
        error("Matrix is not square")
    end

    if sum(abs.(matrix - transpose(matrix))) > tol
        error("Matrix is not symmetric")
    end

    g = DataGraph()

    for j in 1:dim2
        for i in j:dim1
            add_edge!(g, i, j)
            add_edge_data!(g, i, j, matrix[i,j], weight_name)
        end
    end

    return g
end

function mvts_to_graph(mvts, weight_name::String="weight", tol=1e-9)

    mvts_cov  = cov(mvts)
    mvts_prec = inv(mvts_cov)

    g = sym_matrix_to_graph(mvts_prec, weight_name, tol)

    return g
end

function filter_nodes(g::DataGraph, filter_val::Real; attribute::String=g.node_attributes[1])
    node_attributes = g.node_attributes
    edge_attributes = g.edge_attributes
    node_data       = g.node_data
    node_index      = g.node_index
    nodes           = g.nodes
    edge_data    = g.edge_data
    edges_index     = g.edges_index
    edges           = g.edges
    node_positions  = g.node_positions

    if length(node_attributes) == 0
        error("No node weights are defined")
    end

    T = eltype(g)

    new_g = DataGraph()

    am = create_adj_mat(g)

    bool_vec = g.node_data[:, attribute] .< filter_val

    new_am = am[bool_vec, bool_vec]

    new_nodes     = nodes[bool_vec]
    new_node_data = node_data[bool_vec, :]

    if length(node_positions) > 0 && length(node_positions) == length(nodes)
        new_node_pos = node_positions[bool_vec]
        new_g.node_positions  = new_node_pos
    end

    new_node_index = Dict()

    new_edges      = []
    new_edge_index = Dict()
    old_edge_map   = []
    fadjlist       = [Vector{T}() for i in 1:length(new_nodes)]

    for i in 1:length(new_nodes)
        new_node_index[new_nodes[i]] = i
    end

    for j in 1:length(new_nodes)
        for i in (j + 1):length(new_nodes)
            if new_am[i,j]
                new_edge = (i, j)
                push!(new_edges, new_edge)
                new_edge_index[(new_edge)] = length(new_edges)

                @inbounds node_neighbors = fadjlist[i]
                index = searchsortedfirst(node_neighbors, j)
                insert!(node_neighbors, index, j)

                @inbounds node_neighbors = fadjlist[j]
                index = searchsortedfirst(node_neighbors, i)
                insert!(node_neighbors, index, i)

                old_edge = _get_node(nodes_index(new_nodes[j]), nodes_index(new_nodes[i]))
                if old_edge in edges
                    push!(old_edge_map, edges_index[old_edge])
                end
            end
        end
    end

    if length(edge_attributes) > 0
        new_edge_data = edge_data[old_edge_index, :]
        new_g.edge_data    = new_edge_data
    end

    new_g.ne              = length(new_edges)
    new_g.fadjlist        = fadjlist
    new_g.nodes           = new_nodes
    new_g.edges           = new_edges
    new_g.edges_index     = new_edge_index
    new_g.nodes_index     = new_node_index
    new_g.node_attributes = attributes
    new_g.edge_attributes = edge_attributes
    new_g.node_data       = new_node_data
    new_g.node_positions  = new_node_pos
    new_g.adj_mat         = new_am

    return new_g
end

function filter_edges(g::DataGraph, filter_val::Real; attribute::String = g.edge_attributes[1])
    nodes           = g.nodes
    edges           = g.edges
    node_attributes = g.node_attributes
    edge_attributes = g.edge_attributes
    edge_data       = g.edge_data
    node_index      = g.nodes_index

    if length(edge_attributes) == 0
        error("No node weights are defined")
    end

    T = eltype(g)

    bool_vec = g.edge_data[:, attribute] .< filter_val

    new_edges = edges[bool_vec]
    new_edge_data = edge_data[bool_vec, :]

    new_edge_index = Dict()

    fadjlist = [Vector{T}() for i in 1:length(nodes)]
    for i in 1:length(new_edges)
        new_edge_index[new_edges[i]] = i

        node1, node2 = new_edges[i]

        @inbounds node_neighbors = fadjlist[node1]
        index = searchsortedfirst(node_neighbors, node2)
        insert!(node_neighbors, index, node2)

        @inbounds node_neighbors = fadjlist[node2]
        index = searchsortedfirst(node_neighbors, node1)
        insert!(node_neighbors, index, node1)
    end

    new_g = DataGraph()

    new_g.ne              = T(length(new_edges))
    new_g.fadjlist        = fadjlist
    new_g.nodes           = nodes
    new_g.edges           = new_edges
    new_g.edge_data       = new_edge_data
    new_g.nodes_index     = node_index
    new_g.edges_index     = new_edge_index
    new_g.node_attributes = node_attributes
    new_g.edge_attributes = edge_attributes
    new_g.node_positions  = g.node_positions
    new_g.node_data       = g.node_data

    return new_g
end

function run_EC_on_nodes(g::DataGraph, thresh; attribute::String = g.node_attributes[1])
    nodes        = g.nodes
    node_data    = g.node_data

    am = create_adj_mat(g)

    for i in 1:length(nodes)
        if am[i, i] == 1
            am[i, i] = 2
        end
    end

    ECs = zeros(length(thresh))

    for (j,i) in enumerate(thresh)
        bool_vec  = node_data[:, attribute] .< i
        new_am    = am[bool_vec, bool_vec]
        num_nodes = sum(bool_vec)
        num_edges = sum(new_am)/2
        ECs[j]    = num_nodes-num_edges
    end

    return ECs
end

function run_EC_on_edges(g::DataGraph, thresh; attribute::String = g.edge_attributes[1])
    edge_data = g.edge_data
    nodes        = g.nodes

    ECs = zeros(length(thresh))

    num_nodes = length(nodes)

    for (j,i) in enumerate(thresh)
        bool_vec  = edge_data[:, attribute] .< i
        num_edges = sum(bool_vec)
        ECs[j]    = num_nodes - num_edges
    end

    return ECs
end

function aggregate(g::DataGraph, node_set, new_name)
    nodes           = g.nodes
    nodes_index     = g.nodes_index
    node_data       = g.node_data
    node_attributes = g.node_attributes
    node_positions  = g.node_positions

    if !(issubset(node_set, nodes))
        undef_nodes = setdiff(node_set, nodes)
        println()
        for i in undef_nodes
            println("Node $i is not defined in graph")
        end
        error("Node set includes nodes that are not defined")
    end

    if new_name in setdiff(nodes, node_set)
        error("New node name already exists in set of non-aggregated nodes")
    end

    T = eltype(g)

    new_g = DataGraph()

    new_nodes = setdiff(nodes, node_set)
    push!(new_nodes, new_name)

    new_node_dict = Dict()

    for (i,j) in enumerate(new_nodes)
        new_node_dict[j] = i
    end

    # Get indices of old nodes
    old_indices = []
    for i in node_set
        old_index = nodes_index[i]
        push!(old_indices, old_index)
    end

    indices_to_keep = setdiff(1:length(nodes), old_indices)

    if length(node_attributes) > 0
        node_data_to_avg   = node_data[old_indices,:]
        node_weight_avg    = Statistics.mean(node_data_to_avg; dims=1)

        node_data_to_keep = node_data[indices_to_keep, :]
        new_node_data     = vcat(node_data_to_keep, node_weight_avg)
        setnames!(new_node_data, node_attributes, 2)

        new_g.node_attributes = node_attributes
        new_g.node_data       = new_node_data
    end

    if length(node_positions) > 1
        new_node_positions = node_positions[indices_to_keep]
        old_pos            = node_positions[old_indices]

        xvals = 0
        yvals = 0

        for j in 1:length(node_set)
            xvals += old_pos[j][1]
            yvals += old_pos[j][2]
        end

        push!(new_node_positions, Point(xvals/length(node_set), yvals/length(node_set)))
        new_g.node_positions = new_node_positions
    end

    edges           = g.edges
    edge_data       = g.edge_data
    edge_attributes = g.edge_attributes
    edges_index     = g.edges_index

    fadjlist = [Vector{T}() for i in 1:length(new_nodes)]

    node_name_mapping = Dict{T, Any}()
    new_edges    = Vector{Tuple{T, T}}()
    new_edges_index  = Dict{Tuple{T, T}, T}()
    edge_bool_vec = [false for i in 1:length(edges)]
    #edge_bool_vec_avg = [false for i in 1:length(edges)]
    edge_bool_avg_index = Dict{Tuple{T, T}, Vector{T}}()
    new_edge_data = fill(NaN, (0, length(edge_attributes)))

    for i in 1:length(nodes)
        node_name_mapping[nodes_index[nodes[i]]] = nodes[i]
    end

    for (i, edge) in enumerate(edges)
        node1_bool = edge[1] in old_indices
        node2_bool = edge[2] in old_indices

        if !node1_bool && !node2_bool
            new_node1 = new_node_dict[node_name_mapping[edge[1]]]
            new_node2 = new_node_dict[node_name_mapping[edge[2]]]

            push!(new_edges, (new_node1, new_node2))
            new_edges_index[(new_node1, new_node2)] = length(edges)

            @inbounds node_neighbors = fadjlist[new_node1]
            index = searchsortedfirst(node_neighbors, new_node2)
            insert!(node_neighbors, index, new_node2)

            @inbounds node_neighbors = fadjlist[new_node2]
            index = searchsortedfirst(node_neighbors, new_node1)
            insert!(node_neighbors, index, new_node1)

            if length(edge_attributes) > 0
                new_edge_data = vcat(new_edge_data, edge_data[edges_index[edge], :]')
            end

            edge_bool_vec[i] = true
        elseif !node1_bool && node2_bool
            new_node1 = new_node_dict[node_name_mapping[edge[1]]]
            new_node2 = length(new_nodes)

            if !((new_node1, new_node2) in edges)
                push!(new_edges, (new_node1, new_node2))
                new_edges_index[(new_node1, new_node2)] = length(edges)

                @inbounds node_neighbors = fadjlist[new_node1]
                index = searchsortedfirst(node_neighbors, new_node2)
                insert!(node_neighbors, index, new_node2)

                @inbounds node_neighbors = fadjlist[new_node2]
                index = searchsortedfirst(node_neighbors, new_node1)
                insert!(node_neighbors, index, new_node1)

                if length(edge_attributes) > 0
                    new_edge_data = vcat(new_edge_data, fill(NaN, (1, length(edge_attributes))))
                    edge_bool_avg_index[(new_node1, new_node2)] = Vector{T}([edges_index[edge]])
                end
            else
                push!(edge_bool_avg_index[(new_node1, new_node2)], edges_index[edge])
            end
        elseif node1_bool && !node2_bool
            new_node1 = new_node_dict[node_name_mapping[edge[2]]]
            new_node2 = length(new_nodes)

            if !((new_node1, new_node2) in edges)
                push!(new_edges, (new_node1, new_node2))
                new_edges_index[(new_node1, new_node2)] = length(edges)

                @inbounds node_neighbors = fadjlist[new_node1]
                index = searchsortedfirst(node_neighbors, new_node2)
                insert!(node_neighbors, index, new_node2)

                @inbounds node_neighbors = fadjlist[new_node2]
                index = searchsortedfirst(node_neighbors, new_node1)
                insert!(node_neighbors, index, new_node1)

                if length(edge_attributes) > 0
                    new_edge_data = vcat(new_edge_data, fill(NaN, (1, length(edge_attributes))))
                    edge_bool_avg_index[(new_node1, new_node2)] = Vector{T}([edges_index[edge]])
                end
            else
                if length(edge_attributes) > 0
                    push!(edge_bool_avg_index[(new_node1, new_node2)], edges_index[edge])
                end
            end
        end
    end

    if length(edge_attributes) > 0

        for edge in keys(edge_bool_avg_index)
            new_index = new_edges_index[edge]
            edge_data_to_avg = edge_data[edge_bool_avg_index[edge], :]

            edge_data_avg = Statistics.mean(edge_data_to_avg; dims = 1)
            new_edge_data[new_index, :] = edge_data_avg[:]
        end

        new_edge_data = NamedArray(new_edge_data)
        setnames!(new_edge_data, edge_attributes, 2)

        new_g.edge_attributes = edge_attributes
        new_g.edge_data       = new_edge_data
    end

    new_g.ne          = length(new_edges)
    new_g.fadjlist    = fadjlist
    new_g.nodes       = new_nodes
    new_g.nodes_index = new_node_dict
    new_g.edges       = new_edges
    new_g.edges_index = new_edges_index

    return new_g
end
