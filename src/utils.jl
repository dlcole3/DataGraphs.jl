function get_EC(g::Graph)
    nodes = g.nodes
    edges = g.edges

    EC = length(nodes) - length(edges)

    return EC
end

function matrix_to_graph(matrix, weight_name::String="weight")

    dim1, dim2 = size(matrix)
    g = Graph()

    edges        = Vector{Tuple{Int}}()
    nodes        = Vector{Any}()
    nodes_index  = Dict{Any, Int}()
    edges_index  = Dict{Tuple{Int}, Int}()
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
                if i != dim1
                    edge = (i + column_offset, i + column_offset + 1)
                    push!(edges, edge)
                    edges_index[edge] = length(edges)
                end
                if j != dim2
                    edge = (i + column_offset, i + column_offset + dim1)
                    push!(edges, edge)
                    edges_index[edge] = length(edges)
                end
            end
        end
    end

    g.nodes           = nodes
    g.nodes_index     = nodes_index
    g.edges           = unique_edges
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

    g = Graph()

    for j in 1:dim2
        for i in 1:dim1
            add_edge!(g, i, j)
            add_edge_weight!(g, i, j, matrix[i,j], weight_name)
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


function filter_nodes(g::Graph, filter_val::Real; attribute::String=g.node_attributes[1])
    node_attributes = g.node_attributes
    edge_attributes = g.edge_attributes
    node_data       = g.node_data
    node_index      = g.node_index
    nodes           = g.nodes
    edge_weights    = g.edge_data
    edges_index     = g.edges_index
    edges           = g.edges
    node_positions  = g.node_positions

    if length(node_attributes) == 0
        error("No node weights are defined")
    end

    new_g = Graph()

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
    old_edge_map = []

    for i in 1:length(new_nodes)
        new_node_index[new_nodes[i]] = i
    end

    for j in 1:length(new_nodes)
        for i in (j + 1):length(new_nodes)
            if new_am[i,j]
                new_edge = (i, j)
                push!(new_edges, new_edge)
                new_edge_index[(new_edge)] = length(new_edges)

                old_edge = _get_node(nodes_index(new_nodes[j]), nodes_index(new_nodes[i]))
                if old_edge in edges
                    push!(old_edge_map, edges_index[old_edge])
            end
        end
    end

    if length(edge_attributes) > 0
        new_edge_weights = edge_weights[old_edge_index, :]
        new_g.edge_weights    = new_edge_weights
    end

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

function filter_edges(g::Graph, filter_val::Real; attribute::String = g.edge_attributes[1])
    nodes           = g.nodes
    edges           = g.edges
    node_attributes = g.node_attributes
    edge_attributes = g.edge_attributes
    edge_data       = g.edge_data
    node_index      = g.nodes_index

    if length(edge_attributes) == 0
        error("No node weights are defined")
    end

    bool_vec = g.edge_weights[:, attribute] .< filter_val

    new_edges = edges[bool_vec]
    new_edge_weights = edge_weights[bool_vec, :]

    new_edge_index = Dict()

    for i in 1:length(new_edges)
        new_edge_index[new_edges[i]] = i
    end


    new_g = Graph()

    new_g.nodes           = nodes
    new_g.edges           = new_edges
    new_g.edge_weights    = new_edge_weights
    new_g.nodes_index     = node_index
    new_g.edges_index     = new_edge_index
    new_g.node_attributes = node_attributes
    new_g.edge_attributes = edge_attributes
    new_g.node_positions  = g.node_positions
    new_g.node_data       = g.node_data

    return new_g

end

function run_EC_on_nodes(g::Graph, thresh; attribute::String = g.node_attributes[1])
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

function run_EC_on_edges(g::Graph, thresh; attribute::String = g.edge_attributes[1])
    edge_weights = g.edge_weights
    nodes        = g.nodes

    ECs = zeros(length(thresh))

    num_nodes = length(nodes)

    for (j,i) in enumerate(thresh)
        bool_vec  = edge_weights[:, attribute] .< i
        num_edges = sum(bool_vec)
        ECs[j]    = num_nodes - num_edges
    end

    return ECs
end

function aggregate(g::Graph, node_set, new_name)
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

    new_g = Graph()

    new_nodes = setdiff(nodes, node_set)
    push!(new_nodes, new_name)

    new_node_dict = Dict()

    for (i,j) in enumerate(new_nodes)
        new_node_dict[j] = i
    end

    new_g.nodes       = new_nodes
    new_g.nodes_index = new_node_dict

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
    edge_weights    = g.edge_weights
    edge_attributes = g.edge_attributes
    edges_index     = g.edges_index

    if length(edge_attributes) > 0

    end

    am = create_adj_mat(g)

    new_am = am[indices_to_keep, indices_to_keep]

    new_edges      = []
    new_edge_index = Dict()


    if length(node_att) > 0
        node_data_to_avg = node_data[old_indices,:]
        node_weight_avg     = Statistics.mean(node_data_to_avg; dims=1)

        node_data_to_keep = node_data[indices_to_keep, :]
        new_node_data     = vcat(node_data_to_keep, node_weight_avg)
        setnames!(new_node_data, node_att, 2)
    end

    if length(edge_attributes) == 0

        for i in 1:length(edges)
            if setdiff(node_set, edges[i]) == node_set
                push!(new_edges, edges[i])
                new_edge_index[edges[i]] = length(new_edges)
            end
        end


        old_am = am[old_indices, indices_to_keep]
        sum_old_am = sum(old_am, dims=1)

        for i in 1:length(indices_to_keep)
            if sum_old_am[i] >= 1
                new_edge = (nodes[indices_to_keep[i]], new_name)
                push!(new_edges, new_edge)
                new_edge_index[new_edge] = length(new_edges)
            end
        end

        new_g.edges       = new_edges
        new_g.edges_index = new_edge_index
    else
        old_edges = []

        for i in 1:length(edges)
            if setdiff(node_set, edges[i]) == node_set
                push!(new_edges, edges[i])
                new_edge_index[edges[i]] = length(new_edges)
                push!(old_edges, edges_index[edges[i]])
            end
        end

        new_edge_weights = edge_weights[old_edges, :]

        old_am = am[old_indices, indices_to_keep]


        sum_old_am = sum(old_am, dims=1)

        multi_edges = sum_old_am .> 1
        new_edge_weights


        for i in 1:length(indices_to_keep)
            if sum_old_am[i] == 1
                new_edge = (nodes[indices_to_keep[i]], new_name)
                push!(new_edges, new_edge)
                new_edge_index[new_edge] = length(new_edges)

                for j in 1:length(node_set)
                    if (nodes[indices_to_keep[i]], node_set[j]) in edges
                        row_to_add = edge_weights[edges_index[(nodes[indices_to_keep[i]], node_set[j])], :]
                    else
                        row_to_add = edge_weights[edges_index[(node_set[j], nodes[indices_to_keep[i]])], :]
                    end
                end

                new_edge_weights = vcat(new_edge_weights, row_to_add)

            elseif sum_old_am[i] > 1
                new_edge = (nodes[indices_to_keep[i]], new_name)
                push!(new_edges, new_edge)
                new_edge_index[new_edge] = length(new_edges)

                edges_to_replace = []

                for j in 1:length(node_set)
                    if (nodes[indices_to_keep[i]], node_set[j]) in edges
                        push!(edges_to_replace, edges_index[(nodes[indices_to_keep[i]], node_set[j])])
                    end
                    if (node_set[j], nodes[indices_to_keep[i]]) in edges
                        push!(edges_to_replace, edges_index[(node_set[j], nodes[indices_to_keep[i]])])
                    end
                end
                println(edges_to_replace)

                rows_to_avg = edge_weights[edges_to_replace, :]
                avgs        = Statistics.mean(rows_to_avg, dims=1)
                new_edge_weights = vcat(new_edge_weights, avgs)

            end
        end
        setnames!(new_edge_weights, edge_attributes, 2)
        new_g.edges        = new_edges
        new_g.edges_index  = new_edge_index
        new_g.edge_weights = new_edge_weights
        new_g.edge_attributes = edge_attributes

    end

    return new_g
end
