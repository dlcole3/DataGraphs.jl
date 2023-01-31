"""
    filter_nodes(datadigraph, filter_value; attribute_name)

Removes the nodes of the graph whose weight value of `attribute_name` is greater than the given
`filter_value`. If `attribute_name` is not specified, this defaults to the first attribute within
the DataGraph's `NodeData`.
"""
function filter_nodes(dg::DataDiGraph, filter_val::R; attribute::String=dg.node_data.attributes[1]) where {R <: Real}
    node_attributes    = dg.node_data.attributes
    edge_attributes    = dg.edge_data.attributes
    node_attribute_map = dg.node_data.attribute_map
    edge_attribute_map = dg.edge_data.attribute_map
    node_data          = dg.node_data.data
    node_map           = dg.node_map
    nodes              = dg.nodes
    edge_data          = dg.edge_data.data
    edge_map           = dg.edge_map
    edges              = dg.edges
    node_positions     = dg.node_positions

    if length(node_attributes) == 0
        error("No node weights are defined")
    end

    T = eltype(dg)
    T1 = eltype(get_node_data(dg))
    M1 = typeof(get_node_data(dg))
    T2 = eltype(get_edge_data(dg))
    M2 = eltype(get_edge_data(dg))

    new_dg = DataDiGraph{T, T1, T2, M1, M2}()

    am = Graphs.LinAlg.adjacency_matrix(dg.g)

    bool_vec = node_data[:, node_attribute_map[attribute]] .< filter_val

    new_am = am[bool_vec, bool_vec]

    new_nodes     = nodes[bool_vec]
    new_node_data = node_data[bool_vec, :]

    if length(node_positions) > 0 && length(node_positions) == length(nodes)
        new_node_pos = node_positions[bool_vec]
        new_dg.node_positions  = new_node_pos
    else
        new_node_pos = [[0.0 0.0]]
    end

    new_node_map = Dict()

    new_edges      = Vector{Tuple{T, T}}()
    new_edge_map   = Dict{Tuple{T, T}, T}()
    old_edge_index = Vector{Int}()
    fadjlist       = Vector{Vector{T}}([Vector{T}() for i in 1:length(new_nodes)])  ### TODO: if new_nodes is length 0, this is a vector of type any
    badjlist       = Vector{Vector{T}}([Vector{T}() for i in 1:length(new_nodes)])  ### TODO: if new_nodes is length 0, this is a vector of type any

    for i in 1:length(new_nodes)
        new_node_map[new_nodes[i]] = i
    end

    for j in 1:length(new_nodes)
        for i in 1:length(new_nodes)
            if new_am[i,j] == 1
                new_edge = (i, j)
                push!(new_edges, new_edge)
                new_edge_map[new_edge] = length(new_edges)

                @inbounds node_neighbors = fadjlist[i]
                index = searchsortedfirst(node_neighbors, j)
                insert!(node_neighbors, index, j)

                @inbounds node_neighbors = badjlist[j]
                index = searchsortedfirst(node_neighbors, i)
                insert!(node_neighbors, index, i)

                old_edge = _get_edge(node_map[new_nodes[i]], node_map[new_nodes[j]])
                if old_edge in edges
                    push!(old_edge_index, edge_map[old_edge])
                end
            end
        end
    end

    if length(edge_attributes) > 0
        new_edge_data         = edge_data[old_edge_index, :]
        new_dg.edge_data.data = new_edge_data
        new_dg.edge_data.attributes_map = dg.edge_data.attributes_map
    end

    simple_digraph = Graphs.SimpleDiGraph(T(length(new_edges)), fadjlist, badjlist)

    new_dg.g                    = simple_digraph
    new_dg.nodes                = new_nodes
    new_dg.edges                = new_edges
    new_dg.edge_map             = new_edge_map
    new_dg.node_map             = new_node_map
    new_dg.node_data.attributes = node_attributes
    new_dg.edge_data.attributes = edge_attributes
    new_dg.node_data.data       = new_node_data
    new_dg.node_positions       = new_node_pos
    new_dg.node_data.attribute_map = dg.node_data.attribute_map

    return new_dg
end

"""
    filter_edges(datadigraph, filter_value; attribute_name)

Removes the edges of the graph whose weight value of `attribute_name` is greater than the given
`filter_value`. If `attribute_name` is not specified, this defaults to the first attribute within
the DataGraph's `EdgeData`.
"""
function filter_edges(dg::DataDiGraph, filter_val::R; attribute::String = dg.edge_data.attributes[1]) where {R <: Real}
    nodes           = dg.nodes
    edges           = dg.edges
    node_attributes = dg.node_data.attributes
    edge_attributes = dg.edge_data.attributes
    edge_data       = dg.edge_data.data
    node_map        = dg.node_map

    node_attribute_map = dg.node_data.attribute_map
    edge_attribute_map = dg.edge_data.attribute_map

    if length(edge_attributes) == 0
        error("No node weights are defined")
    end

    T = eltype(dg)
    T1 = eltype(get_node_data(dg))
    M1 = typeof(get_node_data(dg))
    T2 = eltype(get_edge_data(dg))
    M2 = eltype(get_edge_data(dg))

    bool_vec = dg.edge_data.data[:, edge_attribute_map[attribute]] .< filter_val

    new_edges = edges[bool_vec]
    new_edge_data = edge_data[bool_vec, :]

    new_edge_map = Dict{Tuple{T, T}, T}()

    fadjlist = Vector{Vector{T}}([Vector{T}() for i in 1:length(nodes)])   ### TODO: if new_nodes is length 0, this is a vector of type any
    badjlist = Vector{Vector{T}}([Vector{T}() for i in 1:length(nodes)])   ### TODO: if new_nodes is length 0, this is a vector of type any

    for i in 1:length(new_edges)
        new_edge_map[new_edges[i]] = i

        node1, node2 = new_edges[i]

        @inbounds node_neighbors = fadjlist[node1]
        index = searchsortedfirst(node_neighbors, node2)
        insert!(node_neighbors, index, node2)

        @inbounds node_neighbors = fadjlist[node2]
        index = searchsortedfirst(node_neighbors, node1)
        insert!(node_neighbors, index, node1)
    end

    new_dg = DataDiGraph{T, T1, T2, M1, M2}()

    simple_digraph = Graphs.SimpleDiGraph(T(length(new_edges)), fadjlist, badjlist)

    new_dg.g                    = simple_digraph
    new_dg.nodes                = nodes
    new_dg.edges                = new_edges
    new_dg.edge_data.data       = new_edge_data
    new_dg.node_map             = node_map
    new_dg.edge_map             = new_edge_map
    new_dg.node_data.attributes = node_attributes
    new_dg.edge_data.attributes = edge_attributes
    new_dg.node_positions       = dg.node_positions
    new_dg.node_data.data       = dg.node_data.data

    new_dg.node_data.attribute_map = dg.node_data.attribute_map
    new_dg.edge_data.attribute_map = dg.edge_data.attribute_map

    return new_dg
end

"""
    aggregate(datadigraph, node_list, aggregated_node_name)

Aggregates all the nodes in `node_list` into a single node which is called `aggregated_node_name`.
If nodes have any weight/attribute values defined, These are averaged across all values in the
`node_list`. Edge weights are also averaged when two or more nodes in the `node_list` are connected
to the same node and these edges have weights defined on them.
"""
function aggregate(dg::DataDiGraph, node_set, new_name)
    nodes              = dg.nodes
    node_map           = dg.node_map
    node_data          = dg.node_data.data
    node_attributes    = dg.node_data.attributes
    node_attribute_map = dg.node_data.attribute_map
    node_positions     = dg.node_positions

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

    T = eltype(dg)
    T1 = eltype(get_node_data(dg))
    M1 = typeof(get_node_data(dg))
    T2 = eltype(get_edge_data(dg))
    M2 = eltype(get_edge_data(dg))

    new_dg = DataDiGraph{T, T1, T2, M1, M2}()

    new_nodes = setdiff(nodes, node_set)
    push!(new_nodes, new_name)

    new_node_dict = Dict()

    for (i,j) in enumerate(new_nodes)
        new_node_dict[j] = i
    end

    # Get indices of old nodes
    agg_node_indices = []
    for i in node_set
        old_index = node_map[i]
        push!(agg_node_indices, old_index)
    end

    indices_to_keep = setdiff(1:length(nodes), agg_node_indices)

    if length(node_attributes) > 0
        node_data_to_avg   = node_data[agg_node_indices, :]
        node_weight_avg    = Statistics.mean(node_data_to_avg; dims=1)

        node_data_to_keep = node_data[indices_to_keep, :]
        new_node_data     = vcat(node_data_to_keep, node_weight_avg)

        new_dg.node_data.attributes    = node_attributes
        new_dg.node_data.attribute_map = node_attribute_map
        new_dg.node_data.data          = new_node_data
    end

    if length(node_positions) > 1
        new_node_positions = node_positions[indices_to_keep]
        old_pos            = node_positions[agg_node_indices]

        xvals = 0
        yvals = 0

        for j in 1:length(node_set)
            xvals += old_pos[j][1]
            yvals += old_pos[j][2]
        end

        push!(new_node_positions, Point(xvals/length(node_set), yvals/length(node_set)))
        new_dg.node_positions = new_node_positions
    end

    edges              = dg.edges
    edge_data          = dg.edge_data.data
    edge_attributes    = dg.edge_data.attributes
    edge_attribute_map = dg.edge_data.attribute_map
    edge_map           = dg.edge_map

    fadjlist = Vector{Vector{T}}([Vector{T}() for i in 1:length(new_nodes)])   ### TODO: if new_nodes is length 0, this is a vector of type any
    badjlist = Vector{Vector{T}}([Vector{T}() for i in 1:length(new_nodes)])   ### TODO: if new_nodes is length 0, this is a vector of type any


    node_name_mapping   = Dict{T, Any}()
    new_edges           = Vector{Tuple{T, T}}()
    new_edge_map        = Dict{Tuple{T, T}, T}()
    edge_bool_vec       = [false for i in 1:length(edges)]
    edge_bool_avg_index = Dict{Tuple{T, T}, Vector{T}}()
    new_edge_data       = fill(0, (0, length(edge_attributes)))

    for i in 1:length(nodes)
        node_name_mapping[node_map[nodes[i]]] = nodes[i]
    end

    for (i, edge) in enumerate(edges)
        node1_bool = edge[1] in agg_node_indices
        node2_bool = edge[2] in agg_node_indices

        if !node1_bool && !node2_bool
            new_node1 = new_node_dict[node_name_mapping[edge[1]]]
            new_node2 = new_node_dict[node_name_mapping[edge[2]]]

            push!(new_edges, (new_node1, new_node2))
            new_edge_map[(new_node1, new_node2)] = length(new_edges)

            @inbounds node_neighbors = fadjlist[new_node1]
            index = searchsortedfirst(node_neighbors, new_node2)
            insert!(node_neighbors, index, new_node2)

            @inbounds node_neighbors = badjlist[new_node2]
            index = searchsortedfirst(node_neighbors, new_node1)
            insert!(node_neighbors, index, new_node1)

            if length(edge_attributes) > 0
                new_edge_data = vcat(new_edge_data, edge_data[edge_map[edge], :]')
            end

            edge_bool_vec[i] = true

        elseif !node1_bool && node2_bool
            new_node1 = new_node_dict[node_name_mapping[edge[1]]]
            new_node2 = length(new_nodes)

            if !((new_node1, new_node2) in new_edges)
                push!(new_edges, (new_node1, new_node2))
                new_edge_map[(new_node1, new_node2)] = length(new_edges)

                @inbounds node_neighbors = fadjlist[new_node1]
                index = searchsortedfirst(node_neighbors, new_node2)
                insert!(node_neighbors, index, new_node2)

                @inbounds node_neighbors = badjlist[new_node2]
                index = searchsortedfirst(node_neighbors, new_node1)
                insert!(node_neighbors, index, new_node1)

                if length(edge_attributes) > 0
                    new_edge_data = vcat(new_edge_data, fill(0, (1, length(edge_attributes))))
                    edge_bool_avg_index[(new_node1, new_node2)] = Vector{T}([edge_map[edge]])
                end
            else
                if length(edge_attributes) > 0
                    if (new_node1, new_node2) in keys(edge_bool_avg_index)
                        push!(edge_bool_avg_index[(new_node1, new_node2)], edge_map[edge])
                    else
                        edge_bool_avg_index[(new_node1, new_node2)] = Vector{T}([edge_map[edge]])
                    end
                end
            end
        elseif node1_bool && !node2_bool
            new_node1 = new_node_dict[node_name_mapping[edge[2]]]
            new_node2 = length(new_nodes)

            if !((new_node1, new_node2) in new_edges)
                push!(new_edges, (new_node1, new_node2))
                new_edge_map[(new_node1, new_node2)] = length(new_edges)

                @inbounds node_neighbors = fadjlist[new_node1]
                index = searchsortedfirst(node_neighbors, new_node2)
                insert!(node_neighbors, index, new_node2)

                @inbounds node_neighbors = badjlist[new_node2]
                index = searchsortedfirst(node_neighbors, new_node1)
                insert!(node_neighbors, index, new_node1)

                if length(edge_attributes) > 0
                    new_edge_data = vcat(new_edge_data, fill(0, (1, length(edge_attributes))))
                    edge_bool_avg_index[(new_node1, new_node2)] = Vector{T}([edge_map[edge]])
                end
            else
                if length(edge_attributes) > 0
                    if (new_node1, new_node2) in keys(edge_bool_avg_index)
                        push!(edge_bool_avg_index[(new_node1, new_node2)], edge_map[edge])
                    else
                        edge_bool_avg_index[(new_node1, new_node2)] = Vector{T}([edge_map[edge]])
                    end
                end
            end
        end
    end

    if length(edge_attributes) > 0

        for edge in keys(edge_bool_avg_index)
            new_index = new_edge_map[edge]
            edge_data_to_avg = edge_data[edge_bool_avg_index[edge], :]

            edge_data_avg = Statistics.mean(edge_data_to_avg; dims = 1)
            new_edge_data[new_index, :] = edge_data_avg[:]
        end

        new_dg.edge_data.attributes    = edge_attributes
        new_dg.edge_data.attribute_map = edge_attirbute_map
        new_g.edge_data.data           = new_edge_data
    end

    simple_digraph = Graphs.SimpleDiGraph(T(length(new_edges)), fadjlist, badjlist)

    new_dg.g        = simple_digraph
    new_dg.nodes    = new_nodes
    new_dg.node_map = new_node_dict
    new_dg.edges    = new_edges
    new_dg.edge_map = new_edge_map

    return new_dg
end
