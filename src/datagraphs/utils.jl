"""
    get_EC(datagraph)

Returns the Euler Characteristic for a DataGraph. The Euler Characteristic is equal to the
number of nodes minus the number of edges or the number of connected components minus the
number of cycles
"""
function get_EC(dg::DataGraph)
    nodes = dg.nodes
    edges = dg.edges

    EC = length(nodes) - length(edges)

    return EC
end

"""
    matrix_to_graph(matrix, diagonal = false, attribute="weight")

Constructs a `DataGraph` object from a matrix and saves the matrix data as node attributes under
the name `attribute`. If `diagonal = false`, the graph has a mesh structure, where each matrix entry is represented by
a node, and each node is connected to the adjacent matrix entries/nodes. If `diagonal = true`, entries of the matrix
are also connected to the nodes diagonal to them (i.e., entry (i,j) is connected to (i-1, j-1), (i + 1, j -1), etc.).
"""
function matrix_to_graph(matrix::AbstractMatrix, diagonal::Bool = false, attribute::String = "weight")

    dim1, dim2 = size(matrix)
    dg = DataGraph()

    fadjlist  = [Vector{Int}() for i in 1:(dim1 * dim2)]
    edges     = dg.edges
    nodes     = dg.nodes
    node_map  = dg.node_map
    edge_map  = dg.edge_map
    node_data = fill(NaN, (dim1*dim2, 1))

    for j in 1:dim2
        for i in 1:dim1
            push!(nodes, (i, j))
            node_map[(i, j)] = length(nodes)
            node_data[length(nodes), 1] = matrix[i, j]
        end
    end

    if dim1 > 1 && dim2 > 1
        for j in 1:dim2
            column_offset = dim1 * (j - 1)
            for i in 1:dim1
                if j != dim2
                    edge = (i + column_offset, i + column_offset + dim1)
                    push!(edges, edge)
                    push!(fadjlist[i + column_offset], i + column_offset + dim1)
                    push!(fadjlist[i + column_offset + dim1], i + column_offset)
                    edge_map[edge] = length(edges)
                end

                if i != dim1
                    edge = (i + column_offset, i + column_offset + 1)
                    push!(fadjlist[i + column_offset], i + column_offset + 1)
                    push!(fadjlist[i + column_offset + 1], i + column_offset)
                    push!(edges, edge)
                    edge_map[edge] = length(edges)
                end

                if !diagonal
                    # Add diagonal from top left to bottom right
                    if i != dim1 && j != dim2
                        edge = (i + column_offset, i + column_offset + dim1 + 1)
                        push!(fadjlist[i + column_offset], i + column_offset + dim1 + 1)
                        push!(fadjlist[i + column_offset + dim1 + 1], i + column_offset)
                        push!(edges, edge)
                        edge_map[edge] = length(edges)
                    end
                    # add diagonal from bottom left to top right
                    if i != 1 && j != dim2
                        edge = (i + column_offset, i + column_offset + dim1 - 1)
                        push!(fadjlist[i + column_offset], i + column_offset + dim1 - 1)
                        push!(fadjlist[i + column_offset + dim1 - 1], i + column_offset)
                        push!(edges, edge)
                        edge_map[edge] = length(edges)
                    end
                end
            end
        end
    end

    simple_graph = Graphs.SimpleGraph(length(edges), fadjlist)
    dg.node_data.attributes                 = [attribute]
    dg.node_data.attribute_map[attribute] = 1

    dg.g               = simple_graph
    dg.nodes           = nodes
    dg.node_map        = node_map
    dg.edges           = edges
    dg.edge_map        = edge_map
    dg.node_data.data  = node_data

    return dg
end

"""
    symmetric_matrix_to_graph(matrix; attribute="weight", tol = 1e-9)

Constructs a `DataGraph` object from a symmetric matrix and saves the values of the matrix to
their corresponding edges. The resulting graph is fully connected (every node is connected to every node)
and the number of nodes is equal to the dimension of the matrix. Matrix values are saved as edge
weights under the name `attribute`. `tol` is the tolerance used when testing that the matrix is symmetric.
"""
function symmetric_matrix_to_graph(matrix::M; attribute::String="weight", tol::R = 1e-9) where {M <: AbstractMatrix, R <: Real}

    dim1, dim2 = size(matrix)

    if dim1 != dim2
        error("Matrix is not square")
    end

    if !all(abs.(matrix - matrix') .<= tol)
        error("Matrix is not symmetric")
    end

    dg = DataGraph()

    for j in 1:dim2
        for i in (j + 1):dim1
            DataGraphs.add_edge!(dg, i, j)
            add_edge_data!(dg, i, j, matrix[i,j], attribute)
        end
    end

    return dg
end

"""
    mvts_to_graph(mvts, attribute)

Converts a multivariate time series to a graph based on the covariance matrix. This first
calculates the covariance of the multivariate time series (`mvts`) and then computes the
covariance. It then forms the precision matrix by taking the inverse of the covariance and
uses the function `symmetric_matrix_to_graph` to form the edge-weighted graph.
"""
function mvts_to_graph(mvts, attribute::String="weight", tol::R=1e-9) where {R <: Real}

    mvts_cov  = cov(mvts)
    mvts_prec = inv(mvts_cov)

    dg = symmetric_matrix_to_graph(mvts_prec, attribute, tol)

    return dg
end

"""
    tensor_to_graph(tensor, attribute)

Constructs a graph from a 3-D array (a tensor). Each entry of the tensor is represented by
a node, and each node is connected to the adjacent nodes in each dimension. This function
creates the graph structure and saves the values of the tensor to their corresponding nodes
as weight values under the name `attribute`.
"""
function tensor_to_graph(tensor::A, attribute::String="weight") where {A <: AbstractArray}

    if length(size(tensor)) != 3
        error("Tensor must have 3 dimensions; given tensor has $(length(size(tensor))) dimensions")
    end

    dim1, dim2, dim3 = size(tensor)
    dg = DataGraph()

    fadjlist  = [Vector{Int}() for i in 1:(dim1 * dim2 * dim3)]
    edges     = dg.edges
    nodes     = dg.nodes
    node_map  = dg.node_map
    edge_map  = dg.edge_map
    node_data = Array{eltype(tensor), 2}(undef, (dim1 * dim2 * dim3, 1)) #fill(NaN, (dim1*dim2, 1))

    for k in 1:dim3
        for j in 1:dim2
            for i in 1:dim1
                push!(nodes, (i, j, k))
                node_map[(i, j, k)] = length(nodes)
                node_data[length(nodes), 1] = tensor[i, j, k]
            end
        end
    end

    if dim1 > 1 && dim2 > 1 && dim3 > 1
        for k in 1:dim3
            dim3_offset = dim1 * dim2 * (k - 1)
            for j in 1:dim2
                column_offset = dim1 * (j - 1)
                for i in 1:dim1
                    if j != dim2
                        edge = (i + column_offset + dim3_offset, i + column_offset + dim3_offset + dim1)
                        push!(edges, edge)
                        push!(fadjlist[i + column_offset + dim3_offset], i + column_offset + dim3_offset + dim1)
                        push!(fadjlist[i + column_offset + dim3_offset + dim1], i + column_offset + dim3_offset)
                        edge_map[edge] = length(edges)
                    end

                    if i != dim1
                        edge = (i + column_offset + dim3_offset, i + column_offset + dim3_offset + 1)
                        push!(fadjlist[i + column_offset + dim3_offset], i + column_offset + dim3_offset + 1)
                        push!(fadjlist[i + column_offset + dim3_offset + 1], i + column_offset + dim3_offset)
                        push!(edges, edge)
                        edge_map[edge] = length(edges)
                    end

                    if k != dim3
                        edge = (i + column_offset + dim3_offset, i + column_offset + dim3_offset + dim1 * dim2)
                        push!(fadjlist[i + column_offset + dim3_offset], i + column_offset + dim3_offset + dim1 * dim2)
                        push!(fadjlist[i + column_offset + dim3_offset + dim1 * dim2], i + column_offset + dim3_offset)
                        push!(edges, edge)
                        edge_map[edge] = length(edges)
                    end
                end
            end
        end
    end

    simple_graph = Graphs.SimpleGraph(length(edges), fadjlist)
    dg.node_data.attributes                 = [attribute]
    dg.node_data.attribute_map[attribute] = 1

    dg.g               = simple_graph
    dg.nodes           = nodes
    dg.node_map        = node_map
    dg.edges           = edges
    dg.edge_map        = edge_map
    dg.node_data.data  = node_data

    return dg
end

"""
    filter_nodes(datagraph, filter_value; attribute_name)

Removes the nodes of the graph whose weight value of `attribute_name` is greater than the given
`filter_value`. If `attribute_name` is not specified, this defaults to the first attribute within
the DataGraph's `NodeData`.
"""
function filter_nodes(dg::DataGraph, filter_val::R; attribute::String=dg.node_data.attributes[1]) where {R <: Real}
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

    new_dg = DataGraph()

    am = adjacency_matrix(dg)

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
    fadjlist       = [Vector{T}() for i in 1:length(new_nodes)]  ### TODO: if new_nodes is length 0, this is a vector of type any

    for i in 1:length(new_nodes)
        new_node_map[new_nodes[i]] = i
    end

    for j in 1:length(new_nodes)
        for i in (j + 1):length(new_nodes)
            if new_am[i,j] == 1
                new_edge = (j, i)
                push!(new_edges, new_edge)
                new_edge_map[new_edge] = length(new_edges)

                @inbounds node_neighbors = fadjlist[i]
                index = searchsortedfirst(node_neighbors, j)
                insert!(node_neighbors, index, j)

                @inbounds node_neighbors = fadjlist[j]
                index = searchsortedfirst(node_neighbors, i)
                insert!(node_neighbors, index, i)

                old_edge = _get_edge(node_map[new_nodes[j]], node_map[new_nodes[i]])
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

    simple_graph = Graphs.SimpleGraph(T(length(new_edges)), fadjlist)

    new_dg.g                    = simple_graph
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
    filter_edges(datagraph, filter_value; attribute_name)

Removes the edges of the graph whose weight value of `attribute_name` is greater than the given
`filter_value`. If `attribute_name` is not specified, this defaults to the first attribute within
the DataGraph's `EdgeData`.
"""
function filter_edges(dg::DataGraph, filter_val::R; attribute::String = dg.edge_data.attributes[1]) where {R <: Real}
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

    bool_vec = dg.edge_data.data[:, edge_attribute_map[attribute]] .< filter_val

    new_edges = edges[bool_vec]
    new_edge_data = edge_data[bool_vec, :]

    new_edge_map = Dict{Tuple{T, T}, T}()

    fadjlist = [Vector{T}() for i in 1:length(nodes)]

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

    new_dg = DataGraph()

    simple_graph = Graphs.SimpleGraph(T(length(new_edges)), fadjlist)

    new_dg.g                    = simple_graph
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
    run_EC_on_nodes(datagraph, threshold_range; attribute_name, scale = false)

Returns the Euler Characteristic Curve by filtering the nodes of the graph at each value in `threshold_range`
and computing the Euler Characteristic after each filtration. If `attribute_name` is not defined, it defaults
to the first attribute in the DataGraph's `NodeData`. `scale` is a Boolean that indicates whether to scale
the Euler Characteristic by the total number of objects (nodes + edges) in the original graph
"""
function run_EC_on_nodes(dg::DataGraph, thresh; attribute::String = dg.node_data.attributes[1], scale::Bool = false)
    nodes        = dg.nodes
    node_data    = dg.node_data.data

    node_attribute_map = dg.node_data.attribute_map

    am = Graphs.LinAlg.adjacency_matrix(dg.g)

    if scale
        scale_factor = length(dg.nodes) + length(dg.edges)
    else
        scale_factor = 1
    end

    for i in 1:length(nodes)
        if am[i, i] == 1
            am[i, i] = 2
        end
    end

    ECs = zeros(length(thresh))

    for (j,i) in enumerate(thresh)
        bool_vec  = node_data[:, node_attribute_map[attribute]] .< i
        new_am    = am[bool_vec, bool_vec]
        num_nodes = sum(bool_vec)
        num_edges = sum(new_am.nzval) / 2
        ECs[j]    = num_nodes - num_edges
    end

    return ECs ./ scale_factor
end

"""
    run_EC_on_edges(datagraph, threshold_range; attribute_name, scale = false)

Returns the Euler Characteristic Curve by filtering the edges of the graph at each value in `threshold_range`
and computing the Euler Characteristic after each filtration. If `attribute_name` is not defined, it defaults
to the first attribute in the DataGraph's `EdgeData`. `scale` is a Boolean that indicates whether to scale
the Euler Characteristic by the total number of objects (nodes + edges) in the original graph
"""
function run_EC_on_edges(dg::DataGraph, thresh; attribute::String = dg.edge_data.attributes[1], scale::Bool = false)
    edge_data = dg.edge_data.data
    nodes     = dg.nodes
    edge_attribute_map = dg.edge_data.attribute_map

    if scale
        scale_factor = length(dg.nodes) + length(dg.edges)
    else
        scale_factor = 1
    end

    ECs = zeros(length(thresh))

    num_nodes = length(nodes)

    for (j,i) in enumerate(thresh)
        bool_vec  = edge_data[:, edge_attribute_map[attribute]].< i
        num_edges = sum(bool_vec)
        ECs[j]    = num_nodes - num_edges
    end

    return ECs ./ scale_factor
end

"""
    aggregate(datagraph, node_list, aggregated_node_name)

Aggregates all the nodes in `node_list` into a single node which is called `aggregated_node_name`.
If nodes have any weight/attribute values defined, These are averaged across all values in the
`node_list`. Edge weights are also averaged when two or more nodes in the `node_list` are connected
to the same node and these edges have weights defined on them.
"""
function aggregate(dg::DataGraph, node_set, new_name)
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

    new_dg = DataGraph()

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

    fadjlist = [Vector{T}() for i in 1:length(new_nodes)]

    node_name_mapping   = Dict{T, Any}()
    new_edges           = Vector{Tuple{T, T}}()
    new_edge_map        = Dict{Tuple{T, T}, T}()
    edge_bool_vec       = [false for i in 1:length(edges)]
    edge_bool_avg_index = Dict{Tuple{T, T}, Vector{T}}()
    new_edge_data       = fill(NaN, (0, length(edge_attributes)))

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

            @inbounds node_neighbors = fadjlist[new_node2]
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

                @inbounds node_neighbors = fadjlist[new_node2]
                index = searchsortedfirst(node_neighbors, new_node1)
                insert!(node_neighbors, index, new_node1)

                if length(edge_attributes) > 0
                    new_edge_data = vcat(new_edge_data, fill(NaN, (1, length(edge_attributes))))
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

                @inbounds node_neighbors = fadjlist[new_node2]
                index = searchsortedfirst(node_neighbors, new_node1)
                insert!(node_neighbors, index, new_node1)

                if length(edge_attributes) > 0
                    new_edge_data = vcat(new_edge_data, fill(NaN, (1, length(edge_attributes))))
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

    simple_graph = Graphs.SimpleGraph(T(length(new_edges)), fadjlist)

    new_dg.g        = simple_graph
    new_dg.nodes    = new_nodes
    new_dg.node_map = new_node_dict
    new_dg.edges    = new_edges
    new_dg.edge_map = new_edge_map

    return new_dg
end
