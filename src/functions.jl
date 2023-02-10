function Graphs.connected_components(dg::D) where {D <: DataGraphUnion}
    connected_components_list = Graphs.connected_components(dg.g)

    nodes  = dg.nodes
    components = []

    for i in connected_components_list
        push!(components, nodes[i])
    end

    return components
end

function Graphs.is_connected(dg::D) where {D <: DataGraphUnion}
    return Graphs.is_connected(dg.g)
end

function Graphs.common_neighbors(dg::D, node1, node2) where {D <: DataGraphUnion}
    node_map = dg.node_map
    return Graphs.common_neighbors(dg.g, node_map[node1], node_map[node2])
end

function Graphs.neighbors(dg::D, node) where {D <: DataGraphUnion}
    node_map = dg.node_map
    return Graphs.neighbors(dg.g, node_map[node])
end

function Graphs.core_number(dg::D) where {D <: DataGraphUnion}
    return Graphs.core_number(dg.g)
end

function Graphs.k_core(dg::D, k = -1) where {D <: DataGraphUnion}
    k_core_index = Graphs.k_core(dg.g, k)
    return dg.nodes[k_core_index]
end

function Graphs.k_shell(dg::D, k = -1) where {D <: DataGraphUnion}
    k_shell_index = Graphs.k_shell(dg.g, k)
    return dg.nodes[k_shell_index]
end

function Graphs.k_crust(dg::D, k = -1) where {D <: DataGraphUnion}
    k_crust_index = Graphs.k_crust(dg.g, k)
    return dg.nodes[k_crust_index]
end

function Graphs.eccentricity(dg::D, node) where {D <: DataGraphUnion}
    node_map = dg.node_map
    return Graphs.eccentricity(dg.g, node_map[node])
end

function Graphs.diameter(dg::D) where {D <: DataGraphUnion}
    return Graphs.diameter(dg.g)
end

function Graphs.periphery(dg::D) where {D <: DataGraphUnion}
    return Graphs.periphery(dg.g)
end

function Graphs.radius(dg::D) where {D <: DataGraphUnion}
    return Graphs.radius(dg.g)
end

function Graphs.center(dg::D) where {D <: DataGraphUnion}
    return Graphs.center(dg.g)
end

function Graphs.complement(dg::D) where {D <: DataGraphUnion}
    return Graphs.complement(dg.g)
end

function Graphs.cycle_basis(dg::D) where {D <: DataGraphUnion}
    numbered_cycles = Graphs.cycle_basis(dg.g)

    nodes  = dg.nodes
    cycles = []

    for i in numbered_cycles
        push!(cycles, nodes[i])
    end

    return cycles
end

function Graphs.indegree(dg::D) where {D <: DataGraphUnion}
    return Graphs.indegree(dg.g)
end

function Graphs.indegree(dg::D, node) where {D <: DataGraphUnion}
    node_map = dg.node_map
    return Graphs.indegree(dg, node_map[node])
end

function Graphs.outdegree(dg::D) where {D <: DataGraphUnion}
    return Graphs.outdegree(dg.g)
end

function Graphs.outdegree(dg::D, node) where {D <: DataGraphUnion}
    node_map = dg.node_map
    return Graphs.outdegree(dg, node_map[node])
end

function Graphs.degree(dg::D) where {D <: DataGraphUnion}
    return Graphs.degree(dg.g)
end

function Graphs.degree(dg::D, node) where {D <: DataGraphUnion}
    node_map = dg.node_map
    return Graphs.degree(dg, node_map[node])
end

function Graphs.degree_histogram(dg::D) where {D <: DataGraphUnion}
    return degree_histogram(dg.g)
end

function Graphs.degree_centrality(dg::D) where {D <: DataGraphUnion}
    return degree_centrality(dg.g)
end

"""
    average_degree(datagraph)

Returns the average degree for `datagraph`
"""
function average_degree(dg::D) where {D <: DataGraphUnion}
    degrees = Graphs.degree(dg)
    return sum(degrees) / length(degrees)
end

"""
    Graphs.has_path(datagraph, src_node, dst_node)

Returns true if a path exists in the `datagraph` between `src_node` to `dst_node`. Else returns false
"""
function has_path(dg::D, src_node::Any, dst_node::Any) where {D <: DataGraphUnion}
    node_map = dg.node_map
    nodes    = dg.nodes

    if !(src_node in nodes && dst_node in nodes)
        error("User has passed a node that does not exist in the DataGraph")
    end

    src_index = node_map[src_node]
    dst_index = node_map[dst_node]

    return Graphs.has_path(dg.g, src_index, dst_index)
end

"""
    Graphs.has_path(datagraph, src_node, intermediate_node, dst_node)

Returns true if a path exists in the `datagraph` between `src_node` and `dst_node` which
passes through the `intermediate node`. Else returns false
"""
function has_path(dg::D, src_node, intermediate_node, dst_node) where {D <: DataGraphUnion}
    node_map = dg.node_map
    nodes    = dg.nodes

    if !(src_node in nodes && intermediate_node in nodes && dst_node in nodes)
        error("User has passed a node that does not exist in the DataGraph")
    end

    src_index = node_map[src_node]
    int_index = node_map[intermediate_node]
    dst_index = node_map[dst_node]

    if Graphs.has_path(dg.g, src_index, int_index) && Graphs.has_path(dg.g, int_index, dst_index)
        return true
    else
        return false
    end
end

"""
    get_path(datagraph, src_node, dst_node; algorithm = "Dijkstra")

Returns the shortest path in the `datagraph` between `src_node` and `dst_node`.
Shortest path is computed by Dijkstra's algorithm

`algorithm` is a string key word. Options are limited to "Dijkstra", "BellmanFord"
"""
function get_path(dg::D, src_node, dst_node; algorithm = "Dijkstra") where {D <: DataGraphUnion}
    node_map = dg.node_map
    nodes    = dg.nodes

    if !(src_node in nodes && dst_node in nodes)
        error("User has passed a node that does not exist in the DataGraph")
    end

    src_index = node_map[src_node]
    dst_index = node_map[dst_node]

    if algorithm == "Dijkstra"
        path_state = Graphs.dijkstra_shortest_paths(dg.g, [src_index])
    elseif algorithm == "BellmanFord"
        path_state = Graphs.bellman_ford_shortest_paths(dg.g, [src_index])
    else
        error("$algorithm is not a supported algorithm option")
    end

    index_path = Graphs.enumerate_paths(path_state, dst_index)

    if length(index_path) == 0
        println("Path between nodes does not exist")
        return []
    end

    path = Vector{Any}(undef, length(index_path))

    for i in 1:length(index_path)
        path[i] = nodes[index_path[i]]
    end

    return path
end

"""
    get_path_with_intermediate(datagraph, src_node, intermediate_node, dst_node; algorithm = "Dijkstra")

Returns the shortest path in the `datagraph` between `src_node` and `dst_node`
which passes through `intermediate node`.

`algorithm` is a string key word. Options are limited to "Dijkstra", "BellmanFord"
"""
function get_path_with_intermediate(
    dg::D,
    src_node::Any,
    intermediate_node::Any,
    dst_node::Any;
    algorithm = "Dijkstra"
) where {D <: DataGraphUnion}

    node_map = dg.node_map
    nodes    = dg.nodes

    if !(src_node in nodes && intermediate_node in nodes && dst_node in nodes)
        error("User has passed a node that does not exist in the DataGraph")
    end

    src_index = node_map[src_node]
    int_index = node_map[intermediate_node]
    dst_index = node_map[dst_node]

    if algorithm == "Dijkstra"
        path_state_src = Graphs.dijkstra_shortest_paths(dg.g, [src_index])
        path_state_int = Graphs.dijkstra_shortest_paths(dg.g, [int_index])
    elseif algorithm == "BellmanFord"
        path_state_src = Graphs.bellman_ford_shortest_paths(dg.g, [src_index])
        path_state_int = Graphs.bellman_ford_shortest_paths(dg.g, [int_index])
    else
        error("$algorithm is not a supported algorithm option")
    end

    index_path_to_int = Graphs.enumerate_paths(path_state_src, int_index)

    index_path_to_dst = Graphs.enumerate_paths(path_state_int, dst_index)

    if length(index_path_to_int) == 0 || length(index_path_to_dst) == 0
        println("Path through intermediate node does not exist")
        return []
    end

    path = Vector{Any}(undef, (length(index_path_to_int) + length(index_path_to_dst) - 1))
    nodes = dg.nodes

    path_to_int_len = length(index_path_to_int)
    path_to_dst_len = length(index_path_to_dst)

    for i in 1:path_to_int_len
        path[i] = nodes[index_path_to_int[i]]
    end

    for i in 2:path_to_dst_len
        path[i + path_to_int_len - 1] = nodes[index_path_to_dst[i]]
    end

    return path
end
