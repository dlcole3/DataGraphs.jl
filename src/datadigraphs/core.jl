

function DataDiGraph{T, T1, T2, M1, M2}() where {T <: Integer, T1, T2,  M1 <: Matrix{T1}, M2 <: Matrix{T2}}
    nodes = Vector{Any}()
    edges = Vector{Tuple{Int, Int}}()

    ne = 0
    fadjlist = Vector{Vector{Int}}()
    badjlist = Vector{Vector{Int}}()

    node_map = Dict{Any, Int}()
    edge_map = Dict{Tuple{Int, Int}, Int}()
    node_attributes = String[]
    edge_attributes = String[]
    node_attribute_map = Dict{String, Int}()
    edge_attribute_map = Dict{String, Int}()
    node_data = Array{Float64}(undef, 0, 0)
    edge_data = Array{Float64}(undef, 0, 0)

    node_positions = [[0.0 0.0]]

    g = SimpleDiGraph(ne, fadjlist, badjlist)

    node_data_struct = NodeData(node_attributes, node_attribute_map, node_data)
    edge_data_struct = EdgeData(edge_attributes, edge_attribute_map, edge_data)

    DataDiGraph{T, T1, T2, M1, M2}(
        g, nodes, edges, node_map, edge_map,
        node_data_struct, edge_data_struct, node_positions
    )
end

DataDiGraph() = DataDiGraph{Int, Float64, Float64, Matrix{Float64}, Matrix{Float64}}()

function DataDiGraph(adj_mat::AbstractMatrix{T}) where {T <: Real}

    dima, dimb = size(adj_mat)
    isequal(dima, dimb) || throw(ArgumentError("Adjacency / distance matrices must be square"))
    LinearAlgebra.issymmetric(adj_mat) || throw(ArgumentError("Adjacency / distance matrices must be symmetric"))

    dg = DataDiGraph()

    maxc = length(adj_mat.colptr)
    @inbounds for c = 1:(maxc - 1)
        for rind = adj_mat.colptr[c]:(adj_mat.colptr[c + 1] - 1)
            isnz = (adj_mat.nzval[rind] != zero(T))
            if isnz
                r = adj_mat.rowval[rind]
                DataGraphs.add_edge!(dg, r, c)
            end
        end
    end
    return dg
end

"""
    add_node!(dg, node_name)

Add the node `node_name` to the DataDiGraph `dg`
"""
function add_node!(
    dg::DataDiGraph, node_name::Any
)
    nodes      = dg.nodes
    attributes = dg.node_data.attributes
    node_map   = dg.node_map

    T = eltype(dg)

    # If new node is not in the list of nodes, add it
    # otherwise, print that the node exists and don't do anything
    if !(node_name in nodes)
        push!(nodes,node_name)
        push!(dg.g.fadjlist, Vector{T}())
        push!(dg.g.badjlist, Vector{T}())

        # If there are data currently defined on the other nodes, add a NaN value to
        # the end of the weight array for the new node
        if length(attributes)>0
            node_data = dg.node_data.data
            row_to_add = fill(NaN, (1, length(attributes)))
            node_data = vcat(node_data, row_to_add)
            dg.node_data.data = node_data
        end

        # Add the new node as a key to the dictionary
        node_map[node_name] = length(nodes)
        dg.node_map = node_map
        return true
    else
       println("Node already exists")
       return false
    end
end


"""
    add_edge!(dg, node_1, node_2)
    add_edge!(dg, (node1, node2))

Add an edge to the DataDiGraph, `dg`. If the nodes are not defined in the graph, they are added to the graph
"""
function add_edge!(dg::DataDiGraph, node1::Any, node2::Any)
    edges      = dg.edges
    nodes      = dg.nodes
    attributes = dg.edge_data.attributes
    edge_map   = dg.edge_map

    if !(node1 in nodes)
        add_node!(dg, node1)
    end
    if !(node2 in nodes)
        add_node!(dg, node2)
    end

    nodes       = dg.nodes
    node_map    = dg.node_map

    node1_index = node_map[node1]
    node2_index = node_map[node2]

    edge = (node1_index, node2_index)

    # If the edge isn't already defined, then add the edge; add to weight arrays too
    if !(edge in edges)
        push!(edges, edge)
        dg.g.ne += 1

        @inbounds node_neighbors = dg.g.fadjlist[node1_index]
        index = searchsortedfirst(node_neighbors, node2_index)
        insert!(node_neighbors, index, node2_index)

        @inbounds node_neighbors = dg.g.badjlist[node2_index]
        index = searchsortedfirst(node_neighbors, node1_index)
        insert!(node_neighbors, index, node1_index)

        if length(attributes)>0
            edge_data  = dg.edge_data.data
            row_to_add = fill(NaN, (1, length(attributes)))
            edge_data  = vcat(edge_data, row_to_add)
            dg.edge_data.data = edge_data
        end

        edge_map[edge] = length(edges)
        return true
    else
        return false
    end
end

function add_edge!(dg::DataDiGraph, edge::Tuple{Any, Any})
    DataGraphs.add_edge!(dg, edge[1], edge[2])
end


function add_node_data!(dg::DataDiGraph, node::Any, node_weight::Number, attribute::String)
    nodes         = dg.nodes
    attributes    = dg.node_data.attributes
    node_map      = dg.node_map
    node_data     = dg.node_data.data
    attribute_map = dg.node_data.attribute_map

    if !(node in nodes)
        error("node does not exist in graph")
    end

    if length(attributes) < 1
        node_data = Array{eltype(dg.node_data.data)}(undef, length(nodes), 0)
        dg.node_data.data = node_data
    end

    if !(attribute in attributes)
        # Add new column to node_weight array
        push!(attributes, attribute)
        attribute_map[attribute] = length(attributes)
        new_col = fill(NaN, (length(nodes), 1))
        node_data = hcat(node_data, new_col)
        node_data[node_map[node], attribute_map[attribute]] = node_weight
        dg.node_data.data = node_data
        return true
    else
        node_data[node_map[node], attribute_map[attribute]] = node_weight
        return true
    end
end

function add_edge_data!(dg::DataDiGraph, node1::Any, node2::Any, edge_weight::T, attribute::String) where {T <: Real}
    edges         = dg.edges
    attributes    = dg.edge_data.attributes
    edge_map      = dg.edge_map
    node_map      = dg.node_map
    attribute_map = dg.edge_data.attribute_map

    node1_index = node_map[node1]
    node2_index = node_map[node2]

    edge = (node1_index, node2_index)

    if !(edge in edges)
        error("edge does not exist in graph")
    end

    if length(attributes) == 0
        edge_data = Array{eltype(dg.edge_data.data)}(undef, length(edges), 0)
        dg.edge_data.data = edge_data
    end

    if !(attribute in attributes)
        edge_data = dg.edge_data.data
        # Add new column to node_weight array
        push!(attributes, attribute)
        attribute_map[attribute] = length(attributes)
        new_col = fill(NaN, (length(edges), 1))
        edge_data = hcat(edge_data, new_col)
        edge_data[edge_map[edge], attribute_map[attribute]] = edge_weight
        dg.edge_data.data = edge_data
        return true
    else
        edge_data = dg.edge_data.data
        edge_data[edge_map[edge], attribute_map[attribute]] = edge_weight
        return true
    end
end

function add_edge_data!(dg::DataDiGraph, edge::Tuple{Any, Any}, edge_weight::T, attribute::String) where {T <: Real}
    add_edge_data!(dg, edge[1], edge[2], edge_weight, attribute)
end
