abstract type AbstractDataGraph{T} <: Graphs.AbstractGraph{T} end

mutable struct DataGraph{T} <: AbstractDataGraph{T}
    ne::Int
    fadjlist::Vector{Vector{T}}

    nodes::Vector{Any}
    edges::Vector{Tuple{T, T}}
    nodes_index::Dict{Any, T}
    edges_index::Dict{Tuple{T, T}, T}

    node_attributes::Vector{String}
    edge_attributes::Vector{String}

    node_data::NamedArray
    edge_data::NamedArray
    node_positions::Array{Union{GeometryBasics.Point{2,Float64}, Array{Float64, 2}},1}
end

function Base.eltype(datagraph::DataGraph)
    return eltype(eltype(datagraph.fadjlist))
end

function DataGraph(nodes::Vector{Any} = Vector{Any}(), edges::Vector{Tuple{T, T}} = Vector{Tuple{Int, Int}}(),
    ne::Int = length(edges),
    fadjlist::Vector{Vector{T}} = [Vector{T} for i in 1:length(nodes)],
    node_attributes::Vector{String} = String[],
    edge_attributes::Vector{String} = String[],
    nodes_index::Dict{Any, Int} = Dict{Any, Int}(),
    edges_index::Dict{Tuple{T}, Int} = Dict{Any, Int}(),
    node_data::NamedArray = [],
    edge_data::NamedArray = [],
    node_positions = [[0.0 0.0]]
) where T <: Int

    if length(edges) != ne
        error("Defined edges do not match ne")
    end
    if ne != length(edges_index)
        error("edges_index does not match the number of edges")
    end
    if length(nodes) != length(nodes_index)
        error("nodes_index does not match the number of nodes")
    end

    DataGraph{T}(
        ne, fadjlist, nodes, edges, nodes_index, edges_index,
        node_attributes, edge_attributes, node_data, edge_data,
        node_positions
    )
end

function DataGraph()
    nodes = Vector{Any}()
    edges = Vector{Tuple{Int, Int}}()

    ne = 0
    fadjlist = Vector{Vector{Int}}()

    node_attributes = String[]
    edge_attributes = String[]
    nodes_index = Dict{Any, Int}()
    edges_index = Dict{Tuple{Int, Int}, Int}()
    node_data = []
    edge_data = []
    node_positions = [[0.0 0.0]]

    DataGraph{Int}(
        ne, fadjlist, nodes, edges, nodes_index, edges_index,
        node_attributes, edge_attributes, node_data, edge_data,
        node_positions
    )
end

function _get_edge(node1_index, node2_index)
    if node2_index > node1_index
        return (node1_index, node2_index)
    else
        return (node2_index, node1_index)
    end
end


"""
    add_node!(g, node_name)

Add the node `node_name` to the graph `g`
"""
function add_node!(g::DataGraph,node_name::Any)
    nodes       = g.nodes
    attributes  = g.node_attributes
    nodes_index = g.nodes_index

    T = eltype(g)

    # If new node is not in the list of nodes, add it
    # otherwise, print that the node exists and don't do anything
    if !(node_name in nodes)
        push!(nodes,node_name)
        push!(g.fadjlist, Vector{T}())

        # If there are data currently defined on the other nodes, add a NaN value to
        # the end of the weight array for the new node
        if length(attributes)>0
            node_data = g.node_data
            row_to_add = NamedArray(fill(NaN, (1, length(attributes))))
            node_data = vcat(node_data, row_to_add)
            setnames!(node_data, attributes, 2)
            g.node_data = node_data
        end

        # Add the new node as a key to the dictionary
        nodes_index[node_name] = length(nodes)
        g.nodes_index = nodes_index
        return true

    else
       println("Node already exists")
       return false
    end
end



"""
    add_edge!(g, node_1, node_2)

Add an edge to the graph, `g`. If the nodes are not defined in the graph, they are added to the graph
"""
function add_edge!(g::DataGraph, node1::Any, node2::Any)
    # TODO: do things differently if edge attributes are defined or not defined;
    edges       = g.edges
    nodes       = g.nodes
    attributes  = g.edge_attributes
    edges_index = g.edges_index

    if !(node1 in nodes)
        add_node!(g, node1)
    end
    if !(node2 in nodes)
        add_node!(g, node2)
    end

    nodes       = g.nodes
    nodes_index = g.nodes_index

    node1_index = nodes_index[node1]
    node2_index = nodes_index[node2]

    edge = _get_edge(node1_index, node2_index)

    # If the edge isn't already defined, then add the edge; add to weight arrays too
    if !(edge in edges)
        push!(edges, edge)
        g.ne += 1

        @inbounds node_neighbors = g.fadjlist[node1_index]
        index = searchsortedfirst(node_neighbors, node2_index)
        insert!(node_neighbors, index, node2_index)

        @inbounds node_neighbors = g.fadjlist[node2_index]
        index = searchsortedfirst(node_neighbors, node1_index)
        insert!(node_neighbors, index, node1_index)


        if length(attributes)>0
            edge_data = g.edge_data
            row_to_add = NamedArray(fill(NaN, (1,length(attributes))))
            edge_data = vcat(edge_data, row_to_add)
            setnames!(edge_data, attributes, 2)
            g.edge_data = edge_data
        end

        edges_index[edge] = length(edges)
        return true
    else
        return false
    end
end

function add_node_data!(g::DataGraph, node::Any, node_weight::Number, attribute::String)
    nodes   = g.nodes
    attributes   = g.node_attributes
    nodes_index  = g.nodes_index
    node_data = g.node_data

    if !(node in nodes)
        error("node does not exist in graph")
    end

    if length(attributes) < 1
        node_data = NamedArray(fill(NaN, (length(nodes), 1)), (1:length(nodes),[attribute]))
        g.node_data = node_data
        push!(attributes, attribute)
    end

    if !(attribute in attributes)
        # Add new column to node_weight array
        push!(attributes, attribute)
        new_col = NamedArray(fill(NaN, length(nodes)), nodes)
        node_data = hcat(node_data, new_col)
        setnames!(node_data, attributes, 2)
        node_data[nodes_index[node], attribute] = node_weight

    else
        node_data[nodes_index[node], attribute] = node_weight
        return true
    end
end


function add_edge_data!(g::DataGraph, node1::Any, node2::Any, edge_weight::Real, attribute::String)
    edge_array   = g.edges
    attributes   = g.edge_attributes
    edges_index  = g.edges_index
    nodes_index  = g.nodes_index

    node1_index = nodes_index[node1]
    node2_index = nodes_index[node2]

    edge = _get_edge(node1_index, node2_index)

    if !(edge in edge_array)
        error("edge does not exist in graph")
    end

    if length(attributes) == 0
        new_array = NamedArray(fill(NaN, (length(edge_array), 1)), (1:length(edge_array),[attribute]))
        g.edge_data = new_array
        push!(attributes, attribute)
    end

    if !(attribute in attributes)
        edge_data = g.edge_data
        # Add new column to node_weight array
        push!(attributes, attribute)
        new_col = NamedArray(fill(NaN, length(edge_array)), edge_array)
        edge_data = hcat(edge_data, new_col)
        setnames!(edge_data, attributes, 2)
        edge_data[edges_index[edge], attribute] = edge_weight
    else
        edge_data = g.edge_data
        edge_data[edges_index[edge], attribute] = edge_weight
    end

end

function create_adj_mat(g::DataGraph; sparse::Bool = true)
    nodes       = g.nodes
    edges       = g.edges
    nodes_index = g.nodes_index

    nn = length(nodes)

    if sparse
        mat = SparseArrays.spzeros(Bool, nn, nn)
    else
        mat = zeros(Bool, nn, nn)
    end

    for i in edges
       mat[i[1], i[2]] = 1
       mat[i[2], i[1]] = 1
    end

    return mat
end
