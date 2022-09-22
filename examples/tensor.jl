using Revise
using DataGraphs, Graphs

abc = rand(10, 5, 4)

tensor_graph = tensor_to_graph(abc)

tensor_graph.edges

include("plots.jl")

plot_graph(tensor_graph)

xyz = rand(128, 48, 48)

tensor_graph1 = tensor_to_graph(abc)
@time tensor_graph1 = tensor_to_graph(abc)
