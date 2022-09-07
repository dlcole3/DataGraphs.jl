using Revise
using Colors, TestImages, Images
using DataGraphs

include("plots.jl")

img = Images.load("./examples/Bucky_Badger.jpg")
#picture from https://www.wikiwand.com/en/Bucky_Badger

imgg = Gray.(img)

mat = convert(Array{Float64}, imgg)

mat_graph = DataGraphs.matrix_to_graph(mat)
mat_graph.node_positions = DataGraphs.set_matrix_node_positions(mat_graph.nodes, mat)

DataGraphs.plot_graph(mat_graph; plot_edges=false, markersize = .5, save_fig = true, fig_name = "full_graph.png")

for i in 1:8
    @time filtered_mat_graph = DataGraphs.filter_nodes(mat_graph, i*.125; attribute = "weight")
    println("done with filter on $i")
    DataGraphs.plot_graph(filtered_mat_graph, plot_edges = false, markersize = .5, save_fig = true, fig_name = "graph_$i.png")
    println("Done with $i")
end

mat3 = channelview(img)

mat_graph_R = DataGraphs.matrix_to_graph(mat3[1, :, :])
mat_graph_R.node_positions = DataGraphs.set_matrix_node_positions(mat_graph_R.nodes, mat3[1, :, :])

mat_graph_G = DataGraphs.matrix_to_graph(mat3[2, :, :])
mat_graph_G.node_positions = DataGraphs.set_matrix_node_positions(mat_graph_G.nodes, mat3[2, :, :])

mat_graph_B = DataGraphs.matrix_to_graph(mat3[3, :, :])
mat_graph_B.node_positions = DataGraphs.set_matrix_node_positions(mat_graph_B.nodes, mat3[3, :, :])

for i in 1:8
    @time filtered_mat_graph = DataGraphs.filter_nodes(mat_graph_R, i*.125; attribute = "weight")
    println("done with filter on $i")
    DataGraphs.plot_graph(filtered_mat_graph, plot_edges = false, markersize = .5, save_fig = true, fig_name = "graph_R$i.png")
    println("Done with $i")
end

for i in 1:8
    @time filtered_mat_graph = DataGraphs.filter_nodes(mat_graph_G, i*.125; attribute = "weight")
    println("done with filter on $i")
    DataGraphs.plot_graph(filtered_mat_graph, plot_edges = false, markersize = .5, save_fig = true, fig_name = "graph_G$i.png")
    println("Done with $i")
end

for i in 1:8
    @time filtered_mat_graph = DataGraphs.filter_nodes(mat_graph_B, i*.125; attribute = "weight")
    println("done with filter on $i")
    DataGraphs.plot_graph(filtered_mat_graph, plot_edges = false, markersize = .5, save_fig = true, fig_name = "graph_B$i.png")
    println("Done with $i")
end

#DataGraphs.plot_graph(mat_graph, color=:gray, C=1, K=.01, xdim = 400, ydim = 400)
