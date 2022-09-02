using Colors, TestImages, Images
include("../src/DataGraphs.jl")
img = Images.load("./examples/Bucky_Badger.jpg")
#picture from https://www.wikiwand.com/en/Bucky_Badger

imgg = Gray.(img)

mat = convert(Array{Float64}, imgg)

mat_graph = DataGraphs.matrix_to_graph(mat)
DataGraphs.plot_graph(mat_graph, color=:gray, C=1, K=.01, xdim = 400, ydim = 400)
