using Revise
using DataGraphs, Graphs
using JLD, LinearAlgebra
using Plots, Statistics

data = JLD.load("examples/brain.jld")["data"]
thresh = 0:.0002:.2

ECs = Array{Any,2}(undef, length(thresh), 30)

for i in 1:30
    mat = abs.(data[i,:,:]) - I
    h = symmetric_matrix_to_graph(mat[:,:])
    println("built graph")
    EC_vals = run_EC_on_edges(h, thresh)
    println("Running EC")
    ECs[:,i] .= EC_vals
    println(i)
end


for i in 1:30
    dg = DataGraph()
    datai = (abs.(data[i,:,:]) - I)

    for j in 1:39
        for k in 1:39
            add_edge!(dg, k, j)
            add_edge_data!(dg, k, j, datai[k,j], "weight")
        end
    end

    println("built graph")
    EC_vals = run_EC_on_edges(dg, thresh)
    ECs[:,i] .= EC_vals
end

adult = mean(ECs[:,1:6], dims=2)
child = mean(ECs[:,7:30], dims=2)

plt = plot(thresh, adult, label="Developed")
plot!(thresh, child, label="Underdeveloped")

include("plots.jl")

h = symmetric_matrix_to_graph(data[1,:,:])
x = DataGraph()
x.nodes = h.nodes
plot_graph(x)
h.node_positions = x.node_positions
plot_graph(h; color=:gray, linealpha=.2, xdim = 400, ydim = 400, save_fig=true, fig_name="full_plot.png")
new_h = filter_edges(h, .00005)
plot_graph(new_h;get_new_positions=false, color=:gray, linealpha=.2, markersize=10)

ah = aggregate(h, [1,2,3,4,5,6,7], "agg")

plot_graph(ah;save_pos=false, get_new_positions=false, color=:gray, linealpha=.2, markersize=10)
