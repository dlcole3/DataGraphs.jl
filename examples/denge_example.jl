using Revise
using DataGraphs, Graphs, DelimitedFiles, Statistics
using DataGraphPlots

# The data and methods in this example are primarily from the paper of
# Barros de Souze at al., 2022. https://doi.org/10.1088/1742-5468/aca0e5
# "The Euler characteristic as a topological marker for outbreaks in
# vector borne disease"

data_any = readdlm((@__DIR__)*"/Recife_data/Dengue_Recife_new_cases_jan_2014_to_dec_2021.csv", ',', skipstart = 1)[:, 2:end]

data = Matrix{Float64}(undef, size(data_any))

data .= data_any

cor_1 = cor(data[1:30, :], dims = 1)

bit_vec = (!).(isnan.(cor_1[:, 1]))

cor1b = cor_1[bit_vec, bit_vec]

cor1g = symmetric_matrix_to_graph(cor1b)

function find_smallest_filter_graph(graph, filter_value)
    all_comps = length(graph.nodes) + length(graph.edges)
    for i in maximum(graph.edge_data.data[:, 1]):(-filter_value):minimum(graph.edge_data.data[:, 1])
        fg = filter_edges(graph, i)
        if length(connected_components(fg)) > 1
            global f_val = i + filter_value
            fg_opt = filter_edges(graph, filter_value)
            global EC_val = get_EC(fg_opt)
            global avg_deg = average_degree(fg_opt)
            max_cl = maximal_cliques(fg_opt.g)
            global len_cliques = length(max_cl)
            global k_clique    = length(max_cl[1])
            global k_val       = length(clique_percolation(fg_opt.g, k = 6))
            break
        end
    end
    return all_comps, f_val, EC_val, avg_deg, len_cliques, k_clique, k_val
end

sols = zeros((2900, 9))

for i in 1:2900
    cor_mat = cor(data[i:(i + 6), :], dims = 1)

    bit_vec = (!).(isnan.(cor_mat[:, 1]))

    sym_mat = cor_mat[bit_vec, bit_vec]

    node_count = size(sym_mat, 1)

    if node_count <= 2
        continue
    end

    sym_graph = symmetric_matrix_to_graph(sym_mat)

    all_comps, f_val, EC_val, avg_deg, lc, mc, kval = find_smallest_filter_graph(sym_graph, .01)

    sols[i, 1] = node_count
    sols[i, 2] = all_comps
    sols[i, 3] = f_val
    sols[i, 4] = EC_val
    sols[i, 5] = avg_deg
    sols[i, 6] = EC_val / all_comps
    sols[i, 7] = lc
    sols[i, 8] = mc
    sols[i, 9] = kval


    if i%10 == 0
        println("DONE WITH $i")
    end
end


new_cases = zeros(2900)

for i in 1:2900
    data_to_avg = data[i:(i + 6), :]

    new_cases[i] = sum(data_to_avg) / 7

end

plot(1:2900, new_cases, legend = false)
xlabel!("Time (day)")
ylabel!("Average New Cases (7 day window)")

plot(1:2900, sols[:, 4], legend = false)
xlabel!("Time (day)")
ylabel!("Euler Characteristic")

plot(1:2900, sols[:, 9], legend = false)
xlabel!("Time (day)")
ylabel!("Number of Communities from Clique Filtration")


#f_val = 0
#for i in 1:-.01:-.25
#    fg = filter_edges(cor1g, i)
#    if length(connected_components(fg)) > 1
#        global f_val = i + .01
#        println(get_EC(filter_edges(cor1g, f_val)))
#        break
#    end
#end
