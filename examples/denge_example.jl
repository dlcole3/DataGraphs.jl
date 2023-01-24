using Revise
using DataGraphs, Graphs, DelimitedFiles, Statistics
using DataGraphPlots

# The data and methods in this example are primarily from the paper of
# Barros de Souze at al., 2022. https://doi.org/10.1088/1742-5468/aca0e5
# "The Euler characteristic as a topological marker for outbreaks in
# vector borne disease"

full_data = readdlm((@__DIR__)*"/Recife_data/Dengue_Recife_new_cases_jan_2014_to_dec_2021.csv", ',', skipstart = 1)[:, 2:end]

# Set the matrix type to be Int
data = Matrix{Int}(undef, size(full_data))
data .= full_data

# Define a function to get the smallest filtered graph that only contains one connected component
function find_smallest_filter_graph(graph, filter_value, k = 5)
    # Define total number of components
    num_comps = length(graph.nodes) + length(graph.edges)

    # Define variables as local
    local f_val, EC_val, avg_deg, max_cl, len_max_cl, k_clique, num_k_comms

    # iterate through filtration levels
    for i in maximum(get_edge_data(graph)[:, 1]):(-filter_value):minimum(get_edge_data(graph)[:, 1])
        filtered_graph = filter_edges(graph, i)
        # If the graph splits into more than 1 connected component at filtration level i,
        # then get the last graph filtration and calculate topological characteristics
        if length(connected_components(filtered_graph)) > 1
            f_val = i + filter_value # Filtration threshold
            filtered_graph_opt = filter_edges(graph, filter_value)
            EC_val = get_EC(filtered_graph_opt) # Euler characteristic
            avg_deg = average_degree(filtered_graph_opt) # Average node degree
            max_cl = maximal_cliques(filtered_graph_opt.g)
            len_max_cl  = length(max_cl) # Number of maximal cliques
            k_clique    = length(max_cl[1]) # Length of maximal cliques (k value)
            num_k_comms = length(clique_percolation(
                filtered_graph_opt.g, k = k)
                ) # Number of communities for k-clique percolation
            break
        end
    end

    return num_comps, f_val, EC_val, avg_deg, len_max_cl, k_clique, num_k_comms
end

# Create array for storing solutions
sols = zeros((2915, 9))

# Run through all 7 day windows in the data
for i in 1:2915
    # Form a correlation matrix based on 7 days of data
    cor_mat = cor(data[i:(i + 6), :], dims = 1)

    # Remove the matrix entries that are NaNs
    bit_vec = (!).(isnan.(cor_mat[:, 1]))
    sym_mat = cor_mat[bit_vec, bit_vec]

    # If there are not more than 2 nodes in the graph, skip this iteration
    node_count = size(sym_mat, 1)

    if node_count <= 2
        continue
    end

    # Build the edge-weighted graph from the correlation matrix
    sym_graph = symmetric_matrix_to_graph(sym_mat)

    # Find the smallest filtration level possible and get the TDA metrics
    num_comps, f_val, EC_val, avg_deg, len_max_cl, k_clique, num_k_comms = find_smallest_filter_graph(sym_graph, .01, 5)

    sols[i, 1] = node_count
    sols[i, 2] = num_comps
    sols[i, 3] = f_val
    sols[i, 4] = EC_val
    sols[i, 5] = avg_deg
    sols[i, 6] = EC_val / num_comps
    sols[i, 7] = len_max_cl
    sols[i, 8] = k_clique
    sols[i, 9] = num_k_comms


    if i%10 == 0
        println("DONE WITH $i")
    end
end

# Get the moving average for a seven day window
new_cases = zeros(2915)

for i in 1:2915
    data_to_avg = data[i:(i + 6), :]

    new_cases[i] = sum(data_to_avg) / 7
end


# Plot results
using Plots
plot(1:2915, new_cases, legend = false)
xlabel!("Time (day)")
ylabel!("Average New Cases (7 day window)")

plot(1:2915, sols[:, 4], legend = false)
xlabel!("Time (day)")
ylabel!("Euler Characteristic")

plot(1:2915, sols[:, 9], legend = false)
xlabel!("Time (day)")
ylabel!("Number of Communities from 5-Clique Filtration")
