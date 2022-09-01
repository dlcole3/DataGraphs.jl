
function plot_graph(g::DataGraph;
    get_new_positions::Bool=true,
    save_pos::Bool=true,
    plot_edges::Bool=true,
    C::Real=1.0,
    K::Real=.1,
    iterations::Int=300,
    tol::Real=.1,
    xdim = 800,
    ydim = 800,
    linewidth = 1,
    linealpha = 1,
    legend::Bool = false,
    color = :blue,
    markercolor = :black,
    markersize = 5
)

    plt_options = Dict(:framestyle => :box, :grid => false, :size => (xdim,ydim), :axis => nothing, :legend => legend)
    line_options = Dict(:linecolor => color, :linewidth => linewidth, :linealpha => linealpha)

    am = create_adj_mat(g)

    if get_new_positions
        pos = NetworkLayout.sfdp(LightGraphs.Graph(am); tol = tol, C = C, K = K, iterations = iterations)
    else
        pos = g.node_positions
    end

    if save_pos
        g.node_positions = pos
    end

    if !get_new_positions && length(g.node_positions) == 0
        error("No positions in graph")
    end


    plt = scatter([i[1] for i in pos], [i[2] for i in pos];markercolor=markercolor, markersize=markersize, plt_options...)

    node_dict = g.nodes_index

    if plot_edges
        for i in g.edges
            from = i[1]
            to   = i[2]

            from_index = node_dict[from]
            to_index   = node_dict[to]

            plot!(plt,[pos[from_index][1], pos[to_index][1]], [pos[from_index][2], pos[to_index][2]]; line_options...)
        end
    end
    display(plt)
end

function set_matrix_node_positions(nodes, mat)
    dim1, dim2 = size(mat)

    positions = []
    for i in 1:length(nodes)
        node_val  = nodes[i]
        node_x    = Float64(node_val[2] * 5)
        node_y    = Float64((dim1 - node_val[1] * 5))
        push!(positions, Point(node_x, node_y))
    end

    return positions
end
