"""
contour_plot.jl
objectives : for plotting contour plot with different
            mesh geometry and precision
"""

using Plots

function plot_mesh(nodes)
    points = vec(nodes)
    x = [p[1] for p in points]
    y = [p[2] for p in points]
    
    # Get min and max values
    x_min, x_max = minimum(x), maximum(x)
    y_min, y_max = minimum(y), maximum(y)
    
    s = scatter(x, y, 
            aspect_ratio=:equal,
            markersize=1,
            xlabel="X",
            ylabel="Y",
            legend=false,
            title="Computational Nodes",
            xticks=[x_min, x_max],
            yticks=[y_min, y_max],
            fontfamily="Computer Modern")
    display(s)
end

"""
TODO  

"""
function plot_temp(data)

    # Get unique coordinates
    x_vals = sort(unique([p[1] for p in data]))
    y_vals = sort(unique([p[2] for p in data]))
    
    # Create temperature matrix
    temp_matrix = zeros(length(y_vals), length(x_vals))
    for point in data
        i = findfirst(y -> y ≈ point[2], y_vals)
        j = findfirst(x -> x ≈ point[1], x_vals)
        temp_matrix[i, j] = point[3]
    end
    
    # contour plot
    c = contourf(x_vals, y_vals, temp_matrix,
            title="Temperature Distribution",
            xlabel="X Coordinate (m)", 
            ylabel="Y Coordinate (m)",
            colorbar_title="Temperature (°C)",
                 clabels = true,
                 linecolor=:white,
            aspect_ratio=:equal,
            levels=25,
            color=:turbo,
            size=(700, 600),
            dpi=300,
            linewidth=0,
            fontfamily="Computer Modern")
    display(c)
end

function plot_velocity(nodes)
    
    # Extract x, y, u, v from nodes
    n_i, n_j = size(nodes)
    xp = [nodes[i,j][1] for i in 1:n_i, j in 1:n_j]
    yp = [nodes[i,j][2] for i in 1:n_i, j in 1:n_j]
    u = [nodes[i,j][4] for i in 1:n_i, j in 1:n_j]
    v = [nodes[i,j][5] for i in 1:n_i, j in 1:n_j]
    
    # Plot velocity field
    vec = 5
    q = quiver(xp, yp, quiver=(u', v', vec), 
           aspect_ratio=:equal, 
           xlabel="x", ylabel="y",
           legend=false)
    display(q)
end

# Plot with log scaling on x-axis (error tolerance)
function plot_error(error_trial, counter_array)
    error_vec = vec(error_trial)
    counter_vec = vec(counter_array)
    plt = plot(
        error_vec, counter_vec;
        xscale = :log10,
        xlabel = "Error tolerance",
        ylabel = "Number of iterations",
        title = "Convergence",
        legend = false,
        lw = 2,
        marker = :circle,
        markersize = 6,
        markercolor = :blue,
        linecolor = :blue,
        grid = :on,
        framestyle = :box,
        background_color = :white,
        guidefont = font(12, "Computer Modern"),
        tickfont = font(10, "Computer Modern"),
        titlefont = font(14, "Computer Modern")
    )

    display(plt)
    return plt
end
