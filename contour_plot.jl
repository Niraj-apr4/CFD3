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

function contour_plot(data)
    """
    TODO  

    """

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
            aspect_ratio=:equal,
            levels=25,
            color=:hot,
            size=(700, 600),
            dpi=300,
            linewidth=0,
            fontfamily="Computer Modern")
    display(c)
end

function plot_convergence(loops_array, tolerance_array)
    plot(loops_array, tolerance_array,
         xlabel="Iteration Number",
         ylabel="Convergence Tolerance",
         title="Convergence History(logarithmic plot)",
         yscale=:log10,                    # Log scale for tolerance
         color=:blue,
         marker=:square,
         grid=:true,
         framestyle=:box,
         markersize=5,
         size=(800, 600),
         fontfamily="Computer Modern",
         titlefontsize=16,
         guidefontsize=14,
         tickfontsize=12,
         legend=false)
    
    # annotation
    hline!([tolerance_array[7]], linestyle=:dash, color=:green, label="Target")
    vline!([loops_array[7]], linestyle=:dash, color=:green, label="Target")
end
