
using DelimitedFiles

"""
# Arguments
- `L::Float`: length of full cv
- `H::Float`: Height of full cv
- `n::Integer`: number of differential cv
                more n more  finer the mesh is 
# Outputs
[points,nodes]
- `points::Array{Array{Float64}}`:
- `nodes::Array{Array{Float64}}`:

"""
function generate_2Dmesh(L,H,n,sx=1,sy=1)
    points = Array{Array{Float64}}(undef,n,n)

    # create the CV grid >>> 
    if sx == 1
        #unifrom
        x_array = range(0, L, n)  
    else
        # Geometric stretching
        x_array = [L * (sx^(i-1) - 1) / (sx^(n-1) - 1) for i in 1:n]
    end
    
    if sy == 1
        #uniform
        y_array = range(0, H, n) 
    else
        # Geometric stretching
        y_array = [H * (sy^(j-1) - 1) / (sy^(n-1) - 1) for j in 1:n]
    end

    x_array = range(0,L,n)
    y_array = range(0,H,n)
    for i = 1:n 
        for j = 1:n 
            points[i,j] = [x_array[i]*sx, y_array[j]*sy ,0]
        end
    end
    # now the points array consists of CV grids with initial temp
    # assigned to be  0 

    # STEP 2 ##### 

    # Compute  computational nodes for each cv >>>
    nodes = Array{Array{Float64}}(undef,n-1,n-1)
    for i = 1:n-1 ,j = 1:n-1
        init_x = points[i,j][1]
        init_y = points[i,j][2]

        final_x = points[i+1,j][1]
        final_y = points[i,j+1][2]

        delta_x = final_x - init_x
        delta_y = final_y - init_y

        nodes[i,j] = [points[i,j][1] + delta_x/2 , points[i,j][2] +
            delta_y/2 , 0]
    end
    # <<<
    return[points,nodes]
end



"""
# Arguments
- `xc_file::String`: path to xc.dat file
- `yc_file::String`: path to yc.dat file

# Outputs
[points,nodes]
- `points::Array{Array{Float64}}`: grid line coordinates
- `nodes::Array{Array{Float64}}`: cell center coordinates
"""
function generate_2Dmesh_from_dat(xc_file="xc.dat", yc_file="yc.dat")
    # Get the directory where this script is located
    script_dir = @__DIR__
    
    # Read grid coordinates from files (using full path)
    x_array = vec(readdlm(joinpath(script_dir, xc_file)))
    y_array = vec(readdlm(joinpath(script_dir, yc_file)))
    
    n = length(x_array)  # number of grid lines
    
    # Create points array (grid line intersections)
    points = Array{Array{Float64}}(undef, n, n)
    for i = 1:n 
        for j = 1:n 
            points[i,j] = [x_array[i], y_array[j], 0]
        end
    end
    
    # Compute computational nodes (cell centers)
    nodes = Array{Array{Float64}}(undef, n-1, n-1)
    for i = 1:n-1, j = 1:n-1
        init_x = points[i,j][1]
        init_y = points[i,j][2]
        final_x = points[i+1,j][1]
        final_y = points[i,j+1][2]
        delta_x = final_x - init_x
        delta_y = final_y - init_y
        nodes[i,j] = [points[i,j][1] + delta_x/2, points[i,j][2] + delta_y/2, 0]
    end
    
    return [points, nodes]
end

points_given, nodes_given = generate_2Dmesh_from_dat("xc.dat", "yc.dat")
