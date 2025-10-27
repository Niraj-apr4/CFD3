using DelimitedFiles
include("plot_functions.jl")
include("mesh_geometry.jl")

# Parameters
rho = 1
k = 1
gamma = 1/50
H = 2


# STEP 1 create mesh ######################### 

# earlier logic of creating generating mesh  >>
# n = 30 
# L = 1 
# H = 0.5
# sx = 1
# sy = 1 

# mesh2 = generate_2Dmesh(L,H,n,sx,sy)
# nodes = mesh2[2]
# points = mesh2[1]
# <<

nodes  = nodes_given
points = points_given
################################################

# STEP 2 add velocity field ###################
# >> add velocity to nodes >>
function add_velocity_to_nodes!(nodes, u_file="u.dat", v_file="v.dat")
    
    script_dir = @__DIR__
    
    # Read velocity data
    u = vec(readdlm(joinpath(script_dir, u_file)))
    v = vec(readdlm(joinpath(script_dir, v_file)))
    
    # Get actual dimensions from data length
    total_points = length(u)
    ni = Int(sqrt(total_points))  # Assuming square grid
    nj = ni
    
    # Reshape to 2D
    u2d = reshape(u, ni, nj)
    v2d = reshape(v, ni, nj)
    
    # Get dimensions from nodes array
    n_i = size(nodes, 1)
    n_j = size(nodes, 2)
    
    # Assign u and v to nodes (taking first n_i×n_j values)
    for i = 1:n_i, j = 1:n_j
        nodes[i,j][4] = u2d[i,j]  # u velocity
        nodes[i,j][5] = v2d[i,j]  # v velocity
    end
    
    return nodes
end

add_velocity_to_nodes!(nodes, "u.dat", "v.dat")
# plot_velocity(nodes)

###########################################################

# STEP 3 Add boundary conditon ############################

# TODO check if there will be points or nodes
n = size(points)[1]

# identify the boundary 
BoundaryS = nodes[:   ,   1] 
BoundaryE = nodes[n-1 ,   :]
BoundaryN = nodes[:   , n-1]
BoundaryW = nodes[1   ,   :]

# identify inlets and  outlets nodes 
IOlength = 0.068 * H 
IOnodes = 11  # approx


# confirm the plot 
# contour_plot(nodes)

# BoundaryW 
for i = 1:length(BoundaryW)
    #fix BoundaryW temp value to 50
    if i <= IOnodes # outlet B
        BoundaryW[i][4] = 0 # conditions
        BoundaryW[i][3] = 20  # conditons 

    elseif i >= size(BoundaryW)[1]-IOnodes # outlet A 
        BoundaryW[i][4] = 1  # conditons 
        BoundaryW[i][3] = 20  # conditons 
    end
end

# boundaryE
for i = 1:length(BoundaryE)
    if i <= IOnodes #  outlets C
        BoundaryE[i][4] = 1  # conditons 
        BoundaryE[i][3] = 50 
    else
        BoundaryE[i][3] = 50 
    end
end

#boundaryN 
for i = 1:length(BoundaryN)
    mid = length(BoundaryN)/2
    mid = Int64(round(Int , mid))

    IOmid = IOnodes/2
    IOmid = Int64(round(Int , IOmid))
    # @show IOmid
    if i >= mid - IOmid && i <= mid + IOmid 
        # @show i
        BoundaryN[i][5] = 0 
    end
end

BoundaryS
for i = 1:length(BoundaryS)
    #fix BoundaryS temp value to 50
    BoundaryS[i][3] = 50
end

# contour_plot(nodes)
############################################################

# {{{
# coeff_array size same as node array 
function gen_coeff_array(nodes)
    coeff_array = Array{Array{Float64}}(undef,n-1,n-1)
    for i=1:n-1 # iteration along rows
        for j = 1:n-1 # iteration along columns

            global rho # global parameters

            # TODO change the boundary conditons
            # parameters

            if i == 1 || j == 1 || i == n-1 || j == n-1 
                coeff_array[i,j] = [0 0 0 0 0]
            else
                P_node = nodes[i,j] # current node
                W_node = nodes[i-1 , j  ] # Adjacent nodes
                S_node = nodes[i   , j-1]
                E_node = nodes[i+1 , j  ]
                N_node = nodes[i   , j+1]

                # node geometry
                delta_xWP = P_node[1] - W_node[1]
                delta_xPE = E_node[1] - P_node[1]  
                delta_ySP = P_node[2] - S_node[2]
                delta_yPN = N_node[2] - P_node[2]

                # TODO Please confirm if this is correct 
                delta_x = delta_xWP/2 + delta_xPE/2
                delta_y = delta_ySP/2 + delta_yPN/2

                Aw = delta_y # face areas in 2D case
                Ae = delta_y
                An = delta_x
                As = delta_x

                # Fw = rho * W_node[4] * Aw
                # Fe = rho * E_node[4] * Ae
                # Fs = rho * S_node[5] * As
                # Fn = rho * N_node[5] * An

                #           OR

                # todo check which implementation is correct
                Fw = rho * (W_node[4] + P_node[4]) / 2 * Aw
                Fe = rho * (E_node[4] + P_node[4]) / 2 * Aw
                Fs = rho * (S_node[5] + P_node[5]) / 2 * Aw
                Fn = rho * (N_node[5] + P_node[5]) / 2 * Aw


                delF = Fe - Fw + Fn - Fs

                Dw = gamma / delta_xWP * Aw
                De = gamma / delta_xPE * Ae
                Ds = gamma / delta_ySP * As
                Dn = gamma / delta_yPN * An

                aW = max( Fw, (Dw + Fw/2),0) 
                aE = max(-Fe, (De - Fe/2),0) 
                aS = max( Fs, (Ds + Fs/2),0) 
                aN = max(-Fn, (Dn - Fn/2),0) 
                aP = aW + aE+ aS + aN + delF

                coeff_array[i,j] = [aW aE aS aN aP]
            end
        end
    end
    return coeff_array
end
# }}}

coeff_array = gen_coeff_array(nodes)

# TODO review needed 
function TDMA!(nodes)
    # global variables
    global coeff_array

    n = size(nodes)[1]
    
    # Solve along j-direction (columns) for each row i
    for i = 2:n-2
        m = n - 3  # size of the system
        
        # Arrays for tridiagonal system
        a = zeros(m)  # lower diagonal
        b = zeros(m)  # main diagonal
        c = zeros(m)  # upper diagonal
        d = zeros(m)  # RHS
        
        # Build the tridiagonal system
        for j = 2:n-2
            loc = j - 1
            
            aW = coeff_array[i,j][1] 
            aE = coeff_array[i,j][2] 
            aS = coeff_array[i,j][3] 
            aN = coeff_array[i,j][4]
            aP = coeff_array[i,j][5]
            
            TW = nodes[i-1, j][3]
            TE = nodes[i+1, j][3]
            
            a[loc] = -aS
            b[loc] = aP
            c[loc] = -aN
            d[loc] = aW * TW + aE * TE
        end
        
        # Forward elimination
        for j = 2:m
            factor = a[j] / b[j-1]
            b[j] = b[j] - factor * c[j-1]
            d[j] = d[j] - factor * d[j-1]
        end
        
        # Back substitution
        T = zeros(m)
        T[m] = d[m] / b[m]
        for j = m-1:-1:1
            T[j] = (d[j] - c[j] * T[j+1]) / b[j]
        end
        
        # Update nodes
        for j = 2:n-2
            nodes[i, j][3] = T[j-1]
        end
    end
end


function Gauss_Siedel!(nodes)
    # global variables
    global coeff_array

    # TODO implement Gauss siedel 
    for i = 2 : n-2, j = 2 : n-2 # internal nodes
        aW = coeff_array[i,j][1] 
        aE = coeff_array[i,j][2] 
        aS = coeff_array[i,j][3] 
        aN = coeff_array[i,j][4]
        aP = coeff_array[i,j][5]

        TW = nodes[i-1, j][3]
        TS = nodes[i, j-1][3]
        TE = nodes[i+1, j][3]
        TN = nodes[i, j+1][3]
        # Calculate new temperature
        T_new = (aW*TW + aE*TE + aS*TS + aN*TN) / aP

        # update the temp 
        nodes[i,j][3] = T_new
    end
end

# tolerance = 0.1
# coeff_array = gen_coeff_array(nodes)

# for i=1:10
#     TDMA!(nodes, coeff_array)
# end




# for i = 2:n-2 , j = 2:n-2
#     # global parameters
#     global tolerance 

#     aW = coeff_array[i,j][1]
#     aE = coeff_array[i,j][2]
#     aS = coeff_array[i,j][3]
#     aN = coeff_array[i,j][4]
#     aP = coeff_array[i,j][5]
    
#     TP = nodes[i,j][3]
#     TW = nodes[i-1 , j  ][3]
#     TS = nodes[i   , j-1][3]
#     TE = nodes[i+1 , j  ][3]
#     TN = nodes[i   , j+1][3]
#     # TODO please change 

#     global epsilon = aE* TE + aW * TW + aN * TN + aS * TS - aP * TP
#     @show epsilon 
#     while epsilon > tolerance
#         TDMA!(nodes, coeff_array)
#     end
# end

# TDMA!(nodes, coeff_array)
# contour_plot(nodes)
