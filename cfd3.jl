using DelimitedFiles
include("plot_functions.jl")
include("mesh_geometry.jl")

# Parameters
rho = 1
k = 1
gamma = 1/50
H = 2


# STEP 1 create mesh ######################### 

# import mesh 
nodes  = nodes_given
points = points_given

n = size(points)[1]
coeff_array = Array{Array{Float64}}(undef,n-1,n-1)
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

# identify the boundary 
BoundaryS = nodes[:   ,   1] 
BoundaryE = nodes[n-1 ,   :]
BoundaryN = nodes[:   , n-1]
BoundaryW = nodes[1   ,   :]

# identify inlets and  outlets nodes 
IOlength = 0.068 * H 
IOnodes = 5  # approx


# confirm the plot 
# contour_plot(nodes)

# BoundaryW 
for i = 1:length(BoundaryW)
    #fix BoundaryW temp value to 50
    if i <= IOnodes # outlet B
        BoundaryW[i][4] = 0 # UB = 0 

    elseif i >= size(BoundaryW)[1]-IOnodes # outlet A 
        BoundaryW[i][4] = 1  #  UA = 1
        BoundaryW[i][3] = 20  # TA = 20 (Dirichlet)
    end
end

# boundaryE
for i = 1:length(BoundaryE)
    if i <= IOnodes #  outlets C
        BoundaryE[i][4] = 1  # UC = 1
    else # other than outlet
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
        BoundaryN[i][5] = 0 #VD = 0 
    end
    # BoundaryN[i][3] = 50
end

#BoundaryS
# for i = 1:length(BoundaryS)
#     BoundaryS[i][3] = 50
# end

# contour_plot(nodes)
############################################################

# {{{
# coeff_array size same as node array 
function coeff_array!(nodes)
    # global parameters
    global coeff_array, rho, IOnodes, n, gamma
    
    for i=1:n-1 # iteration along rows
        for j = 1:n-1 # iteration along columns
            
            # Initialize
            aW = 0.0
            aE = 0.0
            aS = 0.0
            aN = 0.0
            Su = 0.0
            Sp = 0.0
            
            # Check if we're at a boundary that would cause out of bounds
            if i == 1 || i == n-1 || j == 1 || j == n-1
                # For boundary nodes, set default coefficients
                # These will be overridden by boundary conditions below
                P_node = nodes[i,j]
                
                # Get adjacent nodes safely
                W_node = (i > 1) ? nodes[i-1, j] : P_node
                E_node = (i < n-1) ? nodes[i+1, j] : P_node
                S_node = (j > 1) ? nodes[i, j-1] : P_node
                N_node = (j < n-1) ? nodes[i, j+1] : P_node
                
                # Calculate geometry (use safe values)
                if i > 1
                    delta_xWP = P_node[1] - W_node[1]
                else
                    delta_xWP = 0.0
                end
                
                if i < n-1
                    delta_xPE = E_node[1] - P_node[1]
                else
                    delta_xPE = 0.0
                end
                
                if j > 1
                    delta_ySP = P_node[2] - S_node[2]
                else
                    delta_ySP = 0.0
                end
                
                if j < n-1
                    delta_yPN = N_node[2] - P_node[2]
                else
                    delta_yPN = 0.0
                end
                
                delta_x = (delta_xWP + delta_xPE) / 2
                delta_y = (delta_ySP + delta_yPN) / 2
                
                Aw = max(delta_y, 1e-10)
                Ae = max(delta_y, 1e-10)
                An = max(delta_x, 1e-10)
                As = max(delta_x, 1e-10)
                
                # Mass flow rates
                Fw = (i > 1) ? rho * (W_node[4] + P_node[4]) / 2 * Aw : 0.0
                Fe = (i < n-1) ? rho * (E_node[4] + P_node[4]) / 2 * Ae : 0.0
                Fs = (j > 1) ? rho * (S_node[5] + P_node[5]) / 2 * As : 0.0
                Fn = (j < n-1) ? rho * (N_node[5] + P_node[5]) / 2 * An : 0.0
                
                # Diffusion terms
                Dw = (i > 1 && delta_xWP > 0) ? gamma / delta_xWP * Aw : 0.0
                De = (i < n-1 && delta_xPE > 0) ? gamma / delta_xPE * Ae : 0.0
                Ds = (j > 1 && delta_ySP > 0) ? gamma / delta_ySP * As : 0.0
                Dn = (j < n-1 && delta_yPN > 0) ? gamma / delta_yPN * An : 0.0
                
                # Apply boundary conditions based on location
                # (boundary condition logic will follow below)
                
            else
                # Internal node - normal calculation
                P_node = nodes[i,j]
                W_node = nodes[i-1, j]
                S_node = nodes[i, j-1]
                E_node = nodes[i+1, j]
                N_node = nodes[i, j+1]
                
                # node geometry
                delta_xWP = P_node[1] - W_node[1]
                delta_xPE = E_node[1] - P_node[1]  
                delta_ySP = P_node[2] - S_node[2]
                delta_yPN = N_node[2] - P_node[2]
                delta_x = delta_xWP/2 + delta_xPE/2
                delta_y = delta_ySP/2 + delta_yPN/2
                
                Aw = delta_y
                Ae = delta_y
                An = delta_x
                As = delta_x
                
                # Mass flow rates
                Fw = rho * (W_node[4] + P_node[4]) / 2 * Aw
                Fe = rho * (E_node[4] + P_node[4]) / 2 * Ae
                Fs = rho * (S_node[5] + P_node[5]) / 2 * As
                Fn = rho * (N_node[5] + P_node[5]) / 2 * An
                
                delF = Fe - Fw + Fn - Fs
                
                Dw = gamma / delta_xWP * Aw
                De = gamma / delta_xPE * Ae
                Ds = gamma / delta_ySP * As
                Dn = gamma / delta_yPN * An
            end
            
            # Now apply boundary conditions
            # DIRICHLET CONDITION
            # Boundary node at A (left side, inlet with u=1)
            if i == 1 && nodes[i,j][4] == 1
                aE = (i < n-1) ? max(-Fe, De - Fe/2, 0) : 0
                aW = 0
                aN = (j < n-1) ? max(-Fn, Dn - Fn/2, 0) : 0
                aS = (j > 1) ? max(Fs, Ds + Fs/2, 0) : 0
                Sp = -(Dw + Fw)
                Su = (Dw + Fw) * nodes[i,j][3]  # Use its own temp for BC
            
            # Boundary node at right except at C (right side with T=50)
            elseif i == n-1 && j <= n-1-IOnodes
                aE = 0
                aW = (i > 1) ? max(Fw, Dw + Fw/2, 0) : 0
                aN = (j < n-1) ? max(-Fn, Dn - Fn/2, 0) : 0
                aS = (j > 1) ? max(Fs, Ds + Fs/2, 0) : 0
                Sp = -(De - Fe)
                Su = (De - Fe) * nodes[i,j][3]  # T=50 is already set
            
            # NEUMANN CONDITION (insulated)
            # Left side boundaries (except A)
            elseif i == 1
                aE = (i < n-1) ? max(-Fe, De - Fe/2, 0) : 0
                aW = 0
                aN = (j < n-1) ? max(-Fn, Dn - Fn/2, 0) : 0
                aS = (j > 1) ? max(Fs, Ds + Fs/2, 0) : 0
                Sp = 0
                Su = 0
            
            # Boundary node at C (right side outlet)
            elseif i == n-1
                aE = 0
                aW = (i > 1) ? max(Fw, Dw + Fw/2, 0) : 0
                aN = (j < n-1) ? max(-Fn, Dn - Fn/2, 0) : 0
                aS = (j > 1) ? max(Fs, Ds + Fs/2, 0) : 0
                Sp = 0
                Su = 0
            
            # Top Node
            elseif j == n-1
                aE = (i < n-1) ? max(-Fe, De - Fe/2, 0) : 0
                aW = (i > 1) ? max(Fw, Dw + Fw/2, 0) : 0
                aN = 0
                aS = (j > 1) ? max(Fs, Ds + Fs/2, 0) : 0
                Sp = 0
                Su = 0
            
            # Bottom node
            elseif j == 1
                aE = (i < n-1) ? max(-Fe, De - Fe/2, 0) : 0
                aW = (i > 1) ? max(Fw, Dw + Fw/2, 0) : 0
                aN = (j < n-1) ? max(-Fn, Dn - Fn/2, 0) : 0
                aS = 0
                Sp = 0
                Su = 0
            
            # Interior nodes 
            else
                aW = max(Fw, Dw + Fw/2, 0) 
                aE = max(-Fe, De - Fe/2, 0) 
                aS = max(Fs, Ds + Fs/2, 0) 
                aN = max(-Fn, Dn - Fn/2, 0) 
                Su = 0
                Sp = 0
            end
            
            aP = aW + aE + aS + aN - Sp
            coeff_array[i,j] = [aW, aE, aS, aN, aP, Su, Sp]
        end
    end
    return coeff_array
end
coeff_array = coeff_array!(nodes)

# TODO review needed 
function TDMA!(nodes)
    # global variables
    global coeff_array, n
    
    # Solve along j-direction (columns) for each row i
    for i = 1:n-1
        # Skip if this is a boundary row where we shouldn't solve
        if i == 1 || i == n-1
            continue
        end
        
        m = n - 1  # size of the system (all j nodes in coeff_array)
        
        # Arrays for tridiagonal system
        a = zeros(m)  # lower diagonal
        b = zeros(m)  # main diagonal
        c = zeros(m)  # upper diagonal
        d = zeros(m)  # RHS
        
        # Build the tridiagonal system
        for j = 1:n-1
            aW = coeff_array[i,j][1] 
            aE = coeff_array[i,j][2] 
            aS = coeff_array[i,j][3] 
            aN = coeff_array[i,j][4]
            aP = coeff_array[i,j][5]
            Su = coeff_array[i,j][6]
            
            # Safely get adjacent temperatures
            TW = (i > 1) ? nodes[i-1, j][3] : 0.0
            TE = (i < n-1) ? nodes[i+1, j][3] : 0.0
            
            a[j] = -aS
            b[j] = aP
            c[j] = -aN
            d[j] = aW * TW + aE * TE + Su
        end
        
        # Forward elimination
        for j = 2:m
            if b[j-1] != 0
                factor = a[j] / b[j-1]
                b[j]=b[j] - factor * c[j-1]
                d[j] = d[j] - factor * d[j-1]
            end
        end
        
        # Back substitution
        T = zeros(m)
        if b[m] != 0
            T[m] = d[m] / b[m]
        end
        for j = m-1:-1:1
            if b[j] != 0
                T[j]=(d[j] - c[j] * T[j+1]) / b[j]
            end
        end
        
        # Update nodes
        for j = 1:n-1
            nodes[i, j][3] = T[j]
        end
    end
end

# TODO Gauss_Siedel donot converge and solution dont agree
# check this one 

function Gauss_Siedel!(nodes)
    # global variables
    global coeff_array, n

    # Iterate over all nodes in coeff_array
    for i = 1:n-1, j = 1:n-1
        aW = coeff_array[i,j][1] 
        aE = coeff_array[i,j][2] 
        aS = coeff_array[i,j][3] 
        aN = coeff_array[i,j][4]
        aP = coeff_array[i,j][5]
        Su = coeff_array[i,j][6]
        Sp = coeff_array[i,j][7]

        # Safely get adjacent temperatures
        TW = (i > 1) ? nodes[i-1, j][3] : 0.0
        TS = (j > 1) ? nodes[i, j-1][3] : 0.0
        TE = (i < n-1) ? nodes[i+1, j][3] : 0.0
        TN = (j < n-1) ? nodes[i, j+1][3] : 0.0
        
        # Calculate new temperature
        if aP != 0
            T_new = (aW*TW + aE*TE + aS*TS + aN*TN + Su) / aP
            # update the temp 
            nodes[i,j][3] = T_new
        end
    end
end

function error_test(nodes)
    global coeff_array
    # only implemented for internal nodes
    # for simplicity
    for i = 2:n-2 ,j = 2:n-2
        aW = coeff_array[i,j][1] 
        aE = coeff_array[i,j][2] 
        aS = coeff_array[i,j][3] 
        aN = coeff_array[i,j][4]
        aP = coeff_array[i,j][5]

        TW = nodes[i-1, j][3] 
        TS = nodes[i, j-1][3]
        TE = nodes[i+1, j][3]
        TN = nodes[i, j+1][3]

        error = aW * TW + aS * TS + aE * TE + aN * TN - aP * TP
        @show error
    end
end
