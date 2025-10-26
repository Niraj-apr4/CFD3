
include("contour_plot.jl")
include("mesh_geometry.jl")

######################
# STEP 1 create mesh # 
######################

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
