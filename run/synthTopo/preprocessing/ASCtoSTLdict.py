path = '../'
DEM_file = 'synthDEM.asc'

subsample = 2  # subsampling factor of the original DEM
dist0 = 25.0  # semi-width of the fissure/radius of the cylinder
dist_flat = 15.0
enne = 1.0  # shape parameter (1 linear, 2 spherical)
depth = 10.0  # depth of the fissure
nlevels = 4  # levels of refinements of the subsampled grid
RBF_interpolation = True # smooth the topography inside the crater
top_smooth_flag = True # smooth the top of crater volume

lagrangian_layer_depth = 10.0

xb = 0.0  # horizontal x-translation of the base of the fissure/conduit
yb = 0.0  # horizontal y-translation of the base of the fissure/conduit

xc = -140.0  # x of first point of fissure/x-center of cylinder (UTM)
yc = 30.0  # y of first point of fissure/y-center of cylinder (UTM)

# FOR CYLINDRICAL FISSURE
points = [(xc, yc), (xc, yc)]

# FOR LINEAR FISSURE
# points = [(xc, yc), (xc+20, yc+10), (xc+40, yc+15)]

conduit_radius = 5.0
conduit_length = 20.0
conduit_buffer = 3.0

saveDicts_flag = True
z_atm = 2000
offset_mesh = 50.0
delta_mesh = 100.0
