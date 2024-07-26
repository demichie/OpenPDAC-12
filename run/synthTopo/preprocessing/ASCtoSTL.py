from shapely.geometry import LineString, Point
from shapely.affinity import translate
import numpy as np
from scipy.spatial import Delaunay
from stl import mesh
from scipy import interpolate
from linecache import getline
import sys
import os.path
import time
import pandas as pd

import skimage.filters.rank
import skimage.morphology
import skimage.io


def saveDicts(xmin, xmax, ymin, ymax, Zinit, offset_mesh, path):

    try:

        from ASCtoSTLdict import z_atm

    except ImportError:

        print('Missing parameter in dict: z_atm (float)')
        sys.exit(1)

    try:

        from ASCtoSTLdict import delta_mesh

    except ImportError:

        print('Missing parameter in dict: delta_mesh (float>0)')
        sys.exit(1)

    # Define the 3D block and compute the number of cells in each direction
    nx = int(np.ceil((xmax - xmin) / delta_mesh))
    ny = int(np.ceil((ymax - ymin) / delta_mesh))
    zmin = np.min(np.min(Zinit)) - offset_mesh
    zmax = np.max(np.max(Zinit)) + z_atm
    nz = int(np.ceil((zmax - zmin) / delta_mesh))

    # In order to write blockMeshDict dictionary, write the vertices of the
    # block and the discretization
    fid1 = open('blockMesh.mod', 'w')
    fid1.write('vertices\n')
    fid1.write('(\n')
    fid1.write('(%f %f %f) \n' % (xmin, ymin, zmin))
    fid1.write('(%f %f %f) \n' % (xmax, ymin, zmin))
    fid1.write('(%f %f %f) \n' % (xmax, ymax, zmin))
    fid1.write('(%f %f %f) \n' % (xmin, ymax, zmin))
    fid1.write('(%f %f %f) \n' % (xmin, ymin, zmax))
    fid1.write('(%f %f %f) \n' % (xmax, ymin, zmax))
    fid1.write('(%f %f %f) \n' % (xmax, ymax, zmax))
    fid1.write('(%f %f %f) \n' % (xmin, ymax, zmax))
    fid1.write(');\n\n')
    fid1.write('blocks\n')
    fid1.write('(\n')
    fid1.write(
        '    hex (0 1 2 3 4 5 6 7) (%d %d %d) simpleGrading (1 1 1) \n' %
        (nx, ny, nz))
    fid1.write(');\n\n')
    fid1.close()

    # Write blockMeshDict dictionary by concatenating the header with the other
    # parts. Save it as path/system/blockMeshDict
    path_system = path + 'system/'
    command = "cat templates/blockMesh.start blockMesh.mod templates/blockMesh.end > " + \
        path_system + "blockMeshDict"

    # Call the operating system to execute the specified commands
    os.system(command)
    os.system('rm blockMesh.mod')

    fout = open((path_system + 'snappyHexMeshDict'), 'w')
    with open('./templates/snappyHexMeshDict.template') as f:
        content = f.readlines()
        textdata = [x.strip() for x in content]
    for i in range(len(textdata)):
        textdata[i] = textdata[i].replace('xxxxxx', str(zmax - 10.0))
        fout.write('%s\n' % textdata[i])
    fout.close()


# Print iterations progress


def printProgressBar(iteration,
                     total,
                     prefix='',
                     suffix='',
                     decimals=1,
                     bar_length=100):
    """
    Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = 'X' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' %
                     (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def main(argv):

    try:

        from ASCtoSTLdict import DEM_file

    except ImportError:

        print('Missing parameter in dict: DEM_file (str)')
        sys.exit(1)

    print('Reading DEM file: ' + DEM_file)
    # Parse the topography header
    hdr = [getline(DEM_file, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    # values are associated to the centers of the DEM pixels
    xs_DEM = lx + 0.5 * cell + np.linspace(0, (cols - 1) * cell, cols)
    ys_DEM = ly + 0.5 * cell + np.linspace(0, (rows - 1) * cell, rows)

    # Load the topography into a numpy array
    DEM = pd.read_table(DEM_file,
                        delim_whitespace=True,
                        header=None,
                        skiprows=6).astype(float).values
    DEM = np.flipud(DEM)
    DEM[DEM == nd] = 0.0

    try:

        from ASCtoSTLdict import xc

    except ImportError:

        print('Missing parameter in dict: xc (float)')
        sys.exit(1)

    try:

        from ASCtoSTLdict import yc

    except ImportError:

        print('Missing parameter in dict: yc (float)')
        sys.exit(1)

    try:

        from ASCtoSTLdict import xb

    except ImportError:

        xb = 0.0

    try:

        from ASCtoSTLdict import yb

    except ImportError:

        yb = 0.0

    xinit = np.linspace(0, (cols - 1) * cell, cols) - xc
    yinit = np.linspace(0, (rows - 1) * cell, rows) - yc

    xinit = xs_DEM - xc
    yinit = ys_DEM - yc

    try:

        from ASCtoSTLdict import offset_mesh

    except ImportError:

        print('Missing parameter in dict: offset_mesh (float>0)')
        sys.exit(1)

    xmin = np.amin(xinit) + offset_mesh
    xmax = np.amax(xinit) - offset_mesh

    print('xmin,xmax', xmin, xmax)

    ymin = np.amin(yinit) + offset_mesh
    ymax = np.amax(yinit) - offset_mesh

    print('ymin,ymax', ymin, ymax)

    Xinit, Yinit = np.meshgrid(xinit, yinit)
    Zinit = DEM

    try:

        from ASCtoSTLdict import saveDicts_flag

    except ImportError:

        saveDicts_flag = False

    try:

        from ASCtoSTLdict import path

    except ImportError:

        print('Missing parameter in dict: path (str)')
        sys.exit(1)

    if saveDicts_flag:

        saveDicts(xmin, xmax, ymin, ymax, Zinit, offset_mesh, path)

    try:

        from ASCtoSTLdict import points

    except ImportError:

        print('Missing parameter in dict: points (coords)')
        sys.exit(1)

    # initialize the centerline
    line = LineString(points)

    # translate the linestring (relative reference system with (lx,ly)=(0,0))
    line = translate(line, -xc, -yc)

    bb = line.bounds

    print('bb', bb)

    # interpolate with values from original fine grid
    f = interpolate.interp2d(xinit, yinit, Zinit, kind='linear')

    try:

        from ASCtoSTLdict import subsample

    except ImportError:

        print('Missing parameter in dict: subsample (int)')
        sys.exit(1)

    # coarsening of the original grid
    xinit = xinit[::subsample]
    yinit = yinit[::subsample]

    Xinit = Xinit[::subsample, ::subsample]
    Yinit = Yinit[::subsample, ::subsample]
    Zinit = Zinit[::subsample, ::subsample]

    # list of points of DEM for the refined nested grid
    x_check = []
    y_check = []
    z_check = []
    z_check_org = []

    # list of displacements
    dx = []
    dy = []
    dz_total = []
    dz2 = []
    dz_crater = []

    try:

        from ASCtoSTLdict import dist0

    except ImportError:

        print('Missing parameter in dict: dist0 (float)')
        sys.exit(1)

    try:

        from ASCtoSTLdict import dist_flat

    except ImportError:

        print('Missing parameter in dict: dist_flat (float)')
        print('Set to default value: dist_flat=0')
        dist_flat = 0

    if dist_flat >= dist0:
    
        print('dist_flat should be < dist0')
        print('Set to default value: dist_flat=0')
        dist_flat = 0
            
    try:

        from ASCtoSTLdict import nlevels

    except ImportError:

        print('Missing parameter in dict: nlevels (int)')
        sys.exit(1)

    # boundaing box of first refinement level
    xmin_bb = bb[0] - nlevels * dist0
    ymin_bb = bb[1] - nlevels * dist0
    xmax_bb = bb[2] + nlevels * dist0
    ymax_bb = bb[3] + nlevels * dist0

    # x bounding box indexes
    idx_min = np.searchsorted(xinit, xmin_bb, side="left") - 1
    idx_max = np.searchsorted(xinit, xmax_bb, side="left")

    x0 = xinit[idx_min:idx_max + 1]

    # y bounding box indexes
    idy_min = np.searchsorted(yinit, ymin_bb, side="left") - 1
    idy_max = np.searchsorted(yinit, ymax_bb, side="left")

    y0 = yinit[idy_min:idy_max + 1]

    # add original points on the west of first bounding box
    x_check.extend(Xinit[:, 0:idx_min].ravel())
    y_check.extend(Yinit[:, 0:idx_min].ravel())
    z_check.extend(Zinit[:, 0:idx_min].ravel())
    z_check_org.extend(Zinit[:, 0:idx_min].ravel())

    dx.extend(0.0 * Xinit[:, 0:idx_min].ravel())
    dy.extend(0.0 * Xinit[:, 0:idx_min].ravel())
    dz_total.extend(0.0 * Xinit[:, 0:idx_min].ravel())
    dz2.extend(0.0 * Xinit[:, 0:idx_min].ravel())
    dz_crater.extend(0.0 * Xinit[:, 0:idx_min].ravel())

    # add original points on the south of first bounding box
    x_check.extend(Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
    y_check.extend(Yinit[0:idy_min, idx_min:idx_max + 1].ravel())
    z_check.extend(Zinit[0:idy_min, idx_min:idx_max + 1].ravel())
    z_check_org.extend(Zinit[0:idy_min, idx_min:idx_max + 1].ravel())

    dx.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
    dy.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
    dz_total.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
    dz2.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())
    dz_crater.extend(0.0 * Xinit[0:idy_min, idx_min:idx_max + 1].ravel())

    # add original points on the north of first bounding box
    x_check.extend(Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    y_check.extend(Yinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    z_check.extend(Zinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    z_check_org.extend(Zinit[idy_max + 1:, idx_min:idx_max + 1].ravel())

    dx.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    dy.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    dz_total.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    dz2.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())
    dz_crater.extend(0.0 * Xinit[idy_max + 1:, idx_min:idx_max + 1].ravel())

    # add original points on the east of first bounding box
    x_check.extend(Xinit[:, idx_max + 1:].ravel())
    y_check.extend(Yinit[:, idx_max + 1:].ravel())
    z_check.extend(Zinit[:, idx_max + 1:].ravel())
    z_check_org.extend(Zinit[:, idx_max + 1:].ravel())

    dx.extend(0.0 * Xinit[:, idx_max + 1:].ravel())
    dy.extend(0.0 * Xinit[:, idx_max + 1:].ravel())
    dz_total.extend(0.0 * Xinit[:, idx_max + 1:].ravel())
    dz2.extend(0.0 * Xinit[:, idx_max + 1:].ravel())
    dz_crater.extend(0.0 * Xinit[:, idx_max + 1:].ravel())

    dist_lev0 = 1.e10

    # loop over revinement levels
    for ilevel in range(nlevels):

        print('')
        print('Refinement level:', ilevel)

        X0, Y0 = np.meshgrid(x0, y0)

        X0_1d = X0.ravel()
        Y0_1d = Y0.ravel()

        nxy = len(X0_1d)

        dist_lev = dist0 * (nlevels - ilevel)
        print('Distance range:', dist_lev0, dist_lev)
        print('Points to check:', nxy)

        i = 0

        for j, (x, y) in enumerate(zip(X0_1d, Y0_1d)):

            printProgressBar(np.asarray(i) + 1,
                             nxy,
                             prefix='Progress:',
                             suffix='Complete',
                             decimals=1,
                             bar_length=50)

            dist = line.distance(Point(x, y))

            # distance of point (x,y) from line
            if (dist > dist_lev) and (dist <= dist_lev0):

                dx.append(0.0)
                dy.append(0.0)
                dz_total.append(0.0)
                dz2.append(0.0)
                dz_crater.append(0.0)

                x_check.append(x)
                y_check.append(y)

                z = f(x, y)

                z_check.append(float(z))
                z_check_org.append(float(z))

            i += 1

        # bounding box of new refinement level
        xmin_bb = bb[0] - dist_lev
        ymin_bb = bb[1] - dist_lev
        xmax_bb = bb[2] + dist_lev
        ymax_bb = bb[3] + dist_lev

        idx_min = np.searchsorted(x0, xmin_bb, side="left") - 1
        idx_max = np.searchsorted(x0, xmax_bb, side="left")
        nx = 2 * (idx_max - idx_min) + 1

        x1 = np.linspace(x0[idx_min], x0[idx_max], nx)

        idx_min = np.searchsorted(y0, ymin_bb, side="left") - 1
        idx_max = np.searchsorted(y0, ymax_bb, side="left")
        ny = 2 * (idx_max - idx_min) + 1
        y1 = np.linspace(y0[idx_min], y0[idx_max], ny)

        x0 = x1
        y0 = y1
        dist_lev0 = dist_lev

    print('')
    print('Refinement level:', nlevels)

    X1, Y1 = np.meshgrid(x1, y1)

    X1_1d = X1.ravel()
    Y1_1d = Y1.ravel()

    nxy = len(X1_1d)
    print('Distance range:', dist_lev0, 0)
    print('Points to check:', nxy)

    n_out = len(x_check)

    i = 0

    Inside_1d = np.zeros_like(X1_1d, dtype=np.uint8)

    # refinement for the level corresponding to the crater
    for j, (x, y) in enumerate(zip(X1_1d, Y1_1d)):

        # distance of point (x,y) from line
        dist = line.distance(Point(x, y))

        if (dist <= dist_lev):

            Inside_1d[j] = 1

    Inside = np.reshape(Inside_1d, X1.shape)

    # define "next to" - this may be a square, diamond, etc
    selem = skimage.morphology.disk(1)

    # create masks for the two kinds of edges
    Inside_edges = (skimage.filters.rank.minimum(Inside, selem) == 0) & (
        skimage.filters.rank.maximum(Inside, selem) == 255)

    shift = 1
    edgex1 = Inside & (np.roll(~Inside, shift=shift, axis=0)
                       | np.roll(~Inside, shift=-shift, axis=0))
    edgey1 = Inside & (np.roll(~Inside, shift=shift, axis=1)
                       | np.roll(~Inside, shift=-shift, axis=1))
    edge = edgex1 | edgey1

    print(np.sum(edge))

    xedge = X1[edge == True].reshape(-1, 1)
    yedge = Y1[edge == True].reshape(-1, 1)

    zedge = np.zeros_like(xedge)

    for j, (x, y) in enumerate(zip(xedge, yedge)):

        zedge[j] = f(x, y)

    print(zedge.shape)
    xy = np.concatenate([xedge, yedge], axis=1)

    from scipy.interpolate import RBFInterpolator
    rbf = RBFInterpolator(xy, zedge, epsilon=2)

    volume = 0.0

    try:

        from ASCtoSTLdict import conduit_radius

    except ImportError:

        print('Missing parameter in dict: conduit_radius (float)')
        sys.exit(1)

    if conduit_radius > 0:

        try:

            from ASCtoSTLdict import conduit_length

        except ImportError:

            print('Missing parameter in dict: conduit_length (float)')
            sys.exit(1)

        try:

            from ASCtoSTLdict import conduit_buffer

        except ImportError:

            print('Missing parameter in dict: conduit_buffer (float)')
            sys.exit(1)

        x_conduit_top = []
        y_conduit_top = []
        z_conduit_top = []
        x_conduit_bottom = []
        y_conduit_bottom = []
        z_conduit_bottom = []
        conduit_volume = 0.0

    # interpolation for the level corresponding to the crater
    for j, (x, y) in enumerate(zip(X1_1d, Y1_1d)):

        printProgressBar(np.asarray(i) + 1,
                         nxy,
                         prefix='Progress:',
                         suffix='Complete',
                         decimals=1,
                         bar_length=50)

        # distance of point (x,y) from line
        dist = line.distance(Point(x, y))

        try:

            from ASCtoSTLdict import enne

        except ImportError:

            print('Missing parameter in dict: enne (>0)')
            sys.exit(1)

        try:

            from ASCtoSTLdict import depth

        except ImportError:

            print('Missing parameter in dict: depth (>0)')
            sys.exit(1)

        if (dist < dist_lev):
        
            dz_rel = -(1.0 - (np.maximum(0,dist-dist_flat) / (dist0-dist_flat))**enne)**(1.0 / enne)

            if (conduit_radius > 0) and (dist <=
                                         conduit_radius + conduit_buffer):

                alpha_buffer = 1.0 - \
                    np.maximum(0.0, (dist-conduit_radius)/conduit_buffer)
                alpha_buffer = alpha_buffer**4

                conduit_dz_rel = alpha_buffer * conduit_length / depth
                dz_rel -= conduit_dz_rel
                conduit_volume += depth * conduit_dz_rel
                
            else:
            
                conduit_dz_rel = 0.0    

            dx.append(dz_rel * xb)
            dy.append(dz_rel * yb)
            dz_total.append(depth * dz_rel)
            dz_crater.append(depth * (dz_rel + conduit_dz_rel))

            volume += -depth * dz_rel

            try:

                from ASCtoSTLdict import lagrangian_layer_depth

            except ImportError:

                lagrangian_layer_depth = 0.0

            if lagrangian_layer_depth != 0:

                dz2.append(-lagrangian_layer_depth)

            else:

                dz2.append(depth * dz_rel)

            try:

                from ASCtoSTLdict import RBF_interpolation

            except ImportError:

                RBF_interpolation = False

            z_org = f(x, y)

            if RBF_interpolation:

                z = rbf([[x, y]])

            else:
            
                z = z_org


            if (conduit_radius > 0) and (dist <=
                                         conduit_radius + conduit_buffer):

                x_conduit_top.append(x + xb)
                y_conduit_top.append(y + yb)
                z_conduit_top.append(
                    float(z) + depth * (dz_rel + conduit_dz_rel))
                    
                x_conduit_bottom.append(x + dz_rel * xb)
                y_conduit_bottom.append(y + dz_rel * yb)
                z_conduit_bottom.append(float(z) + depth * dz_rel)

            x_check.append(x)
            y_check.append(y)
            z_check.append(float(z))
            z_check_org.append(float(z_org))

        i += 1

    volume *= (X1[0, 1] - X1[0, 0])**2

    print(f"Crater volume [m3] : {float(volume):e}")

    x_org = np.array(x_check)
    y_org = np.array(y_check)
    z_org = np.array(z_check_org)
    z_smooth = np.array(z_check)

    x_new = np.array(x_check) + np.array(dx)
    y_new = np.array(y_check) + np.array(dy)
    z_new = np.array(z_check) + np.array(dz_total)
    z_new2 = np.array(z_check) + np.array(dz2)
    z_new_crater = np.array(z_check) + np.array(dz_crater)



    points = []

    for i, (x, y) in enumerate(zip(x_check, y_check)):

        points.append([x, y])

    points = np.asarray(points)

    print('')
    print('Building Delaunay triangulation')

    tri = Delaunay(points)
    faces = tri.simplices

    print('Delaunay triangulation completed')

    # split the triangulation
    inner_tri_list = []

    tri_simpl = tri.simplices.copy()

    
    vert = range(n_out, len(x_check))

    # check if the points are in triangles (at least one vertex)
    inner_tri_list = np.in1d(tri_simpl[:, 0], vert)
    inner_tri_list = np.logical_or(inner_tri_list,
                                   np.in1d(tri_simpl[:, 1], vert))
    inner_tri_list = np.logical_or(inner_tri_list,
                                   np.in1d(tri_simpl[:, 2], vert))

    inner_tri_list = np.arange(tri_simpl.shape[0])[inner_tri_list]

    tri_lagr = tri_simpl[inner_tri_list, :]
       
    # list of indexes of points outside crater area 
    points_out = list(range(n_out))
    
    

    # check if the points are in triangles (at least one vertex)
    outer_tri_list = np.in1d(tri_simpl[:, 0], points_out)
    outer_tri_list = np.logical_or(outer_tri_list,
                                   np.in1d(tri_simpl[:, 1], points_out))
    outer_tri_list = np.logical_or(outer_tri_list,
                                   np.in1d(tri_simpl[:, 2], points_out))

    outer_tri_list = list(np.arange(tri_simpl.shape[0])[outer_tri_list])

    inner_tri_list = list(set(range(tri_simpl.shape[0])) - set(outer_tri_list))

    tri_out = tri_simpl[outer_tri_list, :]

    tri_in = tri_simpl[inner_tri_list, :]
    
    # list of points belonging to inner triangles
    tri_in1D = list(np.unique(tri_in.ravel()))
    # list of points belonginh to outer triangles
    tri_out1D = list(np.unique(tri_out.ravel()))
    
    # points belonging to both lists (crater area boundary points)
    edge_points = list(set(tri_in1D) & set(tri_out1D))


    try:

        from ASCtoSTLdict import top_smooth_flag

    except ImportError:

        print('Missing parameter in dict: top_smooth_flag (bool)')
        print('Set to default: top_smooth_flag = False')

    # set elevation of points to topography value
    for i in edge_points:
        # print(z_org[i]-z_new[i])

        if top_smooth_flag:

            z_new[i] = z_smooth[i]
            dz_total[i] = 0.0

            z_new_crater[i] = z_smooth[i]
            dz_crater[i] = 0.0

        else:
        
            z_new[i] = z_org[i]
            dz_total[i] = 0.0

            z_new_crater[i] = z_org[i]
            dz_crater[i] = 0.0
        

    for i in inner_tri_list:
    
        dz_tri = dz_total[int(tri_simpl[i,0])]
        dz_tri += dz_total[int(tri_simpl[i,1])]
        dz_tri += dz_total[int(tri_simpl[i,2])]
        
        if dz_tri >= -0.001:
        
            print('QUI',i,dz_tri)
        
            inner_tri_list.remove(i)
            outer_tri_list.append(i)
        
            

    tri_in = tri_simpl[inner_tri_list, :]    
    tri_out = tri_simpl[outer_tri_list, :]



    if top_smooth_flag:
        vertices_org = np.column_stack((x_org, y_org, z_smooth))
    else:
        # original elevation 3D point refined grid
        vertices_org = np.column_stack((x_org, y_org, z_org))

    # modified elevation 3D point refined grid
    vertices = np.column_stack((x_new, y_new, z_new))

    # modified elevation 3D point refined grid
    vertices_crater = np.column_stack((x_new, y_new, z_new_crater))

    # modified elevation 3D point refined grid
    vertices2 = np.column_stack((x_new, y_new, z_new2))

    # Create the full mesh
    print('')
    print('Saving full topography stl')

    surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices[f[j], :]

    output_dir = '../constant/triSurface'
    # Check whether the specified output path exists or not
    isExist = os.path.exists(output_dir)

    if not isExist:

        # Create a new directory because it does not exist
        os.makedirs(output_dir)
        print('The new directory ' + output_dir + ' is created!')

    surface.save('../constant/triSurface/surface.stl')

    # Create the inside mesh
    print('Saving inside crater topography stl')

    faces = tri_in

    surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices[f[j], :]

    surface.save('../constant/triSurface/surface_in.stl')

    surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices_org[f[j], :]

    surface.save('../constant/triSurface/surface_top.stl')


    # Create the outside mesh
    print('Saving out of crater topography stl')
    faces = tri_out

    surface = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices[f[j], :]

    surface.save('../constant/triSurface/surface_out.stl')

    # Create the inside mesh closed on top
    print('Saving closed surface total stl')

    faces = tri_in

    surface = mesh.Mesh(np.zeros(2 * faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices[f[2-j], :]
            surface.vectors[i + faces.shape[0]][j] = vertices_org[f[j], :]
  
    volume, cog, inertia = surface.get_mass_properties()
    print("Closed",surface.is_closed())
    print("Total Volume = {0}".format(volume))

    surface.save('../constant/triSurface/surface_total_closed.stl')


    # Create the inside mesh closed on top
    print('Saving closed surface crater stl')

    surface = mesh.Mesh(np.zeros(2 * faces.shape[0], dtype=mesh.Mesh.dtype))

    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices_crater[f[2-j], :]
            surface.vectors[i + faces.shape[0]][j] = vertices_org[f[j], :]

    volume, cog, inertia = surface.get_mass_properties()
    print("Closed",surface.is_closed())
    print("Total Volume = {0}".format(volume))

    surface.save('../constant/triSurface/surface_crater_closed.stl')



    # Create the inside mesh closed on top
    print('Saving lagrangian stl')

    faces = tri_lagr

    surface = mesh.Mesh(np.zeros(2 * faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            surface.vectors[i][j] = vertices2[f[2-j], :]            
            surface.vectors[i + faces.shape[0]][j] = vertices_org[f[j], :]

    volume, cog, inertia = surface.get_mass_properties()
    print("Closed",surface.is_closed())
    print("Lagrangian Volume = {0}".format(volume))

    surface.save('../constant/triSurface/surface_lagrangian.stl')

    if conduit_radius > 0:

        conduit_volume *= (X1[0, 1] - X1[0, 0])**2

        print(f"Conduit volume [m3] : {float(conduit_volume):e}")
        # Conduit stl

        points = []

        for i, (x, y) in enumerate(zip(x_conduit_top, y_conduit_top)):

            points.append([x, y])

        points = np.asarray(points)

        print('')
        print('Building Delaunay triangulation for conduit')

        tri = Delaunay(points)
        faces = tri.simplices
        
        for count,tri_n in enumerate(tri.neighbors):
        
            check_edge = False
        
            if tri_n[0]==-1:
            
                check_edge = True
                i1 = faces[count][1] 
                i2 = faces[count][2] 
                
            if tri_n[1]==-1:
            
                check_edge = True
                i1 = faces[count][0] 
                i2 = faces[count][2] 

            if tri_n[2]==-1:
            
                check_edge = True
                i1 = faces[count][0] 
                i2 = faces[count][1] 

            if check_edge:
            
                z_conduit_bottom[i1] = z_conduit_top[i1]
                z_conduit_bottom[i2] = z_conduit_top[i2]


        # modified elevation 3D point refined grid
        vertices_top = np.column_stack(
            (x_conduit_top, y_conduit_top, z_conduit_top))

        # modified elevation 3D point refined grid
        vertices_bottom = np.column_stack(
            (x_conduit_bottom, y_conduit_bottom, z_conduit_bottom))

        print('Saving closed surface conduit stl')

        surface = mesh.Mesh(np.zeros(2 * faces.shape[0],
                                     dtype=mesh.Mesh.dtype))
        for i, f in enumerate(faces):
            for j in range(3):
                surface.vectors[i][j] = vertices_top[f[j], :]
                surface.vectors[i +
                                faces.shape[0]][j] = vertices_bottom[f[2-j], :]
        volume, cog, inertia = surface.get_mass_properties()
        print("Closed",surface.is_closed())
        print("conduit Volume = {0}".format(volume))


        surface.save('../constant/triSurface/surface_conduit_closed.stl')
        

if __name__ == "__main__":
    main(sys.argv[1:])
