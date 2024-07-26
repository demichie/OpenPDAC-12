import vtk
import os
import sys
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from linecache import getline
import pandas as pd
sys.path.insert(0,'./preprocessing')
from preprocessing.ASCtoSTLdict import *
from matplotlib.colors import LightSource
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker

toll = 1.0
step_dens = 50.0
    

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def readerVTK(filename):

    reader = vtk.vtkDataSetReader() 

    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.ReadAllFieldsOn()
    reader.Update()

    mydata = dsa.WrapDataObject(reader.GetOutput())
   
    origId = np.array(mydata.PointData['origId']).astype(int)
    ImpCpuId = mydata.PointData['lmpCpuId']
    rho = mydata.PointData['rho']
    d = mydata.PointData['d']
    U = mydata.PointData['U']

    Point_coordinates = reader.GetOutput().GetPoints().GetData()
    position = numpy_support.vtk_to_numpy(Point_coordinates)

    #print(d)
    #print(U)
    #print(position)
    #print(rho)
    #print(origId)

    return origId , d, U, position, rho


def readASC(DEM_file):

    print('Reading DEM file: ' + DEM_file)
    # Parse the topography header
    hdr = [getline(DEM_file, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    xs_DEM = lx + 0.5 * cell + np.linspace(0, (cols - 1) * cell, cols)
    ys_DEM = ly + 0.5 * cell + np.linspace(0, (rows - 1) * cell, rows)

    extent = lx, lx + cols * cell, ly, ly + rows * cell

    # Load the topography into a numpy array
    DEM = pd.read_table(DEM_file, delim_whitespace=True, header=None,
                    skiprows=6).astype(float).values
    DEM = np.flipud(DEM)
    DEM[DEM == nd] = 0.0

    xinit = xs_DEM
    yinit = ys_DEM

    xmin = np.amin(xinit)
    xmax = np.amax(xinit)

    print('xmin,xmax',xmin,xmax)

    ymin = np.amin(yinit)
    ymax = np.amax(yinit)

    print('ymin,ymax',ymin,ymax)

    Xinit, Yinit = np.meshgrid(xinit, yinit)
    Zinit = DEM
   
    return Xinit,Yinit,Zinit,cell,extent


def main():

    foamCommand = 'foamToVTK -fields "()" -useTimeName'
    os.system(foamCommand)

    print('DEM_file',DEM_file)
    filename = './preprocessing/'+DEM_file
    Xinit,Yinit,Zinit,cell,extent = readASC(filename)
 
 
    ls = LightSource(azdeg=45, altdeg=45)

    current_dir = os.getcwd()
    current_dir_name = current_dir.split('/')[-1]

    working_dir = current_dir+'/VTK/lagrangian/cloud'
    files = os.listdir(working_dir)
    n_files = len(files)
    file_idx = []
    for filename in files:
        # print(filename)
        file_idx.append(float(filename.split('_')[1].rsplit('.', 1)[0]))
        
    times = [float(x) for x in file_idx]  
      
    sorted_files = [files[i] for i in np.argsort(file_idx)]
    sorted_times = [times[i] for i in np.argsort(file_idx)]
    
    print('Times read',sorted_times)

    full_filename = working_dir+'/'+sorted_files[1]

    origId,d, U, position, rho = readerVTK(full_filename)

    diams = np.unique(d)
    
    nballistics = d.shape[0]
    n_times = n_files
    print('nballistics',nballistics)
    A = np.zeros((nballistics,3,n_times))
    B = np.zeros((nballistics,3,n_times))
    matr = np.zeros((n_times,8,nballistics))
    
    for i,filename in enumerate(sorted_files[:]):
    
        full_filename = working_dir+'/'+filename
        print(full_filename)
        origId,d, U, position, rho = readerVTK(full_filename)
        A[:,:,i] = U
        B[:,:,i] = position
        for k in range(nballistics):
            matr[i,1:4,k] = B[k,:,i]
            matr[i,4:7,k] = A[k,:,i]
            matr[i,7,k] = LA.norm(matr[i,4:7,k])
    
    
    # Determine impact time for all clasts as the first time index
    # such that velocity (in norm) falls below toll tolerance and, 
    # in the previous time step, velocity along z-axis is negative. 
    # Otherwise, set impact time to be zero  
    
    t_impact = np.zeros((nballistics,1)).astype(int) 
    time_impact = np.zeros((nballistics,1)).astype(float) 

    for s in range(nballistics):
        for l in range(3, n_times):
            if (matr[l,-1,s] < toll and matr[l-1,-2,s] < 0):
                t_impact[s] = int(l)
                time_impact[s] = sorted_times[l]
                print(s,l,time_impact[s])
                break

    ## Calculate mean and maximum velocities along particle trajectory
    velocities = np.zeros((nballistics, 4))
    for k in range(nballistics):
        velocities[k,0] = int(k)
        velocities[k,1] = d[k]
        velocity = matr[:int(t_impact[k])+1, -1, k]
        mean_velocity = np.mean(velocity)
        velocities[k,2] = mean_velocity 
        max_velocity = np.amax(velocity)
        velocities[k,3] = max_velocity
        
    C = ['index','diameter [m]','mean vel [m/s]','max vel [m/s]']
    df = pd.DataFrame(velocities) 
    df[0] = df[0].astype(int)
    df.to_csv("velocities.csv", header=C,index=False) 

    # Determine the mass of each block 
    r = d/2 
    V = 4/3*(np.pi)*(r**3) 
    m = rho*V 
    K = np.zeros((nballistics,1)) 

    mat1 = np.zeros((nballistics, 9))

    for s in range(nballistics):
        mat1[s,0] = s
        if t_impact[s] != 0:
            
            K[s] = 0.5*m[s]*(mat1[s,2]**2) 
            
            mat1[s,1] = d[s]
            mat1[s,2] = rho[s]
            mat1[s,3] = time_impact[s]
            mat1[s,4] = matr[t_impact[s],1,s]
            mat1[s,5] = matr[t_impact[s],2,s]
            mat1[s,6] = matr[t_impact[s],3,s]
            mat1[s,7] = matr[t_impact[s]-1,-1,s]
            mat1[s,8] = K[s]
            
    C = ['index','diameter [m]','density [kg/m3]','[impact time [s]','x [m]','y [m]','z [m]','impact velocity [m/s]','landing energy [J]']

    df = pd.DataFrame(mat1) 
    df[0] = df[0].astype(int)
    df.to_csv("impacts.csv", header=C,index=False) 
            
    x = np.array(position[:,0]) 
    y = np.array(position[:,1]) 
    diam = np.array(d)

    xy = np.vstack([x,y])
    kde = gaussian_kde(xy,bw_method=1.0)
    z = gaussian_kde(xy)(xy)
    
    xmin = np.amin(Xinit) - 0.5*cell   
    xmax = np.amax(Xinit) + 0.5*cell

    ymin = np.amin(Yinit) - 0.5*cell    
    ymax = np.amax(Yinit) + 0.5*cell    
    
    step_dens = 50.0
    x_density = np.arange(xmin,xmax,step=step_dens)
    y_density = np.arange(ymin,ymax,step=step_dens)

    x_dens_min = np.amin(x_density)
    x_dens_max = np.amax(x_density)
    y_dens_min = np.amin(y_density)
    y_dens_max = np.amax(y_density)
    
    extent_density = x_dens_min , x_dens_max , y_dens_min , y_dens_max
    
    print('extent_density',extent_density)

    xx,yy = np.meshgrid(x_density,y_density)
    
    count_ballistic_class = np.zeros((len(diams)+1,xx.shape[0],xx.shape[1]))
    zz = np.zeros_like(xx)
    
    for xi,yi,di in zip(x,y,diam):
    
        i = ( xi + xc - x_dens_min ) / step_dens
        j = ( yi + yc - y_dens_min ) / step_dens

        i = int(np.ceil(i))		
        j = int(np.ceil(j))
                
        count_ballistic_class[-1,j,i] +=1        
        k = int(np.argwhere(di==diams)[0]) 
        count_ballistic_class[k,j,i] +=1        
        
 
    for i in range(len(diams)+1):
     
        fig, ax = plt.subplots()
    
        im = ax.imshow(ls.hillshade(np.flipud(Zinit),vert_exag=1.0,
                     dx=cell,dy=cell),cmap='gray',extent=extent)
        
        ax.set_aspect('equal', 'box')
   
        zz[:,:] = np.squeeze(count_ballistic_class[i,:,:])
        sum_zz = np.nansum(zz)
        print('sum_zz',sum_zz)
        zz = zz / np.sum(zz) * 100.0
        zz = np.log10(zz)   
    
        min_arr = np.amin(zz)
        max_arr = np.amin(zz)
    
        zz_max = np.amax(zz)
        zz_linspace = np.linspace(0,zz_max,num=11)
        ticks = []
        for val in zz_linspace:
            ticks.append(str(val))
            
        if 'domain_size_x' in locals(): 

            # values for blockMeshDict
            xmin = np.maximum(xmin,xc-0.5 * domain_size_x)
            xmax = np.minimum(xmax,xc+0.5 * domain_size_x)

        if 'domain_size_y' in locals(): 

            # values for blockMeshDict
            ymin = np.maximum(ymin,yc-0.5 * domain_size_y)
            ymax = np.minimum(ymax,yc+0.5 * domain_size_y)    
    
        im_ratio = (ymax-ymin) / (xmax-xmin)
                   
        levels = np.linspace(min_arr, max_arr, 11)
        label_str = 'Probabilty [0;1]'

        cmap = plt.get_cmap('terrain_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
       
        # ax.scatter(x+xc,y+yc,c=z,s=3,alpha=0.1,edgecolors='none')
    
        p1 = ax.imshow(np.flipud(zz),cmap=cmap, interpolation='nearest', 
                       extent=extent_density, alpha=0.65)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        clb = plt.colorbar(p1)

        label_str = 'Log of % ballistic'
        clb.set_label(label_str, labelpad=-40, y=1.05, rotation=0)
      
        if i<len(diams):
                
            string = '_d'+str(i)+'_'
            title = 'Diameter = '+str(diams[i])+'m' 
            
        else:
        
            string = '_tot_'
            title = 'All sizes'
      
        plt.title(title)
      
        png_file = current_dir_name+string+'ballistic.png'
             
        plt.savefig(png_file, dpi=200)
        plt.close(fig)
 
        nd = -9999
 
        zz[zz==-np.inf] = nd   		

        asc_file = current_dir_name+string+'ballistic.asc'
    
        header = "ncols     %s\n" % xx.shape[1]
        header += "nrows    %s\n" % xx.shape[0]
        header += "xllcorner " + str(x_dens_min) + "\n"
        header += "yllcorner " + str(y_dens_min) + "\n"
        header += "cellsize " + str(step_dens) + "\n"
        header += "NODATA_value " + str(nd)

        np.savetxt(asc_file,
               np.flipud(zz),
               header=header,
               fmt='%1.5f',
               comments='')

        # plt.show()

    
if __name__ == '__main__':

    main()

