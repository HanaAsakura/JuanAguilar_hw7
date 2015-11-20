from matplotlib import cm
from matplotlib import pyplot
import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import griddata

arr = np.loadtxt(sys.argv[0])
a = arr[0]
b= arr[1]
c= arr[2]
d=arr[3]
l=arr[4]

min_a = np.amin(a)
max_a = np.amax(a)
min_b = np.amin(b)
max_b = np.amax(b)
min_c = np.amin(c)
max_c = np.amax(c)
min_d = np.amin(d)
max_d = np.amax(d)

with PdfPages('mcmc_lotka_volerra_graphs.pdf') as pdf:
    n_points = np.size(a)
    points = np.ones((n_points,2))
    #print shape(points)
    points[:,0] = a
    points[:,1] = b
    grid_a, grid_b = np.mgrid[min_a:max_a:200j, min_b:max_b:200j]
    grid_l = griddata(points, -np.log(l), (grid_a, grid_b), method='cubic')
    imshow(grid_l.T, extent=(min_a,max_a,min_b,max_b), aspect='auto',origin='lower')
    plt.title("a,b")
    pdf.savefig()
    plt.close()
    
    points = np.ones((n_points,2))
    #print shape(points)
    points[:,0] = a
    points[:,1] = c
    grid_a, grid_c = np.mgrid[min_a:max_a:200j, min_c:max_c:200j]
    grid_l = griddata(points, -np.log(l), (grid_a, grid_c), method='cubic')
    imshow(grid_l.T, extent=(min_a,max_a,min_c,max_c), aspect='auto',origin='lower')
    plt.title("a,c")
    pdf.savefig()
    plt.close()
    
    points = np.ones((n_points,2))
    #print shape(points)
    points[:,0] = a
    points[:,1] = d
    grid_a, grid_d = np.mgrid[min_a:max_a:200j, min_d:max_d:200j]
    grid_l = griddata(points, -np.log(l), (grid_a, grid_d), method='cubic')
    imshow(grid_l.T, extent=(min_a,max_a,min_d,max_d), aspect='auto',origin='lower')
    plt.title("a,d")
    pdf.savefig()
    plt.close()
    
    points = np.ones((n_points,2))
    #print shape(points)
    points[:,0] = b
    points[:,1] = c
    grid_b, grid_c = np.mgrid[min_b:max_b:200j, min_c:max_c:200j]
    grid_l = griddata(points, -np.log(l), (grid_b, grid_c), method='cubic')
    imshow(grid_l.T, extent=(min_b,max_b,min_c,max_c), aspect='auto',origin='lower')
    plt.title("b,c")
    pdf.savefig()
    plt.close()
    
    points = np.ones((n_points,2))
    #print shape(points)
    points[:,0] = b
    points[:,1] = d
    grid_b, grid_b = np.mgrid[min_b:max_b:200j, min_d:max_d:200j]
    grid_l = griddata(points, -np.log(l), (grid_b, grid_d), method='cubic')
    imshow(grid_l.T, extent=(min_b,max_b,min_d,max_d), aspect='auto',origin='lower')
    plt.title("b,d")
    pdf.savefig()
    plt.close()
    
    points = np.ones((n_points,2))
    #print shape(points)
    points[:,0] = c
    points[:,1] = d
    grid_c, grid_d = np.mgrid[min_c:max_c:200j, min_d:max_d:200j]
    grid_l = griddata(points, -np.log(l), (grid_c, grid_d), method='cubic')
    imshow(grid_l.T, extent=(min_c,max_c,min_d,max_d), aspect='auto',origin='lower')
    plt.title("c,d")
    pdf.savefig()
    plt.close()
    
    
