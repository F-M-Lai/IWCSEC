# -*- coding: utf-8 -*-
""" 
    This module defines the polygons.
    It requires modules: crystal_class, basis, CrystalLattice,
    ps_setd_mingled, mlab, tvtk and enclose_polyhedron.
"""
import numpy as np
import nlopt
from mayavi import mlab
from tvtk.api import tvtk 
import crystal_class as CC
from crystal_shape import ps_setd_mingled, enclose_polyhedron 
from Inverse_Wulff_construction import crystal_plane, surface_areas, optimal_distance_local
import intervals

    
def crystal_surface(intersection_point):  

    surf_point = [len(intersection_point)] # number of intersections  
    for i in np.arange(len(intersection_point)):
        surf_point.append(i)
        
    return surf_point  
    
def surface_distance(hkl):    

    '''
     	Distance between the center of the polyhedron and the surface 
    '''  
    
    surf_dist = []
    for surf in hkl:
        while True:
            d = input('Please input the distance of %s surface: ' % surf)
            if d < 0:
                print ''
                print 'Distance should be positive.'
                print 'Please input the correct value.'
                print ''            
            else:   
                surf_dist.append(d)
                break
    return surf_dist    

def surface_color(hkl):
  
    surf_color = []
    color_range = intervals.closed(0,1)
    for surf in hkl:
        while True:
            c = input('Please input the color of %s surface: ' % surf)
            if (type(c) is not tuple) or (c[0] not in color_range) or (c[1] not in color_range) or (c[2] not in color_range):
                print ''
                print 'The type of surface color should be tuple,'
                print 'and the elements in it should belong to the interval [0, 1].'
                print 'Please input the correct value.'
                print ''            
            else:   
                surf_color.append(c)
                break
    return surf_color    

def surface_index(): 
    
    while True: 
        hkl = input('Please input surface families: ')
        if (type(hkl) is not list) or (isinstance(hkl[0], int) == True): 
            print ''
            print 'The type of suface families should be list.'
            print 'Please input the correct value.'
            print ''
        else:             
            break   
    return hkl
    
#=========================================================================   
    
def main():   

    print ''      
    print '1. Format of surface families should be [[a1, b1, c1], [a2, b2, c2],...,[an, bn, cn]], n >= 1.'
    print '2. Distance between the center of the polyhedron and the surface should be positive.'
    print '3. Color (r, g, b) is represented by the number from 0 to 1.'
    print ''  
    hkl = surface_index()                 # surface family
    print ''
    surf_dist  = surface_distance(hkl)    # distance
    print ''
    surf_color = surface_color(hkl)       # surface color
    print 'surface_color', surf_color
    print ''     
    planes, number_planes = crystal_plane(hkl) 
    ps_setd_mingled(number_planes, planes, surf_dist) 
    diameter = 20.0      
    volume   = 4.0/3.0 * np.pi * (diameter*0.5)**3     
    qfp      = enclose_polyhedron(number_planes, planes, volume)
    p1       = tvtk.PolyData()  
    n_planes = len(qfp[2]) 
    for i in np.arange(len(number_planes)):
        for surf_index in np.arange(sum(number_planes[0:i]), sum(number_planes[0:i+1])):  
            print 'surf_index', surf_index 
            p1.points = qfp[2][surf_index]                    # the coordinates of the intersection of each plane.
            faces     = crystal_surface(qfp[2][surf_index])
            cells     = tvtk.CellArray()                      # create a new CellArray object to assign the polys property.
            cells.set_cells(1, faces)                         # the first parameter is the number of faces (here is 1), 	                              
            p1.polys  = cells                                 # and the second parameter is an array describing the composition of each face.
            p1.point_data.scalars = np.linspace(0.0, 1.0, len(p1.points))  
            mlab.figure(number_planes, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
            mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color = surf_color[i])        
	axe = tvtk.AxesActor(total_length=(3,3,3))
    mlab.show()
               
if __name__ == "__main__":  

    main()

    
    
    
    
    