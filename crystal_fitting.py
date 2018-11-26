# -*- coding: utf-8 -*-
""" 
    This program is mainly used to calculate the functional relationship between energy ratio and area fraction.
    It mainly requires modules: enclosure, crystal_plane, surface_areas, optimal_distance_local, area_change_trend, 
    bounded_bisecting, and total_surface_energy.
"""
import numpy as np
from cgeometry import enclosure
import nlopt
import matplotlib.pyplot as plt
from Inverse_Wulff_construction import crystal_plane, optimal_distance_local, surface_families_is_2, area_change_trend, bounded_bisecting, total_surface_energy
	  
def target_areas(lb1, ub1, n_point, hkl):
    '''
	area fraction of different surfaces.
    input:
        lb1, ub1: the upper and lower bounds of the area where the area fraction is located.
		n_point: number of divisions of the interval.   
		hkl: the surface couple.
	return: 
		frac_area_couple: the given fractional areas. 
        x_data: the fractional area that need to be fitted to the energy ratios.
    ''' 	
    frac_area = np.linspace(lb1, ub1, n_point)  
    frac_area_couple = []	
    for i in range(n_point):	
        frac_area_couple.extend([[frac_area[i], 1-frac_area[i]]])
        
    while True: 
        fitting_surface = input('Please select the surface that need to be fitted to the energy ratios: ')
        if fitting_surface != hkl[0] and fitting_surface != hkl[1]: 
            print ''
            print ('Please choose the surface in the %s.' % hkl)
            print ('Please input the correct value.')
            print ''
        elif fitting_surface == hkl[0]:
            x_data = frac_area
            break
        elif fitting_surface == hkl[1]:
            x_data = 1 - frac_area 
            break
        else: 
            break     
    
    return frac_area_couple, x_data, fitting_surface
	
def fitted_fun(a, x):
    '''
        this function is used to fit the energy ratios and the fractional area.
    '''
    
    return (1-x)*(a[0]-a[1]*np.sqrt(x/((a[0]-1)*x + 1))) + x*(1.0/(a[2]-a[3]*np.sqrt((1-x)/((a[2]-1)*(1-x) + 1))))
    		   
def _residual(a, x_data, y_data):

    return y_data - fitted_fun(a, x_data)

def _residual_sq(a, grad):

    r = y_data - fitted_fun(a, x_data)
	
    return sum(r*r)


if __name__=="__main__" :

    print ''      
    print '1. Format of surface couple should be [[a1, b1, c1], [a2, b2, c2]].'
    print ''
    # hkl = [[1,2,2], [1,0,0]]
    hkl = surface_families_is_2()
    print ''   
    lb1     = 0.0                       # lower bound of fractional area
    ub1     = 1.0                       # upper bound of fractional area
    n_point = 12                        # number of points that need to be fitted 
    energy_ratio_two = [0.1, 1.9]       # function's increase and decrease      
    frac_area_couple, x_data, fitting_surface = target_areas(lb1, ub1, n_point, hkl)	     
        
    planes, number_planes = crystal_plane(hkl)
    object_func    = total_surface_energy      # objective function
    local_optimize = nlopt.LN_COBYLA           # local optimization method
    lb_ub	       = [0.0, 2.0]                # iteration interval
    Iter_max       = 5000                      # the maximum number of iterations
    re_xtol        = 1e-6                      # relative tolerance on optimization parameters           
    abs_ftol       = 1e-6                      # absolute tolerance on function value
    d0             = [1.0, 1.0]                # initial value    
    dx             = 0.0001                    # step size      
    energy_ratio   = []
    area_trend     = area_change_trend(energy_ratio_two, object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx)      
    for area_target in frac_area_couple:	
        print 'area_trend', area_trend 
        print area_target
        crit_lb_ub, opt_dist = bounded_bisecting(lb_ub, area_trend, area_target, 1e-6, object_func, local_optimize, 
                                                 Iter_max, re_xtol, abs_ftol, d0, dx)        
        energy_ratio.extend([crit_lb_ub])                        
    print energy_ratio       

    print 'area is', x_data   
    y_data   = np.array(energy_ratio)  # the fitted function unsupported operand type(s) for -: 'int' and 'list'.
    Iter_max = 50000                
    re_xtol  = 1e-12
    re_ftol  = 1e-12
    d0       = [2, 1, 2, 1]
    dx       = 0.0001
    # Iter_max, re_xtol, re_ftol, d0, dx are for the fitting function "_residual_sq".
    # The local optimization method is still "COBYLA".
    fp = optimal_distance_local(_residual_sq, local_optimize, Iter_max, re_xtol, re_ftol, d0, dx)	

    print ''            
    print '************************************************ Parameters of fitting function ************************************************************' 
    print ''   
    print 'Coexisting_Surfaces_1: ', hkl[0]
    print 'Coexisting_Surfaces_2: ', hkl[1]
    print ''
    # print 'Fractional area of [100] surface:', frac_area
    print 'Surface energy-ratios of Gamma_%s/Gamma_%s:' % (hkl[0], hkl[1]), energy_ratio 
    print ''
    print 'Parameters [a0, a1, a2, a3] of fitting function: [%s, %s, %s, %s]'  % (fp[0][0], fp[0][1], fp[0][2], fp[0][3])
    print ''
    print '************************************************************* End ***************************************************************************'
    print '' 
    
    x1  = min(x_data)  
    x2  = max(x_data)
    xr  = x2-x1
    y1  = min(y_data)
    y2  = max(y_data)
    yr  = y2-y1
    x   = np.linspace(x1, x2, 500)
    y   = fitted_fun(fp[0], x)    
    plt.figure(figsize=(5,4))	    
    l2, = plt.plot(x, y, 'b-', linewidth=2.0)
    l1, = plt.plot(x_data, y_data, 'ro')

    plt.xlim(0, 1)
    # plt.ylim(0.55, 1.75) 
    plt.yticks([0.55, 0.85, 1.15, 1.45, 1.75])
    # plt.title("$(113)$") 
    S1 = bytes(hkl[0][0]*100+hkl[0][1]*10+hkl[0][2])
    S2 = bytes(hkl[1][0]*100+hkl[1][1]*10+hkl[1][2])
    S3 = bytes(fitting_surface[0]*100+fitting_surface[1]*10+fitting_surface[2])    
    plt.xlabel("Fractional Area, $A_{%s}/(A_{%s}+A_{%s})$" % (S3, S1, S2))
    plt.ylabel("Energy Ratio, $\gamma_{%s}/\gamma_{%s}$" % (S1, S2))	
    plt.legend(handles = [l1,l2], labels = ['Actual Values', 'Fitting Curve'], loc = 'upper right')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.savefig("%s-%s.png" %(S1, S2), dpi=200)
    plt.show()
    
   
    
    