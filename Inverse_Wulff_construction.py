# -*- coding: utf-8 -*-
""" 
    This program is used to calculate the critical free energy ratios. 
    In the optimization part we use the local optimization method.
    It requires modules: enclosure, crystal_class, basis and crystal_shape.
"""
import numpy as np
from cgeometry import enclosure
import nlopt
import crystal_class as CC
import basis
from crystal_lattice import CrystalLattice
from crystal_shape import ps_setd_mingled

def crystal_plane(hkl):
    '''
       return the geometric planes of the polyhedron, where the 32 point group is considered, 
	   and the lattice parameter is taken as 0.4nm.
	   hkl: Miller indices of crystallographic planes.
    ''' 
    global number_planes, planes    
    b1     = basis.cubic(0.4)   
    cc1    = CC.CrystalClass32()
    latt   = CrystalLattice(b1, cc1)
    planes = []
    number_planes = []
    for index in hkl:
        pf = latt.geometric_plane_family(index)
        number_planes.append(len(pf))
        planes.extend(pf)
		
    return planes, number_planes

def surface_areas(d):
    '''
    input:
        origin-to-planes distances.
		
	    volume_calc: the natural or unscaled volume of the polyhedron.
		area_p: areas of all the planes.
		points_p: interception points on the planes, grouped by plane.   
		
	return: 
		area_ps: the area of each plane family (sum of areas of surfaces of a plane family). 
    '''  
    global Gamma_ratio
    ps_setd_mingled(number_planes, planes, d)
    planes_c = [p.as_list() for p in planes]
    pr       = enclosure(planes_c)   
    volume_calc = pr[0]                          # the natural or unscaled volume of the polyhedron.
    area_p      = pr[1]                          # areas of all the planes.
    points_p    = pr[2]                          # interception points on the planes, grouped by plane.                                       
    area_ps     = []                             # area of each plane family: sum of areas of surfaces of a plane family.               
    index_start = 0
    for subset_n in number_planes:
        index_end = index_start + subset_n
        area_ps.append(np.sum(area_p[index_start : index_end]))
        index_start = index_end 
    global volume
    diameter = 20.0 
    volume   = 4.0/3.0*np.pi* (diameter*0.5)**3   # this volume was used to scale the calculated volume.
    area_ps  = np.array(area_ps)*np.power(volume/volume_calc, 2.0/3) 
	
    return area_ps
	
def area_fraction_fun(d,grad):
    '''
     	objective function of area fraction.
    '''  
    areas      = surface_areas(d)
    areas_fra  = areas/sum(areas)                                 # area fraction. 
    areas_diff = np.array(area_target) - np.array(areas_fra)      # the difference between the target area and the optimized area.
    areas_squ  = sum([x*x for x in areas_diff])                   # the square of the area difference.
    # print 'Fractional surface area = ', areas_fra
	
    return areas_squ    
	
def optimal_distance_local(object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx):
    ''' The NLopt (Non-Linear Optimization) library is used for the optimization problem.
	    Local optimization methods include COBYLA, BOBYQA, PRAXIS, etc. Here we use COBYLA.	
    input:
        Iter_max: the maximum number of iterations.
        re_xtol；relative tolerance on optimization parameters.
        abs_ftol：absolute tolerance on function value.	
        d0：initial value.  		
        dx: initial steps.
	return: 
		optd: the optimized values of the optimization parameters.
		enumerated_constant:  enumerated constant which takes on 1, 2, 3, 4, 5, and 6.		
		minf: the function value corresponding to optd.
    '''  

    opt = nlopt.opt(local_optimize, len(d0))  
    opt.set_min_objective(object_func)    
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    opt.set_initial_step(dx)
    optd = opt.optimize(d0)
    minf = opt.last_optimum_value()
    enumerated_constant = opt.last_optimize_result()	
    if enumerated_constant   == 1:
        print 'NLOPT_SUCCESS = 1, Generic success return value.'
    elif enumerated_constant   == 2:
        print 'NLOPT_SUCCESS = 2, stopval (above) was reached.'
    elif enumerated_constant   == 3:
        print 'NLOPT_SUCCESS = 3, ftol_rel or ftol_abs (above) was reached.'
    elif enumerated_constant   == 4:
        print 'NLOPT_SUCCESS = 4, xtol_rel or xtol_abs (above) was reached.'
    elif enumerated_constant   == 5:
        print 'NLOPT_SUCCESS = 5, maxeval (above) was reached.'
    # elif enumerated_constant == 6:
    else:
        print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
 
    return optd, enumerated_constant, minf
	
def optimal_distance_global(object_func, gloabl_optimize, Iter_max, re_xtol, abs_ftol, d0, dx, lb, ub):
    '''
     	This is the global optimization algorithm in NLopt, mainly including
		ESCH, ISRES, MLSL, etc.
    '''
    opt = nlopt.opt(gloabl_optimize, len(d0))  
    opt.set_min_objective(object_func)          
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    #opt.set_ftol_rel(1e-10)
    opt.set_initial_step(dx)
	#a = opt.get_numevals()
    opt.set_lower_bounds(lb)
    opt.set_upper_bounds(ub)
    opt.set_population(100)
    opt.set_local_optimizer(nlopt.opt(nlopt.LN_BOBYQA, len(d0)))
    optd   = opt.optimize(d0)
    minf   = opt.last_optimum_value()
    enumerated_constant = opt.last_optimize_result()
    if enumerated_constant   == 1:
        print 'NLOPT_SUCCESS = 1, Generic success return value.'
    elif enumerated_constant   == 2:
        print 'NLOPT_SUCCESS = 2, stopval (above) was reached.'
    elif enumerated_constant   == 3:
        print 'NLOPT_SUCCESS = 3, ftol_rel or ftol_abs (above) was reached.'
    elif enumerated_constant   == 4:
        print 'NLOPT_SUCCESS = 4, xtol_rel or xtol_abs (above) was reached.'
    elif enumerated_constant   == 5:
        print 'NLOPT_SUCCESS = 5, maxeval (above) was reached.'
    # elif enumerated_constant == 6:
    else:
        print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
    
    return optd, enumerated_constant, minf

def sum_area_is_1():
    '''
     	Determine the sum of the input fractional area is 1.
    '''
    while True: 
        area_target = input('Please input the fractional surface area [f1, f2]: ')
        if len(area_target) != 2 or sum(area_target) != 1 or area_target[0] < 0 or area_target[1] < 0: 
            print ''
            print ('The sum of f1 and f2 should be 1, and 0 <= f1,f2 <= 1.')
            print ('Please input the correct value.')
            print ''
        else: 
            break 

    return area_target
    
def surface_families_is_2():
    '''
     	Determine the number of surface family is 2.
    '''
    while True: 
        hkl = input('Please input the surface couple: ')
        if len(hkl) != 2 or (type(hkl) is not list): 
            print ''
            print ('hkl needs to contain two surface, and the type should be a list.')
            print ('Please input the correct value.')
            print ''
        else: 
            break 

    return hkl 

def total_surface_energy(d, grad):
    '''
     	objective function for surface energy ratio.
    '''

    areas = surface_areas(d)
    # print 'Gamma_ratio', Gamma_ratio
    return Gamma_ratio*areas[0] + areas[1]

def area_change_trend(energy_ratio_two, object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx):
    '''
     	Determine the trend of fractional area as a function of surface energy.
        
        energy_ratio_two: Two numbers used to determine the function's increase and decrease.
        area_trend: Fractional area change trend      
    '''
    area_trend = []
    global Gamma_ratio
    for Gamma_ratio in energy_ratio_two:       
        res    = optimal_distance_local(object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx)
        areas_trend  = surface_areas(res[0]) 
        area_trend.extend([areas_trend[0]])
	
    return area_trend  

def bounded_bisecting(lb_ub, area_trend, area_target, rel_tol, object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx):
    '''
     	Bounded bisecting minimization method for surface energy-ratios.
        
        lb_ub: Upper bound and lower bound of the iteration interval.
        area_target: Fractional area change trend
        rel_tol: Relative tolerance of bounded bisecting        
    '''
    k  = 0         # number of iterations
    lb = lb_ub[0]  # lower bound of iteration interval
    ub = lb_ub[1]  # upper bound of iteration interval
        
    while abs(ub - lb)/ub > rel_tol:		    
        Gamma_ratio = (lb + ub)/2.0 
        # print 'Gamma_ratio-----1', Gamma_ratio        
        k     = k + 1
        global Gamma_ratio 
        res   = optimal_distance_local(object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx)
        # print 'Gamma_ratio-----2', Gamma_ratio
        areas = surface_areas(res[0]) 
        areas = areas/sum(areas)
        # print 'fractional surface area is', areas
        if area_target != [1.0, 0.0] and area_target != [0.0, 1.0]:
            if (area_trend[0] >= area_trend[1]  and areas[0] >= area_target[0]) or (area_trend[0] <= area_trend[1] and areas[0] <= area_target[0]):                 		
                lb = Gamma_ratio
            elif (area_trend[0] > area_trend[1] and areas[0] < area_target[0])  or (area_trend[0] < area_trend[1]  and areas[0] > area_target[0]): 	               	
                ub = Gamma_ratio              
        elif area_target == [1.0, 0.0]:
            if (area_trend[0] > area_trend[1]   and areas[0] < area_target[0])  or (area_trend[0] < area_trend[1]  and areas[0] == area_target[0]):
                ub = Gamma_ratio 
            elif (area_trend[0] < area_trend[1] and areas[0] < area_target[0])  or (area_trend[0] > area_trend[1]  and areas[0] == area_target[0]):
                lb = Gamma_ratio
        elif area_target == [0.0, 1.0]:
            if (area_trend[0] > area_trend[1]   and areas[0] > area_target[0])  or (area_trend[0] < area_trend[1]  and areas[0] == area_target[0]):
                lb = Gamma_ratio  
            elif (area_trend[0] < area_trend[1] and areas[0] > area_target[0])  or (area_trend[0] > area_trend[1]  and areas[0] == area_target[0]):
                ub = Gamma_ratio 
        print 'lb = ', lb, 'ub = ', ub                
        print 'iteration number = ', k 
        crit_lb_ub = round((ub + lb)/2.0, 3)      # retain three decimal places.
        opt_dist   = res[0]	
        print crit_lb_ub  
        print '--------------------------------------------'
        
    return crit_lb_ub, opt_dist    
    
#=========================================================================   
    
def main():

    print ''
    print '******************************************************************************************'
    print '*                                       WELCOME TO                                       *'
    print '*                                         IWCSEC                                         *'
    print '*                 Inverse Wulff Construction for Surface Energy Calculation              *'
    print '*                                   Version IWCSEC.0.1.0                                 *'    
    print '*                                                                                        *'    
    print '*                              Written by Haibo Guo and Fuming Lai                       *'
    print '*                            guohaibo@shu.edu.cn, fuming7710@163.com                     *'
    print '******************************************************************************************'

    print ''      
    print '1. Format of surface couple should be [[a1, b1, c1], [a2, b2, c2]].'
    print '2. For a given fractional area [f1, f2], the sum of f1 and f2 is 1, i.e. f1 + f2 = 1.'
    print '3. Fractional areas should satisfy 0 <= f1,f2 <= 1.'
    print '4. Calculate the case of only one surface familiy, the fractional area should be set to' 
    print '   [1, 0] or [0, 1]. At this point, the obtained energy-ratio is a limit.'  
    print ''        
    
    hkl = surface_families_is_2()
    global planes, number_planes, area_target, Gamma_ratio 
    area_target = sum_area_is_1()
    print '\n'  

    planes, number_planes = crystal_plane(hkl) 
    object_func      = total_surface_energy      # objective function
    local_optimize   = nlopt.LN_COBYLA           # local optimization method
    lb_ub	         = [0.0, 2.0]                # iteration interval
    Iter_max         = 5000                      # the maximum number of iterations
    re_xtol          = 1e-12                     # relative tolerance on optimization parameters           
    abs_ftol         = 1e-12                     # absolute tolerance on function value
    d0               = [1.0, 1.0]                # initial value    
    dx               = 0.0001                    # step size 	
    energy_ratio_two = [0.1, 1.9]                # function's increase and decrease
    rel_tol          = 1e-6                      # Relative tolerance of bounded bisecting
    area_trend       = area_change_trend(energy_ratio_two, object_func, local_optimize, Iter_max, re_xtol, abs_ftol, d0, dx)   
    print 'area_trend', area_trend 
    crit_lb_ub, opt_dist = bounded_bisecting(lb_ub, area_trend, area_target, rel_tol, object_func, local_optimize,
                                             Iter_max, re_xtol, abs_ftol, d0, dx)    
    
    print ''            
    print '******************************* Surface-energy ratio **************************************' 
    print ''   
    print 'Coexisting_Surfaces_1: ', hkl[0], '     Fractional surface area: ', area_target[0]
    print 'Coexisting_Surfaces_2: ', hkl[1], '     Fractional surface area: ', area_target[1] 
    print ''       
    print 'Surface-energy ratio of Gamma_%s/Gamma_%s: '  % (hkl[0], hkl[1]), crit_lb_ub
    print 'Distance between the center of the polyhedron and the surface:', opt_dist
    print 'Calculated fractional area: ', surface_areas(opt_dist)/sum(surface_areas(opt_dist))
    print ''
    print '************************************** End ************************************************'
    print ''
    
if __name__=="__main__" :
	
    main()


