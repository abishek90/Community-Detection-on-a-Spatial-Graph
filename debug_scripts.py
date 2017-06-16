# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 09:08:45 2017

@author: abishek
"""

# This script contains several functions useful in understanding and tuning the parameters of the algorithm contained in simulation.py. The key reason for a seperate debug  script was to see how the pre-defined parameter of graph pruning affects component size. 

# A Script to test the maximum connected component size.
def return_max_T_Good_component(Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
    maxc = 0
    maxnc = 0
    currc = 0
    grids_uc = Z_Grid.keys()
    numc = 1
    
    for zk in Z_Grid.keys():
        grids_uc.remove(zk)
        if Is_T_Good_PP(zk,Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
            curr_grid = zk
            grid_stack = [curr_grid]
            break
            
    
    #curr_grid = grids_uc[0]
    #grid_stack = [curr_grid]
    #grids_uc.remove(curr_grid)
    
    while grids_uc != [] :
        currc = 1
        currnc = len(Z_Grid[curr_grid])
        while grid_stack != []:
            
            for i in Z_Grid_nghbs[curr_grid]:
                if i in Z_Grid.keys() and i in grids_uc:
                    grids_uc.remove(i)
                    if Is_T_Good_PP(i,Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
                        grid_stack = grid_stack + [i]
                        currc += 1
                        currnc += len(Z_Grid[i])
                    
            grid_stack.remove(curr_grid)
            if grid_stack != []:
                curr_grid = grid_stack[0]
        maxc = max(maxc,currc)
        maxnc = max(maxnc,currnc)
        if currc > 1:
            print currc
            
        for zk in grids_uc:
            grids_uc.remove(zk)
            if Is_T_Good_PP(zk,Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
                curr_grid = zk
                grid_stack = [curr_grid]
                break
            
    #    curr_grid = grids_uc[0]
    #    grids_uc.remove(curr_grid)
    #    grid_stack = [curr_grid]
        numc+=1
    
    return maxc
    
def connected_component():
    maxc = 0
    maxnc = 0
    currc = 0
    grids_uc = Z_Grid.keys()
    numc = 1
    
#    for zk in Z_Grid.keys():
#        grids_uc.remove(zk)
#        if Is_T_Good_PP(zk,Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
#            curr_grid = zk
#            grid_stack = [curr_grid]
#            break
            
    
    curr_grid = grids_uc[0]
    grid_stack = [curr_grid]
    grids_uc.remove(curr_grid)
    
    while grids_uc != [] :
        currc = 1
        currnc = len(Z_Grid[curr_grid])
        while grid_stack != []:
            
            for i in Z_Grid_nghbs[curr_grid]:
                if i in Z_Grid.keys() and i in grids_uc:
                    grids_uc.remove(i)
                    grid_stack  = grid_stack + [i]
                    currc +=1
                    currnc += len(Z_Grid[i])
#                    if Is_T_Good_PP(i,Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
#                        grid_stack = grid_stack + [i]
#                        currc += 1
#                        currnc += len(Z_Grid[i])
                    
            grid_stack.remove(curr_grid)
            if grid_stack != []:
                curr_grid = grid_stack[0]
        maxc = max(maxc,currc)
        maxnc = max(maxnc,currnc)
        if currc > 1:
            print currc
            
#        for zk in grids_uc:
#            grids_uc.remove(zk)
#            
#            if Is_T_Good_PP(zk,Z_Grid,Z_Grid_nghbs,graphh_PP,colors):
#                curr_grid = zk
#                grid_stack = [curr_grid]
#                break
            
        curr_grid = grids_uc[0]
        grids_uc.remove(curr_grid)
        grid_stack = [curr_grid]
        numc+=1
    
    return maxc
    

def test_A_goodness():
    nagood = 0
    nagoodcells = 0
    nabad = 0
    nabadcells = 0
    for zk in Z_Grid.keys():
        if Is_Good_PP(zk,Z_Grid,Z_Grid_nghbs,graphh_PP):
            nagood += 1
            nagoodcells += len(Z_Grid[zk])
        else:
            nabad +=1 
            nabadcells += len(Z_Grid[zk])
            
 def rechangeZ():
     Z_Grid = create_Grid_String(2,N,R,locations)
     Z_Grid_nghbs = create_neighbors_Z_Grid(Z_Grid,2)
     est1 = community_Detection(graphh_PP,Z_Grid,Z_Grid_nghbs)
     est2 = community_Detection2(graphh_PP,Z_Grid,Z_Grid_nghbs)
     print np.absolute(sum(est1*colors)/num_nodes)
     print np.absolute(sum(est2*colors)/num_nodes)
    
