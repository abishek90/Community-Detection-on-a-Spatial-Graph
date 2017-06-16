# -*- coding: utf-8 -*-
"""
Created on Sun May 28 16:16:46 2017

@author: abishek
"""
# A Simulation of the Algorithm on synthetic data.

from __future__ import division
import numpy as np
import pandas as pd
import math
import random
from numpy import linalg as LA
import matplotlib.pyplot as plt



def Z_projection(R,d,x):
    z_proj = np.empty((d,1));
    z_proj_str = None
    #mulc = R/(4*math.pow(d,1/d));
    mulc = R*0.6;
    for dims in range(d):
        x_test1 = math.ceil(x[dims]/mulc)
        x_test2 = math.floor(x[dims]/mulc)
        
        if (np.absolute(x_test1*mulc - x[dims]) < 0.5001*mulc) :
            z_proj[dims] = x_test1
            if(z_proj_str == None):
                z_proj_str = str(int(x_test1)) + " "
            else:
                z_proj_str += str(int(x_test1)) + " ";
            
        elif (np.absolute(x_test2*mulc - x[dims]) < 0.5001*mulc) :
            if (z_proj_str == None):
                z_proj_str = str(int(x_test2)) + " "
            else:
                z_proj_str += str(int(x_test2)) + " ";
            z_proj[dims] = x_test2
            
        else:  
            print "Trouble"
           
    return z_proj_str


def create_Grid_String(d,N,R,locations):
    Z_Grid = {}
    num_points = max(locations.shape);
    for i in range(num_points):
        z_c = Z_projection(R,d,locations[i]);
        if z_c in Z_Grid:
            Z_Grid[z_c].append(i);
        else:
            Z_Grid[z_c] = [i];
    return Z_Grid
    
    
def create_neighbors_Z_Grid(Z_Grid,dimensions):
    Z_Grid_nghbs = {}
    for ks in Z_Grid.keys():
        Z_Grid_nghbs[ks] = []
        ks_val_lists = ks.split()
        for dims in range(dimensions):
            if ks in Z_Grid_nghbs :
                ks_val_lists[dims]= str(int(ks_val_lists[dims]) + 1)
                Z_Grid_nghbs[ks].append(" ".join(ks_val_lists) + " ")
                ks_val_lists[dims] = str(int(ks_val_lists[dims]) - 2)
                Z_Grid_nghbs[ks].append(" ".join(ks_val_lists) + " ")
                ks_val_lists[dims] = str(int(ks_val_lists[dims]) + 1)
    return Z_Grid_nghbs
    
        

def return_Threshold(x,y,R,a,b,lamb): # The Model hyperparameter for constant connection probability
    dxy = LA.norm(x-y);
    if dxy > 2*R:
        print dxy
    area_intersection = (2*R*R*math.acos(dxy/(2*R))) - ((dxy/2)*math.sqrt((4*R*R) - (dxy*dxy)));
    return lamb*area_intersection*((a+b)/2)*((a+b)/2)
    
    
############################################################################
    
def Pairwise_Classify(i,j,graphh,locations,N,R,a,b,lamb):
    num_nodes = max(graphh.shape)

    #dxy = LA.norm(locations[i]-locations[j]);
    dxy = return_distance(locations[i],locations[j],N,2)
    area_intersection = (2*R*R*math.acos(dxy/(2*R))) - ((dxy/2)*math.sqrt((4*R*R) - (dxy*dxy)));

    nsame = 0;
    nmisses = 0;
    nperks = 0;
    for k in range(num_nodes):
        #if k != i and k !=j and graphh[i][k] == 1 and graphh[j][k] == 1 :
        dik = LA.norm(locations[i] - locations[k])
        djk = LA.norm(locations[j] - locations[k])
        if dik > R or djk > R:
            continue
        if  graphh[i][k] == 1 and graphh[j][k] == 1 :
            nsame += 1
        if graphh[i][k] == 0 and graphh[j][k] == 0:
            nmisses += 1
        if graphh[i][k] + graphh[j][k] == 1:
            nperks +=1
    nsame = max(1,nsame)
    nmisses = max(1,nmisses)
    nperks = max(1,nperks)
    mlest = -(math.pow(a-b,2))*lamb*area_intersection + nsame*math.log((a*a + b*b)/(2*a*b))  + nmisses*math.log((math.pow(1-a,2) + math.pow(1-b,2))/(2*(1-a)*(1-b))) + (math.pow(a-b,2)/2)*lamb*area_intersection + nperks*math.log((a*(1-a) + b*(1-b))/(a*(1-b) + b*(1-b))); 
    
    if mlest > -0.5 :
        return 1
    else:
        return -1
    
#    lm = -1    
#    if mlest > -0.5:
#        lm = 1
#    #else:
##        #return -1
#    lr = -1
#    if nsame > return_Threshold(locations[i],locations[j],R,a,b,lamb):
#       lr = 1
#    return lm*lr
#    if ret == 0:
        #return -1
    #else:
        #return 1
########################################################################################

def Is_Good(z_key, Z_Grid, Z_Grid_nghbs,graphh,locations,N,R,a,b,lamb):
    if z_key in Z_Grid.keys():
        
        node_cons_list = Z_Grid[z_key];
        for zns in Z_Grid_nghbs[z_key]:
            if zns in Z_Grid.keys():
                node_cons_list = node_cons_list + Z_Grid[zns]
        # Now Test for Triangle inconsistencies
        val = 0
        
        for i in range(len(node_cons_list)):
            for j in range(i+1,len(node_cons_list)):
                for k in range(j+1,len(node_cons_list)):
                    T1 = Pairwise_Classify(node_cons_list[i],node_cons_list[j],graphh,locations,N,R,a,b,lamb);
                    T2 = Pairwise_Classify(node_cons_list[j],node_cons_list[k],graphh,locations,N,R,a,b,lamb);
                    T3 = Pairwise_Classify(node_cons_list[i],node_cons_list[k],graphh,locations,N,R,a,b,lamb);
                    if T1*T2*T3 == -1 :
                        val = -1
                        break
                if val == -1 :
                    break
            if val == -1:
                break
        if val == -1:
            return False
        else:
            return True
                    
        
        
    else:
        return False
        
        
def Is_Good_PP(z_key,Z_Grid,Z_Grid_nghbs,graphh_PP):
    if z_key in Z_Grid.keys():
        node_cons_list = Z_Grid[z_key]
        for zns in Z_Grid_nghbs[z_key]:
            if zns in Z_Grid.keys():
                node_cons_list = node_cons_list + Z_Grid[zns]
        val = 0
        for i in range(len(node_cons_list)):
            for j in range(i+1,len(node_cons_list)):
                for k in range(j+1,len(node_cons_list)):
                    
                    T1 = graphh_PP[node_cons_list[i]][node_cons_list[j]]
                    T2 = graphh_PP[node_cons_list[j]][node_cons_list[k]]
                    T3 = graphh_PP[node_cons_list[k]][node_cons_list[i]]
                    if T1*T2*T3 == -1:
                        val = -1
                        break
                if val == -1:
                    break
            if val == -1:
                break
        if val == -1:
            return False
        else:
            return True
    else:
        return False
                
   
def Is_T_Good_PP(z_key,Z_Grid,Z_Grid_nghbs,graph_PP,colors):
    if z_key in Z_Grid.keys():
        node_cons_list = Z_Grid[z_key]
        for zns in Z_Grid_nghbs[z_key]:
            if zns in Z_Grid.keys():
                node_cons_list = node_cons_list + Z_Grid[zns]
        val = 0
        for i in range(len(node_cons_list)):
            for j in range(i+1,len(node_cons_list)):
                ii = node_cons_list[i]
                jj = node_cons_list[j]
                if graphh_PP[ii][jj] != colors[ii]*colors[jj] :
                    val = -1
                    break
            if val == -1:
                break
            
                
                
        if val == -1: 
            return False
        else:
            return True
    else:
        return False
     
   
def community_Detection(graphh_PP,Z_Grid,Z_Grid_nghbs):
    num_nodes = max(graphh_PP.shape)
    estimates = np.empty((num_nodes,1))
    for i in range(num_nodes):
        estimates[i] = 1
    grids_unconsidered = Z_Grid.keys()
   
    curr_grid = None
    for ks in Z_Grid.keys():
        grids_unconsidered.remove(ks)
        curr_grid = ks
        break
        #if Is_Good(ks,Z_Grid,Z_Grid_nghbs,graphh,locations,R,a,b,lamb):
            #curr_grid = ks
            #break
    if curr_grid == None:
        return estimates
 
    curr_node = Z_Grid[curr_grid][0]
    curr_component = [] # Finished Partitions
    grid_stack = [curr_grid]
    
    max_c_size = 0;
    
    while grids_unconsidered != []:
        
        curr_c_size = 0
        while grid_stack != [] :
            # Produce a complete clustering of the nodes in the connectecd component of the current grid
            # Breadth First Search 
            # Partition points in curr_grid using curr_node
            for pps in Z_Grid[curr_grid]:
                if pps == curr_node:
                    continue
                else:
                    estimates[pps] = graphh_PP[pps][curr_node]*estimates[curr_node]
            #curr_component.append(curr_grid)
            curr_component = curr_component + [curr_grid]
            #print grids_unconsidered
            
            for ggs in Z_Grid_nghbs[curr_grid]:
                
                if ggs in grids_unconsidered  :
                       grids_unconsidered.remove(ggs)
                       grid_stack.append(ggs)
                       
                       #if Is_Good(ggs,Z_Grid,Z_Grid_nghbs,graphh,locations,R,a,b,lamb):
                           #grid_stack.append(ggs)
            
            grid_stack.remove(curr_grid)
            if grid_stack != []:
                curr_grid = grid_stack[0]
            else:
                break
            # Find a new curr_node
            for ggs in curr_component:
                if ggs in Z_Grid_nghbs[curr_grid]:
                    curr_node = Z_Grid[ggs][0]
                    
                    curr_c_size += 1
        # Need to start over with a new component
        if curr_c_size > max_c_size:
            max_c_size = curr_c_size
        curr_grid = None
        for ks in Z_Grid.keys():
            if ks in grids_unconsidered:
                grids_unconsidered.remove(ks)
                curr_grid = ks
                break
                #if Is_Good(ks,Z_Grid,Z_Grid_nghbs,graphh,locations,R,a,b,lamb):
                    #curr_grid = ks
                    #comp_count+=1
                    #break
        if curr_grid != None:
            curr_node = Z_Grid[curr_grid][0]
            curr_component = []
            grid_stack = [curr_grid] 
        else:
            break
    
      
    return estimates
        
   
   
def community_Detection2(graphh_PP,Z_Grid,Z_Grid_nghbs):
    num_nodes = max(graphh_PP.shape)
    estimates = np.empty((num_nodes,1))
    for i in range(num_nodes):
        estimates[i] = 1
    grids_unconsidered = Z_Grid.keys()
   
    curr_grid = None
    for ks in Z_Grid.keys():
        grids_unconsidered.remove(ks)

        if Is_Good_PP(ks,Z_Grid,Z_Grid_nghbs,graphh_PP):
            curr_grid = ks
            break
    if curr_grid == None:
        return estimates
 
    curr_node = Z_Grid[curr_grid][0]
    curr_component = [] # Finished Partitions
    grid_stack = [curr_grid]
    
    max_c_size = 0;
    
    while grids_unconsidered != []:
        
        curr_c_size = 0
        while grid_stack != [] :
            # Produce a complete clustering of the nodes in the connectecd component of the current grid
            # Breadth First Search 
            # Partition points in curr_grid using curr_node
            for pps in Z_Grid[curr_grid]:
                if pps == curr_node:
                    continue
                else:
                    estimates[pps] = graphh_PP[pps][curr_node]*estimates[curr_node]
            #curr_component.append(curr_grid)
            curr_component = curr_component + [curr_grid]
            #print grids_unconsidered
            
            for ggs in Z_Grid_nghbs[curr_grid]:
                
                if ggs in grids_unconsidered  :
                       grids_unconsidered.remove(ggs)
                       #grid_stack.append(ggs)
                       
                       if Is_Good_PP(ggs,Z_Grid,Z_Grid_nghbs,graphh_PP):
                           grid_stack.append(ggs)
            
            grid_stack.remove(curr_grid)
            if grid_stack != []:
                curr_grid = grid_stack[0]
            else:
                break
            # Find a new curr_node
            for ggs in curr_component:
                if ggs in Z_Grid_nghbs[curr_grid]:
                    curr_node = Z_Grid[ggs][0]
                    
                    curr_c_size += 1
        # Need to start over with a new component
        if curr_c_size > max_c_size:
            max_c_size = curr_c_size
        curr_grid = None
        for ks in Z_Grid.keys():
            if ks in grids_unconsidered:
                grids_unconsidered.remove(ks)
                #curr_grid = ks
                #break
                if Is_Good_PP(ks,Z_Grid,Z_Grid_nghbs,graphh_PP):
                    curr_grid = ks
                    #comp_count+=1
                    break
        if curr_grid != None:
            curr_node = Z_Grid[curr_grid][0]
            curr_component = []
            grid_stack = [curr_grid] 
        else:
            break
    
      
    return estimates         
    


def test_pairwise(Z_Grid,Z_Grid_nghbs,graphh,locations,colors,N,R,a,b,lamb):
    tc = 0
    cou = 0
    
    unconsidered = Z_Grid.keys()
    for zk in Z_Grid.keys():
        node_cons_list = []
        if zk in unconsidered:
            unconsidered.remove(zk)
            node_cons_list = Z_Grid[zk]
        
    
            for i in range(len(node_cons_list)):
                for j in range(i+1,len(node_cons_list)):
                    ii = node_cons_list[i]
                    jj = node_cons_list[j]
                    lltest = Pairwise_Classify(node_cons_list[i],node_cons_list[j],graphh,locations,N,R,a,b,lamb)
                    cou+=1
                    if lltest == colors[ii]*colors[jj]:
                        tc += 1
    return tc/cou
    

    
def create_PP(graphh,locations,N,R,a,b,lamb):
    num_nodes = max(locations.shape)
    graphh_PP = np.empty((num_nodes,num_nodes))
    for i in range(num_nodes):
        for j in range(i+1,num_nodes):
            if LA.norm(locations[i] - locations[j]) < 1.95*R :
                lval = Pairwise_Classify(i,j,graphh,locations,N,R,a,b,lamb)
                graphh_PP[i][j] = lval
                graphh_PP[j][i] = lval
            else:
                graphh_PP[i][j] = 0
                graphh_PP[j][i] = 0
    return graphh_PP



def return_distance(x,y,N,dims):
    aa = np.empty((dims,1))
    for i in range(dims):
        val = np.absolute(x[i] - y[i])
        if val > math.sqrt(N)/2 :
            aa[i] = math.sqrt(N) - val
        else:
            aa[i] = val
    return LA.norm(aa,2)
            

N = 1000;
R = 0.8;
lamb = 5;
a = 0.8;
b = 0.2;
dims = 2;
nl = math.pow(N,1/dims)/2;

num_nodes = np.random.poisson(lamb*N);

locations = np.empty((num_nodes,dims));
colors = np.empty((num_nodes,1));
graphh = np.empty((num_nodes,num_nodes));
graphh_PP = np.empty((num_nodes,num_nodes));

for x in range(0,num_nodes):
    colors[x] = 2*(random.random() < 0.5) - 1;
    for xx in range(0,dims):
       locations[x][xx] = np.random.uniform(-nl,nl);
       


        
for x in range(0,num_nodes):
    graphh[x][x] = 0
    for y in range(x+1,num_nodes):
        #distxy = LA.norm(locations[x] - locations[y],2);
        distxy = return_distance(locations[x],locations[y],N,dims)
        edgeran = random.random();
        if ((colors[x] == colors[y]) and (edgeran < a) and (distxy < R) ):
            graphh[x][y] = 1;
            graphh[y][x] = 1;
        elif((colors[x] != colors[y]) and (edgeran < b) and (distxy < R)):
             graphh[x][y] = 1;
             graphh[y][x] = 1;
        else:
            graphh[x][y] = 0;
            graphh[y][x] = 0;
            
            
#graphh_PP = create_PP(locations,R,a,b,lamb)
#num_grids_line = 4*N/R ;
Z_Grid = create_Grid_String(dims,N,R,locations)
Z_Grid_nghbs = create_neighbors_Z_Grid(Z_Grid,dims)
print test_pairwise(Z_Grid,Z_Grid_nghbs,graphh,locations,colors,N,R,a,b,lamb)




#evals = LA.eigvalsh(graphh)
#plt.plot(evals)


# Creating a Dummy Graph
#davg = sum(sum(graphh))/(2*num_nodes);
#
#graphh_dummy = np.empty((num_nodes,num_nodes))
#
#for x in range(num_nodes):
#    graphh_dummy[x][x] = 0
#    for y in range(x+1,num_nodes):
#        edgeran = random.random();
#        if edgeran > (2*davg)/num_nodes :
#            graphh_dummy[x][y] = 0
#            graphh_dummy[y][x] = 0
#            continue
#        else:
#            graphh_dummy[x][y] = 1
#            graphh_dummy[y][x] = 1
#
#evals_d = LA.eigvalsh(graphh_dummy)
#plt.plot(evals_d,'g')

#rint np.dot(colors.reshape(num_nodes,),est2.reshape(num_nodes,))
#print estimates