# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 17:22:36 2017

@author: abishek
"""

from __future__ import division
import numpy as np
import math
import random
from numpy import linalg as LA
import collections



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
            print ("Trouble")

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
    
def create_two_dim_nghbs(z_key):
    nbr_list = []
    z_key_vals_list = z_key.split()
    nbr_list.append(str(int(z_key_vals_list[0]) + 1) + " "+ str(int(z_key_vals_list[1]) + 1) + " ")
    nbr_list.append(str(int(z_key_vals_list[0]) - 1) + " "+ str(int(z_key_vals_list[1]) + 1) + " ")
    nbr_list.append(str(int(z_key_vals_list[0]) - 1) + " "+ str(int(z_key_vals_list[1]) - 1) + " ")
    nbr_list.append(str(int(z_key_vals_list[0]) + 1) + " "+ str(int(z_key_vals_list[1]) - 1) + " ")
    return nbr_list
  

##########################################################################
# CORE ALGORITHM 
##########################################################################

def form_Local_Partition(z_key_list,Z_Grid,Z_Grid_nghbs,graph):
    num_nodes = max(graph.shape)
    local_estimates = np.empty((num_nodes,1))
    for i in range(num_nodes):
        local_estimates[i] = 0
       
    # This will be either a naive spectral method or an improved SDP
       
    node_list = []
    for i in z_key_list:
        if i in Z_Grid.keys():
            node_list = node_list + Z_Grid[i]

    if len(node_list) == 1:
        local_estimates[node_list[0]] = 1
        return local_estimates
    else:
        small_graph = np.empty((len(node_list),len(node_list)))
        for i in range(len(node_list)):
            small_graph[i][i] = 0
        
        for i in range(len(node_list)):
            for j in range(i+1,len(node_list)):
                ii = node_list[i]
                jj = node_list[j]
                small_graph[i][j] = graph[ii][jj]
                small_graph[j][i] = graph[jj][ii]
        
        # NAIVE SPECTRAL ALGORITHM. Can be replaced with SDP if desired
        
        sm_graph_avg = sum(sum(small_graph))/(2*len(node_list))
        [w,v] = np.linalg.eig(small_graph - sm_graph_avg*np.ones((len(node_list),len(node_list)))  )
        idx = np.abs(w).argsort()[::-1]
        est_small = np.sign(v[:,idx[1]])
        for i in range(len(node_list)):
            local_estimates[node_list[i]] = est_small[i] 
        
        # Do a majority clean-up of the local solution
        
        local_estimates_old = local_estimates
        local_estimates_new = np.zeros((num_nodes,1))
        num_total_iters = 0
        for iters in range(num_total_iters):
            lsum = 0
            for i in range(len(node_list)):
                lsum = 0
                for j in range(len(node_list)):
                    if small_graph[i][j] == 1:
                        lsum += local_estimates_old[node_list[j]]
                local_estimates_new[node_list[i]] = np.sign(lsum)
                if lsum == 0:
                    local_estimates_new[node_list[i]] = local_estimates_old[node_list[i]]
            local_estimates_old = local_estimates_new
            
#        local_estimates_new = np.empty((num_nodes,1))
#        for i in range(num_nodes):
#            if i not in node_list:
#                local_estimates_new[i] = 0
#        local_estimates_old = local_estimates
#        for num_updates in range(5):
#            
#            for k in range(len(node_list)):
#                lsum = 0
#                
#                for kk in range(len(node_list)):
#                    if small_graph[k][kk] == 1:
#                        lsum += local_estimates_old[node_list[kk]]
#                    local_estimates_new[node_list[k]] = np.sign(lsum)
#            local_estimates_old = local_estimates_new
        return local_estimates_old

# We can update the above naive algorithm with sophisticated ones like SDP    

# The Key algorithm for getting a rough clustering   
def get_rough_estimates(graph,Z_Grid,Z_Grid_nghbs):
    num_nodes = max(graph.shape)
    num_shuffles = 10;
    estimates = np.empty((num_nodes,num_shuffles))
    estimates_majority = np.empty((num_nodes,num_shuffles))

    for i in range(num_nodes):
        for j in range(num_shuffles):
            estimates[i][j] = 0
            
    for shuffles in range(num_shuffles):
        # Synchronized Majorirty
        node_estimate_list = {}
        for i in range(num_nodes):
            node_estimate_list[i] = []
        Z_keys_shuffle = list(Z_Grid.keys())
        random.shuffle(Z_keys_shuffle)
        grids_estimated = []
        grids_not_estimated = Z_keys_shuffle
        curr_grid = grids_not_estimated[0]
        grid_stack = [curr_grid]
        grids_not_estimated.remove(curr_grid)
            
        while len(grids_not_estimated) > 0 or len(grid_stack) > 0:
            # Curr Grid's 8 neighbors needs to be clustered
            total_list = create_two_dim_nghbs(curr_grid) + Z_Grid_nghbs[curr_grid] + [curr_grid]
            local_estimates = form_Local_Partition(total_list,Z_Grid,Z_Grid_nghbs,graph)
            
            # Test whether the local estimate is aligned
            tsum = sum(estimates[:,shuffles]*local_estimates[:,0])
            if tsum < 0:
                local_estimates = -1*local_estimates
            for i in range(num_nodes):
                if local_estimates[i] !=0:
                    node_estimate_list[i] = node_estimate_list[i] + [local_estimates[i]]                
                
                if estimates[i][shuffles] == 0 and local_estimates[i] != 0:
                    estimates[i][shuffles] = local_estimates[i]
            
            # Update the stack and curr_grid
               
            for i in Z_Grid_nghbs[curr_grid]:
                if i  in grids_not_estimated  :
                    grid_stack = grid_stack + [i]
                    grids_not_estimated.remove(i)

            grid_stack.remove(curr_grid)
            grids_estimated = grids_estimated + [curr_grid]

            if grid_stack != []:
                curr_grid = grid_stack[0]
            elif len(grids_not_estimated) > 0 :
                curr_grid = grids_not_estimated[0]
                grid_stack = [curr_grid]
                grids_not_estimated.remove(curr_grid)
            else:
                break
      
        for ii in range(num_nodes):
            if estimates[ii][shuffles] == 0:
                estimates[ii][shuffles] = 1
            estimates_majority[ii][shuffles] = np.sign(sum(node_estimate_list[ii]))
            if estimates_majority[ii][shuffles] == 0:
                estimates_majority[ii][shuffles] = 1

    return estimates_majority
            
def global_cleanup(graph,estimates,locations,N,R,a,b,lamb,dims):
    num_nodes = max(graph.shape)
    # Run a single (or few) hill-climbing steps to produce an overall clustering
    est_new = np.empty((num_nodes,1))
    est_old = estimates
    num_iters = max(15*int(math.log(N)),5)
    #num_iters = 1
    neighbor_list = graph_neighbors(graph,locations,N,R,dims)
    for iters in range(num_iters):
        
        for i in range(num_nodes):
            N_plus = 0 
            N_minus = 0 
            N_plus_bar = 0
            N_minus_bar = 0
            local_n_list = neighbor_list[i]
            for jj in range(len(local_n_list)):
             
                j = local_n_list[jj]
                if graph[i][j] == 1 and est_old[j] == 1:
                    N_plus += 1
                elif graph[i][j] == 1 and est_old[j] == -1:
                    N_minus += 1
                elif graph[i][j] == 0 and est_old[j] == 1:
                    N_plus_bar +=1
                elif graph[i][j] == 0 and est_old[j] == -1:
                    N_minus_bar += 1
            
            
            dec_var = (N_plus - N_minus)*math.log(a/b) + (N_plus_bar - N_minus_bar)*math.log((1-a)/(1-b))
            if dec_var > 0:
                est_new[i] = 1
            else:
                est_new[i] = -1
        est_old = est_new
                
    
    return est_new
   
# The main control loop to Perform the Clustering in the above steps

def return_clustering(graph,locations,N,R,a,b,lamb,dims,colors):
    # The main routine to produce an estimate
    num_nodes = max(graph.shape)
    dims = locations.shape[1]            
    # Now we will perform the graph clustering based on our algorithm
                
    Z_Grid = create_Grid_String(dims,N,R,locations)
    Z_Grid_nghbs = create_neighbors_Z_Grid(Z_Grid,dims)
    
    rough_estimates = get_rough_estimates(graph,Z_Grid,Z_Grid_nghbs)
    num_shuffles = rough_estimates.shape[1]
    #print ("Rough Estimates")
    #for uy in range(num_shuffles):
        #print (np.abs(sum(colors[:,0]*rough_estimates[:,uy])/num_nodes))
 
    print ("   ")
    print ("Intermediate Clustering Performance by choosing random different Connected Components")
    #better_estimates = global_cleanup(graph,rough_estimates,locations,N,R,a,b,lamb)
    #final_estimate_collections = rough_estimates
    
    final_estimate_collections = np.empty((num_nodes,num_shuffles))
    for uy in range(num_shuffles):
        final_estimate_collections[:,uy] = global_cleanup(graph,rough_estimates[:,uy],locations,N,R,a,b,lamb,dims)[:,0]
 

    rsum = 0
    for j in range(num_shuffles):
        print (np.abs(sum(colors[:,0]*final_estimate_collections[:,j])/num_nodes))
        rsum += np.abs(sum(colors[:,0]*final_estimate_collections[:,j])/num_nodes)
    
    #print ("Answers and Benchmarks From Selection")
    print ("   ")
    print ("The FINAL TRUE OUTPUT")
    
    # Taking a Majority Vote 
    corrected_estimates = np.empty((num_nodes,num_shuffles))
    for i in range(num_shuffles):
        tsum = 0
        for j in range(num_shuffles):
            if j == i:
                continue
            tsum += sum(final_estimate_collections[:,i]*final_estimate_collections[:,j])
        if tsum < 0:
            corrected_estimates[:,i] = -1*final_estimate_collections[:,i]
        else:
            corrected_estimates[:,i] = final_estimate_collections[:,i]
            
    estt_cum = np.empty((num_nodes,1))
    for i in range(num_nodes):
        temp_sum = sum(corrected_estimates[i,:])
        if temp_sum < 0:
            estt_cum[i] = -1
        else:
            estt_cum[i] = 1
    estt_cum_corrected = global_cleanup(graph,estt_cum,locations,N,R,a,b,lamb,dims) 
    return estt_cum_corrected[:,0]

def evaluate_overlap(estimates,colors,graph,locations,N,R,a,b,lamb,dims): 
    num_nodes = colors.shape[0]     
    final_algo_performance = np.abs(sum(colors[:,0]*estimates)/num_nodes)
    print (final_algo_performance)
    
#    # Sanity Check

    [w,v] = np.linalg.eig(graph - ((sum(sum(graph)))/(2*num_nodes))*np.ones((num_nodes,num_nodes)))
    idx = np.abs(w).argsort()[::-1]

    #print ("Spectral Algorithm ignoring spatial labels PLUS GLOBAL")
    sanity = np.sign(v[:,idx[1]])
    updated_sanity = global_cleanup(graph,np.sign(v[:,idx[1]]),locations,N,R,a,b,lamb,dims) 
    sanity_perf = np.abs(sum(updated_sanity[:,0]*colors[:,0] )/num_nodes)
    #print (sanity_perf)
    
     

# Returns the threshold for Exact-Recovery regime
    
def what_is_lambda(a,b):
    return math.pow(math.pi*(1 - math.sqrt(a*b) - math.sqrt((1-a)*(1-b))),-1)


def generate_graph(N,R,a,b,lamb,dims):
    dims = 2
    nl = math.pow(N,1/dims)/2;
    num_nodes = np.random.poisson(lamb*N);
    
    locations = np.empty((num_nodes,dims))
    graph = np.empty((num_nodes,num_nodes))
    colors = np.empty((num_nodes,1))
    for x in range(num_nodes):
        graph[x][x] = 0
        colors[x] = 2*(random.random() < 0.5) - 1
        for xx in range(dims):
            locations[x][xx] = np.random.uniform(-nl,nl)
    
    for x in range(num_nodes):
        for y in range(x+1,num_nodes):
            distxy = return_distance(locations[x],locations[y],N,dims)
            edgeran = random.random()
            if ((colors[x] == colors[y]) and (edgeran < a) and (distxy < R) ):
                graph[x][y] = 1;
                graph[y][x] = 1;
            elif((colors[x] != colors[y]) and (edgeran < b) and (distxy < R)):
                 graph[x][y] = 1;
                 graph[y][x] = 1;
            else:
                graph[x][y] = 0;
                graph[y][x] = 0;

    return [graph, locations, N, R, a, b, lamb, colors]
    
def graph_neighbors(graph,locations,N,R,dims):
    num_nodes = max(graph.shape)
    neighbor_list = {}
    
    for i in range(num_nodes):
        neighbor_list[i] = []
    for i in range(num_nodes):
        for j in range(i+1,num_nodes):
            
            if return_distance(locations[i],locations[j],N,dims) < R:
                neighbor_list[i] = neighbor_list[i] + [j]
                neighbor_list[j] = neighbor_list[j] + [i]
    return neighbor_list
    
def return_distance(x,y,N,dims):
    aa = np.empty((dims,1))
    for i in range(dims):
        val = np.absolute(x[i] - y[i])
        if val > math.sqrt(N)/2 :
            aa[i] = math.sqrt(N) - val
        else:
            aa[i] = val
    return LA.norm(aa,2)

def estimate_parameters_simulation_data(graph,locations,dims): # Testing the learning algorithm
    nl = 0 
    num_nodes = max(graph.shape)
    for dd in range(dims):
        lmax = max(abs(locations[:,dd]))
        if lmax > nl:
            nl = lmax
    N_hat = math.pow(2*nl,dims)
    lamb_hat = max(graph.shape)/N_hat
    d_avg = sum(sum(graph))/(num_nodes)
    delt = 0
    R_hat = 0
        
    R_hat = 0
    for i in range(num_nodes):
        for j in range(i+1,num_nodes):
            l_dist = return_distance(locations[i,:],locations[j,:],N_hat,dims)
            if graph[i][j] == 1 and l_dist > R_hat:
                R_hat = l_dist
            for k in range(j+1,num_nodes):
                if graph[i][j] == 1 and graph[i][k] == 1 and graph[k][j] == 1:
                    delt += 3
    delt = delt/num_nodes

    d_avg_hat = d_avg/(lamb_hat*math.pi*R_hat*R_hat)
    C_R = 18.5 # Dimension dependent which we will insert from an hard-coded estimation
    l_stat = max((8*delt/(C_R*math.pow(R_hat*lamb_hat,dims))) - math.pow(2*d_avg_hat,3),0)
    a_hat = min(d_avg_hat + 0.5*math.pow(l_stat,0.333) , 1)
    b_hat = max(d_avg_hat - 0.5*math.pow(l_stat,0.333) , 0)
    return [N_hat, R_hat, a_hat, b_hat, lamb_hat]



def estimate_parameters_real_data(graph_list,locations,dims):
 
    num_nodes = len(locations.keys())
    location_values = locations.values()
    xaxis_list = [item[0] for item in location_values]
    yaxis_list = [item[1] for item in location_values]
    x_offset = 0.5*(max(xaxis_list) + min(xaxis_list))
    y_offset = 0.5*(max(yaxis_list) + min(yaxis_list))
    scale_factor = (max(yaxis_list) - y_offset)/(max(xaxis_list) - x_offset)
    locations_scaled = {}
    
    fnew = open("User_Locations_Scaled","w")
    for users in locations.keys():
        newx = (locations[users][0] - x_offset)*scale_factor
        newy = locations[users][1] - y_offset
        locations_scaled[int(users)] = [float(newx),float(newy)]
        fnew.write(str(users) + " " + str(newx) + " " + str(newy) + "\n")
    fnew.close()
    locations_scaled.keys().sort()
    
    # Computing Average Degree
    dsum = 0
    nsum = 0
    for users in graph_list:
        dsum += len(graph_list[users])
        nsum +=1
    d_avg = dsum/num_nodes
    N_hat = math.pow(max(yaxis_list) - y_offset,2)
    lamb_hat = num_nodes/N_hat

    num_triangles_data = 494728 # Total Number of Triangles
    delt = num_triangles_data/num_nodes
    R_hat = 0
    list_of_edge_lengths = []
    
    
    for i in graph_list.keys():
        
        if i not in graph_list.keys() or i not in locations_scaled.keys():
            continue
        for jj in range(len(graph_list[i])):
            j = graph_list[i][jj]
            if j < i:
                continue
      
            if j not in graph_list.keys() or j not in locations_scaled.keys():
                continue
            l_dist = LA.norm(np.array(locations_scaled[i])-np.array(locations_scaled[j]))
            list_of_edge_lengths = list_of_edge_lengths + [l_dist]            
            if l_dist > R_hat:
                R_hat = l_dist

    f_longest_edge = open("Longest_Edge_Stats","r")
    f_longest_edge.write("Longest Edge Length " + str(R_hat))
    f_longest_edge.close()
    

    d_avg_hat = d_avg/(lamb_hat*math.pi*R_hat*R_hat)
    C_R = 18.5 # Dimension dependent which we will insert from an hard-coded estimation
    l_stat = max((8*delt/(C_R*math.pow(R_hat*lamb_hat,dims))) - math.pow(2*d_avg_hat,3),0)
    a_hat = min(d_avg_hat + 0.5*math.pow(l_stat,0.333) , 1)
    b_hat = max(d_avg_hat - 0.5*math.pow(l_stat,0.333) , 0)
    return [N_hat, R_hat, a_hat, b_hat, lamb_hat]
    
    
def estimate_from_clean_data(locations_scaled,graph_list,dims=2):
    locations_scaled.keys().sort()
    
    # Computing N_hat
    xaxis_list = []
    yaxis_list = []
    for locs in locations_scaled.keys():
        xaxis_list = xaxis_list + [locations_scaled[locs][0]]
        yaxis_list = yaxis_list + [locations_scaled[locs][1]]
    
    
    num_nodes = min(len(locations_scaled.keys()), len(graph_list.keys()))
    # Computing Average Degree
    dsum = 0
    nsum = 0
    for users in graph_list:
        dsum += len(graph_list[users])
        nsum +=1
    d_avg = dsum/num_nodes
    N_hat = math.pow(max(yaxis_list),2)
    lamb_hat = num_nodes/N_hat

    num_triangles_data = 494728 # Total Number of Triangles
    delt = num_triangles_data/num_nodes
    R_hat = 0
    list_of_edge_lengths = []
    
    
    for i in graph_list.keys():
        
        if i not in graph_list.keys() or i not in locations_scaled.keys():
            continue
        for jj in range(len(graph_list[i])):
            j = graph_list[i][jj]
            if j < i:
                continue
      
            if j not in graph_list.keys() or j not in locations_scaled.keys():
                continue
            l_dist = LA.norm(np.array(locations_scaled[i])-np.array(locations_scaled[j]))
            list_of_edge_lengths = list_of_edge_lengths + [l_dist]            
            if l_dist > R_hat:
                R_hat = l_dist

    f_longest_edge = open("Longest_Edge_Stats","r")
    f_longest_edge.write("Longest Edge Length " + str(R_hat))
    f_longest_edge.close()
    

    d_avg_hat = d_avg/(lamb_hat*math.pi*R_hat*R_hat)
    C_R = 18.5 # Dimension dependent which we will insert from an hard-coded estimation
    l_stat = max((8*delt/(C_R*math.pow(R_hat*lamb_hat,dims))) - math.pow(2*d_avg_hat,3),0)
    a_hat = min(d_avg_hat + 0.5*math.pow(l_stat,0.333) , 1)
    b_hat = max(d_avg_hat - 0.5*math.pow(l_stat,0.333) , 0)
    return [N_hat, R_hat, a_hat, b_hat, lamb_hat]
    
    
def visualize_cleaned_real_data():
    f_locs = open("User_Locations_Scaled","r")
    f_edges = open("loc-brightkite_edges","r")
    locations_scaled = {}
    for lines in f_locs:
        lstring = lines.split(" ")
        locations_scaled[int(lstring[0])] = [float(lstring[1]),float(lstring[2])]
    graph_list = {}
    for lines in f_edges:
        lstring = lines.split('\t')
        from_v = int(lstring[0])
        to_v = int(lstring[1])
        if from_v in graph_list.keys():
            graph_list[from_v] = graph_list[from_v]+[int(to_v)]
        else:
            graph_list[from_v] = [int(to_v)]
    f_locs.close()
    f_edges.close()
    estimate_from_clean_data(locations_scaled,graph_list)
                

def parse_brightkite_rawdata():
    f_locs = open("User_Locations","r")
    f_edges = open("loc-brightkite_edges","r")
    #num_nodes = 58228
    graph_list =  {}
    for lines in f_edges:
        lstring = lines.split('\t')
        from_v = int(lstring[0])
        to_v = int(lstring[1])
        if from_v in graph_list:
            graph_list[from_v] = graph_list[from_v] + [int(to_v)]
        else:
            graph_list[from_v] = [int(to_v)]
    
    locations = {}
    for lines in f_locs:
        lstring = lines.split(" ")
        locations[int(lstring[0])] = [float(lstring[1]),float(lstring[2])]
        
    [N_hat, R_hat, a_hat, b_hat, lamb_hat] = estimate_parameters_real_data(graph_list,locations,2)
    f_locs.close()
    f_edges.close()
    f_parameters = open("BrightKite_Fitted_Parameters","w")
    f_parameters.write("N_hat " + str(N_hat) + "\n")
    f_parameters.write("R_hat " + str(R_hat)+"\n")
    f_parameters.write("a_hat " + str(a_hat) + "\n")
    f_parameters.write("b_hat " + str(b_hat) + "\n")
    f_parameters.write("lamb_hat " + str(lamb_hat) + "\n")
    
    f_parameters.close()

def main():
    N = 50
    R = math.sqrt(math.log(N))
    a=0.7
    b=0.5
    dims = 2
    lamb = what_is_lambda(a,b) 
    #lamb = 2.3
    [graph,locations, N, R, a, b, lamb, colors] = generate_graph(N,R,a,b,lamb,2)
    # Real Data Plug - Use the function to estimate .

    # Read and estimate parameters using    
    
    estimates = return_clustering(graph,locations,N,R,a,b,lamb,dims,colors)
    evaluate_overlap(estimates,colors,graph,locations,N,R,a,b,lamb,dims)
    
    
   
if __name__ == '__main__':
    main()
    
    
