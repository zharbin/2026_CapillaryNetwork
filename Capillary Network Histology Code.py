# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:16:10 2025

@author: Zach Harbin
"""


import matplotlib
#matplotlib.use('TkAgg')  # Set backend to enable an interactive window
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time

import matplotlib as mpl
import matplotlib.patheffects as pe


mpl.rcParams['font.family'] = 'Arial'

start_time = time.time()  # Start the timer

# Parameters for the Capillary Network
Hist_SA_Per = 1.4082
Hist_SA_Per_Bounds = np.array([Hist_SA_Per*0.99, Hist_SA_Per*1.01])


timepoint = 4


cube_length = 500 # length of the cube in um


def random_branch_length():
    # Define bins and counts
    # For Cube
    bin_ranges = [(5, 10), (10, 50), (50, 100), (100, 300)]
    counts = [3218, 5540, 1701, 1310]
    
    # For DePietro Validation
    # bin_ranges = [(5, 10), (10, 50), (50, 100), (100, 500), (500, 1000)]
    # counts = [3218, 5540, 1701, 1264, 46]
    
    # Step 1: Normalize the counts to create a PDF
    total_counts = sum(counts)
    pdf = [c / total_counts for c in counts]  # Normalized probabilities

    # Choose a bin based on the PDF
    selected_bin_index = np.random.choice(len(bin_ranges), p=pdf)
    selected_bin = bin_ranges[selected_bin_index]
    
    # Generate a random length within the selected bin
    sampled_length = np.random.uniform(selected_bin[0], selected_bin[1])
    return sampled_length

def random_branch_order():
    # Define bins and counts
    bin_ranges = [1, 2, 3, 4, 5]
    counts = [10, 46, 32, 9, 3]
    
    # Step 1: Normalize the counts to create a PDF
    total_counts = sum(counts)
    pdf = [c / total_counts for c in counts]  # Normalized probabilities

    # Choose a bin based on the PDF
    selected_order_index = np.random.choice(len(bin_ranges), p=pdf)
    selected_order = bin_ranges[selected_order_index]
    
    return selected_order

def generate_cylinder_points(start, end, diameter, thickness, num_points):
    """
    Generate random points inside a hollow cylinder.

    Parameters:
    - start, end: (x, y, z) coordinates of the cylinder centerline endpoints
    - diameter: Outer diameter of the cylinder
    - thickness: Thickness of the cylinder wall
    - num_points: Number of random points to generate

    Returns:
    - points: (N, 3) array of (x, y, z) points inside the cylinder
    """
    # Compute cylinder properties
    r_outer = diameter / 2
    r_inner = r_outer - thickness
    length_vector = np.array(end) - np.array(start)
    cylinder_length = np.linalg.norm(length_vector)
    unit_vector = length_vector / cylinder_length
    
    # Generate random points along the cylinder axis
    t_vals = np.random.uniform(0, cylinder_length, num_points)
    random_offsets = t_vals[:, None] * unit_vector  # Move along the centerline
    points_on_axis = np.array(start) + random_offsets

    # Generate random radial points
    radii = np.sqrt(np.random.uniform(r_inner**2, r_outer**2, num_points))
    angles = np.random.uniform(0, 2 * np.pi, num_points)
    x_offsets = radii * np.cos(angles)
    y_offsets = radii * np.sin(angles)

    # Create orthonormal basis for cylinder cross-section
    v = np.array([1, 0, 0]) if np.abs(unit_vector[2]) > 0.99 else np.array([0, 0, 1])
    perp1 = np.cross(unit_vector, v)
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(unit_vector, perp1)

    # Offset points from centerline
    points = points_on_axis + x_offsets[:, None] * perp1 + y_offsets[:, None] * perp2
    
    return points

def Capillary_Thickness_Gen(week):
    # Define bins and counts
    # For Cube
    if week == 0:
        counts = [0, 6, 6, 6, 13, 5, 6, 2, 3]
    
    if week == 1:
        counts = [0, 0, 3, 2, 12, 20, 7, 4, 4]
    
    if week == 4:
        counts = [0, 6, 9, 29, 11, 5, 3, 1, 0]
        
    if week == 16:
        counts = [0, 2, 7, 15, 15, 15, 5, 1, 0]
        
    bin_ranges = [(0.5, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 3.5), (3.5, 4), (4, 4.5), (4.5, 5)]
    

    # Step 1: Normalize the counts to create a PDF
    total_counts = sum(counts)
    pdf = [c / total_counts for c in counts]  # Normalized probabilities

    # Choose a bin based on the PDF
    selected_bin_index = np.random.choice(len(bin_ranges), p=pdf)
    selected_bin = bin_ranges[selected_bin_index]
    
    # Generate a random length within the selected bin
    sampled_thickness = np.random.uniform(selected_bin[0], selected_bin[1])
    return sampled_thickness


def Capillary_Feret_Diameter_Gen(week):
    # Define bins and counts
    # For Cube
    if week == 0:
        counts = [0, 15, 13, 10, 4, 2, 0, 0, 2, 1]
        
    if week == 1:
        counts = [0, 0, 11, 15, 12, 6, 5, 1, 2, 1]
    
    if week == 4:
        counts = [0, 6, 25, 14, 11, 3, 2, 0, 1, 3]
        
    if week == 16:
        counts = [0, 8, 17, 11, 8, 13, 1, 1, 1, 0]
        
    bin_ranges = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 28), (28, 32), (32, 36), (36, 40)]
    

    # Step 1: Normalize the counts to create a PDF
    total_counts = sum(counts)
    pdf = [c / total_counts for c in counts]  # Normalized probabilities

    # Choose a bin based on the PDF
    selected_bin_index = np.random.choice(len(bin_ranges), p=pdf)
    selected_bin = bin_ranges[selected_bin_index]
    
    # Generate a random length within the selected bin
    sampled_diameter = np.random.uniform(selected_bin[0], selected_bin[1])
    return sampled_diameter



Total_CD31_slice = []
Total_Capillaries_in_Slice = []
Total_Capillary_Density = []

for i in range(500):
    print(i)
    if len(Total_CD31_slice) == 50:
        break
    # Initialize the NetworkX graph
    G = nx.Graph()
    
    cur_branch_num = 0
    node_id = 0  # Unique identifier for each node
    loop_num = 0
    num_loops = 11660
    CD31_slice_avg = 0 
    
    # Keep track of all nodes to add branches
    capillary_diameter_list = []
    capillary_thickness_list = []
    
    All_branch_entry_point = []
    All_branch_exit_point = []
    
    branches_in_slice = 0
    Histology_Slice_Points = []
    total_volume_in_slice = 0
    
    while abs(CD31_slice_avg - Hist_SA_Per) / Hist_SA_Per > 0.01:
        if CD31_slice_avg > 1.01 * Hist_SA_Per:
            print('SKIP - Restarting for loop')
            break  # Exit the while loop and restart the for loop iteration
            
        branch_order = random_branch_order()
        num_branches_added = 0
        loop_num += 1
        if branch_order >= 1:
            x = np.random.uniform(0, cube_length)
            y = np.random.uniform(0, cube_length)
            z = np.random.uniform(0, cube_length)
            starting_node_position = (x, y, z)
            
            # Add the starting node to the graph
            G.add_node(node_id, pos=starting_node_position)
            starting_node_id = node_id
            node_id += 1
            
                # Generate a new position
            while True:
                # Random direction in spherical coordinates
                theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
                phi = np.random.uniform(0, np.pi)        # Polar angle
                branch_length = random_branch_length()
        
                # Convert spherical to Cartesian coordinates
                dx = branch_length * np.sin(phi) * np.cos(theta)
                dy = branch_length * np.sin(phi) * np.sin(theta)
                dz = branch_length * np.cos(phi)
        
                # Calculate new position
                new_position = (starting_node_position[0] + dx, starting_node_position[1] + dy, starting_node_position[2] + dz)
        
                # Check if the new position is within the cylinder
                if 0 < new_position[0] < cube_length and 0 < new_position[1] < cube_length and 0 < new_position[2] < cube_length:
                    break  # Exit the loop when a valid position is found
                
            # Add the new node to the graph
            G.add_node(node_id, pos=new_position)
            new_node_id = node_id
            node_id += 1
    
            # Add an edge between the starting node and the new node
            G.add_edge(starting_node_id, new_node_id)
    
            # Increment the branch count
            cur_branch_num += 1
            num_branches_added += 1
            capillary_diameter_list.append(Capillary_Feret_Diameter_Gen(timepoint))
            capillary_thickness_list.append(Capillary_Thickness_Gen(timepoint))
            
            
        if branch_order >= 2:
            # Direction vector of the previous edge
            prev_edge_dir = np.array([
                new_position[0] - starting_node_position[0],
                new_position[1] - starting_node_position[1],
                new_position[2] - starting_node_position[2]
            ])
            prev_edge_dir /= np.linalg.norm(prev_edge_dir)  # Normalize the vector
        
            for _ in range(2):  # Create two new branches
                max_attempts = 100  # Prevent infinite loops
                attempts = 0
                while attempts < max_attempts:
                    # Random direction in spherical coordinates
                    theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
                    phi = np.random.uniform(0, np.pi)        # Polar angle
                    branch_length = random_branch_length()
            
                    # Convert spherical to Cartesian coordinates
                    dx = branch_length * np.sin(phi) * np.cos(theta)
                    dy = branch_length * np.sin(phi) * np.sin(theta)
                    dz = branch_length * np.cos(phi)
                    
                    # Calculate the new position
                    # Calculate new position
                    sec_order_branch_pos = (new_position[0] + dx, new_position[1] + dy, new_position[2] + dz)
            
                    new_edge_dir = np.array([
                        sec_order_branch_pos[0] - new_position[0],
                        sec_order_branch_pos[1] - new_position[1],
                        sec_order_branch_pos[2] - new_position[2]
                    ])
                    new_edge_dir /= np.linalg.norm(new_edge_dir)  # Normalize the vector
                    
                    
                    # Ensure the new direction is on the opposite side of the plane
                    if np.sign(prev_edge_dir[0]) == np.sign(new_edge_dir[0]) and np.sign(prev_edge_dir[1]) == np.sign(new_edge_dir[1]) and np.sign(prev_edge_dir[2]) == np.sign(new_edge_dir[2]) and 0 < sec_order_branch_pos[0] < cube_length and 0 < sec_order_branch_pos[1] < cube_length and 0 < sec_order_branch_pos[2] < cube_length:  # Valid direction (opposite side of the plane)
                        break  # Valid direction found
            
                    attempts += 1
            
                if attempts == max_attempts:
                    #print("Warning: Could not find a valid branch direction")
                    continue  # Skip this branch
            
                # Add the new node and edge
                G.add_node(node_id, pos=sec_order_branch_pos)
                G.add_edge(new_node_id, node_id)  # Connect to the previous node
                node_id += 1
                cur_branch_num += 1
                num_branches_added += 1
                capillary_diameter_list.append(Capillary_Feret_Diameter_Gen(timepoint))
                capillary_thickness_list.append(Capillary_Thickness_Gen(timepoint))
                    
                
        if branch_order >= 3:
            # Get the last two nodes that were added
            nodes_in_graph = list(G.nodes())
            node1_id = nodes_in_graph[-2]  # Second last node
            node2_id = nodes_in_graph[-1]  # Last node added
        
            # Get the positions of the last two nodes
            node1_pos = G.nodes[node1_id]['pos']
            node2_pos = G.nodes[node2_id]['pos']
            
            # 50% chance to apply branches to both nodes
            apply_to_both = np.random.rand() < 0.5
            if apply_to_both:
                selected_nodes = [(node1_id, node1_pos), (node2_id, node2_pos)]
                third_order_branches = 4
            else:
                # Randomly select one of the two nodes
                selected_nodes = [np.random.choice([node1_id, node2_id])]
                selected_nodes = [(selected_nodes[0], G.nodes[selected_nodes[0]]['pos'])]
                third_order_branches = 2
                # Apply two branches to the selected node(s)
            for base_node_id, base_pos in selected_nodes:
                connected_edges = list(G.edges(base_node_id))
                _, connected_node_id = connected_edges[0]  # Unpack the single edge
    
                # Calculate the direction of the existing branch
                connected_node_pos = np.array(G.nodes[connected_node_id]['pos'])
                base_node_pos = np.array(base_pos)
            
                # Direction vector for the existing branch
                existing_branch_dir = base_node_pos - connected_node_pos
                existing_branch_dir /= np.linalg.norm(existing_branch_dir)  # Normalize
                
                for _ in range(2):  # Two branches per node
                    max_attempts = 100  # Limit attempts to prevent infinite loops
                    attempts = 0
                    while attempts < max_attempts:
                        # Random direction in spherical coordinates
                        theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
                        phi = np.random.uniform(0, np.pi)        # Polar angle
                        branch_length = random_branch_length()
        
                        # Convert spherical to Cartesian coordinates
                        dx = branch_length * np.sin(phi) * np.cos(theta)
                        dy = branch_length * np.sin(phi) * np.sin(theta)
                        dz = branch_length * np.cos(phi)
        
                        # Calculate the new branch position
                        new_branch_pos = (base_pos[0] + dx, base_pos[1] + dy, base_pos[2] + dz)
                        
                        new_edge_dir = np.array([
                            new_branch_pos[0] - base_pos[0],
                            new_branch_pos[1] - base_pos[1],
                            new_branch_pos[2] - base_pos[2]
                        ])
                        new_edge_dir /= np.linalg.norm(new_edge_dir)  # Normalize the vector
        
                        # Check if the new position is within the cylinder
                        if np.sign(existing_branch_dir[0]) == np.sign(new_edge_dir[0]) and np.sign(existing_branch_dir[1]) == np.sign(new_edge_dir[1]) and np.sign(existing_branch_dir[2]) == np.sign(new_edge_dir[2]) and 0 < new_branch_pos[0] < cube_length and 0 < new_branch_pos[1] < cube_length and 0 < new_branch_pos[2] < cube_length:
                            break  # Valid direction found
        
                        attempts += 1
                
                    if attempts == max_attempts:
                        #print("Warning: Could not find a valid branch direction")
                        continue  # Skip this branch
        
                    G.add_node(node_id, pos=new_branch_pos)
                    G.add_edge(base_node_id, node_id)  # Connect to the current base node
                    node_id += 1
                    cur_branch_num += 1
                    num_branches_added += 1
                    capillary_diameter_list.append(Capillary_Feret_Diameter_Gen(timepoint))
                    capillary_thickness_list.append(Capillary_Thickness_Gen(timepoint))
    
        if branch_order >= 4:
            # Get the nodes created during branch order 3
            nodes_in_graph = list(G.nodes())
            third_order_start = -third_order_branches
            third_order_nodes = nodes_in_graph[third_order_start:]  # Last added nodes
            
            # Ensure at least one node sprouts two branches
            node_with_two_branches = np.random.choice(third_order_nodes)
            sprouting_nodes_4th = [node_with_two_branches]
            for node in third_order_nodes:
                if node != node_with_two_branches and np.random.rand() < 0.5:
                    sprouting_nodes_4th.append(node)
                    
            num_4th_sprouting_nodes = len(sprouting_nodes_4th)
            # Apply branches
            for base_node_id in sprouting_nodes_4th:
                base_pos = np.array(G.nodes[base_node_id]['pos'])
                
                # Find the connected edge and calculate the existing branch direction
                connected_edges = list(G.edges(base_node_id))
                _, connected_node_id = connected_edges[0]  # Assumes one edge for branch orders
                connected_node_pos = np.array(G.nodes[connected_node_id]['pos'])
                
                # Direction vector for the existing branch
                existing_branch_dir = base_pos - connected_node_pos
                existing_branch_dir /= np.linalg.norm(existing_branch_dir)  # Normalize
                # Add new branches
                for _ in range(2):  # Two branches per sprouting node
                    max_attempts = 100
                    attempts = 0
                    while attempts < max_attempts:
                        # Random direction in spherical coordinates
                        theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
                        phi = np.random.uniform(0, np.pi)        # Polar angle
                        branch_length = random_branch_length()
                        
                        # Convert spherical to Cartesian coordinates
                        dx = branch_length * np.sin(phi) * np.cos(theta)
                        dy = branch_length * np.sin(phi) * np.sin(theta)
                        dz = branch_length * np.cos(phi)
                        
                        # Calculate the new branch position
                        new_branch_pos = (base_pos[0] + dx, base_pos[1] + dy, base_pos[2] + dz)
                        new_edge_dir = np.array([
                            new_branch_pos[0] - base_pos[0],
                            new_branch_pos[1] - base_pos[1],
                            new_branch_pos[2] - base_pos[2]
                        ])
                        new_edge_dir /= np.linalg.norm(new_edge_dir)  # Normalize
                        
                        # Check if the new position is within the cylinder
                        if np.sign(existing_branch_dir[0]) == np.sign(new_edge_dir[0]) and np.sign(existing_branch_dir[1]) == np.sign(new_edge_dir[1]) and np.sign(existing_branch_dir[2]) == np.sign(new_edge_dir[2]) and 0 < new_branch_pos[0] < cube_length and 0 < new_branch_pos[1] < cube_length and 0 < new_branch_pos[2] < cube_length:
                            break  # Valid direction found
                        
                        attempts += 1
                    
                    if attempts == max_attempts:
                        #print(f"Warning: Could not find a valid branch direction for node {base_node_id}")
                        continue  # Skip this branch
                    
                    G.add_node(node_id, pos=new_branch_pos)
                    G.add_edge(base_node_id, node_id)  # Connect to the current base node
                    node_id += 1
                    cur_branch_num += 1
                    num_branches_added += 1
                    capillary_diameter_list.append(Capillary_Feret_Diameter_Gen(timepoint))
                    capillary_thickness_list.append(Capillary_Thickness_Gen(timepoint))
    
        if branch_order >= 5:
            # Get the nodes created during branch order 4
            nodes_in_graph = list(G.nodes())
            fourth_order_start = -num_4th_sprouting_nodes
            fourth_order_nodes = nodes_in_graph[fourth_order_start:]  # Last added nodes
            
            # Ensure at least one node sprouts two branches
            node_with_two_branches = np.random.choice(fourth_order_nodes)
            sprouting_nodes_5th = [node_with_two_branches]
            for node in fourth_order_nodes:
                if node != node_with_two_branches and np.random.rand() < 0.5:
                    sprouting_nodes_5th.append(node)
            
            num_sprouting_nodes = len(sprouting_nodes_5th)
            # Apply branches
            for base_node_id in sprouting_nodes_5th:
                base_pos = np.array(G.nodes[base_node_id]['pos'])
                
                # Find the connected edge and calculate the existing branch direction
                connected_edges = list(G.edges(base_node_id))
                _, connected_node_id = connected_edges[0]  # Assumes one edge for branch orders
                connected_node_pos = np.array(G.nodes[connected_node_id]['pos'])
                
                # Direction vector for the existing branch
                existing_branch_dir = base_pos - connected_node_pos
                existing_branch_dir /= np.linalg.norm(existing_branch_dir)  # Normalize
                # Add new branches
                for _ in range(2):  # Two branches per sprouting node
                    max_attempts = 100
                    attempts = 0
                    while attempts < max_attempts:
                        # Random direction in spherical coordinates
                        theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
                        phi = np.random.uniform(0, np.pi)        # Polar angle
                        branch_length = random_branch_length()
                        
                        # Convert spherical to Cartesian coordinates
                        dx = branch_length * np.sin(phi) * np.cos(theta)
                        dy = branch_length * np.sin(phi) * np.sin(theta)
                        dz = branch_length * np.cos(phi)
                        
                        # Calculate the new branch position
                        new_branch_pos = (base_pos[0] + dx, base_pos[1] + dy, base_pos[2] + dz)
                        new_edge_dir = np.array([
                            new_branch_pos[0] - base_pos[0],
                            new_branch_pos[1] - base_pos[1],
                            new_branch_pos[2] - base_pos[2]
                        ])
                        new_edge_dir /= np.linalg.norm(new_edge_dir)  # Normalize
                        
                        # Check if the new position is within the cylinder
                        if np.sign(existing_branch_dir[0]) == np.sign(new_edge_dir[0]) and np.sign(existing_branch_dir[1]) == np.sign(new_edge_dir[1]) and np.sign(existing_branch_dir[2]) == np.sign(new_edge_dir[2]) and 0 < new_branch_pos[0] < cube_length and 0 < new_branch_pos[1] < cube_length and 0 < new_branch_pos[2] < cube_length:
                            break  # Valid direction found
                        
                        attempts += 1
                    
                    if attempts == max_attempts:
                        #print(f"Warning: Could not find a valid branch direction for node {base_node_id}")
                        continue  # Skip this branch
                    
                    G.add_node(node_id, pos=new_branch_pos)
                    G.add_edge(base_node_id, node_id)  # Connect to the current base node
                    node_id += 1
                    cur_branch_num += 1
                    num_branches_added += 1
                    capillary_diameter_list.append(Capillary_Feret_Diameter_Gen(timepoint))
                    capillary_thickness_list.append(Capillary_Thickness_Gen(timepoint))
                    
    
        
        
        new_edges = list(G.edges)[-num_branches_added:]  

        # Extract the corresponding diameter and thickness values
        capillary_diameter_subset = capillary_diameter_list[-num_branches_added:]
        capillary_thickness_subset = capillary_thickness_list[-num_branches_added:]
        
        
        All_Slices_CD31 = []
        All_Slices_Total_Branches = []
        #slice_locations = np.array([25,75,125,175,225,275,325,375,425,475])
        #slice_locations = np.array([25,50,75,100,125,150,175,200,225])
        slice_locations = np.array([250])
        for i in range(len(slice_locations)):
            # Defining the slice
            z_slice = slice_locations[i] #1000  # Center of the slice with respect to the z-axis
            thickness = 4  # Thickness of the slice (From Histology)
            z_min = z_slice - (thickness / 2) - (max(capillary_diameter_list)/2)
            z_max = z_slice + (thickness / 2) + (max(capillary_diameter_list)/2)
            
            z_min_hist = z_slice - (thickness / 2)
            z_max_hist = z_slice + (thickness / 2)
            
            # Count branches in the slice and calculate their lengths within the slice
            
            
            for k, (u, v) in enumerate(new_edges):
                pos_u = np.array(G.nodes[u]['pos'])
                pos_v = np.array(G.nodes[v]['pos'])
                z_u = pos_u[2]
                z_v = pos_v[2]
                Capillary_Diameter = capillary_diameter_subset[k]
                Capillary_Thickness = capillary_thickness_subset[k]
                
                # Check if the branch intersects the slice
                if (z_min <= z_u <= z_max) or (z_min <= z_v <= z_max) or (z_u < z_min and z_v > z_max) or (z_v < z_min and z_u > z_max):
                
                    # Calculate the intersection points (entry and exit)
                    if z_u < z_min:  # Interpolate at z_min
                        t_entry = (z_min - z_u) / (z_v - z_u)
                        entry_point = pos_u + t_entry * (pos_v - pos_u)
                    elif z_u > z_max:  # Interpolate at z_max
                        t_entry = (z_max - z_u) / (z_v - z_u)
                        entry_point = pos_u + t_entry * (pos_v - pos_u)
                    else:
                        entry_point = pos_u
            
                    if z_v < z_min:  # Interpolate at z_min
                        t_exit = (z_min - z_u) / (z_v - z_u)
                        exit_point = pos_u + t_exit * (pos_v - pos_u)
                    elif z_v > z_max:  # Interpolate at z_max
                        t_exit = (z_max - z_u) / (z_v - z_u)
                        exit_point = pos_u + t_exit * (pos_v - pos_u)
                    else:
                        exit_point = pos_v
                    
                    All_branch_entry_point.append(entry_point)
                    All_branch_exit_point.append(exit_point)
                    
                    cut_branch_length = np.linalg.norm(np.array(exit_point) - np.array(entry_point))
                    
                    if CD31_slice_avg < (Hist_SA_Per*0.90):
                        Cylinder_num_points = 50000
                    else:
                        Cylinder_num_points = 400000
                    
                    Capillary_Points = generate_cylinder_points(entry_point, exit_point, Capillary_Diameter, Capillary_Thickness, Cylinder_num_points)
                    points_in_slice = (
                        (Capillary_Points[:, 2] >= z_min_hist) & (Capillary_Points[:, 2] <= z_max_hist) &  # Z-condition
                        (Capillary_Points[:, 0] >= 0) & (Capillary_Points[:, 0] <= cube_length) &  # X-condition
                        (Capillary_Points[:, 1] >= 0) & (Capillary_Points[:, 1] <= cube_length)    # Y-condition
                    )
                    vol_frac = np.sum(points_in_slice) / Cylinder_num_points
                    if vol_frac != 0:
                        branches_in_slice += 1
                        Histology_Slice_Points.append(Capillary_Points[points_in_slice])
                    est_capillary_volume = vol_frac * (np.pi/4)*(Capillary_Diameter**2 - (Capillary_Diameter - (2*Capillary_Thickness))**2)*cut_branch_length 
                    total_volume_in_slice += est_capillary_volume
                    
    
            Est_CD31 = (total_volume_in_slice/(cube_length*cube_length*thickness))*100
            Est_Capillary_Den = (cur_branch_num/(cube_length/1000)**3)
            All_Slices_CD31.append(Est_CD31)
            All_Slices_Total_Branches.append(total_volume_in_slice)
            
        CD31_slice_avg = np.average(All_Slices_CD31)
        if CD31_slice_avg > (Hist_SA_Per*0.90):
            print('Swiched')
            
    if abs(CD31_slice_avg - Hist_SA_Per) / Hist_SA_Per <= 0.01:
        Total_CD31_slice.append(CD31_slice_avg)
        Total_Capillaries_in_Slice.append(branches_in_slice)
        Total_Capillary_Density.append(Est_Capillary_Den)
    


Average_Total_CD31_slice = np.average(Total_CD31_slice)
Average_Capillaries_in_Slice = np.average(Total_Capillaries_in_Slice)
Std_Capillaries_in_Slice = np.std(Total_Capillaries_in_Slice)
Average_Capillary_Density = np.average(Total_Capillary_Density)
Std_Capillary_Density = np.std(Total_Capillary_Density)

# Create a 3D plot
import matplotlib as mpl

# ===============================
# Force Arial Globally
# ===============================
mpl.rcParams['font.family'] = 'Arial'

# ===============================
# Create Square Figure
# ===============================
fig, ax = plt.subplots(figsize=(6, 6))  # square figure

# Plot slice points (all red)
for points_array in Histology_Slice_Points:
    ax.scatter(
        points_array[:, 0],
        points_array[:, 1],
        color='red',      # all red
        s=1
    )

# ===============================
# Axis Limits
# ===============================
ax.set_xlim(0, cube_length)
ax.set_ylim(0, cube_length)

# Force TRUE square aspect ratio
ax.set_aspect('equal', adjustable='box')

# ===============================
# Labels (Arial + Bold)
# ===============================
ax.set_xlabel('X (μm)', fontsize=18, fontweight='bold')
ax.set_ylabel('Y (μm)', fontsize=18, fontweight='bold')


# Title (Arial normal)
ax.set_title(
    f'Simulated Histology Slice \n %CD31+ : {Hist_SA_Per}%',
    fontsize=20,
    fontweight='bold'
)

# then thicken them with a stroke
for lab in [ax.xaxis.label, ax.yaxis.label, ax.title]:
    lab.set_path_effects([pe.withStroke(linewidth=0.6, foreground="black")])

# ===============================
# Blue Borders
# ===============================
for spine in ax.spines.values():
    spine.set_color('blue')
    spine.set_linewidth(3)

# Make tick labels Arial
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontname('Arial')
    label.set_fontsize(16)

plt.tight_layout()
plt.show()

fig.savefig("hist_slice_cap_network.png", dpi=1000, bbox_inches="tight")







end_time = time.time()  # End the timer

print(f"Execution time: {end_time - start_time:.4f} seconds")


# =========================================
# 3D Capillary Network Plot with Cube Frame
# =========================================

from mpl_toolkits.mplot3d.art3d import Line3DCollection

# Extract node positions
pos = nx.get_node_attributes(G, "pos")

# Create list of 3D segments for edges
segments = []
for u, v in G.edges():
    p1 = np.array(pos[u])
    p2 = np.array(pos[v])
    segments.append([p1, p2])

# Create 3D figure
fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(111, projection='3d')
ax.set_anchor('W')

# SHIFT PLOT LEFT (fix whitespace issue)
ax.set_position([0.0, 0.05, 0.80, 0.9])

# Plot edges (red)
line_collection = Line3DCollection(
    segments,
    colors='red',
    linewidths=3,
    alpha=0.9
)
ax.add_collection3d(line_collection)

# Plot nodes (black)
xyz = np.array([pos[n] for n in G.nodes()])
ax.scatter(
    xyz[:, 0],
    xyz[:, 1],
    xyz[:, 2],
    c='black',
    s=15,
    depthshade=True
)

# Axis limits
ax.set_xlim(0, cube_length)
ax.set_ylim(0, cube_length)
ax.set_zlim(0, cube_length)


from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================
# Add histology slice plane
# ==========================
z_slice = float(slice_locations[0])   # uses your existing slice_locations = np.array([250])
L = cube_length

# Corners of the plane (a square spanning the cube in x-y)
plane_verts = [[
    (0, 0, z_slice),
    (L, 0, z_slice),
    (L, L, z_slice),
    (0, L, z_slice)
]]

# Filled translucent plane
plane = Poly3DCollection(
    plane_verts,
    facecolor='royalblue',
    edgecolor='none',   # <-- important
    alpha=0.18
)
ax.add_collection3d(plane)

# Thick outline (like your example)
outline_x = [0, L, L, 0, 0]
outline_y = [0, 0, L, L, 0]
outline_z = [z_slice] * len(outline_x)
ax.plot(
    outline_x,
    outline_y,
    outline_z,
    color='royalblue',     # <-- change here
    linewidth=5
)

# Optional: show slice thickness as a thin "slab" (two planes)
# thickness = 4.0  # your histology thickness
# z1 = z_slice - thickness/2
# z2 = z_slice + thickness/2
# ax.plot(outline_x, outline_y, [z1]*len(outline_x), linewidth=2, alpha=0.6)
# ax.plot(outline_x, outline_y, [z2]*len(outline_x), linewidth=2, alpha=0.6)

# Title (Arial normal)
ax.set_title(
    f'Simulated Capillary Network \n Number of Capillaries : {cur_branch_num}',
    fontsize=20,
    fontweight='bold'
)

# Labels
ax.set_xlabel('X (μm)', fontweight='bold', fontsize=18, labelpad=10)
ax.set_ylabel('Y (μm)', fontweight='bold', fontsize=18, labelpad=10)
ax.set_zlabel('Z (μm)', fontweight='bold', fontsize=18, labelpad=10)

# then thicken them with a stroke
for lab in [ax.xaxis.label, ax.yaxis.label, ax.zaxis.label, ax.title]:
    lab.set_path_effects([pe.withStroke(linewidth=0.6, foreground="black")])

# Make tick numbers larger (and Arial)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.tick_params(axis='z', labelsize=16)

for label in ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels():
    label.set_fontname('Arial')

# =========================================
# Draw Cube Wireframe
# =========================================

L = cube_length

cube_vertices = np.array([
    [0, 0, 0],
    [L, 0, 0],
    [L, L, 0],
    [0, L, 0],
    [0, 0, L],
    [L, 0, L],
    [L, L, L],
    [0, L, L]
])

cube_edges = [
    (0,1),(1,2),(2,3),(3,0),      # bottom
    (4,5),(5,6),(6,7),(7,4),      # top
    (0,4),(1,5),(2,6),(3,7)       # verticals
]

for edge in cube_edges:
    p1 = cube_vertices[edge[0]]
    p2 = cube_vertices[edge[1]]
    ax.plot(
        [p1[0], p2[0]],
        [p1[1], p2[1]],
        [p1[2], p2[2]],
        color='black',
        linewidth=2
    )



ax.grid(False)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Equal cube proportions
try:
    ax.set_box_aspect((1, 1, 1))
except:
    pass

plt.subplots_adjust(left=0, right=0.85)
plt.show()

fig.savefig("capillary_network_v2.png", dpi=1000, bbox_inches="tight", pad_inches=0.5)