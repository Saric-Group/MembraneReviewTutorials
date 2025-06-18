import numpy as np
import random as rnd
from numpy import arange, pi, sin, cos, arccos
import sys

def generate_triangular_lattice(n_rows, n_cols, spacing=1.0):

    """
    Generate a triangular lattice of points.
    
    Parameters:
    - n_rows: Number of rows of points in the lattice.
    - n_cols: Number of columns of points in the lattice.
    - spacing: Distance between adjacent points.
    
    Returns:
    - vertices: Array of point positions.
    - faces: Array of faces (triangles) as indices into the vertices array.
    """

    vertices = []
    faces = []

    top_line = []
    bottom_line = []
    left_line = []
    right_line = []

    # this gives points in a triangular lattice
    assigned_num = 0
    vertices_at_edge = []

    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spacing
            y = row * spacing * np.sqrt(3) / 2
            if row % 2 == 1:
                x += spacing/2
            vertices.append(np.array([x, y, 0, assigned_num]))

            if(col==0) or (col == n_cols-1):
                if assigned_num not in vertices_at_edge:
                    vertices_at_edge.append(assigned_num)

            if (row == 0) or (row == n_rows-1):
                if assigned_num not in vertices_at_edge:
                    vertices_at_edge.append(assigned_num)
            
            if (col==0 and row >0 and row < n_rows-1):
                left_line.append(assigned_num)
            if (col == n_cols-1 and row >0 and row < n_rows-1):
                right_line.append(assigned_num)
            if (row == 0 ):
                bottom_line.append(assigned_num)
            if (row == n_rows-1):
                top_line.append(assigned_num)

            assigned_num+=1
                
    vertices = np.array(vertices)

    # get the faces
    for row in range(n_rows):
        for col in range(n_cols - 1):
            current = row * n_cols + col
            right = current + 1
            top = (row+1) * n_cols + col
            below = (row-1) * n_cols + col
            
            if row%2==1:
                top +=1
                below += 1
            if row==0:
                if [current, right, top] not in faces:
                    faces.append([current, right, top])
            elif row>0 and row<n_rows-1:
                if [current, right, top] not in faces:
                    faces.append([current, right, top])
                if [current, below, right] not in faces:
                    faces.append([current, below, right])
            else:
                if [current, below, right] not in faces:
                    faces.append([current, below, right])
    faces = np.array(faces)
    
    return vertices, faces, np.array(vertices_at_edge), top_line, bottom_line, left_line, right_line

def prepare_config_hexagonal_flatpatch(network_columns, seed, displace_z, scale_lattice = 1.5, radius = 0.8, sigma_nanops = 10):

    # generate the lattice
    network_rows = int(2/np.sqrt(3) * network_columns)
    lattice_vertices, lattice_faces, vertices_at_edge, top_line, bottom_line, left_line, right_line = generate_triangular_lattice(network_rows, network_columns)

    xcm = np.mean(lattice_vertices[:, 0])
    ycm = np.mean(lattice_vertices[:, 1])
    zcm = np.mean(lattice_vertices[:, 2])

    lattice_vertices[:, 0] -= xcm 
    lattice_vertices[:, 1] -= ycm 
    lattice_vertices[:, 2] -= zcm 

    # set the scale of the lattice
    lattice_vertices[:, 0] *=scale_lattice
    lattice_vertices[:, 1] *=scale_lattice
    lattice_vertices[:, 2] *=scale_lattice
    Nbeads = len(lattice_vertices)

    # extent of the patch
    extent_patch = np.max(lattice_vertices[:, 0])
    scaled_extent = radius * extent_patch

    # cut a circle/hexagon from the patch
    distances_patch = np.sqrt(lattice_vertices[:, 0]**2 + lattice_vertices[:, 1]**2)
    distances_in_circle = np.where(distances_patch>scaled_extent)[0]
    edge_vertices = lattice_vertices[distances_in_circle]

    # note that by moving this you will also take the rest down
    lattice_vertices[:, 2] -= displace_z

    # impose external seed
    rnd.seed(seed)

    print("Writing down the file --> "+ "flatpatch_N_"+str(Nbeads)+"_Sigma_"+str(int(sigma_nanops))+"_.dat")
    file = open("flatpatch_N_"+str(Nbeads)+"_Sigma_"+str(int(sigma_nanops))+"_.dat", "w")

    # write down the file
    # write it such that the vertices at edge are ordered
    central_bead_identification=np.zeros((len(lattice_vertices), 3))
    j = 0
    for ver in distances_in_circle:
        file.writelines("{}\t{}\t{}\t{}\n".format(1, lattice_vertices[ver, 0], lattice_vertices[ver, 1], lattice_vertices[ver, 2]))
        central_bead_identification[j, 0] = lattice_vertices[ver, 0]
        central_bead_identification[j, 1] = lattice_vertices[ver, 1]
        central_bead_identification[j, 2] = lattice_vertices[ver, 2]

        j+=1
    
    for i in range(Nbeads):
        if i not in distances_in_circle:
            file.writelines("{}\t{}\t{}\t{}\n".format(1, lattice_vertices[i, 0], lattice_vertices[i, 1], lattice_vertices[i, 2]))
            central_bead_identification[j, 0] = lattice_vertices[i, 0]
            central_bead_identification[j, 1] = lattice_vertices[i, 1]
            central_bead_identification[j, 2] = lattice_vertices[i, 2]
            j+=1
    
    # what is the membrane particle closest to the center?
    distances = np.sqrt(central_bead_identification[:, 0]**2 + central_bead_identification[:, 1]**2)
    index_max_z = np.where(distances==np.min(distances))[0][0]

    file.writelines("{}\t{}\t{}\t{}\n".format(2, central_bead_identification[index_max_z, 0], central_bead_identification[index_max_z, 1], 0.5*(sigma_nanops+1)+central_bead_identification[index_max_z, 2]))
    file.close()

    """
    # OVITO FILE
    file = open("test.dump", "w")

    file.writelines("ITEM: TIMESTEP\n")
    file.writelines("0\n")
    file.writelines("ITEM: NUMBER OF ATOMS\n")
    file.writelines(str(Nbeads+1)+"\n")
    file.writelines("ITEM: BOX BOUNDS pp pp pp\n")
    file.writelines("-100.0 100.0\n")
    file.writelines("-100.0 100.0\n")
    file.writelines("-100.0 100.0\n")
    file.writelines("ITEM: ATOMS id type x y z\n")
    count_vertices = 0

    count_vertices=0
    for i in range(len(edge_vertices)):
        file.writelines("{} {} {} {} {}\n".format(count_vertices+1, 2, edge_vertices[i, 0], edge_vertices[i, 1], edge_vertices[i, 2]))
        count_vertices+=1
    
    for i in range(Nbeads):
        if i not in distances_in_circle:
            file.writelines("{} {} {} {} {}\n".format(count_vertices+1, 1, lattice_vertices[i, 0], lattice_vertices[i, 1], lattice_vertices[i, 2]))
            count_vertices+=1

    file.writelines("{} {} {} {} {}".format(Nbeads+1, 3, central_bead_identification[index_max_z, 0], central_bead_identification[index_max_z, 1], 0.5*(sigma_nanops+1)+central_bead_identification[index_max_z, 2]))
    
    print("COUNT VERTICES : ", count_vertices)
    print("Vertices at edge:", len(distances_in_circle))
    print("INDEX MAX Z + 1: ", index_max_z)
    print(f"EXTENSION BOX: {2*np.max(lattice_vertices[:, 0])}")

    file.close()
    """

    return Nbeads, len(edge_vertices), index_max_z+1, 2*np.max(lattice_vertices[:, 0])
