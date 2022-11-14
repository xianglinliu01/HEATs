import numpy as np
import os
from GenPosition import generate_random_position

def genNeighbors_BBC(n_repeats=[1,1,1], r_max=1.0, file_dir="./", elements=("Mo","Nb","Ta","W")):
    """
    Arguements: 
        n_repeats: the repeats along the 3 directions of the BCC unit cell
        r_max: the cutoff for counting neighbors, using lattice unit (i.e.,lattice_constant = 1.0)
        file_dir: the relative location of the file
    Return:
        a saved file "neighbors.dat" at the designated location
    """
    n_atom = 2 * np.product(n_repeats) # two atoms in a BCC unit cell
    print(f"n_atom={n_atom}")
    print(f"r_max={r_max}")
    
    ##### try to evenly split atoms for each element, also works if n_atom%n_ele != 0
    ele_dict = {}
    for i in elements:
        ele_dict[i]=0
    n_ele = len(elements)
    for i in range(n_atom):
        ele_dict[elements[i%n_ele]] += 1

    lattice_new = generate_random_position(positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],\
                        lattice_vector=np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]),\
                        ele_dict=ele_dict, n_repeat=n_repeats)
    lattice_new.reset_r_max(r_max)
    file_neighbor = os.path.join(os.getcwd(), file_dir, "neighbors.dat")
    lattice_new.write_neighbors(file_neighbor)
    file_dat = os.path.join(os.getcwd(), file_dir, "positions.xyz")
    lattice_new.write_positions_xyz(file_dat)

def genNeighbors_FCC(n_repeats=[1,1,1], r_max=1.0, file_dir="./", elements=("Co","Cr","Fe","Ni")):
    """
    Arguements:
        n_repeats: the repeats along the 3 directions of the FCC unit cell
        r_max: the cutoff for counting neighbors, using lattice unit (i.e.,lattice_constant = 1.0)
        file_dir: the relative location of the file
    Return:
        a saved file "neighbors_FCC.dat" at the designated location
    """
    n_atom = 4 * np.product(n_repeats) # 4 atoms in a FCC unit cell
    print(f"n_atom={n_atom}")
    print(f"r_max={r_max}")

    ##### try to evenly split atoms for each element, also works if n_atom%n_ele != 0
    ele_dict = {}
    for i in elements:
        ele_dict[i]=0
    n_ele = len(elements)
    for i in range(n_atom):
        ele_dict[elements[i%n_ele]] += 1

    lattice_new = generate_random_position(positions=[[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],\
                        lattice_vector=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),\
                        ele_dict={'Va':n_atom}, n_repeat=n_repeats)
    lattice_new.reset_r_max(r_max)
    file_neighbor = os.path.join(os.getcwd(), file_dir, "neighbors_FCC.dat")
    lattice_new.write_neighbors(file_neighbor)

def genNeighbors_SC(n_repeats=[1,1,1], r_max=1.0, file_dir="./", elements=("Va")):
    """
    Arguements:
        n_repeats: the repeats along the 3 directions of the SC unit cell
        r_max: the cutoff for counting neighbors, using lattice unit (i.e.,lattice_constant = 1.0)
        file_dir: the relative location of the file
    Return:
        a saved file "neighbors_SC.dat" at the designated location
    """
    n_atom = 1 * np.product(n_repeats) # 1 atom in a simple cubic unit cell
    print(f"n_atom={n_atom}")
    print(f"r_max={r_max}")

    ##### try to evenly split atoms for each element, also works if n_atom%n_ele != 0
    ele_dict = {}
    for i in elements:
        ele_dict[i]=0
    n_ele = len(elements)
    for i in range(n_atom):
        ele_dict[elements[i%n_ele]] += 1

    lattice_new = generate_random_position(positions=[[0.0, 0.0, 0.0]],\
                        lattice_vector=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),\
                        ele_dict={'Va':n_atom}, n_repeat=n_repeats)
    lattice_new.reset_r_max(r_max)
    file_neighbor = os.path.join(os.getcwd(), file_dir, "neighbors_SC.dat")
    lattice_new.write_neighbors(file_neighbor)

if __name__ == '__main__':
    dir_tmp = "./neighbor_BCC"
    os.makedirs(dir_tmp, exist_ok=True)
    file_dir = os.path.join(os.getcwd(), dir_tmp)
    print(f"save the neighbors data at {file_dir}")
    genNeighbors_BBC(n_repeats=[6,6,6], r_max=1.01, file_dir=dir_tmp)
    #genNeighbors_FCC(n_repeats=[2,2,2], r_max=1.5, file_dir=dir_tmp)
    #genNeighbors_SC(n_repeats=[2,2,2], r_max=1.5, file_dir=dir_tmp)
