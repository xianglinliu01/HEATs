import random
import numpy as np
import os
from Lattice import Lattice

def shuffle_elements(ele_list, n_list, seed=None):
    """
    Arguments:
        ele_list: a list of different elements 
        n_list: number of atoms for each element
    Return:
        elements: a list 
    """
    assert len(ele_list) == len(n_list), "len(ele_list) == len(n_list) not true in shuffle_elements"
    elements = []
    for i in range(len(ele_list)):
        elements += n_list[i]*[ele_list[i]]
    if seed:
        random.seed(seed)
    random.shuffle(elements)
    return elements

def expand_lattice(lattice, ele_list, n_list, n_repeat):
    """
    Expand the small unit cell to large supercell
    """
    assert np.sum(n_list) == lattice.n_atom*np.product(n_repeat), \
    f"np.sum(n_list) {np.sum(n_list)} == lattice.n_atom*np.product(n_repeat) {lattice.n_atom*np.product(n_repeat)} not true in expand_lattice"
    lattice_vector = lattice.lattice_vector
    data_position = np.array([i.position for i in lattice.atom_list])
    data_expand = []
    for k in range(n_repeat[2]):
        for j in range(n_repeat[1]):
            for i in range(n_repeat[0]):
                position_shift =(i*lattice_vector[0] + j*lattice_vector[1] + k*lattice_vector[2])
                for position in data_position + position_shift:
                    data_expand.append(position)
    data_expand = np.array(data_expand)
    elements = shuffle_elements(ele_list, n_list)
    lattice_vector_new = np.array([lattice_vector[0]*n_repeat[0], \
                                   lattice_vector[1]*n_repeat[1],\
                                   lattice_vector[2]*n_repeat[2]])                                                                                             
    return Lattice(elements=elements, positions=data_expand, \
                   lattice_const=lattice.lattice_const, \
                   lattice_vector=lattice_vector_new, cal_neighbor=False)

def generate_random_position(positions=[[0.0, 0.0, 0.0]], \
                             lattice_vector=[[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]],\
                             ele_dict={'Va':1}, \
                             n_repeat=[1,1,1], lattice_const=1.0):
    """
    Generate the position file and save it at file_postion.
    Example:
        1: generate_random_position(file_position)
            Defualt, 1 Vaccum atom cubic cell
        2: generate_random_position(file_position, positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],\
                        lattice_vector=np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]),\
                        ele_dict={'Mo':8, 'V':8}, n_repeat=[2, 2, 2])
            2 by 2 by 2 BCC cell, a total of 16 atoms, with 8 Mo and 8 V atoms.
        3: position_fcc = [[0,0,0],[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
           generate_random_position(file_position, positions = position_fcc, \
                        lattice_vector = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]),\
                        ele_dict={'Mo':16, 'V':16}, n_repeat=[2,2,2])
            2 by 2 by 2 FCC cell, a total of 32 atoms, with 16 Mo and 16 V atoms.
    """
    n_atom_small = len(positions)
    lattice = Lattice(elements=['Va']*n_atom_small, positions=positions, \
                      lattice_const=lattice_const,\
                      lattice_vector=lattice_vector, cal_neighbor=False)
    ele_list = list(ele_dict.keys())
    n_list = list(ele_dict.values())
    lattice_new = expand_lattice(lattice, ele_list, n_list, n_repeat)
    lattice_new.print_info()
    return lattice_new

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Short sample app')
	parser.add_argument('-a', action="store", default='Mo, Nb')
	parser.add_argument('-b', action="store", default="4, 4")
	parser.add_argument('-c', action="store", default="2, 2, 2")
	parser.add_argument('-d', action="store", default="6.2")
	parser.add_argument('-e', action="store", default="")
	results = parser.parse_args()

	elements = [i.strip() for i in results.a.split(sep=',')]
	n_ele = [int(i) for i in results.b.split(sep=',')]
	n_repeat = [int(i) for i in results.c.split(sep=',')]
	lattice_const = float(results.d)
	assert len(elements) == len(n_ele), "len(elements) == len(n_ele) not true"
	ele_dict = {}
	for i in range(len(elements)):
		ele_dict[elements[i]] = n_ele[i]
	out_dir = os.path.join(os.getcwd(), results.e)

	lattice_new = generate_random_position(positions=[[0.0, 0.0, 0.0]],\
	                        lattice_vector=np.array([[0.5, 0.5, -0.5],[0.5, -0.5, 0.5],[-0.5, 0.5, 0.5]]),\
	                        ele_dict=ele_dict, n_repeat=n_repeat, lattice_const=lattice_const)
	file_dat = os.path.join(out_dir, 'position.dat')
	lattice_new.write_positions_LSMS2(file_dat)
	# file_xyz = os.path.join(out_dir, 'positions.xyz')
	# lattice_new.write_positions_xyz(file_xyz)
