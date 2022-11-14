import numpy as np
from scipy import spatial
#import pandas as pd
#import matplotlib.pyplot as plt
#import copy
#import time
#import os
#from os import listdir

class Atom:
    """
    Arguments:
        index(int): the index of the atom
        ele(string): the chemical element of the atom, vaccuum by default
        position(float, len=3 vector): the position of the atom, origin by default
        neighbor(dict{distance(rounded float): [index(int)]}): a dictionary of neighboring atoms 
    """
    def __init__(self, index=0, ele='Va', position=[0.0, 0.0, 0.0], neighbor={}):
        self.index = index
        self.ele = ele
        self.position = np.array(position)
        self.neighbor = neighbor # a list of neighboring positions

    def __repr__(self):
        return f"Atom_{self.index}: {self.ele}"


class Lattice:
    """
    Arguments:
        elements(string, iterable): a list of elements.
        positions(float, iterable): a list of positions.
        lattice_const(float): more like a scaling factor, not a vector.
        lattice_vector(3 by 3 np.array): lattice vectors of the supercell.
        r_max(float): radius cutoff for searching neighbors.
        cal_neighbor(bool): whether calculate the neighboring index of each atom.
    Attributes:
        Lattice.atom_list(Atom, list): a list of Atoms.
    """
    def __init__(self, elements=['Va'], positions=[[0.0, 0.0, 0.0]], lattice_const=1.0, \
                 lattice_vector=np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]), \
                 r_max=0.0, cal_neighbor=True):
        """
        When cal_neighbor=False, do not calculate the atoms' neighbors.
        """
        self.lattice_const = lattice_const
        self.lattice_vector = np.array(lattice_vector)
        self.r_max = r_max
        assert len(elements) == len(positions), "len(elements) == len(positions) not true"
        self.n_atom = len(elements)
        self.atom_list = []
        self.__n_repeat = [0]*3
        for i in range(self.n_atom):
            self.atom_list.append(Atom(index=i, ele=elements[i], position=positions[i]))
        if cal_neighbor:
            self.cal_neighbors()
        
    @staticmethod
    def find_position_vector(position, lattice_vector):
        """
        Auxiliary function to calculate the position vector 
        """
        return position.dot(np.linalg.inv(lattice_vector))
        
    @staticmethod
    def distance(v1, v2):
        """
        Auxiliary function to calculate the distance of two vectors
        """
        v3 = np.array(v1) - np.array(v2)
        return np.linalg.norm(v3)

    def cal_repeats(self, r_max, lattice_vector):
        """
        Calculate the number of repeats needed to expand the lattice
        """
        n_repeat = [0]*3
        norm_vector = np.cross(lattice_vector[0], lattice_vector[1])
        length = norm_vector.dot(lattice_vector[2])/np.linalg.norm(norm_vector)
        n_repeat[2] = int(np.ceil(r_max/abs(length)))

        norm_vector = np.cross(lattice_vector[1], lattice_vector[2])
        length = norm_vector.dot(lattice_vector[0])/np.linalg.norm(norm_vector)
        n_repeat[0] = int(np.ceil(r_max/abs(length)))

        norm_vector = np.cross(lattice_vector[2], lattice_vector[0])
        length = norm_vector.dot(lattice_vector[1])/np.linalg.norm(norm_vector)
        n_repeat[1] = int(np.ceil(r_max/abs(length)))
        
        return n_repeat
    
    def expand_data(self, data):
        """
        Expand data so that the neighboring atoms within r_max is always included in data_expand
        """
        data_expand = []
        lattice_vector = self.lattice_vector
        n_repeat = self.cal_repeats(self.r_max, lattice_vector)
        self.__n_repeat = n_repeat
        for k in range(-n_repeat[2], n_repeat[2]+1):
            for j in range(-n_repeat[1], n_repeat[1]+1):
                for i in range(-n_repeat[0], n_repeat[0]+1):
                    position_shift =(i*lattice_vector[0] + j*lattice_vector[1] + k*lattice_vector[2])
                    data_expand.append(np.array(data) + position_shift)
        return np.array(data_expand).reshape(-1, 3)

    def cal_neighbors(self):
        """
        Calculate the neighboring atom index of each atom in format [{distance: [index]}]
        neighbor_list: a list of neighboring positions for each atom, calculated with k-d tree.
        see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
        """
        data_position = np.array([i.position for i in self.atom_list])
        data_expand = self.expand_data(data_position)
        tree = spatial.KDTree(data_expand)
        neighbor_list = []
        for i in data_position:
            neighbor_list.append([j for j in tree.query_ball_point(i, self.r_max)])            
        for i in range(len(self.atom_list)):
            dict_tmp = {}
            for index in neighbor_list[i]:
                rd = self.distance(data_position[i], data_expand[index])
                key = round(rd, 10) # round off float point number to get the key
                try:
                    #dict_tmp[key].append(index%self.n_atom)
                    dict_tmp[key].append(index)
                except:
                    #dict_tmp[key] = [index%self.n_atom]
                    dict_tmp[key] = [index]
            for key in dict_tmp:
                tmp = sorted(dict_tmp[key]) # sort neighboring atoms according to index (expanded lattice)
                dict_tmp[key] = [i%self.n_atom for i in tmp] # convert expanded index to smaller index within the lattice
            self.atom_list[i].neighbor = dict_tmp
            
    def reset_r_max(self, r_max):
        """
        Reset the value of self.r_max, and call cal_neighbors to calculate the neighbors.
        """
        self.r_max = r_max
        self.cal_neighbors()
        
    def update_elements(self, elements):
        for i in range(self.n_atom):
            self.atom_list[i].ele = elements[i]
            
    def update_positions(self, positions):
        """
        Recalculate the neighboring positions
        """
        for i in range(self.n_atom):
            self.atom_list[i].position = positions[i]
        self.cal_neighbors(self.r_max) 
            
    def print_info(self):
        print("="*50)
        print("Parameters of class Lattice")
        print("="*50)
        print(f"lattice_const: {self.lattice_const}")
        print(f"lattice_vector: \n{self.lattice_vector}")
        print(f"n_atom: {self.n_atom}")
        print(f"r_max: {self.r_max}")
        print(f"atom_list: {self.atom_list}")
        
    def write_positions_LSMS2(self, file_position):
        with open(file_position, "w+") as file:
            file.write("#Lattice constant:\n")
            file.write("%16.12f\n"%self.lattice_const)
            file.write("# Lattice vector:\n")
            for i in range(3):
                file.write("%16.12f\t%16.12f\t%16.12f\n"%tuple(self.lattice_vector[i]))
            file.write("# Number of medium atoms:\t%d\n"%self.n_atom)
            n_species = len(set([atom.ele for atom in self.atom_list]))
            file.write("# Number of medium atom type:\t%d\n"%n_species)
            for atom in self.atom_list:
                file.write("%s\t%16.12f\t%16.12f\t%16.12f\n"%(atom.ele, atom.position[0], atom.position[1], atom.position[2]))
                
    def write_positions_xyz(self, file_position):
        with open(file_position, "w+") as file:
            file.write("%d\n"%self.n_atom)
            n_species = len(set([atom.ele for atom in self.atom_list]))
            file.write("# Number of medium atom type:\t%d\n"%n_species)
            file.write("# Lattice constant:\t%f\n"%self.lattice_const)
            for atom in self.atom_list:
                file.write("%s\t%16.12f\t%16.12f\t%16.12f\n"%(atom.ele, atom.position[0], atom.position[1], atom.position[2]))
                
    def write_neighbors(self, file_position):
        with open(file_position, "w+") as file:
            for key in sorted(self.atom_list[0].neighbor.keys())[1:]: #ignore the self neighbor
                for index in self.atom_list[0].neighbor[key]:
                    file.write(str(key)+',')
            file.write('\n')
                    
            for atom in self.atom_list:
                for key in sorted(atom.neighbor.keys())[1:]: # ignore the self neighbor
                    for index in atom.neighbor[key]:                        
                        file.write("%d,"%index)
                file.write("\n")

    def cal_EPI_Xdata(self):
        """
        Calculate the X data for fitting the EPI model.
        """
        elements = set([i.ele for i in self.atom_list])
        pairs_record = {}
        keys_sorted = sorted(self.atom_list[0].neighbor.keys())[1:]
        n_shell = len(keys_sorted)
        for i_shell in range(n_shell):
            for i in elements:
                for j in elements:
                    k12 = sorted([i, j])            
                    key = str(i_shell)+'_'+k12[0]+k12[1]
                    pairs_record[key] = 0
        for atom in self.atom_list:
            ele = atom.ele
            for i in range(len(keys_sorted)):
                key = keys_sorted[i]
                atoms_tmp = atom.neighbor[key]
                for index in atoms_tmp:
                    ele_2 = self.atom_list[index].ele
                    k12 = sorted([ele, ele_2])
                    pair_key = str(i)+'_'+k12[0]+k12[1]
                    pairs_record[pair_key] += 1
        total_list = [0]*n_shell
        total_list
        for i in pairs_record:
            total_list[int(i[0])] += pairs_record[i]
        total_list    

        data_out = []
        for i in range(n_shell):
            data_out.append([])

        for i in pairs_record:
            i_shell = int(i[0])
            data_out[i_shell].append(pairs_record[i]/total_list[i_shell])
        return data_out
            
    def write_EPI_Xdata(self, file_position):
        data_out = self.cal_EPI_Xdata()
        with open(file_position, "w+") as file:
            for i in data_out:
                for j in i:
                    file.write("%f,"%j)
                file.write("\n")

            


 
