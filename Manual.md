# 3DMonteCarlo 0.5

## C++ Directory
#### Utilitiy functions
* class Timer(): record time elapsed.
  * double getTime()
  * void printTime()
  * void resetTime()
* namespace Variance: calculate the variance of a std::vector list.
  * double mean( const vector<double>& list )
  * double variance(const vector<double>& list )
* Read_csv: Read csv format data.
  * vector< vector<string> > read_csv(string filename, char delim=',')
* namespace Random: Generate random numbers with Mersenne Twister. Can specify see or not with different constructor: Random_mt() or Random_mt(unsigned int seed).
  * class Random_mt: double getRandom()
  * class Random_mt_int: 
* namespace TEST_utilities: functions for convenience
  * void TEST_PRINT_LINE(std::string title)
  * void TEST_ASSERT(bool condition, const char* message)
  
#### Test functions
* void TEST_getRandom()
* void TEST_shuffle()
* void TEST_Neighbor(const string &filename)
* void TEST_calpair(const string & filename)
* void TEST_MonteCarlo(vector<string> elements=vector<string> {"E1", "E2", "E3", "E4"})

#### MonteCarlo.cpp
* struct MonteCarlo: A class for Monte Carlo simulation;
  * MonteCarlo(Configuration config_in, Neighbor nb_in, vector< vector<double> > para_EPI_in)

* struct Configuration:
  * Constructor for equal chemical concentrations.
  * vector< vector<int> > cal_pair(const Neighbor &nb)
  * 
* struct Neighbor: output the neighboring atoms for all atoms in the configuration; 
* ../Python/EPI.dat: the effective pair interaction parameters 
* ../Python/neighbors.dat: the neighboring atoms for all atom; note the first line is the distance between atoms. This file is outputed from HEATS_tutorial.ipynb

* struct Observables: The observables from Monte Carlo simulation, such as magnetization, specific heats, etc.

## Python Directory
#### Lattice.py
* class Atom: return element, position, and its neighbors (optional)
* class Lattice: A relatively large class. The main role is to establish a super cell lattice and find  the neighbors of each atom via the highly efficient [KD tree algorithm](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html). The radius cutoff can be set via the "reset_r_max" method. The configurations can be outputed to LSMS and xyz format via the "write_positions_LSMS2" and "write_positions_xyz" methods.

#### GenPosition.py 
Mainly used to produce different random configurations.
* shuffle_elements(ele_list, n_list, seed=None): if seed is given, always use the same seed; otherwise use random.shuffle, which gives different values everytime the function being called.
* generate_random_position(): the main function needed. See the Docstrings for details.

#### genNeighbors.py
The interface function. This is the only place you are supposed to modify for different supercell size, r_max, n_atom, etc.
Note that BCC structure is used by default. You can use FCC or SC by calling the correspoinding functions.
* genNeighbors_BBC(n_repeats=[1,1,1], r_max=1.0, file_dir="./"): Generate BCC supercells. n_repeats along the BCC unit cell to enlarge the cell size. r_max set the radius cutoff, and file_dir is the output directory.
* genNeighbors_BCC(n_repeats=[1,1,1], r_max=1.0, file_dir="./"): The same as above, but for BCC structure.
* genNeighbors_SC(n_repeats=[1,1,1], r_max=1.0, file_dir="./"): The same as above, but for simple cubic structure.