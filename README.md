# MonteCarlo3D

#### Description
**3-dimensional Monte Carlo simulation of high entropy alloys**
* Cannonical ensemble (fixed temperature).
* Effective pair interaction (EPI) model.
* Temperature parallelization.

For more information about EPI and high entropy alloys, you may check the following papers:
* [Robust data-driven approach for predicting the configurational energy of high entropy alloys](https://www.sciencedirect.com/science/article/pii/S0264127519306859)
* [Monte Carlo simulation of order-disorder transition in refractory high entropy alloys: A data-driven approach](https://www.sciencedirect.com/science/article/abs/pii/S0927025620306261)

#### Software Architecture
Software architecture description

#### Installation

Make sure you have the following packages already installed:
* gcc: support -std=c++11.
* MPI: support mpicxx, OpenMPI is used for testing, but the program should also works for MPICH.
* python: >=3.0, include numpy, scipy, and matplotlib.
* If you want to use the notebook, please install Jupyter Notebook.

#### Instructions
Example:
1. Generate "neighbors.dat" in "./Python/neighbor_tmp" via
   ```bash
   python ./Python/genNeighbors.py
   ```
   * Note that the first line in "neighbors.dat" is the distances of neighboring atoms to the centering atom. The radiu-cutoff r_max can be set in the main function of "genNeighbors.py". 
   * The supercell size can also be set in the main function of "genNeighbors.py" by setting "n_repeats=[4,4,4]" to other values (by default, ).
   * The "EPI.dat" are the energy model parameters for the pair-interaction model. They can be obtained from the linear regression of energy vs EPI. By default, it includes at most 6 neighboring shells.

2. Run Monte Carlo simulation via
   ```bash
   cd ./C++
   mpicxx -O3 -std=c++11 MonteCarlo.cpp
   mpirun -n 4 ./a.out
   ```
    * The supercell size is obtained from the "neighbors.dat" mentioned above.
    * Temperature parallelization is implemented with MPI. The following variables can be set in the main funciton:
      * T_min and T_max: minimum and maximum temperatures.
      * n_T: the number of different temperatures. parallelized with MPI.
      * n_measure: the number of sweeps used for making measurement.
      * n_discard: the number of sweeeps discarded between each measurement.
      * You can save the output to a file with
        ```bash
        mpirun -n 4 ./a.out > ../Output/output.dat
        ``` 
        Note that the outputs start with a bunch of TEST, and the only the last TEST carries out the Monte Carlo simulation 
3. You can also ignore the above steps, and simply run 
   ```bash
   ./run.sh
   ```

#### Contribution

1.  Fork the repository
2.  Create Feat_xxx branch
3.  Commit your code
4.  Create Pull Request

#### Contact Info
xianglinliu01@gmail.com