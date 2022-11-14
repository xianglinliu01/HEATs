#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <sstream> 
#include <stdexcept>
#include <map>
#include "Random.hpp"
#include "Timer.hpp"
#include "TEST_utility.hpp" 
#include "Print_vector.hpp"
#include "Read_csv.hpp"
#include "Variance.hpp"
#include "Write_config.hpp"
#include <mpi.h>

//#define TEST_ALL

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::array;
using Random::Random_mt;
using Random::Random_mt_int;
using TEST_utility::TEST_PRINT_LINE;
using TEST_utility::TEST_ASSERT;
using Print_vector::print_vector;
using Variance::variance;
using Variance::mean;

//Random_mt random_mt;
#define SEED_g 2020
#ifdef SEED_g
    Random_mt random_mt(SEED_g);
#else
    Random_mt random_mt;
#endif 

void TEST_getRandom()
{
    TEST_PRINT_LINE("TEST_getRrandom");
    int total = 1000000;
    int count = 0;
    for (int i=0; i<total; i++)
    {
        double x = random_mt.getRandom();
        double y = random_mt.getRandom();
        if (x*x + y*y < 1.0)
        {
            count+=1;
        }
    }
    printf("Pi from Monte Carlo = %10.8f\n", 4*float(count)/total);
}

void TEST_shuffle()
{
    TEST_PRINT_LINE("TEST_shuffle");
    int n_ele = 4;
    int n_atom = 16;
    vector<int> atoms;
    for(int i=0; i<n_atom; i++)
    {
        atoms.push_back(i/n_ele);
    }
    print_vector(atoms);
    std::random_shuffle(atoms.begin(), atoms.end());
    print_vector(atoms);
}

struct Neighbor
{
public:
    Neighbor(const string &filename)
    {
        vector< vector<string> > data = read_csv(filename);
        shell_list = calShell(data);
        neighbors = calNeighbor(data);
        n_shell = *std::max_element(shell_list.begin(), shell_list.end()) + 1;
        n_atom_shell = vector<int> (n_shell, 0);
        for (const auto &x: shell_list)
        {
            n_atom_shell[x] += 1;
        }

    }

    int n_shell;
    vector<int> n_atom_shell;

    const vector<int> &get_shell_list() const
    {
        return shell_list;
    }

    const vector< vector<int> > &get_neighbors() const
    {
        return neighbors;
    }

private:
    vector<int> shell_list;
    vector< vector<int> > neighbors;

    vector<int> calShell(const vector< vector<string> > &data)
    {
        vector<string> v1 = data[0];
        int index = 0;
        vector<int> result;
        result.push_back(0);
        for (int i=1; i<v1.size(); i++)
        {
            if (v1[i] != v1[i-1]) index += 1;
            result.push_back(index);
        }
        return result;
    }

    vector< vector<int> > calNeighbor(const vector< vector<string> > &data)
    {
        vector< vector<int> > result;
        for (int i=1; i<data.size(); i++)
        {
            vector<int> row;
            for (int j=0; j<data[i].size(); j++)
            {
                row.push_back( std::stoi(data[i][j]) );
            }
            result.push_back(row);
        }
        return result;
    }
};

void TEST_Neighbor(const string &filename)
{
    TEST_PRINT_LINE("TEST_Neighbor");
    Neighbor nb(filename);
    cout << "n_shell= " << nb.n_shell << endl;
    cout << "n_atom_shell= " << endl;
    print_vector(nb.n_atom_shell);
    cout << "shell_list:" << endl;
    print_vector(nb.get_shell_list());
    cout << "neighbors:" << endl;
    for (auto x: nb.get_neighbors())
    {
        print_vector(x);
    }
}


// Get the pair index for element i and j for a n-component system.
int get_pair_index(int ind_i, int ind_j, int n)
{
    int ind_min = std::min(ind_i, ind_j);
    int ind_max = std::max(ind_i, ind_j);
    int index = n*(n + 1)/2 - (n - ind_min)*(n - ind_min + 1)/2
            + (ind_max-ind_min);
    return index;
}

struct Configuration
{
public:
    Configuration(): n_atom(0), n_ele(0), n_pair(0){;}

    // Constructor for equal chemical concentrations.
    Configuration(vector<string> elements_in, int n_atom_in, bool fixSeed=true):
    elements(elements_in), n_atom(n_atom_in), n_ele(elements_in.size()), 
    n_pair(n_ele*(n_ele + 1)/2)
    {
        // Sort the elements to make sure they are ordered.
        std::sort(elements.begin(), elements.end());

        // Make sure the same number of atoms for each element.
        TEST_ASSERT(n_atom%n_ele == 0, 
            "n_atom%n_ele != 0");

        for(int i=0; i<n_atom; i++)
        {
            atoms.push_back(i%n_ele);
        }

        shuffleAtoms(fixSeed);
    }

    // Constructor for nonequal chemical configurations
    Configuration(vector<string> elements_in, vector<int> n_atom_vec, bool fixSeed=true):
    elements(elements_in), 
    n_atom(std::accumulate(n_atom_vec.begin(), n_atom_vec.end(), 0)), 
    n_ele(elements_in.size()), n_pair(n_ele*(n_ele + 1)/2)
    {
        // Sort the elements to make sure they are ordered.
        std::sort(elements.begin(), elements.end());

        // Make sure elements and n_atom_vec have the same length.
        TEST_ASSERT(n_atom_vec.size() == elements.size(), 
            "n_atom_vec.size() != elements.size()");

        for (int i=0; i<elements.size(); i++)
        {
            for (int j=0; j<n_atom_vec[i]; j++)
            {
                atoms.push_back(i);
            }
        }

        shuffleAtoms(fixSeed);
    }
    vector<int> atoms;
    vector<string> elements;
    const int n_atom;
    const int n_ele;
    const int n_pair;

    vector< vector<int> > cal_pair(const Neighbor &nb)
    {
        vector< vector<int> > pair_vec;
        // Initialized pair_vec[n_shell][n_pair] with zeros. 
        for (int i=0; i<(nb.n_shell); i++)
        {
            vector<int> tmp_vec(n_pair, 0);
            pair_vec.push_back(tmp_vec);
        }
        // Make sure the neighbor file has the correct number of atoms.
        TEST_ASSERT(nb.get_neighbors().size() == n_atom, "nb.neighbors != n_atom");
        const vector<int> &shell_vec = nb.get_shell_list();
        for (int i=0; i<n_atom; i++)
        {
            int ind_i = atoms[i];
            const vector<int> &neighbor_i = nb.get_neighbors()[i];
            for (int j=0; j<neighbor_i.size(); j++)
            {
                int ind_j = atoms[neighbor_i[j]];
                int index = get_pair_index(ind_i, ind_j, n_ele);
                int ind_shell = shell_vec[j]; 
                pair_vec[ind_shell][index] +=1;
            }
        }
        return pair_vec;
    }

private:
    // Using a random generator to shuffle the atoms. 
    // Use fixed Seed or not depending on fixSeed.
    void shuffleAtoms(bool fixSeed)
    {
        if (fixSeed)
        {
            std::mt19937 mt(SEED_g);
            std::shuffle(atoms.begin(), atoms.end(), mt);
        }
        else
        {
            std::random_device rd;
            std::mt19937 mt(rd());
            std::shuffle(atoms.begin(), atoms.end(), mt);
        }   
    }

};

// random_mt_int defined for random select. 
string EPI_file = "../Python/neighbor_BCC/EPI.dat";
string neighborFile = "../Python/neighbor_BCC/neighbors.dat";
    Neighbor nb_g(neighborFile);
#ifdef SEED_g
    Random_mt_int random_mt_int(0, nb_g.n_atom_shell[0]-1, SEED_g);
#else
    Random_mt_int random_mt_int(0, nb_g.n_atom_shell[0]-1);
#endif

struct Observables
{
public:
    Observables(){;}

    Observables(int n_shell, int n_SRO):
    SRO_parameters(vector< vector<double> >(n_shell, vector<double>(n_SRO, 0.0)) ){;}

    double temperature;
    double specific_heat;
    vector< vector<double> > SRO_parameters;

    void print_info()
    {
        printf("temperature=%12.10f\n", temperature);
        printf("specific_heat=%12.10f\n", specific_heat);
        printf("SRO_parameters:\n");
        for (auto x: SRO_parameters)
        {
            print_vector(x);
        }
    }
};

struct MonteCarlo
{
public:
    MonteCarlo(Configuration config_in, Neighbor nb_in, vector< vector<double> > para_EPI_in):
    config(config_in), nb(nb_in), para_EPI(para_EPI_in), kb(0.0861733)
    {
        n_shell = nb.n_shell;
        n_pair = config.n_pair;
        n_atom = config.n_atom;
        accept_count = 0;

        energy_total = cal_energy();
        //cout << "energy_total=" << energy_total << endl;

        //cout << "accept_count=" << accept_count << endl;
    }

    Configuration config;
    Neighbor nb;
    vector< vector<double> > para_EPI;
    int n_pair;
    int n_shell;
    int n_atom;

    vector<double> energy_ave;
    vector<double> accept_ratio_ave;

    const double kb;

    void TEST_random_swap()
    {
        TEST_PRINT_LINE("TEST_random_swap");
        printf("Before swap:\n");
        print_vector(config.atoms);
        double T = 100.0;
        int i = 0;
        random_swap(i, T);
        printf("After swap:\n");
        print_vector(config.atoms);     
    }

    void TEST_sweep(int n_sweep=1, double T=2000)
    {
        TEST_PRINT_LINE("TEST_sweep");
        printf("Number of sweeps=%d\n", n_sweep);
        printf("Before %d sweeps:\n", n_sweep);
        print_vector(config.atoms);
        for (int i=0; i<n_sweep; i++)
        {
            energy_total = cal_energy();
            accept_count = 0; 
            for (int j=0; j<n_atom; j++)
            {
                random_swap(j, T); // energy_total and accept_count will be updated
                //printf("swap the %d-th atom, accept_count=%d\n", j, accept_count);
            }
            energy_ave.push_back(energy_total/n_atom);
            accept_ratio_ave.push_back(float(accept_count)/n_atom);
        }
        printf("After %d sweeps\n", n_sweep);
        print_vector(config.atoms);
        printf("Recorded energy_ave for each sweep:\n");
        print_vector(energy_ave);
        printf("Recorded accept_ratio_ave for each sweep:\n");
        print_vector(accept_ratio_ave);

        cout << "mc.energy_total=" << energy_total << endl;
        cout << "mc.cal_energy()" << cal_energy() << endl;

        energy_ave.clear();
        accept_ratio_ave.clear();
    }   

    Observables cal_observables(double T, int n_measure=1000, int n_discard=1000)
    {
        int n_shell_SRO = 2; // number of shells of SRO to be recoreded

        // discard the first n_discard sweeps used to reach thermal equilibrium.
        energy_total = cal_energy();
        for (int i=0; i<n_discard; i++)
        {
            accept_count = 0;
            for (int j=0; j<n_atom; j++)
            {
                random_swap(j, T);
            }
        }

        energy_ave.clear();
        accept_ratio_ave.clear();
        energy_total = cal_energy();

        Observables observables(n_shell_SRO, n_pair); // initialized with zeros
        for (int i=0; i<n_measure; i++)
        {   
            accept_count = 0; 
            for (int j=0; j<n_atom; j++)
            {
                random_swap(j, T);
            }
            energy_ave.push_back(energy_total/n_atom);
            accept_ratio_ave.push_back(float(accept_count)/n_atom);

            vector< vector<int> > pair_vec = config.cal_pair(nb);
            for(int j=0; j<n_shell_SRO; j++)
            {
                for(int k=0; k<n_pair; k++)
                {
                    // Divid pair_vec by the total number of atoms in the shell to get SRO parameters.
                    observables.SRO_parameters[j][k] += 
                    float(pair_vec[j][k])/ float(n_atom * nb.n_atom_shell[j]);
                }
            }
        }
        
        observables.temperature = T;
        observables.specific_heat = n_atom * variance(energy_ave)/(kb*T*T);
        /*
        for (int j=0; j< energy_ave.size(); j++)
        {
          printf("energy_ave: %d  %10.8f\n", j, energy_ave[j]);
        }
        //printf("Observables energy_ave:\n");
        //print_vector(energy_ave);
        cout << "mean(energy_ave)=" << mean(energy_ave) << endl;
        cout << "variance=" << variance(energy_ave) << endl;
        //cout << "kb*T*T=" << kb*T*T << endl;
        */
        for(int j=0; j<n_shell_SRO; j++)
        {
            for(int k=0; k<n_pair; k++)
            {
                observables.SRO_parameters[j][k] /= n_measure;
            }
        }
        return observables;
    }

    double cal_energy()
    {
        double energy = 0.0;
        vector< vector<int> > pair_vec = config.cal_pair(nb);

        for (int i=0; i<n_shell; i++)
            for (int j=0; j<n_pair; j++)
            {
                energy += pair_vec[i][j]*para_EPI[i][j]/float(nb.n_atom_shell[i]);
            }
        return energy;
    }

private:
    int accept_count;
    double energy_total;

    void random_swap(int site, double T)
    {
        const vector<int> &nb_i = nb.get_neighbors()[site];
        int ind_rd = random_mt_int.getRandom(); 
        // randomly select an integer in [0, nb.n_atom_shell[0]-1]
        int site_nb = nb_i[ind_rd]; // randomly select a site。

        if (config.atoms[site] == config.atoms[site_nb])
        {
            return;
        }
        double energy_before = cal_energy_i(site) + cal_energy_i(site_nb);
//cout << "energy_before=" << energy_before << endl;
        // Switch the element of site i and j
        stwitch_ele(site, site_nb);
        double energy_after = cal_energy_i(site) + cal_energy_i(site_nb);
//cout << "energy_after=" << energy_after << endl;
        double energy_diff = energy_after - energy_before;
        // Note the double counting of AB bond is cancelled off by E_after - E_before

        double beta = 1./(kb*T);
        if (exp(-beta*energy_diff) > random_mt.getRandom())
        {
            accept_count += 1;
            energy_total += energy_diff;
        }
        else
        {
            stwitch_ele(site, site_nb); // switch back
        }
        return;
    }

    void stwitch_ele(int i, int j)
    {
        int tmp = config.atoms[i];
        config.atoms[i] = config.atoms[j];
        config.atoms[j] = tmp;
    }

    double cal_energy_i(int site)
    {
        vector< vector<int> > pair_vec;
        // Initialized pair_vec[n_shell][n_pair] with zeros. 
        for (int i=0; i<(n_shell); i++)
        {   
            vector<int> tmp_vec(n_pair, 0);
            pair_vec.push_back(tmp_vec);
        }

        const vector<int> &shell_vec = nb.get_shell_list();
        int ind_i = config.atoms[site];
        const vector<int> &neighbor_i = nb.get_neighbors()[site];
        int n_ele = config.n_ele;

        for (int j=0; j<neighbor_i.size(); j++)
        {   
            int ind_j = config.atoms[neighbor_i[j]]; 
            int index = get_pair_index(ind_i, ind_j, n_ele);
            int ind_shell = shell_vec[j]; 
            pair_vec[ind_shell][index] +=1;
        }
        double energy = 0.0;

        for (int i=0; i<n_shell; i++)
            for (int j=0; j<n_pair; j++)
            {
                energy += 2*pair_vec[i][j]*para_EPI[i][j]/nb.n_atom_shell[i];
            }
        return energy;
    }
};

void TEST_configuration()
{
    TEST_PRINT_LINE("TEST_configuration");
    vector<string> elements{"W", "Mo", "Nb", "Ta"};
    int n_atom = 16;
    vector<int> n_atom_vec{4, 4, 5, 3};
    Configuration config_1(elements, n_atom);
    printf("config_1.atoms:\n");
    print_vector(config_1.elements);
    printf("n_atom = %d\n", config_1.n_atom);
    printf("n_pair = %d\n", config_1.n_pair);
    print_vector(config_1.atoms);

    Configuration config_2(elements, n_atom_vec);
    printf("config_2.atoms:\n");
    print_vector(config_2.elements);
    print_vector(n_atom_vec);
    printf("n_atom = %d\n", config_2.n_atom);
    printf("n_pair = %d\n", config_1.n_pair);
    print_vector(config_2.atoms);
}

void TEST_calpair(const string & filename)
{
    TEST_PRINT_LINE("TEST_cal_pair");
    Neighbor nb(filename);
    vector<string> elements{"W", "Mo", "Nb", "Ta"};
    int n_atom = nb_g.get_neighbors().size();
    Configuration config(elements, n_atom);
    vector< vector<int> > result = config.cal_pair(nb);
    printf("nb.n_shell = %d\n", nb.n_shell);
    printf("nb.n_pair = %d\n", config.n_pair);
    for (auto x: result)
    {
        print_vector(x);
    }
}

void TEST_MonteCarlo(vector<string> elements=vector<string> {"E1", "E2", "E3", "E4"}, int n_config_save=1)
{
    int n_atom = nb_g.get_neighbors().size();
    Configuration config_1(elements, n_atom);

    int n_ele = elements.size();
    int n_pair = n_ele*(n_ele + 1)/2;
    int n_shell = nb_g.n_shell;

    vector< vector<string> > data = read_csv(EPI_file, ','); 
    // EPI_file is a global variable
    vector<double> EPI_raw;
    for (auto x: data)
    {
        for (auto y: x)
        {
            EPI_raw.push_back(atof(y.c_str()));
        }
    }
    cout << "EPI_raw.size()=" << EPI_raw.size() << endl;
    cout << "n_shell=" << n_shell <<  " " << "n_pair="  << n_pair << endl;
    TEST_ASSERT(EPI_raw.size() >= n_shell*n_pair, "EPI_raw.size < n_shell*n_pair");

    vector< vector<double> > para_EPI;

    for (int i=0; i<n_shell; i++)
    {
        vector<double> tmp;
        for (int j=0; j<n_pair; j++)
        {
            tmp.push_back(EPI_raw[i*n_pair+j]);
        }       
        para_EPI.push_back(tmp);
    }

    //cout << "n_shell=" << n_shell << endl;
    //cout << "n_pair=" << n_pair << endl;

    MonteCarlo mc(config_1, nb_g, para_EPI);
    

#ifdef TEST_ALL
    mc.TEST_random_swap();
    mc.TEST_sweep(1000, 500);
#endif
    Observables obs;
    double T_max = 2000.0;
    double T_min = 50;
    int n_T = 50;
    double delta_T = (T_max - T_min)/(n_T -1);

    // dynamical memory allocation
    double *T_record = new double [n_T]();
    double *CV_record = new double [n_T]();
    double *SRO_record = new double [n_T*n_pair]();
    double *T_record_g = new double [n_T]();
    double *CV_record_g = new double [n_T]();
    double *SRO_record_g = new double [n_T*n_pair]();
    
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int id;
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    //printf("id=%d\n", id);
    //printf("np=%d\n", np);
    if (id==0)
    {
        TEST_PRINT_LINE("TEST_MonteCarlo");
        cout << "para_EPI" << endl;
        for (auto x: para_EPI)
        {
            print_vector(x);
        }
        cout << "n_atom=" << n_atom << endl;
        cout << "n_ele=" << n_ele << endl;
        cout << "n_pair=" << n_pair << endl;
        cout << "n_shell=" << n_shell << endl;
        cout << "kb=" << mc.kb << endl;
        cout << "n_T=" << n_T << endl;
    }
    // test Monte Carlo simulation at a fixed temperature
    // double T_test = 1000;
    // Observables obs1 = mc.cal_observables(T_test, 1000);
    // printf("obs.print_info at T= %f\n", T_test);
    // obs1.print_info();
    //printf("id = %d", id);
    // test Monte Carlo simulation at different temperatures
    int n_measure = 40000;
    int n_discard = 10000;
    for (int i=0; i<n_T; i++)
    {
        if (i%np == id)
        {
            double T_i = T_max - i*delta_T;
            cout << "T=" << T_i << endl;
            obs = mc.cal_observables(T_i, n_measure, n_discard);
            T_record[i] = obs.temperature;
            CV_record[i] = obs.specific_heat;
            for (int j=0; j<n_pair; j++)
            {
                SRO_record[i*n_pair+j] = obs.SRO_parameters[0][j]; 
                // store the nearest neighbor SRO
            }
        }
        else
        {
            continue;
        }
    }
    MPI_Reduce(T_record, T_record_g, n_T, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(CV_record, CV_record_g, n_T, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(SRO_record, SRO_record_g, n_T*n_pair, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
        printf("T_record_g:\n");
        for (int i=0; i<n_T; i++){cout << T_record_g[i] << endl;}
        printf("CV_record_g:\n");
        for (int i=0; i<n_T; i++){cout << CV_record_g[i] << endl;}
        printf("SRO_record_g:\n");
        for (int i=0; i<n_T*n_pair; i++){cout << SRO_record_g[i] << endl;}
    }

    // delete dynamically allocated memory
    delete [] T_record;
    delete [] CV_record;
    delete [] SRO_record;
    delete [] T_record_g;
    delete [] CV_record_g;
    delete [] SRO_record_g;
    MPI_Finalize();
}

struct Inputs
{
public:
    int n_measure = 40000;
    int n_discard = 10000;
    int n_skip = 100; // For taking snapshot, skip n_skip MC sweeps between each saved configurations
    double T_max = 2000.0;
    double T_min = 50;
    int n_T = 50;
    int n_config_save = 1;

    Inputs(){;}

    Inputs(string file_input)
    {
        ifstream infile;
        infile.open(file_input);
        if (!infile.is_open()) throw std::runtime_error("Could not open file");
        string line, word;
        vector<string> lines_record;
        if(infile.good())
        {
            while(getline(infile, line))
            {
                if(line[0] == '#')
                {
                    continue;
                }
                else
                {
                    lines_record.push_back(line);
                }
            }
        }
        infile.close();

        n_measure = std::stoi(lines_record[0]);
        n_discard = std::stoi(lines_record[1]);
        n_skip    = std::stoi(lines_record[2]);
        T_max     = std::stod(lines_record[3]);
        T_min     = std::stod(lines_record[4]);
        n_T       = std::stoi(lines_record[5]);
        n_config_save = std::stoi(lines_record[6]);
    }

    void print_inputs()
    {
        cout << "n_measure = " << n_measure << endl;
        cout << "n_discard = " << n_discard << endl;
        cout << "n_skip = " << n_skip << endl;
        cout << "T_max = " << T_max << endl;
        cout << "T_min = " << T_min << endl;
        cout << "n_T = " << n_T << endl;
        cout << "n_config_save = " << n_config_save << endl;
    }
};

Inputs input_g("input.dat");

void MonteCarloSimulation()
{
    //*********************************************************
    // set up the input parameters
    //*********************************************************
    int n_measure = input_g.n_measure;
    int n_discard = input_g.n_discard;
    int n_skip = input_g.n_skip; // For taking snapshot, skip n_skip MC sweeps between each saved configurations
    double T_max = input_g.T_max;
    double T_min = input_g.T_min;
    int n_T = input_g.n_T;
    int n_config_save = input_g.n_config_save;
    int n_atom = nb_g.get_neighbors().size();
    //*********************************************************
    // Read positions.xyz, and use it to initialize a Configuration.
    // Note Config_show is different from Configurations: The
    // purpose of Config_show is to record the lattice structures,
    // while Configuration is just a list of atoms, no lattice info
    //*********************************************************
    Config_show config;
    config.read_position("positions.xyz");
    printf("config.n_atom= %i\n", config.n_atom);
    cout << "n_atom =" << n_atom << endl;
    TEST_ASSERT(config.n_atom == n_atom, "nb_g.get_neighbors().size() != config.n_atom! in MonteCarloSimulation!");
    vector<string> elements;
    vector<int> n_atom_vec;
    for (auto const& w: config.get_atom_count())
    {
        elements.push_back(w.first);
        n_atom_vec.push_back(w.second);
    }
    Configuration config_1(elements, n_atom_vec); // initialize a Configuration
    int n_ele = elements.size();
    int n_pair = n_ele*(n_ele + 1)/2;
    int n_shell = nb_g.n_shell;
    //*********************************************************
    // read the EPI parameters and generate para_EPI
    //*********************************************************
    vector< vector<string> > data = read_csv(EPI_file, ','); 
    // EPI_file is a global variable
    vector<double> EPI_raw;
    for (auto x: data)
    {
        for (auto y: x)
        {
            EPI_raw.push_back(atof(y.c_str()));
        }
    }
    //cout << "EPI_raw.size()=" << EPI_raw.size() << endl;
    cout << "EPI_raw.size()=" << EPI_raw.size() << endl;
    cout << "n_shell=" << n_shell <<  " " << "n_pair="  << n_pair << endl;
    TEST_ASSERT(EPI_raw.size() >= n_shell*n_pair, "EPI_raw.size < n_shell*n_pair");
    vector< vector<double> > para_EPI;
    for (int i=0; i<n_shell; i++)
    {
        vector<double> tmp;
        for (int j=0; j<n_pair; j++)
        {
            tmp.push_back(EPI_raw[i*n_pair+j]);
        }       
        para_EPI.push_back(tmp);
    }

    //*********************************************************
    // Start the Monte Carlo simulation
    //*********************************************************
    MonteCarlo mc(config_1, nb_g, para_EPI);
//#ifdef TEST_ALL
    mc.TEST_random_swap();
    mc.TEST_sweep(10, 2000);
//#endif
    Observables obs;
    double delta_T = (T_max - T_min)/(n_T -1);

    // dynamical memory allocation to store observables
    double *T_record = new double [n_T]();
    double *CV_record = new double [n_T]();
    double *SRO_record = new double [n_T*n_pair]();
    double *T_record_g = new double [n_T]();
    double *CV_record_g = new double [n_T]();
    double *SRO_record_g = new double [n_T*n_pair]();
    
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int id;
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    //printf("id=%d\n", id);
    //printf("np=%d\n", np);
    if (id==0)
    {
        TEST_PRINT_LINE("Read input.dat");
        input_g.print_inputs();
        TEST_PRINT_LINE("MonteCarloSimulation");
        cout << "para_EPI" << endl;
        for (auto x: para_EPI)
        {
            print_vector(x);
        }
        cout << "n_atom=" << n_atom << endl;
        cout << "n_ele=" << n_ele << endl;
        cout << "n_pair=" << n_pair << endl;
        cout << "n_shell=" << n_shell << endl;
        cout << "kb=" << mc.kb << endl;
        cout << "n_T=" << n_T << endl;
    }

    // file_config for storing the renewed configuration, along with the lattice info
    TEST_ASSERT(n_skip >= 0, "n_skip need to >= 0!");
    for (int i=0; i<n_T; i++)
    {  
        if (i%np == id)
        {
            double T_i = T_max - i*delta_T;
            cout << "T=" << T_i << endl;
            string file_config = "configurations" + to_string(int(T_i)) + ".xyz";
            obs = mc.cal_observables(T_i, n_measure, n_discard);
            T_record[i] = obs.temperature;
            CV_record[i] = obs.specific_heat;
            for (int j=0; j<n_pair; j++)
            {
                SRO_record[i*n_pair+j] = obs.SRO_parameters[0][j]; 
                // store the nearest neighbor SRO
            }
            //*********************************************************
            // do extra simulations of n_skip*n_config_save at each temperature
            // and save n_config_save configurations, skipping n_skip sweeps in between.
            //*********************************************************
            for (int j=0; j<n_config_save; j++)
            {
                mc.cal_observables(T_i, 0, n_skip+1);   
                vector<int> atoms = mc.config.atoms;
                vector<string> elements = mc.config.elements;
                double energy_ave = mc.cal_energy()/mc.n_atom;
                if (j == 0)
                {
                    config.write_config(file_config, atoms, elements, false, energy_ave);
                }
                else
                {
                    config.write_config(file_config, atoms, elements, true, energy_ave);
                    //the above "true" arg means saving configurations using the append mode
                }
                // print_vector(elements);
                // print_vector(atoms);             
            }
        }
        else
        {
            continue;
        }
    }
    MPI_Reduce(T_record, T_record_g, n_T, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(CV_record, CV_record_g, n_T, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(SRO_record, SRO_record_g, n_T*n_pair, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
        printf("T_record_g:\n");
        for (int i=0; i<n_T; i++){cout << T_record_g[i] << endl;}
        printf("CV_record_g:\n");
        for (int i=0; i<n_T; i++){cout << CV_record_g[i] << endl;}
        //printf("SRO_record_g:\n");
        //for (int i=0; i<n_T*n_pair; i++){cout << SRO_record_g[i] << endl;}
    }

    // delete dynamically allocated memory
    delete [] T_record;
    delete [] CV_record;
    delete [] SRO_record;
    delete [] T_record_g;
    delete [] CV_record_g;
    delete [] SRO_record_g;
    MPI_Finalize();
}

int main(int argc, const char * argv[])
{
    Timer time; //record time
        
    #ifdef TEST_ALL
        TEST_getRandom();
        time.printTime();

        time.resetTime();
        TEST_shuffle();
        time.printTime();

        time.resetTime();
        TEST_configuration();
        time.printTime();
        time.resetTime();
        TEST_Neighbor("../Python/neighbors.dat");
        time.printTime();

        time.resetTime();
        TEST_calpair("../Python/neighbors.dat");
        time.printTime();
    
        time.resetTime();
        TEST_Config_show();
        time.printTime();
        time.resetTime();
        TEST_MonteCarlo();
        time.printTime();
    #endif
    time.resetTime();
    MonteCarloSimulation();
    time.printTime();

    return 0;
}
