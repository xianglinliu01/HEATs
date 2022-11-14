#ifndef Write_config
#define Write_config

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <map>
#include "TEST_utility.hpp"
#include "Print_vector.hpp"
using namespace std;
using TEST_utility::TEST_ASSERT;
using TEST_utility::TEST_PRINT_LINE;
using Print_vector::print_vector;

struct Config_show
{
public:
    int n_atom;
    vector< string> atoms_record;
    vector< vector<string> > positions;

    void read_position(string file_name)
    {
        ifstream infile;
        string file_loc = "../Python/neighbor_BCC/" + file_name;
        infile.open(file_loc);
        if (!infile.is_open()) throw std::runtime_error("Could not open file");
        string line, word;
        getline(infile, line);
        n_atom = std::stoi(line);
        
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
                    stringstream ss (line);
                    vector<string> row;
                    ss >> word;
                    atoms_record.push_back(word); // record the atoms in the xyz template
                    while (ss >>  word) // record the positions of each atom
                    {
                        row.push_back(word);
                    }
                    positions.push_back(row); 
                }
            }
        }
        infile.close();
    }

    const map<string, int> get_atom_count(bool print = false)
    {
        map<string, int> count;
        string word;
        TEST_ASSERT(atoms_record.size() > 0, "length of atoms_record less than 1. Check whether atoms_record is initialized!");
        for (auto const& word: atoms_record)
        {
            count[word]++;
        }
        if (print)
        {
            for (auto const& w: count)
            {
                cout<< "count atoms in Config_show: " << w.first << ":" << w.second << endl;  
            }
        }
        return count;
    }

    // write new configurations (i.e. change atoms from atom_record to elements[atoms[i]])
    void write_config(string file_name, vector<int> atoms, vector<string> elements,  bool append=false, double energy_ave=0.0)
    {
        ofstream outfile;
        string file_loc = "../Output/Configurations/" + file_name;
        if (append)
        {
            outfile.open(file_loc, ios::out | ios::app);
        }
        else
        {
            outfile.open(file_loc, ios::out);
        }
        if (!outfile.is_open()) throw std::runtime_error("Could not open file");
        outfile << n_atom << endl;
        outfile << "# energy_ave = " << energy_ave << " Ry" << endl;
        TEST_ASSERT(atoms.size() == n_atom, "atoms.size() == n_atom not true!");
        for (int i=0; i<n_atom; i++)
        {
            outfile << elements[atoms[i]] << " " << positions[i][0] << " " << positions[i][1] << " " <<  positions[i][2] << endl;
        }
        outfile.close();
    }
};

void TEST_Config_show()
{
    TEST_PRINT_LINE("TEST_Config_show");
    Config_show config;
    config.read_position("positions.xyz");
    cout << "n_atom: " << config.n_atom << endl;
    vector<int> atoms(config.n_atom, 0);
    vector<string> elements(1, "Va");
    config.write_config("test.dat", atoms, elements);
    config.get_atom_count(true);
    config.write_config("test.dat", atoms, elements, true);
}
#endif
