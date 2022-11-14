#ifndef Read_csv_h
#define Read_csv_h

#include <vector>
#include <string>
#include <fstream>
#include <sstream> 
#include <stdexcept>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
using std::string;

vector< vector<string> > read_csv(string filename, char delim=',')
{
    vector< vector<string> > result;
    std::ifstream myFile(filename);
    if (!myFile.is_open()) throw std::runtime_error("Could not open file");

    string line, word;
    if(myFile.good())
    {
        while(std::getline(myFile, line))
        {
            std::stringstream ss(line);
            vector<string> row;
            while (std::getline(ss, word, delim))
            {
                row.push_back(word);
            }
            result.push_back(row);
        }
    }

    myFile.close();
    return result;
}
#endif
