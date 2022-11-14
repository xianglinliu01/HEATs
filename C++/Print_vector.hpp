#ifndef print_vector_h
#define print_vector_h

#include<iostream>
#include<cstdio>
#include<vector>
#include<string>

using std::cout;
using std::endl;
using std::vector;
using std::string;

namespace Print_vector
{
    void print_vector(const vector<int> &v1)
    {
        for(int i=0; i<v1.size(); i++)
        {
            printf("%d ", v1[i]);
        }
        printf("\n");
    }

    void print_vector(const vector<double> &v1)
    {
        for(int i=0; i<v1.size(); i++)
        {
            printf("%16.8f", v1[i]);
        }
        printf("\n");
    }

    void print_vector(const vector<string> &v1)
    {
        for(int i=0; i<v1.size(); i++)
        {
            cout << v1[i] << " ";
        }
        cout << endl;
    }
}
#endif
