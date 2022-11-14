#ifndef Variance
#define Variance

#include<iostream>
#include<vector>
#include<string>

using std::cout;
using std::endl;
using std::vector;
using std::string;

namespace Variance
{
    double mean( const vector<double>& list )
    {
        double sum = 0;
        for (int i=0; i<list.size(); i++){ sum += list[i]; }
        return sum/double(list.size());
    }

    double variance(const vector<double>& list )
    {
        double var = 0;
        double average = mean(list);
        for (int i=0; i<list.size(); i++)
        {
            var += (list[i]-average)*(list[i]-average);
        }
        var = var/double(list.size());
        return var;
    }
}

#endif
