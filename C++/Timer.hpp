#ifndef Timer_h
#define Timer_h

#include<ctime>
#include<iostream>

using std::cout;
using std::endl;

//  A Timer class to record the calculation time.
class Timer
{
public:
    Timer(){c_start=std::clock();}
    
    double getTime()
    {
        c_end = std::clock();
        time_elapsed_s = double(c_end-c_start) / CLOCKS_PER_SEC;
        return time_elapsed_s;
    }
    
    void printTime()
    {
        double time = getTime();
        if (time < 1.)
        {
            cout << "CPU time used: " << 1000.*time << " ms\n";
        }
        else
        {
            cout << "CPU time used: " << time << " s\n";
        }
        
    }
    
    void resetTime()
    {
        c_start = std::clock();
    }
    
private:
    std::clock_t c_start;
    std::clock_t c_end;
    double time_elapsed_s;
};
#endif
