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
#include "Write_config.hpp"//#define TEST_ALL

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
    Timer time;
    for (int i=0; i<total; i++)
    {
        double x = random_mt.getRandom();
        double y = random_mt.getRandom();
        if (x*x + y*y < 1.0)
        {
            count+=1;
        }
    }
    time.printTime();
    printf("Pi from Monte Carlo = %10.8f\n", 4*float(count)/total);
}

int main()
{
    TEST_getRandom();
    return 0;
}
