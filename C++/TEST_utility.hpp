#ifndef TEST_utility_h
#define TEST_utility_h

#include<cstdio>
#include<iostream>
#include <stdexcept>

using std::cout;
using std::endl;

namespace TEST_utility
{
    void TEST_PRINT_LINE(std::string title)
    {
        printf("============================================\n");
        cout << title << endl;
        printf("============================================\n");
    }

    void TEST_ASSERT(bool condition, const char* message)
    {
        try
        {
            if (!condition)
            {
                throw std::runtime_error(message);
            }
        }
        catch(std::runtime_error& exception)
        {
            std::cerr << "Error: " << exception.what() << endl;
            abort();
        }
    }
}
#endif
