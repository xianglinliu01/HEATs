#ifndef Random_mt_hpp
#define Random_mt_hpp

namespace Random
{
    class Random_mt
    {
    public:
        // random seed
        Random_mt():
        mt(rd()),dist(0.0, 1.0)
        {;}

        // specify seed
        Random_mt(unsigned int seed):mt(seed), dist(0.0, 1.0) {;}
        
        double getRandom() { return dist(mt); }
        
    private:
        std::random_device rd;
        std::mt19937 mt;
        std::uniform_real_distribution<double> dist;
    };

    class Random_mt_int
    {
    public:
        // random seed

        Random_mt_int(int min, int max):
        mt(rd()),dist(min, max)
        {;}

        // specify seed
        Random_mt_int(int min, int max, unsigned int seed):mt(seed), dist(min, max) {;}
        
        int getRandom() { return dist(mt); }
        
    private:
        std::random_device rd;
        std::mt19937 mt;
        std::uniform_int_distribution<int> dist;
    };
}
#endif
