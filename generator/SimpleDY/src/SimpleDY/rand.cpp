#include "rand.h"

#include <random>

namespace SimpleDY
{
    double rand()
    {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<double> dist(0.0, 1.0);

        return dist(gen);
    }

    double rand(double min, double max)
    {
        return rand() * (max - min) + min;
    }
}
