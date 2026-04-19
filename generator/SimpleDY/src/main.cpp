#include "./SimpleDY/process.h"

constexpr double sqrtS = 8.0e3; // GeV
constexpr double mMin = 20; // GeV
constexpr double mMax = 120; // GeV

int main()
{
    SimpleDY::Process process(sqrtS, mMin, mMax);

    process.init();
    process.run();
    process.writeToFile("/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/events/events.csv");
}
