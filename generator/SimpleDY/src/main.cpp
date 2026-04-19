#include "./SimpleDY/process.h"

constexpr double sqrtS = 8.0e3; // GeV
constexpr double mMin = 20; // GeV
constexpr double mMax = 120; // GeV
const std::string exportFilePath = "/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/events/events.csv";
const std::string pdfDataLocation = "/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/lhapdf";
const std::string pdfSet = "NNPDF40_lo_as_01180";

int main()
{
    SimpleDY::Process process(sqrtS, mMin, mMax);

    process.init(pdfDataLocation, pdfSet);
    process.run();
    process.writeToFile(exportFilePath);
}
