#pragma once

#include "base.h"
#include "event.h"

#include <memory>
#include <vector>
#include <string>

#include <LHAPDF/LHAPDF.h>

namespace SimpleDY
{
    class Process
    {
    public:
        Process() = delete; 

        Process(double sqrtS, double mMin, double mMax)
            : m_sqrtS(sqrtS), m_mMin(mMin), m_mMax(mMax)
        {
        }

        ~Process() = default;

        void init(const std::string& pdfDataLocation, const std::string& pdfSet);
        void run();
        void writeToFile(const std::string& filePath);

    private:
        void _clear();
        Event _sampleNextEventKinematics();
        double _computePInverse(const Event& event);
        void _determineMaxWeight();
        void _computeTotalCrossSection();

    private:
        const double m_sqrtS;   
        const double m_mMin, m_mMax;
        int m_nEventTrials = 0;
        double m_MaxWeight = 0.0;
        double m_TotalCrossSection = 0.0;

        std::unique_ptr<LHAPDF::PDF> m_pdfs;
        std::vector<Event> m_events;
    };  

} // namespace SimpleDY
