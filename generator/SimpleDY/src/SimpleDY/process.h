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

        void init();
        void run();
        void writeToFile(const std::string& filePath);

    private:
        Event _sampleNextEventKinematics();

    private:
        const double m_sqrtS;   
        const double m_mMin, m_mMax;

        std::unique_ptr<LHAPDF::PDF> m_pdfs;
        std::vector<Event> m_events;
    };  

} // namespace SimpleDY
