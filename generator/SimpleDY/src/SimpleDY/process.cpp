#include "process.h"

#include "isr.h"
#include "rand.h"
#include "born.h"
#include "file.h"

#include <iostream>

namespace SimpleDY
{
    static constexpr int N_TRIAL_EVENTS = 10e4;
    static constexpr int N_ACCEPTED_EVENTS = 10e3;
    static constexpr double SECURITY_FACTOR = 1.1;
    
    namespace
    {
        void __calcP(Event& event)
        {
            double sinTh = std::sqrt(1.0 - event.cosTh * event.cosTh);

            double e = event.m / 2;
            double pz = event.m / 2. * event.cosTh;

            event.p1.e = e * std::cosh(event.y) + pz * std::sinh(event.y);
            event.p1.x = event.m / 2. * sinTh * std::cos(event.phi);
            event.p1.y = event.m / 2. * sinTh * std::sin(event.phi);
            event.p1.z = e * std::sinh(event.y) + pz * std::cosh(event.y);

            event.p2.e = e * std::cosh(event.y) - pz * std::sinh(event.y);
            event.p2.x = - event.m / 2. * sinTh * std::cos(event.phi);
            event.p2.y = - event.m / 2. * sinTh * std::sin(event.phi);
            event.p2.z = e * std::sinh(event.y) - pz * std::cosh(event.y);
        }

        std::string __toString(const Event& event)
        {
            return std::to_string(event.m) + ", " 
                + std::to_string(event.s) + ", "
                + std::to_string(event.y) + ", "
                + std::to_string(event.cosTh);
        }

    } // namespace
    
    void Process::init(const std::string& pdfDataLocation, const std::string& pdfSet)
    {
        LHAPDF::setPaths(pdfDataLocation);
        m_pdfs = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfSet, 0));
    }

    void Process::run()
    {
        _clear();

        _determineMaxWeight();
        
        while (m_events.size() < N_ACCEPTED_EVENTS)
        {
            Event event = _sampleNextEventKinematics();
            BornEvent bornEvent = sampleBornEvent(event, m_sqrtS, m_pdfs);
            event.iParton = bornEvent.iParton;

            double weight = bornEvent.dSigma * _computePInverse(event);
        
            ASSERT(weight < m_MaxWeight, "Weight exceeds maximum weight.");

            double r = rand(0, 1);
            if (r < weight / m_MaxWeight)
            {
                auto emission = generateFirstISR(event, m_pdfs);
                if (emission.hasValue())
                {
                    std::cout << "Generated emission with t = " << emission.getValue().t << " GeV^2 "
                        << "and z = " << emission.getValue().z << "." << std::endl;
                }
                else 
                    std::cout << "No emission generated." << std::endl;

                m_events.push_back(event); // accept event
            }

            m_nEventTrials++;
        }

        std::cout << "Acceptance ratio: " << double(N_ACCEPTED_EVENTS) / m_nEventTrials << std::endl;
        
        _computeTotalCrossSection();
        
        std::cout << "Total cross section: " << m_TotalCrossSection  * Physics::GEV2_TO_MB << " mb." << std::endl;
    }

    void Process::writeToFile(const std::string& filePath)
    {
        std::string fileContent;
            
        for (const Event& event: m_events)
            fileContent.append(__toString(event) + '\n');
        
        File file = File(filePath);
        file.write(fileContent);
    }

    void Process::_clear()
    {
        m_events.clear();
        m_events.reserve(N_ACCEPTED_EVENTS);
        m_nEventTrials = 0;
        m_TotalCrossSection = 0.0;
        m_MaxWeight = 0.0;
    }

    Event Process::_sampleNextEventKinematics()
    {
        Event event;

        event.m = rand(m_mMin, m_mMax);
        event.s = event.m*event.m;
        
        event.yMax = std::log(m_sqrtS / event.m);
        event.y = rand(-event.yMax, event.yMax);

        event.phi = rand(0, 2*Math::PI);

        // Sample cos(theta) from p(c) = 3(1+c^2)/8
        double u = rand(0.0, 1.0);
        event.cosTh = 2.0 * std::sinh(std::asinh(4.0 * u - 2.0) / 3.0);

        event.x1 = (event.m / m_sqrtS) * std::exp(event.y);
        event.x2 = (event.m / m_sqrtS) * std::exp(-event.y);

        __calcP(event);

        return event;
    }

    // The inverse sampling factor
    double Process::_computePInverse(const Event& event)
    {
        double pCos = 3.0 / 8.0 * (1.0 + event.cosTh * event.cosTh);
        return 2.0 * Math::PI * (m_mMax-m_mMin) * 2.0 * event.yMax / pCos;            
    }

    void Process::_determineMaxWeight()
    {
        double maxWeight = 0.0;
        
        for (int i = 0; i < N_TRIAL_EVENTS; i++)
        {   
            Event event = _sampleNextEventKinematics();
            BornEvent bornEvent = sampleBornEvent(event, m_sqrtS, m_pdfs);

            double weight = bornEvent.dSigma * _computePInverse(event);
            
            if (weight > maxWeight)
                maxWeight = weight;
        }

        m_MaxWeight = SECURITY_FACTOR * maxWeight;
    }

    void Process::_computeTotalCrossSection()
    {
        m_TotalCrossSection = m_MaxWeight * N_ACCEPTED_EVENTS / m_nEventTrials;
    }

} // namespace SimpleDY
