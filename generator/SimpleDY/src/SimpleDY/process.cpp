#include "process.h"

#include "rand.h"
#include "born.h"
#include "file.h"

namespace SimpleDY
{
    static constexpr int nTrialEvents = 10e4;
    static constexpr int nAcceptedEvents = 10e4;
    static constexpr double securityFactor = 1.2;
    
    namespace
    {
        void _calcP(Event& event)
        {
            auto sin_th = std::sin(std::acos(event.cos_th));

            double e = event.m / 2;
            double pz = event.m / 2. * event.cos_th;

            event.p1.e = e * std::cosh(event.y) + pz * std::sinh(event.y);
            event.p1.x = event.m / 2. * sin_th * std::cos(event.phi);
            event.p1.y = event.m / 2. * sin_th * std::sin(event.phi);
            event.p1.z = e * std::sinh(event.y) + pz * std::cosh(event.y);

            event.p2.e = e * std::cosh(event.y) - pz * std::sinh(event.y);
            event.p2.x = - event.m / 2. * sin_th * std::cos(event.phi);
            event.p2.y = - event.m / 2. * sin_th * std::sin(event.phi);
            event.p2.z = e * std::sinh(event.y) - pz * std::cosh(event.y);
        }

        std::string _toString(const Event& event)
        {
            return std::to_string(event.m) + ", " 
                + std::to_string(event.s) + ", "
                + std::to_string(event.y) + ", "
                + std::to_string(event.cos_th);
        }

    } // namespace
    
    void Process::init(const std::string& pdfDataLocation, const std::string& pdfSet)
    {
        LHAPDF::setPaths(pdfDataLocation);
        m_pdfs = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfSet, 0));
    }

    void Process::run()
    {
        _determineMaxWeight();
        
        while (m_events.size() < nAcceptedEvents)
        {
            Event event = _sampleNextEventKinematics();
            double weight = _computeEventWeight(event);
        
            double r = rand(0, 1);
            if (r < weight / m_MaxWeight)
                m_events.push_back(event); // accept event
            
            m_nEventTrials++;
        }

        std::cout << "Event Trials: " << m_nEventTrials << std::endl;
        
        _computeTotalCrossSection();
        
        std::cout << "Total cross section: " << m_TotalCrossSection  * GEV2_TO_MB << " mb." << std::endl;
    }

    void Process::writeToFile(const std::string& filePath)
    {
        std::string fileContent;
            
        for (const Event& event: m_events)
            fileContent.append(_toString(event) + '\n');
        
        File file = File(filePath);
        file.write(fileContent);
    }

    Event Process::_sampleNextEventKinematics()
    {
        Event event;

        event.m = rand(m_mMin, m_mMax);
        event.s = event.m*event.m;
        
        event.y_max = std::log(m_sqrtS / event.m);
        event.y = rand(-event.y_max, event.y_max);

        event.cos_th = rand(-1, 1);
        event.phi = rand(0, 2*PI);

        event.x1 = (event.m / m_sqrtS) * std::exp(event.y);
        event.x2 = (event.m / m_sqrtS) * std::exp(-event.y);

        _calcP(event);

        return event;
    }

    double Process::_computeEventWeight(const Event& event)
    {
        double dsigma = computeDSigma(event, m_sqrtS, m_pdfs);
        
        // The inverse sampling factor
        double p_inv = 8.0*PI * (m_mMax-m_mMin) * event.y_max;            
        
        return dsigma * p_inv;
    }

    void Process::_determineMaxWeight()
    {
        double maxWeight = 0.0;
        
        for (int i = 0; i < nTrialEvents; i++)
        {   
            Event event = _sampleNextEventKinematics();
            double weight = _computeEventWeight(event);
            
            if (weight > maxWeight)
                maxWeight = weight;
        }

        m_MaxWeight = securityFactor * maxWeight;
    }

    void Process::_computeTotalCrossSection()
    {
        m_TotalCrossSection = m_MaxWeight * nAcceptedEvents / m_nEventTrials;
    }

} // namespace SimpleDY
