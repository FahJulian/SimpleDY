#include "process.h"

#include "rand.h"
#include "born.h"
#include "file.h"

namespace SimpleDY
{
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
                + std::to_string(event.cos_th) + ", "
                + std::to_string(event.weight);
        }

    } // namespace
    
    void Process::init()
    {
        LHAPDF::setPaths("/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/lhapdf");
        m_pdfs = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF("NNPDF40_lo_as_01180", 0));
    }

    void Process::run()
    {
        for (int i = 0; i < 1000000; i++)
        {   
            Event event = _sampleNextEventKinematics();
            
            double prefactor = ALPHA * ALPHA / 2.0 / NC / m_sqrtS / m_sqrtS / event.m;
            double kernel = computeBornKernel(event, m_pdfs);
            
            // dσ / dM dy dcosθ dφ
            double dsigma = prefactor * kernel;
            
            // The inverse sampling factor
            double p_inv = 8.0*PI * (m_mMax-m_mMin) * event.y_max;            
            
            event.weight = dsigma * p_inv;
            
            m_events.push_back(event);
        }
        
        std::cout << "Total cross section: " << computeSigma(m_events) << " mb." << std::endl;
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

} // namespace SimpleDY
