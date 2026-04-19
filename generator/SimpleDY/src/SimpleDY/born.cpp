#include "born.h"

#include "rand.h"
#include "file.h"

#include <math.h>

#include <iostream>

#include <LHAPDF/LHAPDF.h>

namespace SimpleDY
{   
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double ALPHA = 1.0 / 137.035999084;
    static constexpr double NC = 3.0;
    static constexpr double GEV2_TO_MB = 0.389379338;

    static const int sqrt_S = 8e3;
    static const int m_min = 20;
    static const int m_max = 200;

    void sampleKinematics(Event& event)
    {
        event.m = rand(m_min, m_max);
        event.s = event.m*event.m;
        
        event.y_max = std::log(sqrt_S / event.m);
        event.y = rand(-event.y_max, event.y_max);

        event.cos_th = rand(-1, 1);
        event.phi = rand(0, 2*PI);

        event.x1 = (event.m / sqrt_S) * std::exp(event.y);
        event.x2 = (event.m / sqrt_S) * std::exp(-event.y);
    }

    void calc_p(Event& event)
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

    double computeLuminosity(const Event& event, const std::unique_ptr<LHAPDF::PDF>& pdf)
    {
        double luminosity = 0;
        for (int iFlavour = 1; iFlavour <= 5; iFlavour++)
        {
            double f11 = pdf->xfxQ2(iFlavour, event.x1, event.s) / event.x1;
            double f12 = pdf->xfxQ2(-iFlavour, event.x2, event.s) / event.x2;
            double f21 = pdf->xfxQ2(-iFlavour, event.x1, event.s) / event.x1;
            double f22 = pdf->xfxQ2(iFlavour, event.x2, event.s) / event.x2;
            
            double qcharge_sq = (iFlavour == 2 || iFlavour == 4) ? 4.0/9.0 : 1.0/9.0;

            luminosity += qcharge_sq * (f11*f12 + f21*f22);
        }
    }

    double computeSigma(const std::vector<Event>& events)
    {
        double sigma = 0.0;
        for (const Event& event: events)
            sigma += event.weight;
        double prefactor = ALPHA * ALPHA / 2.0 / NC / sqrt_S / sqrt_S;
        sigma *= prefactor;
        sigma /= events.size();

        return sigma * GEV2_TO_MB;
    }

    std::string toString(const Event& event)
    {
        return std::to_string(event.m) + ", " 
            + std::to_string(event.s) + ", "
            + std::to_string(event.y) + ", "
            + std::to_string(event.cos_th) + ", "
            + std::to_string(event.weight);
    }

    void writeToFile(const std::vector<Event>& events, const std::string& filePath)
    {
        std::string fileContent;
            
        for (const Event& event: events)
            fileContent.append(toString(event) + '\n');
        
        File file = File(filePath);
        file.write(fileContent);
    }

    void doStuff()
    {
        LHAPDF::setPaths("/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/lhapdf");
        std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF("NNPDF40_lo_as_01180", 0));

        std::vector<Event> events;
        for (int i = 0; i < 100000; i++)
        {   
            Event event;

            sampleKinematics(event);
            calc_p(event);

            double luminosity = computeLuminosity(event, pdf);

            double born_factor = (1.0 + event.cos_th * event.cos_th) / event.m;
            double sampling_factor_inv = 8.0*PI * (m_max-m_min) * event.y_max;            

            event.weight = luminosity * born_factor * sampling_factor_inv;

            events.push_back(event);
        }

        std::cout << "Total cross section: " << computeSigma(events) << " mb." << std::endl;

        writeToFile(events, "/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/events/events.csv");
    }
}
