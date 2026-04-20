#include "isr.h"

#include "rand.h"

#include <math.h>

namespace SimpleDY
{
    static constexpr double DELTA = 1.0e-4;
    static constexpr double B = 2.0;
    static constexpr double T_CUTOFF = 2; // Cutoff in GeV^2

    // TODO include running coupling
    static constexpr double ALPHA_S = 0.12;

    namespace
    {
        double __computeZNormalization(double zMin, double zMax)
        {
            return log((1 - zMin) / (1 - zMax));
        }

        double __computeTrialTFromSudakov(double sudakovTrial, double tMax, double zNormalization)
        {
            return tMax * exp(log(sudakovTrial) / B / zNormalization);
        }

        double __computeAcceptanceRatio(double tTrial, double zTrial, double x, int iParton, 
            const std::unique_ptr<LHAPDF::PDF>& pdf)
        {
            double pdfRatio = pdf->xfxQ2(iParton, x / zTrial, tTrial) / pdf->xfxQ2(iParton, x, tTrial);
            double prefactor = ALPHA_S / 2 / Math::PI * Physics::C_F / B;

            return prefactor * (1 + zTrial * zTrial) * pdfRatio;
        }

    } // namespace 

    Optional<ShowerEmission> generateFirstISR(const Event& event, const std::unique_ptr<LHAPDF::PDF>& pdf)
    {
        ShowerEmission emission;

        // Decide which leg; TODO: Weight by luminosity
        emission.leg = (rand() >= 0.5) ? 1 : -1;
        double x = emission.leg == 1 ? event.x1 : event.x2;

        double tMax = event.s;
        
        double zMin = x;
        double zMax = 1 - DELTA;
        double zNormalization = __computeZNormalization(zMin, zMax);

        while (tMax > T_CUTOFF)
        {
            // Sample the trial Sudakov factor and compute trial t from it
            double sudakovTrial = rand(0.0, 1.0);
            emission.t = __computeTrialTFromSudakov(sudakovTrial, tMax, zNormalization);

            // If the emission is unresolvable, do no emission
            if (emission.t <= T_CUTOFF)
                return {};

            // Sample another uniform value and compute z from it
            double u = rand(0.0, 1.0);
            emission.z = 1.0 - (1.0 - zMin) * exp(-u * zNormalization);

            // Sample phi uniformly
            emission.phi = rand(0.0, 2.0*Math::PI);

            // Decide whether to accept or reject the emission
            if (rand() < __computeAcceptanceRatio(emission.t, emission.z, x, emission.leg * event.iParton, pdf))
                return emission;
            else
                tMax = emission.t;
        }

        // No resolvable emission was generated
        return {};
    }
    
} // namespace SimpleDY
