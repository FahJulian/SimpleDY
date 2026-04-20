#include "born.h"

#include <math.h>

namespace SimpleDY
{   
    namespace
    {
        struct __EWFactors
        {
            double hU, hF;
        };

        __EWFactors __neutralCurrentFactors(int iFlavour, double m_sq)
        {
            // quark charge and axial and vector couplings
            double qQ = (iFlavour % 2 == 0) ? 2.0 /3.0 : -1.0 / 3.0;
            double aQ = (iFlavour % 2 == 0) ? 0.5 : -0.5;
            double vQ = (iFlavour % 2 == 0) ? 0.5 - (4.0 / 3.0) * Physics::S_W_SQ : -0.5 + (2.0 / 3.0) * Physics::S_W_SQ;

            double dZ = (m_sq - Physics::M_Z * Physics::M_Z) * (m_sq - Physics::M_Z * Physics::M_Z) 
                + Physics::M_Z * Physics::M_Z * Physics::GAMMA_Z * Physics::GAMMA_Z;

            double ReChi   = Physics::KAPPA * m_sq * (m_sq - Physics::M_Z * Physics::M_Z) / dZ;
            double AbsChiSq = Physics::KAPPA * Physics::KAPPA * m_sq * m_sq / dZ;

            double hU = qQ * qQ
                    - 2.0 * qQ * Physics::V_L * vQ * ReChi
                    + (Physics::V_L * Physics::V_L + Physics::A_L * Physics::A_L) * (vQ * vQ + aQ * aQ) * AbsChiSq;

            double hF = -2.0 * qQ * Physics::A_L * aQ * ReChi
                    + 4.0 * Physics::V_L * Physics::A_L * vQ * aQ * AbsChiSq;

            return {hU, hF};
        }

    } // namespace

    double computeDSigma(const Event& event, double sqrtS, const std::unique_ptr<LHAPDF::PDF>& pdf)
    {
        double kernel = 0.0;

        for (int iFlavour = 1; iFlavour <= 5; iFlavour++)
        {
            double q1  = pdf->xfxQ2( iFlavour, event.x1, event.s) / event.x1;
            double qb2 = pdf->xfxQ2(-iFlavour, event.x2, event.s) / event.x2;
            double qb1 = pdf->xfxQ2(-iFlavour, event.x1, event.s) / event.x1;
            double q2  = pdf->xfxQ2( iFlavour, event.x2, event.s) / event.x2;

            double Lplus  = q1 * qb2 + qb1 * q2;
            double Lminus = q1 * qb2 - qb1 * q2;

            __EWFactors ew = __neutralCurrentFactors(iFlavour, event.s);

            kernel += Lplus  * ew.hU * (1.0 + event.cosTh * event.cosTh)
                    + Lminus * (2.0 * ew.hF * event.cosTh);
        }

        double prefactor = Physics::ALPHA * Physics::ALPHA / 2.0 / Physics::NC / sqrtS / sqrtS / event.m;
        return prefactor * kernel;
    }

} // namespace SimpleDY
