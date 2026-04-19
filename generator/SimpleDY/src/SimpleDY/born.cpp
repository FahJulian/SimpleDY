#include "born.h"

#include <math.h>

namespace SimpleDY
{   
    namespace
    {
        struct _EWFactors
        {
            double H_U, H_F;
        };

        _EWFactors _neutralCurrentFactors(int iFlavour, double m_sq)
        {
            // charged lepton axial and vector couplings
            static constexpr double a_l = -0.5;
            static constexpr double v_l = -0.5 + 2.0 * s_W_sq;

            // quark charge and axial and vector couplings
            double q_q = (iFlavour % 2 == 0) ? 2.0 /3.0 : -1.0 / 3.0;
            double a_q = (iFlavour % 2 == 0) ? 0.5 : -0.5;
            double v_q = (iFlavour % 2 == 0) ? 0.5 - (4.0 / 3.0) * s_W_sq : -0.5 + (2.0 / 3.0) * s_W_sq;

            double D_Z = (m_sq - m_Z * m_Z) * (m_sq - m_Z * m_Z) + m_Z * m_Z * gamma_Z * gamma_Z;

            double Re_chi   = kappa * m_sq * (m_sq - m_Z * m_Z) / D_Z;
            double abs_chi_sq = kappa * kappa * m_sq * m_sq / D_Z;

            double H_U = q_q * q_q
                    - 2.0 * q_q * v_l * v_q * Re_chi
                    + (v_l * v_l + a_l * a_l) * (v_q * v_q + a_q * a_q) * abs_chi_sq;

            double H_F = -2.0 * q_q * a_l * a_q * Re_chi
                    + 4.0 * v_l * a_l * v_q * a_q * abs_chi_sq;

            return {H_U, H_F};
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

            _EWFactors ew = _neutralCurrentFactors(iFlavour, event.s);

            kernel += Lplus  * ew.H_U * (1.0 + event.cos_th * event.cos_th)
                    + Lminus * (2.0 * ew.H_F * event.cos_th);
        }

        double prefactor = ALPHA * ALPHA / 2.0 / NC / sqrtS / sqrtS / event.m;
        return prefactor * kernel;
    }

} // namespace SimpleDY
