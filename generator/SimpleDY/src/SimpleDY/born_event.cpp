#include "born_event.h"

#include "SimpleDY/rand.h"
#include "SimpleDY/process.h"

#include <math.h>

namespace SimpleDY
{
    namespace
    {
        // TODO: This still looks ugly
        std::tuple<double, double> __neutralCurrentCouplings(bool upType, double mSq)
        {
            // quark charge and axial and vector couplings
            double qQ = upType ? 2.0 / 3.0 : -1.0 / 3.0;
            double aQ = upType ? 0.5 : -0.5;
            double vQ = upType ? 0.5 - (4.0 / 3.0) * Physics::S_W_SQ : -0.5 + (2.0 / 3.0) * Physics::S_W_SQ;

            double dZ = (mSq - Physics::M_Z * Physics::M_Z) * (mSq - Physics::M_Z * Physics::M_Z) 
                + Physics::M_Z * Physics::M_Z * Physics::GAMMA_Z * Physics::GAMMA_Z;

            double ReChi   = Physics::KAPPA * mSq * (mSq - Physics::M_Z * Physics::M_Z) / dZ;
            double AbsChiSq = Physics::KAPPA * Physics::KAPPA * mSq * mSq / dZ;

            double hU = qQ * qQ
                    - 2.0 * qQ * Physics::V_L * vQ * ReChi
                    + (Physics::V_L * Physics::V_L + Physics::A_L * Physics::A_L) * (vQ * vQ + aQ * aQ) * AbsChiSq;

            double hF = -2.0 * qQ * Physics::A_L * aQ * ReChi
                    + 4.0 * Physics::V_L * Physics::A_L * vQ * aQ * AbsChiSq;

            return {hU, hF};
        }

    } // namespace

    void BornEvent::sampleKinematics()
    {
        // Sample the invariant mass of the event uniformly from the allowed range
        m_mBoson = rand(m_process.getMMin(), m_process.getMMax());

        // Sample boson rapidity uniformly from the kinematically allowed region (x_1,2 \leq 1)
        double yBosonMax = std::log(m_process.getSqrtS() / m_mBoson);
        m_yBoson = rand(-yBosonMax, yBosonMax);

        // Sample \cos(\theta) from the probability distr. p(c) = 3(1+c^2)/8
        double u = rand(0.0, 1.0);
        m_cosTh = 2.0 * std::sinh(std::asinh(4.0 * u - 2.0) / 3.0);

        // Sample \phi uniformly
        m_phi = rand(0, 2.0 * Math::PI);

        // Calculate x1 and x2 from the sampled variables
        m_x1 = m_mBoson / m_process.getSqrtS() * std::exp(m_yBoson);
        m_x2 = m_mBoson / m_process.getSqrtS() * std::exp(-m_yBoson);
    }

    void BornEvent::computeWeightAndSampleParton()
    {
        // For each quark flavour, compute the individual contribution to the cross section
        auto channels = _computePartonChannelContributions();

        m_dSigma = 0.0;
        for (auto [partonId, dSigma] : channels)
            m_dSigma += dSigma;

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, m_dSigma);
        for (auto [partonId, dSigma] : channels)
        {
            m_partonId = partonId;
            if (u < dSigma)
                break;

            u -= dSigma;
        }
    }

    double BornEvent::_computeInverseSamplingDensity()
    {
        double jacobianM = m_process.getMMax() - m_process.getMMin();
        double jacobianY = 2.0 * std::log(m_process.getSqrtS() / m_mBoson);
        double jacobianCosTh = 8.0 / 3.0 / (1.0 + m_cosTh * m_cosTh);
        double jacobianPhi = 2.0 * Math::PI;

        return jacobianM * jacobianY * jacobianCosTh * jacobianPhi;  
    }

    std::vector<std::tuple<int, double>> BornEvent::_computePartonChannelContributions()
    {
        double mSq = m_mBoson * m_mBoson;
        double physicsPrefactor = Physics::ALPHA * Physics::ALPHA / 2.0 / Physics::NC 
            / m_process.getSqrtS() / m_process.getSqrtS() / m_mBoson;
        double samplingPrefactor  = _computeInverseSamplingDensity();

        std::vector<std::tuple<int, double>> channels;
        channels.reserve(10);
        
        for (int partonId = 1; partonId <= 5; partonId++)
        {
            // The luminosity factors if the quark is on leg 1 and the antiquark on leg 2
            double f1  = m_process.getPdfs()->xfxQ2( partonId, m_x1, mSq) / m_x1;
            double fb1 = m_process.getPdfs()->xfxQ2(-partonId, m_x2, mSq) / m_x2;
            
            // The luminosity factors if the quark is on leg 2 and the antiquark on leg 1
            double f2  = m_process.getPdfs()->xfxQ2( partonId, m_x2, mSq) / m_x2;
            double fb2 = m_process.getPdfs()->xfxQ2(-partonId, m_x1, mSq) / m_x1;

            // Compute the couplings 
            bool upType = partonId % 2 == 0;
            auto [c1, c2] = __neutralCurrentCouplings(upType, mSq);

            // Compute the weight if the quark is in leg 1 or leg 2, respectively
            double weight1 = f1 * fb1 * (c1 * (1.0 + m_cosTh * m_cosTh) + 2.0 * c2 * m_cosTh);
            double weight2 = f2 * fb2 * (c1 * (1.0 + m_cosTh * m_cosTh) - 2.0 * c2 * m_cosTh);

            // Store the two contributions to the differential cross section with 
            // the respective id of the parton on leg 1
            channels.push_back({  partonId, samplingPrefactor * physicsPrefactor * weight1 });
            channels.push_back({ -partonId, samplingPrefactor * physicsPrefactor * weight2 });
        }       

        return channels;
    }

} // namespace SimpleDY
