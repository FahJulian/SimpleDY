#include "event.h"

#include "SimpleDY/process.h"
#include "SimpleDY/born_event.h"

#include <math.h>

namespace SimpleDY
{
    namespace 
    {
        static constexpr double __ALLOWED_KINEMATIC_MISMATCH = 1.0e-10;

    } // namespace

    void Event::reconstructMomenta()
    {
        double pT = std::sqrt(m_emission.getT());
        double mT = std::sqrt(m_bornEvent.getS() + m_emission.getT());

        double x1PreEm = m_bornEvent.getX1() / (m_emission.getLeg() == 1 ? m_emission.getZ() : 1.0);
        double x2PreEm = m_bornEvent.getX2() / (m_emission.getLeg() == 1 ? 1.0 : m_emission.getZ());

        auto yBoson = _calculateBosonRapidity(x1PreEm, x2PreEm, mT);
        if (!yBoson.has_value())
        {
            m_emission.reject();
            reconstructMomenta();
            return;
        }

        m_pBoson = {
            mT * std::cosh(yBoson.value()),
            -pT * std::cos(m_emission.getPhi()),
            -pT * std::sin(m_emission.getPhi()),
            mT * std::sinh(yBoson.value())
        };

        m_p1In = {
            0.5 * x1PreEm * m_process.getSqrtS(),
            0.0,
            0.0,
            0.5 * x1PreEm * m_process.getSqrtS()
        };

        m_p2In = {
            0.5 * x2PreEm * m_process.getSqrtS(),
            0.0,
            0.0,
            -0.5 * x2PreEm * m_process.getSqrtS()
        };

        double p = std::sqrt(m_bornEvent.getS()) / 2.0;
        double sinTh = std::sqrt(1.0 - m_bornEvent.getCosTh() * m_bornEvent.getCosTh());

        FourVector p1Rest = {
            p,
            p * sinTh * std::cos(m_bornEvent.getPhi()),
            p * sinTh * std::sin(m_bornEvent.getPhi()),
            p * m_bornEvent.getCosTh()
        };

        FourVector p2Rest = { p1Rest.e, -p1Rest.getThreeVec() };

        m_p1Out = p1Rest.boost(m_pBoson);
        m_p2Out = p2Rest.boost(m_pBoson);
        m_pGluon = m_p1In + m_p2In - m_pBoson;

        ASSERT(_hasValidKinematics(), "Invalid kinematics");
    }

    std::string Event::toString() const
    {
        FourVector sum = m_p1Out + m_p2Out;
        double pT = std::sqrt(sum.pX*sum.pX + sum.pY*sum.pY);

        return std::to_string(m_bornEvent.getMBoson()) + ", "
            + std::to_string(m_bornEvent.getCosTh()) + ", "
            + std::to_string(pT);
    }

    std::optional<double> Event::_calculateBosonRapidity(double x1PreEm, double x2PreEm, double mT) const
    {
        if (m_emission.isRejected())
            return m_bornEvent.getYBoson();

        double a = x1PreEm * m_process.getSqrtS();
        double b = x2PreEm * m_process.getSqrtS();

        double C = (a * b + m_bornEvent.getS()) / mT;
        double disc = C * C - 4.0 * a * b;
        
        if (disc < 0) return { };   // unphysical kinematics, reject the emission

        double sqrtDisc = std::sqrt(disc);

        double u1 = (C + sqrtDisc) / (2.0 * b);
        double u2 = (C - sqrtDisc) / (2.0 * b);

        double y1 = std::log(u1);
        double y2 = std::log(u2);

        // pick the solution closer to the underlying Born rapidity
        return (std::abs(y1 - m_bornEvent.getYBoson()) < std::abs(y2 - m_bornEvent.getYBoson())) ? y1 : y2;
    }

    bool Event::_hasValidKinematics() const
    {
        FourVector totalIn = m_p1In + m_p2In;
        FourVector totalOut = m_p1Out + m_p2Out + m_pGluon;
        double mismatch = (totalIn - totalOut) * (totalIn - totalOut) / m_bornEvent.getS();

        return std::abs(mismatch) < __ALLOWED_KINEMATIC_MISMATCH;
    }

} // namespace SimpleDY
