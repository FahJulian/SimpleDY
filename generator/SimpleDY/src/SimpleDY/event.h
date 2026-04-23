#pragma once

#include "SimpleDY/base.h"
#include "SimpleDY/emission.h"
#include "SimpleDY/born_event.h"
#include "SimpleDY/four_vector.h"

#include <string>
#include <optional>

namespace SimpleDY
{
    class Process;

    class Event
    {
    public:
        Event(const Process& process, const BornEvent& bornEvent, const Emission& emission)
            : m_process(process), m_bornEvent(bornEvent), m_emission(emission)
        {
        }

        void reconstructMomenta();
        std::string toString() const;

    private:
        std::optional<double> _calculateBosonRapidity(double x1PreEm, double x2PreEm, double mT) const;
        bool _hasValidKinematics() const;

    private:
        const Process& m_process;
        const BornEvent m_bornEvent;
        Emission m_emission;
        FourVector m_p1In;
        FourVector m_p2In;
        FourVector m_p1Out;
        FourVector m_p2Out;
        FourVector m_pBoson;
        FourVector m_pGluon;
    };
    
} // namespace SimpleDY
