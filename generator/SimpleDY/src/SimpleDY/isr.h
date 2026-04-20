#pragma once

#include "base.h"
#include "event.h"
#include "optional.h"

#include <memory>

#include <LHAPDF/LHAPDF.h>

namespace SimpleDY
{
    struct ShowerEmission
    {
        int leg;    // +1 or -1
        double t;   // Ordering variable = pT^2
        double z;   // Momentum fraction
        double phi; 
    };
    
    Optional<ShowerEmission> generateFirstISR(const Event& event, const std::unique_ptr<LHAPDF::PDF>& pdf);
    
} // namespace SimpleDY
