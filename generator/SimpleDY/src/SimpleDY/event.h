#pragma once

#include "base.h"

namespace SimpleDY
{
    struct Event
    {
        double m, s, y, cosTh, phi, x1, x2, yMax;
        FourMomentum p1, p2;
    };
    
} // namespace SimpleDY
