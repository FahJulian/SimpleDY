#pragma once

#include "base.h"

namespace SimpleDY
{
    struct Event
    {
        double m, s, y, cos_th, phi, x1, x2, y_max;
        FourMomentum p1, p2;
        double weight;
    };
    
} // namespace SimpleDY
