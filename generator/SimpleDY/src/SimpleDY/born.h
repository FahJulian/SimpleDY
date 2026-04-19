#pragma once

namespace SimpleDY
{
    struct FourMomentum
    {
        double e, x, y, z;
    };

    struct Event
    {
        double m, s, y, cos_th, phi, x1, x2, y_max;
        FourMomentum p1, p2;
        double weight;
    };

    void doStuff();
}
