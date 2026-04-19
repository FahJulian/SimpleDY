#pragma once

namespace SimpleDY
{
    struct Event
    {
        double m, s, y, cos_th, phi, x1, x2;
        double e1, p1x, p1y, p1z;
        double e2, p2x, p2y, p2z;
        double weight;
    };

    void doStuff();
}
