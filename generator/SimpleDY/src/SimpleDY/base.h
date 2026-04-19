#pragma once

namespace SimpleDY
{
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double ALPHA = 1.0 / 137.035999084;
    static constexpr double NC = 3.0;
    static constexpr double GEV2_TO_MB = 0.389379338;
    
    struct FourMomentum
    {
        double e, x, y, z;
    };

} // namespace SimpleDY
