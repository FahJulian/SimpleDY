#pragma once

namespace SimpleDY
{
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double ALPHA = 1.0 / 137.035999084;
    static constexpr double NC = 3.0;
    static constexpr double GEV2_TO_MB = 0.389379338;
    static constexpr double s_W_sq  = 0.23126;
    static constexpr double c_W_sq  = 1.0 - s_W_sq;
    static constexpr double m_Z     = 91.1876;   
    static constexpr double gamma_Z = 2.4952;    
    static constexpr double kappa   = 1.0 / (4.0 * s_W_sq * c_W_sq);
    
    struct FourMomentum
    {
        double e, x, y, z;
    };

} // namespace SimpleDY
