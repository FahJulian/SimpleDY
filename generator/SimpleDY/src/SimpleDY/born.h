#pragma once

#include "base.h"
#include "event.h"

#include <memory>

#include <LHAPDF/LHAPDF.h>

namespace SimpleDY
{
    double computeDSigma(const Event& event, double sqrtS, const std::unique_ptr<LHAPDF::PDF>& pdf);

} // namespace SimpleDY
