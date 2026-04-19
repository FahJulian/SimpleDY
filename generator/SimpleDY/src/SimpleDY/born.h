#pragma once

#include "base.h"
#include "event.h"

#include <memory>
#include <vector>

#include <LHAPDF/LHAPDF.h>

namespace SimpleDY
{
    double computeDSigma(const Event& event, double sqrtS, const std::unique_ptr<LHAPDF::PDF>& pdf);
    double computeDSigma(const std::vector<Event>& events);

} // namespace SimpleDY
