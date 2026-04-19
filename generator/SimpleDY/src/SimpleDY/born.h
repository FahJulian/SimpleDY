#pragma once

#include "base.h"
#include "event.h"

#include <memory>
#include <vector>

#include <LHAPDF/LHAPDF.h>

namespace SimpleDY
{
    double computeBornKernel(const Event& event, const std::unique_ptr<LHAPDF::PDF>& pdf);
    double computeSigma(const std::vector<Event>& events);

} // namespace SimpleDY
