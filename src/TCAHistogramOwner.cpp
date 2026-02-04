// Standard C++ includes

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

TCAHistogramOwner::TCAHistogramOwner(const char* name, const char* title)
    : TNamed(name, title), fID(fgIDCounter++), fHistograms()
{
}