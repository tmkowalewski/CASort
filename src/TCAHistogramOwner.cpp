// Standard C++ includes

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

TCAHistogramOwner::TCAHistogramOwner(const char* name, const char* title)
    : TNamed(name, title), fOwnerID(fgOwnerIDCounter++), fHistograms()
{
}

TCAHistogramOwner::~TCAHistogramOwner()
{
}

// std::vector<std::shared_ptr<TH1>> TCAHistogramOwner::CreateThreadLocalPtrs()
// {
//     std::vector<std::shared_ptr<TH1>> ptrs;
//     // Unfinished
//     return ptrs;
// }

// void TCAHistogramOwner::PrintInfo() const
// {
//     /// Unfinished
// }