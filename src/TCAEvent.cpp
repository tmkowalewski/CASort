// C++ standard includes

// ROOT includes

// Project includes
#include "TCAEvent.hpp"

TCAEvent::TCAEvent(TCAExperiment *experiment)
    : fExperiment(experiment)
{
}

TCAEvent::~TCAEvent()
{
    for (auto &ptr : fData)
    {
        delete ptr;
    }
}