#ifndef TCAEVENT_HPP
#define TCAEVENT_HPP

// Standard C++ includes
#include <array>
#include <cstddef>

// ROOT includes
#include <TTreeReaderArray.h>

// Project includes

// Forward declarations
class TCAExperiment;

class TCAEvent
{
public:
    static inline constexpr int kNModules = 4;

    enum FilterID
    {
        kAmplitude = 0,
        kChannelTime = 1,
        kPileUp = 2,
        kModuleTime = 3,
        kTriggerTime = 4,
        kIntLong = 5,
        kIntShort = 6,
        kNFilters
    };

    typedef std::array<TTreeReaderArray<double>*, kNModules * kNFilters> EventDataArray;

    // Constructors
    TCAEvent() = delete;
    TCAEvent(const TCAEvent&) = delete;
    TCAEvent(TCAExperiment* experiment);

    // Destructor
    ~TCAEvent();

    // Getters

    TCAExperiment* GetExperiment() const { return fExperiment; }

    // Setters

    void SetExperiment(TCAExperiment* experiment) { fExperiment = experiment; }

    // Operators
    inline double operator()(size_t moduleID, size_t filterID, size_t idx = 0) const
    {
        return (*fData[moduleID * kNFilters + filterID])[idx];
    }
    TTreeReaderArray<double>& operator[](size_t idx) const
    {
        return *fData[idx];
    }

private:
    TCAExperiment* fExperiment = nullptr;
    EventDataArray fData;
};

#endif // TCAEVENT_HPP