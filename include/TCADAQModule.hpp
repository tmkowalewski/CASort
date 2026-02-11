#ifndef TCADAQMODULE_HPP
#define TCADAQMODULE_HPP

// Standard C++ includes
#include <vector>

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

// Forward declarations
class TCADetector;

class TCADAQModule : public TCAHistogramOwner
{
public:
    // Constructor
    TCADAQModule();
    TCADAQModule(const char* name, const char* title, const char* type);

    // Destructor
    virtual ~TCADAQModule();

    // Getters
    inline size_t GetModuleID() const { return fModuleID; }
    inline const char* GetType() const { return fType.c_str(); }

    // Setters

    // Methods
    virtual void PrintInfo() const;

protected:
    inline static size_t fgModuleIDCounter = 0; // Static counter to assign unique module IDs

    std::string fType;                   // Type of DAQ module (e.g. "MDPP16SCP")
    size_t fModuleID;                    // Unique ID for this DAQ module
    const size_t fChannelCount;          // Number of channels in this DAQ module
    std::vector<TCADetector> fDetectors; // Detectors associated with this DAQ module
};

#endif // TCADAQMODULE_HPP