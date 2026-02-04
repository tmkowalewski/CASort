#ifndef TCADAQMODULE_HPP
#define TCADAQMODULE_HPP

// Standard C++ includes

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

class TCADAQModule : public TCAHistogramOwner
{
public:
    // Constructor
    TCADAQModule();
    TCADAQModule(const char* name, const char* title, const char* type);
    // Destructor
    virtual ~TCADAQModule();

    // Methods
    virtual void PrintInfo() const override;

protected:
    std::string fType; // Type of DAQ module (e.g. "MDPP16SCP")
    const size_t fChannelCount; // Number of channels in this DAQ module

};

#endif // TCADAQMODULE_HPP