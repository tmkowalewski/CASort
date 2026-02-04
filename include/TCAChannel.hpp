#ifndef TCACHANNEL_HPP
#define TCACHANNEL_HPP

// Standard C++ includes

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

class TCAChannel : public TCAHistogramOwner
{
public:
    // Constructor
    TCAChannel();
    TCAChannel(const char* name, const char* title, const char* type);
    // Destructor
    virtual ~TCAChannel();

    // Methods
    virtual void PrintInfo() const override;
protected:

};

#endif // TCACHANNEL_HPP
