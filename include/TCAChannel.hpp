#ifndef TCACHANNEL_HPP
#define TCACHANNEL_HPP

// Standard C++ includes

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

class TCAChannel : public TCAHistogramOwner
{
public:
    // Constructors
    TCAChannel();
    TCAChannel(const char* name, const char* title, const char* type);

    // Destructors
    virtual ~TCAChannel();

    // Getters
    inline size_t GetChannelID() const { return fChannelID; }

    // Setters

    // Methods
    virtual void PrintInfo() const;

protected:
    inline static size_t fgChannelIDCounter = 0; // Static counter to assign unique IDs
    const size_t fChannelID;                     // Unique ID for the channel
};

#endif // TCACHANNEL_HPP
