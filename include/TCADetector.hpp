#ifndef TCADETECTOR_HPP
#define TCADETECTOR_HPP

// Standard C++ includes

// ROOT includes

// Project includes
#include "TCAHistogramOwner.hpp"

class TCADetector : public TCAHistogramOwner
{
public:
    // Constructors
    TCADetector();
    TCADetector(const char* name);
    TCADetector(const char* name, const char* title);

    // Destructors
    virtual ~TCADetector();

    // Getters

    // Setters

    // Methods
    virtual void PrintInfo() const;

protected:
    inline static size_t fgDetectorIDCounter = 0; // Static counter to assign unique IDs
    size_t fDetectorID;                           // Unique ID for the detector
};

#endif // TCADETECTOR_HPP