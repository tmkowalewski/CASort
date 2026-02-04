#ifndef TCAHISTOGRAMOWNER_HPP
#define TCAHISTOGRAMOWNER_HPP

// Standard C++ includes

// ROOT includes
#include <TNamed.h>
#include <TObjArray.h>

// Project includes

class TCAHistogramOwner : public TNamed
{
public:
    // Constructor
    TCAHistogramOwner(const char* name, const char* title);
    // Destructor
    virtual ~TCAHistogramOwner();

    // Methods
    inline static size_t GetOwnerCount() { return fgIDCounter; }
    inline size_t GetID() const { return fID; }



    virtual void PrintInfo() const;

private:
    inline static size_t fgIDCounter = 0;            // Static counter to assign unique

    const size_t fID;                                // Unique ID for the histogram owner
    TObjArray fHistograms;                     // Array of histograms owned by this owner



};

#endif // TCAHISTOGRAMOWNER_HPP
