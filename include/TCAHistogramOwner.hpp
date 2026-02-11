#ifndef TCAHISTOGRAMOWNER_HPP
#define TCAHISTOGRAMOWNER_HPP

// Standard C++ includes

// ROOT includes
#include <TH1.h>
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

    // Getters
    inline static size_t GetOwnerCount() { return fgOwnerIDCounter; }
    inline size_t GetOwnerID() const { return fOwnerID; }
    inline const TObjArray& GetHistograms() const { return fHistograms; }

    template <typename T>
    inline T* GetHistogramAt(const size_t index) const { return static_cast<T*>(fHistograms.At(index)); }

    // Setters

    // Methods
    // std::vector<std::shared_ptr<TH1>> CreateThreadLocalPtrs();

    template <typename T, typename... Args>
    void AddHistogram(Args&&... args)
    {
        auto hist = new T(std::forward<Args>(args)...);
        fHistograms.Add(hist);
    }

    // virtual void PrintInfo() const;

protected:
    inline static size_t fgOwnerIDCounter = 0; // Static counter to assign unique IDs

    const size_t fOwnerID; // Unique ID for the histogram owner
    TObjArray fHistograms; // Array of histograms owned by this owner
};

#endif // TCAHISTOGRAMOWNER_HPP
