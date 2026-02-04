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

    // Methods
    inline static size_t GetOwnerCount() { return fgIDCounter; }
    inline size_t GetID() const { return fID; }
    inline const TObjArray& GetHistograms() const { return fHistograms; }
    std::vector<std::shared_ptr<TH1>> GetThreadLocalHistograms();
    inline TH1D* GetHistogramAt(const size_t index) const { return static_cast<TH1D*>(fHistograms.At(index)); }

    template <typename T, typename... Args>
    void AddHistogram(Args&&... args)
    {
        auto hist = std::make_shared<T>(std::forward<Args>(args)...);
        fHistograms.Add(hist.get());
    }

    virtual void PrintInfo() const;

protected:
    inline static size_t fgIDCounter = 0;            // Static counter to assign unique

    const size_t fID;                                // Unique ID for the histogram owner
    TObjArray fHistograms;                     // Array of histograms owned by this owner
};

#endif // TCAHISTOGRAMOWNER_HPP
