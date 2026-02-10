#ifndef TCAHISTOGRAM_HPP
#define TCAHISTOGRAM_HPP

// Standard C++ includes

// ROOT includes
#include <ROOT/TThreadedObject.hxx>

// Project includes

// Forward declarations
class TCAEvent;

template <typename T>
class TCAHistogram : public TNamed
{
public:
    template <typename... Args>
    explicit TCAHistogram(Args&&... args)
        : fHistogram(ROOT::TThreadedObject<T>(std::forward<Args>(args)...))
    {
        fName = fHistogram.Get()->GetName();
        fTitle = fHistogram.Get()->GetTitle();
    }
    virtual ~TCAHistogram() = default;

    template <typename... Args>
    inline void Fill(Args&&... args) { fFillFunction(std::forward<Args>(args)...); }

    void SetFillFunction(const std::function<void(std::shared_ptr<T>, TCAEvent* event)>& func) { fFillFunction = func; }

    auto GetPtr() { return fHistogram.Get(); }
    auto GetRawPtr() { return fHistogram.Get().get(); }
    auto GetThreadLocalPtr() { return fHistogram.Get(); }
    auto Merge() { return fHistogram.Merge(); }
    auto Write() { return this->Merge()->Write(); }

protected:
    ROOT::TThreadedObject<T> fHistogram;
    std::function<void(std::shared_ptr<T>, TCAEvent* event)> fFillFunction = [this](std::shared_ptr<T> thisHist, TCAEvent* event) {}; // Default does nothing
};

#endif // TCAHISTOGRAM_HPP