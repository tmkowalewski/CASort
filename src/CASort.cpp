// C++ Includes
#include <string>
#include <thread>

// ROOT Includes
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TTreeReaderArray.h>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>

// Project Includes

static const int kMaxThreads = std::thread::hardware_concurrency(); // Number of threads to use for processing, defaults to system max
static const constexpr int kDigitizerBins = 1 << 16;                // Number of bins in the digitizer (16-bit)
static const constexpr int kDigitizerChannels = 16;                 // Number of channels in digitizer
static const constexpr double kNsPerBin = 0.098;                    // Conversion factor from bin to nanoseconds

// Histograms Namespace
namespace Histograms
{
    auto cc_amp = ROOT::TThreadedObject<TH2D>("cc_amp", "Clover Cross Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_cht = ROOT::TThreadedObject<TH2D>("cc_cht", "Clover Cross Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_plu = ROOT::TThreadedObject<TH2D>("cc_plu", "Clover Cross Pile-Up;Pile-Up Multiplicity;Channel;Counts/Bin", 4, 0, 4, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_mdt = ROOT::TThreadedObject<TH1D>("cc_mdt", "Clover Cross Module Time;Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto cc_trt = ROOT::TThreadedObject<TH2D>("cc_trt", "Clover Cross Trigger Time;Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

}

// Main function
int main(int argc, char *argv[])
{
    // Introduction
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input filename> <output filename>" << std::endl;
        return 1;
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

    std::cout << "===== Welcome to CASort! =====" << std::endl;
    std::cout << "----- Current Configuration -----" << std::endl;
    std::cout << "Input file: " << input_filename << std::endl;
    std::cout << "Output file: " << output_filename << std::endl;
    std::cout << "Max Threads: " << kMaxThreads << std::endl;
    std::cout << "---------------------------------" << std::endl;

    // Set Max Threads
    ROOT::EnableImplicitMT(kMaxThreads);
    // ROOT::EnableThreadSafety();

    // Read the input ROOT file
    auto infile = TFile::Open(input_filename.c_str());
    if (!infile || infile->IsZombie())
    {
        std::cerr << "Error opening input file" << std::endl;
        return 1;
    }
    std::cout << "Opened file: " << input_filename << std::endl;

    // Create a TTreeReaderMT to read the TTree
    std::cout << "Processing events..." << std::endl;
    ROOT::TTreeProcessorMT EventProcessor(input_filename.c_str(), "clover");

    // Fill Function
    auto fillHistograms = [&](TTreeReader &eventReader)
    {
        TTreeReaderArray<double> cc_amp_val(eventReader, "clover_cross.amplitude");
        TTreeReaderArray<double> cc_cht_val(eventReader, "clover_cross.channel_time");
        TTreeReaderArray<double> cc_mdt_val(eventReader, "clover_cross.module_timestamp");
        TTreeReaderArray<double> cc_plu_val(eventReader, "clover_cross.pileup");
        TTreeReaderArray<double> cc_trt_val(eventReader, "clover_cross.trigger_time");

        TTreeReaderArray<double> cb_amp_val(eventReader, "clover_back.amplitude");
        TTreeReaderArray<double> cb_cht_val(eventReader, "clover_back.channel_time");
        TTreeReaderArray<double> cb_mdt_val(eventReader, "clover_back.module_timestamp");
        TTreeReaderArray<double> cb_plu_val(eventReader, "clover_back.pileup");
        TTreeReaderArray<double> cb_trt_val(eventReader, "clover_back.trigger_time");

        TTreeReaderArray<double> ps_amp_val(eventReader, "pos_sig.amplitude");
        TTreeReaderArray<double> ps_cht_val(eventReader, "pos_sig.channel_time");
        TTreeReaderArray<double> ps_mdt_val(eventReader, "pos_sig.module_timestamp");
        TTreeReaderArray<double> ps_plu_val(eventReader, "pos_sig.pileup");
        TTreeReaderArray<double> ps_trt_val(eventReader, "pos_sig.trigger_time");

        TTreeReaderArray<double> ce_inl_val(eventReader, "cebr_all.integration_long");
        TTreeReaderArray<double> ce_cht_val(eventReader, "cebr_all.channel_time");
        TTreeReaderArray<double> ce_mdt_val(eventReader, "cebr_all.module_timestamp");
        TTreeReaderArray<double> ce_ins_val(eventReader, "cebr_all.integration_short");
        TTreeReaderArray<double> ce_trt_val(eventReader, "cebr_all.trigger_time");

        auto cc_amp = Histograms::cc_amp.Get();
        auto cc_cht = Histograms::cc_cht.Get();
        auto cc_plu = Histograms::cc_plu.Get();
        auto cc_trt = Histograms::cc_trt.Get();
        auto cc_mdt = Histograms::cc_mdt.Get();

        // Loop over the entries in the tree
        while (eventReader.Next())
        {

            // Module Time
            cc_mdt->Fill(cc_mdt_val[0] * kNsPerBin);

            // Trigger Times
            cc_trt->Fill(cc_trt_val[0] * kNsPerBin, 0);
            cc_trt->Fill(cc_trt_val[1] * kNsPerBin, 1);
            // Main Loop

            // Detector Loop
            for (size_t det = 0; det < 4; det++)
            {
                // Crystal Loop
                for (size_t xtal = 0; xtal < 4; xtal++)
                {
                    auto ch = det * 4 + xtal; // Channel number 0-15

                    // Raw Histograms
                    cc_amp->Fill(cc_amp_val[ch], ch);
                    cc_cht->Fill(cc_cht_val[ch], ch);
                    cc_plu->Fill(cc_plu_val[ch], ch);
                }
            }
        }
    };

    // Loop over the entries in the TTree and fill the histograms appropriately
    EventProcessor.Process(fillHistograms);

    // Save the histograms to a new ROOT file
    auto outfile = std::make_unique<TFile>(output_filename.c_str(), "RECREATE");
    if (!outfile || outfile->IsZombie())
    {
        std::cerr << "Error creating output file" << std::endl;
        return 1;
    }

    Histograms::cc_amp.Merge()->Write();
    Histograms::cc_cht.Merge()->Write();
    Histograms::cc_plu.Merge()->Write();
    Histograms::cc_trt.Merge()->Write();
    Histograms::cc_mdt.Merge()->Write();

    std::cout << "Saved histograms to file: " << outfile->GetName() << std::endl;

    outfile->Close();

    std::cout << "Done!" << std::endl;

    return 0;
}