#ifndef CACONFIGURATION_HPP
#define CACONFIGURATION_HPP

// C++ Includes
#include <algorithm>
#include <thread>

// ROOT Includes
#include <ROOT/RTaskArena.hxx>
#include <TROOT.h>

// Number of hardware threads to use in processing
const unsigned int kMaxThreads = std::thread::hardware_concurrency();

// Initialize ROOT's thread pool before any TThreadedObject is created
namespace
{
    struct ROOTThreadPoolInitializer
    {
        ROOTThreadPoolInitializer()
        {
            ROOT::EnableImplicitMT(kMaxThreads);
            ROOT::EnableThreadSafety();
        }
    };
    static ROOTThreadPoolInitializer rootThreadPoolInit;
}

// Debug Mode
#define DEBUG 1

// Which modules to process
#define PROCESS_CLOVER_CROSS true
#define PROCESS_CLOVER_BACK true
#define PROCESS_POS_SIG false
#define PROCESS_CEBR_ALL false

// Calibration File Name Templates
#define RUN_FILE_NAME_TEMPLATE "root_data_70Ge_run%03d.mvmelst.bin_tree.root"

#endif // CACONFIGURATION_HPP