#ifndef CACONFIGURATION_HPP
#define CACONFIGURATION_HPP

// C++ Includes
#include <thread>
#include <algorithm>

// Number of hardware threads to use in processing
extern const int kMaxThreads = std::thread::hardware_concurrency(); // Number of threads to use for processing, defaults to system max

// Which modules to process
#define PROCESS_CLOVER_CROSS true
#define PROCESS_CLOVER_BACK false
#define PROCESS_POS_SIG false
#define PROCESS_CEBR_ALL false

// Calibration File Name Templates
#define RUN_FILE_NAME_TEMPLATE "root_data_70Ge_run%s.mvmelst.bin_tree.root"

// Config Constants
#define ENERGY_PER_BIN 0.25    // Energy per bin (keV) for clover detectors for 16-bit digitizer
#define MAX_ENERGY 10000       // Maximum energy (keV) for clover detectors
#define CLOVER_COIN_WINDOW 100 // Time window (ns) for clover coincidence

#define POS_SIG_ENERGY_THRESHOLD 150 // Energy threshold (keV) for positive signal detectors

#define CEBR_ALL_ENERGY_THRESHOLD 150 // Energy threshold (keV) for CeBr (also LaBr and MPAD) detectors

#endif // CACONFIGURATION_HPP