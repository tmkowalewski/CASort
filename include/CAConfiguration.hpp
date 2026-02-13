#ifndef CACONFIGURATION_HPP
#define CACONFIGURATION_HPP

// C++ Includes
#include <algorithm>
#include <thread>

// Which modules to process
#define PROCESS_POS_SIG false
#define PROCESS_CEBR_ALL false

// Number of hardware threads to use in processing
const unsigned int kMaxThreads = std::min(20U, std::thread::hardware_concurrency()); // Number of threads to use for processing, defaults to system max

// Debug Mode
#define DEBUG 1

// Calibration File Name Templates
#define RUN_FILE_NAME_TEMPLATE "root_data_70Ge_run%03d.mvmelst.bin_tree.root"

#endif // CACONFIGURATION_HPP