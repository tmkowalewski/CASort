#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

// Number of hardware threads to use in processing
#define N_THREADS 1 // std::thread::hardware_concurrency() // Number of threads to use for processing, defaults to system max

// Which modules to process
#define PROCESS_CLOVER_CROSS true
#define PROCESS_CLOVER_BACK false
#define PROCESS_POS_SIG false
#define PROCESS_CEBR_ALL false

// Calibration File Directories
#define ENERGY_CAL_DIR "/home/tylermk/TUNL/Data/NRF/70Ge/energy_calibration/calibrations"
#define GAIN_MATCH_DIR "energy_calibration/gain_matching/"

// Config Constants
#define CLOVER_ENERGY_THRESHOLD 150  // Energy threshold (keV) for clover detectors
#define CLOVER_ADDBACK_THRESHOLD 150 // Energy threshold (keV) for add-back
#define CLOVER_ADDBACK_WINDOW 150    // Time window (ns) around primary hit for add-back
#define ENERGY_PER_BIN 0.25          // Energy per bin (keV) for clover detectors for 16-bit digitizer
#define MAX_ENERGY 10000             // Maximum energy (keV) for clover detectors
#define CLOVER_COIN_WINDOW 100       // Time window (ns) for clover coincidence

#define POS_SIG_ENERGY_THRESHOLD 150 // Energy threshold (keV) for positive signal detectors

#define CEBR_ALL_ENERGY_THRESHOLD 150 // Energy threshold (keV) for CeBr (also LaBr and MPAD) detectors

#endif // CONFIGURATION_HPP