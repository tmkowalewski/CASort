#ifndef CAGAINDRIFT_HPP
#define CAGAINDRIFT_HPP

// C++ Includes
#include <string>

// ROOT Includes

// Project Includes

namespace CAGainDrift
{

    void ApplyGainDriftCorrections(const std::string& gain_shift_dir,
        const std::string& run_file_dir);

} // namespace CAGainDrift

#endif // CAGAINDRIFT_HPP