#ifndef TCAEXPERIMENT_HPP
#define TCAEXPERIMENT_HPP

// Standard C++ includes
#include <string> // std::string
#include <vector> // std::vector

// ROOT includes
#include <TNamed.h>

// Forward declarations to avoid circular dependencies
class TCADAQModule; // Represents a data acquisition module
class TCARun;       // Represents an experimental run

// Experiment class
// Encapsulates the name and configuration of an experiment, along with lists
// of DAQ modules and runs. Provides methods to add, remove, and query these.
class TCAExperiment : public TNamed
{
public:
    // Constructors

    TCAExperiment();


    TCAExperiment(const std::string& file_name);


    ~TCAExperiment();

    // Accessors (Getters)

    // Returns the name (path) of the configuration file.
    const std::string& GetFileName() const { return file_name_; }

    // Finds and returns a pointer to the DAQModule with the given name.
    // @param module_name: Name of the desired DAQ module.
    // @return Pointer to the module if found, nullptr otherwise.
    const TCADAQModule* GetDAQModule(const std::string& module_name) const;

    // Returns a const reference to the internal list of DAQModule pointers.
    const std::vector<TCADAQModule*>& GetDAQModules() const { return daq_modules_; }




    void AddDAQModule(TCADAQModule* module);


    void RemoveDAQModule(const std::string& module_name);


    void AddRun(TCARun* run);


    void RemoveRun(unsigned int run_number);


    void PrintInfo() const;

private:
    std::string file_name_;                // Configuration file path
    std::vector<TCADAQModule*> daq_modules_; // Owned DAQ modules
    std::vector<TCARun*> runs_;              // Owned runs
};

#endif // TCAEXPERIMENT_HPP