/* Configuration class for the program using JSON */

#pragma once 

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class MutoConfig {

public:
    MutoConfig() {}
    MutoConfig(std::string filename) { 
        if (!loadFile(filename)) {
            std::cerr << "Configuration file load error for: " << filename << std::endl;
        } 
    }
    ~MutoConfig() {}

    bool loadFile (std::string filename) {
        std::ifstream infile(filename.c_str());
        if (!infile.fail()) {
            infile >> fConfig;
            infile.close();
            return true;
        } else {
            return false;
        }
    }

    std::string getConfigString() {
        return fConfig.dump();
    }

    json::value_type get(std::string attr) {
        return fConfig[attr];
    }

    bool exist(std::string key) {
        return fConfig.find(key) != fConfig.end();
    }

private:
    json fConfig;
};