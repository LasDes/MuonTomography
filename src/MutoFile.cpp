#include "MutoFile.h"

// get current working directory 
#ifdef _WIN32
    #include <windows.h>
    #include <direct.h>
    char winpath[ MAX_PATH ]; 
    std::string MutoFile::fCwd = std::string( _getcwd( winpath, MAX_PATH ) );
#endif

#ifdef __unix__
    #include <unistd.h>
    char unixpath[ PATH_MAX ]; 
    std::string MutoFile::fCwd = std::string( getcwd( unixpath, PATH_MAX ) );
#endif



std::string MutoFile::getCurrentPath() {
    return fCwd;
}




