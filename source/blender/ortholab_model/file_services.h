#ifndef FILE_SERVICES_H_
#define FILE_SERVICES_H_

#include <iostream>

namespace OrthoLab
{
    std::string ensure_directory_exists(const std::string &relateive_dir);

    bool file_exists(const std::string &relateive_dir);    
}
#endif // FILE_SERVICES_H_

