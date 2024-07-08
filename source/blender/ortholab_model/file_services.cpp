#include "file_services.h"
#include <spdlog/spdlog.h>
#include <filesystem>


namespace fs = std::filesystem;

namespace OrthoLab
{
    std::string ensure_directory_exists(const std::string &relateive_dir)
    {
        fs::path current = fs::current_path();
        fs::path new_path = current / relateive_dir;
        //std::cout << "Check for Directory: " << new_path << std::endl;
        if (!fs::exists(new_path))
        { // Check if the directory exists
            spdlog::trace("Directory not exist, creating: {}" , relateive_dir);
            if (fs::create_directory(new_path))
            { // Try to create the directory 
                spdlog::trace("Directory created {}", relateive_dir);
            }
            else
            {
                spdlog::error("Failed to create directory {}", relateive_dir);
            }
        }
        else
        {
            spdlog::trace("Directory already exists: {}", relateive_dir); 
        }
        return new_path.string();
    }
    bool file_exists(const std::string &relateive_dir)
    {
        fs::path current = fs::current_path();
        fs::path new_path = current / relateive_dir;
        return fs::exists(new_path);
    }
}