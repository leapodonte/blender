#include "base_helpers.h"

#ifdef _WIN32
#include <windows.h>
BOOL setEnv(const std::string &name, const std::string &value) {
    return SetEnvironmentVariable(name.c_str(), value.c_str());
}
#else
#include <cstdlib>
int setEnv(const std::string &name, const std::string &value, int overwrite = 1) {
    return setenv(name.c_str(), value.c_str(), overwrite);
}
#endif


namespace OrthoLab {

    float get_variable_float(const char* v_name, float default_value)
    {
        float value = 0;
        const char *env_var = getenv(v_name);

        if (env_var != nullptr)
        {
            try
            {
                value = std::stof(env_var);
                std::cout << v_name<< " the integer value is: " << value << std::endl;
            }
            catch (const std::invalid_argument &e)
            {
                std::cerr << v_name << " invalid number: " << env_var << std::endl;
            }
            catch (const std::out_of_range &e)
            {
                std::cerr << v_name << " number out of range: " << env_var << std::endl;
            }
        }
        else
        {
            std::cout << "No value defined for: " << v_name << std::endl;
            return default_value;
        }
        return value;
    }

    bool set_env_variable(const std::string &name, const std::string &v)
    {
        // Set the environment variable
        std::cout << "set_env_variable: " << name <<" = " << v << std::endl;
        bool success = _putenv((name+"="+v).c_str());
        return success;
    }

} // namespace OrthoLab
