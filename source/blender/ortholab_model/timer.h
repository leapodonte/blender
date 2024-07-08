#ifndef ORTHO_TIMER_H
#define ORTHO_TIMER_H

#include <chrono>
namespace OrthoLab
{
    class Timer
    {
    public:
        Timer() : start_time_(std::chrono::steady_clock::now()) {}

        void reset()
        {
            start_time_ = std::chrono::steady_clock::now();
        }

        double elapsed() const
        {
            return std::chrono::duration_cast<std::chrono::duration<double>>(
                       std::chrono::steady_clock::now() - start_time_)
                .count();
        }

    private:
        std::chrono::time_point<std::chrono::steady_clock> start_time_;
    };
} // namespace OrthoLab
#endif
