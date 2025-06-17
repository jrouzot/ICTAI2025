//
// Created by Julien Rouzot on 20/05/24.
//

#pragma once

#include <iostream>
#include <sstream>
#include <string>

namespace logger {

    enum class Verbosity {
        Error,
        Warn,
        Debug
    };

    class Logger {
    public:
        explicit Logger(Verbosity verbosity) : mVerbosity(verbosity) {}

        template<typename... Args>
        void Debug(Args... args) {
            if (mVerbosity >= Verbosity::Debug) {
                Log("DEBUG", args...);
            }
        }

        template<typename... Args>
        void Warn(Args... args) {
            if (mVerbosity >= Verbosity::Warn) {
                Log("WARN", args...);
            }
        }

        template<typename... Args>
        void Error(Args... args) {
            if (mVerbosity >= Verbosity::Error) {
                Log("ERROR", args...);
            }
        }

        template<typename T>
        void DebugVector(const std::vector<T>& vec) {
            if (mVerbosity >= Verbosity::Debug) {
                std::ostringstream oss;
                oss << "[";
                for (int i(0); i < vec.size(); ++i) {
                    oss << vec[i];
                    if (i < vec.size() - 1) {
                        oss << ", ";
                    }
                }
                oss << "]";
                std::cout << "[DEBUG]: " << oss.str() << std::endl;
            }
        }

    private:
        Verbosity mVerbosity;

        template<typename... Args>
        void Log(const char* level, Args... args) {
            std::ostringstream oss;
            (oss << ... << args);  // Fold expression to unpack and append all arguments
            std::cout << "[" << level << "]: " << oss.str() << std::endl;
        }
    };

} // namespace logger


