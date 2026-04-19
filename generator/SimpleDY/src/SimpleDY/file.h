#pragma once

#include "base.h"

#include <string>

namespace SimpleDY
{
    class File
    {
    public:
        File(const std::string& filePath);

        std::string read() const;

        // std::vector<std::string> readLines() const;

        void write(const std::string& text);

    private:
        std::string mFilePath;
    };
}
