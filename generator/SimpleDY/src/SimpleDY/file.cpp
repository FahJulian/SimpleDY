#include "file.h"

#include <sstream>
#include <fstream>

namespace SimpleDY
{
    File::File(const std::string& filePath)
        : mFilePath(filePath)
    {
    }

    std::string File::read() const
    {
        std::ifstream stream = std::ifstream(mFilePath.c_str());

        std::stringstream buffer;
        buffer << stream.rdbuf();

        buffer.seekg(0, buffer.end);
        int size = buffer.tellg();

        char* content = new char[size];

        buffer.seekg(std::ifstream::beg);
        buffer.read(content, size);

        stream.close();
        return content;
    }

    void File::write(const std::string& text)
    {
        std::ofstream stream = std::ofstream(mFilePath.c_str());

        stream << text;

        stream.close();
    }
}
