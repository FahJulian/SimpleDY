#pragma once

#include "SimpleDY/base.h"

#include <string>
#include <sstream>

namespace SimpleDY
{
    class Process;

    class LesHouchesSerializer
    {
    public:
        LesHouchesSerializer(const Process& process)
            : m_process(process)
        {
        }

        void serialize(const std::string& filePath);

    private:
        const Process& m_process;
    };

} // namespace SimpleDY
