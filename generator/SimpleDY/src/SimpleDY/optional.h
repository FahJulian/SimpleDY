#pragma once

#include "base.h"

namespace SimpleDY
{
    struct _Nullopt
    {
    };

    template<typename T>
    class Optional
    {
    public:
        constexpr Optional()
            : m_hasValue(false), _()
        {
        }

        constexpr Optional(_Nullopt _)
            : m_hasValue(false), _()
        {
        }

        Optional(const T& value)
            : m_hasValue(true), m_value(value)
        {
        }

        ~Optional()
        {
            if (m_hasValue)
                m_value.~T();
        }

        T* operator->()
        {
            return &m_value;
        }

        const T* operator->() const
        {
            return &m_value;
        }

        bool hasValue() const
        {
            return m_hasValue;
        }

        T& getValue()
        {
            return m_value;
        }

        const T& getValue() const
        {
            return m_value;
        }

    private:
        bool m_hasValue;

        union
        {
            T m_value;
            char _;
        };
    };

    constexpr _Nullopt NULLOPT = { };

} // namespace SimpleDY
