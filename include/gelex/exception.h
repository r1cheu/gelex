#pragma once

#include <filesystem>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>

namespace gelex
{

class GelexException : public std::runtime_error
{
   public:
    using std::runtime_error::runtime_error;
};

class FileNotFoundException : public GelexException
{
   public:
    explicit FileNotFoundException(const std::filesystem::path& path)
        : GelexException(std::format("File not found: '{}'", path.string()))
    {
    }
};

class FileIOException : public GelexException
{
   public:
    explicit FileIOException(std::string_view message)
        : GelexException(std::string(message))
    {
    }
};

class InvalidDataException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class OutputFileExistsException : public GelexException
{
   public:
    explicit OutputFileExistsException(const std::filesystem::path& path)
        : GelexException(
              std::format("Output file already exists: '{}'", path.string()))
    {
    }
};

class InvalidOperationException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class NotNumberException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class InvalidFileException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class InvalidRangeException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class InvalidArgumentException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class WrongHeaderException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class InconsistentColumnCountException : public GelexException
{
   public:
    using GelexException::GelexException;
};

// Helper functions for exception message formatting

inline std::string enrich_with_line_info(
    std::string_view message,
    int line_number)
{
    return std::format("{} (line {})", message, line_number);
}

inline std::string enrich_with_file_info(
    std::string_view message,
    std::string_view path)
{
    return std::format("{} (file [{}])", message, path);
}

inline std::string enrich_with_file_info(
    std::string_view message,
    const std::filesystem::path& path)
{
    return std::format("{} (file [{}])", message, path.string());
}

}  // namespace gelex
