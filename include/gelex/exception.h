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

// ----------------------------------------------------------------
// ----------------- File-related exceptions ----------------------
// ----------------------------------------------------------------
class FileNotFoundException : public GelexException
{
   public:
    explicit FileNotFoundException(const std::filesystem::path& path)
        : GelexException(std::format("File not found: '{}'", path.string()))
    {
    }
};

class FileOpenException : public GelexException
{
   public:
    explicit FileOpenException(std::string_view message)
        : GelexException(std::string(message))
    {
    }
};

class FileWriteException : public GelexException
{
   public:
    explicit FileWriteException(std::string_view message)
        : GelexException(std::string(message))
    {
    }
};

class FileExistsException : public GelexException
{
   public:
    explicit FileExistsException(const std::filesystem::path& path)
        : GelexException(
              std::format("Output file already exists: '{}'", path.string()))
    {
    }
};

class FileFormatException : public GelexException
{
   public:
    using GelexException::GelexException;
};

// ----------------------------------------------------------------
// ------------------ Parse-related exceptions ---------------------
// ----------------------------------------------------------------

class DataParseException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class NumberParseException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class HeaderFormatException : public GelexException
{
   public:
    using GelexException::GelexException;
};

// ---------------------------------------------------------------
// -------------- Argument-related exceptions --------------------
// ---------------------------------------------------------------

class ArgumentValidationException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class ColumnRangeException : public GelexException
{
   public:
    using GelexException::GelexException;
};

class InvalidOperationException : public GelexException
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
