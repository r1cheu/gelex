#pragma once

#include <filesystem>
#include <format>
#include <string>

namespace gelex
{

enum class ErrorCode : uint8_t
{
    FileNotFound,
    FileIOError,
    InvalidData,
    OutputFileExists,

    NotNumber,
    InvalidFile,
    InvalidRange,
    InvalidArgument,

    WrongHeader,
    InconsistColumnCount,

    Unknown,
};

struct Error
{
    ErrorCode code;
    std::string message;
};

constexpr std::string_view to_string_view(ErrorCode code) noexcept
{
    using namespace std::string_view_literals;

    switch (code)
    {
        case ErrorCode::FileNotFound:
            return "File not found"sv;
        case ErrorCode::FileIOError:
            return "File io error"sv;
        case ErrorCode::NotNumber:
            return "Not a number"sv;
        case ErrorCode::OutputFileExists:
            return "OutputFileExists"sv;
        case ErrorCode::InvalidFile:
            return "Invalid file format"sv;
        case ErrorCode::InvalidRange:
            return "Invalid Range specified"sv;
        case ErrorCode::InvalidArgument:
            return "Invalid argument"sv;
        case ErrorCode::WrongHeader:
            return "Incorrect header"sv;
        case ErrorCode::InconsistColumnCount:
            return "Inconsistent column count"sv;
        case ErrorCode::InvalidData:
            return "Invalid data input"sv;
        case ErrorCode::Unknown:
            return "An unknown error occurred"sv;
        default:
            break;
    }
    return "An unknown error occurred"sv;
}

template <typename E>
E enrich_with_line_info(E&& error, int line_number)
{
    std::string original_message = std::move(error.message);
    error.message = std::format("{} (line {})", original_message, line_number);
    return std::forward<E>(error);
}

template <typename E>
E enrich_with_file_info(E&& error, std::string_view path)
{
    std::string original_message = std::move(error.message);
    error.message = std::format("{} (file [{}])", original_message, path);
    return std::forward<E>(error);
}

template <typename E>
E enrich_with_file_info(E&& error, const std::filesystem::path& path)
{
    std::string original_message = std::move(error.message);
    error.message
        = std::format("{} (file [{}])", original_message, path.string());
    return std::forward<E>(error);
}

}  // namespace gelex
