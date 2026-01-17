#ifndef GELEX_EXCEPTION_H_
#define GELEX_EXCEPTION_H_

#include <Eigen/Core>
#include <stdexcept>

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
    using GelexException::GelexException;
};

class FileOpenException : public GelexException
{
    using GelexException::GelexException;
};

class FileWriteException : public GelexException
{
    using GelexException::GelexException;
};

class FileExistsException : public GelexException
{
    using GelexException::GelexException;
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

class InvalidInputException : public GelexException
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
}  // namespace gelex

#endif  // GELEX_EXCEPTION_H_
