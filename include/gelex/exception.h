/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
