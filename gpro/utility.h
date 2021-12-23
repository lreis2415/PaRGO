/*!
 * \file utility.h
 * \brief Utility functions derived from the Common Cross-platform Geographic Library (CCGL)
 *        https://github.com/crazyzlj/CCGL
 *
 * Changelog:
 *   - 1. 2019-07-18  Transport from CCGL (master branch of commit '72460fc573e526ef3edbae40d8e85a32a8aefd6c')
 *
 * \author Liangjun Zhu (crazyzlj)
 */

#ifndef UTILITY_H
#define UTILITY_H

/*! `NDEBUG` or `_DEBUG` mean not build on `DEBUG` mode. */
#ifndef NDEBUG
#ifndef _DEBUG
#define _DEBUG
#endif /* _DEBUG */
#endif /* NDEBUG */

/*! A reference to x64 architecture */
#if defined(_WIN64) || defined(__x86_64) || defined(__LP64__)
#define CPP_64
#endif

/*! A reference to MSVC environment */
#if defined _MSC_VER
#define CPP_MSVC
#endif /* _MSC_VER */

/*! A reference to Intel C++ compiler */
#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC)
#define CPP_ICC
#endif /* __INTEL_COMPILER */

/*! A reference to GCC compiler */
#if defined(__GNUC__)
#define CPP_GCC
/*! A reference to GCC compiler on macOS */
#if defined(__APPLE__)
#define CPP_APPLE
#endif /* __APPLE__ */
#endif /* __GNUC__ */

#include <memory>
#include <stdexcept>
#include <cfloat>
#include <map>
#include <string>
#include <vector>
#include <cstring> // strcasecmp in GCC

/// platform
#if defined WINDOWS
// For MSVC and MINGW64 in Windows OS
// #define _WINSOCKAPI_    // In order to stop windows.h including winsock.h
// _WINSOCKAPI_ is defined by <winsock2.h>
#include <winsock2.h>
#include <windows.h>
#endif /* WINDOWS */

#if defined CPP_GCC
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <cerrno>
#endif /* CPP_GCC */

// define some macro for string related built-in functions
#ifdef CPP_MSVC
#define stringcat strcat_s
#define stringcpy strcpy_s
#define strprintf sprintf_s
#define strtok strtok_s
#define stringscanf sscanf_s
#else
#define stringcat strcat
#define stringcpy strcpy
#define strprintf snprintf
#define strtok strtok_r
#define stringscanf sscanf
#endif /* CPP_MSVC */

#if defined(__MINGW32_MAJOR_VERSION) || defined(__MINGW64_VERSION_MAJOR)
#define MINGW
#endif

#if defined(MINGW) || defined(_MSC_VER)
#define strcasecmp _stricmp
#endif /* MINGW or MSVC */

#if defined(__clang__)
// Clang
#if ((__clang_major__ * 100) + __clang_minor__) >= 400
#if __has_feature(cxx_noexcept)
#define HAS_NOEXCEPT
#endif /* NOEXCEPT */
#if __has_feature(cxx_override_control)
#define HAS_OVERRIDE
#endif /* OVERRIDE */
#if __has_feature(cxx_variadic_templates)
#define HAS_VARIADIC_TEMPLATES
#endif /* VARIADIC_TEMPLATES */
#if __has_feature(cxx_defaulted_functions)
#define HAS_DEFAULT_FUNC
#endif /* #define HAS_DEFAULT_FUNC */
#endif /* Clang */
#elif defined(CPP_ICC)
// Intel C++
#if ((__INTEL_COMPILER >= 1400) && (__INTEL_COMPILER != 9999)) || (__ICL >= 1400)
#define HAS_NOEXCEPT
#define HAS_OVERRIDE
#define HAS_VARIADIC_TEMPLATES
#define HAS_DEFAULT_FUNC
#endif /* Intel C++ */
#elif defined(CPP_GCC)
// GNU GCC
#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 406 && (__cplusplus >= 201103L || (defined(__GXX_EXPERIMENTAL_CXX0X__) && __GXX_EXPERIMENTAL_CXX0X__))
#define HAS_NOEXCEPT
#define HAS_OVERRIDE
#define HAS_VARIADIC_TEMPLATES
#define HAS_DEFAULT_FUNC
#endif /* GCC */
#elif defined(_MSC_VER)
// MS Visual C++
#if _MSC_VER >= 1900
#define HAS_NOEXCEPT
#endif /* Visual Studio 2015 or later */
#if _MSC_VER >= 1800
#define HAS_VARIADIC_TEMPLATES
#define HAS_DEFAULT_FUNC
#endif /* Visual Studio 2013 or later */
#if _MSC_VER>= 1600
#define HAS_OVERRIDE
#endif /* Visual Studio 2010 or later */
#endif /* Figure out HAS_NOEXCEPT, HAS_VARIADIC_TEMPLATES, and HAS_OVERRIDE or not */

/*! A compatible reference to `noexcept` or `throw()` if not supported by the compiler. */
#ifdef HAS_NOEXCEPT
#define NOEXCEPT noexcept
#else
#define NOEXCEPT throw()
#endif /* HAS_NOEXCEPT */

/*! A compatible reference to `override` or blank if not supported by the compiler. */
#ifdef HAS_OVERRIDE
#define OVERRIDE override
#else
#define OVERRIDE
#endif /* HAS_OVERRIDE */

/*! A compatible reference to `override` or blank if not supported by the compiler. */
#ifdef HAS_DEFAULT_FUNC
#define DEFAULT =default
#else
#define DEFAULT {}
#endif /* HAS_OVERRIDE */

/*
* Avoid the compile error on MSVC like this:
*   warning C4251: 'CLASS_TEST::m_structs':
*           class 'std::vector<_Ty>' needs to have dll-interface to be used by clients of class
* refers to http://www.cnblogs.com/duboway/p/3332057.html
*/
#ifdef MSVC
#define DLL_STL_LIST(STL_API, STL_TYPE) \
    template class STL_API std::allocator< STL_TYPE >; \
    template class STL_API std::vector<STL_TYPE, std::allocator< STL_TYPE > >;
#endif /* MSVC */

#ifdef USE_GDAL
/* Ignore warning on Windows MSVC compiler caused by GDAL.
* refers to http://blog.csdn.net/liminlu0314/article/details/8227518
*/
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#pragma warning(disable: 4100 4190 4251 4275 4305 4309 4819 4996)
#endif /* Ignore warnings of GDAL */
#endif /* USE_GDAL */

#ifdef WINDOWS
#define SEP             "\\"
#else
#define SEP             "/"
#endif /* WINDOWS */

/*! Default NoData value for raster data etc. */
#ifndef NODATA_VALUE
#define NODATA_VALUE    (-9999.0f)
#endif /* NODATA_VALUE */

/*! Missing float value */
#ifndef MISSINGFLOAT
#define MISSINGFLOAT    (-1 * FLT_MAX)
#endif /* MISSINGFLOAT */

/*! Maximum float value */
#ifndef MAXIMUMFLOAT
#define MAXIMUMFLOAT    FLT_MAX
#endif /* MAXIMUMFLOAT */

/*! Maximum length of full file path */
#ifndef PATH_MAX
#define PATH_MAX        1024
#endif /* PATH_MAX */

/*! A approximation of Zero */
#ifndef UTIL_ZERO
#define UTIL_ZERO       1.0e-6f
#endif /* UTIL_ZERO */

/*! A approximation of PI */
#ifndef PI
#define PI              3.14159265358979323846f
#endif /* PI */

using std::vector;
using std::string;

/*!
 * \brief return the rank of this process
 * \return rank
 */
int GetRank();

/*!
 * \brief Trim given string's heading and tailing by "<space>,\n,\t,\r"
 * \sa TrimSpaces
 * \param[in] s \a string information
 * \return Trimmed string
 */
string& Trim(string& s);

/*!
 * \brief Get Uppercase of given string
 * \param[in] str
 * \return Uppercase string
 */
string GetUpper(const string& str);

/*!
 * \brief Match \a char ignore cases
 * \param[in] a, b \a char*
 * \return true or false
 * \sa StringMatch()
 */
bool StringMatch(const char* a, const char* b);

/*!
 * \brief Match Strings in UPPERCASE manner
 * \param[in] text1, text2
 * \return true or false
 */
bool StringMatch(const string& text1, const string& text2);

/*!
 * \brief Splits the given string by spaces
 * \param[in] item \a string information
 * \return The split strings vector
 */
vector<string> SplitString(const string& item);

/*!
 * \brief Splits the given string based on the given delimiter
 * \param[in] item \a string information
 * \param[in] delimiter \a char
 * \return The split strings vector
 */
vector<string> SplitString(const string& item, char delimiter);

/*!
 * \brief Check the given directory path is exists or not.
 */
bool DirectoryExists(const string& dirpath);

/*!
 * \brief Clean a directory if exists, otherwise create it.
 */
bool CleanDirectory(const string& dirpath);

/*!
 * \brief Delete a directory if exists.
 *
 * Reference:
 *   - 1. Windows: https://stackoverflow.com/questions/734717/how-to-delete-a-folder-in-c
 *   - 2. Linux: https://www.linuxquestions.org/questions/programming-9/deleting-a-directory-using-c-in-linux-248696/
 */
bool DeleteDirectory(const string& dirpath, bool del_subdirs = true);

/*!
 * \brief Get the root path of the current executable file
 * \return \a string root path
 */
string GetAppPath();

/*!
 * \brief Return the absolute file path from a given file path
 * \param[in] full_filename
 * \return absolutePath
 * \sa GetPathFromFullName
 */
string GetAbsolutePath(string const& full_filename);

/*!
 * \brief Return the file name from a given file's path
 * \param[in] full_filename
 * \return CoreFileName
 * \sa GetPathFromFullName
 */
string GetCoreFileName(string const& full_filename);

/*!
 * \brief Return the suffix of a given file's path
 * \param[in] full_filename
 * \return Suffix
 * \sa GetPathFromFullName
 */
string GetSuffix(string const& full_filename);

/*!
 * \brief Replace the suffix by a given suffix
 * \param[in] full_filename
 * \param[in] new_suffix
 * \return new full_filename
 */
string ReplaceSuffix(string const& full_filename, string const& new_suffix);

/*!
 * \brief Get Path From full file path string
 * \param[in] full_filename \a string
 * \return filePath string
 * \sa GetCoreFileName
 */
string GetPathFromFullName(string const& full_filename);

/*!
 * \brief Return a flag indicating if the given file exists
 * \param[in] filename String path of file
 * \return True if Exists, and false if not.
 */
bool FileExists(string const& filename);

/*!
 * \brief Return a flag indicating if the given path exists
 * \param[in] path String path
 * \return True if Exists, and false if not.
 */
bool PathExists(string const& path);

/*!
 * \brief Delete the given file if existed.
 * \param[in] filepath \a string File path, full path or relative path
 * \return 0 if deleted successful, else return nonzero value, e.g. -1.
 */
int DeleteExistedFile(const string& filepath);

/*!
 * \brief Find files in given paths
 * \param[in] lp_path, expression
 * \param[out] vec_files
 * \return 0 means success
 */
int FindFiles(const char* lp_path, const char* expression, vector<string>& vec_files);

/*!
 * \brief Load short plain text file as string vector, ignore comments begin with '#' and empty lines
 * \param[in] filepath Plain text file path
 * \param[out] content_strs Each line without CRLF or LF stored in vector
 * \return True when read successfully, and false with empty content_strs when failed
 */
bool LoadPlainTextFile(const string& filepath, vector<string>& content_strs);

#endif
