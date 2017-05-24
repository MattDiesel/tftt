
#ifndef FORMATSTRING_H
#define FORMATSTRING_H

#include <string>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <iomanip>


namespace util {

namespace fmt {


/** Prints the Nth argument of a packed parameter list (specialization for length = 1)
 *
 * @author Diesel
 */
template<typename T>
inline bool print(std::ostream& os, int N, T v)
{
    if (!N) {
        os << v;
        return true;
    }
    return false;
}

/** Prints the Nth argument of a packed parameter list (specialization for length = 1)
 *
 * @author Diesel
 */
template<>
inline bool print(std::ostream& os, int N, int v)
{
    if (!N) {
        os << v;
        return true;
    }
    return false;
}

/** Prints the Nth argument of a packed parameter list.
 *
 * @author Diesel
 */
template<typename T, typename... Args>
inline bool print(std::ostream& os, int N, T first, Args... args)
{
    if (!N) {
        os << first;
        return true;
    }
    return print(os, N-1, args...);
}


} // namespace fmt


/** Formats a string according to a format pattern.
 *
 * @author Diesel
 */
template<class ...TA>
inline std::string formatString(std::string f, TA... pargs)
{
    std::ostringstream oss;

    std::string::size_type c = 0, p, n;
    int w;
    for (; c < f.length();) {
        p = f.find('{', c);

        if (p == std::string::npos) {
            oss << f.substr(c);
            break;
        }

        oss << f.substr(c, p-c);

        n = 0;
        while (std::isdigit(f[++p])) {
            n = n*10 + (f[p]-'0');
        }

        if (f[p] == ',') {
            w = 0;
            while (std::isdigit(f[++p]))
                w = w*10 + (f[p]-'0');

            oss << std::setw(w);

            if (f[p] == '<') { // Left align
                oss << std::left;
                p++;
            }
            else if (f[p] == '>') {
                oss << std::right;
                p++;
            }
            else if (f[p] == '|') {
                oss << std::internal;
                p++;
            }

            if (f[p] != '}') {
                oss << std::setfill(f[p++]);
            }
        }

        if (f[p] != '}')
            throw std::invalid_argument("Format string bad.");
        p++;

        if (!fmt::print(oss, n, pargs...)) {
            throw std::out_of_range("Index was out of range of arguments");
        }

        c = p;
    }

    return oss.str();
}


} // namespace util


#endif
