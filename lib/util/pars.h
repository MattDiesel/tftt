
#ifndef UTIL_PARS_H
#define UTIL_PARS_H


#include <string>
#include <map>
#include <sstream>


namespace util {


extern std::map<std::string, std::string> pars;

void getpars(std::string parfile);
void getpars(int argc, char* argv[]);

bool sfetch(std::string s, std::string& ret);


template<typename T>
inline bool tfetch(std::string s, T& ret) {
    std::string val;
    if (!sfetch(s, val))
        return false;

    std::istringstream in(val);
    in >> ret;
    return !in.fail(); 
}


} // namespace util


#endif
