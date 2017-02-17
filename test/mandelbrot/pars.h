
#ifndef DELTAWING_PARS_H
#define DELTAWING_PARS_H


#include <string>
#include <map>


extern std::map<std::string, std::string> pars;


void getpars(std::string parfile);
bool sfetch(std::string s, std::string& ret);

template<typename T>
inline bool tfetch(std::string s, T& ret) {
    std::string val;
    if (!sfetch(s, val))
        return false;

    std::istringstream in(val);
    in >> ret;
    return !in.fail(); // Todo: Check conversion?
}

#endif
