
#define PARS_DEBUG


#include <string>
#include <map>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <algorithm>

#ifdef PARS_DEBUG
#include <iostream>
#endif

#include "pars.h"


namespace util {


std::map<std::string, std::string> pars;


std::string& rtrim(std::string& str)
{
    auto it1 = std::find_if(str.rbegin(), str.rend(), [](char ch){
        return !std::isspace<char>(ch , std::locale::classic());
    });
    str.erase( it1.base() , str.end() );
    return str;
}


void getpars(std::string parfile) {
    std::ifstream ifs(parfile);

    ifs.setf(std::ios_base::skipws);

    int i = 0;
    std::string key;
    std::string val;
    for (;ifs;) {
        // Ignore WS
        ifs >> std::ws;

        if (ifs.peek() == '#') {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        ifs >> std::ws; std::getline(ifs, key, '=');
        ifs >> std::ws; std::getline(ifs, val, '\n');

        pars[rtrim(key)] = val;

#ifdef PARS_DEBUG
        std::cout << "   '" << key << "' = '" << val << "'\n";
#endif
    }
}


bool sfetch(std::string s, std::string& ret) {
    auto search = pars.find(s);
    if(search != pars.end()) {
        ret = search->second;
        return true;
    }
    return false;
}


} // namespace util
