

#include <string>
#include <map>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include <iostream>

#include "pars.h"


std::map<std::string, std::string> pars;


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

        std::getline(ifs, key, '=');
        std::getline(ifs, val, '\n');

        pars[key] = val;

        // std::cout << "   " << key << " = " << val << "\n";
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

template<typename T>
bool tfetch(std::string s, T& ret) {
    std::string val;
    if (!sfetch(s, val))
        return false;

    std::istringstream in(val);
    in >> ret;
    return !in.fail(); // Todo: Check conversion.
}


bool dfetch(std::string s, double& ret) {
    return tfetch(s, ret);
}
bool ifetch(std::string s, int& ret) {
    return tfetch(s, ret);
}
