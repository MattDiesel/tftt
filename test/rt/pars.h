
#ifndef DELTAWING_PARS_H
#define DELTAWING_PARS_H


#include <string>
#include <map>


extern std::map<std::string, std::string> pars;


void getpars(std::string parfile);

bool dfetch(std::string s, double& ret);
bool ifetch(std::string s, int& ret);
bool sfetch(std::string s, int& ret);


#endif
