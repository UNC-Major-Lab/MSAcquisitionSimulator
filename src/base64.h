//
// Created by Dennis Goldfarb on 9/4/15.
//

#ifndef MSACQUISITIONSIMULATOR_BASE64_H
#define MSACQUISITIONSIMULATOR_BASE64_H


#include <string>

std::string base64_encode(unsigned char const* , unsigned int len);
std::string base64_decode(std::string const& s);

#endif //MSACQUISITIONSIMULATOR_BASE64_H
