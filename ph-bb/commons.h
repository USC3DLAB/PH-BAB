//
//  commons.h
//  PH-BB
//
//  Created by Semih Atakan on 7/5/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef phbab_commons_h
#define phbab_commons_h

#include <iostream>
#include <string>
#include <istream>

//#include <string>


const std::string currentDateTime();

double get_cpu_time();
double get_wall_time();

std::istream& safeGetline(std::istream& is, std::string& t);

#endif
