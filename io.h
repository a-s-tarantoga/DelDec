/*
 * (c) 30.07.2014 Martin Hünniger
 */
#ifndef IO_H
#define IO_H

#include "configure.h"
#include "point.h"

#include <vector>
#include <istream>
#include <fstream>
#include <strstream>

template<int dim, typename Float>
std::ifstream& operator>>(std::ifstream& ifs, std::vector<DD::Point<dim,double>>) {

    return ifs;
}

#endif // IO_H
