/*
 * (c) 30.07.2014 Martin Hünniger
 *
 *   This file is part of DelDec.
 *
 *   DelDec is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   DelDec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
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
