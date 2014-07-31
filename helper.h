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
#ifndef HELPER_H
#define HELPER_H

#include "configure.h"

#include "point.h"
#include <iostream>

namespace NAMESPACE {

template<typename Float>
Float sqr( Float a ) {
    return a*a;
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
int pos_neg(T val) {
    return (T(0) <= val) ? 1 : -1;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[ ";
    for( auto e : v ) std::cout << e << " ";
    std::cout << "]";
    return os;
}

}

#endif // HELPER_H
