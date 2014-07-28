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

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[ ";
    for( auto e : v ) std::cout << e << " ";
    std::cout << "]";
    return os;
}

}

#endif // HELPER_H
