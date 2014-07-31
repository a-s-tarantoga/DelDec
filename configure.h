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

#ifndef FC_CONFIGURE_H
#define FC_CONFIGURE_H

#include <iostream>

#ifdef NDEBUG
  #undef NDEBUG
#endif

#include <assert.h>
#include <iostream>
#include <iomanip>

#define NAMESPACE DD         // Namespace of the library
//#define ASSERTION_MODE       // enables (cheap) assertions
//#define DEBUG_MODE           // enables debugging checks(cheap and
                                // expensive ones)
#define SAVE_MEMORY          // enables compression of the subspan
//#define CLEAR_WORKER_MARKS   // clears the local mak list after a worker
                                // ended a computation
//#define USE_MPI
#define USE_SSE
//#define USE_OPENCL
//#define USE_THREADS
//#define USE_OPT
//#define SAVE_MORE_MEMORY     // lets the workers compress the subspan
//#define SLIM_INFINITY        // enables the deletion of redundant
                                // representants of the infinite maximum

// Assertion mode:
#ifdef ASSERTION_MODE
    #define ASSERT(condition) { \
        if(!(condition)){ \
            std::cerr << "ASSERT FAILED: " << #condition << " @ " \
                      << __FILE__ << " (" << __LINE__ << ")" \
                      << std::endl; \
        } \
    }
    //#define ASSERT(expr) std::cerr << __FILE__ << " " << __LINE__ << std::endl; assert(expr)
#else
  #define ASSERT(expr)
#endif // ASSERTION_MODE

#define TIME_SECTION(expr) { \
    time_t begin = clock();  \
    expr;                    \
    time_t end   = clock();  \
    std::cout << "Time for " << #expr << " " \
              << double(end-begin)/CLOCKS_PER_SEC << "s\n"; \
}

#ifdef DEBUG_MODE
    #define DEBUG_OUTPUT(expr) std::cout expr
    #define DEBUG_CALL(call) call
    #define ON_DEBUG(expr) expr
#else
    #define DEBUG_OUTPUT(expr)
    #define DEBUG_CALL(call)
    #define ON_DEBUG(expr)
#endif

#endif
