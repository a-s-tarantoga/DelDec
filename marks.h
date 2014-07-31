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
#ifndef MARKS_H
#define MARKS_H

#include "ball.h"
#include <vector>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

namespace NAMESPACE {

struct Marks {

    typedef boost::dynamic_bitset<> mark_type;

    static mark_type compute(const std::vector<unsigned long> &vec, int max_size) {

        std::vector<long> sorted_vec(vec.size());
        std::copy(vec.begin(),
                  vec.end(),
                  sorted_vec.begin());
        std::sort(sorted_vec.begin(),sorted_vec.end());
        int log_max_size = 0;
        int c = 1;
        while(max_size-1 >= c) {
            ++log_max_size; c *= 2;
        }
        mark_type res(vec.size()*log_max_size);
        for(auto i=0;i<sorted_vec.size();++i ) {
            mark_type bits( log_max_size, sorted_vec[i] );
            for( auto j=0;j<log_max_size;++j) {
                res[log_max_size*i+j] = bits[j];
            }
        }
        return res;
    }

    static mark_type compute_without_index(const std::vector<unsigned long> &vec, unsigned long index, int max_size) {
        int log_max_size = 0;
        int c = 1;
        while(c <= max_size-1) {
            ++log_max_size; c *= 2;
        }
#ifndef BLUBB
        unsigned long sorted[vec.size()-1];
        int i = 0;
        for( unsigned long el : vec ) {
            if(el != index) {
                sorted[i] = el;
                ++i;
            }
        }
        std::sort(sorted, sorted+vec.size()-1);
        mark_type res(vec.size()*log_max_size);
        for(auto k=0;k<vec.size()-1;++k ) {
            mark_type bits( log_max_size, sorted[k] );
            for( auto j=0;j<log_max_size;++j) {
                res[log_max_size*k+j] = bits[j];
            }
        }
        return res;
#else
        std::vector<long> sorted_vec(0);
        for( auto el : vec) {
            if(el != index) {
                sorted_vec.push_back(el);
            }
        }
        std::sort(sorted_vec.begin(),sorted_vec.end());
        mark_type res(vec.size()*log_max_size);
        for(auto i=0;i<sorted_vec.size();++i ) {
            mark_type bits( log_max_size, sorted_vec[i] );
            for( auto j=0;j<log_max_size;++j) {
                res[log_max_size*i+j] = bits[j];
            }
        }
        return res;
#endif
    }

    static mark_type compute_exceptional(int dim, int max_size) {
        int log_max_size = 0;
        int c = 1;
        while(max_size-1 >= log_max_size) {
            ++log_max_size; c *= 2;
        }
        mark_type res((dim+1)*log_max_size);
        for(auto i=0;i<(dim+1)*log_max_size;++i ) {
            res[i] = 1;
        }
        return res;
    }

    struct Cmp {

        bool operator()( const mark_type& a,
                         const mark_type& b ) const
        {
            if( a.size() < b.size() ) return true;
            else if( a.size() > b.size() ) return false;

            for( int i = a.size()-1; i >= 0; --i ) {
                if( a[i] < b[i] ) return true;
                else if( a[i] > b[i] ) return false;
            }
            return false;
        }
    };

};
} // NAMESPACE

#endif // MARKS_H
