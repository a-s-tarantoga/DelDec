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

#ifndef FC_POINT_H
#define FC_POINT_H

#include "configure.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <initializer_list>

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <strstream>
#include <sstream>
#include <iterator>


#include <boost/serialization/vector.hpp>
#ifdef USE_MPI
#include <boost/mpi.hpp>
#endif

namespace NAMESPACE {

template<int dim=3, typename Float=double>
class Point
{
public: // types
    typedef typename std::vector<Float>::const_iterator Const_iterator;
    typedef typename std::vector<Float>::iterator       Iterator;
    typedef typename std::vector<Float>::size_type      size_type;

public: // construction and destruction

    /**
     * @brief Point Constructs a d-dimensional point with zero coordinates
     */
    Point()
        : c(dim)
        , flag(false)
    {}

    /**
     * @brief Point constructs the point from a initializer list
     * @param list the initializer list
     */
    Point(const std::initializer_list<Float> list)
        : c(dim)
        , flag(false)
    {
        ASSERT(list.size() == dim);
        std::copy(list.begin(), list.end(), c.begin());
    }

    /**
     * @brief Point Constructs a d-dimensional point with Cartesian center coordinates [first,first+d).
     * @param d the dimension of the point
     * @param first Iterator to the first data element used to initialize the point.
     */
    template<typename InputIterator>
    Point(const InputIterator first )
        : c(first,first+dim)
        , flag(false)
    {}

    Point(Point && other)
       : c(std::move(other.c))
       , flag(std::move(other.flag))
       , squared_norm_value(std::move(other.squared_norm_value))
    {}

    Point(const Point& other)
        : c(other.begin(),other.end())
        , flag(other.flag)
        , squared_norm_value(other.squared_norm_value)
    {}

    Point& operator=(Point other)
    {
        std::swap(c,other.c);
        std::swap(flag,other.flag);
        std::swap(squared_norm_value,other.squared_norm_value);
        return *this;
    }

public: // access

    /**
     * @brief operator [] returns a const-reference to the i-th coordinate.
     * @param i the index of the coordinate to return
     * @return the value at the i-th coordinate
     */
    const Float& operator[]( const size_type i ) const
    {
        ASSERT( 0 <= i && i < c.size() );
        return c[i];
    }
    Float& operator[]( const size_type i )
    {
        ASSERT( 0 <= i && i < c.size() );
        return c[i];
    }

    const Float* data() const { return c.data(); }
    Float* pointer() const { return const_cast<Float*>( &c.at(0) ); }

    Const_iterator begin() const { return c.begin(); }
    Const_iterator end() const   { return c.end(); }
    Iterator begin() { return c.begin(); }
    Iterator end()   {  return c.end(); }

    size_type size() const { return c.size(); }

    /**
     * @brief perturb Perturbs the point in a ball with radius max
     * remember to initialize the random generator!!
     * @param max
     */
    void perturb( const Float max )
    {
        typename std::vector<Float>::iterator it;
        for( it = c.begin(); it != c.end(); ++it ) {
            *it += static_cast<Float>( 2*max*rand()/RAND_MAX-max );
        }
        flag = false;
    }

    const Float squared_norm() const
    {
        if(!flag) {
            squared_norm_value = std::inner_product(c.begin(),c.end(),c.begin(),Float(0.0));
            flag = true;
        }
        return squared_norm_value;
    }

private: // Serialization

    friend class boost::serialization::access;

    template<class Archive>
    void serialize( Archive &ar, const unsigned int )
    {
        ar & BOOST_SERIALIZATION_NVP( c );
    }


    friend std::ifstream& operator>>(std::ifstream& ifs, Point& point) {
        std::string str;
        std::getline(ifs, str);
        std::stringstream ss(str);
        std::copy(std::istream_iterator<Float>(ss),
                  std::istream_iterator<Float>(),
                  point.begin());
        return ifs;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& point ) {
        os << "[ ";
        for(auto it=point.c.begin();it!=point.c.end();++it) {
            os << *it << " ";
        }
        os << "]";
        return os;
    }

    friend std::ofstream& operator<<(std::ofstream& ofs, const Point& point ) {
        for(auto it=point.c.begin();it!=point.c.end();++it) {
            ofs << *it << " ";
        }
        return ofs;
    }

private: // member fields
    std::vector<Float> c;          // Cartesian center coordinates.
    mutable Float squared_norm_value;
    mutable bool flag;
};

} // namespace FC_NAMESPACE

#ifdef USE_MPI
namespace boost {
namespace mpi {
template<>
template<int dim, typename Float>
struct is_mpi_datatype<NAMESPACE::Point<dim,Float> > : mpl::true_ { };
}
}
#endif

#endif
