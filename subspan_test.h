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
#ifndef SUBSPAN_TEST_H
#define SUBSPAN_TEST_H

#include "subspan.h"
#include <vector>
#include <cstdlib>
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;


template<int dim>
void test_function() {

    typedef DD::Subspan<dim,double>                subspan_type;
    typedef typename subspan_type::point_type      point_type;
    typedef typename subspan_type::point_list_type point_list_type;

    point_list_type list;
    for(int i=0;i<dim;++i) {
        std::vector<double> v(dim);
        v[i] = 1.0;
        point_type p(v.begin());
        list.push_back(p);
    }
    std::vector<double> v(dim);
    for(int i=0;i<dim;++i) {
        v[i] = 1.0;
    }
    point_type p(v.begin());
    list.push_back(p);

    subspan_type s(list, list[0], 0);

    BOOST_CHECK(s.size() == 1);
    BOOST_CHECK(s.global_index(0) == 0);
    BOOST_CHECK(s.any_member() == 0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        BOOST_CHECK_CLOSE(lambdas[0],1,0.001);
    }

    s.add_point(1);

    BOOST_CHECK(s.size() == 2);
    BOOST_CHECK(s.global_index(0) == 1);
    BOOST_CHECK(s.global_index(1) == 0);
    BOOST_CHECK(s.any_member() == 0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        BOOST_CHECK_CLOSE(lambdas[0],0.5,0.001);
        BOOST_CHECK_CLOSE(lambdas[1],0.5,0.001);
    }

    s.add_point(2);

    BOOST_CHECK(s.size() == 3);
    BOOST_CHECK(s.global_index(0) == 1);
    BOOST_CHECK(s.global_index(1) == 2);
    BOOST_CHECK(s.global_index(2) == 0);
    BOOST_CHECK(s.any_member() == 0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        if(dim == 2) {
            BOOST_CHECK_CLOSE(lambdas[0],0.5,0.001);
            BOOST_CHECK_CLOSE(lambdas[1],0.0,0.001);
            BOOST_CHECK_CLOSE(lambdas[2],0.5,0.001);
        } else {
            BOOST_CHECK_CLOSE(lambdas[0],1.0/3,0.001);
            BOOST_CHECK_CLOSE(lambdas[1],1.0/3,0.001);
            BOOST_CHECK_CLOSE(lambdas[2],1.0/3,0.001);
        }
    }

    if(dim == 2) goto remove;

    s.add_point(3);

    BOOST_CHECK(s.size() == 4);
    BOOST_CHECK(s.global_index(0) == 1);
    BOOST_CHECK(s.global_index(1) == 2);
    BOOST_CHECK(s.global_index(2) == 3);
    BOOST_CHECK(s.global_index(3) == 0);
    BOOST_CHECK(s.any_member() == 0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        BOOST_CHECK_CLOSE(lambdas[0],0.25,0.001);
        BOOST_CHECK_CLOSE(lambdas[1],0.25,0.001);
        BOOST_CHECK_CLOSE(lambdas[2],0.25,0.001);
        BOOST_CHECK_CLOSE(lambdas[3],0.25,0.001);
    }

    if(dim == 3) goto remove;

    s.add_point(4);

    BOOST_CHECK(s.size() == 5);
    BOOST_CHECK(s.global_index(0) == 1);
    BOOST_CHECK(s.global_index(1) == 2);
    BOOST_CHECK(s.global_index(2) == 3);
    BOOST_CHECK(s.global_index(3) == 4);
    BOOST_CHECK(s.global_index(4) == 0);
    BOOST_CHECK(s.any_member() == 0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        BOOST_CHECK_CLOSE(lambdas[0],0.2,0.001);
        BOOST_CHECK_CLOSE(lambdas[1],0.2,0.001);
        BOOST_CHECK_CLOSE(lambdas[2],0.2,0.001);
        BOOST_CHECK_CLOSE(lambdas[3],0.2,0.001);
        BOOST_CHECK_CLOSE(lambdas[4],0.2,0.001);
    }

    s.add_point(5);

    BOOST_CHECK(s.size() == 6);
    BOOST_CHECK(s.global_index(0) == 1);
    BOOST_CHECK(s.global_index(1) == 2);
    BOOST_CHECK(s.global_index(2) == 3);
    BOOST_CHECK(s.global_index(3) == 4);
    BOOST_CHECK(s.global_index(4) == 5);
    BOOST_CHECK(s.global_index(5) == 0);
    BOOST_CHECK(s.any_member() == 0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        BOOST_CHECK_CLOSE(lambdas[0],0.125,0.001);
        BOOST_CHECK_CLOSE(lambdas[1],0.125,0.001);
        BOOST_CHECK_CLOSE(lambdas[2],0.125,0.001);
        BOOST_CHECK_CLOSE(lambdas[3],0.125,0.001);
        BOOST_CHECK_CLOSE(lambdas[4],0.375,0.001);
        BOOST_CHECK_CLOSE(lambdas[5],0.125,0.001);
    }

remove:

    s.remove_point(0);

    s.update_miniball_center();
    {
        double lambdas[s.size()];
        s.compute_center_coefficients(lambdas);
        switch(dim) {
        case 5:
            BOOST_CHECK_CLOSE(lambdas[0],2.0/13,0.001);
            BOOST_CHECK_CLOSE(lambdas[1],2.0/13,0.001);
            BOOST_CHECK_CLOSE(lambdas[2],2.0/13,0.001);
            BOOST_CHECK_CLOSE(lambdas[3],5.0/13,0.001);
            BOOST_CHECK_CLOSE(lambdas[4],2.0/13,0.001);
            break;
        case 3:
            BOOST_CHECK_CLOSE(lambdas[0],1.0/3,0.001);
            BOOST_CHECK_CLOSE(lambdas[1],1.0/3,0.001);
            BOOST_CHECK_CLOSE(lambdas[2],1.0/3,0.001);
            break;
        case 2:
            BOOST_CHECK_CLOSE(lambdas[0],0.5,0.001);
            BOOST_CHECK_CLOSE(lambdas[1],0.5,0.001);
        }
    }
}

BOOST_AUTO_TEST_CASE( subspan_test )
{
    test_function<2>();
    test_function<3>();
    test_function<5>();
}

#endif // SUBSPAN_TEST_H
