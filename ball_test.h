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
#ifndef BALL_TEST_H
#define BALL_TEST_H

#include "ball.h"
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE( ball_test )
{
    typedef DD::Ball<2,double>         ball_type;
    typedef ball_type::point_type      point_type;
    typedef ball_type::point_list_type point_list_type;
    typedef ball_type::float_list_type float_list_type;

    point_type p1({1.0, 0.0});
    point_type p2({0.0, 1.0});
    point_type p3({1.0, 1.0});

    point_list_type point_list({p1, p2, p3});

    ball_type test_ball( point_list,
                         point_list[0],
                         0 );

    BOOST_CHECK( test_ball.size() == 1 );
}

#endif // BALL_TEST_H
