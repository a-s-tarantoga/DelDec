/*
 * (c) 30.07.2014 Martin Hünniger
 */
#ifndef POINT_TEST_H
#define POINT_TEST_H

#include "point.h"
#include <vector>
#include <cstdlib>
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE( point_test )
{
    typedef DD::Point<2,double> point_type_2D;
    typedef DD::Point<3,double> point_type_3D;
    typedef DD::Point<5,double> point_type_5D;

    point_type_2D p1({1.0, 0.0});
    point_type_2D p2({0.0, 1.0});
    point_type_2D p3({1.0, 1.0});
    point_type_2D p4;
    point_type_5D p5({-1.0, -2.0, -3.0, -4.0, -5.0});

    BOOST_CHECK( p1[0] == 1.0 && p1[1] == 0.0 );
    BOOST_CHECK( p2[0] == 0.0 && p2[1] == 1.0 );
    BOOST_CHECK( p5[0] == -1.0 &&
                 p5[1] == -2.0 &&
                 p5[2] == -3.0 &&
                 p5[3] == -4.0 &&
                 p5[4] == -5.0 );
    BOOST_CHECK( p1.squared_norm() ==  1.0 );
    BOOST_CHECK( p2.squared_norm() ==  1.0 );
    BOOST_CHECK( p3.squared_norm() ==  2.0 );
    BOOST_CHECK( p5.squared_norm() == 55.0 );
    p4 = std::move(p3);
    BOOST_CHECK( p4[0] == 1.0 && p4[1] == 1.0 );
    p3 = p4;
    BOOST_CHECK( p4[0] == p3[0] && p4[1] == p3[1] );
    point_type_2D p6(p4);
    BOOST_CHECK( p6[0] == p3[0] && p6[1] == p3[1] );
    std::vector<double> v({2.0, 3.0, 4.0});
    point_type_3D p7(v.begin());
    BOOST_CHECK( p7[0] == v[0] && p7[1] == v[1] && p7[2] == v[2] );
    srand(clock());
    p7.perturb(0.001);
    BOOST_CHECK_CLOSE(p7[0], v[0], 0.1);
    BOOST_CHECK_CLOSE(p7[1], v[1], 0.1);
    BOOST_CHECK_CLOSE(p7[2], v[2], 0.1);
}
/*
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( point_failure_test, 1 )

BOOST_AUTO_TEST_CASE( point_failure_test ) {

    BOOST_CHECK();
}
*/
#endif // POINT_TEST_H
