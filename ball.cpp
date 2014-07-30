/*
 * (c) 30.07.2014 Martin Hünniger
 */

#ifndef BALL_CPP
#define BALL_CPP

#include "point.h"
#include <algorithm>
#include <boost/serialization/vector.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#ifdef USE_MPI
namespace mpi = boost::mpi;
#endif

namespace NAMESPACE {

template<int dim,typename Float>
Ball<dim,Float>::Ball(const point_list_type& S,
                      const point_type& x,
                      const int index )
    : point_list(S)
    , membership(S.size())
    , members(dim+1)
{
    is_infinity = false;

    //QR.reset(new Subspan<dim,Float>(point_list, x, index));
    QR = std::make_shared<Subspan<dim,Float>>(point_list, x, index);

    point_to_add = -1;
    point_to_remove = -1;

    r = 0;
    members[r] = index;
    membership[index] = true;

    up_to_date = true;
    hash_value();
    hash_value_computed = true;
}


//template<int dim,typename Float>
//Ball<dim,Float>::Ball(const point_list_type& S,
//                      const point_type& x,
//                      const int index
//                      const float_list_type &squared_norm )
//    : point_list(S)
//    , squared_norm(squared_norm)
//    , membership(S.size())
//    , members(dim+1)
//{
//    is_infinity = false;
//    QR.reset( new Subspan<dim,Float>(S, x, index ) );

//    point_to_add = -1;
//    point_to_remove = -1;

//    r = 0;
//    members[r] = index;
//    membership[index] = true;

//    up_to_date = true;
//    hash_value();
//    hash_value_computed = true;
//}

template<int dim,typename Float>
Ball<dim,Float>::~Ball()
{
}

template<int dim,typename Float>
Ball<dim,Float>::Ball( const Ball& copy )
    : point_list(copy.point_list)
    , membership(copy.point_list.size())
    , members(dim+1)
    , r(copy.r)
{
    hash_value_computed = copy.hash_value_computed;;
    is_infinity = copy.is_infinity;

    std::copy( copy.membership.begin(),
               copy.membership.end(),
               membership.begin() );
    std::copy( copy.members.begin(),
               copy.members.end(),
               members.begin() );

    point_to_add    = -1;
    point_to_remove = -1;

    QR = copy.QR;

    hash_value_computed = copy.hash_value_computed;
    the_hash_value      = copy.the_hash_value;
    up_to_date = true;
}

template<int dim,typename Float>
Ball<dim,Float>& Ball<dim,Float>::operator=( const Ball& other )
{
    hash_value_computed = other.hash_value_computed;
    is_infinity = other.is_infinity;

    point_list = other.point_list;
    membership = std::vector<bool>(other.point_list.size());
    members = std::vector<int>(other.dim+1);
    r = other.r;

    std::copy( other.membership.begin(),
               other.membership.end(),
               membership );
    std::copy( other.members.begin(),
               other.members.end(),
               members );

    point_to_add = -1;
    point_to_remove = -1;

    QR  = other.QR;

    hash_value_computed = other.hash_value_computed;
    the_hash_value      = other.the_hash_value;
    up_to_date = true;
}

template<int dim,typename Float>
void Ball<dim,Float>::add_point( int index ) {
    ASSERT( !is_member( index ) );

    hash_value_computed = false;

    membership[index] = true;
    members[r+1] = members[r];
    members[r] = index;
    ++r;
    QR = std::make_shared<Subspan<dim,Float>>( *QR.get() );
    QR->add_point( index );
}

template<int dim,typename Float>
void Ball<dim,Float>::add_point_plain( int index ) {
    ASSERT( !is_member( index ) );

    hash_value_computed = false;

    membership[index] = true;
    members[r+1] = members[r];
    members[r] = index;
    ++r;
    //QR = std::make_shared<Subspan<dim,Float>>( *QR.get() );
    QR->add_point( index );
}

template<int dim,typename Float>
void Ball<dim,Float>::remove_point( int local_index ) {

    ASSERT( is_member( global_index( local_index ) ) );
    ASSERT( size() > 1 );

    //up_to_date = false;
    hash_value_computed = false;
    is_infinity = false;

    membership[global_index( local_index )] = false;
    if( local_index != r ) {
        for( int j = local_index + 1; j<r; ++j ) {
            members[j-1] = members[j];
        }
        members[r-1] = members[r];
    }
    --r;
    QR = std::make_shared<Subspan<dim,Float>>( *QR.get() );
    QR->remove_point( local_index );
}

template<int dim,typename Float>
void Ball<dim,Float>::remove_point_plain( int local_index ) {

    ASSERT( is_member( global_index( local_index ) ) );
    ASSERT( size() > 1 );

    //up_to_date = false;
    hash_value_computed = false;
    is_infinity = false;

    membership[global_index( local_index )] = false;
    if( local_index != r ) {
        for( int j = local_index + 1; j<r; ++j ) {
            members[j-1] = members[j];
        }
        members[r-1] = members[r];
    }
    --r;
    //QR = std::make_shared<Subspan<dim,Float>>( *QR.get() );
    QR->remove_point( local_index );
}

template<int dim,typename Float>
Float Ball<dim,Float>::update_driver() {
    ASSERT( QR != NULL );

    //if( !up_to_date )
    //update();

    return QR->update_driver();
}

template<int dim, typename Float>
void Ball<dim,Float>::update_center( Float t ) {
    // We don't need to have the QR updated, since if the ball
    // is grown to t it will mit the newly added point. !!
    ASSERT( QR != NULL );

    //if( point_to_add > -1 )
    //update();

    QR->update_center( t );
}

template<int dim,typename Float>
void Ball<dim,Float>::compute_center_coefficients( Float* lambdas ) {
    ASSERT( up_to_date );
    ASSERT( QR != NULL );
    QR->compute_center_coefficients( lambdas );
}

template<int dim,typename Float>
void Ball<dim,Float>::recalculate() {
    Point<dim,Float> x;
    QR = std::make_shared<Subspan<dim,Float>>(point_list, x, members[r]);
    for( int i=0; i<r; ++i ) QR->add_point( members[i] );
    QR->update_miniball_center();
}

template<int dim, typename Float>
void Ball<dim,Float>::update() {

    if( false ) {
        point_to_remove = -1;
        point_to_add = -1;
        up_to_date = true;
        return;
    }
    QR = std::make_shared<Subspan<dim,Float>>( *QR.get() );

    if( point_to_remove > -1 ) QR->remove_point( point_to_remove );

    point_to_remove = -1;

    if( point_to_add > -1 ) QR->add_point( point_to_add );

    point_to_add = -1;

    up_to_date = true;
}

template<int dim,typename Float>
template<class Archive>
void Ball<dim,Float>::serialize( Archive & ar, const unsigned int )
// Store all the relevant members of the ball in Archive. This serialization
// is needed to transfer the data over a network as used by Boost.MPI.
{
    ar & BOOST_SERIALIZATION_NVP( membership );
    ar & BOOST_SERIALIZATION_NVP( members );
    ar & BOOST_SERIALIZATION_NVP( r );
    ar & BOOST_SERIALIZATION_NVP( is_infinity );
    ar & BOOST_SERIALIZATION_NVP( up_to_date );
    ar & BOOST_SERIALIZATION_NVP( *QR );

    /*
    for( int i=0; i<members.size(); ++i ) {
      DEBUG_OUTPUT( << members[i] << " " );
    }
    DEBUG_OUTPUT( << std::endl );
    for( int i=0; i<membership.size(); ++i ) {
      DEBUG_OUTPUT( << membership[i] << " " );
    }
    DEBUG_OUTPUT( << std::endl );
    DEBUG_OUTPUT( << r << std::endl );
    DEBUG_OUTPUT( << is_infinity << std::endl );
    */
}



template<int dim,typename Float>
auto Ball<dim,Float>::infinite_hash_value() -> hash_value_type {

    unsigned long length = 0;
    unsigned long two_to_the_length = 1;
    while( membership.size()-1 >= two_to_the_length ) {
        ++length;
        two_to_the_length *= 2;
    }

    typename Ball<dim,Float>::hash_value_type res((dim+1)*length);

    for( u_long i=0; i<(dim+1)*length; ++i ) {
        res[i] = 1;
    }
    return res;
}

template<int dim,typename Float>
auto Ball<dim,Float>::hash_value() -> hash_value_type {
    return QR->hash_value();
}

} // NAMESPACE

#endif
