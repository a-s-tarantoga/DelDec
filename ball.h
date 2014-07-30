/*
 * (c) 30.07.2014 Martin Hünniger
 */

#ifndef BALL_H
#define BALL_H

#include "configure.h"
#include "point.h"
#include "subspan.h"

#include <memory>

#include <boost/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#ifdef USE_MPI
#include <boost/mpi.hpp>
#endif

namespace NAMESPACE {

/**
 * @brief The Ball class
 * An instance of this class represents the ball that is defined by a set
 * M of points.
 * To keep track of the balls it has an instance of Subspan<Float>.
 */
template<int dim=3, typename Float=double>
class Ball {
public: // Public types

    typedef Point<dim, Float>       point_type;
    typedef std::vector<point_type> point_list_type;
    typedef std::vector<Float>      float_list_type;
    typedef boost::dynamic_bitset<> hash_value_type;

public: // Construction and destruction

    Ball( const point_list_type& S,
          const point_type& x,
          const int i );
    //    Ball( const point_list_type& S,
    //          const point_type& x,
    //          const int i,
    //          const float_list_type &squared_norm );
    ~Ball();

public: // Copy routines

    Ball( const Ball& copy );
    Ball& operator=( const Ball& other );

public: // Modification

    void add_point( int global_index );    // Add a point
    void add_point_plain( int global_index );    // Add a point without copying the Subspan
    void remove_point( int global_index ); // Remove a point
    void remove_point_plain( int global_index ); // Remove a point without copying the Subspan
    Float update_driver();                 // Compute the new driver
    void update_center( Float t );         // Compute the center for time t
    void update_miniball_center()          // Compute the center of the miniball
    {
        QR->update_miniball_center();
    }
    void compute_center_coefficients( Float* lambdas );
    // Compute the affine coefficients for
    // the center
    void recalculate();                    // Recompute the QR-decomposition by
    // successively inserting the points
    void compute_qr_complete()             // Recompute the QR-decomposition
    {
        QR->compute_qr_complete();
    }

    bool test_criticality() const
    {
        return QR->test_criticality();
    }

    void compress()                        // Compresses the QR-decomposition
    {
        //QR.reset(new Subspan<Float>( *QR.get() ) );
        QR->compress();
    }

    void expand()                          // Restores the QR-decomposition
    {
        QR->expand();
    }

    std::vector<u_long> get_members() {
        return members;
    }

public: // Access

    int size() const
    {
        return r+1;
    }

    bool is_member( int i ) const
    {
        return membership[i];
    }

    int global_index( int i ) const
    {
        return members[i];
    }

    int any_member() const
    {
        return members[r];
    }

    template <class RandomAccessIterator>
    void get_center( RandomAccessIterator p )
    {
        QR->get_center( p );
    }

    template <class RandomAccessIterator>
    void get_driver( RandomAccessIterator p )
    {
        QR->get_driver( p );
    }

    const Float get_radius() const
    {
        return QR->get_radius();
    }

    void print_state() const
    {
        std::cout << "members: ";
        for( int i=0; i<r+1; ++i ) {
            std::cout << members[i] << " ";
        }
        std::cout << "\n";
        std::cout << "membership: ";
        for( int i=0; i<membership.size(); ++i ) {
            std::cout << membership[i] << " ";
        }
        std::cout << "\n";
        std::cout << "r = " << r << "\n";
        std::cout << "is_infinity = " << is_infinity << "\n";
        std::cout << "In QR:\n";
        QR->print_state();
    }

    void print_center() const {
        QR->print_center();
    }

    void print_driver() const {
        QR->print_driver();
    }

    void print_matrices( int i ) const
    // Prints all the matrices that represent the ball
    // 0 - print the center
    // 1 - print center and driver
    // 2 - print Q, R, A, center and driver
    {
        QR->print_matrices( i );
    }

    void print_sorted_members() const
    // Prints the indices of the members of the ball in ascending order
    {
        std::vector<u_long> sorted_members( r+1 );
        std::copy( members.begin(), members.begin()+r+1, sorted_members.begin() );
        std::sort( sorted_members.begin(), sorted_members.end() );
        for( long i=sorted_members.size()-1; i>=0; --i )
            std::cout << sorted_members[i] << " ";
    }

    std::vector<u_long> simplex() const {
        std::vector<u_long> sorted_members(r+1);
        std::copy( members.begin(),
                   members.end(),
                   sorted_members.begin());
        std::sort( sorted_members.begin(),
                   sorted_members.end());
        return sorted_members;
    }

    template<typename RandomAccessIterator1,
             typename RandomAccessIterator2>
    Float time_to_enter_ball( RandomAccessIterator1 p,
                              RandomAccessIterator2 q )
    // Computes the point of time at which the point p enters the ball
    // around the point specified by x in direction v
    {
#ifdef USE_SSE
        return QR->time_to_enter_ball_sse( p, q );
#else
        return QR->time_to_enter_ball_opt( p, q );
#endif
    }

    template<typename RandomAccessIterator>
    Float distance_to_center( RandomAccessIterator p )
    // Compute the distance of the point p to the center
    {
        return QR->distance_to_center( p );
    }

    Float time_to_miniball()
    // Computes the time necessary to reach the miniball from the actual
    // center without catching another point
    {
        return QR->time_to_miniball();
    }

    bool get_infinity()
    // returns true if the ball represents the infinite maximum
    {
        return is_infinity;
    }

    void set_infinity( bool val )
    // guess what?
    {
        is_infinity = val;
    }

public: // Hash value generation

    hash_value_type hash_value();          // Compute a binary identifier
    hash_value_type infinite_hash_value(); // return the hash_value for infty

public: // private helper routines

    void update();

private: // serialization stuff

    friend class boost::serialization::access;

    template<class Archive>
    void serialize( Archive & ar, const unsigned int version );

private: // Member fields

    std::shared_ptr<Subspan<dim,Float>> QR; // the underlying subspan

    const point_list_type& point_list;      // the set points
    std::vector<bool> membership;           // s[i] in the ball iff membership
    std::vector<u_long> members;               // the indices of the points defining
                                            // the ball
    int r;                                  // the size of the base of the ball
    bool up_to_date;                        // Flag to store modification
    bool is_infinity;                       // Tag infinite maximum

    hash_value_type the_hash_value;         // store the hash_value
    bool hash_value_computed;               // Flag to mark its computation

    int point_to_add;                       // the next point erntering the ball
    int point_to_remove;                    // the next point leaving the ball
};

} // FC_NAMESPACE

#ifdef USE_MPI
namespace boost {
namespace mpi {
template<>
template<int dim, typename Float>
struct is_mpi_datatype<NAMESPACE::Ball<dim,Float>> : mpl::true_ { };
}
}
#endif

#include "ball.cpp"

#endif

