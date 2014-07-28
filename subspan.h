/* class representing the affine hull of a point set
 *
 * ripped from Martin Kutz and Kaspar Fischer
 *
 * Martin Huenniger 07.04.2010
 */

#ifndef SUBSPAN_H
#define SUBSPAN_H

#include "point.h"
#include "helper.h"
#include "marks.h"

#include <boost/dynamic_bitset.hpp>
#include <vector>

#ifdef USE_MPI
#include <boost/mpi.hpp>
#endif

namespace NAMESPACE {

/**
 * @brief Subspan
 * An instance of this class represents the affine hull of a non-empty
 * set M of affinely independent points. The set M is not represented
 * explicitly; when an instance of this class is created, it takes a list
 * of points S to it (which the instance will never change and which is
 * assumed to stay fixed during the lifetime of this instance): The set M
 * is the a subset of S, and its members are indentified by their (zero-
 * based) indices in S.
 */
template<int dim=3,typename Float=double>
class Subspan
{
private: // types

public: // public types

    typedef Point<dim,Float>                    point_type;
    typedef std::vector<point_type>             point_list_type;
    typedef typename Marks::mark_type           hash_value_type;
//    typedef boost::dynamic_bitset<>           Hash_value_type;
    typedef typename point_list_type::size_type size_type;

public: // Construction and destruction
    Subspan( const point_list_type &S,
             const point_type &x,
             const size_type i );
    ~Subspan();

public: // Copy routines

    Subspan( const Subspan& copy );
    Subspan& operator=( const Subspan& other );

public: // Modification

    void add_point( size_type global_index );   // inserts a new point
    void remove_point( size_type global_index); // removes a point
    Float update_driver();                      // computes the driver of the balls center x
                                                // and returns their
                                                // squared distance
    void update_center( Float t );              // computes the new center of the ball
    void update_miniball_center();              // compute the center directly
    void compute_center_coefficients( Float *lambdas );
                                                // computes the affine coordinates of
                                                // the center in the ball
                                                // lambdas is (r+1)-dimensional
    void compute_qr_complete();                 // compute the whole qr-decomposition.

public: // access

    size_type size() const
    {
        return r+1;
    }

    bool is_member( size_type i ) const
    {
        ASSERT( i >= 0 && i < point_list.size() );
        return membership[i];
    }

    size_type global_index( size_type i ) const
    {
        ASSERT( i >= 0 && i < size() );
        return members[i];
    }

    size_type any_member() const
    {
        ASSERT( size() > 0 );
        return members[r];
    }

    template <class RandomAccessIterator>
    void get_center( RandomAccessIterator p ) const
    {
        for( size_type i=0; i<dim; ++i )
            *(p+i) = x[i];
    }

    template <class RandomAccessIterator>
    void get_driver( RandomAccessIterator p ) const
    {
        for( size_type i=0; i<dim; ++i )
            *(p+i) = d[i];
    }

    template <class RandomAccessIterator>
    void get_center( RandomAccessIterator p, Float t ) const
    {
        for( size_type i=0; i<dim; ++i )
            *(p+i) = x[i]+t*d[i];
    }

    const Float get_radius() const
    {
        return distance_to_center( point_list[any_member()].begin() );
    }

    void print_state() const
    {
        std::cout << "members:    " << members << std::endl;
        std::cout << "membership: " << membership << std::endl;
        std::cout << "r =         " << r << "\n";
    }

    /**
     * @brief print_matrices
     * Debug helper function that prints out the internal state of the
     * qr-decomposition Q and R are storing the base vectors as columns
     * instead of rows. This might be necessary to fulfill the requirement
     * m>=n for the represented matrix A = QR, but it leads to confusion
     * when inspecting the values of this matrices.
     * @param i loglevel: 0,1,2
     */
    void print_matrices( int i ) const;

    void print_center() const {
        std::cout << "center: [ ";
        for(auto i=0; i<dim; ++i) {
            std::cout << x[i] << " ";
        }
        std::cout << "]\n";
    }

    template<typename RandomAccessIterator1,
             typename RandomAccessIterator2>
    const Float project_point_onto_span( RandomAccessIterator1 p,
                                         RandomAccessIterator2 w ) const;
    // Computes the projection of the point on the span
    // we need the vector from x to the span (orthogonal projection
    // and compute the corresponding formula (computes the driver...)

    template<typename RandomAccessIterator1,
             typename RandomAccessIterator2>
    void find_affine_coefficients( RandomAccessIterator1 x,
                                   RandomAccessIterator2 coeffs );

    Float time_to_enter_ball( size_type i, size_type j ) const;
    // Computes the point of time at which the point p enters the ball
    // around the point specified by x in direction v
    Float time_to_enter_ball_opt( size_type i, size_type j ) const;
    Float time_to_enter_ball_sse( size_type i, size_type j ) const;

    template<typename RandomAccessIterator>
    const Float distance_to_center( RandomAccessIterator p ) const;
    // Computes the distance of point p to center x

    const Float time_to_miniball() const
    // Computes the time it takes to reach the miniball from the actual
    // center without catching another ball
    {
        return sqrt( project_point_onto_span( x, w ) );
    }

public: // Hashvalue generation

    hash_value_type hash_value() const;          // Computes a binary identifier
    hash_value_type infinite_hash_value() const; // return the hash_value for infty

private: // private helper routines

    void append_column();
    // Append the new column u to the right of A it assumes r to be still the old
    // value, i.e. the index of the colum used now for insertion. r is not alterd
    // and should be changed by the caller afterwards.
    // Precondition: r<dim

    void hessenberg_clear( size_type start );
    // Given R in lower Hessenberg form with subdiagonal entries 0 top pos-1
    // already set to 0, clears all the remaining subdiagonal entries via
    // Givens rotations

    void special_update();
    // Update current A to A + u [1,...,1]

private: // Serialization stuff...
    friend class boost::serialization::access;

    template<class Archive>
    void serialize( Archive & ar, const unsigned int version );

private: // member fields
    const point_list_type& point_list;    // a const-reference to the set
    std::vector<bool> membership;         // point_list[i] in M iff membership
    std::vector<size_type> members;       // Entry i of members contains the index
                                          // into S of the i-th point in ;
                                          // members[r] is the origin
                                          //    const float_list_type &squared_norm;
                                          // contains the squared norms of all the
                                          // input points

private: // memberfields for maintaining the affine space
    Float **Q, **R;                       // (dim x dim)-matrices representing the
                                          // affine space as a QR decomposition.
    Float *u,*w;                          // vector from x to driver and helper
    Float *x,*d;                          // x the actual point, d the driver of x
    size_type r;                          // rank of A (#points - 1)
};

} // NAMESPACE

#ifdef USE_MPI
namespace boost {
    namespace mpi {
        template<>
        template<int dim, typename Float>
        struct is_mpi_datatype<NAMESPACE::Subspan<dim,Float>> : mpl::true_ { };
    }
}
#endif

#include "subspan.cpp"

#endif
