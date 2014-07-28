/* class representing the affine hull of a point set
 *
 * ripped from Martin Kutz and Kaspar Fischer
 *
 * Martin Huenniger 07.04.2010
 */

#ifndef SUBSPAN_CPP
#define SUBSPAN_CPP

#include <numeric>
#include <algorithm>
#include <cmath>

#include "helper.h"

// Serialization support for Boost.MPI
#include <boost/serialization/vector.hpp>

#ifdef USE_MPI
namespace mpi = boost::mpi;
#endif

#define AFFINE_ORIGIN point_list[members[r]]

namespace NAMESPACE {

/**
   * @brief givens
   * Determine the Givens coefficients (c,s) satisfying
   *
   *   c*a + s*b = +/-( a^2 + b^2 )
   *   c*b - s*b = 0
   *
   * We don't care about the signs here for efficiency,
   * so make sure not to rely on them anywhere.
   *
   * Update: We now care about the signs and hope, the
   * copysign function won't crash the template stuff...
   *
   * Source: Golub et.al. "Matrix Computations"
   *
   * @param c
   * @param s
   * @param a
   * @param b
   */
template<typename Float>
void givens( Float& c, Float& s, const Float a, const Float b )
{
    using std::abs;
    using std::sqrt;
    //using std::copysign;

    if( b == 0.0 ) {
        c=1.0;
        s=0.0;
    } else if( abs(b) > abs(a) ) {
        const Float t = a/b;
        const Float u = copysign( 1.0, b );
        s = u / sqrt(1.0 + t*t);
        c = s*t;
    } else {
        const Float t = b/a;
        const Float u = copysign( 1.0, a );
        c = u / sqrt(1.0 + t*t);
        s = c*t;
    }
}

template<int dim, typename Float>
Subspan<dim,Float>::Subspan( const point_list_type &S,
                             const point_type &p,
                             const size_type index )
    : point_list(S)
    , membership(point_list.size())
    , members(dim+1)
{
    // allocate storage for Q,R and the vectors
    Q = new Float *[dim];
    R = new Float *[dim];

    for( size_type i=0; i<dim; ++i ) {
        Q[i] = new Float[dim];
        R[i] = new Float[dim];
    }

    w = new Float[dim];  // helper
    u = new Float[dim];  // helper
    x = new Float[dim];  // center
    d = new Float[dim];  // driver

    // initialize Q to the identity matrix
    for( size_type i=0; i<dim; ++i )
        for( size_type j=0; j<dim; ++j )
            Q[i][j] = (i==j)? 1.0 : 0.0;

    for( size_type i=0; i<dim; ++i )
        for( size_type j=0; j<dim; ++j )
            R[i][j] = 0.0;

    for( size_type i=0; i<dim; ++i )
        x[i] = p[i];

    r = 0;
    members[r] = index;
    membership[index] = true;
}

template<int dim, typename Float>
Subspan<dim, Float>::~Subspan()
{
    for( size_type i=0; i<dim; ++i ) {
        delete[] Q[i];
        delete[] R[i];
    }
    delete[] Q;
    delete[] R;

    delete[] w;
    delete[] u;
    delete[] x;
    delete[] d;
}

template<int dim, typename Float>
Subspan<dim,Float>::Subspan(const Subspan<dim,Float>& copy )
    : point_list(copy.point_list)
    , membership(copy.point_list.size())
    , members(dim+1)
    , d(copy.d)
    , r(copy.r)
{
    std::copy( copy.membership.begin(),
               copy.membership.end(),
               membership.begin() );

    std::copy( copy.members.begin(),
               copy.members.end(),
               members.begin() );

    // create the new matrices
    Q = new Float *[dim];
    R = new Float *[dim];
    for( size_type i=0; i<dim; ++i ) {
        Q[i] = new Float[dim];
        R[i] = new Float[dim];
    }

    // create the new vectors
    w = new Float[dim];  // helper
    u = new Float[dim];  // helper
    x = new Float[dim];  // center
    d = new Float[dim];  // driver

    // copy the matrices
    for( size_type j=0; j<dim; ++j )
        for( size_type i=0; i<dim; ++i ) {
            Q[j][i] = copy.Q[j][i];
            R[j][i] = copy.R[j][i];
        }

    // copy the vectors
    for( size_type i=0; i<dim; ++i ) {
        w[i] = copy.w[i];
        u[i] = copy.u[i];
        x[i] = copy.x[i];
        d[i] = copy.d[i];
    }
}

template<int dim,typename Float>
Subspan<dim,Float>& Subspan<dim,Float>::operator=( const Subspan<dim,Float>& other ){

    // Initialize all the stuff
    point_list = other.point_list;
    membership = std::vector<bool>(other.S.size());
    members = std::vector<int>(other.dim+1);
    r = other.r;

    // copy all the stuff
    std::copy( other.membership.begin(),
               other.membership.end(),
               membership );
    std::copy( other.members.begin(),
               other.members.end(),
               members );

    // copy the matrices
    for( int j=0; j<dim; ++j )
        for( int i=0; i<dim; ++i ) {
            Q[j][i] = other.Q[j][i];
            R[j][i] = other.R[j][i];
        }

    // copy the vectors
    for( int i=0; i<dim; ++i ) {
        w[i] = other.w[i];
        u[i] = other.u[i];
        x[i] = other.x[i];
        d[i] = other.d[i];
    }
}

template<int dim,typename Float>
void Subspan<dim,Float>::add_point( size_type index ) {
    ASSERT( !is_member( index ) );

    using std::inner_product;

    //compute S[i] - origin in u
    for( size_type i=0; i<dim; ++i ) {
        u[i] = point_list[index][i] - AFFINE_ORIGIN[i];
    }

    // appends new column u to R and updates QR-decomposition
    // routine works with old r:
    append_column();

    // move origin index and insert new index
    membership[index] = true;
    members[r+1] = members[r];
    members[r] = index;

    ++r;
}

template<int dim,typename Float>
void Subspan<dim,Float>::remove_point( const size_type local_index ) {
    ASSERT( is_member( global_index( local_index ) ) );
    ASSERT( size() > 1 );

    membership[global_index( local_index )] = false;

    if( local_index == r ) {
        // origin is to be deleted

        // choose the right-most entry as new origin
        // all the relative vectors have to be updated
        // by u := old origin - S[global_index(r-1)]
        for( size_type i=0; i<dim; ++i ) {
            u[i] = AFFINE_ORIGIN[i] - point_list[global_index(r-1)][i];
        }
        --r;
        special_update();
    } else {
        // general case: delete column from A = QR;

        // shift higher columns of R and T one step to the left
        Float *dummy = R[local_index];
        for( size_type j = local_index + 1; j<r; ++j ) {
            R[j-1] = R[j];
            members[j-1] = members[j];
        }
        members[r-1] = members[r];      // shift down origin
        R[--r] = dummy;                 // relink trash column

        // zero out all subdiagonal columns
        hessenberg_clear( local_index );
    }
}

template<int dim,typename Float>
Float Subspan<dim,Float>::update_driver() {
    ON_DEBUG(std::cout << "Old driver: [ ";
            for(auto i=0; i<dim;++i)
                std::cout << d[i] << " ";
            std::cout << "]\n");

    Float length_squared = project_point_onto_span( x, d );

    ON_DEBUG(std::cout << "New driver: [ ";
            for(auto i=0; i<dim;++i)
                std::cout << d[i] << " ";
            std::cout << "]\n");

    Float length = sqrt( length_squared );
    for( size_type i=0; i<dim; ++i )
        d[i] = d[i] / length;

    ON_DEBUG(std::cout << "New driver (normalized): [ ";
            for(auto i=0; i<dim;++i)
                std::cout << d[i] << " ";
            std::cout << "]\n");

    return length_squared;

}

template<int dim,typename Float>
void Subspan<dim,Float>::update_center( Float time ) {

    for( size_type i=0; i<dim; ++i )
        x[i] += time*d[i];
}

/**
 * @brief update_miniball_center
 * Solve the linear equation 2*(R^T * R) * w = c by using forward substitution
 * to compute R^T * w' = c/2 and then using backward substitution to compute
 * R * w = w'
 * This method is correct even if r <= dim, since we compute the miniball
 * _inside_ the affine space given by our qr-decomposition.
 */
template<int dim,typename Float>
void Subspan<dim,Float>::update_miniball_center()
{
    ASSERT( r <= dim );

    ON_DEBUG(std::cout << "== computing new center ==\n");

    for( size_type j=0; j<r; ++j ) {
        for( size_type i=0; i<dim; ++i ) {
            u[i] = point_list[global_index(j)][i] - AFFINE_ORIGIN[i];
        }
        w[j] = std::inner_product( u, u+dim, u, Float(0.0) )/2.0;
    }

    // forward substitution
    for( size_type k=0; k<r; ++k ) {
        for( size_type i=0; i < k; ++i ) {
            w[k] -= R[k][i] * w[i];        // we're acting on R^T
        }
        w[k] = w[k] / R[k][k];
    }

    // backward substitution
    for( int k = int(r)-1; k>=0; --k ) {
        for( size_type i = k+1; i<r; ++i )
            w[k] -= R[i][k] * w[i];       // now we're acting on R
        w[k] = w[k] / R[k][k];
    }

    // compute the new center of the ball in global coordinates
    for( size_type i=0; i<dim; ++i ) {
        u[i] = 0.0;
        for( size_type j=0; j<r; ++j ) {
            u[i] +=  w[j] * (point_list[global_index(j)][i] - AFFINE_ORIGIN[i]);
        }
        u[i] = u[i] + AFFINE_ORIGIN[i];
    }

    for( size_type i=0; i<dim; ++i ) {
        x[i] = u[i];
    }
}

template<int dim,typename Float>
void Subspan<dim,Float>::compute_center_coefficients( Float *lambdas ) {
    find_affine_coefficients( x, lambdas );
}


template<int dim,typename Float>
template<typename RandomAccessIterator1,
         typename RandomAccessIterator2>
const Float Subspan<dim,Float>::project_point_onto_span( RandomAccessIterator1 p,
                                                         RandomAccessIterator2 w ) const
{
    using std::inner_product;

    // compute vector from p to origin
    for( size_type i=0; i<dim; ++i )
        w[i] = AFFINE_ORIGIN[i] - p[i];

    // remove projections of w onto the affine hull
    for( size_type j=0; j<r; ++j ) {
        const Float scale = inner_product( w, w+dim, Q[j], Float(0.0) );
        for( size_type i=0; i<dim; ++i )
            w[i] -= scale * Q[j][i];
    }

    return inner_product( w, w+dim, w, Float(0.0) );
}

template<int dim,typename Float>
template<typename RandomAccessIterator1,
         typename RandomAccessIterator2>
void Subspan<dim,Float>::find_affine_coefficients( RandomAccessIterator1 p,
                                                   RandomAccessIterator2 lambdas )
{
    // compute relative positions of p i.e. u = p - origin
    for( size_type i=0; i<dim; ++i ) {
        u[i] = p[i] - AFFINE_ORIGIN[i];
    }

    // calculate Q^T * u into w
    for( size_type i=0; i<dim; ++i ) {
        w[i] = 0.0;
        for( size_type k=0; k<dim; ++k ) {
            w[i] += Q[i][k] * u[k];
        }
    }

    // We compute the coefficients by backsubstitution. Notice that
    //
    //   c = \sum_{ i\in M } \lambda_i * (point_list[i] - origin[i])
    //     = \sum_{ i\in M } \lambda_i * point_list[i] + (1-s) origin[i]
    // where
    //
    //   s = \sum_{ i\in M } \lambda_i
    //
    // We compute the coefficient (1-s) of the origin in the value origin_lambda:

    Float origin_lambda = 1.0;
    for( int j = int(r)-1; j>=0; --j ) {
        for( size_type k = j+1; k<r; ++k ) {
            w[j] -= *(lambdas + k) * R[k][j];
        }
        *(lambdas + j) = w[j] / R[j][j];
        origin_lambda -= *(lambdas + j);
    }
    // the r-th coefficient corresponds to the origin (cf. remove_point()):
    *(lambdas+r) = origin_lambda;
}


template<int dim, typename Float>
Float Subspan<dim,Float>::time_to_enter_ball( size_type p, size_type q ) const
// Compute the follwing formula:
//
//       ||p||^2 - ||q||^2 + 2< x, q-p >
// t_p = -------------------------------
//                  2< d, q-p >
{
    using std::inner_product;

    Float p_squared = point_list[p].squared_norm(); // inner_product( p, p+dim, p, Float(0.0) );
    Float q_squared = point_list[q].squared_norm(); // inner_product( q, q+dim, q, Float(0.0) );


    Float q_minus_p[dim];
#pragma unroll
    for( size_type k=0; k<dim; ++k ) {
        q_minus_p[k] = point_list[q][k] - point_list[p][k];
    }

    Float x_times_qmp = inner_product( x, x+dim, q_minus_p, Float(0.0) );
    Float d_times_qmp = inner_product( d, d+dim, q_minus_p, Float(0.0) );

    Float res = (p_squared - q_squared + 2*x_times_qmp)/(2*d_times_qmp);

    return res;
}

template<int dim, typename Float>
/**
 * @brief time_to_enter_ball_opt
 * Compute the following formula:
 *
 *       ||p_i||^2 - ||p_j||^2 + 2< x, p_j-p_i >
 * t_p = ---------------------------------------
 *                 2< d, p_j-p_i >
 * Uses simple loop unrolling for efficiency reasons
 * @param i
 * @param j
 * @return
 */
Float Subspan<dim,Float>::time_to_enter_ball_opt( size_type i, size_type j ) const
{
    Float p_squared = point_list[i].squared_norm(); // inner_product( p, p+dim, p, Float(0.0) );
    Float q_squared = point_list[j].squared_norm(); // inner_product( q, q+dim, q, Float(0.0) );
    Float q_minus_p[dim];
    Float x_times_qmp = 0;
    Float d_times_qmp = 0;
    size_type k=0;
    for( k=0; k<dim-1; ++(++k) ) {
        q_minus_p[k]   = point_list[j][k]   - point_list[i][k];
        q_minus_p[k+1] = point_list[j][k+1] - point_list[i][k+1];
        x_times_qmp += x[k]*q_minus_p[k] + x[k+1]*q_minus_p[k+1];
        d_times_qmp += d[k]*q_minus_p[k] + d[k+1]*q_minus_p[k+1];
    }
    if( k==(dim-1) ) {
        q_minus_p[k] = point_list[j][k] - point_list[i][k];
        x_times_qmp += x[k]*q_minus_p[k];
        d_times_qmp += d[k]*q_minus_p[k];
    }

    Float res = (p_squared - q_squared + 2.0*x_times_qmp)/(2.0*d_times_qmp);

    return res;
}

template<int dim, typename Float>
/**
 * @brief time_to_enter_ball_opt
 * Compute the following formula:
 *
 *       ||p_i||^2 - ||p_j||^2 + 2< x, p_j-p_i >
 * t_p = ---------------------------------------
 *                 2< d, p_i-p_j >
 * Uses asm code and sse for efficiency reasons
 * @param i
 * @param j
 * @return
 */
Float Subspan<dim,Float>::time_to_enter_ball_sse( size_type p, size_type q ) const
{
    long dimension = dim;
    Float p_squared = point_list[p].squared_norm(); // inner_product( p, p+dim, p, Float(0.0) );
    Float q_squared = point_list[q].squared_norm(); // inner_product( q, q+dim, q, Float(0.0) );
    Float x_times_qmp;
    Float d_times_qmp;
    Float* pd = const_cast<Float*>(point_list[p].data());
    Float* qd = const_cast<Float*>(point_list[q].data());

    __asm__ __volatile__(
                "movq   %2, %%r8                  \n\t" // store address of point_list[p] in %r8
                "movq   %3, %%r9                  \n\t" // store address of point_list[q] in %r9
                "movq   %4, %%r10                 \n\t" // store address of x in %r10
                "movq   %5, %%r11                 \n\t" // store address of d in %r11
                "movq   %6, %%rbx                 \n\t"
                "decq   %%rbx                     \n\t" // save dim-1 in %rbx
                "xorq   %%rcx, %%rcx              \n\t" // set counter %rcx to 0
                "xorpd  %%xmm2, %%xmm2            \n\t" // set the inner product <x,q-p> to 0
                "xorpd  %%xmm3, %%xmm3            \n\t" // set the inner product <d,q-p> to 0
                "1:                               \n\t"
                "movupd (%%r9,%%rcx,8), %%xmm0    \n\t" // compute point_list[q]-point_list[p] (p-q) into xmm0
                "subpd  (%%r8,%%rcx,8), %%xmm0    \n\t" // by processing 2 coords in parallel
                "movupd %%xmm0, %%xmm1            \n\t"
                "mulpd  (%%r10,%%rcx,8), %%xmm1   \n\t" // compute 2 components of <x,q-p>
                "mulpd  (%%r11,%%rcx,8), %%xmm0   \n\t" // compute 2 components of <d,q-p>
                "addpd  %%xmm0, %%xmm2            \n\t" // x_times_qmp += x[k] * q_minus_p[k]
                "addpd  %%xmm1, %%xmm3            \n\t" // d_times_qmp += d[k] * q_minus_p[k]
                "incq   %%rcx                     \n\t" // increase counter by 1
                "incq   %%rcx                     \n\t" // increase counter by 1
                "cmpq   %%rbx,%%rcx               \n\t" // loop if there are two more components
                "jl     1b                        \n\t" // to process
                "cmpq   %%rbx, %%rcx              \n\t"
                "jne    2f                        \n\t" // jump if dimension is even
                "movq   (%%r9,%%rcx,8), %%xmm0    \n\t" // process in the odd case the remaining components
                "subsd  (%%r8,%%rcx,8), %%xmm0    \n\t" // compute point_list[q]-point_list[p] (p-q) into xmm0
                "movq   %%xmm0, %%xmm1            \n\t" //
                "mulsd  (%%r10,%%rcx,8), %%xmm1   \n\t" //
                "mulsd  (%%r11,%%rcx,8), %%xmm0   \n\t" //
                "addsd  %%xmm0, %%xmm2            \n\t" //
                "addsd  %%xmm1, %%xmm3            \n\t" //
                "2:                               \n\t" // merge the computed results
                "movhlps %%xmm2, %%xmm0           \n\t" // add the lower and the higher
                "addsd   %%xmm0, %%xmm2           \n\t" // double in %xmm2
                "movhlps %%xmm3, %%xmm1           \n\t" // like above for
                "addsd   %%xmm1, %%xmm3           \n\t" // %xmm3
                "movq  %%xmm2, %1                 \n\t" // store the results in the output
                "movq  %%xmm3, %0                 \n\t" // variables
                : /* outputs */ "=m"(x_times_qmp),
                                "=m"(d_times_qmp)
                : /* inputs */  "m"(pd), //(point_list[p].data()),
                                "m"(qd), //(point_list[q].data()),
                                "m"(x),
                                "m"(d),
                                "m"(dimension)
                : /* clobbered */ "xmm0", "xmm1", "xmm2", "xmm3",
                "r8", "r9", "r10", "r11", "rcx", "rbx"  );

    Float res = (p_squared - q_squared + 2*x_times_qmp)/(2*d_times_qmp);

    return res;
}


template<int dim, typename Float>
template<typename RandomAccessIterator>
const Float Subspan<dim,Float>::distance_to_center( RandomAccessIterator p ) const
{
    using std::sqrt;
    using std::inner_product;

    for( int i=0; i<dim; ++i )
        u[i] = x[i] - *(p+i);

    return sqrt( inner_product( u, u+dim, u, Float(0.0) ) );
}

template<int dim, typename Float>
void Subspan<dim,Float>::append_column() {
    ASSERT( r<dim );

    // compute new column R[r] = Q^T * u
    for( size_type i=0; i<dim; ++i ) {
        R[r][i] = 0.0;
        for( size_type k=0; k<dim; ++k ) {
            R[r][i] += Q[i][k] * u[k];
        }
    }

    // zero all entries R[r][dim-1] down to R[r][r+1]
    for( size_type j=dim-1; j>r; --j ) {
        // j is the index of the entry to be cleared with the help of entry j-1

        // compute the Givens-coefficients c,s
        Float c, s;
        givens( c, s, R[r][j-1], R[r][j] );

        // rotate one R-entry ( the other one is an implizit zero )
        R[r][j-1] = c*R[r][j-1] + s*R[r][j];
        R[r][j]   = 0.0;//c*R[r][j]   - s*R[r][j-1];


        // rotate two Q-columns
        for( size_type i=0; i<dim; ++i ) {
            const Float a = Q[j-1][i];
            const Float b = Q[j][i];
            Q[j-1][i] = c*a + s*b;
            Q[j][i]   = c*b - s*a;
        }
    }
}

template<int dim, typename Float>
void Subspan<dim,Float>::print_matrices( int i ) const
{
    std::cout << "members: [ ";
    for(auto i=0; i<r+1; ++i) std::cout << members[i] << " ";
    std::cout << "]\n";
    switch(i) {
    case 0: {
        Float B[r][r];
        std::cout << "Q^TQ:\n";
        for( int i=0; i<r; ++i ) {
            for( int j=0; j<r; ++j ) {
                B[i][j] = 0.0;
                for( int k=0; k<dim; ++k ) {
                    B[i][j] += Q[i][k]*Q[j][k];
                }
                std::cout << B[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Q:\n";
        for( int i=0; i<dim; ++i ) {
            for( int j=0; j<dim; ++j ) std::cout << Q[j][i] << " ";
            std::cout << "\n";
        }
        std::cout << "R:\n";
        for( int i=0; i<dim; ++i ) {
            for( int j=0; j<r; ++j ) std::cout << R[j][i] << " ";
            std::cout << "\n";
        }
        std::cout << "A:\n";
        Float A[dim][dim];
        for( int i=0; i<r; ++i ) {
            for( int j=0; j<dim; ++j ) {
                A[i][j] = 0.0;
                for( int k=0; k<dim; ++k ) {
                    A[i][j] += Q[k][j] * R[i][k];
                }
                std::cout << A[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "A':\n";
        for( int i=0; i<r; ++i ) {
            for( int j=0; j<dim; ++j ) {
                A[i][j] = point_list[members[i]][j]-point_list[members[r]][j];
                std::cout << A[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
    case 1:
        std::cout << "d: [ ";
        for( int i=0; i<dim; ++i ) std::cout << d[i] << " ";
        std::cout << "]\n";
    case 2:
        std::cout << "x: [ ";
        for( int i=0; i<dim; ++i ) std::cout << x[i] << " ";
        std::cout << "]\n";
    }
}


template<int dim,typename Float>
void Subspan<dim,Float>::hessenberg_clear( size_type pos )
{

    // clear new subdiagonal entries
    for( ; pos<r; ++pos ) {
        // pos is the index of the entry to be cleared

        //compute Givens coefficients c,s
        Float c, s;
        givens( c, s, R[pos][pos], R[pos][pos+1] );

        // rotate partial R-rows (of the first pair, only one entry is needed,
        // the other one is an implicit zero )
        R[pos][pos]   = c*R[pos][pos]   + s*R[pos][pos+1];
        R[pos][pos+1] = 0.0;//c*R[pos][pos+1] - s*R[pos][pos];

        // then begin at column pos+1
        for( size_type j = pos+1; j<r; ++j ) {
            const Float a = R[j][pos];
            const Float b = R[j][pos+1];
            R[j][pos]   = c*a + s*b;
            R[j][pos+1] = c*b - s*a;
        }

        // rotate Q-colums
        for( size_type i=0; i<dim; ++i ) {
            const Float a = Q[pos][i];
            const Float b = Q[pos+1][i];
            Q[pos][i]   = c*a + s*b;
            Q[pos+1][i] = c*b - s*a;
        }
    }
}


template<int dim,typename Float>
void Subspan<dim,Float>::special_update()
// update current QR-decomposition "A=QR" to A + u*[1,...,1] = Q'R'
{

    // compute w=Q^T * u
    for( size_type i=0; i<dim; ++i ) {
        w[i] = 0.0;
        for( size_type k=0; k<dim; ++k ) {
            w[i] += Q[i][k] * u[k];
        }
    }

    // rotate w down to a multiple of the first unit-vector; record these
    // ops in R and Q
    for( int k = dim-1; k>0; --k ) {
        // k is the index of the entry to be cleared with the help of entry k-1

        // compute givens coefficients
        Float c, s;
        givens( c, s, w[k-1], w[k] );

        // rotate w-entry
        w[k-1] = c*w[k-1] + s*w[k];

        // rotate two R-rows
        // the first column has to be treated seperately in order to
        // account fo the implicit zero in R[k-1][k]
        R[k-1][k]    = -s * R[k-1][k-1];
        R[k-1][k-1] *=  c;
        for( size_type j=k; j<r; ++j ) {
            const Float a = R[j][k-1];
            const Float b = R[j][k];
            R[j][k-1] = c*a + s*b;
            R[j][k]   = c*b - s*a;
        }

        // rotate two Q-colums
        for( size_type i=0; i<dim; ++i ) {
            const Float a = Q[k-1][i];
            const Float b = Q[k][i];
            Q[k-1][i] = c*a + s*b;
            Q[k][i]   = c*b - s*a;
        }
    }

    // add w*(1,...,1)^T to new R
    // which means simply to add u[0] to each column
    // since the other entries of u have just been eliminated
    for( size_type j=0; j<r; ++j ) {
        R[j][0] += w[0];
    }
    // clear subdiagonal entries
    hessenberg_clear(0);
}

template<int dim,typename Float>
void Subspan<dim,Float>::compute_qr_complete()
{

    DEBUG_OUTPUT( << "Computing the QR-decomposition from scratch...\n" );

    // Initialize Q and R
    for( int i=0; i<dim; ++i )
        for( int j=0; j<dim; ++j ) {
            Q[i][j]  = i == j ? 1.0 : 0.0;
        }

    // Insert the base-vectors in R
    for( int i=0; i<r; ++i ) {
        for( int j=0; j<dim; ++j ) {
            R[i][j] = point_list[members[i]][j] - AFFINE_ORIGIN[j];
        }
    }
    // Compute the QR-decomposition
    for( int pos=0; pos<r; ++pos ) {

        // zero all entries R[pos][dim-1] down to R[pos][r+1]
        for( int j=dim-1; j>pos; --j ) {
            // j is the index of the entry to be cleared with the help of entry j-1

            // compute the Givens-coefficients c,s
            Float c, s;
            givens( c, s, R[pos][j-1], R[pos][j] );

            // rotate one R-entry ( the other one is an implizit zero )
            R[pos][j-1] = c*R[pos][j-1] + s*R[pos][j];
            R[pos][j]   = Float(0.0);//c*R[pos][j]   - s*R[pos][j-1];

            for( int k=pos+1; k<r; ++k ) {
                const Float a = R[k][j-1];
                const Float b = R[k][j];
                R[k][j-1] = c*a + s*b;
                R[k][j]   = c*b - s*a;
            }

            // rotate two Q-columns
            for( int i=0; i<dim; ++i ) {
                const Float a = Q[j-1][i];
                const Float b = Q[j][i];
                Q[j-1][i] = c*a + s*b;
                Q[j][i]   = c*b - s*a;
            }
        }
    }
}

#include <iostream>
#include <iomanip>

template<int dim, typename Float>
template<class Archive>
void Subspan<dim,Float>::serialize( Archive & ar, const unsigned int )
// Store the data in aspecified Archive. This routine is needed to transfer
// The data of an object over the network as needed by Boost.MPI.
{

    ar & BOOST_SERIALIZATION_NVP( membership );

    ar & BOOST_SERIALIZATION_NVP( members );

    for( int i=0; i<dim; ++i ) {
        for( int j=0; j<dim; ++j ) {
            ar & BOOST_SERIALIZATION_NVP( Q[i][j] );
        }
    }
    for( int i=0; i<dim; ++i ) {
        for( int j=0; j<dim; ++j ) {
            ar & BOOST_SERIALIZATION_NVP( R[i][j] );
        }
    }
    for( int i=0; i<dim; ++i ) {
        ar & BOOST_SERIALIZATION_NVP( x[i] );
    }
    for( int i=0; i<dim; ++i ) {
        ar & BOOST_SERIALIZATION_NVP( d[i] );
    }
    ar & BOOST_SERIALIZATION_NVP( r );
}

template<int dim,typename Float>
auto Subspan<dim, Float>::infinite_hash_value() const -> hash_value_type {
    return Marks::compute_exceptional(r,point_list.size());
}


template<int dim, typename Float>
auto Subspan<dim,Float>::hash_value() const -> hash_value_type {
    return Marks::compute(members,point_list.size());
}

} // NAMESPACE

#endif
