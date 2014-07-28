#ifndef DELAUNAYDECOMPOSITION_CPP
#define DELAUNAYDECOMPOSITION_CPP

#include "configure.h"
#include "delaunaydecomposition.h"

#include <future>

namespace NAMESPACE {


template<int dim, typename Float>
void DelaunayDecomposition<dim,Float>::thread_function(std::promise<ball_list_type> &ball_list,
                     ball_type *ball)
{
    ball_list.set_value(explore(ball));
}

template<int dim, typename Float>
void DelaunayDecomposition<dim,Float>::find_all_simplices(ball_type *start_ball)
{
    agenda_type agenda;
    agenda.push_front(start_ball);
    marks[start_ball->hash_value()] += 1;
    simplex_list.push_back(start_ball->simplex());
    while(!agenda.empty()) {
#ifdef USE_THREADS
        //auto ball_list = explore(ball);
        u_long num_threads = 32;
        std::vector<std::future<ball_list_type>> ball_list_futures;
        std::vector<ball_type*> ball_list;
        for(auto i=0;i<num_threads && !agenda.empty();++i) {
            auto ball = agenda.front();
            agenda.pop_front();
            ball_list_futures.push_back(std::async(&DD::DelaunayDecomposition<dim,Float>::explore,this, ball));
        }

        for(auto i=0;i<ball_list_futures.size();++i) {
            for(auto ball : ball_list_futures[i].get()) {
                if(!explored(ball)) {
                    agenda.push_front(ball);
                    simplex_list.push_back(ball->simplex());
                } else {
                    delete ball;
                }
            }
            //delete ball; // the incidence information is in simplex_list
        }
#else
        auto ball = agenda.front();
        agenda.pop_front();

        auto ball_list = explore(ball);

        for(auto ball : ball_list) {
            if(!explored(ball)) {
                agenda.push_front(ball);
                simplex_list.push_back(ball->simplex());
            } else {
                delete ball;
            }
        }
        delete ball; // the incidence information is in simplex_list
#endif
    }
}


template<int dim,typename Float>
auto DelaunayDecomposition<dim,Float>::explore(ball_type *ball) -> ball_list_type
{
    //std::cout << "Exploring " << ball->simplex() << std::endl;
    ball_list_type ball_list;
    if(ball->size() == 1) return ball_list;

    Float lambdas[ball->size()];

    ball->update_miniball_center();
    ball->compute_center_coefficients(lambdas);

    ON_DEBUG( std::cout << "Lambdas: [ "; for(auto i=0;i<ball->size();++i) std::cout << lambdas[i] <<" "; std::cout << "]\n" );

    for(auto i=0;i<ball->size();++i) {
        int direction      = pos_neg(lambdas[i]);
        auto removed_index = ball->global_index(i);
        mark_type hash = Marks::compute_without_index(ball->get_members(),
                                                      removed_index,
                                                      point_list.size());

        // Every (d-1)-face of a Delaunay d-simplex is common two exactly to Delaunay d-simplices
        // once we cross it, we don't ever need to cross it again
#ifdef USE_THREADS
        mutex.lock();
#endif
        ++marks[hash];
        int num = marks.at(hash);
#ifdef USE_THREADS
        mutex.unlock();
#endif
        if(num < 2) {

        auto clone         = new ball_type(*ball);
        clone->remove_point(i);

        ON_DEBUG( std::cout << "Exploring further: " << ball->simplex() << std::endl);
        ON_DEBUG( std::cout << "Removing point: " << removed_index << std::endl);
        ON_DEBUG( clone->print_center() );
        clone->update_driver();

        long  next_point_index = -1;
        Float next_point_time  = Float(0.0);
#ifdef USE_OPENCL
        next_point_index = find_next_point_opencl( clone,
                                                   removed_index,
                                                   direction,
                                                   next_point_time );
#else
        next_point_index = find_next_point( clone,
                                            removed_index,
                                            direction,
                                            next_point_time );
#endif
        if(next_point_index > -1) {
            // At this stead, clone represents a face of ball that
            // is common to a neighboring Delaunay simplex
            // We may enter the new Delaunay simplex by adding
            // next_point_index
            ON_DEBUG(std::cout << "|| Found new vertex: " << next_point_index << std::endl);
            clone->add_point_plain( next_point_index );
            clone->update_miniball_center();
            ON_DEBUG( clone->print_center() );
            clone->update_driver();
            ON_DEBUG( std::cout << "==========> Found simplex: " << clone->simplex() << std::endl);
            ball_list.push_back(clone);
        } else { // next_point_index == 1
            // Here, clone represents a face of "ball" that is a face
            // of the convex hull of point_list. Thus, the neighboring
            // simplex contains the point at infinity. We don't
            // explore it further.
            ON_DEBUG(std::cout << "No adaequate vertex!\n");
            delete clone;
        }
    }
    }
    ON_DEBUG( std::cout << "Next balls to explore:\n";
            for(auto ball : ball_list) {
              std::cout << "[ ";
              for(auto vertex : ball->simplex())
                std::cout << vertex << " ";
              std::cout << "]\n";});
    return ball_list;
}





template<int dim, typename Float>
template<int number>
int DelaunayDecomposition<dim,Float>::find_next_point_helper( ball_type           *ball,
                                                              const unsigned long  last_index,
                                                              Float               &min_time )
{
    switch(number) {
    case  1: min_time = -std::numeric_limits<Float>::max(); break;
    case -1: min_time =  std::numeric_limits<Float>::max(); break;
    case  0: min_time =  std::numeric_limits<Float>::max(); break;
    }
    int min = -1;
    for( u_long i=0; i < point_list.size(); ++i ) {
        if( !ball->is_member( i ) && (i != last_index) ) {
            Float t_p = time_to_enter_ball( ball, i );
            switch(number) {
            case  1:
                if( t_p < 0.0 && t_p > min_time ) {
                    min_time = t_p;
                    min = i;
                }
                break;
            case  -1:
                if( t_p > 0.0 && t_p < min_time ) {
                    min_time = t_p;
                    min = i;
                }
                break;
            case  0:
                if( std::abs(t_p) < std::abs(min_time) ) {
                    min_time = -t_p; // correct direction of the driver
                    min = i;
                }
            }
        }
    }
    return min;
}

template<int dim, typename Float>
/**
 * @brief DelaunayDecomposition<dim, Float>::find_next_point
 * Helper routine to find the next point in S that enters the ball if we let
 * it shrink or grow. The shrinking/growing is determined by the parameter
 * direction.
 * The reference paramter min_time is used to return the time at which the
 * returned point enters the ball.
 * @param ball
 * @param last_index
 * @param direction
 * @param min_time
 * @return
 */
int DelaunayDecomposition<dim,Float>::find_next_point( ball_type           *ball,
                                                       const unsigned long  last_index,
                                                       const int            direction,
                                                       Float               &min_time )
{
    // The driver vector points from perpendicular of the old balls center on the affine hull to the
    // center we now have to walk in the right direction
    switch( direction ) {
    case  1:
        // the driver points inside the simplex. we have to move in the direction
        // of the negative driver to find the next simplex
#ifdef USE_OPT
        return find_next_point_helper_opt<1>( ball, last_index, min_time);
#else
        return find_next_point_helper<1>( ball, last_index, min_time );
#endif
    case -1:
        // the driver points outside the simplex. we have to move into the direction
        // of the driver to find the next simplex
#ifdef USE_OPT
        return find_next_point_helper_opt<-1>( ball, last_index, min_time);
#else
        return find_next_point_helper<-1>( ball, last_index, min_time );
#endif
    case  0:
        // this case is for finding the initial simplex.
#ifdef USE_OPT
        return find_next_point_helper_opt<0>( ball, last_index, min_time);
#else
        return find_next_point_helper<0>( ball, last_index, min_time );
#endif
    default:
        std::cerr << "Wrong direction in function 'find_next_point'\n";
        exit(1);
    }
}

/******************************************************************************
 * Optimized Routines
 ******************************************************************************/



/**
 * @brief find_next_point_helper_opt
 * Helper routine for the heuristical approach to the search of the next point
 * We choose a initial point that satisfies our needs for the direction in which
 * we will walk and then select all the points inside the ball determined by this
 * first point. Among those candidates the next point that enters the ball must
 * be found.
 * @param ball
 * @param last_index
 * @param min_time
 * @return
 */
template<int dim, typename Float>
template<int number>
long DelaunayDecomposition<dim,Float>::find_next_point_helper_opt(ball_type* ball,
                                                                  const unsigned long last_index,
                                                                  Float &min_time )
{
    int min = -1;
    switch(number) {
    case  1: min_time = -std::numeric_limits<Float>::max(); break;
    case -1: min_time = std::numeric_limits<Float>::max(); break;
    case  0: min_time = std::numeric_limits<Float>::max(); break;
    default:
        std::cerr << "Wrong template parameter in function 'find_next_point_opt_helper': "
                  << number << '\n';
        exit(1);
    }

    for( auto i=0; i < point_list.size(); ++i ) {
        if( ball->is_member( i ) || (i == last_index) ) continue;

        Float t_p = time_to_enter_ball( ball, i );
        // Optimization: We choose a first point "min" that satisfies our needs
        // After that, we only examine the points of the point_list that lie
        // within the circle defined by ball and "min"
        switch(number) {
        case  1:
            if( t_p < 0.0 ) {
                min_time = t_p;
                min = i;
                goto outside;
            }
            break;
        case -1:
            if( t_p > 0.0 ) {
                min_time = t_p;
                min = i;
                goto outside;
            }
            break;
        case  0:
            min_time = t_p;
            min = i;
            goto outside;
            break;
        }
    }

outside:

    // We are heading for infinity if no valid index could be determined
    if( min == -1 ) return min;

    // Compute the center of the ball defined by the points in "ball" and
    // point_list[index] and get its radius
    Float center[dim];
    ball->get_center( center );
    Float actual_radius_squared = squared_distance( center, min );

    // Find all points that lie inside the ball
    std::vector<u_long> points_inside_ball;

    for( auto i=min; i<point_list.size(); ++i ) {
        if( ball->is_member( i ) || (i == last_index) ) continue;
        Float dist_squared = squared_distance( center, i );
        if( dist_squared < actual_radius_squared )
            points_inside_ball.push_back(i);

    }

    // If no other points of point_list lie inside the circle, we found the
    // next point to enter the ball: it is point_list[index].
    if( points_inside_ball.empty() ) {
        return min;
    }

    // Now we search the next point to enter the ball among all points inside the
    // ball
    for(u_long i=0; i < points_inside_ball.size(); ++i ) {

        Float t = time_to_enter_ball( ball, points_inside_ball[i] );

        // Set the point that firstly enters the ball according to the direction
        // in that we walk (template parameter)
        switch(number) {
        case  1:
            if( t > min_time && t < Float(0.0) ) {
                min_time = t;
                min = points_inside_ball[i];
            }
            break;
        case -1:
            if( t < min_time  && t > Float(0.0) ) {
                min_time = t;
                min = points_inside_ball[i];
            }
            break;
        case  0:
            if( std::abs(t) < std::abs(min_time) ) {
                min_time = t;
                min = points_inside_ball[i];
            }
            break;
        }
    }
    return min;
}

/**
 * @brief find_next_point
 * Helper routine to find the next point in S that enters the ball if we let
 * it shrink or grow. The shrinking/growing is determined by the parameter
 * direction.
 * The reference paramter min_time is used to return the time at which the
 * returned point enters the ball.
 * @param ball
 * @param last_index
 * @param direction
 * @param min_time
 * @return
 */
template<int dim, typename Float>
long DelaunayDecomposition<dim,Float>::find_next_point_opt(ball_type           *ball,
                                                           const unsigned long  last_index,
                                                           const int            direction,
                                                           Float               &min_time )
{
    // The direction vector points from the driver to the center of the old ball
    // we now have to walk in the right direction
    switch(direction) {
    case  1: return find_next_point_helper_opt<1>( ball, last_index, min_time );
    case -1: return find_next_point_helper_opt<-1>( ball, last_index, min_time );
    case  0: return find_next_point_helper_opt<0>( ball, last_index, min_time );
    default:
        std::cerr << "Wrong direction given in function 'find_next_point_opt'\n";
        exit(1);
    }
}

template<int dim,typename Float>
Float DelaunayDecomposition<dim,Float>::squared_distance( const Float         *center,
                                                          const unsigned long  index ) const
{
#ifdef USE_SSE
        long dimension = dim;
        Float res = 0.0;
        __asm__ __volatile__ ( "movq    %1, %%r8                \n\t" // center   -> r8
                               "movq    %2, %%r9                \n\t" // point_list[index] -> r9
                               "movq    %3, %%rbx               \n\t" // dim      -> rbx
                               "xorq    %%rcx, %%rcx            \n\t" // 0        -> rcx
                               "decq    %%rbx                   \n\t" // rbx = dim-1
                               "xorpd   %%xmm1, %%xmm1          \n\t" // 0        -> xmm1
                               "1:                              \n\t" //
                               "movupd  (%%r8,%%rcx,8), %%xmm0  \n\t" // center[i] - point_list[index][i]
                               "subpd   (%%r9,%%rcx,8), %%xmm0  \n\t" // into xmm0
                               "mulpd   %%xmm0, %%xmm0          \n\t" // xmm0 = xmm0 * xmm0
                               "addpd   %%xmm0, %%xmm1          \n\t" // xmm1 accumulates squares
                               "incq    %%rcx                   \n\t" //
                               "incq    %%rcx                   \n\t" // i+=2
                               "cmpq    %%rbx, %%rcx            \n\t" // for-loop
                               "jl      1b                      \n\t" //
                               "cmpq    %%rbx, %%rcx            \n\t" // if the dimension is odd
                               "jne     2f                      \n\t" // compute the last coordinate
                               "movq    (%%r8,%%rcx,8), %%xmm0  \n\t"
                               "subsd   (%%r9,%%rcx,8), %%xmm0  \n\t"
                               "mulsd   %%xmm0, %%xmm0          \n\t"
                               "addsd   %%xmm0, %%xmm1          \n\t"
                               "2:                              \n\t"
                               "movhlps %%xmm1, %%xmm0          \n\t" // accumulate the high and low
                               "addsd   %%xmm0, %%xmm1          \n\t" // words of xmm1
                               "movq    %%xmm1, %0              \n\t" // return result in res
                               : /* out */ "=m"(res)
                               : /* in  */ "rm"(center),
                               "rm"(point_list[index].pointer()),
                               "m" (dimension)
                               : /* brk */ "xmm0", "xmm1", "rbx", "rcx", "r8", "r9"
                               );
        return res;
#else
        Float res = 0.0;
        for( auto i=0; i<dim-1; i+=2 ) {
            res += (center[i]   - point_list[index][i])*(center[i] - point_list[index][i]);
            res += (center[i+1] - point_list[index][i+1])*(center[i+1] - point_list[index][i+1]);
        }
        if( dim%2 == 1 ) {
            res += (center[dim-1] - point_list[index][dim-1])*(center[dim-1] - point_list[index][dim-1]);
        }
        return res;

#endif
}

} // NAMESPACE

#endif // DELAUNAYDECOMPOSITION_CPP
