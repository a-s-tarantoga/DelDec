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
#ifndef DELAUNAYDECOMPOSITION_CPP
#define DELAUNAYDECOMPOSITION_CPP

#include "configure.h"
#include "delaunaydecomposition.h"

#include <future>

namespace NAMESPACE {

template<int dim,typename Float>
auto DelaunayDecomposition<dim,Float>::find_initial_simplex(uint start_index) -> ball_type*{

    ASSERT(start_index < point_list.size());

    // find two points of maximal distance
    Float max_distance = Float(0.0);
    long  max = -1;
    for(auto j=0; j<point_list.size(); ++j) {
        Float dist_squared = Float(0.0);
        for(auto i=0; i<dim; ++i) {
            dist_squared += sqr(point_list[start_index][i] - point_list[j][i]);
        }
        if(max_distance < dist_squared) {
            max = j;
            max_distance = dist_squared;
        }
    }

    // set up the starting point
    point_type start_center;
    srand(clock());
    for(auto i=0; i<dim; ++i) {
        start_center[i] = static_cast<Float>((point_list[start_index][i] + point_list[max][i])/2.0);
    }
    start_center.perturb(0.001*sqrt(max_distance));
    ON_DEBUG( std::cout << "Initial center:    " << start_center << std::endl);

    // compute the driver for start_point
    Float min_distance = std::numeric_limits<Float>::max();
    long closest_point = -1;
    for(auto j=0;j<point_list.size();++j) {
        Float dist_squared = Float(0.0);
        for(auto i=0;i<dim;++i) {
            dist_squared += sqr(start_center[i]-point_list[j][i]);
        }
        if(min_distance > dist_squared) {
            closest_point = j;
            min_distance = dist_squared;
        }
    }
    ON_DEBUG( std::cout << "Closest point:     " << closest_point << std::endl );

    // Define the first ball
    ball_type *ball = new ball_type( point_list,
                                     start_center,
                                     closest_point);
    ball->update_driver();

    ON_DEBUG( std::cout << "Simplex:           " << ball->simplex() << std::endl );
    ON_DEBUG( std::cout << "Center:            "; ball->print_center() );
    ON_DEBUG( std::cout << "Driver:            "; ball->print_driver() );
    ON_DEBUG( std::cout << "Start inserting points\n" );

    while((ball->size()<point_list.size()) &&
          (ball->size()<=dim)) {
        // find the point that enters the ball
        Float min_time = Float(0.0);
        auto next_point = find_next_point( ball,
                                           -1,
                                           0,
                                           min_time );
        ball->update_center( min_time );
        ball->add_point( next_point );

        ON_DEBUG( std::cout << "New vertex:        " << next_point << std::endl );
        ON_DEBUG( std::cout << "Simplex:           " << ball->simplex() << std::endl );
        ON_DEBUG( std::cout << "Center:            "; ball->print_center() );
        ON_DEBUG( std::cout << "Driver:            "; ball->print_driver() );

        Float dist_center_to_driver = ball->update_driver();
        // Test if the ball cannot be grown any further
        if((dist_center_to_driver < 1e-100) &&
           (ball->size() <= dim) &&
           (ball->size() < point_list.size())) {
            ON_DEBUG(std::cout << "Stuck in a dead end!\n");
            delete ball;
            return NULL;
        }
    }

    // At this stead, the center of the ball is a Voronoi vertex
    // and thus the points defining the ball form a
    // Delaunay simplex
    ball->update_miniball_center();
    ON_DEBUG(std::cout << "==> Initial simplex found: " << ball->simplex() << std::endl);
    return ball;
}

template<int dim, typename Float>
void DelaunayDecomposition<dim,Float>::find_all_simplices(ball_type *start_ball)
{
    agenda_type agenda;
    agenda.push_front(start_ball);
    marks[start_ball->hash_value()] += 1;
    simplex_list.push_back(start_ball->simplex());

    while(!agenda.empty()) {

#ifndef USE_THREADS
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
#else
        // Asynchronus compute explore the simplices
        u_long num_threads = 2;
        std::vector<std::future<ball_list_type>> ball_list_futures;
        std::vector<ball_type*> ball_list;
        for(auto i=0;i<num_threads && !agenda.empty();++i) {
            auto ball = agenda.front();
            agenda.pop_front();
            ball_list.push_back(ball);
            ball_list_futures.push_back(std::async(&DD::DelaunayDecomposition<dim,Float>::explore,this, ball));
        }
        // Join the asynchronous threads and get their results
        for(auto i=0;i<ball_list_futures.size();++i) {
            for(auto ball : ball_list_futures[i].get()) {
                if(!explored(ball)) {
                    agenda.push_front(ball);
                    simplex_list.push_back(ball->simplex());
                } else {
                    delete ball;
                }
            }
        }
        // Don't forget to delete the explored balls
        for(auto b : ball_list) {
            delete b;
        }
#endif
    }
}


template<int dim,typename Float>
auto DelaunayDecomposition<dim,Float>::explore(ball_type *ball) -> ball_list_type
{
    ball_list_type ball_list;
    if(ball->size() == 1) return ball_list;

    Float lambdas[ball->size()];

    ball->update_miniball_center();
    ball->compute_center_coefficients(lambdas);

    ON_DEBUG( std::cout << "Center lambdas:    [ "; for(auto i=0;i<ball->size();++i) std::cout << lambdas[i] <<" "; std::cout << "]\n" );

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
#ifdef USE_OPENCL
            ball->get_center(io_center);
            cl_queue->enqueueWriteBuffer(cl_device_center, CL_TRUE, 0, dim*sizeof(Float), io_center);
#endif
            auto clone         = new ball_type(*ball);
            clone->remove_point(i);
            clone->update_driver();

            ON_DEBUG( std::cout << "Exploring further: " << ball->simplex() << std::endl);
            ON_DEBUG( std::cout << "Removing point:    " << removed_index << std::endl);
            ON_DEBUG( std::cout << "Center:            "; ball->print_center() );
            ON_DEBUG( std::cout << "Old driver:        "; ball->print_driver() );
            ON_DEBUG( std::cout << "New driver:        "; clone->print_driver() );

            long  next_point_index = -1;
            Float next_point_time  = Float(0.0);
#ifdef USE_OPENCL
            next_point_index = find_next_point_opencl( clone,
                                                       removed_index,
                                                       direction,
                                                       next_point_time );
            //std::cout << "OpenCL: " << next_point_index << " - " << next_point_time << std::endl;
#else
            next_point_index = find_next_point( clone,
                                                removed_index,
                                                direction,
                                                next_point_time );
            //std::cout << "SSE:    " << next_point_index << " - " << next_point_time << std::endl;
#endif
            if(next_point_index > -1) {
                // At this stead, clone represents a face of ball that
                // is common to a neighboring Delaunay simplex
                // We may enter the new Delaunay simplex by adding
                // next_point_index
                clone->add_point_plain( next_point_index );
                clone->update_driver();

                ON_DEBUG( std::cout << "==> Found simplex: " << clone->simplex() << std::endl);

                ball_list.push_back(clone);
            } else { // next_point_index == 1
                // Here, clone represents a face of "ball" that is a face
                // of the convex hull of point_list. Thus, the neighboring
                // simplex contains the point at infinity. We don't
                // explore it further.
                ON_DEBUG(std::cout << "==> No vertex found!\n");
                delete clone;
            }
        }
    }
    ON_DEBUG( std::cout << "Explored Simpl.:   "; for(auto b : ball_list) std::cout << b->simplex() << " "; std::cout << std::endl; );
    return ball_list;
}


template<int dim, typename Float>
bool DelaunayDecomposition<dim,Float>::explored(ball_type* ball)
{
#ifdef USE_THREADS
    mutex.lock();
#endif
    marks[ball->hash_value()] += 1;
    int num = marks.at(ball->hash_value());
#ifdef USE_THREADS
    mutex.unlock();
#endif
    return (num>1);
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

#ifdef USE_OPENCL
template<int dim, typename Float>
void DelaunayDecomposition<dim,Float>::init_OpenCL() {
    typedef std::vector<cl::Platform> platform_list_type;
    typedef std::vector<cl::Device>   device_list_type;

    std::stringstream kernel_source;
    kernel_source <<
           "/*                                                                                          \n"
           " * Compute the following formula:                                                           \n"
           " *                                                                                          \n"
           " *       ||p_i||^2 - ||p_j||^2 + 2< x, p_j-p_i >                                            \n"
           " * t_p = ---------------------------------------                                            \n"
           " *                 2< d, p_j-p_i >                                                          \n"
           " * i     is the global index                                                                \n"
           " * index refers to point p_j                                                                \n"
           " */                                                                                         \n"
           "#ifdef cl_khr_fp64                                                                          \n"
           "    #pragma OPENCL EXTENSION cl_khr_fp64 : enable                                           \n"
           "#elif defined(cl_amd_fp64)                                                                  \n"
           "    #pragma OPENCL EXTENSION cl_amd_fp64 : enable                                           \n"
           "#else                                                                                       \n"
           "    #error \"Double precision floating point not supported by OpenCL implementation.\"      \n"
           "#endif                                                                                      \n"
           "#define Float double                                                                        \n"
           "#define dim " << dim << "                                                                   \n"
           "__kernel void kernel_compute( __constant Float *point_list,                                 \n"
           "                              __constant Float *squared_norm_list,                          \n"
           "                              __global Float *center,                                       \n"
           "                              __global Float *driver,                                       \n"
           "                              const unsigned long index,                                    \n"
           "                              __global Float *result)                                       \n"
           "{                                                                                           \n"
           "    int point_i_base  = get_global_id(0) * dim;             // point to test                \n"
           "    int point_j_base  = index * dim;                        // any point from the ball      \n"
           "    int id            = get_global_id(0);                                                   \n"
           "                                                                                            \n"
           "    Float p_minus_q[dim];                                                                   \n"
           "#pragma unroll                                                                              \n"
           "    for(int i=0;i< dim;++i) {                                                               \n"
           "        p_minus_q[i] = point_list[point_j_base + i] - point_list[point_i_base + i];         \n"
           "    }                                                                                       \n"
           "                                                                                            \n"
           "    Float x_times_pmq = 0.0;                                                                \n"
           "    Float d_times_pmq = 0.0;                                                                \n"
           "#pragma unroll                                                                              \n"
           "    for(int i=0; i<dim; ++i) {                                                              \n"
           "        x_times_pmq += center[i] * p_minus_q[i];                                            \n"
           "        d_times_pmq += driver[i] * p_minus_q[i];                                            \n"
           "    }                                                                                       \n"
           "    x_times_pmq *= 2.0;                                                                     \n"
           "    d_times_pmq *= 2.0;                                                                     \n"
           "                                                                                            \n"
           "    result[id] = (squared_norm_list[id] - squared_norm_list[index] + x_times_pmq)/d_times_pmq;\n"
           "                                                                                            \n"
           "}                                                                                           \n";


    // try to get all supported platforms
    platform_list_type all_platforms;
    cl::Platform::get(&all_platforms);

    if(all_platforms.size()==0) {
        std::cout << "No OpenCL Platforms!" << std::endl;
        exit(1);
    }

    // print all platforms
    cl::Platform default_platform = all_platforms[0];
    for(platform_list_type::iterator it=all_platforms.begin();
        it!=all_platforms.end();++it) {
        std::cout << "Platform: " << it->getInfo<CL_PLATFORM_NAME>() << std::endl;
    }

    // Find all available OpenCL devices
    device_list_type all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices); // We want GPUs

    if(all_devices.size()==0) {
        std::cout << "No OpenCL Devices!" << std::endl;
        exit(1);
    }

    cl::Device default_device = all_devices[all_devices.size()-1];
    for(device_list_type::iterator it=all_devices.begin();it!=all_devices.end();++it) {
        std::cout << "Devices: " << it->getInfo<CL_DEVICE_NAME>() << std::endl;
    }
    std::cout << "Default Device: " << default_device.getInfo<CL_DEVICE_NAME>() << std::endl;

    device_list_type default_device_vector;
    default_device_vector.push_back(default_device);
    cl_context = new cl::Context(default_device_vector);

    // Setup the program on the OpenCL device
    typedef typename std::pair<const char*,unsigned long> src_type;
    cl::Program::Sources sources;
    sources.push_back(src_type(kernel_source.str().c_str(), kernel_source.str().length()));
    cl_program = new cl::Program(*cl_context,sources);
    if(cl_program->build(default_device_vector)!=CL_SUCCESS) {
        std::cout << "Program did not build!" << std::endl;
        std::cout << cl_program->getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    // Prepare for sending data between Host and OpenCL device
    cl_queue = new cl::CommandQueue(*cl_context,default_device);

    // Store the point_list on the OpenCL device
    unsigned long point_list_buffer_size = point_list.size() * dim * sizeof(Float);
    Float point_buffer[point_list.size()*dim];
    cl_point_list_buffer = cl::Buffer(*cl_context, CL_MEM_READ_ONLY, point_list_buffer_size);
    unsigned long counter = 0;
    for(u_long i=0; i<point_list.size(); ++i) {
        for(u_long j=0; j<dim; ++j) {
            point_buffer[counter] = point_list[i][j];
            ++counter;
        }
    }
    cl_queue->enqueueWriteBuffer(cl_point_list_buffer, CL_TRUE, 0, point_list_buffer_size, point_buffer);

    // Store the squared norms for all points on the OpenCL device
    unsigned long squared_norm_buffer_size = point_list.size() * sizeof(Float);
    cl_squared_norm_buffer = cl::Buffer(*cl_context, CL_MEM_READ_ONLY, squared_norm_buffer_size);
    Float norm_buffer[point_list.size()];
    for(u_long i=0; i<point_list.size(); ++i) {
        norm_buffer[i] = point_list[i].squared_norm();
    }
    cl_queue->enqueueWriteBuffer(cl_squared_norm_buffer, CL_TRUE, 0, squared_norm_buffer_size, norm_buffer);

    // Create pinned memory for the center:
    cl_center = cl::Buffer(*cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, dim*sizeof(Float));
    cl_device_center = cl::Buffer(*cl_context, CL_MEM_READ_ONLY, dim*sizeof(Float));
    io_center = (Float*)cl_queue->enqueueMapBuffer(cl_center, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, dim*sizeof(Float));

    // Create pinned memory for the driver:
    cl_driver = cl::Buffer(*cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, dim*sizeof(Float));
    cl_device_driver = cl::Buffer(*cl_context, CL_MEM_READ_ONLY, dim*sizeof(Float));
    io_driver = (Float*)cl_queue->enqueueMapBuffer(cl_driver, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, dim*sizeof(Float));

    //Create pinned memory for the result:
    cl_result = cl::Buffer(*cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, point_list.size()*sizeof(Float));
    cl_device_result = cl::Buffer(*cl_context, CL_MEM_WRITE_ONLY, point_list.size()*sizeof(Float));
    io_result = (Float*)cl_queue->enqueueMapBuffer(cl_result, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, point_list.size()*sizeof(Float));

    // Register OpenCl function
    kernel_compute = cl::KernelFunctor( cl::Kernel(*cl_program,"kernel_compute"),
                                       *cl_queue,
                                        cl::NullRange,
                                        cl::NDRange(point_list.size()),
                                        cl::NullRange );
}

template<int dim, typename Float>
int DelaunayDecomposition<dim,Float>::find_next_point_opencl(ball_type           *ball,
                                                             const unsigned long  last_index,
                                                             const int            direction,
                                                             Float                &min_time )
{
    // Setup memory for the driver of the ball
    ball->get_driver(io_driver);
    cl_queue->enqueueWriteBuffer(cl_device_driver, CL_TRUE, 0, dim*sizeof(Float), io_driver);

    unsigned long index = ball->any_member();

    // Call OpenCl function
    kernel_compute( cl_point_list_buffer,
                    cl_squared_norm_buffer,
                    cl_center,
                    cl_driver,
                    index,
                    cl_result );

    // Fetch the memory
    cl_queue->enqueueReadBuffer(cl_device_result, CL_TRUE, 0, point_list.size()*sizeof(Float), io_result);

    // Try to find the point that enters the ball next
    long   min = -1;
    switch(direction) {
    case  1: min_time = -std::numeric_limits<Float>::max(); break;
    case -1: min_time =  std::numeric_limits<Float>::max(); break;
    case  0: min_time =  std::numeric_limits<Float>::max(); break;
    }
    for( u_long i=0; i<point_list.size(); ++i ) {
        if(i == last_index || ball->is_member(i)) continue;

        switch(direction) {
        case  1:
            if(io_result[i] < 0.0 && io_result[i] > min_time) {
                min_time = io_result[i];
                min      = i;
            }
            break;
        case -1:
            if(io_result[i] > 0.0 && io_result[i] < min_time) {
                min_time = io_result[i];
                min      = i;
            }
            break;
        case  0:
            if(std::abs(io_result[i]) < std::abs(min_time)) {
                min_time = -io_result[i]; // correct direction of the driver
                min      = i;
            }
            break;
        }
    }
    return min;
}

template<int dim, typename Float>
void DelaunayDecomposition<dim,Float>::shutdown_OpenCL() {
    //cl_queue->enqueueUnmapMemObject(cl_center,io_center);
    //cl_queue->enqueueUnmapMemObject(cl_driver,io_driver);
    //cl_queue->enqueueUnmapMemObject(cl_result,io_result);
    cl_queue->finish();
    delete cl_context;
    delete cl_program;
    delete cl_queue;
}

#endif


} // NAMESPACE

#endif // DELAUNAYDECOMPOSITION_CPP
