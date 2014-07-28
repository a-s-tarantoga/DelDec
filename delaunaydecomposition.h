#ifndef DELAUNAYDECOMPOSITION_H
#define DELAUNAYDECOMPOSITION_H

//#include "configure.h"
#include "point.h"
#include "ball.h"
#include "marks.h"
#include "helper.h"

#include <vector>
#include <queue>
#include <map>
#include <iostream>
#include <initializer_list>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <limits>
#include <string>
#include <fstream>

#include <thread>
#include <future>

#include <CL/cl.hpp>

#ifdef USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

namespace NAMESPACE {

template<int dim=3,typename Float=double>
class DelaunayDecomposition {
public: // typedefs
    typedef Point<dim,Float>           point_type;
    typedef std::vector<point_type>    point_list_type;
    typedef std::vector<Float>         float_list_type;
    typedef Ball<dim,Float>            ball_type;
    typedef std::vector<ball_type*>    ball_list_type;
    typedef Marks::mark_type           mark_type;
    //typedef std::map<mark_type, unsigned long, Marks::Cmp>
    typedef std::map<mark_type, unsigned long, std::less<mark_type>>
    mark_map_type;
    typedef std::deque<ball_type*>     agenda_type;
    typedef std::vector<u_long>        simplex_type;
    typedef std::vector<simplex_type>  simplex_list_type;

public: // construction & destruction
    DelaunayDecomposition()
        : point_list(0)
        , needs_update(true)
    {
        ASSERT( dim > 1 );
    }

    DelaunayDecomposition(const std::vector<point_type> &list)
        : point_list(list.size())
        , needs_update(true)
    {
        std::copy(list.begin(),
                  list.end(),
                  point_list.begin());
    }

    DelaunayDecomposition(const std::string& filename)
        : point_list(0)
        , needs_update(true)
    {
        std::ifstream file(filename);
        if(!file.good()) {
            std::cerr << "Error opening " << filename << std::endl;
            exit(-1);
        }
        srand(clock());
        while(!file.eof()) {
            if(file.peek() != -1) {
                point_type point;
                file >> point;
                //point.perturb(0.000001);
                point_list.push_back(point);
            }
        }
        file.close();
    }

    ~DelaunayDecomposition()
    {}

public: // access
    bool is_empty()
    {
        return point_list.size() == 0;
    }

    void compute() {
        ASSERT(point_list.size() > dim);
#ifdef USE_OPENCL
        init_openCL();
#endif
        marks.clear();
        uint start_index = 0;
        ball_type *initial_ball;
        do {
            initial_ball = find_initial_simplex(start_index);
            ++start_index;
        } while(!initial_ball && start_index < point_list.size());

        ASSERT(initial_ball);

        if(!initial_ball) {
            std::cerr << "An initial ball could not be found ! :-(\n";
            exit(-2);
        }

        find_all_simplices(initial_ball);
        needs_update = false;
    }

    unsigned long num_simplices() const {
        return simplex_list.size();
    }

private:

    ball_type *find_initial_simplex(uint start_index) {
        ASSERT(start_index < point_list.size());
        if(!needs_update) return NULL;

        // find two points of maximal distance
        Float max_distance = Float(0.0);
        long  max = -1;
        for(auto j=0; j<point_list.size(); ++j) {
            Float dist_squared = Float(0.0);
            for(auto i=0; i<dim; ++i) {
                dist_squared += sqr(point_list[start_index][i] - point_list[j][i]);
            }
            ON_DEBUG(std::cout << "point: " << j << " dist_squared: " << dist_squared << std::endl);
            if(max_distance < dist_squared) {
                max = j;
                max_distance = dist_squared;
            }
        }
        ON_DEBUG(std::cout << "max: " << max << " max_distance: " << max_distance << std::endl);
        srand(clock());

        // set up the starting point
        point_type start_center;
        for(auto i=0; i<dim; ++i) {
            start_center[i] = static_cast<Float>((point_list[start_index][i] + point_list[max][i])/2.0);
        }
        start_center.perturb(0.001*sqrt(max_distance));
        ON_DEBUG( std::cout << "start_center: " << start_center << std::endl);

        // compute the driver for start_point
        Float min_distance = std::numeric_limits<Float>::max();
        long closest_point = -1;
        for(auto j=0;j<point_list.size();++j) {
            Float dist_squared = Float(0.0);
            for(auto i=0;i<dim;++i) {
                dist_squared += sqr(start_center[i]-point_list[j][i]);
            }
            ON_DEBUG(std::cout << "point: " << j << " dist_squared: " << dist_squared << std::endl);
            if(min_distance > dist_squared) {
                closest_point = j;
                min_distance = dist_squared;
            }
        }
        ON_DEBUG(std::cout << "Selected closest_point: " << closest_point << std::endl);

        // Define the first ball
        ball_type *ball = new ball_type( point_list,
                                         start_center,
                                         closest_point);
        ball->update_driver();
        ON_DEBUG(ball->print_matrices(1));
        ON_DEBUG(std::cout << "Inserting points...\n");

        while((ball->size()<point_list.size()) &&
              (ball->size()<=dim)) {
            // find the point that enters the ball
            Float min_time = Float(0.0);
            auto next_point = find_next_point( ball,
                                               -1,
                                               0,
                                               min_time );
            ON_DEBUG( std::cout << "Inserting point: " << next_point
                      << " with t_p " << min_time << std::endl );
            ball->update_center( min_time );
            ball->add_point( next_point );

            ON_DEBUG( ball->print_center() );
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
        ON_DEBUG(
                simplex_type simplex = ball->simplex();
                std::cout << " => [ ";
                for(auto vertex : simplex) {
                    std::cout << vertex << " ";
                }
                std::cout << "]\n");
        ON_DEBUG(std::cout << "\nInitial simplex found!\n\n");
        return ball;
    }

    /**
     * @brief find_all_simplices  This routine explores the whole space in depth first manner to find all
     * simplices which circumsenters are critical points.
     * It also determines the infinite maximum.
     * It explores the whole pointset for maxima in a depth first manner
     * and stores them for later exploration.
     * @param start_ball the initial ball representing an initial simplex
     * of the Delaunay decomposition.
     */
    void find_all_simplices(ball_type *start_ball);

    void thread_function(std::promise<ball_list_type> & ball_list,
                         ball_type *ball);

    Float sqr(Float a) { return a*a; }

    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template <typename T>
    int pos_neg(T val) {
        return (T(0) <= val) ? 1 : -1;
    }

    /**
     * @brief explore
     * finds all the neighboring simplices to a given ball
     * @param ball
     * @return
     */
    ball_list_type explore(ball_type *ball);

    bool explored(ball_type* ball)
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

    /**
     * @brief squared_distance
     * Computes the squared distance between center and the point with index in the point list
     * @param center
     * @param index
     * @return
     */
    Float squared_distance( const Float   *center,
                            unsigned long  index ) const;


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
    int find_next_point( ball_type           *ball,
                         const unsigned long  last_index,
                         const int            direction,
                         Float               &min_time );

    /**
     * @brief find_next_point_helper
     * Does the actual computeation for find_next_point.
     * @param ball
     * @param last_index
     * @param min_time
     * @return
     */
    template<int number>
    int find_next_point_helper( ball_type           *ball,
                                const unsigned long  last_index,
                                Float               &min_time );

    /**************************************************************************
     * Optimized Routines
     **************************************************************************/

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
    long find_next_point_opt(ball_type *ball,
                             const unsigned long   last_index,
                             const int    direction,
                             Float &min_time );

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
    template<int number>
    long find_next_point_helper_opt(ball_type* ball,
                                    const unsigned long last_index,
                                    Float &min_time );

    Float time_to_enter_ball(ball_type *ball,
                             const unsigned long index)
    {
        Float res = ball->time_to_enter_ball(index, ball->any_member());
        ON_DEBUG( std::cout << "For index: " << index << " time to enter ball: " << res << std::endl );
        return res;
    }

    /**************************************************************************
     * OpenCL stuff
     **************************************************************************/
#ifdef USE_OPENCL
    void init_openCL() {
        typedef std::vector<cl::Platform> platform_list_type;
        typedef std::vector<cl::Device>   device_list_type;


        //const std::string kernel_source =
        std::ostringstream kernel_source;
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
               "                                                                                            \n"
//               "__kernel void store_value( const Float value,                                               \n"
//               "                           const int   index,                                               \n"
//               "                            __global Float* dest)                                           \n"
//               "{                                                                                           \n"
//               "    dest[index] = value;                                                                    \n"
//               "}                                                                                           \n"
               "                                                                                            \n"
               "__kernel void kernel_compute( __constant Float *point_list,                                 \n"
               "                              __constant Float *squared_norm_list,                          \n"
               "                              __global Float *center,                                       \n"
               "                              __global Float *driver,                                       \n"
               "                              const unsigned long index,                                    \n"
               "                              const int dim,                                                \n"
               "                              __global Float *result)                                       \n"
               "{                                                                                           \n"
               "    int point_i_base  = get_global_id(0) * dim;             // point to test                \n"
               "    int point_j_base  = index * dim;                        // any point from the ball      \n"
               "    int id            = get_global_id(0);                                                   \n"
               "                                                                                            \n"
               "    Float p_minus_q[dim];                                                                   \n"
               "#pragma unroll                                                                              \n"
               "    for(int i=0;i< " << dim << ";++i) {                                                     \n"
               "        p_minus_q[i] = point_list[point_j_base + i] - point_list[point_i_base + i];         \n"
               "    }                                                                                       \n"
               "                                                                                            \n"
               "    Float x_times_pmq = 0.0;                                                                \n"
               "    Float d_times_pmq = 0.0;                                                                \n"
               "#pragma unroll                                                                              \n"
               "    for(int i=0; i<" << dim <<"; ++i) {                                                     \n"
               "        x_times_pmq += center[i] * p_minus_q[i];                                            \n"
               "        d_times_pmq += driver[i] * p_minus_q[i];                                            \n"
               "    }                                                                                       \n"
               "    x_times_pmq *= 2.0;                                                                     \n"
               "    d_times_pmq *= 2.0;                                                                     \n"
               "                                                                                            \n"
               "    result[id] = (squared_norm_list[id] - squared_norm_list[index] + x_times_pmq)/d_times_pmq;\n"
               "    //result[id] = p_minus_q[0];                                                            \n"
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

        //try to get all supported decies
        device_list_type all_devices;
        default_platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);

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

        typedef typename std::pair<const char*,unsigned long> src_type;
        cl::Program::Sources sources;
        sources.push_back(src_type(kernel_source.str().c_str(), kernel_source.str().length()));

        cl_program = new cl::Program(*cl_context,sources);

        if(cl_program->build(default_device_vector)!=CL_SUCCESS) {
            std::cout << "Program did not build!" << std::endl;
            std::cout << cl_program->getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
            exit(1);
        }

        // point_list_buffer stores all points on the OpenCL device
        unsigned long point_list_buffer_size = point_list.size() * dim * sizeof(Float);
        Float point_buffer[point_list.size()*dim];
        cl_point_list_buffer = new cl::Buffer(*cl_context, CL_MEM_READ_ONLY, point_list_buffer_size);
        unsigned long counter = 0;
        for(u_long i=0; i<point_list.size(); ++i) {
            for(u_long j=0; j<dim; ++j) {
                point_buffer[counter] = point_list[i][j];
                ++counter;
            }
        }

        // squared_norm_buffer stores all the squared norms of all points
        unsigned long squared_norm_buffer_size = point_list.size() * sizeof(Float);
        cl_squared_norm_buffer = new cl::Buffer(*cl_context, CL_MEM_READ_ONLY, squared_norm_buffer_size);
        Float norm_buffer[point_list.size()];
        for(u_long i=0; i<point_list.size(); ++i) {
            norm_buffer[i] = point_list[i].squared_norm();
        }

        cl_queue = new cl::CommandQueue(*cl_context,default_device);
        cl_queue->enqueueWriteBuffer(*cl_point_list_buffer, CL_TRUE, 0, point_list_buffer_size, point_buffer);
        cl_queue->enqueueWriteBuffer(*cl_squared_norm_buffer, CL_TRUE, 0, squared_norm_buffer_size, norm_buffer);
    }


    int find_next_point_opencl(ball_type           *ball,
                               const unsigned long  last_index,
                               const int            direction,
                               Float                &min_time )
    {
        Float  center[dim];
        Float  driver[dim];
        //Float  *result;
        Float  result[point_list.size()];
        u_long index = ball->any_member();
        int    dimension = dim;

        ball->get_center(center);
        ball->get_driver(driver);
        u_long buffer_size = dim * sizeof(Float);
        u_long result_size = point_list.size() * sizeof(Float);

        cl::Buffer cl_center(*cl_context, CL_MEM_READ_WRITE, buffer_size);
        cl::Buffer cl_driver(*cl_context, CL_MEM_READ_WRITE, buffer_size);
        cl::Buffer cl_result(*cl_context, CL_MEM_READ_WRITE, result_size);

//TIME_SECTION(
        cl_queue->enqueueWriteBuffer(cl_center, CL_TRUE, 0, buffer_size, center);
        cl_queue->enqueueWriteBuffer(cl_driver, CL_TRUE, 0, buffer_size, driver);
//);
        // Call OpenCl stuff
        cl::KernelFunctor kernel_compute(cl::Kernel(*cl_program,"kernel_compute"),
                                         *cl_queue,
                                         cl::NullRange,
                                         cl::NDRange(point_list.size()),
                                         cl::NullRange);
//        cl::KernelFunctor store_value(cl::Kernel(*cl_program,"store_value"),
//                                      *cl_queue,
//                                      cl::NullRange,
//                                      cl::NDRange(1),
//                                      cl::NullRange);

//TIME_SECTION(
//#pragma unroll
//        for(int i=0;i<dim;++i) {
//            store_value(center[i],i,cl_center);
//            store_value(driver[i],i,cl_driver);
//        }
//)
        kernel_compute(*cl_point_list_buffer,
                       *cl_squared_norm_buffer,
                        cl_center,
                        cl_driver,
                        index,
                        dim,
                        cl_result);

        //result = (Float*)cl_queue->enqueueMapBuffer(cl_result, CL_TRUE, 0, 0, result_size);
        cl_queue->enqueueReadBuffer(cl_result, CL_TRUE, 0, result_size, result);


        // Find the right value:
        long   min = -1;
        switch(direction) {
        case  1: min_time = -std::numeric_limits<Float>::max(); break;
        case -1: min_time =  std::numeric_limits<Float>::max(); break;
        case  0: min_time =  std::numeric_limits<Float>::max(); break;
        }
        for( u_long i=0; i<point_list.size(); ++i ) {
            if(i == last_index || ball->is_member(i)) continue;

            ON_DEBUG( std::cout << "r[" << i << "] = " << result[i] << ", " );

            switch(direction) {
            case  1:
                if(result[i] < 0.0 && result[i] > min_time) {
                    min_time = result[i];
                    min      = i;
                }
                break;
            case -1:
                if(result[i] > 0.0 && result[i] < min_time) {
                    min_time = result[i];
                    min      = i;
                }
                break;
            case  0:
                if(std::abs(result[i]) < std::abs(min_time)) {
                    min_time = -result[i]; // correct direction of the driver
                    min      = i;
                }
                break;
            }
        }
        ON_DEBUG( std::cout << std::endl);
        return min;
    }
#endif







    /**************************************************************************
     * IO
     **************************************************************************/

    void sort_simplices() {
        std::sort(simplex_list.begin(), simplex_list.end());
    }

    friend std::ostream& operator<<(std::ostream &os, DelaunayDecomposition& D ) {
        for(auto point : D.point_list) {
            os << point << std::endl;
        }
        D.sort_simplices();
        for(auto simplex : D.simplex_list) {
//            for(auto vertex : simplex) {
//                os << vertex << " ";
//            }
            os << simplex << std::endl;
        }
        return os;
    }

    friend std::ofstream& operator<<( std::ofstream& ofs, DelaunayDecomposition& D) {
        for(auto point : D.point_list) {
            ofs << point << std::endl;
        }
        D.sort_simplices();
        for(auto simplex : D.simplex_list) {
            ofs << simplex << " ";
        }
        return ofs;
    }

    friend std::ofstream& print_point_list( std::ofstream& ofs, DelaunayDecomposition& D) {
        for(auto point : D.point_list) {
            ofs << point << std::endl;
        }
        return ofs;
    }

private:
    DelaunayDecomposition( const DelaunayDecomposition& );
    DelaunayDecomposition& operator=( const DelaunayDecomposition& );

private:
    point_list_type   point_list;
    float_list_type   squared_norm_list;
    simplex_list_type simplex_list;
    mark_map_type     marks;
    bool              needs_update;

#ifdef USE_THREADS
    std::mutex        mutex;
#endif

#ifdef USE_OPENCL
    // OpenCL stuff
    cl::Context      *cl_context;
    cl::CommandQueue *cl_queue;
    cl::Program      *cl_program;
    cl::Buffer       *cl_point_list_buffer;
    cl::Buffer       *cl_squared_norm_buffer;
#endif
};

} // NAMESPACE

#include "delaunaydecomposition.cpp"

#endif // DELAUNAYDECOMPOSITION_H
