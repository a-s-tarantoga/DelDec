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
        if( point_list.empty() ) {
            return;
        }
#ifdef USE_OPENCL
        init_OpenCL();
#endif
        marks.clear();
        simplex_list.clear();

        uint start_index = 0;
        ball_type *initial_ball;

        if(point_list.size() == 1) {
            initial_ball->add_point(0);
            simplex_list.push_back(initial_ball->simplex());
            return;
        }

        do {
            initial_ball = find_initial_simplex(start_index);
            ++start_index;
        } while(!initial_ball && start_index < point_list.size());

        if(!initial_ball) {
            std::cerr << "An initial ball could not be found !\n";
            exit(-2);
        }

        find_all_simplices(initial_ball);
#ifdef USE_OPENCL
        shutdown_OpenCL();
#endif
    }

    unsigned long num_simplices() const {
        return simplex_list.size();
    }

private:

    ball_type *find_initial_simplex(uint start_index);

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

    /**
     * @brief explore
     * finds all the neighboring simplices to a given ball
     * @param ball
     * @return
     */
    ball_list_type explore(ball_type *ball);

    /**
     * @brief explored
     * Finds out whether the ball has been already explored or not
     * @param ball the ball under consideration
     * @return true if already explored false otherwise
     */
    bool explored(ball_type* ball);

    /**
     * @brief time_to_enter_ball
     * Computes the "time" or distance at which the point enters the ball as its center is moved
     * in direction of the driver and keeping the vertices of the subspan of the ball fixed
     * @param ball the ball under consideration
     * @param index the index of a point in point_list
     * @return a Float representing the "time" or the "distance" the center of the ball has to be moved
     */
    Float time_to_enter_ball(ball_type *ball,
                             const unsigned long index)
    {
        Float res = ball->time_to_enter_ball(index, ball->any_member());
        return res;
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
     * @param ball the ball under consideration
     * @param last_index the point that was removed from the ball as last
     * @param direction the sign of the driver determined by the coefficient
     * of last_index in the conve combination of the center
     * @param min_time the "time" or "distance" at which the new point enters the ball
     * @return the index of the new point
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
     * @param ball the ball under consideration
     * @param last_index the point that was removed from the ball as last
     * @param min_time the "time" or "distance" at which the new point enters the ball
     * @return the index of the new point
     */
    template<int number>
    long find_next_point_helper_opt(ball_type* ball,
                                    const unsigned long last_index,
                                    Float &min_time );

    /**************************************************************************
     * OpenCL stuff
     **************************************************************************/
#ifdef USE_OPENCL
    /**
     * @brief init_openCL
     * Procedure to set up the context for computing on accelerated devices
     */
    void init_OpenCL();

    /**
     * @brief find_next_point_opencl
     * Perform the computation of the point that enters the ball next in direction of the balls
     * driver by using the accelerated hardware as GPUs
     * @param ball represents the actual ball
     * @param last_index the point removed lastly from the ball
     * @param direction the sign of the driver defined by the coefficient in the conve combination of the center
     * @param min_time the "time" or "distance" at which the point enters the bal
     * @return the index of the point that enters the ball next
     */
    int find_next_point_opencl(ball_type           *ball,
                               const unsigned long  last_index,
                               const int            direction,
                               Float                &min_time );

    void shutdown_OpenCL();
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
    cl::Buffer        cl_point_list_buffer;
    cl::Buffer        cl_squared_norm_buffer;
    cl::Buffer        cl_center;
    cl::Buffer        cl_device_center;
    Float            *io_center;
    cl::Buffer        cl_driver;
    cl::Buffer        cl_device_driver;
    Float            *io_driver;
    cl::Buffer        cl_result;
    cl::Buffer        cl_device_result;
    Float            *io_result;
    cl::KernelFunctor kernel_compute;
#endif
};

} // NAMESPACE

#include "delaunaydecomposition.cpp"

#endif // DELAUNAYDECOMPOSITION_H
