/*
 * (c) 30.07.2014 Martin Hünniger
 */

#include <iostream>
#include "delaunaydecomposition.h"
#include "io.h"

#include <fstream>
#include <iostream>

#include <ctime>

static const int dim = 16;

typedef DD::DelaunayDecomposition<dim,double> delaunay_type;
typedef delaunay_type::point_type           point_type;
typedef delaunay_type::point_list_type      point_list_type;

point_type generate_point(int dim) {
    point_type point;
    for(int i=0;i<dim;++i) {
        point[i] = static_cast<double>( 20.0*rand()/RAND_MAX-10.0 );//*1e+100;
    }
    return point;
}

point_list_type generate_points(int dim, int num) {
    srand(clock());
    point_list_type points(num);
    for(int i=0; i<num; ++i) {
        points[i] = generate_point(dim);
    }
    return points;
}

void print_usage( int argc, char** argv ) {
    std::cout << "Usage:\n"
              << argv[0] << " <num_points>\n"
              << "num_points is the number of points randomly generated in " << dim <<"D ambient space.\n";
    exit(-1);
}

int main(int argc, char** argv)
{
    if(argc != 2) print_usage(argc,argv);

    unsigned long num_points = atoi(argv[1]);
    point_list_type points = generate_points(dim, num_points);

//    std::string filename("/Users/pirx/Desktop/datapoints2d.txt");
//    std::string testfilename("/Users/pirx/data_for_delaunay.txt");
//    delaunay_type D(filename);

    delaunay_type d(points);

    std::cout << "Generated the points" << std::endl;
    std::cout << "Now computing the Delaunay decomposition" << std::endl;
    clock_t begin = std::clock();
    d.compute();
    clock_t end = std::clock();
    std::cout << d;
    std::cout << "Number of simplices: " << d.num_simplices() << std::endl;
    std::cout << "Computed in: " << double(end-begin) / CLOCKS_PER_SEC << "s\n";

    return 0;
}

