/*
 * (c) 30.07.2014 Martin Hünniger
 */


/*
 * Compute the following formula:
 *
 *       ||p_i||^2 - ||p_j||^2 + 2< x, p_j-p_i >
 * t_p = ---------------------------------------
 *                 2< d, p_j-p_i >
 * i     is the global index
 * index refers to point p_j
 */
__kernel void kernel_compute( __global double *point_list,
                              __global double *squared_norm_list,
                              __global double *center,
                              __global double *driver,
                              const unsigned long index,
                              const int dim
                              __global double *result)
{
    int point_i_base  = index * dim;         // any point from the ball
    int point_j_base  = global_id(0) * dim;  // point to test
    int id            = global_id(0);

    double p_minus_q[dim];
    for(int i=0; i<dim; ++i) {
        p_minus_q[i] = point_list[point_j_base + i] - point_list[point_i_base + i];
    }

    double x_times_pmq = 0.0;
    double d_times_pmq = 0.0;
    for(int i=0; i<dim; ++i) {
        x_times_pmq += center[i] * p_minus_q[i];
        d_times_pmq += driver[i] * p_minus_q[i];
    }
    x_times_pmq *= 2.0;
    d_times_pmq *= 2.0;

    result[id] = (squared_norm[id] - squared_norm[index] + x_times_pmq)/d_times_pmq;
}
