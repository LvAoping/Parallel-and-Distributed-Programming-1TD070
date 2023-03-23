#ifndef _BUFFON_LAPLACE_H_
#define _BUFFON_LAPLACE_H_
#include<stdlib.h>

/**
 * This program estimates the value of PI using a Buffon-Laplace simulation.
 * @author Malin Kallen 2019-2020
 */

#define PI 3.14159265358979323846 // Needed for generation of random angles

struct point {
	double x;
	double y;
};
typedef struct point point_t;

/**
 * Decide whether an needle with the eye inside the rectangle [0 A]x[0 B]
 * crosses the boundary of the rectangle.
 * @param point Location of the point of the needle
 * @param A X value of the upper boundary of the rectangle
 * @param B Y value of the upper boundary of the rectangle
 */
int crossing(point_t eye, point_t point, double A, double B, int Nx, int Ny);

/**
 * Repeatedly drop a needle of length L at a random position on a surface with a
 * grid of equally spaced lines at distances A and B in the x and y direction
 * respectively. Count the number of times the needle crosses a grid line.
 * @param A Grid line space in x
 * @param B Grid line space in y
 * @param L Length of the needle
 * @param N Number of times to drop the needle
 * @return Number of times the needle crosses a grid line
 */
int run_simulation(double A, double B, double L, int N, int Nx, int Ny);

/**
 * Throw a needle of length L such that the end with an eye is within the
 * rectangle [0 A]x[0 B] and the angle of the needle is completely random.
 * @param A Width of the rectangle in which the eye will land
 * @param B Height of the rectangle in which the point will land
 * @param L Length of the needle
 * @param eye Position of the eye after throw
 * @param point Position of the point after throw
 */
void throw_needle(double L, point_t *eye, point_t *point, int Nx, int Ny);

#endif
