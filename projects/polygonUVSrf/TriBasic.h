#ifndef _TRIBASIC_H
#define _TRIBASIC_H
#define	ZERO_VALUE	(1.0e-15)

#define DOT(a,b)	(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#ifdef __cplusplus
extern "C" {
#endif

inline int	LineToPoint(
	double pnt[3],		// (i) Coordinate value of point
	double bound1[3],	// (i) Line segment end point 1
	double bound2[3],	// (i) Line segment end point 2
	double dir[3],		// (i) Unit direction vector of line segment
	double len,			// (i) Line segment length
	double npnt[3],		// (o) nearest point
	double * dist_p		// (o) nearest distance squared)
)
{
	if (len > ZERO_VALUE) 
	{
		double	inner = (pnt[0] - bound1[0]) * dir[0] + (pnt[1] - bound1[1]) * dir[1] + (pnt[2] - bound1[2]) * dir[2];
		if (inner <= 0.0)
		{
			npnt[0] = bound1[0];
			npnt[1] = bound1[1];
			npnt[2] = bound1[2];
		}
		else if (inner >= len)
		{
			npnt[0] = bound2[0];
			npnt[1] = bound2[1];
			npnt[2] = bound2[2];
		}
		else {
			npnt[0] = bound1[0] + inner * dir[0];
			npnt[1] = bound1[1] + inner * dir[1];
			npnt[2] = bound1[2] + inner * dir[2];
		}
	}
	else {
		npnt[0] = bound1[0];
		npnt[1] = bound1[1];
		npnt[2] = bound1[2];
	}

	*dist_p = (pnt[0] - npnt[0]) *(pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) *(pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) *(pnt[2] - npnt[2]);
	return 0;
}

inline int	TriToPoint(
	double pnt[3],		// (i) Coordinate value of point
	double vertex1[3],		// (i) vertex 1
	double vertex2[3],		// (i) vertex 2
	double vertex3[3],		// (i) vertex 3
	double npnt[3],		// (o) nearest point
	double * dist_p		// (o) output the square of the nearest distance)
)
{
	double dir12[3];
	double dir23[3];
	double dir31[3];

	dir12[0] = vertex2[0] - vertex1[0];
	dir12[1] = vertex2[1] - vertex1[1];	
	dir12[2] = vertex2[2] - vertex1[2];
	
	dir23[0] = vertex3[0] - vertex2[0];
	dir23[1] = vertex3[1] - vertex2[1];
	dir23[2] = vertex3[2] - vertex2[2];
	
	dir31[0] = vertex1[0] - vertex3[0];
	dir31[1] = vertex1[1] - vertex3[1];
	dir31[2] = vertex1[2] - vertex3[2];

	const double length12 = sqrt(DOT(dir12, dir12));
	if (length12 > ZERO_VALUE)
	{
		dir12[0] = dir12[0] / length12;
		dir12[1] = dir12[1] / length12;
		dir12[2] = dir12[2] / length12;
	}
	const double length23 = sqrt(DOT(dir23, dir23));
	if (length23 > ZERO_VALUE)
	{
		dir23[0] = dir23[0] / length23;
		dir23[1] = dir23[1] / length23;
		dir23[2] = dir23[2] / length23;
	}
	const double length31 = sqrt(DOT(dir31, dir31));
	if (length31 > ZERO_VALUE)
	{
		dir31[0] = dir31[0] / length31;
		dir31[1] = dir31[1] / length31;
		dir31[2] = dir31[2] / length31;
	}

	double	normal[3];
	normal[0] = dir12[2] * dir31[1] - dir12[1] * dir31[2];
	normal[1] = dir12[0] * dir31[2] - dir12[2] * dir31[0];
	normal[2] = dir12[1] * dir31[0] - dir12[0] * dir31[1];
	const double length = sqrt(DOT(normal, normal));
	if (length > ZERO_VALUE)
	{
		normal[0] = normal[0] / length;
		normal[1] = normal[1] / length;
		normal[2] = normal[2] / length;
	}

	if (length12 < ZERO_VALUE && length23 < ZERO_VALUE)
	{
		npnt[0] = vertex1[0];
		npnt[1] = vertex1[1];
		npnt[2] = vertex1[2];
		*dist_p = 
			(pnt[0] - npnt[0]) * (pnt[0] - npnt[0]) +
			(pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) +
			(pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	if (length12 < ZERO_VALUE)
	{
		LineToPoint(pnt, vertex1, vertex3, dir23, length23, npnt, dist_p);
		return 0;
	}
	if (length23 < ZERO_VALUE) {
		LineToPoint(pnt, vertex1, vertex2, dir12, length12, npnt, dist_p);
		return 0;
	}
	if (DOT(normal, normal) < ZERO_VALUE)
	{
		if (DOT(dir12,dir31) > 0.0) {
			LineToPoint(pnt, vertex2, vertex3, dir23, length23, npnt, dist_p);
		}
		else {
			if (length12 > length31) {
				LineToPoint(pnt, vertex1, vertex2, dir12, length12, npnt, dist_p);
			}
			else {
				LineToPoint(pnt, vertex3, vertex1, dir31, length31, npnt, dist_p);
			}
		}
		return 0;
	}

	double	vec_diff[3];

	vec_diff[0] = pnt[0] - vertex1[0];
	vec_diff[1] = pnt[1] - vertex1[1];
	vec_diff[2] = pnt[2] - vertex1[2];

	const double dotprd = DOT(normal, vec_diff);
	double	plane[3];
	plane[0] = vec_diff[0] - dotprd * normal[0];
	plane[1] = vec_diff[1] - dotprd * normal[1];
	plane[2] = vec_diff[2] - dotprd * normal[2];
	double	plane_on[3];
	plane_on[0] = vertex1[0] + plane[0];
	plane_on[1] = vertex1[1] + plane[1];
	plane_on[2] = vertex1[2] + plane[2];

	const double dotprd1 = DOT(plane, dir12);
	double	outer[3];
	outer[0] = plane[1] * dir12[2] - plane[2] * dir12[1];
	outer[1] = plane[2] * dir12[0] - plane[0] * dir12[2];
	outer[2] = plane[0] * dir12[1] - plane[1] * dir12[0];
	const double RL1 = DOT(normal, outer);

	plane[0] = plane_on[0] - vertex2[0];
	plane[1] = plane_on[1] - vertex2[1];
	plane[2] = plane_on[2] - vertex2[2];
	const double dotprd2 = DOT(plane, dir23);
	outer[0] = plane[1] * dir23[2] - plane[2] * dir23[1];
	outer[1] = plane[2] * dir23[0] - plane[0] * dir23[2];
	outer[2] = plane[0] * dir23[1] - plane[1] * dir23[0];
	const double RL2 = DOT(normal, outer);

	plane[0] = plane_on[0] - vertex3[0];
	plane[1] = plane_on[1] - vertex3[1];
	plane[2] = plane_on[2] - vertex3[2];
	const double dotprd3 = DOT(plane, dir31);
	outer[0] = plane[1] * dir31[2] - plane[2] * dir31[1];
	outer[1] = plane[2] * dir31[0] - plane[0] * dir31[2];
	outer[2] = plane[0] * dir31[1] - plane[1] * dir31[0];
	const double RL3 = DOT(normal, outer);

	if (RL1 <= 0.0 && RL2 <= 0.0 && RL3 <= 0.0) 
	{
		npnt[0] = plane_on[0];
		npnt[1] = plane_on[1];
		npnt[2] = plane_on[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	else if (dotprd1 >= 0.0 && dotprd1 <= length12 && RL1 >= 0.0) 
	{
		LineToPoint(pnt, vertex1, vertex2, dir12, length12, npnt, dist_p);
		return 0;
	}
	else if (dotprd2 >= 0.0 && dotprd2 <= length23 && RL2 >= 0.0) 
	{
		LineToPoint(pnt, vertex2, vertex3, dir23, length23, npnt, dist_p);
		return 0;
	}
	else if (dotprd3 >= 0.0 && dotprd3 <= length31 && RL3 >= 0.0) 
	{
		LineToPoint(pnt, vertex3, vertex1, dir31, length31, npnt, dist_p);
		return 0;
	}
	else if (dotprd1 <= 0.0 && dotprd3 >= length31) 
	{
		npnt[0] = vertex1[0];
		npnt[1] = vertex1[1];
		npnt[2] = vertex1[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	else if (dotprd2 <= 0.0 && dotprd1 >= length12) 
	{
		npnt[0] = vertex2[0];
		npnt[1] = vertex2[1];
		npnt[2] = vertex2[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	else if (dotprd3 <= 0.0 && dotprd2 >= length23) 
	{
		npnt[0] = vertex3[0];
		npnt[1] = vertex3[1];
		npnt[2] = vertex3[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	return -1;
}


#ifdef __cplusplus
};
#endif

#endif

