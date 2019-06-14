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

inline int	TriToPointSub(
	double pnt[3],		// (i) Coordinate value of point
	double vtx1[3],		// (i) vertex 1
	double vtx2[3],		// (i) vertex 2
	double vtx3[3],		// (i) vertex 3
	double norml[3],	// (i) Unit normal vector
	double dir12[3],	// (i) Unit vector from vertices 1 to 2
	double leng12,		// (i) distance from vertex 1 to 2
	double dir23[3],	// (i) Unit vector from vertices 2 to 3
	double leng23,		// (i) distance from vertex 2 to 3
	double dir31[3],	// (i) Unit vector from vertex 3 to 1
	double leng31,		// (i) distance from vertex 3 to 1
	double npnt[3],		// (o) nearest point
	double * dist_p		// (o) output the square of the nearest distance)
)
{

	if (leng12 < ZERO_VALUE && leng23 < ZERO_VALUE)
	{
		npnt[0] = vtx1[0];
		npnt[1] = vtx1[1];
		npnt[2] = vtx1[2];
		*dist_p = 
			(pnt[0] - npnt[0]) * (pnt[0] - npnt[0]) +
			(pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) +
			(pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	if (leng12 < ZERO_VALUE)
	{
		LineToPoint(pnt, vtx1, vtx3, dir23, leng23, npnt, dist_p);
		return 0;
	}
	if (leng23 < ZERO_VALUE) {
		LineToPoint(pnt, vtx1, vtx2, dir12, leng12, npnt, dist_p);
		return 0;
	}
	if (DOT(norml, norml) < ZERO_VALUE) 
	{
		if (DOT(dir12,dir31) > 0.0) {
			LineToPoint(pnt, vtx2, vtx3, dir23, leng23, npnt, dist_p);
		}
		else {
			if (leng12 > leng31) {
				LineToPoint(pnt, vtx1, vtx2, dir12, leng12, npnt, dist_p);
			}
			else {
				LineToPoint(pnt, vtx3, vtx1, dir31, leng31, npnt, dist_p);
			}
		}
		return 0;
	}

	double	vec_diff[3];

	vec_diff[0] = pnt[0] - vtx1[0];
	vec_diff[1] = pnt[1] - vtx1[1];
	vec_diff[2] = pnt[2] - vtx1[2];

	const double dotprd = DOT(norml, vec_diff);
	double	plane[3];
	plane[0] = vec_diff[0] - dotprd * norml[0];
	plane[1] = vec_diff[1] - dotprd * norml[1];
	plane[2] = vec_diff[2] - dotprd * norml[2];
	double	plane_on[3];
	plane_on[0] = vtx1[0] + plane[0];
	plane_on[1] = vtx1[1] + plane[1];
	plane_on[2] = vtx1[2] + plane[2];

	const double dotprd1 = DOT(plane, dir12);
	double	outer[3];
	outer[0] = plane[1] * dir12[2] - plane[2] * dir12[1];
	outer[1] = plane[2] * dir12[0] - plane[0] * dir12[2];
	outer[2] = plane[0] * dir12[1] - plane[1] * dir12[0];
	const double RL1 = DOT(norml, outer);

	plane[0] = plane_on[0] - vtx2[0];
	plane[1] = plane_on[1] - vtx2[1];
	plane[2] = plane_on[2] - vtx2[2];
	const double dotprd2 = DOT(plane, dir23);
	outer[0] = plane[1] * dir23[2] - plane[2] * dir23[1];
	outer[1] = plane[2] * dir23[0] - plane[0] * dir23[2];
	outer[2] = plane[0] * dir23[1] - plane[1] * dir23[0];
	const double RL2 = DOT(norml, outer);

	plane[0] = plane_on[0] - vtx3[0];
	plane[1] = plane_on[1] - vtx3[1];
	plane[2] = plane_on[2] - vtx3[2];
	const double dotprd3 = DOT(plane, dir31);
	outer[0] = plane[1] * dir31[2] - plane[2] * dir31[1];
	outer[1] = plane[2] * dir31[0] - plane[0] * dir31[2];
	outer[2] = plane[0] * dir31[1] - plane[1] * dir31[0];
	const double RL3 = DOT(norml, outer);

	if (RL1 <= 0.0 && RL2 <= 0.0 && RL3 <= 0.0) 
	{
		npnt[0] = plane_on[0];
		npnt[1] = plane_on[1];
		npnt[2] = plane_on[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	else if (dotprd1 >= 0.0 && dotprd1 <= leng12 && RL1 >= 0.0) 
	{
		LineToPoint(pnt, vtx1, vtx2, dir12, leng12, npnt, dist_p);
		return 0;
	}
	else if (dotprd2 >= 0.0 && dotprd2 <= leng23 && RL2 >= 0.0) 
	{
		LineToPoint(pnt, vtx2, vtx3, dir23, leng23, npnt, dist_p);
		return 0;
	}
	else if (dotprd3 >= 0.0 && dotprd3 <= leng31 && RL3 >= 0.0) 
	{
		LineToPoint(pnt, vtx3, vtx1, dir31, leng31, npnt, dist_p);
		return 0;
	}
	else if (dotprd1 <= 0.0 && dotprd3 >= leng31) 
	{
		npnt[0] = vtx1[0];
		npnt[1] = vtx1[1];
		npnt[2] = vtx1[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	else if (dotprd2 <= 0.0 && dotprd1 >= leng12) 
	{
		npnt[0] = vtx2[0];
		npnt[1] = vtx2[1];
		npnt[2] = vtx2[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	else if (dotprd3 <= 0.0 && dotprd2 >= leng23) 
	{
		npnt[0] = vtx3[0];
		npnt[1] = vtx3[1];
		npnt[2] = vtx3[2];
		*dist_p = (pnt[0] - npnt[0]) *  (pnt[0] - npnt[0]) + (pnt[1] - npnt[1]) * (pnt[1] - npnt[1]) + (pnt[2] - npnt[2]) * (pnt[2] - npnt[2]);
		return 0;
	}
	return -1;
}

inline int	TriToPoint(
	double pnt[3],		// (i) Coordinate value of point
	double vertex1[3],	// (i) vertex 1
	double vertex2[3],	// (i) vertex 2
	double vertex3[3],	// (i) vertex 3
	double npnt[3],		// (o) nearest point
	double * dist_p		// (o) nearest distance squared)
)
{

	double	dir1[3];
	double	dir2[3];
	double	dir3[3];
	dir1[0] = vertex2[0] - vertex1[0];	dir1[1] = vertex2[1] - vertex1[1];	dir1[2] = vertex2[2] - vertex1[2];
	dir2[0] = vertex3[0] - vertex2[0];	dir2[1] = vertex3[1] - vertex2[1];	dir2[2] = vertex3[2] - vertex2[2];
	dir3[0] = vertex1[0] - vertex3[0];	dir3[1] = vertex1[1] - vertex3[1];	dir3[2] = vertex1[2] - vertex3[2];

	double	normal[3];
	normal[0] = dir1[2] * dir3[1] - dir1[1] * dir3[2];
	normal[1] = dir1[0] * dir3[2] - dir1[2] * dir3[0];
	normal[2] = dir1[1] * dir3[0] - dir1[0] * dir3[1];
	const double length = sqrt(DOT(normal, normal));
	if (length > ZERO_VALUE)
	{
		normal[0] = normal[0] / length;
		normal[1] = normal[1] / length;
		normal[2] = normal[2] / length;
	}

	double	udir1[3];
	double	udir2[3];
	double	udir3[3];

	const double length1 = sqrt(DOT(dir1, dir1));
	if (length1 > ZERO_VALUE)
	{
		udir1[0] = dir1[0] / length1;
		udir1[1] = dir1[1] / length1;
		udir1[2] = dir1[2] / length1;
	}
	const double length2 = sqrt(DOT(dir2, dir2));
	if (length2 > ZERO_VALUE)
	{
		udir2[0] = dir2[0] / length2;
		udir2[1] = dir2[1] / length2;
		udir2[2] = dir2[2] / length2;
	}
	const double length3 = sqrt(DOT(dir3, dir3));
	if (length3 > ZERO_VALUE)
	{
		udir3[0] = dir3[0] / length3;
		udir3[1] = dir3[1] / length3;
		udir3[2] = dir3[2] / length3;
	}

	TriToPointSub(pnt, vertex1, vertex2, vertex3, normal,
		udir1, length1, udir2, length2,
		udir3, length3, npnt, dist_p);
	return 0;
}

#ifdef __cplusplus
};
#endif

#endif

