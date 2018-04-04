/* Pointd
* author: Clement Laroche
* SINTEF, Oslo, Norway
*
* Just a few utility functions related to vectors
* and d-dimensional points
* Consider using the vectors (dx1 matrices) from NEWMAT instead
*/

#ifndef POINTD_H
#define POINTD_H

#include <vector>
#include <assert.h>
#include <iostream>
#include "GoTools/utils/config.h"

using std::vector;
using std::ostream;

namespace Go {

// An utility function for vectors
template<typename T>
bool IsElementOfVector(T elem, vector<T> vec) {
	for (int i=0;i<vec.size();i++)
		if (vec[i]==elem)
			return true;
	return false;
}

typedef double PointdReal;
typedef vector<PointdReal> Pointd;

inline int GetPointdDimension(Pointd p) {
	return p.size();
}

// Compute ||p||^2
inline PointdReal GetPointdSquaredNorm(Pointd p) {
	PointdReal res = 0;
	for (int i=0;i<p.size();i++)
		res += p[i]*p[i];
	return res;
}

// Create a point with each coordinate uniformly randomized in [-maxcoef, maxcoef]
inline Pointd GetRandomPointd(int dim, double maxcoef = 1.0, bool intcoef = false) {
	Pointd res;
	int sign;
	for (int i=0;i<dim;i++) {
		sign = (rand()%2)*2-1;
		res.push_back(sign*(double)rand()/(double)RAND_MAX*maxcoef);
		if (intcoef) res[i] = floor(res[i]);
	}
	return res;
}

// Round each coordinate of the point
inline void RoundPointd(Pointd& p) {
	for (int i=0;i<p.size();i++)
		p[i] = round(p[i]);
}

inline Pointd operator+(Pointd l, Pointd r) {
	assert(GetPointdDimension(l)==GetPointdDimension(r));
	Pointd res(GetPointdDimension(l));
	for (int i=0;i<GetPointdDimension(l);i++)
		res[i] = l[i]+r[i];
	return res;
}

inline Pointd& operator+=(Pointd& l, Pointd r) {
	assert(GetPointdDimension(l)==GetPointdDimension(r));
	for (int i=0;i<GetPointdDimension(l);i++)
		l[i] += r[i];
	return l;
}

inline Pointd operator-(Pointd l, Pointd r) {
	assert(GetPointdDimension(l)==GetPointdDimension(r));
	Pointd res(GetPointdDimension(l));
	for (int i=0;i<GetPointdDimension(l);i++)
		res[i] = l[i]-r[i];
	return res;
}

inline Pointd operator-(Pointd p) {
	Pointd res(GetPointdDimension(p));
	for (int i=0;i<GetPointdDimension(p);i++)
		res[i] = -p[i];
	return res;
}

inline Pointd& operator-=(Pointd& l, Pointd r) {
	assert(GetPointdDimension(l)==GetPointdDimension(r));
	for (int i=0;i<GetPointdDimension(l);i++)
		l[i] -= r[i];
	return l;
}

inline Pointd operator*(PointdReal factor, Pointd p) {
	Pointd res(GetPointdDimension(p));
	for (int i=0;i<GetPointdDimension(p);i++)
		res[i] = factor*p[i];
	return res;
}

inline Pointd operator*(Pointd p, PointdReal factor) {
	return factor*p;
}

inline Pointd operator/(Pointd p, PointdReal div) {
	return (1.0f/div)*p;
}

inline Pointd& operator*=(Pointd& p, PointdReal factor) {
	for (int i=0;i<GetPointdDimension(p);i++)
		p[i] *= factor;
	return p;
}

inline Pointd& operator/=(Pointd& p, PointdReal div) {
	return p *= 1.0f/div;
}

inline PointdReal operator*(Pointd l, Pointd r) {
	assert(GetPointdDimension(l)==GetPointdDimension(r));
	PointdReal res = 0;
	for (int i=0;i<GetPointdDimension(l);i++)
		res += l[i]*r[i];
	return res;
}

inline bool operator==(Pointd l, Pointd r) {
	if (GetPointdDimension(l)!=GetPointdDimension(r))
		return false;
	for (int i=0;i<GetPointdDimension(l);i++)
		if (l[i]!=r[i])
			return false;
	return true;
}

inline bool operator!=(Pointd l, Pointd r) {
	return !(l==r);
}

// Lexicographical order
inline bool operator>(Pointd l, Pointd r) {
	int diml = GetPointdDimension(l);
	int dimr = GetPointdDimension(r);
	for (int i=0;i<diml && i<dimr;i++)
		if (l[i]<r[i])
			return true;
		else if (l[i]>r[i])
			return false;
	return diml<dimr;
}

// Lexicographical order
inline bool operator<(Pointd l, Pointd r) {
	return r>l;
}

// Print the coordinates of the point on a line
inline ostream& operator<<(ostream& os, const Pointd& p) {
	for (int i=0;i<GetPointdDimension(p);i++)
		(i==0 ? os : os << " ") << p[i];
	return os;
}

} // namespace Go

#endif
