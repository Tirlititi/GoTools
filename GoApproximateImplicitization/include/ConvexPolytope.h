/* Convex Polytope
* author: Clement Laroche
* SINTEF, Oslo, Norway
*
* n-dimensional polytope representation
* Convex hull algorithm
*/

#ifndef CONVEXPOLYTOPE_H
#define CONVEXPOLYTOPE_H

#include <assert.h>
#include <vector>
#include "Pointd.h"
#include "GoTools/utils/config.h"

using std::vector;

namespace Go {

struct ConvexPolytope; // Base structure for polytope representations
struct HRepConvexPolytope; // Hyperplane-based representation of a polytope \! Not implemented !/
struct VRepConvexPolytope; // Vertex-based representation of a polytope

enum ConvexPolytopeRepresentation {
	CONVEX_POLYTOPE_HREP,
	CONVEX_POLYTOPE_VREP
};

struct ConvexPolytope {


	ConvexPolytopeRepresentation rep_;
};

struct HRepConvexPolytope : public ConvexPolytope {

};

struct VRepConvexPolytope : public ConvexPolytope {

	VRepConvexPolytope(vector<Pointd> pointcloud); // Compute the convex hull of the point cloud

	Pointd GetCenterOfGravity();
	int GetExtremalPointCount() {
		return extremal_points_.size();
	}
	Pointd GetExtremalPoint(int index) {
		assert(index>=0 && index<extremal_points_.size());
		return extremal_points_[index];
	}
	vector<Pointd>& GetExtremalPoints() {
		return extremal_points_;
	}

	vector<Pointd> extremal_points_;
};

}

#endif
