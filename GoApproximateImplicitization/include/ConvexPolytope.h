/*
* Copyright (C) 2018 SINTEF ICT,
* Applied Mathematics, Norway.
*
* Contact information: E-mail: tor.dokken@sintef.no
* SINTEF ICT, Department of Applied Mathematics,
* P.O. Box 124 Blindern,
* 0314 Oslo, Norway.
*
* This file is part of GoTools.
*
* GoTools is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version. 
*
* GoTools is distributed in the hope that it will be useful,        
* but WITHOUT ANY WARRANTY; without even the implied warranty of         
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public
* License along with GoTools. If not, see
* <http://www.gnu.org/licenses/>.
*
* In accordance with Section 7(b) of the GNU Affero General Public
* License, a covered work must retain the producer line in every data
* file that is created or manipulated using GoTools.
*
* Other Usage
* You can be released from the requirements of the license by purchasing
* a commercial license. Buying such a license is mandatory as soon as you
* develop commercial activities involving the GoTools library without
* disclosing the source code of your own applications.
*
* This file may be used in accordance with the terms contained in a
* written agreement between you and SINTEF ICT. 
*/

/* Convex Polytope
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
