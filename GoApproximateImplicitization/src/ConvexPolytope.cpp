#include "ConvexPolytope.h"

#include "newmat.h"

using std::vector;
using std::abs;
using NEWMAT::Matrix;

namespace Go {

Pointd GetHyperplaneNormal(vector<Pointd>& pointcloud, vector<int>& hyperplane) {
	Pointd& basepoint = pointcloud[hyperplane[0]];
	int dim = GetPointdDimension(basepoint);
	Pointd res(dim);
	Matrix m(dim-1,dim-1);
	int i,j,k;
	int sign = 1;
	for (i=0;i<dim;i++) {
		for (j=1;j<dim;j++)
			for (k=0;k+1<dim;k++)
				m(j,k+1) = pointcloud[hyperplane[j]][k<i ? k : k+1]-basepoint[k<i ? k : k+1];
		res[i] = sign*m.Determinant();
		sign = -sign;
	}
	return res;
}

void QuickHullProceedHyperplane(vector<int>& convexhull, vector<Pointd>& pointcloud, vector<int>& proceedindex, vector<int> hyperplane, bool hyperplaneside) {
	if (proceedindex.size()==0)
		return;
	int dim = GetPointdDimension(pointcloud[0]);
	assert(dim>0 && dim==hyperplane.size());
	PointdReal maxdist = 0, curdist;
	vector<int> hyperplanehalfpoint;
	int i, j, maxdistindex = -1;
	Matrix hyperplanematrix(dim,dim);
	Pointd& basepoint = pointcloud[hyperplane[0]];
	assert(GetPointdDimension(basepoint)==dim);

	for (i=1;i<dim;i++) {
		assert(GetPointdDimension(pointcloud[hyperplane[i]])==dim);
		for (j=0;j<dim;j++)
			hyperplanematrix(i+1,j+1) = pointcloud[hyperplane[i]][j]-basepoint[j];
	}
	Pointd hyperplanenormal = GetHyperplaneNormal(pointcloud,hyperplane);
	PointdReal hyperplaneconstant = -(hyperplanenormal*basepoint);
	vector<int> proceededpoints;
	for (i=0;i<proceedindex.size();i++) {
		assert(GetPointdDimension(pointcloud[proceedindex[i]])==dim);
		for (j=0;j<dim;j++)
			hyperplanematrix(1,j+1) = pointcloud[proceedindex[i]][j]-basepoint[j];
		PointdReal det = hyperplanematrix.Determinant();
		if ((hyperplaneside && det>0) || (!hyperplaneside && det<0)) {
			hyperplanehalfpoint.push_back(proceedindex[i]);
			proceededpoints.push_back(proceedindex[i]);
			curdist = abs(hyperplanenormal*pointcloud[proceedindex[i]] + hyperplaneconstant);
			if (curdist>maxdist || (curdist==maxdist && maxdistindex>=0 && GetPointdSquaredNorm(pointcloud[proceedindex[i]]-basepoint)>GetPointdSquaredNorm(pointcloud[maxdistindex]-basepoint))) {
				maxdist = curdist;
				maxdistindex = proceedindex[i];
			}
		}
	}
	if (maxdistindex==-1)
		return;

	convexhull.push_back(maxdistindex);
	for (i=0;i<hyperplanehalfpoint.size();i++)
		if (hyperplanehalfpoint[i]==maxdistindex) {
			hyperplanehalfpoint.erase(hyperplanehalfpoint.begin()+i);
			break;
		}
	vector<int> newhyperplane = hyperplane;
	for (i=0;i<dim;i++) {
		newhyperplane[i] = maxdistindex;
		QuickHullProceedHyperplane(convexhull,pointcloud,hyperplanehalfpoint,newhyperplane,hyperplaneside);
		newhyperplane[i] = hyperplane[i];
	}
	for (i=0;i<proceedindex.size();i++)
		if (IsElementOfVector(proceedindex[i],proceededpoints))
			proceedindex.erase(proceedindex.begin()+(i--));
}

VRepConvexPolytope::VRepConvexPolytope(vector<Pointd> pointcloud) {
	if (pointcloud.size()==0)
		return;
	int dim = GetPointdDimension(pointcloud[0]);
	int i,j,k;
	assert(dim>0);

	// If pointcloud.size()<=dim, use another algorithm or project the point cloud on a smaller space
	assert(pointcloud.size()>dim);

	vector<int> convexhull;

	// First find [dim] extremal points defining an hyperplane
	Pointd mincoord(pointcloud[0]), maxcoord(pointcloud[0]);
	vector<int> minindex(dim,0), maxindex(dim,0);
	for (i=1;i<pointcloud.size();i++) {
		assert(GetPointdDimension(pointcloud[i])==dim);
		for (j=0;j<dim;j++) {
			if (pointcloud[i][j]<mincoord[j]) {
				mincoord[j] = pointcloud[i][j];
				minindex[j] = i;
			} else if (pointcloud[i][j]==mincoord[j]) { // Handle degeneracy (vertices inside a supporting hyperplane)
				k = 0;
				while (k<dim && pointcloud[i][k]==pointcloud[minindex[j]][k])
					k++;
				if (k<dim && pointcloud[i][k]<pointcloud[minindex[j]][k])
					minindex[j] = i;
			}
			if (pointcloud[i][j]>maxcoord[j]) {
				maxcoord[j] = pointcloud[i][j];
				maxindex[j] = i;
			} else if (pointcloud[i][j]==maxcoord[j]) {
				k = 0;
				while (k<dim && pointcloud[i][k]==pointcloud[maxindex[j]][k])
					k++;
				if (k<dim && pointcloud[i][k]>pointcloud[maxindex[j]][k])
					maxindex[j] = i;
			}
		}
	}
	for (j=0;j<dim && convexhull.size()<dim;j++) {
		if (!IsElementOfVector(minindex[j],convexhull))
			convexhull.push_back(minindex[j]);
		if (!IsElementOfVector(maxindex[j],convexhull))
			convexhull.push_back(maxindex[j]);
	}

	// Verify that enough points were found to form a starting hyperplane
	assert(convexhull.size()>=dim);

	// Run the QuickHull algorithm
	vector<int> hyperplane(dim);
	for (i=0;i<dim;i++)
		hyperplane[i] = convexhull[i];

	vector<int> proceedindex;
	proceedindex.reserve(pointcloud.size());
	for (i=0;i<pointcloud.size();i++)
		if (!IsElementOfVector(i,hyperplane))
			proceedindex.push_back(i);

	QuickHullProceedHyperplane(convexhull,pointcloud,proceedindex,hyperplane,true);
	QuickHullProceedHyperplane(convexhull,pointcloud,proceedindex,hyperplane,false);

	extremal_points_.resize(convexhull.size());
	for (i=0;i<convexhull.size();i++)
		extremal_points_[i] = pointcloud[convexhull[i]];
}

Pointd VRepConvexPolytope::GetCenterOfGravity() {
	assert(extremal_points_.size()>0);
	Pointd res = extremal_points_[0];
	for (int i=1;i<extremal_points_.size();i++)
		res += extremal_points_[i];
	res /= extremal_points_.size();
	return res;
}

} // namespace Go
