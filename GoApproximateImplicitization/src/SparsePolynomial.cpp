#include "SparsePolynomial.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include "newmatio.h"

#define SP_DEBUG_PRINT	false

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::pow;

// Doesn't really belong there...
Matrix GetRandomMatrix(int rowdim, int coldim) {
	Matrix res(rowdim, coldim);
	int i, j;
	for (i=1;i<=rowdim;i++)
		for (j=1;j<=coldim;j++)
			res(i,j) = (rand() % 1000)+1;
	return res;
}

template<typename T>
T Power(T base, int exponent) {
	T res = 1;
	for (int i=0;i<exponent;i++)
		res *= base;
	return res;
}

template<> float Power<float>(float base, int exponent) {
	return pow(base,exponent);
}

template<> double Power<double>(double base, int exponent) {
	return pow(base,exponent);
}

template<> int Power<int>(int base, int exponent) {
	return pow(base,exponent);
}

namespace Go {

template<typename T>
SparsePolynomial<T> _Determinant_Rec(vector< vector< SparsePolynomial<T> > >& matrix, int row, vector<int>& col) {
	if (col.size()==1)
		return matrix[row][col[0]];
	int i, j, colelem;
	bool neg = true;
	SparsePolynomial<T> result = 0;
	for (i=0;i<col.size();i++) {
		neg = !neg;
		if (matrix[row][col[i]]==0)
			continue;
		colelem = col[i];
		for (j=i;j+1<col.size();j++)
			col[j] = col[j+1];
		col.pop_back();
		if (neg)
			result += (-matrix[row][colelem])*_Determinant_Rec(matrix,row+1,col);
		else
			result += matrix[row][colelem]*_Determinant_Rec(matrix,row+1,col);
		col.push_back(col[col.size()-1]);
		for (j=col.size()-2;j>i;j--)
			col[j] = col[j-1];
		col[i] = colelem;
	}
	return result;
}

template<typename T>
SparsePolynomial<T> Determinant(vector< vector< SparsePolynomial<T> > >& matrix) {
	int s = matrix.size();
	vector<int> col(s);
	for (int i=0;i<s;i++) {
		assert(matrix[i].size()==s);
		col[i] = i;
	}
	if (s==0)
		return 0;
	if (s==1)
		return matrix[0][0];
	return _Determinant_Rec(matrix,0,col);
}

struct _SparseResultantRowResult {
	Pointd vert;
	int imax;
	int jmax;

	_SparseResultantRowResult(Pointd v, int i, int j) : vert(v), imax(i), jmax(j) {}
	bool operator>(_SparseResultantRowResult& other) {
		return vert>other.vert;
	}
	bool operator<(_SparseResultantRowResult& other) {
		return vert<other.vert;
	}
};

PowerProduct PowerProduct::operator*(PowerProduct factor) const {
	int minnumvar = min(GetVariableCount(),factor.GetVariableCount());
	int numvar = max(GetVariableCount(),factor.GetVariableCount());
	vector<int> resultpow(numvar);
	int i;
	for (i=0;i<minnumvar;i++)
		resultpow[i] = variablePower_[i]+factor.variablePower_[i];
	for (;i<GetVariableCount();i++)
		resultpow[i] = variablePower_[i];
	for (;i<factor.GetVariableCount();i++)
		resultpow[i] = factor.variablePower_[i];
	return PowerProduct(resultpow);
}

bool PowerProduct::operator==(PowerProduct other) const {
	return (*this==other.variablePower_);
}

bool PowerProduct::operator==(vector<int> other) const {
	int i;
	for (i=0;i<GetVariableCount() && i<other.size();i++)
		if (variablePower_[i]!=other[i])
			return false;
	while (i<GetVariableCount())
		if (variablePower_[i++]!=0)
			return false;
	while (i<other.size())
		if (other[i++]!=0)
			return false;
	return true;
}

bool PowerProduct::operator<(PowerProduct other) const {
	for (int i=0;i<other.variablePower_.size();i++)
		if (GetDegree(i)<other.variablePower_[i])
			return true;
		else if (GetDegree(i)>other.variablePower_[i])
			return false;
	return false;
}

template<typename T>
bool Monomial<T>::operator==(Monomial<T> comp) const {
	if (coefficient_!=comp.coefficient_)
		return false;
	int i, minvar = min(GetVariableCount(),comp.GetVariableCount());
	for (i=0;i<minvar;i++)
		if (variablePower_[i]!=comp.variablePower_[i])
			return false;
	while (i<GetVariableCount())
		if (variablePower_[i++]!=0)
			return false;
	while (i<comp.GetVariableCount())
		if (comp.variablePower_[i++]!=0)
			return false;
	return true;
}

template<typename T>
T Monomial<T>::Evaluate(vector<T> varval) const {
	assert(varval.size()>=GetVariableCount());
	T res = coefficient_;
	for (int i=0;i<GetVariableCount();i++)
		res *= Power(varval[i],GetDegree(i));
	return res;
}

template<typename T>
SparsePolynomial<T> SparsePolynomial<T>::operator+(SparsePolynomial<T> summand) const {
	int resvarcount = max(GetVariableCount(),summand.GetVariableCount());
	SparsePolynomial<T> result(resvarcount,0);
	vector<bool> newfromsum(summand.terms_.size(),true);
	result.terms_.pop_back();
	int i, j;
	for (i=0;i<terms_.size();i++) {
		if (terms_[i].coefficient_==0)
			continue;
		for (j=0;j<summand.terms_.size();j++)
			if (static_cast<PowerProduct>(terms_[i])==static_cast<PowerProduct>(summand.terms_[j])) {
				result.terms_.push_back(Monomial<T>(terms_[i].variablePower_,terms_[i].coefficient_+summand.terms_[j].coefficient_));
				newfromsum[j] = false;
				break;
			}
		if (j>=summand.terms_.size())
			result.terms_.push_back(Monomial<T>(terms_[i].variablePower_,terms_[i].coefficient_));
	}
	for (i=0;i<summand.terms_.size();i++) {
		if (summand.terms_[i].coefficient_==0)
			continue;
		if (newfromsum[i])
			result.terms_.push_back(Monomial<T>(summand.terms_[i].variablePower_,summand.terms_[i].coefficient_));
	}
	if (result.terms_.size()==0)
		result.terms_.push_back(Monomial<T>(resvarcount,0));
	for (i=0;i<result.terms_.size();i++)
		result.terms_[i].variablePower_.resize(resvarcount,0);
	return result;
}

template<typename T>
SparsePolynomial<T> SparsePolynomial<T>::operator*(SparsePolynomial<T> factor) const {
	int resvarcount = max(GetVariableCount(),factor.GetVariableCount());
	SparsePolynomial<T> result(resvarcount,0);
	Monomial<T> curmon(resvarcount,0);
	result.terms_.pop_back();
	int j,k;
	for (int i=0;i<terms_.size();i++) {
		if (terms_[i].coefficient_==0)
			continue;
		for (j=0;j<factor.terms_.size();j++) {
			if (factor.terms_[j].coefficient_==0)
				continue;
			curmon = terms_[i]*factor.terms_[j];
			for (k=0;k<result.terms_.size();k++)
				if (static_cast<PowerProduct>(curmon)==static_cast<PowerProduct>(result.terms_[k])) {
					result.terms_[k].coefficient_ += curmon.coefficient_;
					break;
				}
			if (k>=result.terms_.size())
				result.terms_.push_back(curmon);
		}
	}
	if (result.terms_.size()==0)
		result.terms_.push_back(Monomial<T>(resvarcount,0));
	return result;
}

template<typename T>
SparsePolynomial<T> SparsePolynomial<T>::operator*(T factor) const {
	if (factor==0)
		return SparsePolynomial<T>(GetVariableCount(),0);
	SparsePolynomial<T> result(*this);
	for (int i=0;i<terms_.size();i++)
		result.terms_[i].coefficient_ *= factor;
	return result;
}

template<typename T>
bool SparsePolynomial<T>::operator==(SparsePolynomial<T> comp) const {
	vector<bool> comptermok(comp.terms_.size(),false);
	int i,j;
	for (i=0;i<terms_.size();i++) {
		if (terms_[i].coefficient_==T(0))
			continue;
		for (j=0;j<comp.terms_.size();j++)
			if (terms_[i]==comp.terms_[j]) {
				comptermok[j] = true;
				break;
			}
		if (j>=comp.terms_.size())
			return false;
	}
	for (j=0;j<comp.terms_.size();j++)
		if (!comptermok[j] && comp.terms_[i].coefficient_!=0)
			return false;
	return true;
}

template<typename T>
SparsePolynomial<T> SparsePolynomial<T>::operator-() const {
	SparsePolynomial<T> result;
	for (int i=0;i<terms_.size();i++)
		result.terms_.push_back(Monomial<T>(terms_[i].variablePower_,-terms_[i].coefficient_));
	return result;
}

SparsePolynomial<double> GetPolynomialPrimitive(const SparsePolynomial<double>& poly, vector<int> primitivedegree) {
	int varcount = poly.GetVariableCount();
	vector<int> varpow(varcount);
	SparsePolynomial<double> result;
	int i,j,k;
	assert(primitivedegree.size()>=varcount);
	for (i=0;i<poly.terms_.size();i++) {
		int invcoef = 1;
		for (j=0;j<varcount;j++) {
			varpow[j] = poly.terms_[i].GetDegree(j)+primitivedegree[j];
			for (k=0;k<primitivedegree[j];k++)
				invcoef *= poly.terms_[i].GetDegree(j)+k+1;
		}
		result += Monomial<double>(varpow,poly.terms_[i].coefficient_*(1/(float)invcoef));
	}
	return result;
}

SparsePolynomial<double> GetPolynomialDerivative(const SparsePolynomial<double>& poly, vector<int> partialdiffdegree) {
	int varcount = poly.GetVariableCount();
	vector<int> varpow(varcount);
	SparsePolynomial<double> result;
	bool zeroterm;
	int i,j,k;
	assert(partialdiffdegree.size()>=varcount);
	for (i=0;i<poly.terms_.size();i++) {
		int coef = 1;
		zeroterm = false;
		for (j=0;j<varcount;j++) {
			varpow[j] = poly.terms_[i].GetDegree(j)-partialdiffdegree[j];
			if (varpow[j]<0) {
				zeroterm = true;
				break;
			}
			for (k=0;k<partialdiffdegree[j];k++)
				coef *= poly.terms_[i].GetDegree(j)-k;
		}
		if (zeroterm)
			continue;
		result += Monomial<double>(varpow,poly.terms_[i].coefficient_*coef);
	}
	if (result.terms_.size()==0)
		result.terms_.push_back(Monomial<double>(varcount,0));
	return result;
}

SparsePolynomial<double> Approximate_L2Projection(const SparsePolynomial<double>& poly, vector<vector<int>> newpolybasis) {
	int matsize = newpolybasis.size();
	int varcount = poly.GetVariableCount();
	vector<int> integrandvect(varcount,1);
	vector<double> upperboundvect(varcount,1);
	Matrix sysmat(matsize,matsize);
	Matrix sysmatrhs(matsize,1);
	int i,j,k,lhsdenom;
	// Solve the linear system in a_i :
	// for 0<=i<matsize,
	// integral_[0,1] x^newpolybasis[i]*this(x) dx = sum_j a_j * integral_[0,1] x^newpolybasis[i]*x^newpolybasis[j] dx
	for (i=0;i<matsize;i++) {
		assert(newpolybasis[i].size()<=varcount);
		SparsePolynomial<double> rhsintegral = GetPolynomialPrimitive(SparsePolynomial<double>(Monomial<double>(newpolybasis[i],1))*(poly),integrandvect);
		sysmatrhs(i+1,1) = rhsintegral.Evaluate(upperboundvect);
		for (j=0;j<matsize;j++) {
			lhsdenom = 1;
			for (k=0;k<varcount;k++)
				lhsdenom *= newpolybasis[i][k]+newpolybasis[j][k]+1;
			sysmat(i+1,j+1) = 1/(float)lhsdenom;
		}
	}
	Matrix resvec(matsize,1);
	resvec = sysmat.i()*sysmatrhs;
	SparsePolynomial<Real> res(varcount,0);
	for (i=0;i<matsize;i++)
		res += Monomial<Real>(newpolybasis[i],resvec(i+1,1));
	return res;
}

template<typename T>
VRepConvexPolytope SparsePolynomial<T>::GetNewtonPolytope(bool sortnodes, vector<Pointd>* out_powernode) const {
	vector<Pointd> nodes;
	int i,j;
	for (i=0;i<terms_.size();i++) {
		if (terms_[i].coefficient_==0)
			continue;
		Pointd node;
		for (j=0;j<terms_[i].GetVariableCount();j++)
			node.push_back(terms_[i].GetDegree(j));
		nodes.push_back(node);
	}
	if (nodes.size()==0)
		nodes.push_back(Pointd(GetVariableCount(),0));
	if (sortnodes)
		std::sort(nodes.begin(),nodes.end()); // descending order
	if (out_powernode)
		*out_powernode = nodes;
	if (sortnodes) {
		VRepConvexPolytope res(nodes);
		std::sort(res.GetExtremalPoints().begin(),res.GetExtremalPoints().end()); // descending order
		return res;
	}
	return VRepConvexPolytope(nodes);
}

template<typename T>
string SparsePolynomial<T>::to_string(char varname) const {
	bool constantonecoef, zeropolynomial = true;
	bool printtimes, printadd = false, printindex = GetVariableCount()>1;
	std::stringstream ss;
	int i, j;
	for (i=0;i<terms_.size();i++) {
		if (terms_[i].coefficient_==0)
			continue;
		zeropolynomial = false;
		printtimes = true;
		if (printadd) {
			std::stringstream coefss;
			if (terms_[i].coefficient_==1 || terms_[i].coefficient_==-1)
				printtimes = false;
			coefss << terms_[i].coefficient_;
			string coefstr = coefss.str();
			if (coefstr.length()>0 && coefstr[0]=='-') {
				ss << " - ";
				if (printtimes)
					ss << coefstr.substr(1);
			} else {
				ss << " + ";
				if (printtimes)
					ss << coefstr;
			}
		} else {
			if (terms_[i].coefficient_==1) {
				printtimes = false;
			} else if (terms_[i].coefficient_==-1) {
				printtimes = false;
				ss << "-";
			} else {
				ss << terms_[i].coefficient_;
			}
		}
		printadd = true;
		constantonecoef = !printtimes;
		for (j=0;j<terms_[i].variablePower_.size();j++) {
			if (terms_[i].variablePower_[j]==0)
				continue;
			if (printtimes)
				ss << "*";
			printtimes = true;
			ss << varname;
			if (printindex)
				ss << "_" << j;
			if (terms_[i].variablePower_[j]==1)
				continue;
			if (terms_[i].variablePower_[j]>0 && terms_[i].variablePower_[j]<10)
				ss << "^" << terms_[i].variablePower_[j];
			else 
				ss << "^{" << terms_[i].variablePower_[j] << "}";
		}
		if (!printtimes && constantonecoef)
			ss << "1";
	}
	if (zeropolynomial)
		return "0";
	return ss.str();
}

template<typename T>
void SparsePolynomialSystem<T>::GetSparseResultantMatrix(vector< vector<T> >* out_matrix, vector<Pointd>* out_minknodes, \
														 Pointd* opt_delta, Matrix* opt_lifting) const {
	Matrix outnum, outdenom;
	Pointd delta;
	Matrix lifting;
	int d = GetSystemSize();
	int dim = GetVariableCount();
	int i,j,k;

	assert(d==dim+1);

	if (opt_delta) delta = *opt_delta;
	else delta = GetRandomPointd(dim, 0.1);
	if (opt_lifting) lifting = *opt_lifting;
	else lifting = GetRandomMatrix(d, d);

	// Simplify the equations by their lower degree
	vector< SparsePolynomial<T> > simplifiedpoly = polynomial_;
	int lowdegree;
	for (i=0;i<d;i++)
		for (j=0;j<dim;j++) {
			lowdegree = simplifiedpoly[i].GetLowDegree(i);
			for (k=0;k<simplifiedpoly[i].terms_.size();k++)
				simplifiedpoly[i].terms_[k].variablePower_[j] -= lowdegree;
		}

	vector<VRepConvexPolytope> newtonpol;
	vector< vector<Pointd> > polypowernodevec(d);
	if (SP_DEBUG_PRINT) cout << "===Newton Polys" << endl;
	for (i=0;i<d;i++) {
		newtonpol.push_back(simplifiedpoly[i].GetNewtonPolytope(true, &polypowernodevec[i]));
		if (SP_DEBUG_PRINT) {
			for (j=0;j<polypowernodevec[i].size();j++)
				cout << "polypowernodevec " << polypowernodevec[i][j] << endl;
			for (j=0;j<newtonpol[i].GetExtremalPointCount();j++)
				cout << "extremal point " << newtonpol[i].GetExtremalPoint(j) << endl;
		}
	}
	if (SP_DEBUG_PRINT) cout << "===" << endl;

	SparsePolynomialSystem<Real> lhs, conveq;
	SparsePolynomial<Real> objective;
	vector<int> pointnum;
	make_eqns(&lhs,&conveq,&objective,&pointnum,newtonpol,lifting);

	Pointd startvert;
	startvert = start_vertex(newtonpol,lhs,conveq,objective,pointnum,delta);

	vector<Pointd> rowvert;
	vector<int> rowimax, rowjmax;
	compute_rows(&rowvert,&rowimax,&rowjmax,startvert,polypowernodevec,newtonpol,lhs,conveq,objective,pointnum,delta);

	// Warning: variables of symbolicmatrix are not the same as lhs, conveq or objective
	int pointnumsum = 0;
	vector<int> incrpointnumsum(polypowernodevec.size(),0);
	for (i=0;i<polypowernodevec.size();i++) {
		incrpointnumsum[i] = pointnumsum;
		pointnumsum += polypowernodevec[i].size();
	}

	Pointd mixednode, mixedpolynode;
	vector<int> varpow(pointnumsum,0);
	int matsize = rowvert.size();
	vector< vector< SparsePolynomial<T> > > symbolicmatrix(matsize,vector< SparsePolynomial<T> >(matsize,Monomial<T>(varpow,0)));

	int pnum;
	for (i=0;i<matsize;i++) {
		vector<Pointd>& minknewtonsummand = newtonpol[rowimax[i]].GetExtremalPoints();
		vector<Pointd>& minksummand = polypowernodevec[rowimax[i]];
		pnum = minksummand.size();
		mixednode = rowvert[i] - minknewtonsummand[rowjmax[i]];

		j = 0;
		for (k=0;k<pnum;k++) {
			mixedpolynode = mixednode + minksummand[k];
			while (j<matsize && mixedpolynode!=rowvert[j])
				j++;
			assert(j<matsize);
			varpow[incrpointnumsum[rowimax[i]]+k] = 1; // C_rowimax[i] X_k
			symbolicmatrix[i][j] = Monomial<T>(varpow,1);
			varpow[incrpointnumsum[rowimax[i]]+k] = 0;
			j++;
		}
	}
	if (SP_DEBUG_PRINT) {
		cout << "===Symbolic Matrix" << endl << endl;
		for (i=0;i<matsize;i++) {
			for (j=0;j<matsize;j++)
				cout << symbolicmatrix[i][j].to_string() << "\t";
			cout << endl;
		}
		cout << "===" << endl;
	}

	vector<T> symbolicval(pointnumsum,0);
	T coefval;
	varpow.assign(dim,0);
	pointnumsum = 0;
	for (i=0;i<d;i++)
		for (j=0;j<polypowernodevec[i].size();j++) {
			for (k=0;k<dim;k++)
				varpow[k] = polypowernodevec[i][j][k];
			coefval = simplifiedpoly[i].GetCoefficient(varpow);
//			coefval = round(coefval);
			symbolicval[pointnumsum++] = coefval; // C_i X_j
		}

	vector< vector<T> > numericmatrix(matsize,vector<T>(matsize));
	for (i=0;i<matsize;i++)
		for (j=0;j<matsize;j++)
			numericmatrix[i][j] = symbolicmatrix[i][j].Evaluate(symbolicval);

	*out_matrix = numericmatrix;
	*out_minknodes = rowvert;
}


// Private

template<typename T>
void SparsePolynomialSystem<T>::make_eqns(SparsePolynomialSystem<Real>* out_lhs, SparsePolynomialSystem<Real>* out_conv, SparsePolynomial<Real>* out_obj, vector<int>* out_cd, \
	vector<VRepConvexPolytope> newtonpoly, Matrix lifting) const {

	int d = GetSystemSize();
	int dim = d-1;
	int cdsum = 0;
	int i, j, k;
	for (i=0;i<d;i++)
		cdsum += newtonpoly[i].GetExtremalPointCount();

	assert(lifting.Nrows()==d);
	assert(lifting.Ncols()==d);

	SparsePolynomialSystem<Real> VertLhsConvEq(dim,cdsum);
	SparsePolynomialSystem<Real> ConvEq(d,cdsum);
	SparsePolynomial<Real> Objective(cdsum,0);
	vector<int> cd(d,0);

	vector<int> lampow(cdsum,0);
	int cdincrsum = 0;
	for (i=0;i<d;i++) {
		vector<Pointd>& ipoly = newtonpoly[i].GetExtremalPoints();
		cd[i] = newtonpoly[i].GetExtremalPointCount();
		SparsePolynomial<Real> lambda(cdsum,0);

		for (j=0;j<cd[i];j++) {
			assert(GetPointdDimension(ipoly[j])==dim);
			lampow[cdincrsum] = 1;
			lambda += Monomial<Real>(lampow);
			for (k=0;k<dim;k++) {
				VertLhsConvEq[k] += ((T)ipoly[j][k])*SparsePolynomial<Real>(lampow);
				Objective += ((T)(lifting(i+1,k+1)*ipoly[j][k]))*SparsePolynomial<Real>(lampow);
			}
			lampow[cdincrsum] = 0;
			cdincrsum++;
		}
		ConvEq[i] = lambda - Monomial<Real>(lampow); // sum lambda_j = 1

		Objective += Monomial<Real>(cdsum,lifting(i+1,d));
	}

	*out_lhs = VertLhsConvEq;
	*out_conv = ConvEq;
	*out_obj = Objective;
	*out_cd = cd;

}

inline
bool IsNextToZero(Real x) {
	return abs(x)<=0.001;
}

template<typename T>
Pointd SparsePolynomialSystem<T>::start_vertex(vector<VRepConvexPolytope> newtonpoly, SparsePolynomialSystem<Real> lhs, SparsePolynomialSystem<Real> conveq, SparsePolynomial<Real> objective, vector<int> cd, Pointd delta) const {

	int dim = GetVariableCount();
	int d = GetSystemSize();
	int i;

	assert(d==dim+1 && newtonpoly.size()==d && GetPointdDimension(delta)==dim);

	Pointd gravity_center(dim,0);
	Pointd current_vert(dim,0);
	Pointd small_shift(dim,-1);
	vector<int> opt_vert_res;
	small_shift[0] = -2;

	for (i=0;i<d;i++)
		gravity_center += newtonpoly[i].GetCenterOfGravity();

	RoundPointd(gravity_center);

	bool foundstartvert;
	current_vert = gravity_center;
	while (true) {

		opt_vert_res = opt_verts(current_vert-delta, lhs, conveq, objective, cd);

		foundstartvert = false;
		for (i=0;i<d;i++)
			if (opt_vert_res[i]>=0) {
				foundstartvert = true;
				break;
			}
		if (foundstartvert)
			break;

		// Update the small_shift (it must go through [-1, 1]^dim)
		for (i=0;i<dim && small_shift[i]==1;i++) {}
		assert(i<dim);
		small_shift[i--] += 1;
		while (i>=0)
			small_shift[i--] = -1;
		current_vert = gravity_center+small_shift;
	}

	return current_vert;
}

vector<Real> MinimizeLinearSystem(SparsePolynomial<Real> obj, SparsePolynomialSystem<Real> cond, vector<int>* feasiblebase = NULL) {
	//  Based on "Introduction to Linear Optimization"
	// of D. Bertsimas and J.N. Tsitsiklis,
	// Chapter 3
	int varcount = obj.GetVariableCount();
	int syssize = cond.GetSystemSize();
	int conddim = varcount-syssize;
	int i, j, k;

//	assert(varcount==cond.GetVariableCount());
	vector<Real> result(varcount,0);

	Matrix linsystemrhs(syssize,1);
	Matrix linsystem(syssize,varcount);
	Matrix basicsystem(syssize,syssize);
	Matrix objsystem(1,varcount);
	Matrix objbasicsystem(1,syssize);
	for (i=1;i<=syssize;i++) {
		for (j=1;j<=varcount;j++)
			linsystem(i,j) = 0;
		linsystemrhs(i,1) = 0;
		objbasicsystem(1,i) = 0;
	}
	for (i=1;i<=varcount;i++)
		objsystem(1,i) = 0;
	int syspow;
	for (i=0;i<syssize;i++) {
		for (j=0;j<cond[i].terms_.size();j++) {
			syspow = -1;
			for (k=0;k<varcount;k++) {
				if (cond[i][j].coefficient_!=0 && cond[i][j].GetDegree(k)>0) {
					assert(cond[i][j].GetDegree(k)==1 && syspow==-1);
					syspow = k;
				}
			}
			if (syspow>=0)
				linsystem(i+1,syspow+1) = cond[i][j].coefficient_;
			else
				linsystemrhs(i+1,1) = cond[i][j].coefficient_;
		}
	}
	for (i=0;i<syssize;i++)
		if (linsystemrhs(i+1,1)<0) {
			linsystemrhs(i+1,1) = -linsystemrhs(i+1,1);
			for (j=0;j<varcount;j++)
				linsystem(i+1,j+1) = -linsystem(i+1,j+1);
		}
	for (j=0;j<obj.terms_.size();j++) {
		syspow = -1;
		for (k=0;k<varcount;k++) {
			if (obj[j].coefficient_!=0 && obj[j].GetDegree(k)>0) {
				assert(obj[j].GetDegree(k)==1 && syspow==-1);
				syspow = k;
			}
		}
		if (syspow>=0)
			objsystem(1,syspow+1) = obj[j].coefficient_;
	}

	vector<int> basicbase;
	vector<int> nonbasicbase;
	vector<int> auxiliarybase(syssize);
	if (feasiblebase==NULL) {
		// No basic solution is given, find one
		vector<int> valpow(varcount+syssize,0);
		SparsePolynomial<Real> auxobj(varcount+syssize,0);
		SparsePolynomialSystem<Real> auxcond(syssize,varcount+syssize);
		for (i=0;i<syssize;i++) {
			auxiliarybase[i] = varcount+i;
			valpow[varcount+i] = 1;
			auxobj += Monomial<Real>(valpow);
			auxcond[i] = cond[i] + Monomial<Real>(valpow);
			valpow[varcount+i] = 0;
		}
		vector<Real> auxresult = MinimizeLinearSystem(auxobj, auxcond, &auxiliarybase);
		if (auxresult.size()==0)
			return vector<Real>();
		for (i=0;i<syssize;i++)
			if (auxresult[varcount+i]>0) {
				return vector<Real>();
			}
		for (i=0;i<syssize;i++)
			if (auxiliarybase[i]>=varcount) {
				int auxpivot = -1;
				Matrix auxbasicsystem(syssize,syssize);
				for (j=1;j<=syssize;j++)
					for (k=1;k<=syssize;k++)
						auxbasicsystem(j,k) = auxiliarybase[k-1]<varcount ? linsystem(j,auxiliarybase[k-1]+1) : (auxiliarybase[k-1]-varcount==j ? 1 : 0);
				for (j=0;j<varcount;j++) {
					Matrix auxdirection = auxbasicsystem.i() * linsystem.SubMatrix(1,syssize,j+1,j+1);
					if (!IsNextToZero(auxdirection(i+1,1)) && !IsElementOfVector(j,auxiliarybase)) {
						auxpivot = j;
						break;
					}
				}
				if (auxpivot<0)
					for (j=0;j<varcount;j++)
						if (!IsElementOfVector(j,auxiliarybase)) {
							auxpivot = j;
							break;
						}
				auxiliarybase[i] = auxpivot;
			}
		feasiblebase = &auxiliarybase;
	}
	for (i=0;i<varcount;i++)
		if (IsElementOfVector(i,*feasiblebase))
			basicbase.push_back(i);
		else
			nonbasicbase.push_back(i);

	for (i=0;i<syssize;i++) {
		for (j=0;j<syssize;j++)
			basicsystem(i+1,j+1) = linsystem(i+1,basicbase[j]+1);
		objbasicsystem(1,i+1) = objsystem(1,basicbase[i]+1);
	}
	Matrix basicsolution = basicsystem.i() * linsystemrhs;
	vector<Matrix> basicdirection(conddim);
	vector<Real> basicreducedcost(conddim);
	while (true) {
		int negativeindex = -1;
		for (i=0;i<conddim;i++) {
			basicdirection[i] = basicsystem.i() * linsystem.SubMatrix(1,syssize,nonbasicbase[i]+1,nonbasicbase[i]+1);
			basicreducedcost[i] = objsystem(1,nonbasicbase[i]+1) - (objbasicsystem * basicdirection[i]).AsScalar();
			if (basicreducedcost[i]<0)
				negativeindex = i;
		}
		if (negativeindex<0) {
			for (j=0;j<syssize;j++)
				result[basicbase[j]] = -basicsolution(j+1,1);
			*feasiblebase = basicbase;
			return result;
		}
		Real theta;
		int basicswapindex = -1;
		for (i=0;i<syssize;i++)
			if (basicdirection[negativeindex](i+1,1)>0) {
				if (basicswapindex<0 || theta<basicsolution(i+1,1)/basicdirection[negativeindex](i+1,1)) {
					theta = basicsolution(i+1,1)/basicdirection[negativeindex](i+1,1);
					basicswapindex = i;
				}
			}
		if (basicswapindex<0) {
			return vector<Real>();
		}
		// Setup next basic data, ie. replace the "basicswapindex" coordinate by the "negativeindex" coordinate
		objbasicsystem(1,basicswapindex+1) = objsystem(1,nonbasicbase[negativeindex]+1);
		for (i=1;i<=syssize;i++)
			basicsystem(i,basicswapindex+1) = linsystem(i,nonbasicbase[negativeindex]+1);
		basicsolution = basicsystem.i() * linsystemrhs;
		int tmpnbb = nonbasicbase[negativeindex];
		nonbasicbase[negativeindex] = basicbase[basicswapindex];
		basicbase[basicswapindex] = tmpnbb;
	}

	return vector<Real>();
}

template<typename T>
vector<int> SparsePolynomialSystem<T>::opt_verts(Pointd vert, SparsePolynomialSystem<Real> lhs, SparsePolynomialSystem<Real> conveq, SparsePolynomial<Real> objective, vector<int> cd) const {

	int dim = GetVariableCount();
	int d = GetSystemSize();
	int i, j;

	assert(GetPointdDimension(vert)==dim && dim==lhs.GetSystemSize());

	SparsePolynomialSystem<Real> lhseqvert = lhs;
	for (i=0;i<dim;i++)
		lhseqvert[i] -= Monomial<Real>(lhseqvert.GetVariableCount(),vert[i]);

	lhseqvert.AddPolynomials(conveq);

	Pointd minimalvert = MinimizeLinearSystem(objective, lhseqvert);

	vector<int> result(d,-1);
	if (minimalvert.size()==0)
		return result;
	
	assert(cd.size()==d);

	int cdsum = 0;
	for (i=0;i<d;i++) {
		// Variables are (X_k = lam_i x_j) with 0<=k<sum(cd), 0<=i<d and 0<=j<cd[i]
		for (j=0;j<cd[i];j++) {
			assert(minimalvert.size()>cdsum);
			if (IsNextToZero(minimalvert[cdsum]-1))
				result[i] = j;
			cdsum++;
		}
	}

	return result;
}

template<typename T>
void SparsePolynomialSystem<T>::compute_rows(vector<Pointd>* out_rowvert, vector<int>* out_rowimax, vector<int>* out_rowjmax, \
		Pointd startvert, vector< vector<Pointd> > polsupport, vector<VRepConvexPolytope> newtonpoly, SparsePolynomialSystem<Real> lhs, SparsePolynomialSystem<Real> conveq, SparsePolynomial<Real> objective, vector<int> cd, Pointd delta) const {

	vector<Pointd> todo(1,startvert);
	vector<Pointd> processlist(1,startvert);
	Pointd processed, vt1, vt2;
	int dim = GetVariableCount();
	int d = GetSystemSize();

	vector<_SparseResultantRowResult> row;
	vector<Pointd> rowvert;
	vector<int> rowimax;
	vector<int> rowjmax;
	int i, j, imax, jmax;

	while (todo.size()>0) {
		processed = todo.back();
		todo.pop_back();

		assert(GetPointdDimension(processed)==dim);

		vector<int> opt_vert_res = opt_verts(processed-delta, lhs, conveq, objective, cd);

		for (imax=d-1;imax>=0 && opt_vert_res[imax]<0;imax--) {}
		assert(imax>=0);
		jmax = opt_vert_res[imax];

		row.push_back(_SparseResultantRowResult(processed,imax,jmax));

		vt2 = processed - newtonpoly[imax].GetExtremalPoint(jmax);

		vector<Pointd>& ipoly = polsupport[imax];
		
		for (i=0;i<ipoly.size();i++) {
			vt1 = vt2 + ipoly[i];
			if (!IsElementOfVector(vt1,processlist)) {
				processlist.push_back(vt1);
				todo.push_back(vt1);
			}
		}
	}

	std::sort(row.rbegin(),row.rend()); // descending order
	rowvert.resize(row.size());
	rowimax.resize(row.size());
	rowjmax.resize(row.size());
	for (i=0;i<row.size();i++) {
		rowvert[i] = row[i].vert;
		rowimax[i] = row[i].imax;
		rowjmax[i] = row[i].jmax;
	}

	*out_rowvert = rowvert;
	*out_rowimax = rowimax;
	*out_rowjmax = rowjmax;
}

} // namespace GoTools
