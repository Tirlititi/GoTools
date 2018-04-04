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

/* Sparse Polynomial
 *
 * Polynomial structures and manipulations
 * Newton Polytope and Sparse implicitization algorithm
 * 
 * The algorithm of Sparse Resultant computation is a port
 * of the Maple code multires
 * <http://www-sop.inria.fr/galaad/software/multires/>
 */

#ifndef SPARSEPOLYNOMIAL_H
#define SPARSEPOLYNOMIAL_H

#include <assert.h>
#include <vector>
#include <map>
#include <iostream>
#include "GoTools/utils/config.h"
#include "newmat.h"
#include "ConvexPolytope.h"

using std::vector;
using std::map;
using std::string;
using NEWMAT::Matrix;
using NEWMAT::Real;

namespace Go {

struct PowerProduct; // A container of a multi-indexed power product
// Note: T must respect the following conditions:
// - Support the operators *, + (binary) and - (binary and unary)
// - Support an implicit conversion from (int)0 to (T)0 and (int)1 to (T)1
template<typename T> struct Monomial; // Child of PowerProduct, with a coefficient of type T on top of it
template<typename T> struct SparsePolynomial; // A polynomial represented as a list of Monomials
template<typename T> struct SparsePolynomialSystem; // A list of SparsePolynomial representing a polynomial system

struct PowerProduct {

	PowerProduct(vector<int> varpow) : variablePower_(varpow) {}
	PowerProduct(const int varpow[]) : variablePower_() {
		variablePower_.assign(varpow,varpow+sizeof(varpow)/sizeof(int));
	}

	PowerProduct operator*(PowerProduct factor) const;
	bool operator==(PowerProduct other) const;
	bool operator==(vector<int> other) const;
	bool operator<(PowerProduct other) const;
	bool operator<=(PowerProduct other) const {
		return (*this==other) || (*this<other);
	}
	bool operator>(PowerProduct other) const {
		return !(*this<=other);
	}
	bool operator>=(PowerProduct other) const {
		return !(*this<other);
	}

	int GetVariableCount() const {
		return variablePower_.size();
	}
	int GetTotalDegree() const {
		int res = 0;
		for (int i = 0 ; i < GetVariableCount() ; i++)
			res += variablePower_[i];
		return res;
	}
	int GetDegree(int varindex) const {
		if (varindex<0 || varindex>=GetVariableCount())
			return 0;
		return variablePower_[varindex];
	}

	vector<int> variablePower_;
};

template<typename T>
struct Monomial: public PowerProduct {

	Monomial(int numvar, T coef) : PowerProduct(vector<int>(numvar,0)), coefficient_(coef) {}
	Monomial(T coef) : Monomial<T>(0,coef) {}
	Monomial(vector<int> varpow, T coef = 1) : PowerProduct(varpow), coefficient_(coef) {}
	Monomial(PowerProduct varpow, T coef = 1) : PowerProduct(varpow), coefficient_(coef) {}
	Monomial(const int varpow[], T coef = 1) : PowerProduct(varpow), coefficient_(coef) {}

	Monomial<T> operator*(Monomial<T> factor) const {
		return Monomial<T>(static_cast<PowerProduct>(*this) * static_cast<PowerProduct>(factor), coefficient_ * factor.coefficient_);
	}
	Monomial<T>& operator*=(Monomial<T> factor) {
		return *this = *this * factor;
	}
	bool operator==(Monomial<T> comp) const;
	bool operator!=(Monomial<T> comp) const {
		return !(*this==comp);
	}

	// Note: negative powers are not supported, especially by this method
	// The type T does not recquire to support inversion
	T Evaluate(vector<T> varval) const;

	T coefficient_;
	friend Monomial<T> operator*(T factor, Monomial<T> monom) {
		return Monomial<T>(static_cast<PowerProduct>(monom),factor*monom.coefficient_);
	}
};

template<typename T>
struct SparsePolynomial {

	// Even the 0 polynomial should have at least one (constant) term with coefficient 0
	SparsePolynomial(int numvar, T constvalue) : terms_(1,Monomial<T>(numvar,constvalue)) {}
	SparsePolynomial(T constvalue) : terms_(1,Monomial<T>(0,constvalue)) {}
	SparsePolynomial(Monomial<T> uniqueterm) : terms_(1,uniqueterm) {}
	SparsePolynomial() {}

	SparsePolynomial<T> operator+(SparsePolynomial<T> summand) const;
	SparsePolynomial<T> operator-() const;
	SparsePolynomial<T> operator*(SparsePolynomial<T> factor) const;
	SparsePolynomial<T> operator*(T factor) const;
	SparsePolynomial<T> operator+(Monomial<T> summand) const {
		return *this + SparsePolynomial<T>(summand);
	}
	SparsePolynomial<T> operator-(SparsePolynomial<T> summand) const {
		return *this + (-summand);
	}
	SparsePolynomial<T> operator-(Monomial<T> summand) const {
		return *this + (-SparsePolynomial<T>(summand));
	}
	SparsePolynomial<T>& operator+=(SparsePolynomial<T> summand) {
		return *this = *this + summand;
	}
	SparsePolynomial<T>& operator+=(Monomial<T> summand) {
		return *this = *this + SparsePolynomial<T>(summand);
	}
	SparsePolynomial<T>& operator-=(SparsePolynomial<T> summand) {
		return *this = *this - summand;
	}
	SparsePolynomial<T>& operator-=(Monomial<T> summand) {
		return *this = *this - SparsePolynomial<T>(summand);
	}
	SparsePolynomial<T>& operator*=(SparsePolynomial<T> factor) {
		return *this = *this * factor;
	}
	SparsePolynomial<T>& operator*=(T factor) {
		return *this = *this * factor;
	}
	bool operator==(SparsePolynomial<T> comp) const;
	bool operator==(T comp) const {
		bool result = false;
		bool compiszero = comp==0;
		int i, j;
		for (i=0;i<terms_.size();i++) {
			if (terms_[i].coefficient_==0)
				continue;
			if (compiszero)
				return false;
			for (j=0;j<terms_[i].variablePower_.size();j++)
				if (terms_[i].variablePower_[j]!=0)
					return false;
			return terms_[i].coefficient_==comp;
		}
		return compiszero;
	}
	bool operator!=(SparsePolynomial<T> comp) const {
		return !(*this==comp);
	}
	bool operator!=(T comp) const {
		return !(*this==comp);
	}
	Monomial<T>& operator[](int index) {
		assert(index>=0 && index<terms_.size());
		return *(terms_.begin()+index);
	}

	int GetTermCount() const {
		return terms_.size();
	}
	int GetVariableCount() const {
		assert(terms_.size()>0);
		int res = terms_[0].GetVariableCount();
		for (int i=1;i<terms_.size();i++)
			if (res<terms_[i].GetVariableCount())
				res = terms_[i].GetVariableCount();
		/*int res = terms_.begin()->GetVariableCount();
		for (auto& t: terms_)
			if (res<t.GetVariableCount())
				res = t.GetVariableCount();*/
		return res;
	}

	T Evaluate(vector<T> varval) const {
		T res = 0;
		for (int i=0;i<terms_.size();i++)
			res += terms_[i].Evaluate(varval);
		/*for (auto& t: terms_)
			res += t.Evaluate(varval);*/
		return res;
	}

	T GetCoefficient(PowerProduct powervalue) const {
		for (int i=0;i<terms_.size();i++)
			if (static_cast<PowerProduct>(terms_[i])==powervalue)
				return terms_[i].coefficient_;
		return 0;
		//return terms_.at(powervalue);
	}

	int GetDegree(int varindex) const {
		assert(terms_.size()>0);
		int res = terms_[0].GetDegree(varindex);
		for (int i=1;i<terms_.size();i++)
			if (res<terms_[i].GetDegree(varindex))
				res = terms_[i].GetDegree(varindex);
		/*int res = terms_.begin()->GetDegree(varindex);
		for (auto& t: terms_)
			if (res<t.GetDegree(varindex))
				res = t.GetDegree(varindex);*/
		return res;
	}
	int GetLowDegree(int varindex) const {
		assert(terms_.size()>0);
		int res = terms_[0].GetDegree(varindex);
		for (int i=1;i<terms_.size();i++)
			if (res>terms_[i].GetDegree(varindex))
				res = terms_[i].GetDegree(varindex);
		/*int res = terms_.begin()->GetDegree(varindex);
		for (auto& t: terms_)
			if (res>t.GetDegree(varindex))
				res = t.GetDegree(varindex);*/
		return res;
	}
	int GetTotalDegree() const {
		assert(terms_.size()>0);
		int res = terms_[0].GetTotalDegree();
		for (int i=1;i<terms_.size();i++)
			if (res<terms_[i].GetTotalDegree())
				res = terms_[i].GetTotalDegree();
		/*int res = terms_.begin()->GetTotalDegree();
		for (auto& t: terms_)
			if (res<t.GetTotalDegree(varindex))
				res = t.GetTotalDegree(varindex);*/
		return res;
	}

	// The Newton Polytope of P(x_1,...,x_d) = sum_{i = i_1,...,i_d in I} a_i x_1^{i_1} ... x_d^{i_d}, with a_i != 0
	// is the convex hull of I, in the d-dimensional space
	VRepConvexPolytope GetNewtonPolytope(bool sortnodes = false, vector<Pointd>* out_powernode = NULL) const;

	string to_string(char varname = 'x') const;

	// Using map instead of vector may improve the speed of polynomial operations
	vector< Monomial<T> > terms_;
//	map<PowerProduct, T> terms_;

	friend SparsePolynomial<T> operator*(T factor, SparsePolynomial<T> monom) {
		return monom*factor;
	}
	friend ostream& operator<<(ostream& os, const SparsePolynomial<T>& p) {
		return os << p.to_string('y');
	}
};

// TODO: Use monomial basis for computing matrix_symbolic instead of doing polynomial manip
template<typename T> SparsePolynomial<T> Determinant(vector< vector< SparsePolynomial<T> > >& matrix);

template<typename T>
struct SparsePolynomialSystem {

	SparsePolynomialSystem(int numpol, int numvar) : polynomial_(numpol, SparsePolynomial<T>(numvar,0)) {}
	SparsePolynomialSystem() : SparsePolynomialSystem(0,0) {}

	SparsePolynomial<T>& operator[](int index) {
		assert(index>=0 && index<polynomial_.size());
		return polynomial_[index];
	}

	int GetSystemSize() const {
		return polynomial_.size();
	}

	int GetVariableCount() const {
		if (polynomial_.size()==0)
			return 0;
		int res = polynomial_[0].GetVariableCount();
		for (int i=1;i<GetSystemSize();i++)
			if (res<polynomial_[i].GetVariableCount())
				res = polynomial_[i].GetVariableCount();
		return res;
	}

	void AddPolynomial(SparsePolynomial<T> pol) {
		polynomial_.push_back(pol);
	}
	void AddPolynomials(SparsePolynomialSystem<T> polsys) {
		for (int i=0;i<polsys.polynomial_.size();i++)
			polynomial_.push_back(polsys[i]);
	}

	// Compute the Sparse Resultant matrix
	void GetSparseResultantMatrix(vector< vector<T> >* out_matrix, vector<Pointd>* out_minknodes, \
								  Pointd* opt_delta = NULL, Matrix* opt_lifting = NULL) const;

	vector< SparsePolynomial<T> > polynomial_;

private:
	void make_eqns(SparsePolynomialSystem<Real>* out_lhs, SparsePolynomialSystem<Real>* out_conv, SparsePolynomial<Real>* out_obj, vector<int>* out_cd, \
		vector<VRepConvexPolytope> newtonpoly, Matrix lifting) const;
	Pointd start_vertex(vector<VRepConvexPolytope> newtonpoly, SparsePolynomialSystem<Real> lhs, SparsePolynomialSystem<Real> conveq, SparsePolynomial<Real> objective, vector<int> cd, Pointd delta) const;
	vector<int> opt_verts(Pointd vert, SparsePolynomialSystem<Real> lhs, SparsePolynomialSystem<Real> conveq, SparsePolynomial<Real> objective, vector<int> cd) const;
	void compute_rows(vector<Pointd>* out_rowvert, vector<int>* out_rowimax, vector<int>* out_rowjmax, \
		Pointd startvert, vector< vector<Pointd> > polsupport, vector<VRepConvexPolytope> newtonpoly, SparsePolynomialSystem<Real> lhs, SparsePolynomialSystem<Real> conveq, SparsePolynomial<Real> objective, vector<int> cd, Pointd delta) const;
};


template struct Monomial<double>;
template struct SparsePolynomial<double>;
template struct SparsePolynomialSystem<double>;
/*template struct Monomial<int>;
template struct SparsePolynomial<int>;
template struct SparsePolynomialSystem<int>;*/
template struct Monomial< SparsePolynomial<double> >;
template struct SparsePolynomial< SparsePolynomial<double> >;
template struct SparsePolynomialSystem< SparsePolynomial<double> >;
template SparsePolynomial<double> Determinant<double>(vector< vector< SparsePolynomial<double> > >& matrix);

SparsePolynomial<double> GetPolynomialPrimitive(const SparsePolynomial<double>& poly, vector<int> primitivedegree);
SparsePolynomial<double> GetPolynomialDerivative(const SparsePolynomial<double>& poly, vector<int> partialdiffdegree);

// Approximate by another polynomial with only terms belonging to the newpolybasis
// and minimizing the error integral_[0, 1] (P(x) - Q(x))^2 dx
SparsePolynomial<double> Approximate_L2Projection(const SparsePolynomial<double>& poly, vector< vector<int> > newpolybasis);

} // end namespace Go

#endif
