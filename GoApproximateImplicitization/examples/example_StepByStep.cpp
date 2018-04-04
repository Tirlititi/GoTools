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

/* Example of Sparse Resultant computation and implicitization */

#include "SparsePolynomial.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

using std::fstream;
using std::ios;
using std::cin;
using std::cout;
using std::endl;
using std::vector;
using namespace Go;

int main(int argc, const char* argv[]) {
	int i, j;

	cout << "FIRST EXAMPLE: SPARSE RESULTANT\n\n";
	SparsePolynomialSystem<double> polysys;
	SparsePolynomial<double> pol;
	vector<int> powval(2);
	// Fill the system with polynomials
	powval[0] = 1;	powval[1] = 0;	pol = Monomial<double>(powval,1.0);		// pol =  x
	powval[0] = 0;	powval[1] = 1;	pol += Monomial<double>(powval,3.0);	// pol += 3y
	polysys.AddPolynomial(pol);
	powval[0] = 1;	powval[1] = 0;	pol = Monomial<double>(powval,2.0);		// pol =  2x
	powval[0] = 0;	powval[1] = 1;	pol += Monomial<double>(powval,-1.0);	// pol += -y
	powval[0] = 2;	powval[1] = 2;	pol += Monomial<double>(powval,5.0);	// pol += 5 x^2 y^2
	powval[0] = 1;	powval[1] = 1;	pol += Monomial<double>(powval,1.0);	// pol += x y
	polysys.AddPolynomial(pol);
	// You may also avoid using a vector for constructing the monomials
	pol = Monomial<double>({0,0},1.0);										// pol = 1
	pol += Monomial<double>({1,0},1.0);										// pol += x
	pol += Monomial<double>({0,1},1.0);										// pol += y
	polysys.AddPolynomial(pol);
	vector< vector<double> > det;
	vector<Pointd> minknodes;
	polysys.GetSparseResultantMatrix(&det, &minknodes);
	cout << "Minkovski Nodes:" << endl;
	for (i=0;i<minknodes.size();i++) {
		for (j=0;j<GetPointdDimension(minknodes[i]);j++)
			cout << minknodes[i][j] << " ";
		cout << endl;
	}
	cout << "Sparse Matrix:" << endl;
	for (i=0;i<det.size();i++) {
		for (j=0;j<det.size();j++)
			cout << det[i][j] << " ";
		cout << endl;
	}
	cout << "FIRST EXAMPLE DONE!" << endl;

	char whatever;
	cin >> whatever;

	cout << "SECOND EXAMPLE: IMPLICITIZATION\n\n";
	SparsePolynomialSystem< SparsePolynomial<double> > paramsys;
	SparsePolynomial< SparsePolynomial<double> > parameq;
	//  A polynomial whose coefficients are polynomials allow to split the variables in two groups:
	// the ones used for Resultant computation and the ones that remain symbolic parameters
	parameq = Monomial< SparsePolynomial<double> >({0,0}, Monomial<double>({1,0,0}));		// x = ...
	parameq -= Monomial< SparsePolynomial<double> >({1,0}, Monomial<double>({0,0,0}, 1.0));	// x = s
	parameq -= Monomial< SparsePolynomial<double> >({0,1}, Monomial<double>({0,0,0}, 3.0));	// x = s + 3t
	paramsys.AddPolynomial(parameq);
	parameq = Monomial< SparsePolynomial<double> >({0,0}, Monomial<double>({0,1,0}));		// y = ...
	parameq -= Monomial< SparsePolynomial<double> >({1,0}, Monomial<double>({0,0,0}, 2.0));	// y = 2s
	parameq -= Monomial< SparsePolynomial<double> >({0,1}, Monomial<double>({0,0,0},-1.0));	// y = 2s - t
	parameq -= Monomial< SparsePolynomial<double> >({2,2}, Monomial<double>({0,0,0}, 5.0));	// y = 2s - t + 5 s^2 t^2
	parameq -= Monomial< SparsePolynomial<double> >({1,1}, Monomial<double>({0,0,0}, 1.0));	// y = 2s - t + 5 s^2 t^2 + s t
	paramsys.AddPolynomial(parameq);
	parameq = Monomial< SparsePolynomial<double> >({0,0}, Monomial<double>({0,0,1}));		// z = ...
	parameq -= Monomial< SparsePolynomial<double> >({0,0}, Monomial<double>({0,0,0}, 1.0));	// z = 1
	parameq -= Monomial< SparsePolynomial<double> >({1,0}, Monomial<double>({0,0,0}, 1.0));	// z = 1 + s
	parameq -= Monomial< SparsePolynomial<double> >({0,1}, Monomial<double>({0,0,0}, 1.0));	// z = 1 + s + t
	paramsys.AddPolynomial(parameq);
	cout << "Equations:" << endl;
	cout << paramsys[0].to_string('t') << endl;
	cout << paramsys[1].to_string('t') << endl;
	cout << paramsys[2].to_string('t') << endl << endl;

	vector< vector< SparsePolynomial<double> > > spmatrix;
	paramsys.GetSparseResultantMatrix(&spmatrix, &minknodes);
	cout << "Minkovski Nodes:" << endl;
	for (i=0;i<minknodes.size();i++) {
		for (j=0;j<GetPointdDimension(minknodes[i]);j++)
			cout << minknodes[i][j] << " ";
		cout << endl;
	}
	cout << "Sparse Implicit Matrix:" << endl;
	for (i=0;i<spmatrix.size();i++) {
		for (j=0;j<spmatrix.size();j++)
			cout << spmatrix[i][j].to_string() << " ";
		cout << endl;
	}

	cin >> whatever;

	cout << "Implicit Equation: " << Determinant(spmatrix) << endl;
	cout << "SECOND EXAMPLE DONE!" << endl;

	cin >> whatever;
    return 0;
}


