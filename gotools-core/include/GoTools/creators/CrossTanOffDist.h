/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
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

#ifndef _CROSSTANOFFDIST_
#define _CROSSTANOFFDIST_

#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>
#include <memory>


namespace Go
{

/// This class defines an evaluator-based offset-curve.
/// We're blending between a set of curves seen as an evaluator
/// based curve.
class CrossTanOffDist : public EvalCurve
{
public:

    /// Constructor
    /// \param poscurve the curve to offset from.
    /// \param tangcv1 tangent curve along poscurve.
    /// \param tangcv2 cross tangent curve along poscurve.
    /// \param blend1 1-dimensional blending function for tangcv1.
    /// \param blend2 1-dimensional blending function for tangcv2.
    /// \param opposite1 space curve corresponding to poscurve,
    ///                  parametrized in the opposite direction.
    /// \param opposite2 space curve corresponding to poscurve,
    ///                  parametrized in the opposite direction.
    ///                  May be equal to opposite1.
    /// \param factor if != 1.0 the poscurve, opposite1 & opposite2
    ///               will be used to define the opposite curve, so
    ///               as to minimize differences in parametrization
    ///               along the cv.
    CrossTanOffDist(shared_ptr<SplineCurve>& poscurve,
		    shared_ptr<SplineCurve>& tangcv1,
		    shared_ptr<SplineCurve>& tangcv2,
		    shared_ptr<SplineCurve>& blend1,
		    shared_ptr<SplineCurve>& blend2,
		    shared_ptr<SplineCurve>& opposite1,
		    shared_ptr<SplineCurve>& opposite2,
		    double factor);


    /// Destructor.
    virtual ~CrossTanOffDist();

    /// Evaluate a point on the curve for a given parameter
    /// \param t the parameter for which to evaluate the curve.
    /// \return the evaluated point
    virtual Point eval( double t) const;

    /// Evaluate a point and a certain number of derivatives 
    /// on the curve for a given parameter.
    /// \param t the parameter for which to evaluate the curve.
    /// \param n the number of derivatives (0 or more)
    /// \retval der pointer to an array of Points where the 
    ///         result will be written.  The position will be stored
    ///         first, then the first derivative (tangent), then the
    ///         second, etc..
    ///         \b NB: For most (all) derived classes of 'EvalCurve', 
    ///         the implementation actually only supports the computation of 
    ///         one derivative, i.e. if n > 1, only one derivative will be 
    ///         computed anyway.
    virtual void eval( double t, int n, Point der[]) const; // n = order of diff

    /// Get the start parameter of the curve.
    /// \return the start parameter of the curve.
    virtual double start() const;

    /// Get the end parameter of the curve.
    /// \return  the end parameter of the curve.
    virtual double end() const;

    /// Get the dimension of the space in which the curve lies.
    /// \return the space dimension of the curve.
    virtual int dim() const;

    /// Check if the curve, evaluated at a given parameter, approximates
    /// a given position within a given tolerance.
    /// \param par the parameter at which to check the curve
    /// \param approxpos the position we want to check whether or not the curve
    ///                  approximates for parameter 'par'.
    /// \param tol1 approximation tolerance,
    /// \param tol2 another approximation tolerance (its use is defined by some of
    ///             the derived classes.
    /// \return 'true' if the curve approximates the point at the parameter, 'false'
    ///         otherwise.
    virtual bool approximationOK(double par, Point approxpos,
				 double tol1, double tol2) const;

private:
    const shared_ptr<SplineCurve> poscurve_;
    std::vector<shared_ptr<SplineCurve> > tangcurves_;
    std::vector<shared_ptr<SplineCurve> > blends_;
    shared_ptr<SplineCurve> oppositepos_;
    shared_ptr<SplineCurve> lengthfac_;
    shared_ptr<SplineCurve> avcross_;
    const double epstol_;

    /// Evaluate up to n derivatives of th offset tangent in the input parameter.
    /// \param t parameter in which to evaluate.
    /// \param derblend the derivatives in the offset direction.
    /// \param derproj the projection of the difference vector between the
    ///                derivatives in the offset direction when using the blending
    ///                functions and when using a linear blend between end offset
    ///                tangents.
    void evalcrtan(double t, Point& blend, Point& projdiff) const ;

    /// Evaluate up to n derivatives of th offset tangent in the input parameter.
    /// \param t parameter in which to evaluate.
    /// \param n the number of derivatives to compute.
    /// \param derblend array of size 'n + 1' containing the derivatives in the
    ///                 offset direction.
    /// \param derproj array of the same size, containing the projection of the
    ///                difference vector between the derivatives in the offset
    ///                direction when using the blending functions and when using
    ///                a linear blend between end offset tangents.
    void evalcrtan(double t, int n, Point derblend[], Point derproj[]) const;

    /// Evaluate the offset point using the blending functions to create a linear
    /// combination of the input tangent curves.
    /// \param t the parameter in which to evaluate.
    /// \return the blended offset point.
    Point evalblend(double t) const;

    /// Evaluate the difference between poscurve_ & oppositepos_ in input parameter.
    /// \param t parameter in which to evaluate.
    /// \return the difference vector between poscurve_ & oppositepos_ in t.
    Point evaldiff(double t) const ;

    /// Evaluate the difference between poscurve_ & oppositepos_ in input parameter,
    /// with the given number of derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n the number of derivatives to compute.
    /// \param der the difference vector between poscurve_ & oppositepos_ in t, up
    ///            to the n'th derivative. Size of vector is 'n + 1'.
    void evaldiff(double t, int n, Point der[]) const ;
};


} // namespace Go

#endif
