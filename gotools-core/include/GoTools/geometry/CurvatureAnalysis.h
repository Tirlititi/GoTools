//===========================================================================
//                                                                           
// File: CurvatureAnalysis.h                                                 
//                                                                           
// Created: Wed Mar 16 14:33:35 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurvatureAnalysis.h,v 1.9 2008-12-15 13:42:07 kfp Exp $
//                                                                           
//===========================================================================

#ifndef _CURVATUREANALYSIS_H
#define _CURVATUREANALYSIS_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/Point.h"

/// \file CurvatureAnalysis.h
///Functions computing fundamental forms and curvature.

namespace Go
{
    /// Computes the coefficients of the first fundamental form.
    /// The returned values are stored 
    /// { E F G [Eu Fu Gu Ev Fv Gv
    ///   [Euu Fuu Guu Euv Fuv Guv Evv Fvv Gvv [...] ] ] }.
    /// Only one derivative is implemented so far.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to compute the first 
    ///          fundamental form
    /// \param v second parameter of point where we want to compute the first
    ///          fundamental form
    /// \param derivs number of (partial) derivatives of the first fundamental
    ///               form that we want to include in the result.  So far, 
    ///               only the computation of the first derivative is actually 
    ///               implemented.
    /// \param form the computed values, on the format described above
    void computeFirstFundamentalForm(const ParamSurface& sf,
				     double u, double v, int derivs,
				     std::vector<double>& form);

    /// Computes the coefficients of the first and second fundamental
    /// forms.
    /// The returned values are stored 
    /// { E F G } in form1 and { e f g } in form2.
    /// This function cannot compute derivatives.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to compute the fundamental
    ///          forms.
    /// \param v second parameter of point where we want to compute the fundamental
    ///          forms
    /// \param form1 the values associated with the first fundamental form will be 
    ///              returned here.
    /// \param form2 the values associated with the second fundamental form will be
    ///              returned here.
    void computeSecondFundamentalForm(const ParamSurface& sf,
				      double u, double v,
				      double form1[3],
				      double form2[3]);

    /// Computes the Gaussian (K) and mean (H) curvatures.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to carry out computation
    /// \param v second parameter of point where weh want to carry out computation
    /// \param K value of Gaussian curvature returned here
    /// \param H value of mean curvature returned here
    void curvatures(const ParamSurface& sf,
		    double u, double v,
		    double& K, double& H);

    void principalCurvatures(const ParamSurface& sf,
			     double u, double v,
			     double& k1, Point& d1,
			     double& k2, Point& d2);

    void minimalCurvatureRadius(const ParamSurface& sf,
				double tolerance,
				double& mincurv,
				double& pos_u,
				double& pos_v,
				double degenerate_eps,
				double curv_tol = 1.0e-3);


    void evaluateMinCurvatureRadius(const ParamSurface& sf,
				    double star_u, double end_u, double start_v, double end_v,
				    double tolerance,
				    std::vector<double>& param_u, std::vector<double>& param_v,
				    std::vector<std::vector<double> >& curvs,
				    double& mincurv,
				    double& minpos_u,
				    double& minpos_v,
				    bool initialize);




} // namespace Go

#endif // _CURVATUREANALYSIS_H
