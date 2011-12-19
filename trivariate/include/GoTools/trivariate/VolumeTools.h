//===========================================================================
//
// File : VolumeTools.h
//
// Created: Tue Nov 25 10:58:13 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: VolumeTools.h,v 1.2 2008-11-27 12:59:19 kfp Exp $
//
// Description:
//
//===========================================================================


#ifndef _VOLUMETOOLS_H
#define _VOLUMETOOLS_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/utils/Array.h"
#include <memory>
#include <vector>



namespace Go
{

  class SurfaceOnVolume;

    /// Analyze periodicity of volume based on number of repeating
    /// knots and control points. The return value is -1 if the volume
    /// edges are disjoint, otherwise k if sf is C^k continuous across the
    /// seam. These are sufficient but not necessary conditions for periodicity,
    /// so it is possible that a call to analyzePeriodicityDerivs() will yield a
    /// higher degree of periodicity.
    /// The current implementation is quite slow, and not optimized for speed.
    /// \param sf reference to the SplineVolume to be analyzed
    /// \param direction specify 'direction' to be '0' to check for periodicity in
    ///                  the first parameter direction, or '1' to check
    ///                  the second parameter direction, or 2 for the third.
    /// \param knot_tol the tolerance used when comparing knot intervals
    /// \return -1 if the volume edges are disjoint, otherwise k if the 
    ///         
    int GO_API
    analyzePeriodicity(const SplineVolume& sf, int direction,
                       double knot_tol = 1e-12);

  /// Describe a volume as a high-dimensional curve in a given direction.
  /// If the volume is rational, the curve will be non-rational
  /// and living in the homogenous space.
  /// \param volume the volume to express as a curve
  /// \param cv_dir the parameter direction that will be kept when defining 
  ///               the curve (the other two will disappear, as the control
  ///               points in this direction will be lumped together and expressed
  ///               as single control points in a higher-dimensional space.
  ///               'cv_dir' takes the values 0, 1 or 2 for keeping the first,
  ///               second or third parameter direction respectively
  /// \return shared pointer to a new SplineCurve, expressing the volume
  ///         as a curve in a high-dimensional space.
  shared_ptr<SplineCurve>
    representVolumeAsCurve(const SplineVolume& volume,
			   int cv_dir);

  /// Describe a curve as a lower-dimensional volume in a given direction.
  /// \param curve the curve that we want to express as a volume
  /// \param cv_dir If this variable is set to 0, then the curve's parameter
  ///               will become the \em first parameter in the generated volume.  
  ///               If it is set to 1, the curve's parameter will become the
  ///               \em second parameter in the generated volume.
  ///               If it is set to 2, the curve's parameter will become the
  ///               \em third parameter in the generated volume.  Other values
  ///               are illegal.
  /// \param other_bas1 the first BsplineBasis for the additional parameter directions.
  ///                   If cv_dir is 0, this will be the second parameter direction
  ///                   on the volume. Otherwise, it wil be the first parameter direction
  /// \param other_bas2 the second BsplineBasis for the additional parameter directions.
  ///                   If cv_dir is 2, this will be the second parameter direction
  ///                   on the volume. Otherwise, it wil be the third parameter direction
  /// \param rational define whether the generated volume shall be specified as 
  ///                 \em rational or not.
  /// \return a shared pointer to a new SplineVolume, expressing the curve
  ///         in a space of lower dimensionality.
  shared_ptr<SplineVolume>
    representCurveAsVolume(const SplineCurve& curve,
			   int cv_dir,
			   const BsplineBasis& other_bas1,
			   const BsplineBasis& other_bas2,
			   bool rational);

  /// Describe a volume as a high-dimensional surface in two given directions.
  /// If the volume is rational, the surface will be non-rational
  /// and living in the homogenous space.
  /// \param volume the volume to express as a surface
  /// \param sf_dir1 the volume parameter direction to become the first parameter
  ///                direction of the surface, either 0, 1 or 2.
  /// \param sf_dir2 the volume parameter direction to become the second parameter
  ///                direction of the surface, either 0, 1 or 2, and different from
  ///                sf_dir1. The control points in the direction different from
  ///                sf_dir1 and sf_dir2 will be lumped together and expressed as
  ///                single control points in a higher-dimensional space.
  /// \return shared pointer to a new SplineSurface, expressing the volume
  ///         as a surface in a high-dimensional space.
  shared_ptr<SplineSurface>
    representVolumeAsSurface(const SplineVolume& volume,
			     int sf_dir1,
			     int sf_dir2);

  /// Describe a surface as a lower-dimensional volume in given directions.
  /// \param surface the surface that we want to express as a volume
  /// \param sf_dir1 The volume parameter direction from the first surface parameter
  ///                direction. The value of sf_dir1 must be either 0, 1 or 2. Other values
  ///                are illegal.
  /// \param sf_dir2 The volume parameter direction from the second surface parameter
  ///                direction. The value of sf_dir2 must be either 0, 1 or 2, and different
  ///                from sf_dir1. Other values are illegal.
  /// \param other_bas the BsplineBasis for the additional volume parameter direction, different
  ///                  sf_dir1 and sfdir_2.
  /// \param rational define whether the generated volume shall be specified as 
  ///                 \em rational or not.
  /// \return a shared pointer to a new SplineVolume, expressing the surface
  ///         in a space of lower dimensionality.
  shared_ptr<SplineVolume>
    representSurfaceAsVolume(const SplineSurface& surface,
			     int sf_dir1,
			     int sf_dir2,
			     const BsplineBasis& other_bas,
			     bool rational);

  bool getVolAdjacencyInfo(shared_ptr<ParamVolume> vol1,
			   shared_ptr<SurfaceOnVolume> vol_sf1,
			   shared_ptr<ParamVolume> vol2,
			   shared_ptr<SurfaceOnVolume> vol_sf2,
			   double tol,
			   int& bd1, int& bd2, int& orientation,
			   bool& same_seq);

  bool getCorrCoefVolEnum(shared_ptr<SplineVolume> vol1,
			  shared_ptr<SplineVolume> vol2,
			  int bd1, int bd2, int orientation,
			  bool same_seq, 
			  std::vector<std::pair<int, int> >& enumeration);

  bool getVolCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			     std::vector<int>& enumeration);

  bool getVolBdCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			       int bd_cv, std::vector<int>& enumeration);

 std::vector<shared_ptr<ParamSurface> > 
    getBoundarySurfaces(shared_ptr<ParamVolume> vol);

 std::vector<shared_ptr<ParamSurface> > 
    getOrientedBoundarySurfaces(shared_ptr<ParamVolume> vol);

 shared_ptr<SurfaceOnVolume> 
   getBoundarySurface(shared_ptr<SplineVolume> vol, int idx);

 shared_ptr<SurfaceOnVolume> 
   getOrientedBoundarySurface(shared_ptr<SplineVolume> vol, int idx);

 void volCommonSplineSpace(shared_ptr<SplineVolume> vol1, int bd1,
			   shared_ptr<SplineVolume> vol2, int bd2,
			   int orientation, bool same_seq);

 // Approximate a parameter curve in the parameter space of a given volume
 // by a space curve
 shared_ptr<SplineCurve> 
   liftVolParamCurve(shared_ptr<ParamCurve> pcurve, 
		     shared_ptr<ParamVolume> vol,
		     double tol);

 // Approximate a space curve by a curve in the parameter space of a 
 // given volume
  
 shared_ptr<SplineCurve> 
   projectVolParamCurve(shared_ptr<ParamCurve> spacecurve, 
		     shared_ptr<ParamVolume> vol,
		     double tol);
} // namespace Go





#endif // _VOLUMETOOLS_H
