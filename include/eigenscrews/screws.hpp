/*
 * C++ Port by jovaldivieso, 2025 (MIT License).
 *
 * ==========================================================================
 * ORIGINAL HEADER (MUST BE RETAINED)
 * ==========================================================================
 *
 *	      Screws.m - a screw package for mathematica
 *
 *		    Richard Murray and Sudipto Sur
 *
 *	     Division of Engineering and Applied Science
 *		  California Institute of Technology
 *			  Pasadena, CA 91125
 *
 * Screws.m is a Mathematica package for performing screw calculus.
 * It follows the treatment described in _A Mathematical Introduction to
 * Robotic Manipulation_, by R. M. Murray, Z. Li, and S. S. Sastry
 * (CRC Press, 1994).  This package implements screw theory in
 * 3-dimensional Euclidean space (some functions work in n dimensions).
 *
 * Copyright (c) 1994, Richard Murray.
 * Permission is granted to copy, modify and redistribute this file,
 * provided that this header message is retained.
 *
 * Revision history:
 *
 * Version 1.1 (7 Dec 1991)
 * Initial implementation (as RobotLinks.m)
 *
 * Version 1.2 (10 Dec 1992)
 * Removed robot specific code, leaving screw theory operations only
 * Changed name to Screws.m to reflect new emphasis
 * Rewrote some functions to use features of Mma more effectively
 *
 * Version 1.2a (30 Jan 1995)
 * Modified TwistExp, TwistMagnitude and TwistAxis, 10/6/94, Jeff Wendlandt
 * ...Changed If[w.w != 0, ... to If[(MatchQ[w,{0,0,0}] || .., ...
 * If[w.w != 0,... does not work when w is a symbolic expression.
 *
 */

#ifndef SCREWS_HPP
#define SCREWS_HPP

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace eigenscrews
{

using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;

/* AxisToSkew[w] generates skew symmetric matrix given 3 vector w */
Eigen::Matrix3d skew(const Eigen::Vector3d& w);

/* SkewToAxis[S] extracts vector w from skewsymmetric matrix S */
Eigen::Vector3d unSkew(const Eigen::Matrix3d& S);

/* SkewExp[w, theta] gives the matrix exponential of an axis w.
 * Handles both unit w (with explicit theta) and scaled w (where norm is theta).
 */
Eigen::Matrix3d skewExp(const Eigen::Vector3d& w, double theta = 1.0);

/* TwistToHomogeneous[xi] converts xi from a 6 vector to a 4X4 matrix */
Eigen::Matrix4d twistToHomogeneous(const Vector6d& xi);

/* HomogeneousToTwist[xi] converts xi from a 4x4 matrix to a 6 vector */
Vector6d homogeneousToTwist(const Eigen::Matrix4d& T);

/* TwistExp[xi, theta] gives the matrix exponential of a twist xi.
 * Handles both unit twists and scaled twists (exponential coordinates). */
Eigen::Isometry3d twistExp(const Vector6d& xi, double theta = 1.0);

/* RigidTwist[g] extracts 6 vector xi from g (Matrix Logarithm) */
Vector6d rigidTwist(const Eigen::Isometry3d& g);

/* ScrewToTwist[h, q, w] returns the twist coordinates of a screw */
Vector6d screwToTwist(double h, const Eigen::Vector3d& q, const Eigen::Vector3d& w,
                      bool infinitePitch = false);

double twistPitch(const Vector6d& xi);
double twistMagnitude(const Vector6d& xi);
std::pair<Eigen::Vector3d, Eigen::Vector3d> twistAxis(const Vector6d& xi);

Matrix6d rigidAdjoint(const Eigen::Isometry3d& g);
Eigen::Isometry3d rigidInverse(const Eigen::Isometry3d& g);

}  // namespace eigenscrews

#endif
