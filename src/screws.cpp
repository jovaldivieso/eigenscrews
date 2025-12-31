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

#include "eigenscrews/screws.hpp"
#include <cmath>

namespace eigenscrews
{

/* AxisToSkew[w] generates skew symmetric matrix given 3 vector w */
Eigen::Matrix3d skew(const Eigen::Vector3d& w)
{
  Eigen::Matrix3d S;
  S << 0.0, -w.z(), w.y(), w.z(), 0.0, -w.x(), -w.y(), w.x(), 0.0;
  return S;
}

/* SkewToAxis[S] extracts vector w from skewsymmetric matrix S */
Eigen::Vector3d unSkew(const Eigen::Matrix3d& S)
{
  return Eigen::Vector3d(S(2, 1), S(0, 2), S(1, 0));
}

/* SkewExp[w,(theta)] gives the matrix exponential of an axis w. */
Eigen::Matrix3d skewExp(const Eigen::Vector3d& w, double theta)
{
  // Robustly handle both unit vectors and scaled exponential coordinates
  double mag = w.norm();
  if (mag < 1e-9)
  {
    return Eigen::Matrix3d::Identity();
  }

  // If w is unit, angle is theta. If w is scaled, angle is theta * mag.
  double angle = mag * theta;
  Eigen::Matrix3d S = skew(w / mag);  // Use normalized axis for Rodrigues formula

  return Eigen::Matrix3d::Identity() + std::sin(angle) * S + (1.0 - std::cos(angle)) * S * S;
}

/* TwistToHomogeneous[xi] converts xi from a 6 vector to a 4X4 matrix */
Eigen::Matrix4d twistToHomogeneous(const Vector6d& xi)
{
  Eigen::Vector3d v = xi.head<3>();
  Eigen::Vector3d w = xi.tail<3>();

  Eigen::Matrix4d mat = Eigen::Matrix4d::Zero();
  mat.block<3, 3>(0, 0) = skew(w);
  mat.block<3, 1>(0, 3) = v;
  return mat;
}

/* HomogeneousToTwist[xi] converts xi from a 4x4 matrix to a 6 vector */
Vector6d homogeneousToTwist(const Eigen::Matrix4d& T)
{
  Vector6d xi;
  Eigen::Vector3d v = T.block<3, 1>(0, 3);
  Eigen::Vector3d w = unSkew(T.block<3, 3>(0, 0));
  xi << v, w;
  return xi;
}

/* TwistExp[xi,(Theta)] gives the matrix exponential of a twist xi. */
Eigen::Isometry3d twistExp(const Vector6d& xi, double theta)
{
  Eigen::Vector3d v_in = xi.head<3>();
  Eigen::Vector3d w_in = xi.tail<3>();
  Eigen::Isometry3d T = Eigen::Isometry3d::Identity();

  double w_norm = w_in.norm();

  if (w_norm < 1e-9)
  {
    // Pure translation
    T.linear() = Eigen::Matrix3d::Identity();
    T.translation() = v_in * theta;
  }
  else
  {
    // Rotation + Translation
    // Normalize to get axis and screw parameters
    double angle = w_norm * theta;
    Eigen::Vector3d w_unit = w_in / w_norm;

    // The linear component v in the formula corresponds to the screw parameter,
    // which is the input v scaled down by the same factor as w.
    Eigen::Vector3d v_param = v_in / w_norm;

    Eigen::Matrix3d R = skewExp(w_unit, angle);
    Eigen::Matrix3d S = skew(w_unit);

    // Formula: (I - e^[w]theta) * (w x v) + w * (w^T v) * theta
    // Using normalized w_unit and v_param ensuring units are correct.
    Eigen::Vector3d p =
        (Eigen::Matrix3d::Identity() - R) * (S * v_param) + w_unit * (w_unit.dot(v_param)) * angle;

    T.linear() = R;
    T.translation() = p;
  }
  return T;
}

/* RigidTwist[g] extracts 6 vector xi from g */
Vector6d rigidTwist(const Eigen::Isometry3d& g)
{
  Eigen::Matrix3d R = g.linear();
  Eigen::Vector3d p = g.translation();

  Eigen::AngleAxisd angleAxis(R);
  double theta = angleAxis.angle();
  Eigen::Vector3d w = angleAxis.axis();
  Eigen::Vector3d v;

  if (std::abs(theta) < 1e-6)
  {
    double mag = p.norm();
    if (mag < 1e-9)
      return Vector6d::Zero();

    theta = mag;
    w = Eigen::Vector3d::Zero();
    v = p;  // Pure translation
  }
  else
  {
    Eigen::Matrix3d S = skew(w);
    Eigen::Matrix3d wwT = w * w.transpose();

    // Solve (I - e^[w]theta) * (w x v) + ... = p for v
    Eigen::Matrix3d A = (Eigen::Matrix3d::Identity() - wwT) * std::sin(theta) +
                        S * (1.0 - std::cos(theta)) + wwT * theta;

    v = A.colPivHouseholderQr().solve(p);

    // Scale results by theta to return exponential coordinates
    v = v * theta;
    w = w * theta;
  }

  Vector6d xi;
  xi << v, w;
  return xi;
}

/* ScrewToTwist[h, q, w] returns the twist coordinates of a screw */
Vector6d screwToTwist(double h, const Eigen::Vector3d& q, const Eigen::Vector3d& w,
                      bool infinitePitch)
{
  Vector6d xi;
  if (infinitePitch)
  {
    xi << w, Eigen::Vector3d::Zero();
  }
  else
  {
    Eigen::Vector3d v = -skew(w) * q + h * w;
    xi << v, w;
  }
  return xi;
}

/* RigidAdjoint[g] gives the adjoint matrix corresponding to g */
Matrix6d rigidAdjoint(const Eigen::Isometry3d& g)
{
  Eigen::Matrix3d R = g.linear();
  Eigen::Vector3d p = g.translation();
  Eigen::Matrix3d pR = skew(p) * R;

  Matrix6d adj;
  adj.block<3, 3>(0, 0) = R;
  adj.block<3, 3>(0, 3) = pR;
  adj.block<3, 3>(3, 0) = Eigen::Matrix3d::Zero();
  adj.block<3, 3>(3, 3) = R;

  return adj;
}

/* RigidInverse[g] gives the inverse transformation of g */
Eigen::Isometry3d rigidInverse(const Eigen::Isometry3d& g)
{
  return g.inverse();
}

/* TwistPitch[xi] gives pitch of screw corresponding to a twist */
double twistPitch(const Vector6d& xi)
{
  Eigen::Vector3d v = xi.head<3>();
  Eigen::Vector3d w = xi.tail<3>();
  double w_sq = w.dot(w);
  if (w_sq < 1e-9)
    return 0.0;
  return w.dot(v) / w_sq;
}

/* TwistMagnitude[xi] gives magnitude of screw corresponding to a twist */
double twistMagnitude(const Vector6d& xi)
{
  Eigen::Vector3d v = xi.head<3>();
  Eigen::Vector3d w = xi.tail<3>();
  if (w.isZero(1e-9))
    return v.norm();
  return w.norm();
}

/* TwistAxis[xi] gives axis of screw corresponding to a twist */
std::pair<Eigen::Vector3d, Eigen::Vector3d> twistAxis(const Vector6d& xi)
{
  Eigen::Vector3d v = xi.head<3>();
  Eigen::Vector3d w = xi.tail<3>();

  if (w.isZero(1e-9))
  {
    return { Eigen::Vector3d::Zero(), v.normalized() };
  }
  else
  {
    double w_sq = w.dot(w);
    Eigen::Vector3d axis_dir = w / w_sq;
    Eigen::Vector3d point = skew(w) * v / w_sq;
    return { point, axis_dir };
  }
}

}  // namespace eigenscrews
