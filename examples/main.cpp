#include "eigenscrews/screws.hpp"
#include <iostream>

int main()
{
  using namespace eigenscrews;
  using Eigen::Vector3d;

  std::cout << "--- Screws Library Test ---\n\n";

  Vector3d q(0, 1, 0);
  Vector3d w(0, 0, 1);
  double pitch = 0.5;

  Vector6d xi = screwToTwist(pitch, q, w);
  std::cout << "1. Calculated Twist (xi) [v; w]:\n" << xi << "\n\n";

  double theta = EIGEN_PI / 2.0;
  Eigen::Isometry3d g = twistExp(xi, theta);

  std::cout << "2. Rigid Body Transformation (g) after theta=PI/2:\n";
  std::cout << "Rotation Matrix (R):\n" << g.rotation() << "\n";
  std::cout << "Translation Vector (p):\n" << g.translation().transpose() << "\n\n";

  Vector6d xi_recovered = rigidTwist(g);

  std::cout << "3. Recovered Twist from g (scaled by theta):\n" << xi_recovered << "\n\n";

  Matrix6d Ad_g = rigidAdjoint(g);
  Vector6d xi_transformed = Ad_g * xi;
  std::cout << "4. Adjoint transformed twist (Ad_g * xi):\n" << xi_transformed << std::endl;

  return 0;
}
