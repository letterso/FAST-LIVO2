/*
 * math_utils.h
 *
 *  Created on: Jul 20, 2012
 *      Author: cforster
 */

#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

namespace vk
{
Eigen::Vector3d triangulateFeatureNonLin(
    const Eigen::Matrix3d& R,
    const Eigen::Vector3d& t,
    const Eigen::Vector3d& feature1,
    const Eigen::Vector3d& feature2);

/// Assumes the bearing vectors f_c and f_r are on the epipolar plane, i.e.
/// perfect triangulation without noise!
bool depthFromTriangulationExact(
    const Eigen::Matrix3d& R_r_c,
    const Eigen::Vector3d& t_r_c,
    const Eigen::Vector3d& f_r,
    const Eigen::Vector3d& f_c,
    double& depth_in_r,
    double& depth_in_c);

double reprojError(
    const Eigen::Vector3d& f1,
    const Eigen::Vector3d& f2,
    double error_multiplier2);

double computeInliers(
    const std::vector<Eigen::Vector3d>& features1,
    const std::vector<Eigen::Vector3d>& features2,
    const Eigen::Matrix3d& R,
    const Eigen::Vector3d& t,
    const double reproj_thresh,
    double error_multiplier2,
    std::vector<Eigen::Vector3d>& xyz_vec,
    std::vector<int>& inliers,
    std::vector<int>& outliers);

void computeInliersOneView(
    const std::vector<Eigen::Vector3d> & feature_sphere_vec,
    const std::vector<Eigen::Vector3d> & xyz_vec,
    const Eigen::Matrix3d &R,
    const Eigen::Vector3d &t,
    const double reproj_thresh,
    const double error_multiplier2,
    std::vector<int>& inliers,
    std::vector<int>& outliers);

//! Direct Cosine Matrix to Roll Pitch Yaw
Eigen::Vector3d dcm2rpy(const Eigen::Matrix3d &R);

//! Roll Pitch Yaw to Direct Cosine Matrix
Eigen::Matrix3d rpy2dcm(const Eigen::Vector3d &rpy);

//! Angle Axis parametrization to Quaternion
Eigen::Quaterniond angax2quat(const Eigen::Vector3d& n, const double& angle);

//! Angle Axis parametrization to Matrix representation
Eigen::Matrix3d angax2dcm(const Eigen::Vector3d& n, const double& angle);

double sampsonusError(
    const Eigen::Vector2d &v2Dash,
    const Eigen::Matrix3d& m3Essential,
    const Eigen::Vector2d& v2);

inline Eigen::Matrix3d sqew(const Eigen::Vector3d& v)
{
  Eigen::Matrix3d v_sqew;
  v_sqew << 0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0;
  return v_sqew;
}

inline double norm_max(const Eigen::VectorXd & v)
{
  double max = -1;
  for (int i=0; i<v.size(); i++)
  {
    double abs = fabs(v[i]);
    if(abs>max){
      max = abs;
    }
  }
  return max;
}

inline Eigen::Vector2d project2d(const Eigen::Vector3d& v)
{
  return v.head<2>()/v[2];
}

inline Eigen::Vector3d unproject2d(const Eigen::Vector2d& v)
{
  return Eigen::Vector3d(v[0], v[1], 1.0);
}

inline Eigen::Vector3d project3d(const Eigen::Vector4d& v)
{
  return v.head<3>()/v[3];
}

inline Eigen::Vector4d unproject3d(const Eigen::Vector3d& v)
{
  return Eigen::Vector4d(v[0], v[1], v[2], 1.0);
}

template<class T>
T getMedian(std::vector<T>& data_vec)
{
  assert(!data_vec.empty());
  typename std::vector<T>::iterator it = data_vec.begin()+floor(data_vec.size()/2);
  nth_element(data_vec.begin(), it, data_vec.end());
  return *it;
}

inline double pyrFromZero_d(double x_0, int level)
{
  return x_0/(1<<level); // = 1 / 2^level
}

inline Eigen::Vector2d pyrFromZero_2d(const Eigen::Vector2d& uv_0, int level)
{
  return Eigen::Vector2d(pyrFromZero_d(uv_0[0], level),
                  pyrFromZero_d(uv_0[1], level));
}

inline void
frameJac_xyz2uv(const Eigen::Vector3d & xyz,
                 const double & focal_length,
                 Eigen::Matrix<double,2,6> & frame_jac)
{
  const double x = xyz[0];
  const double y = xyz[1];
  const double z = xyz[2];
  const double z_2 = z*z;

  frame_jac(0,0) = -1./z *focal_length;
  frame_jac(0,1) = 0;
  frame_jac(0,2) = x/z_2 *focal_length;
  frame_jac(0,3) =  x*y/z_2 * focal_length;
  frame_jac(0,4) = -(1+(x*x/z_2)) *focal_length;
  frame_jac(0,5) = y/z *focal_length;

  frame_jac(1,0) = 0;
  frame_jac(1,1) = -1./z *focal_length;
  frame_jac(1,2) = y/z_2 *focal_length;
  frame_jac(1,3) = (1+y*y/z_2) *focal_length;
  frame_jac(1,4) = -x*y/z_2 *focal_length;
  frame_jac(1,5) = -x/z *focal_length;
}

} // end namespace vk

#endif /* MATH_UTILS_H_ */
