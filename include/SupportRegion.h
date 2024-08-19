#pragma once
#include <vector>
#include <eigen3/Eigen/Dense>
#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/Contact.h>

struct ContactProjection
{
    ContactProjection(){};

    Eigen::Vector3d plane_n;
    Eigen::Vector3d plane_p;
    Eigen::Matrix3d R_frame_0;
    std::vector<Eigen::Vector3d *> pos_pts;
    std::vector<Eigen::Vector3d *> neg_pts;
    std::vector<Eigen::Vector3d *> pos_f;
    std::vector<Eigen::Vector3d *> neg_f;
    std::vector<Eigen::Vector3d> pts;
    std::vector<Eigen::Vector3d> f;
    double max_f;
    Eigen::Vector3d contact_n;

};



std::vector<std::vector<Eigen::Vector3d>> static_zmp_region(const mc_rbdyn::Robot & robot,
                                                              const std::vector<mc_rbdyn::Contact> & contacts,
                                                              const Eigen::Vector3d & plane_p,
                                                              const Eigen::Vector3d & plane_n);

std::vector<Eigen::Vector3d> momentless_force_cone(const mc_rbdyn::Robot & robot,
                                                              const std::vector<mc_rbdyn::Contact> & contacts,
                                                              const Eigen::Vector3d & com,
                                                              const Eigen::Vector3d & plane_n);