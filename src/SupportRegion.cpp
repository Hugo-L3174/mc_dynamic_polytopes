#include "../include/SupportRegion.h"
#include <thread>         // std::thread
#include <chrono>
#include <eigen-quadprog/QuadProg.h>
#include <eigen-quadprog/eigen_quadprog_api.h>
#include <type_traits>
#include "WrenchCones.h"

std::vector<Eigen::Vector3d> concatenate_vectors(const std::vector<Eigen::Vector3d> & a,const std::vector<Eigen::Vector3d> & b)
{
    std::vector<Eigen::Vector3d> pts = a;
    pts.insert( pts.end(), b.begin(), b.end() );
    return pts;
}

/**
 * @brief 
 * 
 * @param QPsuccess 
 * @param pz projection of contacts through the span representation on the plane
 * @param fs span representation of contacts
 * @param v ray of LP
 * @param plane_n plane normal
 * @param target_force commonly a vertical force
 * @return const Eigen::VectorXd 
 */
const Eigen::VectorXd compute_force_limit_LP(bool & QPsuccess,
                                          const std::vector<Eigen::Vector3d> & pz,
                                          const std::vector<Eigen::Vector3d> & fs,
                                          const Eigen::Vector3d & v,
                                          const Eigen::Vector3d & plane_n,
                                          const Eigen::Vector3d & target_force)
{
  const auto t_s_ilp = std::chrono::high_resolution_clock::now();
  assert(pz.size() == fs.size());
  
  const int n_lbd = fs.size(); //nb of variables



  /**
   * Inequality cstr are that decision var are positive
   * 
   *
   * Equality cstr are
   * The resulting force is target_force
   * The moment around plane_n is 0
   */
  Eigen::Matrix4Xd Aeq(4,n_lbd); 

  Eigen::Vector4d beq = {target_force.x(),target_force.y(),target_force.z() ,0};

  /**
   * Mpz is the matrix thatt computes the ZMP w.r.t the decision variables
   * 
   */
  Eigen::Matrix3Xd Mpz(3,n_lbd);

  //In this case, the total force is known, we pre compute the projection w.r.t plane_n
  const double sigma = plane_n.dot(target_force) ;
  for (int i = 0 ; i < n_lbd ; i++)
  {

    Aeq.col(i) << fs[i].x() , fs[i].y() , fs[i].z(), plane_n.dot( pz[i].cross(fs[i]) );
    
    Mpz.col(i) = (pz[i] * plane_n.dot(fs[i]) / sigma);

  }
  
  Eigen::QuadProgDense qp;
  qp.problem(n_lbd,beq.rows(),n_lbd);

  
  const auto t_e_ilp = std::chrono::high_resolution_clock::now();
  const auto dur_i = std::chrono::duration_cast<std::chrono::microseconds>(t_e_ilp - t_s_ilp).count();
  // mc_rtc::log::info("lp init {} us",dur_i);


  const auto t_s_lp = std::chrono::high_resolution_clock::now();
  QPsuccess = qp.solve(1e-6 * Eigen::MatrixXd::Identity(n_lbd,n_lbd),
                        Mpz.transpose() * v,
                        Aeq,beq,
                        -Eigen::MatrixXd::Identity(n_lbd,n_lbd),Eigen::VectorXd::Zero(n_lbd));
  const auto t_e_lp = std::chrono::high_resolution_clock::now();
  const auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t_e_lp - t_s_lp).count();
  // mc_rtc::log::success("lp sol {} us",dur);

  const auto res = qp.result();
  
  if (!QPsuccess)
  {
    mc_rtc::log::error("[compute_force_limit]LP failed");
    return Eigen::VectorXd::Zero(n_lbd);
  }
  else
  {
    return res;
  }
}


const Eigen::Vector3d compute_cmp_cone_LP(bool & QPsuccess,
                                          const std::vector<Eigen::Vector3d> & pz,
                                          const std::vector<Eigen::Vector3d> & fs,
                                          const Eigen::Vector3d & v,
                                          const Eigen::Vector3d & plane_n,
                                          const Eigen::Vector3d & zmp)
{
  assert(pz.size() == fs.size());
  
  const int n_lbd = fs.size();
  Eigen::Matrix3Xd Aeq(3,n_lbd);
  Eigen::Matrix3Xd ray_mat(3,n_lbd);
  Eigen::MatrixXd cstr_normals(6,3);
  Eigen::VectorXd sigma_vec(n_lbd);
  cstr_normals << 1, 0, 0,
                     -1, 0, 0,
                      0, 1, 0,
                      0, -1, 0,
                      0, 0, 1,
                      0, 0, -1;
  
  Eigen::Vector6d cstr_off = Eigen::Vector6d::Ones() * 400;

  Eigen::Vector3d beq = Eigen::Vector3d::Zero();
  for (int i = 0 ; i < n_lbd ; i++)
  {
    sigma_vec(i) = plane_n.dot(fs[i]);
    Aeq.block(0,i,3,1) = sigma_vec(i) * (pz[i] - zmp);
    ray_mat.block(0,i,3,1) = fs[i];
  }

  Eigen::MatrixXd Aineq(n_lbd + cstr_normals.rows(),n_lbd);
  Aineq.block(0,0,n_lbd,n_lbd) = -Eigen::MatrixXd::Identity(n_lbd,n_lbd);
  Aineq.block(n_lbd,0,cstr_normals.rows(),n_lbd) = cstr_normals * ray_mat;
  
  Eigen::VectorXd bineq(Aineq.rows());
  bineq.segment(0,n_lbd).setZero();
  bineq.segment(n_lbd,cstr_off.rows()) = cstr_off;

  Eigen::QuadProgDense qp;
  qp.problem(n_lbd,beq.rows(),bineq.rows());

  QPsuccess = qp.solve(1e-8 * Eigen::MatrixXd::Identity(n_lbd,n_lbd),
                      ray_mat.transpose() * v ,
                      Aeq,beq,
                      Aineq,bineq);

  const auto res = qp.result();
  
  if (!QPsuccess)
  {
    mc_rtc::log::error("[compute_cmp_cone_LP] LP failed");
    return Eigen::Vector3d::Zero();
  }
  else
  {
    return ray_mat * res;
  }
}

Eigen::Matrix3d compute_plane_frame(const std::vector<Eigen::Vector3d> & hull,const Eigen::Vector3d & plane_p,const Eigen::Vector3d & plane_n)
{
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  R.col(2) << plane_n / plane_n.norm();

  size_t i = 0;
  while (i < hull.size())
  {
      auto & p1 = hull[i];
      
      if ((p1 - plane_p).norm() > 1e-4)
      {
        R.col(0) << (p1 - plane_p).cross(plane_n);
        R.col(0) /= R.col(0).norm();
        R.col(1) << (p1 - plane_p) / (p1 - plane_p).norm() ;
      }
      i+=1;
  }
  return R;
}

void project_ray( ContactProjection & c_proj,
                  const Eigen::Vector3d & cone_pt,
                  const Eigen::Vector3d & ray)
{
  const double n_dot_f = c_proj.plane_n.dot(ray);
  //If the ray is parallel to the plane, stop
  if( n_dot_f == 0)
  {
    return;
  }
  const Eigen::Vector3d proj = c_proj.plane_p + (c_proj.plane_n.cross( (cone_pt - c_proj.plane_p).cross(ray) )/n_dot_f);

  c_proj.pts.push_back(proj);
  c_proj.f.push_back(ray);

  if(n_dot_f >=0)
  {
    c_proj.pos_pts.push_back(&c_proj.pts.back());
    c_proj.pos_f.push_back(&c_proj.f.back());
  }
  else
  {
    c_proj.neg_pts.push_back(&c_proj.pts.back());
    c_proj.neg_f.push_back(&c_proj.f.back());
  }  

}

void project_contact_cone(const mc_rbdyn::Robot & r,
                          ContactProjection & c_proj,
                          const mc_rbdyn::Contact & contact)
{
  const auto & s = *(contact.r1Surface());
  std::vector<sva::PTransformd> contact_pts = s.points();
  sva::PTransformd X_0_b;

  c_proj.max_f = computeMaxNormalforce(r,contact);

  //Check is robot surface is surface 1 or 2
  if(r.hasSurface(contact.r1Surface()->name()))
  {
    X_0_b = r.bodyPosW(contact.r1Surface()->bodyName()) ;
  }
  else
  {
    const auto & s2 = *(contact.r2Surface());
    contact_pts = s2.points();
    X_0_b = r.bodyPosW(contact.r2Surface()->bodyName());
  }

  //Rotation matrix of contact frame in world frame
  const Eigen::Matrix3d R_c_0 = (contact.X_b_s() * X_0_b).rotation().transpose() ;
  c_proj.contact_n = R_c_0.col(2).normalized();

  //vertice is a PTransform from parent body frame to vertice
  for (auto & vertice : contact_pts)
  {
    
    //Express vertice in world frame
    const Eigen::Vector3d vertice_0 = (vertice * X_0_b).translation();

    std::vector<Eigen::Vector3d> force_span;
    // if(std::is_same_v<decltype(&mc_rbdyn::Contact::force_span), void (mc_rbdyn::Contact::*)()>)
    // {
    //   force_span = contact.force_span();
    // }
    // else
    // {
    //compute span rays in contact frame
    for (int i = 0 ; i < contact.nrConeGen ; i++)
    {
      const auto R = sva::RotZ(2 * 3.14 * ( ((double) i) / 4.));
      force_span.push_back( (R * Eigen::Vector3d{contact.friction()/sqrt(2),contact.friction()/sqrt(2),1}).normalized() );
    }
    // }

    for (auto & fc : force_span)
    {

      const Eigen::Vector3d f = (R_c_0 * fc);
      project_ray(c_proj,vertice_0,f);
    
    }
  }
}

const ContactProjection contact_projection(const mc_rbdyn::Robot & robot,
                                      const std::vector<mc_rbdyn::Contact> & contacts,
                                      const Eigen::Vector3d & plane_p,
                                      const Eigen::Vector3d & plane_n)
{
  ContactProjection c_proj;
  c_proj.plane_n = plane_n.normalized();
  c_proj.plane_p = plane_p;

  for (auto & c : contacts)
  {
    project_contact_cone(robot,c_proj,c);
  }

  return c_proj;
  
}

std::vector<std::vector<Eigen::Vector3d>> static_zmp_region(const mc_rbdyn::Robot & robot,
                                                              const std::vector<mc_rbdyn::Contact> & contacts,
                                                              const Eigen::Vector3d & plane_p,
                                                              const Eigen::Vector3d & plane_n)
{
  const auto t_s = std::chrono::high_resolution_clock::now();

  const Eigen::Vector3d n = plane_n.normalized();

  const ContactProjection c_proj = contact_projection(robot,contacts,plane_p,plane_n);


  std::vector<Eigen::Vector3d> pts_statics; // Hull of the static zmp region
  const int n_ray = 15; //nb of LP to be done
  pts_statics.resize(n_ray);
  const auto & pts = c_proj.pts;
  const auto & f = c_proj.f;
  auto get_lbd = [&pts_statics,&pts,&f,&n](const int i) -> bool
  {
    bool res;
    const double theta = 2 * M_PI * ( static_cast<double>(i) / static_cast<double>(n_ray) );
    const Eigen::Matrix3d R = sva::RotZ(theta).transpose();
    const Eigen::Vector3d & v = R * Eigen::Vector3d{1.,0.,0.}; 
    const auto l = compute_force_limit_LP(res,pts,f,v,n,Eigen::Vector3d{0,0,1});

    if(res)
    {
      Eigen::Vector3d p = Eigen::Vector3d::Zero();
      for(size_t j = 0 ; j < l.size(); j++)
      {
        const double p_j = n.transpose() * f[j];
        p += l(j) * p_j * pts[j];
      }

      pts_statics[i] = p;
      
      return true;
    }
  
    return false;
    
  };
  
  std::vector<std::thread> th_list;
  if(get_lbd(0))
  {
    for(int i = 1 ; i < n_ray ; i++)
    {

      th_list.push_back(std::thread(get_lbd,i));

    }
    for(auto & th: th_list){th.join();}
  }

  return {pts_statics,c_proj.pts};



}


std::vector<Eigen::Vector3d> momentless_force_cone(const mc_rbdyn::Robot & robot,
                                                              const std::vector<mc_rbdyn::Contact> & contacts,
                                                              const Eigen::Vector3d & com,
                                                              const Eigen::Vector3d & plane_n)
{
  const auto t_s = std::chrono::high_resolution_clock::now();

  const Eigen::Vector3d n = plane_n.normalized();

  const ContactProjection c_proj = contact_projection(robot,contacts,com,plane_n);

  std::vector<ContactProjection> contacts_proj;
  for (auto & c : contacts)
  {
    contacts_proj.push_back(contact_projection(robot,{c},com,plane_n));
  }
    
  std::vector<Eigen::Vector3d> rays;
  const int n_ray = 5;
  rays.resize(n_ray);
  auto get_ray = [&rays,&n,&com,&contacts_proj,&c_proj](const int i) -> bool
  {
    bool res;
    const double theta = 2 * M_PI * ( static_cast<double>(i) / static_cast<double>(n_ray) );
    const Eigen::Matrix3d R = sva::RotZ(theta).transpose();
    const Eigen::Vector3d & v = R * Eigen::Vector3d{1.,0.,0.}; 
    const auto r = compute_cmp_cone_LP(res,c_proj.pts,c_proj.f,v,n,com);

    if(res)
    {

      rays[i] = r;
      
      return true;
    }
  
    return false;
    
  };
  
  std::vector<std::thread> th_list;
  if(get_ray(0))
  {
    for(int i = 1 ; i < n_ray ; i++)
    {

      th_list.push_back(std::thread(get_ray,i));

    }
    for(auto & th: th_list){th.join();}
  }

  return rays;


}