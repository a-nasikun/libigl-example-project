#include <igl/opengl/glfw/Viewer.h>

// Bartels
#include <linear_tet_dphi_dX.h>
#include <linear_tri2d_dphi_dX.h>
#include <linear_tri2dmesh_dphi_dX.h>
#include <linear_tri2dmesh_corotational_dq.h>
#include <linear_tetmesh_neohookean_dq2.h>
#include <assemble.h>
#include <flatten.h>
#include <flatten_multiply.h>
#include <eval_at_point.h>
#include <linear_tetmesh_dphi_dX.h>
#include <simple_psd_fix.h>

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  /*const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    0,6,4,
    0,2,6,
    0,3,2,
    0,1,3,
    2,7,6,
    2,3,7,
    4,6,7,
    4,7,5,
    0,4,5,
    0,5,1,
    1,5,7,
    1,7,3).finished();*/

    Eigen::MatrixXd V, dphidX, X;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    Eigen::VectorXd v; //volumes

    igl::readMESH("../Bartels/data/meshes_mesh/coarser_bunny.mesh", V, E, F);
    std::cout << "V" << V.cols() << "x" << V.rows() << std::endl;
    std::cout << "E" << E.cols() << "x" << E.rows() << std::endl;
    std::cout << "F" << F.cols() << "x" << F.rows() << std::endl;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, E);
  viewer.data().set_face_based(true);
  viewer.launch();
}
