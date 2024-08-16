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
    Eigen::MatrixXi E, Tet;
    Eigen::MatrixXi F;
    Eigen::VectorXd v; //volumes
    Eigen::VectorXi TriTag, TetTag; 

    igl::readMESH("../Bartels/data/meshes_mesh/coarse_bunny.mesh", V, Tet, F);
    //igl::readMSH("../Bartels/data/meshes_mesh/hand.mesh", V, F, Tet, TriTag, TetTag);
    //igl::readMSH("../Bartels/data/meshes_mesh/bunny.msh", V, F, Tet, TriTag, TetTag);

    // Construct the faces for visualization
    F.resize(4 * Tet.rows(), 3);
    for (int i = 0; i < Tet.rows(); ++i) 
    {
        F.row(4 * i + 0) << Tet(i, 0), Tet(i, 1), Tet(i, 2);
        F.row(4 * i + 1) << Tet(i, 0), Tet(i, 2), Tet(i, 3);
        F.row(4 * i + 2) << Tet(i, 0), Tet(i, 1), Tet(i, 3);
        F.row(4 * i + 3) << Tet(i, 2), Tet(i, 1), Tet(i, 3);
    }

    // Energy model
    double YM = 6e5; //young's modulus
    double mu = 0.4; //poissons ratio
    double D = 0.5 * (YM * mu) / ((1.0 + mu) * (1.0 - 2.0 * mu));
    double C = 0.5 * YM / (2.0 * (1.0 + mu));

    igl::volume(V, Tet, v);
    Eigen::MatrixXd Vt = V.transpose();
    Eigen::VectorXd q = Eigen::Map<Eigen::VectorXd>(Vt.data(), V.size(), 1);
    Vt.row(0) *= 2;
    Vt.row(1) *= 4;
    Vt.row(2) *= 3;

    /*
    Eigen::MatrixXd params;
    params.resize(Tet.rows(), 2);
    Eigen::SparseMatrixd H;
    params.col(0) = Eigen::VectorXd::Constant(Tet.rows(), C);
    params.col(1) = Eigen::VectorXd::Constant(Tet.rows(), D);
    sim::linear_tetmesh_neohookean_dq2(H, V, Tet, q, dphidX, v, params);


    std::cout << "V" << V.rows() << "x" << V.cols() << std::endl;
    std::cout << "Tet" << Tet.rows() << "x" << Tet.cols() << std::endl;
    std::cout << "F" << F.rows() << "x" << F.cols() << std::endl;

    std::cout << "H size: " << H.rows() << "x" << H.cols() << std::endl; 
    std::cout << H.topLeftCorner(15, 15) << std::endl; */

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
