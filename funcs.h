#include <fstream>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Sparse>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstdlib>
using namespace Eigen;
using namespace std;

#define PI            3.14159265358979323846
#define G             6.67300E-11
#define POISSON       0.25
#define YOUNGS        1.0e7
#define PERIOD        15.92
#define FRICTION      35.0
#define DENSITY       500.0

class structure
{
typedef SparseMatrix<double, ColMajor> matcm;
typedef SparseMatrix<double, RowMajor> matrm;
typedef Matrix<double, Dynamic, Dynamic> matrixC;
typedef Matrix<int, Dynamic, Dynamic> matrixI;
typedef Matrix<double, 12, 1> Vector12d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 3, 3> matrix3d;
typedef Matrix<double, 6, 6> matrix6d;
typedef Matrix<double, 12, 12> matrix12d;
typedef Matrix<double, 6, 12> matrix612d;


public:
static matrixC node;
static matrixI ele;
static matrix6d matE;
static matrixC DPYieldStr, nodalDPY;
static matrixC mNode;
static matrixC Den;
static int NFELEM, NFNODE, origin;
static int DIM;


void readInput();
void calmatE();
MatrixXd calmatDTDN(int en);
double calVolume(int en);
MatrixXd bodyForce();
int findOrigin();
void allocMass(Eigen::MatrixXd& mNode);
void applyBC(Eigen::SparseMatrix<double> &matM1, Eigen::MatrixXd& matB);
void removeRowDenseM(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumnDenseM(Eigen::MatrixXd& matrix, unsigned int colToRemove);
void removeRowSparseM(Eigen::SparseMatrix<double> &matrix, unsigned int rowToRemove);
void removeColumnSparseM(Eigen::SparseMatrix<double,ColMajor> &matrix, unsigned int colToRemove);
void addBC(Eigen::MatrixXd& matU, Eigen::MatrixXd tempmatU);
void exportOutput(Eigen::SparseMatrix<double> nodalSol, Eigen::MatrixXd matU);
void calStress(Eigen::MatrixXd& mStress, Eigen::MatrixXd matU);
void defEquilbrium();
double DPcriterion(Eigen::MatrixXd Stress);
MatrixXd calPrincipalStress(Eigen::MatrixXd Stress);
void CheckMesh();
void obtainNodal(Eigen::SparseMatrix<double> &nodalSol, Eigen::MatrixXd& mStress);
void allocDen(Eigen::MatrixXd& Den);
};
