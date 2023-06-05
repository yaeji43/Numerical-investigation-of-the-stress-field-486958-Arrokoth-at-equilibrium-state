#include "funcs.h"

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

matrixC structure::node;
matrixI structure::ele;
matrix6d structure::matE;
matrixC structure::DPYieldStr, structure::nodalDPY;
matrixC structure::mNode;
matrixC structure::Den;
Matrix3d structure::matMOI;
int structure::NFELEM, structure::NFNODE, structure::origin;
int structure::DIM;
MatrixXd structure::MFGRV;
VectorXd structure::OM0;
VectorXd structure::OMD0;
VectorXd structure::QU0;
VectorXd structure::XV0;
double structure::SIZ;
double structure::TOT;
double structure::GM;
double structure::POISSON;
double structure::YOUNGS;
double structure::PERIOD;
double structure::FRICTION;
double structure::DENSITY;
string structure::NDFILE;
string structure::ELFILE;
string structure::DESWITCH; // Dynamic effect switch

void structure::allocDen(Eigen::MatrixXd& Den){
        Den.setConstant(DENSITY);
}

void structure::exportOutput(Eigen::SparseMatrix<double> nodalSol, Eigen::MatrixXd matU){

        ofstream fout1("THL_Output_.txt");
        fout1 << "%% +++++++++++++++++++ SIMULATION INFO +++++++++++++++++++++" << endl;
        fout1 << "%% Node    number:= " << NFNODE << endl;
        fout1 << "%% Element number:= " << NFELEM << endl;
        fout1 << "%% Spin period:= " << PERIOD  << endl;
        fout1.setf(ios::scientific);
        fout1.precision(15);

        fout1 << "%% Node Info =================================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                fout1 << setw(15) << node(i,0) << "  ,  " << node(i,1) << "  ,  " << node(i,2) << endl;
        }

        fout1 << "%% Stress Info =================================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                fout1 << setw(15) << nodalSol.coeffRef(6*i+0,0) << "  ,  " << nodalSol.coeffRef(6*i+1,0) << "  ,  " << nodalSol.coeffRef(6*i+2,0) << "  ,  " << nodalSol.coeffRef(6*i+3,0) << "  ,  " << nodalSol.coeffRef(6*i+4,0) << "  ,  " << nodalSol.coeffRef(6*i+5,0) << endl;
        }

        fout1 << "%% Displacement Info =============================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                for (int j = 0; j < 3; j++) {
                        fout1 << setw(15) << matU(i*3 + j) <<" ";
                }
                fout1 << endl;
        }

        fout1 << "%% Mass distribution =============================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                fout1 << setw(15) << mNode(i,0) << endl;
        }

        fout1.close();
}

// Read infile.txt
// Read initial parameters/ constants/ shapemodel locations
void structure::readInitials(){
        ifstream fin("infile.txt");
        if (!fin) {
                cout << "Infile.txt file is could not opened... (check the file)" << endl;
                exit(1);
        }
        else
                cout << "Infile.txt file is successfully opened." << endl;

        XV0 = VectorXd::Zero(6);
        OM0 = VectorXd::Zero(3);
        OMD0 = VectorXd::Zero(3);
        QU0 = VectorXd::Zero(4);
        string strIn;
        int ii = 0;
        while (getline(fin, strIn)) {
                ii++;
                switch (ii) {
                case 1:
                        XV0(0)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 2:
                        XV0(1)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 3:
                        XV0(2)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 4:
                        XV0(3)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 5:
                        XV0(4)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 6:
                        XV0(5)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 7:
                        OM0(0)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 8:
                        OM0(1)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 9:
                        OM0(2)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 10:
                        OMD0(0)  = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 11:
                        OMD0(1)  = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 12:
                        OMD0(2)  = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 13:
                        QU0(0)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 14:
                        QU0(1)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 15:
                        QU0(2)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 16:
                        QU0(3)   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 17:
                        SIZ      = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 18:
                        TOT      = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 19:
                        GM       = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 20:
                        YOUNGS   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 21:
                        POISSON  = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 22:
                        DENSITY  = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 23:
                        FRICTION = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 24:
                        PERIOD   = atof((char*)strIn.substr(38, strIn.length() - 38).c_str());
                case 25:
                        NDFILE   = strIn.substr(38, strIn.length() - 38);
                case 26:
                        ELFILE   = strIn.substr(38, strIn.length() - 38);
                case 27:
                        DESWITCH = strIn.substr(38, strIn.length() - 38);
                }
        }
}



void structure::readInput(){
        DIM = 3;
        // Element
        ifstream fin1(ELFILE);
        if (!fin1) {
                cout << "ERROR: ELEMENT file cannot be opened..." << endl;
                exit(1);
        }
        string stringInput; //Format of Edge
        int Rows_num = 0;
        while(getline(fin1,stringInput)) {
                ++Rows_num;
        }
        fin1.close();
        NFELEM = Rows_num;
        ele.resize(Rows_num, 4); //define the size of ele matrix
        fin1.open(ELFILE);
        int ii = 0;
        while(getline(fin1,stringInput)) {
                ele(ii,0) = atoi((char*)stringInput.substr(0, 8).c_str());
                ele(ii,1) = atoi((char*)stringInput.substr(8, 8).c_str());
                ele(ii,2) = atoi((char*)stringInput.substr(16, 8).c_str());
                ele(ii,3) = atoi((char*)stringInput.substr(24, 8).c_str());
                ii++;
        }
        fin1.close();

        // Node
        ifstream fin2(NDFILE);
        if (!fin2) {
                cout << "ERROR: Node file cannot be opened..." << endl;
                exit(1);
        }

        Rows_num = 0;
        while(getline(fin2,stringInput)) {
                ++Rows_num;
        }
        fin2.close();
        NFNODE = Rows_num;
        node.resize(Rows_num, 3); //define the size of ele matrix
        fin2.open(NDFILE);
        ii = 0;
        while(getline(fin2,stringInput)) {
                node(ii,0) = 1.e+3*atof((char*)stringInput.substr(0, 28).c_str());
                node(ii,1) = 1.e+3*atof((char*)stringInput.substr(28, 28).c_str());
                node(ii,2) = 1.e+3*atof((char*)stringInput.substr(56, 28).c_str());
                ii++;
        }
        fin2.close();

        cout << "Node Numbers= " << NFNODE << endl;
        cout << "Element Numbers= " << NFELEM << endl;
}

void structure::CheckMesh(){
        MatrixXd moiele = MatrixXd::Zero(6,1);
        MatrixXd matMOI = MatrixXd::Zero(3,3);

        // Density Allocation
        Den.resize(NFNODE,1);
        Den.setOnes();
        allocDen(Den);

        // Mass allocation
        mNode.resize(NFNODE,1);
        mNode.setZero();
        allocMass(mNode);

        //Center of Mass
        Vector3d dx;
        Vector3d lm = Vector3d::Zero(3);
        double totalMass = 0.;
        for (int i = 0; i < NFNODE; i++) {
                lm(0) += mNode(i)*node(i,0);
                lm(1) += mNode(i)*node(i,1);
                lm(2) += mNode(i)*node(i,2);
                totalMass += mNode(i);
        }

        dx = lm/totalMass;

        cout << "Original COM:= " << dx(0) << dx(1) << dx(2) << endl;

        for (int i = 0; i < NFNODE; i++) {
                node(i,0) -= dx(0);
                node(i,1) -= dx(1);
                node(i,2) -= dx(2);
        }

        lm << 0., 0., 0.;
        totalMass = 0.;
        dx << 0., 0., 0.;
        for (int i = 0; i < NFNODE; i++) {
                lm(0) += mNode(i)*node(i,0);
                lm(1) += mNode(i)*node(i,1);
                lm(2) += mNode(i)*node(i,2);
                totalMass += mNode(i);
        }

        dx = lm/totalMass;


        cout << "Renewed COM:= "  << dx(0) << dx(1) << dx(2) << endl;

        //Find origin
        origin = findOrigin();
        node(origin,0) = dx(0);
        node(origin,1) = dx(1);
        node(origin,2) = dx(2);
        cout << "Origin:= " << origin << endl;
        cout << node(origin,0) << ", " << node(origin,1) << ", " << node(origin,2) <<  endl;

        // check MOI
        for (int i = 0; i < NFNODE; i++) {
                moiele(0) += mNode(i)*(node(i,1)*node(i,1) + node(i,2)*node(i,2)); //m(y2+z2)
                moiele(1) += mNode(i)*(node(i,0)*node(i,0) + node(i,2)*node(i,2)); //m(x2+z2)
                moiele(2) += mNode(i)*(node(i,0)*node(i,0) + node(i,1)*node(i,1)); //m(x2+y2)
                moiele(3) += -mNode(i)*node(i,0)*node(i,1); //mxy
                moiele(4) += -mNode(i)*node(i,0)*node(i,2); //mxz
                moiele(5) += -mNode(i)*node(i,1)*node(i,2); //myz
        }
        matMOI << moiele(0), moiele(3), moiele(4),
                moiele(3), moiele(1), moiele(5),
                moiele(4), moiele(5), moiele(2);
        cout << "Inertia tensor:=" << endl;
        cout << matMOI << endl;

        // Rotate the mesh to have principal spin pole
        matrix3d eigvec, trans;
        EigenSolver<MatrixXd> es(matMOI);
        eigvec = es.pseudoEigenvectors();
        trans = eigvec.inverse();

        Vector3d lin;
        for (int i = 0; i < NFNODE; i++) {
                lin << node(i,0), node(i,1), node(i,2);
                node.row(i) = trans*lin;
        }
        cout << "Transforming the spin axes for the principal rotation..." << endl;

        moiele = MatrixXd::Zero(6,1);
        matMOI = MatrixXd::Zero(3,3);
        for (int i = 0; i < NFNODE; i++) {
                moiele(0) += mNode(i)*(node(i,1)*node(i,1) + node(i,2)*node(i,2)); //m(y2+z2)
                moiele(1) += mNode(i)*(node(i,0)*node(i,0) + node(i,2)*node(i,2)); //m(x2+z2)
                moiele(2) += mNode(i)*(node(i,0)*node(i,0) + node(i,1)*node(i,1)); //m(x2+y2)
                moiele(3) += -mNode(i)*node(i,0)*node(i,1); //mxy
                moiele(4) += -mNode(i)*node(i,0)*node(i,2); //mxz
                moiele(5) += -mNode(i)*node(i,1)*node(i,2); //myz
        }
        matMOI << moiele(0), moiele(3), moiele(4),
                moiele(3), moiele(1), moiele(5),
                moiele(4), moiele(5), moiele(2);
        cout << "Inertia tensor:=" << endl;
        cout << matMOI << endl;                                                                                                                                                                                                                                                                                                                                                                                                          cout << matMOI << endl;
}

void structure::calStress(Eigen::MatrixXd& mStress, Eigen::MatrixXd matU){
        matrix612d matDTDN;
        Vector12d matUj;
        matrixI order;
        order.resize(1,4);

        for (int i = 0; i < NFELEM; i++) {
                matDTDN = calmatDTDN(i);
                order = ele.row(i);
                matUj<< matU(3*order(0)+0,0), matU(3*order(0)+1,0), matU(3*order(0)+2,0),
                        matU(3*order(1)+0,0), matU(3*order(1)+1,0), matU(3*order(1)+2,0),
                        matU(3*order(2)+0,0), matU(3*order(2)+1,0), matU(3*order(2)+2,0),
                        matU(3*order(3)+0,0), matU(3*order(3)+1,0), matU(3*order(3)+2,0);
                mStress.row(i) = matE*matDTDN*matUj; //6 x 1
        }
}


double structure::DPcriterion(Eigen::MatrixXd Stress){
        Vector3d Ps;
        double alpha, i1, j2, ystar, angle;

        Ps = calPrincipalStress(Stress);

        angle = FRICTION*PI/180.;
        i1 = (Ps(0) + Ps(1) + Ps(2));
        j2 = (1./6.)*(pow(Ps(0)-Ps(1),2) + pow(Ps(1)-Ps(2),2) + pow(Ps(2)-Ps(0),2));
        alpha = 2.*sin(angle)/(sqrt(3)*(3.-sin(angle)));
        ystar = sqrt(3.)*(3.-sin(angle))*(alpha*i1 + sqrt(j2))/(6.*cos(angle));

        return ystar;
}

MatrixXd structure::calPrincipalStress(Eigen::MatrixXd Stress){
        matrix3d Stensor;
        Vector3d Pstress;

        Stensor << Stress(0), Stress(3), Stress(4),
                Stress(3), Stress(1), Stress(5),
                Stress(4), Stress(5), Stress(2);

        Pstress = Stensor.eigenvalues().col(0).real();

        return Pstress;
}

void structure::calmatE(){
        matrix6d tempMatE;

        tempMatE << (1.-POISSON), POISSON, POISSON, 0., 0., 0.,
                POISSON, (1.-POISSON), POISSON, 0., 0., 0,
                POISSON, POISSON,(1.-POISSON), 0., 0., 0,
                0., 0., 0., (1.-2.*POISSON)/2., 0., 0.,
                0., 0., 0., 0., (1.-2.*POISSON)/2., 0.,
                0., 0., 0., 0., 0., (1.-2.*POISSON)/2.;
        matE = (YOUNGS/((1.+POISSON)*(1.-2.*POISSON)))*tempMatE; //verified
}

double structure::calVolume(int en){
        Vector3d a, b, c, d, edge1, edge2, edge3;
        double n, vol;

        a << node(ele(en,0),0), node(ele(en,0),1), node(ele(en,0),2);
        b << node(ele(en,1),0), node(ele(en,1),1), node(ele(en,1),2);
        c << node(ele(en,2),0), node(ele(en,2),1), node(ele(en,2),2);
        d << node(ele(en,3),0), node(ele(en,3),1), node(ele(en,3),2);

        edge1 = a-d;
        edge2 = b-d;
        edge3 = c-d;

        n = edge1.dot(edge2.cross(edge3));
        vol = abs(n)/6.;

        return vol;
}

MatrixXd structure::calmatDTDN(int en){
        matrix612d matDTDN;
        matrix3d tempMatM, matM;
        Vector3d dx1, dx2, dx3;

        dx1 = node.row(ele(en,1)) - node.row(ele(en,0));
        dx2 = node.row(ele(en,2)) - node.row(ele(en,0));
        dx3 = node.row(ele(en,3)) - node.row(ele(en,0));

        tempMatM << dx1(0), dx2(0), dx3(0),
                dx1(1), dx2(1), dx3(1),
                dx1(2), dx2(2), dx3(2);
        matM = tempMatM.inverse();

        matDTDN << -(matM(0,0) + matM(1,0) + matM(2,0)), 0., 0., matM(0,0), 0., 0., matM(1,0), 0., 0., matM(2,0), 0., 0.,
                0., -(matM(0,1)+ matM(1,1) + matM(2,1)), 0., 0., matM(0,1), 0., 0., matM(1,1), 0., 0., matM(2,1), 0.,
                0., 0., -(matM(0,2) + matM(1,2) + matM(2,2)), 0., 0., matM(0,2), 0., 0., matM(1,2), 0., 0., matM(2,2),
                -(matM(0,1) + matM(1,1) + matM(2,1)), -(matM(0,0) + matM(1,0) + matM(2,0)), 0., matM(0,1), matM(0,0), 0., matM(1,1), matM(1,0), 0., matM(2,1), matM(2,0), 0.,
                -(matM(0,2) + matM(1,2) + matM(2,2)), 0., -(matM(0,0) + matM(1,0) + matM(2,0)), matM(0,2), 0., matM(0,0), matM(1,2), 0., matM(1,0), matM(2,2), 0., matM(2,0),
                0., -(matM(0,2) + matM(1,2) + matM(2,2)), -(matM(0,1)+ matM(1,1) + matM(2,1)), 0., matM(0,2), matM(0,1), 0., matM(1,2), matM(1,1), 0., matM(2,2), matM(2,1);

        return matDTDN;
}

MatrixXd structure::bodyForce(){
        // Gravity + centrifugal force
        MatrixXd matB(DIM*NFNODE, 1); // size = 2*NODE x 1
        double omega;
        MFGRV.resize(DIM*NFNODE,1);

        omega = 2.*PI/PERIOD/3600.;

        Vector3d forceCE, dr, xij, forceT;
        for (int i = 0; i < NFNODE; i++) {
                dr = node.row(i) - node.row(origin);
                forceCE(0) = omega*omega*dr(0); //centrifugal force
                forceCE(1) = omega*omega*dr(1);
                forceCE(2) = 0.;

                MatrixXd forceGE = MatrixXd::Zero(3,1);
                for (int j = 0; j < NFNODE; j++) {
                        if (j != i) {
                                xij = node.row(i) - node.row(j);

                                forceGE += -G*mNode(j)/pow(xij.norm(),3)*xij;
                        }
                }
                MFGRV(3*i,0) = forceGE(0); //save gravitational force for each element
                MFGRV(3*i+1,0) = forceGE(1);
                MFGRV(3*i+2,0) = forceGE(2);

                forceT = forceCE + forceGE;
                matB(3*i,0) = forceT(0);
                matB(3*i+1,0) = forceT(1);
                matB(3*i+2,0) = forceT(2);
        }
        return matB;
}

int structure::findOrigin(){
        int nodenumber, n;
        double value, len, tol;

        tol = 1.e+2;
        value = node.row(0).norm();
        n = 0;
        for (int i = 0; i < NFNODE-1; i++) {
                len = node.row(i+1).norm();
                if (len < value) {
                        value = len;
                        n = i+1;
                }
        }

        if (value <= tol) {
                nodenumber = n;
        }

        if (value > tol) {
                cout << "The origin is too far from [0, 0, 0]" << endl << "Simulation stopped !" << endl;
                exit(1);
        }

        return nodenumber;
}

void structure::allocMass(Eigen::MatrixXd& mNode){
        double dV, tV;
        matrixI order;
        order.resize(1,4);

        for (int i = 0; i < NFELEM; i++) {
                dV = calVolume(i); // voume
                order = ele.row(i);

                mNode(order(0)) += Den(order(0))*dV/4.0;
                mNode(order(1)) += Den(order(1))*dV/4.0;
                mNode(order(2)) += Den(order(2))*dV/4.0;
                mNode(order(3)) += Den(order(3))*dV/4.0;
                tV += dV;
        }
}

void structure::applyBCVec(Eigen::MatrixXd& matB){
        int rownum;

        for (int i = 0; i < DIM; i++) {
                rownum = 3*origin;
                removeRowDenseM(matB,rownum);
        }
}

void structure::applyBCMat(Eigen::SparseMatrix<double> &matM1){
        int rownum, colnum;

        for (int i = 0; i < DIM; i++) {
                rownum = 3*origin;
                removeRowSparseM(matM1,rownum);
        }

        for (int i = 0; i < DIM; i++) {
                colnum = 3*origin;
                removeColumnSparseM(matM1, colnum);
        }
}

void structure::removeRowSparseM(Eigen::SparseMatrix<double> &matrix, unsigned int rowToRemove)
{
        unsigned int numRows = matrix.rows()-1;
        unsigned int numCols = matrix.cols();
        // Create another row-major order sparse matrix with the size of (numRows,numCols)
        // Delete the specific row
        if( rowToRemove < numRows ) {
                matrm mat1(numRows,numCols);
                mat1.middleRows(0,rowToRemove) = matrix.middleRows(0,rowToRemove);
                mat1.middleRows(rowToRemove,numRows-rowToRemove) = matrix.middleRows(rowToRemove+1,numRows-rowToRemove);

                matrix = mat1;
        }
        matrix.conservativeResize(numRows,numCols);
}

void structure::removeColumnSparseM(Eigen::SparseMatrix<double> &matrix, unsigned int colToRemove)
{
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols()-1;

        if( colToRemove < numRows ) {
                matcm mat1(numRows,numCols);
                mat1.middleCols(0,colToRemove) = matrix.middleCols(0,colToRemove);
                mat1.middleCols(colToRemove,numCols-colToRemove) = matrix.middleCols(colToRemove+1,numCols-colToRemove);

                matrix = mat1;
        }
        matrix.conservativeResize(numRows,numCols);
}

void structure::removeRowDenseM(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
        unsigned int numRows = matrix.rows()-1;
        unsigned int numCols = matrix.cols();

        if( rowToRemove < numRows )
                matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

        matrix.conservativeResize(numRows,numCols);
}

void structure::removeColumnDenseM(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols()-1;

        if( colToRemove < numCols )
                matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

        matrix.conservativeResize(numRows,numCols);
}

void structure::addBC(Eigen::MatrixXd& matU, Eigen::MatrixXd tempmatU){

        for (int i = 0; i < NFNODE; i++) {
                if (i != origin) {
                        if (i < origin) {
                                matU(DIM*i + 0) = tempmatU(DIM*i + 0);
                                matU(DIM*i + 1) = tempmatU(DIM*i + 1);
                                matU(DIM*i + 2) = tempmatU(DIM*i + 2);
                        }
                        if (i > origin) {
                                matU(DIM*i + 0) = tempmatU(DIM*(i-1) + 0);
                                matU(DIM*i + 1) = tempmatU(DIM*(i-1) + 1);
                                matU(DIM*i + 2) = tempmatU(DIM*(i-1) + 2);
                        }
                }
        }
}

void structure::obtainNodal(Eigen::SparseMatrix<double> &nodalSol, Eigen::MatrixXd mStress){
        matcm eles(6*NFELEM,1);
        matcm matA(6*NFELEM,6*NFNODE);
        Matrix<int, 4, 1> order;

        for (int i = 0; i < NFELEM; i++) {
                eles.coeffRef(6*i+0,0) = mStress(i,0);
                eles.coeffRef(6*i+1,0) = mStress(i,1);
                eles.coeffRef(6*i+2,0) = mStress(i,2);
                eles.coeffRef(6*i+3,0) = mStress(i,3);
                eles.coeffRef(6*i+4,0) = mStress(i,4);
                eles.coeffRef(6*i+5,0) = mStress(i,5);
        }

        for (int en = 0; en < NFELEM; en++) {
                order = ele.row(en);

                for (int k = 0; k < 6; k++) {
                        matA.coeffRef(6*en+k,6*order(0)+k) = 0.25;
                        matA.coeffRef(6*en+k,6*order(1)+k) = 0.25;
                        matA.coeffRef(6*en+k,6*order(2)+k) = 0.25;
                        matA.coeffRef(6*en+k,6*order(3)+k) = 0.25;
                }
        }

        LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
        // lscg.setTolerance(1.e-16);
        lscg.compute(matA);
        nodalSol = lscg.solve(eles);
        std::cout << "NodalSol #iterations:     " << lscg.iterations() << "  estimated error: " << lscg.error() << std::endl;
}

void structure::defEquilbrium(){
        matrix12d tempMatB, matB, matMA, matP;
        matrix612d matDTDN;
        Vector12d matMB;
        double vol, density;
        matrixI order;
        MatrixXd matI = MatrixXd::Identity(12, 12);
        MatrixXd subX(4*DIM,1);
        MatrixXd subP(4*DIM,1);
        matcm matM1(DIM*NFNODE, DIM*NFNODE);
        matcm matM2(DIM*NFNODE, DIM*NFNODE);
        MatrixXd matBF(DIM*NFNODE,1);
        MatrixXd tempMatM2(DIM*NFNODE,1);
        MatrixXd tempmatU(DIM*(NFNODE-1),1);
        MatrixXd matU = MatrixXd::Zero(DIM*NFNODE,1);
        MatrixXd mStress = MatrixXd::Zero(NFELEM,6);

        calmatE();
        tempMatB << 2., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.,
                0., 2., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0.,
                0., 0., 2., 0., 0., 1., 0., 0., 1., 0., 0., 1.,
                1., 0., 0., 2., 0., 0., 1., 0., 0., 1., 0., 0.,
                0., 1., 0., 0., 2., 0., 0., 1., 0., 0., 1., 0.,
                0., 0., 1., 0., 0., 2., 0., 0., 1., 0., 0., 1.,
                1., 0., 0., 1., 0., 0., 2., 0., 0., 1., 0., 0.,
                0., 1., 0., 0., 1., 0., 0., 2., 0., 0., 1., 0.,
                0., 0., 1., 0., 0., 1., 0., 0., 2., 0., 0., 1.,
                1., 0., 0., 1., 0., 0., 1., 0., 0., 2., 0., 0.,
                0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 2., 0.,
                0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 2.;
        matB = (1./120.)*tempMatB; //verified

        order.resize(1,4);
        int col, row, quot, rem;
        for (int i = 0; i < NFELEM; i++) {
                // Output: matM1 and matM2 (M1*U = M2*B)
                order = ele.row(i);
                matDTDN = calmatDTDN(i); //define M1 (M1*U = M2*b)
                vol = calVolume(i); // volume
                density = (Den(order(0))+Den(order(1))+Den(order(2))+Den(order(3)))/4.;
                matMA = (vol/density)*matDTDN.transpose()*matE*matDTDN; //volume //check DEN

                matP = 6.*matB*vol*matI;

                for (int j = 0; j < 4*DIM; j++) {
                        subX = matMA.col(j);
                        subP = matP.col(j);
                        quot = j/3;         // make a function
                        rem = j%3;

                        col = 3*order(0,quot) + rem;

                        for (int k = 0; k < 4*DIM; k++) {
                                quot = k/3; // make a function
                                rem = k%3;

                                row = 3*order(0,quot) + rem;
                                matM1.coeffRef(row,col) += subX(k,0);
                                matM2.coeffRef(row,col) += subP(k,0);
                        }
                }
                if (i%1000 == 0) {
                        cout << "Calculated vertex number " << i << "..." << endl;
                }
        }
        matBF = bodyForce();         //define body force vector
        tempMatM2 = matM2*matBF;

        // Applying BC
        // applyBC(matM1, tempMatM2); //reduced the matrix matM1:matrix tempMatM2:vector
        applyBCMat(matM1);
        applyBCVec(tempMatM2);
        applyBCMat(matM2);

        //Solve matrix by using LeastSquaresConjugateGradient
        LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
        // lscg.setMaxIterations(1.e+6);
        // lscg.setTolerance(1.e-16);
        lscg.compute(matM1);
        tempmatU = lscg.solve(tempMatM2);
        std::cout << "StaticCase #iterations:     " << lscg.iterations() << "  estimated error: " << lscg.error() << std::endl;

        //should update displacement (Future work)
        addBC(matU, tempmatU);

        // Calculate stress
        calStress(mStress, matU);

        // Convert it to nodal stress
        matcm nodalSol(6*NFNODE,1);
        obtainNodal(nodalSol,mStress);

        // Drucker-Prager yield condition
        DPYieldStr.resize(NFELEM,1);
        for (int i = 0; i < NFELEM; i++) {
                DPYieldStr(i) = DPcriterion(mStress.row(i));
        }

        nodalDPY.resize(NFNODE,1);
        Vector6d stress;
        for (int i = 0; i < NFNODE; i++) {
                stress << nodalSol.coeffRef(6*i+0,0),
                        nodalSol.coeffRef(6*i+1,0),
                        nodalSol.coeffRef(6*i+2,0),
                        nodalSol.coeffRef(6*i+3,0),
                        nodalSol.coeffRef(6*i+4,0),
                        nodalSol.coeffRef(6*i+5,0);
                nodalDPY(i) = DPcriterion(stress);
        }

        exportOutput(nodalSol,matU); // Export the data
        cout << "Static case is finished ................ " << endl;

        // Dynamic effect
        // Input: Mat1, Mat2, MatU, timestep (h)
        // Note that input parameters were applied BC
        // Matrix size: (DIM*NFNODE - NBC, DIM*NFNODE - NBC);
        // Vector size: (DIM*NFNODE - NBC, 1);
        string cond1 = "ON";
        string cond2 = "OFF";
        if (cond2.compare(DESWITCH) != 0) {
                if (cond1.compare(DESWITCH) == 0) {
                        cout << "Dynamic effect is 'ON '..." << endl;
                        DynamicEffect(matM1, matM2, tempmatU, SIZ);
                }
                else{
                        cout << "!!!!!!!!!!!!!!!!!!!!!!!! Not proper input condition !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                        cout << "You can choose either of 'ON ' or 'OFF '" << endl;
                }
        }
        else
                cout << "Dynamic effect is 'OFF '..." << endl;
}

void structure::DynamicEffect(Eigen::SparseMatrix<double> &MatM1_bc, Eigen::SparseMatrix<double> &MatM2_bc, Eigen::MatrixXd& u_bc, double h){

        //Initialization
        matcm nodalSol(6*NFNODE,1);
        MatrixXd mStress(NFELEM,6);
        MatrixXd matU(DIM*NFNODE,1);
        MatrixXd u_bc_state = MatrixXd::Zero(DIM*(NFNODE-1),1);
        MatrixXd v_bc_state = MatrixXd::Zero(DIM*(NFNODE-1),1); //u_dot = 0
        VectorXd xstate = VectorXd::Zero(14);
        u_bc_state = u_bc;

        // Initialization
        xstate << XV0, OM0, QU0, 0.;
        xstate(8) = 2.*PI/PERIOD/3600.; // updated with the period

        // Outputfile settings
        // Output
        ofstream fout1("state.txt");
        fout1.setf(ios::scientific);
        fout1.precision(7);
        ofstream fout2("stress.txt");
        fout2.setf(ios::scientific);
        fout2.precision(7);
        ofstream fout3("disp.txt");
        fout3.setf(ios::scientific);
        fout3.precision(7);

        int iimax = (int) (TOT/h);
        cout << "Total iterations will be " << iimax << endl;
        for (int ii = 0; ii < iimax; ii++) {
                tie(u_bc_state, v_bc_state, xstate) = rk4(MatM1_bc, MatM2_bc, u_bc_state, v_bc_state, xstate, h);

                matU = MatrixXd::Zero(DIM*NFNODE,1);
                mStress = MatrixXd::Zero(NFELEM,6);
                addBC(matU, u_bc_state);
                calStress(mStress, matU);     // Calculate stress
                obtainNodal(nodalSol,mStress);     // Convert it to nodal stress

                // Export output data
                // Statevec
                for (int kk = 0; kk < 14; kk++)
                        fout1 << setw(15) << xstate(kk);
                fout1 << endl;

                // Nodal Stress
                for (int jj = 0; jj < NFNODE; jj++) {
                        for (int kk = 0; kk < 6; kk++) {
                                fout2 << setw(15) << nodalSol.coeffRef(6*jj + kk,0);
                        }
                        fout2 << endl;
                }

                // displacement (u)
                for (int kk = 0; kk < NFNODE; kk++) {
                        for (int jj = 0; jj < 3; jj++) {
                                fout3 << setw(15) << matU(3*kk +jj,0);
                        }
                        fout3 << endl;
                }
                cout << "iteration #" << ii << " is done ..." << endl;
        }
        fout1.close();
        fout2.close();
        fout3.close();
}

tuple<MatrixXd, MatrixXd, VectorXd> structure::rk4(Eigen::SparseMatrix<double> MatM1_bc, Eigen::SparseMatrix<double> MatM2_bc, Eigen::MatrixXd u, Eigen::MatrixXd v, Eigen::VectorXd x, double h){

        int dem = DIM*(NFNODE-1);
        MatrixXd k1u(dem,1), k2u(dem,1), k3u(dem,1), k4u(dem,1);
        MatrixXd k1v(dem,1), k2v(dem,1), k3v(dem,1), k4v(dem,1);
        VectorXd k1x(14), k2x(14), k3x(14), k4x(14); //why 6 size?

        VectorXd ua(dem,1), va(dem,1);
        VectorXd ub(dem,1), vb(dem,1);
        VectorXd xa(14), xb(14);

        tie(k1u, k1v, k1x) = force(MatM1_bc, MatM2_bc, u, v, x);
        ua = u + 0.5 * k1u * h;
        va = v + 0.5 * k1v * h;
        xa = x + 0.5 * k1x * h;

        tie(k2u, k2v, k2x) = force(MatM1_bc, MatM2_bc, ua, va, xa);
        ua = u + 0.5 * k2u * h;
        va = v + 0.5 * k2v * h;
        xa = x + 0.5 * k2x * h;

        tie(k3u, k3v, k3x) = force(MatM1_bc, MatM2_bc, ua, va, xa);
        ua = u + k3u * h;
        va = v + k3v * h;
        xa = x + k3x * h;

        tie(k4u, k4v, k4x) = force(MatM1_bc, MatM2_bc, ua, va, xa);
        ub = u + 1. / 6. * (k1u + 2. * k2u + 2. * k3u + k4u) * h;
        vb = v + 1. / 6. * (k1v + 2. * k2v + 2. * k3v + k4v) * h;
        xb = x + 1. / 6. * (k1x + 2. * k2x + 2. * k3x + k4x) * h;

        return {ub, vb, xb};
}

tuple<MatrixXd, MatrixXd, VectorXd> structure::force(Eigen::SparseMatrix<double> MatM1_bc, Eigen::SparseMatrix<double> MatM2_bc, Eigen::MatrixXd u, Eigen::MatrixXd v, Eigen::VectorXd x){
        int dem = DIM*(NFNODE-1);
        MatrixXd a(dem,1);   // Deformation acceleration
        MatrixXd f(DIM * NFNODE,1); // Force calculation
        MatrixXd f_bc(dem,1); // Acceleration withoug bcs
        VectorXd sx(6), dsx(6); // COM state
        Vector3d om, dom;  // Angular velocity
        VectorXd qu(4), dqu(4); // Euler angle
        VectorXd xa(14);   // New state of rigid body.

        // Create state vector x
        // Combine all sx, om, qu
        // COM state
        for (int ii = 0; ii < 6; ii++)
                sx(ii) = x(ii);
        // Angular velocity
        for (int ii = 0; ii < 3; ii++)
                om(ii) = x(ii + 6);
        // Euler angle
        for (int ii = 0; ii < 4; ii++)
                qu(ii) = x(ii + 9);

        Vector3d Rc; // Vector from the Earth centor to asteroid
        double rRc;  // Distance from the Earth center to asteroid
        Rc(0) = sx(0);
        Rc(1) = sx(1);
        Rc(2) = sx(2);
        rRc = Rc.norm();

        // Translation of COM
        dsx(0) = sx(3); // Vel x
        dsx(1) = sx(4); // Vel y
        dsx(2) = sx(5); // Vel z
        dsx(3) = -GM / rRc / rRc / rRc * Rc(0); // Acc x
        dsx(4) = -GM / rRc / rRc / rRc * Rc(1); // Acc y
        dsx(5) = -GM / rRc / rRc / rRc * Rc(2); // Acc z

        // Get angular velocity
        Vector3d b;
        ColPivHouseholderQR<Matrix3d> dec(matMOI);
        b   = -om.cross(matMOI * om);// Rotational condition
        dom = dec.solve(b);

        // Get Euler angle
        dqu = dEuler(qu, om);

        // Convert Euler to DCM
        Matrix3d DCM;
        DCM = Euler_to_DCM(qu);

        // Rc in the body frame
        Vector3d RcB;
        double rRcB;
        RcB  = DCM * Rc;
        rRcB = RcB.norm();

        // Computation of deformation
        // Acceleration for self-rotation and gravity
        f = get_acc(om, dom);
        applyBCVec(f);
        f_bc = f;

        // Update u_ddot from eq26 in Hirabayashi et al., 2021
        // Using fopenmp
        MatrixXd tempMatM2(dem,1);
        MatrixXd tempMata(dem,1);
        tempMatM2 = MatM1_bc*u;

        LeastSquaresConjugateGradient<SparseMatrix<double> > lscg1;
        // lscg1.setMaxIterations(1.e+6);
        // lscg1.setTolerance(1.e-16);
        lscg1.compute(MatM2_bc);
        tempMata = lscg1.solve(tempMatM2);
        std::cout << "Dynamic case #iterations:     " << lscg1.iterations() << "  estimated error: " << lscg1.error() << std::endl;

        a = -tempMata + f_bc;


        // int tidmax = omp_get_max_threads();
        // int nd     = dim / tidmax;
        // int ntid, tid, nd2;
        //
        // #pragma omp parallel private(tid, ntid, nd2)
        // {
        //         tid  = omp_get_thread_num();
        //         ntid = tid * nd;
        //         nd2  = dim - ntid;
        //         if (tid == tidmax - 1)
        //                 a.segment(ntid, nd2) = -nn_nen.block(ntid, 0, nd2, dim) * u + f_bc.segment(ntid, nd2);
        //         else
        //                 a.segment(ntid, nd) = -nn_nen.block(ntid, 0, nd, dim) * u + f_bc.segment(ntid, nd);
        // }

        // Update: X_dot
        // COM state
        for (int ii = 0; ii < 6; ii++)
                xa(ii) = dsx(ii);
        // Angular velocity
        for (int ii = 0; ii < 3; ii++)
                xa(ii + 6) = dom(ii);
        // Euler angle
        for (int ii = 0; ii < 4; ii++)
                xa(ii + 9) = dqu(ii);

        // Time
        xa(13) = 1.;


        return {v,a,xa};
}

VectorXd structure::dEuler(Eigen::VectorXd& qu,Eigen::Vector3d& om)
{
        MatrixXd omM(4, 4);
        VectorXd dEuler(4);

        omM <<     0., -om(0), -om(1), -om(2),
                om(0),     0.,  om(2), -om(1),
                om(1), -om(2),     0.,  om(0),
                om(2),  om(1), -om(0),     0.;

        dEuler = 0.5 * omM * qu;

        return dEuler;
}

Matrix3d structure::Euler_to_DCM(Eigen::VectorXd& qu)
{
        double be0, be1, be2, be3;
        Matrix3d DCM;

        // Interface
        be0 = qu(0);
        be1 = qu(1);
        be2 = qu(2);
        be3 = qu(3);

        DCM(0) = be0 * be0 + be1 * be1 - be2 * be2 - be3 * be3;
        DCM(1) = 2. * (be1 * be2 + be0 * be3);
        DCM(2) = 2. * (be1 * be3 - be0 * be2);
        DCM(3) = 2. * (be1 * be2 - be0 * be3);
        DCM(4) = be0 * be0 - be1 * be1 + be2 * be2 - be3 * be3;
        DCM(5) = 2. * (be2 * be3 + be0 * be1);
        DCM(6) = 2. * (be1 * be3 + be0 * be2);
        DCM(7) = 2. * (be2 * be3 - be0 * be1);
        DCM(8) = be0 * be0 - be1 * be1 - be2 * be2 + be3 * be3;

        return DCM;
}

MatrixXd structure::get_acc(Eigen::Vector3d& om, Eigen::Vector3d& om_dot)
{
        MatrixXd f(DIM * NFNODE,1);
        Vector3d df, r, ttorque;
        ttorque(0) = 0.;
        ttorque(1) = 0.;
        ttorque(2) = 0.;

        // Compute the acceleration acting on nodes
        // MFGRV is the gravity acceleration.
        // df is the total acceleration.
        for (int ii = 0; ii < NFNODE; ii++)
        {
                df(0) = MFGRV(3*ii,0); //gravity
                df(1) = MFGRV(3*ii+1,0);
                df(2) = MFGRV(3*ii+2,0);

                r(0)  = node(ii,0);
                r(1)  = node(ii,1);
                r(2)  = node(ii,2);

                df += -2. * om_dot.cross(r) - om.cross(om.cross(r));

                f(3 * ii + 0, 0) = df(0);
                f(3 * ii + 1, 0) = df(1);
                f(3 * ii + 2, 0) = df(2);
        }
        return f;
}
