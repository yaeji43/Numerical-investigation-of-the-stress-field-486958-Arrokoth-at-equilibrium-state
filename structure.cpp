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
int structure::NFELEM, structure::NFNODE, structure::origin;
int structure::DIM;

void structure::allocDen(Eigen::MatrixXd& Den){
        Den.setConstant(DENSITY);
}

void structure::exportOutput(Eigen::SparseMatrix<double> nodalSol, Eigen::MatrixXd matU){

        ofstream fout1("THL_Output_Low_3.txt");
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

        fout1 << "%% Density Info =============================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                fout1 << Den(i) << endl;
        }

        fout1 << "%% Nodal Solution =============================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                fout1 << setw(15) << nodalDPY(i,0) << endl;
        }

        fout1 << "%% Mass distribution =============================" << endl;
        for (int i = 0; i < NFNODE; i++) {
                fout1 << setw(15) << mNode(i,0) << endl;
        }

        fout1.close();
}

void structure::readInput(){
        DIM = 3;
        // Element
        ifstream fin1("elem.txt");
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
        fin1.open("elem.txt");
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
        ifstream fin2("node.txt");
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
        fin2.open("node.txt");
        ii = 0;
        while(getline(fin2,stringInput)) {
                node(ii,0) = 1.e+3*atof((char*)stringInput.substr(0, 28).c_str());
                node(ii,1) = 1.e+3*atof((char*)stringInput.substr(28, 28).c_str());
                node(ii,2) = 1.e+3*atof((char*)stringInput.substr(56, 28).c_str());
                ii++;
        }
        fin2.close();
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
        cout << matMOI << endl;
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

        // for (int i = 0; i < NFNODE; i++) {
        //         cout << "M"<<i<< " := " << mNode(i) << endl;
        // }

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

void structure::applyBC(Eigen::SparseMatrix<double> &matM1, Eigen::MatrixXd& matB){
        int rownum, colnum;

        for (int i = 0; i < DIM; i++) {
                rownum = 3*origin;
                removeRowSparseM(matM1,rownum);
                removeRowDenseM(matB,rownum);
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

void structure::removeColumnSparseM(Eigen::SparseMatrix<double,ColMajor> &matrix, unsigned int colToRemove)
{
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols()-1;

        if( colToRemove < numRows )
                matrix.middleCols(colToRemove,numCols-colToRemove) = matrix.middleCols(colToRemove+1,numCols-colToRemove);

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

void structure::obtainNodal(Eigen::SparseMatrix<double> &nodalSol, Eigen::MatrixXd& mStress){
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
        lscg.setTolerance(1.e-16);
        lscg.compute(matA);
        nodalSol = lscg.solve(eles);
        std::cout << "+++++++++++++++++++ Nodal solution Info ++++++++++++++++++" << std::endl;
        std::cout << "#iterations:     " << lscg.iterations() << std::endl;
        std::cout << "estimated error: " << lscg.error()      << std::endl;
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
        applyBC(matM1, tempMatM2); //reduced the matrix

        //Solve matrix by using LeastSquaresConjugateGradient
        LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
        lscg.setMaxIterations(1.e+6);
        lscg.setTolerance(1.e-16);
        lscg.compute(matM1);
        tempmatU = lscg.solve(tempMatM2);
        std::cout << "+++++++++++++++++++ Displacement Info ++++++++++++++++++" << std::endl;
        std::cout << "#iterations:     " << lscg.iterations() << std::endl;
        std::cout << "estimated error: " << lscg.error()      << std::endl;

        //should update displacement
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
}
