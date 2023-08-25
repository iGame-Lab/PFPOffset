//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_OSQP_H
#define THICKEN2_OSQP_H

MeshKernel::iGameVertex do_quadratic_error_metric_check(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list,bool &is_succ,double& dist){
    MeshKernel::iGameVertex v = mesh->fast_iGameVertex[vh];
    int m = neighbor_face_list.size();
    Eigen::SparseMatrix<double> hessian(3, 3);     //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian.setZero();
    Eigen::VectorXd gradient(3);                  //Q: n*1向量
    gradient.setZero();
    Eigen::SparseMatrix<double> linearMatrix(m, 3); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    linearMatrix.setZero();
    Eigen::VectorXd lowerBound(m);                  //L: m*1下限向量
    lowerBound.setZero();
    Eigen::VectorXd upperBound(m);                  //U: m*1上限向量
    upperBound.setZero();
    int cnt = 0;
    dist = 0;

    double avg_move_dist=0;
    double max_move_dist = 0;
    MeshKernel::iGameVertex avg_move_vertex(0,0,0);
    for (auto f : neighbor_face_list){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;
        max_move_dist += mesh->fast_iGameFace[f].move_dist;
    }
    avg_move_dist /= (1.0*neighbor_face_list.size());

    for (auto f : neighbor_face_list) {
        MeshKernel::iGameVertex normal
                = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();

        normal = normal * -1;
//        MeshKernel::iGameVertex new_v = v + normal * mesh->fast_iGameFace[f].move_dist;
//
//        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist*1.35;
//        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist*0.8;
        MeshKernel::iGameVertex new_v = v + normal * avg_move_dist;
        avg_move_vertex += normal * avg_move_dist;

        // MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * (1.25+0.1*depth);
        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist * max_distance_limit;
        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist;

        double d = -(normal.x() * new_v.x() + normal.y() * new_v.y() +  normal.z() * new_v.z());

        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d D2(2*d*normal.x(),2*d*normal.y(),2*d*normal.z());
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower = normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper = normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



//        for(int i=0;i<3;i++)
//            for(int j=0;j<3;j++)
//                hessian.coeffRef(i,j)+=A.coeff(i,j)*2;
//        for(int i=0;i<3;i++)
//            gradient.coeffRef(i) +=  D2.coeffRef(i);


        for(int i=0;i<3;i++)
            linearMatrix.coeffRef(cnt,i) = LimitV[i];
        lowerBound.coeffRef(cnt) = lower;
        upperBound.coeffRef(cnt) = upper;

        cnt++;
    }
    avg_move_vertex /= neighbor_face_list.size();//mesh->FastNeighborFhOfVertex_[vh].size();

    avg_move_vertex = v + avg_move_vertex;;
    hessian.coeffRef(0,0) += (2.0);
    hessian.coeffRef(1,1) += (2.0);
    hessian.coeffRef(2,2) += (2.0);

    gradient.coeffRef(0) -= (2.0) * v.x();
    gradient.coeffRef(1) -= (2.0) * v.y();
    gradient.coeffRef(2) -= (2.0) * v.z();

    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(3);   //变量数n
    solver.data()->setNumberOfConstraints(m); //约束数m
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);
    solver.initSolver();


    if(solver.solve()) {
        // if 长度过短
        is_succ = true;
        Eigen::VectorXd QPSolution = solver.getSolution();
        //dist
        MeshKernel::iGameVertex ret {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
        dist = 0;
        for(auto f : neighbor_face_list){
            MeshKernel::iGameVertex normal
                    = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();
            dist += abs((ret - v) * normal);
        }
        return ret;
    }
    else{
        is_succ = false;
        return {0,0,0};
    }
}


#endif //THICKEN2_OSQP_H
