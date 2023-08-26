//
// Created by rainbowwing on 2023/2/17.
//

#include "remeshing.h"

#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>

using namespace vcg;
using namespace std;

class MyEdge;

class MyFace;

class MyVertex;

struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType,
        Use<MyEdge>::AsEdgeType,
        Use<MyFace>::AsFaceType> {
};

class MyVertex
        : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::BitFlags, vertex::Mark> {
};

class MyFace
        : public Face<MyUsedTypes, face::Mark, face::VertexRef, face::VFAdj, face::FFAdj, face::Normal3f, face::BitFlags> {
};

class MyEdge : public Edge<MyUsedTypes> {
};

class MyMesh : public tri::TriMesh<vector<MyVertex>, vector<MyFace>, vector<MyEdge> > {
};

void Remeshing::run(string file_name) {
    MyMesh original, toremesh;

    int loadmask;
    if (tri::io::ImporterOBJ<MyMesh>::Open(original, file_name.c_str(), loadmask) != 0) {
        printf("Error reading file  %s\n", file_name.c_str());
        exit(0);
    }
    float targetLenPerc = 0.5;
    int iterNum = 5;
    float creaseAngle = 30.f;
    float maxSurfDistPerc = 1.0;





    // Mesh cleaning
    tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);

    tri::Clean<MyMesh>::RemoveDuplicateVertex(original);
    tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
    tri::Allocator<MyMesh>::CompactEveryVector(original);
    vcg::tri::Allocator<MyMesh>::CompactEveryVector(original);


    tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(original);
    tri::UpdateBounding<MyMesh>::Box(original);

    vcg::tri::Append<MyMesh, MyMesh>::MeshCopy(toremesh, original);
    tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(toremesh);
    tri::UpdateBounding<MyMesh>::Box(toremesh);

    tri::UpdateTopology<MyMesh>::FaceFace(toremesh);
    float lengthThr = targetLenPerc * (original.bbox.Diag() / 100.f);
    float maxSurfDist = maxSurfDistPerc * (original.bbox.Diag() / 100.f);
    printf("Length Thr: %8.3f ~ %4.2f %% on %5.3f\n", lengthThr, targetLenPerc, original.bbox.Diag());

    vcg::tri::IsotropicRemeshing<MyMesh>::Params params;
    params.SetTargetLen(lengthThr);
    params.SetFeatureAngleDeg(creaseAngle);
    params.iter = iterNum;

    if (maxSurfDistPerc != 0) {
        params.surfDistCheck = false;
        params.maxSurfDist = maxSurfDist;
    } else {
        params.surfDistCheck = false;
    }
    params.splitFlag = true;
    params.collapseFlag = true;
    params.swapFlag = true;
    params.smoothFlag = true;
    params.projectFlag = true;

    params.cleanFlag = true;
    params.userSelectedCreases = false;
    params.surfDistCheck = false;


    printf(" Input mesh %8i v %8i f\n", toremesh.VN(), toremesh.FN());
    vcg::tri::IsotropicRemeshing<MyMesh>::Do(toremesh, original, params);
    cout << "remeshing result " << (file_name).c_str() << endl;
    //vcg::tri::io::ExporterOBJ<MyMesh>::Save(toremesh, (file_name).c_str(),loadmask);
    vcg::tri::io::ExporterOBJ<MyMesh>::Save(toremesh, (file_name).c_str(), loadmask);
    printf("Output mesh %8i v %8i f\n", toremesh.VN(), toremesh.FN());

    return;
}