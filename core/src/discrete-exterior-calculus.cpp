// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

#include<bits/stdc++.h>

using namespace std;
namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    double nvertex = this->mesh.nVertices();
    vector<Eigen::Triplet<double>> triplets;
    for(Vertex v: mesh.vertices()){
        double dual_area = this->barycentricDualArea(v);
        triplets.emplace_back(v.getIndex(), v.getIndex(), dual_area);
    }
    SparseMatrix<double> A;
    A.resize(nvertex, nvertex);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A; // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    double nedge = this->mesh.nEdges();
    vector<Eigen::Triplet<double>> triplets;
    for(Edge ed: mesh.edges()){
        Halfedge he1 = ed.halfedge();
        Halfedge he2 = he1.twin();
        double cot1 = this->cotan(he1);
        double cot2 = this->cotan(he2);
        double dual_area = (1.0/2.0)*(cot1+ cot2);
        triplets.emplace_back(ed.getIndex(), ed.getIndex(), dual_area);
    }
    SparseMatrix<double> A;
    A.resize(nedge, nedge);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A; // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    double nface = this->mesh.nFaces();
    vector<Eigen::Triplet<double>> triplets;
    for(Face f: mesh.faces()){
        double area = (1.0/this->faceArea(f));
        triplets.emplace_back(f.getIndex(), f.getIndex(), area);
    }
    SparseMatrix<double> A;
    A.resize(nface, nface);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    double nedge = this->mesh.nEdges();
    double nvertex = this->mesh.nVertices();
    vector<Eigen::Triplet<double>> triplets;
    for(Edge e:mesh.edges()){
        Vertex v1 = e.halfedge().tailVertex();
        Vertex v2 = e.halfedge().tipVertex();

        triplets.emplace_back(e.getIndex(), v1.getIndex(), -1.0);
        triplets.emplace_back(e.getIndex(), v2.getIndex(), 1.0);
    }
    SparseMatrix<double> A;
    A.resize(nedge, nvertex);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A; // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    double nface = this->mesh.nFaces();
    double nedge = this->mesh.nEdges();
    vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * nface);
    for(Face f1:mesh.faces()){
        for(Halfedge ed2:f1.adjacentHalfedges()){
            Edge e = ed2.edge();
            double sign = (ed2 == e.halfedge()) ? 1.0: -1.0;
            triplets.emplace_back(f1.getIndex(), e.getIndex(), sign);
        }
    }
    SparseMatrix<double> A;
    A.resize(nface, nedge);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A; // placeholder
}

} // namespace surface
} // namespace geometrycentral