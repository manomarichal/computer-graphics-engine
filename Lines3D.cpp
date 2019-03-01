//============================================================================
// @name        : Lines3D.h
// @author      : Mano Marichal
// @date        :
// @version     :
// @copyright   : Computer Graphics - BA1 Informatica - Mano Marichal - University of Antwerp
// @description : Class describing lines draw in 3D
//============================================================================
#include "Lines3D.h"

Lines3D::Lines3D(ini::Configuration &conf) {

}

void Lines3D::rotateAroundX(Matrix &m, const double angle) {
    Matrix s;
    s(2,2) = std::cos(angle);
    s(1,1) = std::cos(angle);
    s(2,1) = std::sin(angle);
    s(1,2) = -1 * std::sin(angle);
    m *= s;
}

void Lines3D::rotateAroundY(Matrix &m, const double angle) {
    Matrix s;
    s(0,0) = std::cos(angle);
    s(2,2) = std::cos(angle);
    s(0,2) = std::sin(angle);
    s(2,0) = -1 * std::sin(angle);
    m *= s;
}

void Lines3D::rotateAroundZ(Matrix &m, const double angle) {
    Matrix s;
    s(0,0) = std::cos(angle);
    s(1,1) = std::cos(angle);
    s(1,0) = std::sin(angle);
    s(0,1) = -1 * std::sin(angle);
    m *= s;
}

void Lines3D::scaleMatrix(Matrix &m, const double scale) {
    Matrix s;
    s(0,0) = scale;
    s(1,1) = scale;
    s(2,2) = scale;
    m *= s;
}

void Lines3D::translateMatrix(Matrix &m, const Vector3D &v){
    Matrix s;
    s(0,3) = v.x;
    s(1,3) = v.y;
    s(2,3) = v.z;
    m *= s;

}

void Lines3D::applyTransformations(Lines3D &f, const Matrix &m) {
    for (auto p:f.getPoints()) p *= m;
}

const std::vector<Vector3D> & Lines3D::getPoints() const {
    return points;
}

void Lines3D::toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt((point.x*point.x) + (point.y*point.y) + (point.z*point.z));
    theta = std::atan2(point.y, point.x);
    phi = std::acos(r);
}


Matrix Lines3D::eyePointTrans(const Vector3D &eyepoint) {
    double theta,phi,r;
    toPolar(eyepoint,theta,phi,r);
    Vector3D v = Vector3D::vector(0,0,-r);
    Matrix m;
    translateMatrix(m,v);
    rotateAroundZ(m, -M_PI/2 - theta);
    rotateAroundX(m, phi);
    return m;
}

Point2D Lines3D::doProjection(const Vector3D &point, const double d) {
    Point2D newPoint;
    newPoint.x = (d*point.x) / (-1*point.z);
    newPoint.y = (d*point.y) / (-1*point.z);
    return newPoint;
}