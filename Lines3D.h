//============================================================================
// @name        : Lines3D.h
// @author      : Mano Marichal
// @date        :
// @version     :
// @copyright   : Computer Graphics - BA1 Informatica - Mano Marichal - University of Antwerp
// @description : Class describing lines draw in 3D
//============================================================================
#ifndef FUCKINGWERKCG_LINES3D_H
#define FUCKINGWERKCG_LINES3D_H

#include <vector>
#include <cmath>
#include <forward_list>
#include <string>
#include <iostream>
#include <map>
#include "vector3d.h"
#include "easy_image.h"
#include "ini_configuration.h"
#include "Lines2D.h"

struct Point3D {
    double x,y,z;
    void iniPoint3D(const std::vector<double> &vec) {
        x = vec[0];
        z = vec[2];
        y = vec[1];
    }

};

struct Face {
    std::vector<int> pointIndexes;
};

class Figure3D {
private:
    Color color;
    std::vector<Vector3D> points;
    std::vector<Point2D> points2D;
    listWithLines lines2D;
    double rotateX, rotateY, rotateZ, scale;
    Vector3D center, eye;
    int nrOfPoints, nrOfLines;

public:
    const std::vector<Vector3D>& getPoints() const;

    Figure3D(const std::string &name, const ini::Configuration &conf);

    void rotateAroundX(Matrix &m, const double angle);

    void rotateAroundY(Matrix &m, const double angle);

    void rotateAroundZ(Matrix &m, const double angle);

    void scaleMatrix(Matrix &m, const double scale);

    void translateMatrix(Matrix &m, const Vector3D &v);

    void applyTransformations(const Matrix &m);

    Matrix eyePointTrans(const Vector3D &eyepoint);

    void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

    void  doProjection(const Vector3D &point, const double d);

    void addLines2D(listWithLines &list);
};


class Wireframe {
    Color backgroundcolor;
    int imageSize, nrOfFigures;
    std::forward_list<Figure3D> figures;
    listWithLines lines;
    const img::EasyImage drawLines2D();
public:
    img::EasyImage drawWireFrame(const ini::Configuration &conf);
};
#endif //FUCKINGWERKCG_LINES3D_H
