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
#include "vector3d.h"
#include "ini_configuration.h"
#include "Lines2D.h"

class Lines3D {
private:
    struct Point3D {
        double x,y,z;
    };

    struct Face {
        std::vector<int> pointIndexes;
    };

    Color color;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

public:
    const std::vector<Vector3D>& getPoints() const;

    Lines3D(ini::Configuration &conf);

    void rotateAroundX(Matrix &m, const double angle);

    void rotateAroundY(Matrix &m, const double angle);

    void rotateAroundZ(Matrix &m, const double angle);

    void scaleMatrix(Matrix &m, const double scale);

    void translateMatrix(Matrix &m, const Vector3D &v);

    void applyTransformations(Lines3D &f, const Matrix &m);

    Matrix eyePointTrans(const Vector3D &eyepoint);

    void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

    Point2D doProjection(const Vector3D &point, const double d);
};

typedef std::forward_list<Lines3D> Figures3D;
#endif //FUCKINGWERKCG_LINES3D_H
