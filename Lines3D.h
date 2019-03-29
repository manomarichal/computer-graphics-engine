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
#include <limits>
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
    Face()= default;
    Face(std::vector<int> indexes) {
        for (auto i:indexes) {
            pointIndexes.emplace_back(i);
        }
    }
    Face(int a, int b, int c) {
        pointIndexes.emplace_back(a);
        pointIndexes.emplace_back(b);
        pointIndexes.emplace_back(c);
    }

    Face(int a, int b, int c, int d) {
        pointIndexes.emplace_back(a);
        pointIndexes.emplace_back(b);
        pointIndexes.emplace_back(c);
        pointIndexes.emplace_back(d);
    }
};

struct FaceWithPoints {
    std::vector<Vector3D> pointIndexes;
};

struct stackPoint3D {
    Vector3D H, L, U, cPos;
    int index;
};

class ZBuffer: public std::vector<std::vector<double>>{
    std::vector<std::vector<double>> zVals;
public:
    ZBuffer(const int width, const int height) {

        double posInf = std::numeric_limits<double>::infinity();
        double negInf = -std::numeric_limits<double>::infinity();

        for (int x=0;x<width;x++) {
            std::vector<double> temp;
            for (int y=0;y<height;y++) {
                temp.emplace_back(posInf);
            }
            zVals.emplace_back(temp);
        }
    };

    double getZVal(int x, int y) {
        return zVals[x][y];
    }

    void setVal(int x, int y, double val) {
        zVals[x][y] = val;
    }
};


class Figure3D {
private:
    // variables used by all 3D figures
    Color color;
    std::vector<Vector3D> points;
    std::vector<Point2D> points2D;
    listWithLines lines2D;
    double rotateX, rotateY, rotateZ, scale;
    Vector3D center, eye;
    int nrOfPoints, nrOfLines;
    std::vector<Face> faces;

    //helper functions
    void createDodecahedronPoint(int a, int b, int c, std::vector<Vector3D> &tempPoints);
    void createDodecahedronFace(int a, int b, int c, int d, int e);

    // variables and functions used by a 3DLSystem
    Vector3D H, L, U;
    Vector3D currentPos;
    double currentPhi, currentTheta, currentAngle, delta;
    int recursionDepth=0; int maxRecursionDepth;
    LParser::LSystem3D system;
    std::forward_list<stackPoint3D> stackPoint;

    // helper functions for 3DLSystems
    void calculateLines(const std::string &input);


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

    // create functions
    void createCube(std::string name, const ini::Configuration &conf);

    void createLineDrawing(std::string name, const ini::Configuration &conf);

    void createTetrahedron(std::string name, const ini::Configuration &conf);

    void createOctahedron(std::string name, const ini::Configuration &conf);

    void createIsocahedron(std::string name, const ini::Configuration &conf);

    void createDodecahedron(std::string name, const ini::Configuration &conf);

    void createSphere(std::string name, const ini::Configuration &conf);

    void createCone(std::string name, const ini::Configuration &conf);

    void createCylinder(std::string name, const ini::Configuration &conf);

    void createTorus(std::string name, const ini::Configuration &conf);

    void create3DLSystem(std::string name, const ini::Configuration &conf);
};


class Wireframe {
    Color backgroundcolor;
    int imageSize, nrOfFigures;
    listWithLines lines;
    const img::EasyImage drawLines2D();
    void drawZbuffLine(ZBuffer &zBuf, img::EasyImage &img, const unsigned int x0, const unsigned int y0,
                      const double z0,
                      const unsigned int x1, const unsigned int y1,
                      const double z1,
                      const Color &color);
public:
    img::EasyImage drawWireFrame(const ini::Configuration &conf);
};


#endif //FUCKINGWERKCG_LINES3D_H
