//============================================================================
// @name        : Lines3D.h
// @author      : Mano Marichal
// @date        :
// @version     :
// @copyright   : Computer Graphics - BA1 Informatica - Mano Marichal - University of Antwerp
// @description : Class describing lines draw in 3D
//============================================================================
#include "Lines3D.h"

// transformation functions
void Figure3D::rotateAroundX(Matrix &m, const double angle) {
    Matrix s;
    s(3, 3) = std::cos(angle);
    s(2, 2) = std::cos(angle);
    s(3, 2) = -std::sin(angle);
    s(2, 3) = std::sin(angle);
    m *= s;
}

void Figure3D::rotateAroundY(Matrix &m, const double angle) {
    Matrix s;
    s(1, 1) = std::cos(angle);
    s(3, 3) = std::cos(angle);
    s(1, 3) = -std::sin(angle);
    s(3, 1) = std::sin(angle);
    //std::cout << "rotation aroud Y s:" << std::endl; s.print(std::cout);
    m *= s;
}

void Figure3D::rotateAroundZ(Matrix &m, const double angle) {
    Matrix s;
    //std::cout << angle << std::endl;
    s(1, 1) = std::cos(angle);
    s(2, 2) = std::cos(angle);
    s(2, 1) = -std::sin(angle);
    s(1, 2) = std::sin(angle);
    //std::cout << "rotation aroud Z s:" << std::endl; s.print(std::cout);
    m *= s;
}

void Figure3D::scaleMatrix(Matrix &m, const double scale) {
    Matrix s;
    s(1, 1) = scale;
    s(2, 2) = scale;
    s(3, 3) = scale;
    m *= s;
}

void Figure3D::translateMatrix(Matrix &m, const Vector3D &v) {
    Matrix s;
    s(4, 1) = v.x;
    s(4, 2) = v.y;
    s(4, 3) = v.z;
    m *= s;

}

void Figure3D::applyTransformations(const Matrix &m) {
    for (auto &p:points) {
        p *= m;
    }
}

const std::vector <Vector3D> &Figure3D::getPoints() const {
    return points;
}

void Figure3D::toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt((point.x * point.x) + (point.y * point.y) + (point.z * point.z));
    theta = std::atan2(point.y, point.x);
    phi = std::acos(point.z / r);
}

Matrix Figure3D::eyePointTrans(const Vector3D &eyepoint) {
    double theta, phi, r;
    toPolar(eyepoint, theta, phi, r);
    Vector3D v = Vector3D::vector(0, 0, -r);
    Matrix m;

    rotateAroundZ(m, (-M_PI / 2) - theta);
    rotateAroundX(m, -phi);
    translateMatrix(m, v);

    /*
   m(1,1) = -std::sin(theta);
   m(1,2) = -std::cos(theta) * std::cos(phi);
   m(1,3) = std::cos(theta) * std::sin(phi);
   m(2,1) = cos(theta);
   m(2,2) = -std::sin(theta) * cos(phi);
   m(2,3) = std::sin(theta) * std::sin(phi);
   m(3,2) = std::sin(phi);
   m(3,3) = std::cos(phi);
   m(4,3) = -r;
    */
    return m;
}

void Figure3D::doProjection(const Vector3D &point, const double d) {
    //std::cout << "old point: " << point.x << "|" << point.y << "|" << point.z << std::endl;
    Point2D newPoint;
    newPoint.x = (d * point.x) / (-point.z);
    newPoint.y = (d * point.y) / (-point.z);
    //std::cout << "newpoint:" << newPoint.x << "|" << newPoint.y << std::endl;
    points2D.emplace_back(newPoint);
}

// helper functions
void Figure3D::addLines2D(listWithLines &list) {
    for (const Line2D &line:lines2D) {
        list.emplace_back(line);
    }
}

void Figure3D::createDodecahedronPoint(int a, int b, int c, std::vector <Vector3D> &tempPoints) {
    points.emplace_back(Vector3D::point((tempPoints[a - 1].x + tempPoints[b - 1].x + tempPoints[c - 1].x) / 3,
                                        (tempPoints[a - 1].y + tempPoints[b - 1].y + tempPoints[c - 1].y) / 3,
                                        (tempPoints[a - 1].z + tempPoints[b - 1].z + tempPoints[c - 1].z) / 3));
}

void Figure3D::createDodecahedronFace(int a, int b, int c, int d, int e) {
    Face temp;
    temp.pointIndexes.emplace_back(a - 1);
    temp.pointIndexes.emplace_back(b - 1);
    temp.pointIndexes.emplace_back(c - 1);
    temp.pointIndexes.emplace_back(d - 1);
    temp.pointIndexes.emplace_back(e - 1);
    faces.emplace_back(temp);
}

// create figure functions
void Figure3D::createLineDrawing(std::string name, const ini::Configuration &conf) {

    nrOfPoints = conf[name]["nrPoints"].as_int_or_die();
    // read in points
    for (int k = 0; k < nrOfPoints; k++) {
        Vector3D temp = Vector3D::point(conf[name]["point" + std::to_string(k)].as_double_tuple_or_die()[0],
                                        conf[name]["point" + std::to_string(k)].as_double_tuple_or_die()[1],
                                        conf[name]["point" + std::to_string(k)].as_double_tuple_or_die()[2]);
        //temp.print(std::cout);
        points.emplace_back(temp);
        //std::cout << k << std::endl;
    }

    nrOfLines = conf[name]["nrLines"].as_int_or_die();
    // read in faces
    for (int k = 0; k < nrOfLines; k++) {
        Face temp;
        temp.pointIndexes.emplace_back(conf[name]["line" + std::to_string(k)].as_int_tuple_or_die()[0]);
        temp.pointIndexes.emplace_back(conf[name]["line" + std::to_string(k)].as_int_tuple_or_die()[1]);
        faces.emplace_back(temp);
    }

}

void Figure3D::createCube(std::string name, const ini::Configuration &conf) {

    // read in points
    points.emplace_back(Vector3D::point(0, 0, 0)); // dummy
    points.emplace_back(Vector3D::point(1, -1, -1));
    points.emplace_back(Vector3D::point(-1, 1, -1));
    points.emplace_back(Vector3D::point(1, 1, 1));
    points.emplace_back(Vector3D::point(-1, -1, 1));
    points.emplace_back(Vector3D::point(1, 1, -1));
    points.emplace_back(Vector3D::point(-1, -1, -1));
    points.emplace_back(Vector3D::point(1, -1, 1));
    points.emplace_back(Vector3D::point(-1, 1, 1));

    // read in faces
    faces.emplace_back(Face(1,5,3,7));
    faces.emplace_back(Face(5,2,8,3));
    faces.emplace_back(Face(2,6,4,8));
    faces.emplace_back(Face(6,1,7,4));
    faces.emplace_back(Face(7,3,8,4));
    faces.emplace_back(Face(1,6,2,5));
}

void Figure3D::createTetrahedron(std::string name, const ini::Configuration &conf) {
    // read in points
    points.emplace_back(Vector3D::point(0, 0, 0)); // dummy
    points.emplace_back(Vector3D::point(1, -1, -1));
    points.emplace_back(Vector3D::point(-1, 1, -1));
    points.emplace_back(Vector3D::point(1, 1, 1));
    points.emplace_back(Vector3D::point(-1, -1, 1));


    // read in faces
    faces.emplace_back(Face(1,2,3));
    faces.emplace_back(Face(2,4,3));
    faces.emplace_back(Face(1,4,2));
    faces.emplace_back(Face(1,3,4));
}

void Figure3D::createOctahedron(std::string name, const ini::Configuration &conf) {

    // read in points
    points.emplace_back(Vector3D::point(0, 0, 0)); // dummy
    points.emplace_back(Vector3D::point(1, 0, 0));
    points.emplace_back(Vector3D::point(0, 1, 0));
    points.emplace_back(Vector3D::point(-1, 0, 0));
    points.emplace_back(Vector3D::point(0, -1, 0));
    points.emplace_back(Vector3D::point(0, 0, -1));
    points.emplace_back(Vector3D::point(0, 0, 1));

    // read in faces
    Face temp;
    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(6);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(6);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(6);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(6);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(5);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(5);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(5);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(5);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

}

void Figure3D::createIsocahedron(std::string name, const ini::Configuration &conf) {
    // read in points
    points.emplace_back(Vector3D::point(0, 0, sqrt(5) / 2));
    for (int i = 2; i <= 6; i++) {
        double x = std::cos(((i - 2) * 2 * M_PI) / 5);
        double y = std::sin(((i - 2) * 2 * M_PI) / 5);
        points.emplace_back(Vector3D::point(x, y, 0.5));
    }
    for (int i = 7; i <= 11; i++) {
        double x = std::cos(M_PI / 5 + ((i - 7) * 2 * M_PI) / 5);
        double y = std::sin(M_PI / 5 + ((i - 7) * 2 * M_PI) / 5);
        points.emplace_back(Vector3D::point(x, y, -0.5));
    }
    points.emplace_back(Vector3D::point(0, 0, -sqrt(5) / 2));

    // read in faces
    Face temp;
    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(3);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(4);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(5);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(5);
    temp.pointIndexes.emplace_back(6);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(1);
    temp.pointIndexes.emplace_back(6);
    temp.pointIndexes.emplace_back(2);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(7);
    temp.pointIndexes.emplace_back(3);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(7);
    temp.pointIndexes.emplace_back(8);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(3);
    temp.pointIndexes.emplace_back(8);
    temp.pointIndexes.emplace_back(4);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();


    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(8);
    temp.pointIndexes.emplace_back(9);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(4);
    temp.pointIndexes.emplace_back(9);
    temp.pointIndexes.emplace_back(5);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(5);
    temp.pointIndexes.emplace_back(9);
    temp.pointIndexes.emplace_back(10);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(5);
    temp.pointIndexes.emplace_back(10);
    temp.pointIndexes.emplace_back(6);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(6);
    temp.pointIndexes.emplace_back(10);
    temp.pointIndexes.emplace_back(11);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(6);
    temp.pointIndexes.emplace_back(11);
    temp.pointIndexes.emplace_back(2);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(2);
    temp.pointIndexes.emplace_back(11);
    temp.pointIndexes.emplace_back(7);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(12);
    temp.pointIndexes.emplace_back(8);
    temp.pointIndexes.emplace_back(7);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(12);
    temp.pointIndexes.emplace_back(9);
    temp.pointIndexes.emplace_back(8);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(12);
    temp.pointIndexes.emplace_back(10);
    temp.pointIndexes.emplace_back(9);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(12);
    temp.pointIndexes.emplace_back(11);
    temp.pointIndexes.emplace_back(10);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    temp.pointIndexes.emplace_back(12);
    temp.pointIndexes.emplace_back(7);
    temp.pointIndexes.emplace_back(11);
    faces.emplace_back(temp);
    temp.pointIndexes.clear();

    for (auto &face:faces) {
        for (auto &i:face.pointIndexes) {
            i = i - 1;
        }
    }

}

void Figure3D::createDodecahedron(std::string name, const ini::Configuration &conf) {

    // read in points
    // generate isocahedron points
    std::vector <Vector3D> tempPoints;
    tempPoints.emplace_back(Vector3D::point(0, 0, sqrt(5) / 2));
    for (int i = 2; i <= 6; i++) {
        double x = std::cos(((i - 2) * 2 * M_PI) / 5);
        double y = std::sin(((i - 2) * 2 * M_PI) / 5);
        tempPoints.emplace_back(Vector3D::point(x, y, 0.5));
    }
    for (int i = 7; i <= 11; i++) {
        double x = std::cos(M_PI / 5 + ((i - 7) * 2 * M_PI) / 5);
        double y = std::sin(M_PI / 5 + ((i - 7) * 2 * M_PI) / 5);
        tempPoints.emplace_back(Vector3D::point(x, y, -0.5));
    }
    tempPoints.emplace_back(Vector3D::point(0, 0, -sqrt(5) / 2));

    // generate dodecahedron points
    createDodecahedronPoint(1, 2, 3, tempPoints);
    createDodecahedronPoint(1, 3, 4, tempPoints);
    createDodecahedronPoint(1, 4, 5, tempPoints);
    createDodecahedronPoint(1, 5, 6, tempPoints);
    createDodecahedronPoint(1, 6, 2, tempPoints);
    createDodecahedronPoint(2, 7, 3, tempPoints);
    createDodecahedronPoint(3, 7, 8, tempPoints);
    createDodecahedronPoint(3, 8, 4, tempPoints);
    createDodecahedronPoint(4, 8, 9, tempPoints);
    createDodecahedronPoint(4, 9, 5, tempPoints);
    createDodecahedronPoint(5, 9, 10, tempPoints);
    createDodecahedronPoint(5, 10, 6, tempPoints);
    createDodecahedronPoint(6, 10, 11, tempPoints);
    createDodecahedronPoint(6, 11, 2, tempPoints);
    createDodecahedronPoint(2, 11, 7, tempPoints);
    createDodecahedronPoint(12, 8, 7, tempPoints);
    createDodecahedronPoint(12, 9, 8, tempPoints);
    createDodecahedronPoint(12, 10, 9, tempPoints);
    createDodecahedronPoint(12, 11, 10, tempPoints);
    createDodecahedronPoint(12, 7, 11, tempPoints);

    // generate faces
    createDodecahedronFace(1, 2, 3, 4, 5);
    createDodecahedronFace(1, 6, 7, 8, 2);
    createDodecahedronFace(2, 8, 9, 10, 3);
    createDodecahedronFace(3, 10, 11, 12, 4);
    createDodecahedronFace(4, 12, 13, 14, 5);
    createDodecahedronFace(5, 14, 15, 6, 1);
    createDodecahedronFace(20, 19, 18, 17, 16);
    createDodecahedronFace(20, 15, 14, 13, 19);
    createDodecahedronFace(19, 13, 12, 11, 18);
    createDodecahedronFace(18, 11, 10, 9, 17);
    createDodecahedronFace(17, 9, 8, 7, 16);
    createDodecahedronFace(16, 7, 6, 15, 20);
}

void Figure3D::createSphere(std::string name, const ini::Configuration &conf) {
    int n = conf[name]["n"].as_int_or_die();
    createIsocahedron(name, conf);

    for (int m = 0; m < n; m++) {

        std::vector <Face> temp;

        for (int i = 0; i < faces.size(); i++) {

            int a, b, c, d, e, f;

            a = faces[i].pointIndexes[0];
            b = faces[i].pointIndexes[1];
            c = faces[i].pointIndexes[2];
            points.emplace_back(Vector3D::point((points[a].x + points[b].x) / 2,
                                                (points[a].y + points[b].y) / 2,
                                                (points[a].z + points[b].z) / 2));
            d = points.size() - 1;
            points.emplace_back(Vector3D::point((points[b].x + points[c].x) / 2,
                                                (points[b].y + points[c].y) / 2,
                                                (points[b].z + points[c].z) / 2));
            f = points.size() - 1;
            points.emplace_back(Vector3D::point((points[a].x + points[c].x) / 2,
                                                (points[a].y + points[c].y) / 2,
                                                (points[a].z + points[c].z) / 2));
            e = points.size() - 1;

            temp.emplace_back(Face(a, d, e));
            temp.emplace_back(Face(b, f, d));
            temp.emplace_back(Face(c, e, f));
            temp.emplace_back(Face(d, f, e));
        }

        faces = temp;
    }

    for (auto &point:points) {
        point.normalise();
    }

}

void Figure3D::createCone(std::string name, const ini::Configuration &conf) {
    double h = conf[name]["height"].as_double_or_die();
    int n = conf[name]["n"].as_int_or_die();

    for (int i = 0; i < n; i++) {
        points.emplace_back(Vector3D::point(std::cos((2 * i * M_PI) / n),
                                            std::sin((2 * i * M_PI) / n), 0));
    }

    points.emplace_back(Vector3D::point(0, 0, h));

    for (int i = 0; i < n; i++) {
        faces.emplace_back(Face(i, (i + 1) % n, n));
    }

    Face temp;
    for (int i = n - 1; i >= 0; i--) temp.pointIndexes.emplace_back(i);
    faces.emplace_back(temp);
}

void Figure3D::createCylinder(std::string name, const ini::Configuration &conf) {
    double h = conf[name]["height"].as_double_or_die();
    int n = conf[name]["n"].as_int_or_die();

    for (int i = 0; i < n; i++) {
        points.emplace_back(Vector3D::point(std::cos((2 * i * M_PI) / n),
                                            std::sin((2 * i * M_PI) / n), 0));
    }

    for (int i = 0; i < n; i++) {
        points.emplace_back(Vector3D::point(std::cos((2 * i * M_PI) / n),
                                            std::sin((2 * i * M_PI) / n), h));
    }

    for (int i = 0; i < n; i++) {
        faces.emplace_back(Face(i, (i + 1) % n, n+((i+1)%n), n+i));
    }

    Face temp;
    for (int i = n - 1; i >= 0; i--) temp.pointIndexes.emplace_back(i);
    faces.emplace_back(temp);

    temp.pointIndexes.clear();
    for (int i = 2*n - 1; i >= n; i--) temp.pointIndexes.emplace_back(i);
    faces.emplace_back(temp);
}

void Figure3D::createTorus(std::string name, const ini::Configuration &conf) {
    double r = conf[name]["r"].as_double_or_die();
    double R = conf[name]["R"].as_double_or_die();
    int n = conf[name]["n"].as_int_or_die();
    int m = conf[name]["m"].as_int_or_die();

    for (int i=0;i<n;i++) {
        for (int j = 0; j < m; j++) {

            double u = (2*i*M_PI) / n;
            double v = (2*j*M_PI) / n;

            points.emplace_back(Vector3D::point((R + r*std::cos(v)) * std::cos(u),
                                                (R + r*std::cos(v)) * std::sin(u),
                                                 R + r*std::sin(v)));
        }
    }

    for (int i=0;i<n;i++) {
        for (int j = 0; j < m; j++) {
            faces.emplace_back(Face(i*m + j, ((i+1)%n)*m + j, ((i+1)%n)*m  + (j + 1)%m, i*m + (j+1)%m));
        }
    }
}

// 3DLSystem functions
void Figure3D::create3DLSystem(std::string name, const ini::Configuration &conf) {
    // parse Lsystem file
    std::string input = conf[name]["inputfile"].as_string_or_die();
    std::ifstream input_stream(input);
    input_stream >> system;
    input_stream.close();

    // initialize variables
    H = Vector3D::vector(1, 0, 0);
    L = Vector3D::vector(0, 1, 0);
    U = Vector3D::vector(0, 0, 1);
    points.emplace_back(Vector3D::point(0, 0, 0));
    maxRecursionDepth = system.get_nr_iterations();
    currentAngle = convertToRad(system.get_angle());
    delta = convertToRad(system.get_angle());
    calculateLines(system.get_initiator());
}

void Figure3D::calculateLines(const std::string &input) {
    //std::cout << input << std::endl;
    for (char c:input) {
        // check if the alphabet contains the symbol, if so replace it
        if (system.get_alphabet().find(c) != system.get_alphabet().end() and recursionDepth < maxRecursionDepth) {
            recursionDepth++;
            calculateLines(system.get_replacement(c));
            continue;
        }
            // max recursion depth reached
        else if (c == '+') {
            Vector3D temp = H * std::cos(delta) + L * std::sin(delta);
            L = -H * std::sin(delta) + L * std::cos(delta);
            H = temp;
        } else if (c == '-') {
            Vector3D temp = H * std::cos(-delta) + L * std::sin(-delta);
            L = -H * std::sin(-delta) + L * std::cos(-delta);
            H = temp;
        } else if (c == '^') {
            Vector3D temp = H * std::cos(delta) + U * std::sin(delta);
            U = -H * std::sin(delta) + U * std::cos(delta);
            H = temp;
        } else if (c == '&') {
            Vector3D temp = H * std::cos(-delta) + U * std::sin(-delta);
            U = -H * std::sin(-delta) + U * std::cos(-delta);
            H = temp;
        } else if (c == '\\') {
            Vector3D temp = L * std::cos(delta) - U * std::sin(delta);
            U = L * std::sin(delta) + U * std::cos(delta);
            L = temp;
        } else if (c == '/') {
            Vector3D temp = L * std::cos(-delta) - U * std::sin(-delta);
            U = L * std::sin(-delta) + U * std::cos(-delta);
            L = temp;
        } else if (c == '|') {
            H = -H;
            L = -L;
        } else if (c == '(') {
            stackPoint3D temp;
            temp.H = H;
            temp.L = L;
            temp.U = U;
            temp.index = points.size() - 1;
            stackPoint.emplace_front(temp);
        } else if (c == ')') {
            H = stackPoint.front().H;
            L = stackPoint.front().L;
            U = stackPoint.front().U;
            points.emplace_back(Vector3D::point(points[stackPoint.front().index]));
            stackPoint.pop_front();
        } else if (system.get_alphabet().find(c) != system.get_alphabet().end()) {
            points.emplace_back(Vector3D::point(points[points.size() - 1] + H));
            if (system.draw(c)) {
                Face temp;
                temp.pointIndexes.emplace_back(points.size() - 2);
                temp.pointIndexes.emplace_back(points.size() - 1);
                faces.emplace_back(temp);
            }
        }
    }
    recursionDepth--;
}

// constructor
Figure3D::Figure3D(const std::string &name, const ini::Configuration &conf) {
    // read information from configuration file
    rotateX = conf[name]["rotateX"].as_double_or_die();
    rotateY = conf[name]["rotateY"].as_double_or_die();
    rotateZ = conf[name]["rotateZ"].as_double_or_die();
    scale = conf[name]["scale"].as_double_or_die();
    center = Vector3D::point(conf[name]["center"].as_double_tuple_or_die()[0],
                             conf[name]["center"].as_double_tuple_or_die()[1],
                             conf[name]["center"].as_double_tuple_or_die()[2]);
    eye = Vector3D::point(conf["General"]["eye"].as_double_tuple_or_die()[0],
                          conf["General"]["eye"].as_double_tuple_or_die()[1],
                          conf["General"]["eye"].as_double_tuple_or_die()[2]);

    // read in faces
    if (conf[name]["type"].as_string_or_die() == "LineDrawing") createLineDrawing(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Cube") createCube(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Octahedron") createOctahedron(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Tetrahedron") createTetrahedron(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Icosahedron") createIsocahedron(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Dodecahedron") createDodecahedron(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Sphere") createSphere(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Cone") createCone(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Cylinder") createCylinder(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "Torus") createTorus(name, conf);

    else if (conf[name]["type"].as_string_or_die() == "3DLSystem") create3DLSystem(name, conf);

    else std::cerr << "unknown figure type" << std::endl;

    // generate transformation matrix
    Matrix m;
    scaleMatrix(m, scale);
    rotateAroundX(m, convertToRad(rotateX));
    rotateAroundY(m, convertToRad(rotateY));
    rotateAroundZ(m, convertToRad(rotateZ));
    translateMatrix(m, center);
    m *= eyePointTrans(eye);
    applyTransformations(m);

    // apply transfaormations on points
    for (Vector3D &point:points) {
        doProjection(point, 1);
    }

    std::cout << points.size() << " | " << points2D.size() << "\n";
    // create lines
    for (const Face &face:faces) {
        for (uint index = 0; index < face.pointIndexes.size(); index++) {
            Line2D lineTemp;
            lineTemp.p1.x = points2D[face.pointIndexes[index]].x;
            lineTemp.p1.y = points2D[face.pointIndexes[index]].y;
            lineTemp.z1 = points[face.pointIndexes[index]].z;


            int n = face.pointIndexes[(index + 1) % face.pointIndexes.size()];
            lineTemp.p2.x = points2D[n].x;
            lineTemp.p2.y = points2D[n].y;
            lineTemp.z2 = points[n].z;


            lineTemp.color.ini(conf[name]["color"].as_double_tuple_or_die());
            lines2D.emplace_back(lineTemp);
        }
    }

}

// draw functions
const img::EasyImage Wireframe::drawLines2D() {


    double xmin = lines.front().p1.x;
    double xmax = lines.front().p1.y;
    double ymin = lines.front().p1.x;
    double ymax = lines.front().p1.y;
    for (const Line2D &line: lines) {
        // calculating xmin and xmax
        if (line.p1.x > xmax) xmax = line.p1.x;
        if (line.p2.x > xmax) xmax = line.p2.x;
        if (line.p1.x < xmin) xmin = line.p1.x;
        if (line.p2.x < xmin) xmin = line.p2.x;

        // calculating ymin and ymax
        if (line.p1.y > ymax) ymax = line.p1.y;
        if (line.p2.y > ymax) ymax = line.p2.y;
        if (line.p1.y < ymin) ymin = line.p1.y;
        if (line.p2.y < ymin) ymin = line.p2.y;
    }
    // calculating the imageSize of the image
    double imagex = 0;
    double imagey = 0;
    double rangex = 0;
    double rangey = 0;

    rangex = xmax - xmin;
    rangey = ymax - ymin;
    imagex = imageSize * (rangex / std::max(rangex, rangey));
    imagey = imageSize * (rangey / std::max(rangex, rangey));
    double d = 0.95 * (imagex / rangex);


    // move line drawing
    double dx, dy;
    dx = (imagex / 2) - (d * ((xmin + xmax) / 2));
    dy = (imagey / 2) - (d * ((ymin + ymax) / 2));

    // draw the lines
    img::EasyImage image(roundToInt(imagex), roundToInt(imagey));
    image.clear(img::Color(backgroundcolor.red * 255, backgroundcolor.green * 255, backgroundcolor.blue * 255));
    for (const Line2D &line: lines) {
        image.draw_line(roundToInt((line.p1.x * d) + dx), roundToInt((line.p1.y * d) + dy),
                        roundToInt((line.p2.x * d) + dx), roundToInt((line.p2.y * d) + dy),
                        img::Color(line.color.red * 255, line.color.green * 255, line.color.blue * 255));
    }
    return image;
}

img::EasyImage Wireframe::drawWireFrame(const ini::Configuration &conf) {
    // read information from configuration file
    imageSize = conf["General"]["size"].as_int_or_die();
    nrOfFigures = conf["General"]["nrFigures"].as_int_or_die();
    backgroundcolor.ini(conf["General"]["backgroundcolor"].as_double_tuple_or_die());

    for (int k = 0; k < nrOfFigures; k++) {
        Figure3D temp("Figure" + std::to_string(k), conf);
        temp.addLines2D(lines);
    }

    return drawLines2D();
}
