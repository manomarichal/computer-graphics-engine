//============================================================================
// @name        : Lines3D.h
// @author      : Mano Marichal
// @date        :
// @version     :
// @copyright   : Computer Graphics - BA1 Informatica - Mano Marichal - University of Antwerp
// @description : Class describing lines draw in 3D
//============================================================================
#include "Lines3D.h"

void Figure3D::rotateAroundX(Matrix &m, const double angle) {
    Matrix s;
    s(3,3) = std::cos(angle);
    s(2,2) = std::cos(angle);
    s(3,2) = std::sin(angle);
    s(2,3) = -1 * std::sin(angle);
    //std::cout << "rotation aroud x s:" << std::endl; s.print(std::cout);
    m *= s;
}

void Figure3D::rotateAroundY(Matrix &m, const double angle) {
    Matrix s;
    s(1,1) = std::cos(angle);
    s(3,3) = std::cos(angle);
    s(1,3) = std::sin(angle);
    s(3,1) = -1 * std::sin(angle);
    //std::cout << "rotation aroud Y s:" << std::endl; s.print(std::cout);
    m *= s;
}

void Figure3D::rotateAroundZ(Matrix &m, const double angle) {
    Matrix s;
    std::cout << angle << std::endl;
    s(1,1) = std::cos(angle);
    s(2,2) = std::cos(angle);
    s(2,1) = std::sin(angle);
    s(1,2) = -1 * std::sin(angle);
    //std::cout << "rotation aroud Z s:" << std::endl; s.print(std::cout);
    m *= s;
}

void Figure3D::scaleMatrix(Matrix &m, const double scale) {
    Matrix s;
    s(1,1) = scale;
    s(2,2) = scale;
    s(3,3) = scale;
    m *= s;
}

void Figure3D::translateMatrix(Matrix &m, const Vector3D &v){
    Matrix s;
    s(4,1) = v.x;
    s(4,2) = v.y;
    s(4,3) = v.z;
    m *= s;

}

void Figure3D::applyTransformations(const Matrix &m) {
    for (auto &p:points) {
        //std::cout << "normal" << std::endl;
        //std::cout << p.x << "|" << p.y << std::endl;
        //std::cout << "transformations" << std::endl;
        p *= m;
        //std::cout << p.x << "|" << p.y << std::endl;
    }
    //std::cout << "all points have been calculated" << std::endl;
}

const std::vector<Vector3D> & Figure3D::getPoints() const {
    return points;
}

void Figure3D::toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt((point.x*point.x) + (point.y*point.y) + (point.z*point.z));
    theta = std::atan2(point.y,point.x);
    phi = std::acos(point.z/r);
    //std::cout << "theta: " << theta << "phi: " << phi << "r: " << r << std::endl;
}


Matrix Figure3D::eyePointTrans(const Vector3D &eyepoint) {
    double theta,phi,r;
    toPolar(eyepoint,theta,phi,r);
    Vector3D v = Vector3D::vector(0,0,-r);
    Matrix m;
    rotateAroundZ(m, (M_PI/2) + theta);
    m.print(std::cout);
    rotateAroundX(m, phi);
    m.print(std::cout);
    translateMatrix(m,v);
    m.print(std::cout);
    return m;
}

void Figure3D::doProjection(const Vector3D &point, const double d) {
    Point2D newPoint;
    newPoint.x = (d*point.x) / (-point.z);
    newPoint.y = (d*point.y) / (-point.z);
    //std::cout << "newpoint:" << newPoint.x << "|" << newPoint.y << std::endl;
    points2D.emplace_back(newPoint);
}

void Figure3D::addLines2D(listWithLines &list) {
    for (const Line2D &line:lines2D) {
        list.emplace_front(line);
    }
}

Figure3D::Figure3D(const std::string &name, const ini::Configuration &conf) {
    // read information from configuration file
    //std::cout << name << std::endl;
    if (conf[name]["type"].as_string_or_die() == "LineDrawing") {
        rotateX = conf[name]["rotateX"].as_double_or_die();
        rotateY = conf[name]["rotateY"].as_double_or_die();
        rotateZ = conf[name]["rotateZ"].as_double_or_die();
        center = Vector3D::point(conf[name]["center"].as_double_tuple_or_die()[0],
                                 conf[name]["center"].as_double_tuple_or_die()[1],
                                 conf[name]["center"].as_double_tuple_or_die()[2]);
        eye = Vector3D::point(conf["General"]["eye"].as_double_tuple_or_die()[0],
                              conf["General"]["eye"].as_double_tuple_or_die()[1],
                              conf["General"]["eye"].as_double_tuple_or_die()[2]);
        scale = conf[name]["scale"].as_double_or_die();
        nrOfPoints = conf[name]["nrPoints"].as_int_or_die();
        nrOfLines = conf[name]["nrLines"].as_int_or_die();

        // read in points
        for (int k=0;k<nrOfPoints;k++) {
            // for each line
            Vector3D temp = Vector3D::point(conf[name]["point" + std::to_string(k)].as_double_tuple_or_die()[0],
                                   conf[name]["point" + std::to_string(k)].as_double_tuple_or_die()[1],
                                   conf[name]["point" + std::to_string(k)].as_double_tuple_or_die()[2]);
            temp.print(std::cout);
            points.emplace_back(temp);
            std::cout << k << std::endl;
        }

        // generate transformation matrix
        Matrix m;
        //std::cout << "eye"<< std::endl; m.print(std::cout); std::cout << std::endl;
        scaleMatrix(m, scale);
        //std::cout << "scale" << std::endl; m.print(std::cout); std::cout << std::endl;
        rotateAroundX(m, convertToRad(rotateX));
        //std::cout << "rx"<< std::endl; m.print(std::cout); std::cout << std::endl;
        rotateAroundY(m, convertToRad(rotateY));
        //std::cout << "ry"<< std::endl; m.print(std::cout); std::cout << std::endl;
        rotateAroundZ(m, convertToRad(rotateZ));
        //std::cout << "rz"<< std::endl; m.print(std::cout); std::cout << std::endl;
        translateMatrix(m, center);
        //std::cout << "trans"<< std::endl; m.print(std::cout); std::cout << std::endl;
        m*=eyePointTrans(eye);
        applyTransformations(m);
        // apply transformations on points
        for (Vector3D &point:points) {
           doProjection(point, 1);
        }

        // create lines
        for (int k=0; k<nrOfLines;k++) {
            Line2D temp;
            temp.p1.x = points2D[conf[name]["line" + std::to_string(k)].as_double_tuple_or_die()[0]].x;
            temp.p1.y = points2D[conf[name]["line" + std::to_string(k)].as_double_tuple_or_die()[0]].y;
            temp.p2.x = points2D[conf[name]["line" + std::to_string(k)].as_double_tuple_or_die()[1]].x;
            temp.p2.y = points2D[conf[name]["line" + std::to_string(k)].as_double_tuple_or_die()[1]].y;
            //std::cout << "newline:" << temp.p1.x << "|" << temp.p1.y << "|" << temp.p2.x << "|" << temp.p2.y << std::endl;
            temp.color.ini(conf[name]["color"].as_double_tuple_or_die());
            lines2D.emplace_front(temp);
        }
    }
}

const img::EasyImage Wireframe::drawLines2D() {
    double xmin = INT64_MAX;
    double xmax = INT64_MIN;
    double ymin = INT64_MAX;
    double ymax = INT64_MIN;
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
    dx = imagex / 2 - (d * ((xmin + xmax) / 2));
    dy = imagey / 2 - (d * ((ymin + ymax) / 2));

    // draw the lines
    img::EasyImage image(roundToInt(imagex), roundToInt(imagey));
    image.clear(img::Color(backgroundcolor.red*255, backgroundcolor.green*255, backgroundcolor.blue*255));
    for (const Line2D &line: lines) {
        //std::cout << "scaled:\n" << roundToInt((line.p1.x * d) + dx) << "|" << roundToInt((line.p1.y * d) + dy)<< "|" << roundToInt((line.p2.x * d) + dx)<< "|" << roundToInt((line.p2.y * d) + dy) << std::endl;
        image.draw_line(roundToInt((line.p1.x * d) + dx), roundToInt((line.p1.y * d) + dy),
                        roundToInt((line.p2.x * d) + dx), roundToInt((line.p2.y * d) + dy),
                        img::Color(line.color.red*255, line.color.green*255, line.color.blue*255));
    }
    return image;
}

img::EasyImage Wireframe::drawWireFrame(const ini::Configuration &conf) {
    // read information from configuration file
    imageSize = conf["General"]["size"].as_int_or_die();
    nrOfFigures = conf["General"]["nrFigures"].as_int_or_die();
    backgroundcolor.ini(conf["General"]["backgroundcolor"].as_double_tuple_or_die());


    for (int k=0;k<nrOfFigures;k++) {
        Figure3D temp("Figure" + std::to_string(k), conf);
        temp.addLines2D(lines);
        figures.emplace_front(temp);
    }

    return drawLines2D();
}