//
// Created by mano on 22.02.19.
//

#ifndef CGENGINE_LINE2D_H
#define CGENGINE_LINE2D_H

#include "l_parser.h"
#include "ini_configuration.h"
#include "easy_image.h"

#include <forward_list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stack>

struct Point2D {
    double x,y;
};

struct Color {
    double red, green, blue;
};

struct Line2D {
    Point2D p1;
    Point2D p2;
    Color color;
};

typedef std::forward_list<Line2D> listWithLines;

inline int roundToInt(double d){
    return static_cast<int>(std::round(d));
}

inline double convertToRad(double a) {
    return ((a * M_PI) / 180.0);
}


class Lines2D {
private:
    LParser::LSystem2D system;
    std::vector<double> backgroundcolor;
    std::vector<double> color;
    Point2D prev, cur;
    int recursionDepth=0; int maxRecursionDepth;
    double currentAngle;
    int imageSize;
    listWithLines lines;
    std::forward_list<std::pair<double, Point2D>> stackPoint;

    // calculates all lines for the 2D l-system and stores them in the member 'lines'
    void calculateLines(const std::string &input);


public:
    // constructor
    Lines2D(const ini::Configuration &conf);

    // draws the current set of lines in the member 'lines' on an image with size
    const img::EasyImage drawLines2D(const listWithLines &list, int size);

    // calculates all lines for the 2D l-system and draws them on an image
    img::EasyImage drawLSystem2D();


};
#endif //CGENGINE_LINE2D_H
