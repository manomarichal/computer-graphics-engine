//
// Created by mano on 22.02.19.
//

#ifndef CGENGINE_LINE2D_H
#define CGENGINE_LINE2D_H

#include "easy_image.h"
#include "l_parser.h"
#include "ini_configuration.h"

#include <forward_list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stack>

inline int roundToInt(double d){
    return static_cast<int>(std::round(d));
}

inline double convertToRad(double a) {
    return ((a * M_PI) / 180.0);
}

struct Point2D {
    double x,y;

};

struct Line2D {
    Point2D p1;
    Point2D p2;
    img::Color color;
};

typedef std::forward_list<Line2D> listWithLines;


class Lines2D {
private:
    LParser::LSystem2D system;                                  // the parsed system
    std::vector<double> backgroundcolor;                        // color used for the background
    std::vector<double> color;                                  // color used for the lines
    Point2D prev, cur;                                          // points used to make the lines
    int recursionDepth=0; int maxRecursionDepth;                // recursion depth for replacement
    double currentAngle;                                        // the current angle
    int imageSize;                                              // size of the image
    listWithLines lines;                                        // std::forward_list containing all lines
    std::forward_list<std::pair<double, Point2D>> stackPoint;   // stack for remembering return points when encountering '(' or ')'

    // calculates all lines for the 2D l-system and stores them in the member 'lines'
    void calculateLines(const std::string &input);


public:
    // constructors
    Lines2D(const ini::Configuration &conf);

    // draws the current set of lines in the member 'lines' on an image with size
    const img::EasyImage drawLines2D(const listWithLines &list, int size);

    // calculates all lines for the 2D l-system and draws them on an image
    img::EasyImage drawLSystem2D();


};
#endif //CGENGINE_LINE2D_H
/*
        // read information from configuration file
        int size = conf["General"]["size"].as_int_or_die();
        std::vector<double> backgroundcolor = conf["General"]["backgroundcolor"].as_double_tuple_or_die();
        std::string input = conf["2DLSystem"]["inputfile"].as_string_or_die();
        std::vector<double> color = conf["2DLSystem"]["color"].as_double_tuple_or_die();

        // parse Lsystem file and make a Lines2D with the needed information
        LParser::LSystem2D system;
        std::ifstream input_stream(input);
        input_stream >> system;
        input_stream.close();
        */