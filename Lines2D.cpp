//
// Created by mano on 22.02.19.
//


#include "Lines2D.h"

// Lines2D
const img::EasyImage Lines2D::drawLines2D(const listWithLines &list, int size) {
    double xmin = INT64_MAX;
    double xmax = INT64_MIN;
    double ymin = INT64_MAX;
    double ymax = INT64_MIN;
    for (const Line2D &line: list) {
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
    // calculating the size of the image
    double imagex = 0;
    double imagey = 0;
    double rangex = 0;
    double rangey = 0;

    rangex = xmax - xmin;
    rangey = ymax - ymin;
    imagex = size * (rangex / std::max(rangex, rangey));
    imagey = size * (rangey / std::max(rangex, rangey));
    double d = 0.95 * (imagex / rangex);

    // move line drawing
    double dx, dy;
    dx = imagex / 2 - (d * ((xmin + xmax) / 2));
    dy = imagey / 2 - (d * ((ymin + ymax) / 2));

    // draw the lines
    img::EasyImage image(roundToInt(imagex), roundToInt(imagey));
    image.clear(img::Color(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]));
    for (const Line2D line: list) {
        image.draw_line(roundToInt((line.p1.x * d) + dx), roundToInt((line.p1.y * d) + dy),
                        roundToInt((line.p2.x * d) + dx), roundToInt((line.p2.y * d) + dy),
                        img::Color(line.color.red, line.color.blue, line.color.green));
    }
    return image;
}

// LSystem2D
Lines2D::Lines2D(const ini::Configuration &conf) {
    // read information from configuration file
    imageSize = conf["General"]["size"].as_int_or_die();
    backgroundcolor = conf["General"]["backgroundcolor"].as_double_tuple_or_die();
    color = conf["2DLSystem"]["color"].as_double_tuple_or_die();
    for (u_long i = 0; i < color.size(); i++) color[i] = color[i] * 255.99;
    for (u_long i = 0; i < backgroundcolor.size(); i++) backgroundcolor[i] = backgroundcolor[i] * 255.99;

    // parse Lsystem file
    std::string input = conf["2DLSystem"]["inputfile"].as_string_or_die();
    std::ifstream input_stream(input);
    input_stream >> system;
    input_stream.close();

    maxRecursionDepth = system.get_nr_iterations();
    cur.x = 0;
    cur.y = 0;
    prev.x = 0;
    prev.y = 0;
    currentAngle = convertToRad(system.get_starting_angle());
}

img::EasyImage Lines2D::drawLSystem2D() {
    this->calculateLines(system.get_initiator());
    return drawLines2D(lines, imageSize);
}

void Lines2D::calculateLines(const std::string &input) {
    for (char c:input) {
        // check if the alphabet contains the symbol, if so replace it
        if (system.get_alphabet().find(c) != system.get_alphabet().end() and recursionDepth < maxRecursionDepth) {
            recursionDepth++;
            calculateLines(system.get_replacement(c));
            continue;
        }
            // max recursion depth reached
        else if (c == '+') currentAngle = currentAngle + convertToRad(system.get_angle());
        else if (c == '-') currentAngle = currentAngle - convertToRad(system.get_angle());
        else if (c == '(') stackPoint.emplace_front(std::pair < double, Point2D > {currentAngle, cur});
        else if (c == ')') {
            cur = stackPoint.front().second;
            prev = cur;
            currentAngle = stackPoint.front().first;
            stackPoint.pop_front();
        } else if (system.get_alphabet().find(c) != system.get_alphabet().end()) {
            cur.x = cur.x + std::cos(currentAngle);
            cur.y = cur.y + std::sin(currentAngle);
            if (system.draw(c)) {
                Line2D temp;
                temp.p1 = cur;
                temp.p2 = prev;
                temp.color = {color[0], color[1], color[2]};
                lines.emplace_front(temp);
            }
            prev = cur;
        }
    }
    recursionDepth--;
}

