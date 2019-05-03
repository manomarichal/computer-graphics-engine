//============================================================================
// @name        : 
// @author      : Mano Marichal
// @date        : 
// @version     : 
// @copyright   : BA1 Informatica - Mano Marichal - University of Antwerp
// @description : 
//============================================================================

//
// Created by mano on 03.05.19.
//

#include "Wireframe.h"


// draw functions
const img::EasyImage Wireframe::drawLines2D(bool zBuffered)
{

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

    double rangex = xmax - xmin;
    double rangey = ymax - ymin;
    double imagex = imageSize * (rangex / std::max(rangex, rangey));
    double imagey = imageSize * (rangey / std::max(rangex, rangey));
    double d = 0.95 * (imagex / rangex);


    // move line drawing
    double dx, dy;
    dx = (imagex / 2) - (d * ((xmin + xmax) / 2));
    dy = (imagey / 2) - (d * ((ymin + ymax) / 2));

    // draw the lines
    img::EasyImage image(roundToInt(imagex), roundToInt(imagey));
    image.clear(img::Color(backgroundcolor.red * 255, backgroundcolor.green * 255, backgroundcolor.blue * 255));

    ZBuffer zBuf(roundToInt(imagex), roundToInt(imagey));

    for (const Line2D &line: lines) {

        if (!zBuffered) {
            image.draw_line(roundToInt((line.p1.x * d) + dx), roundToInt((line.p1.y * d) + dy),
                            roundToInt((line.p2.x * d) + dx), roundToInt((line.p2.y * d) + dy),
                            img::Color(line.color.red * 255, line.color.green * 255, line.color.blue * 255));
        }
        else image.draw_zbuf_line(zBuf,
                                  roundToInt((line.p1.x * d) + dx), roundToInt((line.p1.y * d) + dy), line.z1,
                                  roundToInt((line.p2.x * d) + dx), roundToInt((line.p2.y * d) + dy), line.z2,
                                  img::Color(line.color.red * 255, line.color.green * 255, line.color.blue * 255));

    }
    return image;
}

const img::EasyImage Wireframe::drawZBufferedTriangles(const ini::Configuration &conf)
{
    double xmin = INT64_MAX;
    double xmax = INT64_MIN;
    double ymin = INT64_MAX;
    double ymax = INT64_MIN;

    for (auto &figure:figures) {

        figure.createTriangles();

        for (auto &point:figure.points2D) {

            if (point.x < xmin) xmin = point.x;
            if (point.x > xmax) xmax = point.x;
            if (point.y < ymin) ymin = point.y;
            if (point.y > ymax) ymax = point.y;

        }
    }

    double imagex; double imagey; double rangex; double rangey;

    rangex = xmax - xmin;
    rangey = ymax - ymin;
    imagex = imageSize * (rangex / std::max(rangex, rangey));
    imagey = imageSize * (rangey / std::max(rangex, rangey));

    double d = 0.95 * (imagex / rangex);


    // move line drawing
    double dx, dy;
    dx = (imagex / 2) - (d * ((xmin + xmax) / 2));
    dy = (imagey / 2) - (d * ((ymin + ymax) / 2));

    // make image and zbuffer
    img::EasyImage image(roundToInt(imagex), roundToInt(imagey));
    image.clear(img::Color(backgroundcolor.red * 255, backgroundcolor.green * 255, backgroundcolor.blue * 255));

    ZBuffer zBuf(roundToInt(imagex), roundToInt(imagey));

    for (auto &figure:figures) {

        for (auto &face:figure.faces) {

            image.img::EasyImage::draw_zbuf_triangle(zBuf,
                                                     figure.points[face.pointIndexes[0]],
                                                     figure.points[face.pointIndexes[1]],
                                                     figure.points[face.pointIndexes[2]],
                                                     d, dx, dy,
                                                     img::Color(figure.color.red * 255, figure.color.green * 255, figure.color.blue * 255));

        }
    }
    return image;
}

void Wireframe::isFractal(std::string name, const ini::Configuration &conf, bool zBuf, Figure3D &figure)
{
    int nrIterations = conf[name]["nrIterations"].as_int_or_die();

    std::vector<Figure3D> temp;
    std::vector<Figure3D> temp2;

    temp.emplace_back(figure);

    for (int i=0;i<nrIterations;++i)
    {

        for (auto &figure:temp)
        {
            for(auto &newFigure:createFractal(name, conf, figure))
            {
                temp2.emplace_back(newFigure);
            }
        }

        temp = temp2;
        temp2.clear();
    }

    for (auto &f:temp)
    {
        for (Vector3D &point:f.points)
        {
            f.doProjection(point, 1);
        }

        if (!zBuf) f.createLinesOutOfFaces(name,conf);

        figures.emplace_back(f);
    }
}

std::vector<Figure3D> Wireframe::createFractal(std::string name, const ini::Configuration &conf, Figure3D &fig)
{
    double scale = conf[name]["fractalScale"].as_double_or_die();

    std::vector <Figure3D> fractal;

    for (int i = 0; i < fig.points.size(); ++i) {
        Figure3D tempFig;
        tempFig.points = fig.points;
        tempFig.faces = fig.faces;

        Matrix m;
        fig.scaleMatrix(m, 1 / scale);
        tempFig.applyTransformations(m);

        Matrix p;
        fig.translateMatrix(p, Vector3D::vector(fig.points[i] - tempFig.points[i]));
        tempFig.applyTransformations(p);

        fractal.emplace_back(tempFig);
    }

    return fractal;
}

void Wireframe::createSmallCube(Figure3D &tempFig, int a, int b)
{
}
std::vector<Figure3D> Wireframe::splitSponge(Figure3D &root)
{

    std::vector<Figure3D> splittedCubes;

    double scale = 3;

    for (int i = 0; i < root.points.size(); ++i) {
        Figure3D tempFig;
        tempFig.points = root.points;
        tempFig.faces = root.faces;

        Matrix m;
        root.scaleMatrix(m, 1 / scale);
        tempFig.applyTransformations(m);

        Matrix p;
        root.translateMatrix(p, Vector3D::vector(root.points[i] - tempFig.points[i]));
        tempFig.applyTransformations(p);

        splittedCubes.emplace_back(tempFig);
    }

    for (int i = 0; i < root.points.size(); ++i) {

        Figure3D tempFig;
        tempFig.points = root.points;
        tempFig.faces = root.faces;

        Matrix m;
        root.scaleMatrix(m, 1 / scale);
        tempFig.applyTransformations(m);

        int index = i + 1;
        if (i == 3) index = 0;
        if (i == 7) index = 4;


        Matrix k;
        root.translateMatrix(k, Vector3D::vector(root.points[i] - tempFig.points[i]));
        tempFig.applyTransformations(k);

        Matrix p;
        root.translateMatrix(p, Vector3D::vector( Vector3D::point(root.points[i].x + (root.points[index].x - root.points[i].x)/3,
                                                                  root.points[i].y + (root.points[index].y - root.points[i].y)/3,
                                                                  root.points[i].z + (root.points[index].z - root.points[i].z)/3) - root.points[i]));
        tempFig.applyTransformations(p);

        splittedCubes.emplace_back(tempFig);
    }


    for (int i=0;i<4;i++)
    {
        Figure3D tempFig;
        tempFig.points = root.points;
        tempFig.faces = root.faces;

        Matrix m;
        root.scaleMatrix(m, 1 / scale);
        tempFig.applyTransformations(m);

        int index;

        if (i == 0) index = 5;
        else if (i == 1) index = 4;
        else if (i == 2) index = 7;
        else if (i == 3) index = 6;

        Matrix k;
        root.translateMatrix(k, Vector3D::vector(root.points[i] - tempFig.points[i]));
        tempFig.applyTransformations(k);

        Matrix p;
        root.translateMatrix(p, Vector3D::vector( Vector3D::point(root.points[i].x + (root.points[index].x - root.points[i].x)/3,
                                                                  root.points[i].y + (root.points[index].y - root.points[i].y)/3,
                                                                  root.points[i].z + (root.points[index].z - root.points[i].z)/3) - root.points[i]));
        tempFig.applyTransformations(p);

        splittedCubes.emplace_back(tempFig);
    }

    return splittedCubes;

}
void Wireframe::createMengerSponge(std::string name, const ini::Configuration &conf, Figure3D &root)
{
    int nrIterations = conf[name]["nrIterations"].as_int_or_die();

    std::vector<Figure3D> temp;
    std::vector<Figure3D> temp2;

    temp.emplace_back(root);

    for (int i=0;i<nrIterations;++i)
    {

        for (auto &figure:temp)
        {
            for(auto &newFigure:splitSponge(figure))
            {
                temp2.emplace_back(newFigure);
            }
        }

        temp = temp2;
        temp2.clear();
    }

    for (auto &f:temp)
    {
        for (Vector3D &point:f.points)
        {
            f.doProjection(point, 1);
        }

        f.createLinesOutOfFaces(name,conf);

        figures.emplace_back(f);
    }
}


img::EasyImage Wireframe::drawWireFrame(const ini::Configuration &conf, bool zBuffered, bool zBuffTriangle)
{
    // read information from configuration file
    imageSize = conf["General"]["size"].as_int_or_die();
    nrOfFigures = conf["General"]["nrFigures"].as_int_or_die();
    backgroundcolor.ini(conf["General"]["backgroundcolor"].as_double_tuple_or_die());
    std::string type = conf["General"]["type"].as_string_or_die();


    for (int k = 0; k < nrOfFigures; k++)
    {
        std::string name = "Figure" + std::to_string(k);

        Figure3D temp(name, conf, zBuffTriangle);

        if (conf[name]["type"].as_string_or_die().substr(0, 7) == "Fractal") // FRACTAL
        {
            isFractal(name, conf, zBuffTriangle, temp);
        }
        else if (conf[name]["type"].as_string_or_die() == "MengerSponge")
        {
            createMengerSponge(name, conf, temp);
        }
        else figures.emplace_back(temp);

    }

    if (!zBuffTriangle)
    {
        for (auto &f:figures)
        {
            f.addLines2D(lines);
        }
    }
    else return drawZBufferedTriangles(conf);

    return drawLines2D(zBuffered);
}