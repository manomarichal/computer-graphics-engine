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

#ifndef CGENGINE_WIREFRAME_H
#define CGENGINE_WIREFRAME_H

#include "Lines3D.h"
#include <cmath>
#include <forward_list>
#include <string>
#include <iostream>
#include <map>
#include "vector3d.h"
#include "easy_image.h"
#include "ini_configuration.h"

class Wireframe {
    Color backgroundcolor;
    int imageSize, nrOfFigures;
    listWithLines lines;
    std::vector<std::vector<Figure3D>> allFigures;
    std::vector<Figure3D> figures;

    // draw figure functions
    const img::EasyImage drawLines2D(bool zBuffered);

    const img::EasyImage drawZBufferedTriangles();

    std::vector<Figure3D> createFractal(std::string name, const ini::Configuration &conf, Figure3D &fig);

    void isFractal(std::string name, const ini::Configuration &conf, Figure3D &figure);

    void createMengerSponge(std::string name, const ini::Configuration &conf, Figure3D &fig);

    // menger sponge helper functions
    std::vector<Figure3D> splitSponge(Figure3D &root);

    void createSmallCube(Figure3D &tempFig, int a, int b);

    // light functions
    void initLights(Figure3D &fig);


public:
    img::EasyImage drawWireFrame(const ini::Configuration &conf, bool zBuffered, bool zBuffTriangle, bool light);

};



#endif //CGENGINE_WIREFRAME_H
