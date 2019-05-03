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
    std::vector<Figure3D> figures;

    // draw figure functions
    const img::EasyImage drawLines2D(bool zBuffered);

    const img::EasyImage drawZBufferedTriangles(const ini::Configuration &conf);

public:
    img::EasyImage drawWireFrame(const ini::Configuration &conf, bool zBuffered, bool zBuffTriangle);

    std::vector<Figure3D> createFractal(std::string name, const ini::Configuration &conf, Figure3D &fig);

    void isFractal(std::string name, const ini::Configuration &conf, bool zBuf, Figure3D &figure);

};



#endif //CGENGINE_WIREFRAME_H
