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

const img::EasyImage Wireframe::drawLines2D(bool zBuffered)
{

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
            //std::cout << "\nfrom" << line.p1.x << " | " << line.p1.y << "to" << line.p2.x << " | " << line.p2.y <<std::endl;
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
void Wireframe::createLightZBuffer(std::vector<Figure3D> &figs, Light &light)
{
    double xmin = INT64_MAX;
    double xmax = INT64_MIN;
    double ymin = INT64_MAX;
    double ymax = INT64_MIN;

    for (auto &figure:figs) {
        for (auto &point:figure.points2D) {

            if (point.x < xmin) xmin = point.x;
            if (point.x > xmax) xmax = point.x;
            if (point.y < ymin) ymin = point.y;
            if (point.y > ymax) ymax = point.y;

        }
    }

    double rangex = xmax - xmin;
    double rangey = ymax - ymin;

    double imagex = shadowSize * (rangex / std::max(rangex, rangey));
    double imagey = shadowSize * (rangey / std::max(rangex, rangey));

    light.d = 0.95 * (imagex / rangex);

    // move line drawing
    light.dx = (imagex / 2) - (light.d * ((xmin + xmax) / 2));
    light.dy = (imagey / 2) - (light.d * ((ymin + ymax) / 2));

    light.shadowMask = ZBuffer(roundToInt(imagex), roundToInt(imagey));

    for (auto &figure:figs) {
        for (auto &face:figure.faces) {
            img::EasyImage::draw_zbuf_triangle_colorless(light.shadowMask,
                                                         figure.points[face.pointIndexes[0]],
                                                         figure.points[face.pointIndexes[1]],
                                                         figure.points[face.pointIndexes[2]],
                                                         light.d, light.dx, light.dy);
        }
    }
    //createImageZBuffer(light.shadowMask);

}
void Wireframe::createImageZBuffer(const ZBuffer &zbuffer)
{
    img::EasyImage img(zbuffer.xW, zbuffer.yH);
    for (int i=0; i < zbuffer.xW;i++)
    {
        for (int j=0; j < zbuffer.yH;j++)
        {
            if (zbuffer.getZVal(i,j) == std::numeric_limits<double>::infinity()) continue;

            img.operator()(i, j) = img::Color(255, 0, 0);
        }
    }
    std::ofstream f_out("shadowing264_ZBUFFER.bmp");
    f_out << img;
}
const img::EasyImage Wireframe::drawZBufferedTriangles()
{

    double xmin = INT64_MAX;
    double xmax = INT64_MIN;
    double ymin = INT64_MAX;
    double ymax = INT64_MIN;

    for (auto &figure:figures) {

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

    Matrix m = Matrix::inv(Figure3D::eyePointTrans(eye));


    for (auto &figure:figures)
    {
        std::vector<Light> *lights;
        if (wireframeLights.size() != 0) lights = &wireframeLights;
        else lights = &figure.lights;
            for (auto &face:figure.faces)
        {
            std::vector<double> ambient = figure.ambientReflection.asVector();
            std::vector<double> diffuse = figure.diffuseReflection.asVector();
            std::vector<double> specular = figure.specularReflection.asVector();

            image.img::EasyImage::draw_zbuf_triangle(zBuf,
                                                     figure.points[face.pointIndexes[0]],
                                                     figure.points[face.pointIndexes[1]],
                                                     figure.points[face.pointIndexes[2]],
                                                     d, dx, dy,
                                                     ambient, diffuse, specular,
                                                     figure.reflectionCoefficient, *lights,
                                                     m, applyShadows);
        }
    }
    return image;
}
void Wireframe::initLights(const ini::Configuration &conf)
{
    std::vector<Light> light;
    std::vector<double> defaultTuple = {-1, -1, -1};

    for (int k=0;k<conf["General"]["nrLights"].as_int_or_die();k++)
    {
        Color temp;
        Light tempLight;
        std::string name = "Light" + std::to_string(k);

        // ambient light
        if (conf[name]["ambientLight"].as_double_tuple_if_exists(defaultTuple))
        {
            temp.ini(conf[name]["ambientLight"].as_double_tuple_or_die());
            tempLight.ambientLight.iniColor(temp);
            tempLight.amLight = true;
        }

        // diffuse light
        if (conf[name]["diffuseLight"].as_double_tuple_if_exists(defaultTuple))
        {
            temp.ini(conf[name]["diffuseLight"].as_double_tuple_or_default(defaultTuple));
            tempLight.diffuseLight.iniColor(temp);
            tempLight.infinity = conf[name]["infinity"].as_bool_or_default(false);
            tempLight.difLight = true;

            if (tempLight.infinity)
            {
                tempLight.ldVector = Vector3D::vector(conf[name]["direction"].as_double_tuple_or_die()[0],
                                                      conf[name]["direction"].as_double_tuple_or_die()[1],
                                                      conf[name]["direction"].as_double_tuple_or_die()[2]);
            }
            else tempLight.ldVector = Vector3D::point(conf[name]["location"].as_double_tuple_or_die()[0],
                                                       conf[name]["location"].as_double_tuple_or_die()[1],
                                                       conf[name]["location"].as_double_tuple_or_die()[2]);

            if (applyShadows)
            {
                tempLight.eye = Figure3D::eyePointTrans(tempLight.ldVector);
            }

            tempLight.ldVector *= Figure3D::eyePointTrans(eye);

            if (conf[name]["specularLight"].as_double_tuple_if_exists(defaultTuple))
            {
                tempLight.specularLight.ini(conf[name]["specularLight"].as_double_tuple_or_die());
                tempLight.specLight = true;
            }
        }


        light.emplace_back(tempLight);

    }
    wireframeLights = light;
}
img::EasyImage Wireframe::drawWireFrame(const ini::Configuration &conf, bool zBuffered, bool zBuffTriangle, bool light)
{
    // read information from configuration file
    imageSize = conf["General"]["size"].as_int_or_die();
    nrOfFigures = conf["General"]["nrFigures"].as_int_or_die();
    backgroundcolor.ini(conf["General"]["backgroundcolor"].as_double_tuple_or_die());
    std::string type = conf["General"]["type"].as_string_or_die();
    applyShadows = conf["General"]["shadowEnabled"].as_bool_or_default(false);

    eye = Vector3D::point(conf["General"]["eye"].as_double_tuple_or_die()[0],
                          conf["General"]["eye"].as_double_tuple_or_die()[1],
                          conf["General"]["eye"].as_double_tuple_or_die()[2]);

    // belichting
    if (conf["General"]["type"].as_string_or_die().substr(0, 7) == "Lighted")
    {
        initLights(conf);
    }

    for (int k = 0; k < nrOfFigures; k++)
    {
        std::string name = "Figure" + std::to_string(k);

        Figure3D temp(name, conf, zBuffTriangle, light);

        if (conf[name]["type"].as_string_or_die().substr(0, 7) == "Fractal") // FRACTAL
        {
            if (conf[name]["type"].as_string_or_die() == "FractalBuckyBall") return drawLines2D(zBuffered);

            isFractal(name, conf, temp);
        }
        else if (conf[name]["type"].as_string_or_die() == "MengerSponge")
        {
            createMengerSponge(name, conf, temp);
        }
        else
        {
            std::vector<Figure3D> f;
            f.emplace_back(temp);
            allFigures.emplace_back(f);
        }
    }

    // do thick stuff
    for (unsigned int n = 0; n < allFigures.size(); n++)
    {
        std::string name = "Figure" + std::to_string(n);

        if (conf[name]["type"].as_string_or_die().substr(0, 5) != "Thick" ) continue;

        std::vector<Figure3D> temp;

        for (auto &figure:allFigures[n])
        {
            for (auto &thiccFigure:generateThiccFigure(name, conf, figure))
            {
                temp.emplace_back(thiccFigure);
            }
        }
        allFigures[n] = temp;
    }

    // CREATING SHADOWMASKS
    if (applyShadows)
    {
        shadowSize = conf["General"]["shadowMask"].as_int_or_die();
        for (auto &wlight:wireframeLights)
        {
            if (wlight.infinity) continue;

            std::vector<Figure3D> figs;

            for (unsigned int n = 0; n < allFigures.size(); n++)
            {
                for (auto figure:allFigures[n])
                {
                    figure.applyTransformations(wlight.eye);

                    figure.createTriangles();

                    for (Vector3D &point:figure.points)
                    {
                        figure.doProjection(point, 1);
                    }

                    figs.emplace_back(figure);
                }
            }
            createLightZBuffer(figs, wlight);
        }
    }

    // FIGURES
    Matrix m = Figure3D::eyePointTrans(eye);

    for (unsigned int n = 0; n < allFigures.size(); n++)
    {
        std::string name = "Figure" + std::to_string(n);

        for (auto &figure:allFigures[n])
        {
            figure.applyTransformations(m);

            // colors
            if ((conf["General"]["type"].as_string_or_die().substr(0, 7) == "Lighted"))
            {
                figure.readLights(name, conf);
            }

            else {
                Light tempLight;
                tempLight.ambientLight.ini(conf[name]["color"].as_double_tuple_or_die());
                figure.color = tempLight.ambientLight;
                tempLight.amLight = true;
                figure.ambientReflection.ini({1, 1, 1});
                figure.lights = {tempLight};
            }


            for (Vector3D &point:figure.points)
            {
                figure.doProjection(point, 1);
            }

            if (!zBuffTriangle)
            {
                figure.createLinesOutOfFaces();

                figure.addLines2D(lines);
            }
            else figure.createTriangles();

            figures.emplace_back(figure);
        }
    }

    std::cout << "\n";

    if(zBuffTriangle) return drawZBufferedTriangles();

    return drawLines2D(zBuffered);
}

img::EasyImage Wireframe::drawTextureSphere(const ini::Configuration &conf)
{
    // read information from configuration file
    imageSize = conf["General"]["size"].as_int_or_die();
    nrOfFigures = conf["General"]["nrFigures"].as_int_or_die();
    backgroundcolor.ini(conf["General"]["backgroundcolor"].as_double_tuple_or_die());

    eye = Vector3D::point(conf["General"]["eye"].as_double_tuple_or_die()[0],
                          conf["General"]["eye"].as_double_tuple_or_die()[1],
                          conf["General"]["eye"].as_double_tuple_or_die()[2]);

    for (int k = 0; k < nrOfFigures; k++)
    {
        std::string name = "Figure" + std::to_string(k);

        Figure3D temp(name, conf, false, false);

        Matrix m = Figure3D::eyePointTrans(eye);

        temp.applyTransformations(m);

        temp.createTriangles();

        std::ifstream fin(conf[name]["filename"].as_string_or_die());
        fin >> temp.texture;
        fin.close();

        for (auto &point:temp.points)
        {
            temp.doProjection(point, 1);
        }
        figures.emplace_back(temp);
    }

    double xmin = INT64_MAX;
    double xmax = INT64_MIN;
    double ymin = INT64_MAX;
    double ymax = INT64_MIN;

    for (auto &figure:figures) {

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

    Matrix m = Figure3D::eyePointTrans(eye);

    for (auto &figure:figures)
    {
        for (auto &face:figure.faces)
        {
            image.img::EasyImage::draw_zbuf_triangle_textures(zBuf,
                                                     figure.points[face.pointIndexes[0]],
                                                     figure.points[face.pointIndexes[1]],
                                                     figure.points[face.pointIndexes[2]],
                                                     figure.center*m,
                                                     d, dx, dy, m,
                                                     figure.texture);
        }
    }

    std::cout << "\n";
    return image;
}

// -------------------------------------FRACTALS-----------------------------
void Wireframe::isFractal(std::string name, const ini::Configuration &conf, Figure3D &figure)
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

    allFigures.emplace_back(temp);

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

    allFigures.emplace_back(temp);
}
std::vector<Figure3D> Wireframe::generateThiccFigure(std::string name, const ini::Configuration &conf, Figure3D &root)
{
    std::vector<Figure3D> thiccFigures;

    double r = conf[name]["radius"].as_double_or_die();

    int nC = conf[name]["n"].as_int_or_die();

    int m = conf[name]["m"].as_int_or_die();


    for (int i = 0; i < root.points.size(); ++i)
    {
        Figure3D sphere;
        sphere.createSphere(name, conf, m);

        Matrix m;
        Figure3D::scaleMatrix(m, r);
        Figure3D::translateMatrix(m, Vector3D::vector(root.points[i]));
        sphere.applyTransformations(m);

        thiccFigures.emplace_back(sphere);
    }

    for (const Face &face:root.faces) {

        for (uint index = 0; index < face.pointIndexes.size(); index++) {

            int n = face.pointIndexes[(index + 1) % (face.pointIndexes.size())];

            Vector3D point1 = root.points[face.pointIndexes[index]];
            Vector3D point2 = root.points[n];
            
            double height = sqrt
                    (std::pow(point2.x - point1.x, 2)+
                     std::pow(point2.y - point1.y, 2)+
                     std::pow(point2.z - point1.z, 2));

            Figure3D cylinder;
            cylinder.createCylinder(name, conf, height*1/r, nC);

            Vector3D p1p2 = Vector3D::vector(point2 - point1);
            p1p2.normalise();

            double theta; double phi; double rP;
            cylinder.toPolar(p1p2, theta, phi, rP);

            Matrix m;
            Figure3D::scaleMatrix(m, r);
            Figure3D::rotateAroundY(m, phi);
            Figure3D::rotateAroundZ(m, theta);
            Figure3D::translateMatrix(m, Vector3D::vector(point1));

            cylinder.applyTransformations(m);
            thiccFigures.emplace_back(cylinder);
        }
    }

    return thiccFigures;
}