//
// Created by mano on 29.03.19.
//

#ifndef CGENGINE_ZBUFFER_H
#define CGENGINE_ZBUFFER_H

#include <vector>
#include <limits>

class ZBuffer
{
        std::vector<std::vector<double>> zVals;
    public:
        ZBuffer(const int width, const int height);

        double getZVal(int x, int y);

        void setVal(int x, int y, double val);
};


#endif //CGENGINE_ZBUFFER_H
