#ifndef TRICONVERTER_H
#define TRICONVERTER_H

#include <fstream>
#include <iostream>

#include "DirectedEdgeSurface.h"
#include "Cartesian3.h"

class TriConverter
{
public:
    TriConverter();
    DirectedEdgeSurface convertTriFile(std::istream &in);

private:
    //Values required for Directed Edge Structures
    std::vector<Cartesian3> vertices;
    std::vector<Cartesian3> normals;
    std::vector<int> degrees;
    std::vector<unsigned int> faceVertices;
    std::vector<unsigned int> firstDirectedEdge;
    std::vector<unsigned int> otherHalf;

    //Get info from the tri file - read vertices and faces
    void readTriFile();
    void getFDEandOtherHalves();
    void addFDE(int index);
    int addVertex(Cartesian3 coords);

    int getOtherHalf(int vertexFrom, int vertexTo);

    void calculateNormals();
    void testManifold();

    Cartesian3 getCentreOfGravity();
    float getObjectSize(Cartesian3 centreOfGravity);
};
#endif // TRICONVERTER_H
