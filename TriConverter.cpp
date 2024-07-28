#include "TriConverter.h"

TriConverter::TriConverter()
{
}

DirectedEdgeSurface TriConverter::convertTriFile(std::istream &in)
{
    //Read the .tri file

   //Get the number of faces, in other words the number of times to loop
    int numFaces;
    in >> numFaces;

    //Get the vertices and face vertices
    for(int i = 0; i < numFaces * 3; i++)
    {
        float X, Y, Z;
        in >> X; in >> Y; in >> Z;

        int vertIndex = addVertex(Cartesian3{X,Y,Z});
        faceVertices.emplace_back(vertIndex);
    }

    //Next, get the first directed edges
    getFDEandOtherHalves();

    //Before we continue, we must check whether the mesh is manifold, otherwise loop subdivision won't work
    testManifold();

    //If we reach this point, we've passed the test
    //Next, calculate normals
    calculateNormals();

    //Calculate centre of gravity
    Cartesian3 centreOfGravity = getCentreOfGravity();

    float size = getObjectSize(centreOfGravity);

    //We've now calculate everything needed, so return the DirectedEdgeSurface
    DirectedEdgeSurface convertedTri;
    convertedTri.vertices = vertices;
    convertedTri.faceVertices = faceVertices;
    convertedTri.firstDirectedEdge = firstDirectedEdge;
    convertedTri.otherHalf = otherHalf;
    convertedTri.normals = normals;
    convertedTri.centreOfGravity = centreOfGravity;
    convertedTri.objectSize = size;

    return convertedTri;
}

//Returns the vertex's index
int TriConverter::addVertex(Cartesian3 coords)
{
    //Check whether vertex is already in the vector
    for(int i = 0; i < vertices.size(); i++)
    {
        if(vertices[i] == coords)
        {
            degrees[i]++;
            return i;
        }

    }

    //Otherwise, add new vertex
    int index = vertices.size();
    vertices.emplace_back(coords);
    degrees.emplace_back(1);

    return index;
}

//Calculate the list of first directed edges and other halves
void TriConverter::getFDEandOtherHalves()
{
    //First we initialise the firstDirectedEdges and otherHalves to have -1 values
    firstDirectedEdge.assign(vertices.size(), -1);
    otherHalf.assign(faceVertices.size(), -1);

    //Now, loop through each face to ensure that each vertex has an FDE
    //Then, find the other half of each edge
    for(int i = 0; i < faceVertices.size(); i+= 3)
    {
        //Order of edges on a face
        // 2->0, 0->1, 1->2

        int vert0Index = faceVertices[i];
        int vert1Index = faceVertices[i+1];
        int vert2Index = faceVertices[i+2];

        //First of all, check for FDE
        if(firstDirectedEdge[vert0Index] == -1)
            addFDE(i);
        if(firstDirectedEdge[vert1Index] == -1)
            addFDE(i+1);
        if(firstDirectedEdge[vert2Index] == -1)
            addFDE(i+2);

        //Now, get other halves
        otherHalf[i] = getOtherHalf(vert0Index, vert2Index);
        otherHalf[i+1] = getOtherHalf(vert1Index, vert0Index);
        otherHalf[i+2] = getOtherHalf(vert2Index, vert1Index);
    }

}

void TriConverter::addFDE(int index)
{
    int faceNumber = index / 3;
    int facePosition = index % 3;

    int edgeNumber = faceNumber*3 + ((facePosition + 1) % 3);
    firstDirectedEdge[faceVertices[index]] = edgeNumber;
}

int TriConverter::getOtherHalf(int vertexFrom, int vertexTo)
{
    for(int i = 0; i < faceVertices.size(); i+= 3)
    {
        int vertexA = faceVertices[i];
        int vertexB = faceVertices[i+1];
        int vertexC = faceVertices[i+2];

        //Check for the edge that goes from and to our desired vertices
        if(vertexC == vertexFrom && vertexA == vertexTo)
            return i;
        if(vertexA == vertexFrom && vertexB == vertexTo)
            return i+1;
        if(vertexB == vertexFrom && vertexC == vertexTo)
            return i+2;
    }

    //If not found, then return -1
    return -1;
}


//2 conditions to be manifold:
// 1 - All edges have 2 faces
// 2 - Each vertex has a single cycle
void TriConverter::testManifold()
{
    //First test - ensure that all edges have other halves
    for(int i = 0; i < otherHalf.size(); i++)
    {
        if(otherHalf[i] == -1)
        {
            std::cout << "Error: this mesh is not manifold since edge " << i << " has no other half" << std::endl;
            exit(0);
        }
    }

    //Second test - check whether each vertex has a single cycle
    for(int i = 0; i < vertices.size(); i++)
    {
        //We start with the FDE of the vertex and walk our way around
        int firstEdge = firstDirectedEdge[i];
        int currentEdge = firstEdge;
        int cycleLength = 0;
        bool endCycle = false;

        while (endCycle == false)
        {
            //Get the other half
            int other = otherHalf[currentEdge];

            //Find out which face it's in to discover the next edge
            int face = other / 3;
            int index = other % 3;

            int nextEdge = face*3 + ((index + 1) % 3);
            cycleLength ++;

            //Have we finished the cycle?
            if(nextEdge == firstEdge)
            {
                endCycle = true;
                continue;
            }

            //Otherwise, the next edge becomes the current edge
            currentEdge = nextEdge;
        }


        //Have we walked through all neighbours?
        if(cycleLength != degrees[i])
        {
            std::cout << "Error: this mesh is not manifold since there is not a single cycle for vertex " << i << std::endl;
            exit(0);
        }
    }
}

void TriConverter::calculateNormals()
{

    //We want the normals for each vertex
    for(int i = 0; i < vertices.size(); i++)
    {

        Cartesian3 normal = {0,0,0};

        //To do this, we need to sum the face normals
        //We get the faces by walking around the vertex
        int firstEdge = firstDirectedEdge[i];
        int currentEdge = firstEdge;

        do
        {
            //Get the face containing the FDE
            int edgeFace = currentEdge / 3;
            int edgeIndex = currentEdge % 3;

            //Get the 3 vertices of the face
            int vertex0 = faceVertices[3*edgeFace + ((edgeIndex + 2) % 3)];
            int vertex1 = faceVertices[3*edgeFace + edgeIndex];
            int vertex2 = faceVertices[3*edgeFace + ((edgeIndex + 1) % 3)];

            //Calculate edges
            Cartesian3 edge1 = vertices[vertex1] - vertices[vertex0];
            Cartesian3 edge2 = vertices[vertex2] - vertices[vertex0];

            //Accumulate normal value
            normal = normal + edge1.cross(edge2).unit();

            //Get the other half
            int other = otherHalf[currentEdge];

            //Find out which face it's in to discover the next edge
            int face = other / 3;
            int index = other % 3;

            currentEdge = face*3 + ((index + 1) % 3);
        } while (currentEdge != firstEdge);

       //Normalize the vertex normals
        normals.emplace_back(normal.unit());
    }
}

Cartesian3 TriConverter::getCentreOfGravity()
{
    Cartesian3 centreOfGravity = {0.0, 0.0, 0.0};

   for (int i = 0; i < vertices.size(); i++)
   {
       centreOfGravity = centreOfGravity + vertices[i];
   }

   //Divide through by the number to get the average position
   centreOfGravity = centreOfGravity / vertices.size();

   return centreOfGravity;

}

float TriConverter::getObjectSize(Cartesian3 centreOfGravity)
{
    float size = 0;

    for(int i = 0; i < vertices.size(); i++)
    {
        float distance = (vertices[i] - centreOfGravity).length();

        if(distance > size)
            size = distance;
    }

    return size;
}
