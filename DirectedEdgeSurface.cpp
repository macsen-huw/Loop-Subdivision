///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face " << face << " ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << " " << otherHalf[dirEdge] << std::endl;
    } // WriteObjectStream()

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()     

DirectedEdgeSurface DirectedEdgeSurface::loopSubdivision()
{
    //Initialise a new DirectedEdgeSurface, storing the mesh after the subdivision
    DirectedEdgeSurface next;

    //Include the previous faceVertices, so we can add the centre face to the correct place
    next.faceVertices = faceVertices;

    //First of all, update the positions of the new vertices using the previous mesh
    //This also preserves the indices
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        next.vertices.push_back(updateVertexPosition(i, getNeighbours(i)));
    }

    //Loop through faces to calculate new edge vertices
    for (unsigned int i = 0; i < faceVertices.size(); i+=3)
    {
        int faceVertex0 = faceVertices[i];
        int faceVertex1 = faceVertices[i+1];
        int faceVertex2 = faceVertices[i+2];

        //Calculate the edge vertex for each edge
        Cartesian3 edgeVertex0 = calculateEdgeVertex(faceVertex2, faceVertex0, faceVertex1, i);
        Cartesian3 edgeVertex1 = calculateEdgeVertex(faceVertex0, faceVertex1, faceVertex2, i+1);
        Cartesian3 edgeVertex2 = calculateEdgeVertex(faceVertex1, faceVertex2, faceVertex0, i+2);

        bool edgeVertex0Exists = false;
        bool edgeVertex1Exists = false;
        bool edgeVertex2Exists = false;

        int edgeVertex0ID, edgeVertex1ID, edgeVertex2ID;

        //Check whether these vertices are already in the vertex array
        //If they are present, assign them their correct IDs
        for (int i = 0; i < next.vertices.size(); i++)
        {
            if (next.vertices[i] == edgeVertex0)
            {
                edgeVertex0Exists = true;
                edgeVertex0ID = i;
            }


            if (next.vertices[i] == edgeVertex1)
            {
                  edgeVertex1Exists = true;
                  edgeVertex1ID = i;
            }

            if(next.vertices[i] == edgeVertex2)
            {
                edgeVertex2Exists = true;
                edgeVertex2ID = i;
            }

        }

        //Include these new vertices in the new vertex vector (if they aren't already included)
        if(!edgeVertex0Exists)
        {
            edgeVertex0ID = next.vertices.size();
            next.vertices.push_back(edgeVertex0);
        }

        if(!edgeVertex1Exists)
        {
            edgeVertex1ID = next.vertices.size();
            next.vertices.push_back(edgeVertex1);
        }

        if(!edgeVertex2Exists)
        {
            edgeVertex2ID = next.vertices.size();
            next.vertices.push_back(edgeVertex2);
        }


        //We now have 4 new faces
        //Update the current face to be the new centre face
        next.faceVertices[i] = edgeVertex0ID;
        next.faceVertices[i+1] = edgeVertex1ID;
        next.faceVertices[i+2] = edgeVertex2ID;

        //Add the other 3 faces at the end of the vector
        next.faceVertices.push_back(faceVertex0);
        next.faceVertices.push_back(edgeVertex1ID);
        next.faceVertices.push_back(edgeVertex0ID);

        next.faceVertices.push_back(edgeVertex1ID);
        next.faceVertices.push_back(faceVertex1);
        next.faceVertices.push_back(edgeVertex2ID);

        next.faceVertices.push_back(edgeVertex0ID);
        next.faceVertices.push_back(edgeVertex2ID);
        next.faceVertices.push_back(faceVertex2);
    }

    //We now have updated vertex and face data
    //Set all FDE values to -1 (a value that can't be a vertex ID) and ensure that there are enough for each vertex
    for(int i = 0; i < next.vertices.size(); i++)
    {
        next.firstDirectedEdge.push_back(-1);
    }
    //Set all otherHalf values to -1 and ensure there are enough to fill all edges
    for(int i = 0; i < next.faceVertices.size(); i++)
        next.otherHalf.push_back(-1);

    //Iterate through the new face vertices
    //Whenever an edge reaches a certain vertex, we update the value of the FDE
    for(unsigned i = 0; i < next.faceVertices.size(); i++)
    {
        //We set the edge to be i+1 (since edges are showing the vertex they point to
        if(next.firstDirectedEdge[next.faceVertices[i]] == -1)
        {
            //Find out what face we're in
            int faceIndex = i / 3;
            int faceLocalIndex = i % 3;
            next.firstDirectedEdge[next.faceVertices[i]] = 3*faceIndex + ((faceLocalIndex + 1) % 3);

        }

        //Now find the edge's other half
        //First. we need the origin and destination of the current edge
        int edgeFace = i / 3;
        int edgeFaceIndex = i % 3;
        int vertexFrom = next.faceVertices[3*edgeFace + ((edgeFaceIndex + 2) % 3)];
        int vertexTo = next.faceVertices[3*edgeFace + (edgeFaceIndex % 3)];

        //Iterate per face, and find the other half
        for(int j = 0; j < next.faceVertices.size(); j+=3)
        {

            if(next.faceVertices[j] == vertexTo && next.faceVertices[j+1] == vertexFrom)
            {
                next.otherHalf[i] = j+1;
                break;
            }

            if(next.faceVertices[j+1] == vertexTo && next.faceVertices[j+2] == vertexFrom)
            {
                next.otherHalf[i] = j+2;
                break;
            }

            if(next.faceVertices[j+2] == vertexTo && next.faceVertices[j] == vertexFrom)
            {
                next.otherHalf[i] = j;
                break;
            }
        }

    }

    //Code to show the mesh on screen
    //Compute the normals
    for(int vertex = 0; vertex < next.vertices.size(); vertex++)
    {

        int firstEdge = next.firstDirectedEdge[vertex];
        int currentEdge = firstEdge;

        Cartesian3 normal = {0,0,0};

        do
        {

            //Given an edge number, we calculate face and face index
            int edgeFace = currentEdge / 3;
            int edgeFaceIndex = currentEdge % 3;

            //Get the three vertices in the face and use to calculate normals
            int vertex0 = next.faceVertices[3*edgeFace + ((edgeFaceIndex + 2) % 3)];
            int vertex1 = next.faceVertices[3*edgeFace + ((edgeFaceIndex) % 3)];
            int vertex2 = next.faceVertices[3*edgeFace + ((edgeFaceIndex + 1) % 3)];

            //Calculate the cross product and add to the normal value
            Cartesian3 edge1 = next.vertices[vertex1] - next.vertices[vertex0];
            Cartesian3 edge2 = next.vertices[vertex2] - next.vertices[vertex0];
            normal = normal + edge1.cross(edge2).unit();

            //Now, we want the next face -> we do this by getting the other half of the current edge
            //Get the other half of the edge
            int other = next.otherHalf[currentEdge];

            //Get face and face index for other half edge
            edgeFace = other / 3;
            edgeFaceIndex = other % 3;

            //Go to the next edge
            currentEdge = (edgeFace * 3) + ((edgeFaceIndex + 1) % 3);

        //Keep repeating until we're back at the first edge (completed the cycle)
        } while (currentEdge != firstEdge);

        //Take the average of the normals
        next.normals.push_back(normal.unit());
    }

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    next.centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (next.vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < next.vertices.size(); vertex++)
            next.centreOfGravity = next.centreOfGravity + next.vertices[vertex];

        // and divide through by the number to get the average position
        // also known as the barycentre
        next.centreOfGravity = next.centreOfGravity / next.vertices.size();

        // start with 0 radius
        next.objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < next.vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (next.vertices[vertex] - next.centreOfGravity).length();

            // now test for maximality
            if (distance > next.objectSize)
                next.objectSize = distance;
            } // per vertex
        } // non-empty vertex set


    return next;
}

std::vector<int> DirectedEdgeSurface::getNeighbours(int vertexID)
{
    std::vector<int> neighbours;
    //We want to walk around the vertex in order to get its neighbours
    //Start with the first face vertex
    int firstEdge = firstDirectedEdge[vertexID];
    int currentEdge = firstEdge;
    do
    {

        //Given an edge number, we calculate face and face index
        int edgeFace = currentEdge / 3;
        int edgeFaceIndex = currentEdge % 3;

        //Given face and face index, we can get the vertex and add to the vector
        neighbours.push_back(faceVertices[3*edgeFace + edgeFaceIndex]);

        //Get the other half of the edge
        int other = otherHalf[currentEdge];

        //Get face and face index for other half edge
        edgeFace = other / 3;
        edgeFaceIndex = other % 3;

        //Go to the next edge
        currentEdge = (edgeFace * 3) + ((edgeFaceIndex + 1) % 3);

    //Keep repeating until we're back at the first edge (completed the cycle)
    } while (currentEdge != firstEdge);


    return neighbours;
}

Cartesian3 DirectedEdgeSurface::calculateEdgeVertex(int vertex1, int vertex2, int vertex3, int edgeNumber)
{
    /*

    Edge vertex calculation requires 4 vertices
    The 2 vertices that connect the edge, and the other 2 vertices in the faces (i.e. the vertices in the face not connected to the edge)

    We know the first 3 edges from the faceVertices, but the fourth vertex requires additional calculation:
    We want the vertex that's in the other face that the edge is connected to
    To do this, we have to find the edge's other half and use that to discover the vertex
    */

    int otherEdge = otherHalf[edgeNumber];
    int otherEdgeFace = otherEdge / 3;
    int otherEdgeFaceNumber = otherEdge % 3;

    int vertex4;
    //We have face of the other half edge, now we just want to find the other vertex (i.e. the one that isn't already in the calculation)
    int otherVertex = faceVertices[3*otherEdgeFace];

    if(otherVertex != vertex1 && otherVertex != vertex2)
        vertex4 = otherVertex;

    otherVertex = faceVertices[3*otherEdgeFace + 1];
    if(otherVertex != vertex1 && otherVertex != vertex2)
        vertex4 = otherVertex;

    otherVertex = faceVertices[3*otherEdgeFace + 2];
    if(otherVertex != vertex1 && otherVertex != vertex2)
        vertex4 = otherVertex;

    //Get the positions of all vertices
    Cartesian3 vertex1Position = vertices[vertex1];
    Cartesian3 vertex2Position = vertices[vertex2];
    Cartesian3 vertex3Position = vertices[vertex3];
    Cartesian3 vertex4Position = vertices[vertex4];

    //Calcuate the new edge vertex based on the positions of the 4 vertices
    Cartesian3 edgeVertex = ((3.0/8.0) * (vertex1Position + vertex2Position)) + ((1.0/8.0) * (vertex3Position + vertex4Position));

    return edgeVertex;

}

Cartesian3 DirectedEdgeSurface::updateVertexPosition(int vertex, std::vector<int> vertexNeighbours)
{
    float n = vertexNeighbours.size();
    float beta = 0;

    if (n == 3)
        beta = 9.0/16.0;
    else if (n > 3)
        beta = 5.0/8.0 - std::pow((3.0/8.0 + (1.0/4.0 * std::cos(2*M_PI / n))), 2);

    //Find the sum of all neighbouring vertices
    Cartesian3 sumofNeighbours = {0,0,0};
    for(int neighbour: vertexNeighbours)
        sumofNeighbours = sumofNeighbours + vertices[neighbour];

    //Divide this total by n
    sumofNeighbours = sumofNeighbours / n;

    //Finally, calculate the updated position of the vertex
    Cartesian3 updatedVertex = ((1 - beta) * vertices[vertex]) + (beta * sumofNeighbours);

    return updatedVertex;
}

