# Loop Subdivision

## Introduction
This application subdivides meshes using the Loop Subdivision algorithm. This is used to improve the quality and detail of the mesh.

## Compilation
This application was made using Linux and QT 5.15.3, which can be downloaded at  
https://www.qt.io/download-dev  

If using a Windows machine, the Linux terminal can be accessed by using WSL
Learn more here: https://learn.microsoft.com/en-us/windows/wsl/


To compile the program, enter the following commands in a terminal:  
    
    qmake -project QT+=opengl LIBS+=-lGLU
    qmake
    make

## Usage
Run the program using
`./Loop-Subdivision filename`

### Arguments
`filename` 
- The input mesh to be read and subdivided
- Accepts `.tri` files and `.diredgenormal` files, the latter being a custom file format that stores the mesh as a half edge structure (see the `diredgenormals/` directory for examples)

 **Note** - this algorithm relies on the half edge structure, therefore the input file must be 2-manifold. If it is not, the file will be rejected from the application with one of the following messages:

    Error: this mesh is not manifold since edge X has no other half
    
    Error: this mesh is not manifold since there is not a single cycle for vertex Y

where *X* represents the edge number and *Y* represents the vertex number.

### Interface
The interface contains settings to change how the mesh is viewed including:
- An arcball to change light direction
- An arcball to rotate the model
- Check boxes to visualize vertices and to use flat normals
- A slider that changes vertex size (when visualized)
- X and Y axis sliders to move the model 
- A slider that changes zoom

![Image of Application](/assets/loop%20subdivision.PNG)

### Subdividing
Subdivide the mesh by pressing *Next Subdivision Level*. The result will be rendered upon completion of the subdivision algorithm. The mesh can be subdivided further by simply pressing the button again.

Pressing *Previous Subdivision Level* removes the subdivision and returns the mesh to its previous state.

### Exporting
The mesh currently being rendered can be exported to a file at any point by pressing the *Export File* button. A combo box above this button determines which file type the mesh is to be exported to. The possible file types are:
- `.tri` - File format storing faces and vertices as "triangle soup"
- `.diredgenormal` - Custom file format storing meshes as half edge structures
- `.obj` - Widely recognized object file format

The exported file will be called `surface`, and will be in whichever file format is chosen.

