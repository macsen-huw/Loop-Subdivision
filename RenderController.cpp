/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Controller
//  -----------------------------
//  
//  We're using the Model-View-Controller pattern
//  so most of the control logic goes here
//  which means we need a slot for substantially
//  every possible UI manipulation
//
/////////////////////////////////////////////////////////////////

#include "RenderController.h"
#include <stdio.h>

// constructor
RenderController::RenderController
        (
        // the geometric object to show
        std::vector<DirectedEdgeSurface>      *newDirectedEdgeSurface,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // the render window that it controls
        RenderWindow        *newRenderWindow
        )
    :
    directedEdgeSurface  (newDirectedEdgeSurface),
    renderParameters(newRenderParameters),
    renderWindow    (newRenderWindow),
    dragButton      (Qt::NoButton),
    startx(-1), starty(-1),
    vertexID(-1)
    { // RenderController::RenderController()
    
    // connect up signals to slots

    // signals for arcballs
    QObject::connect(   renderWindow->modelRotator,                 SIGNAL(RotationChanged()),
                        this,                                       SLOT(objectRotationChanged()));
    QObject::connect(   renderWindow->lightRotator,                 SIGNAL(RotationChanged()),
                        this,                                       SLOT(lightRotationChanged()));

    // signals for main widget to control arcball
    QObject::connect(   renderWindow->renderWidget,                 SIGNAL(MouseDown(int, int, int)),
                        this,                                       SLOT(MouseDown(int, int, int)));
    QObject::connect(   renderWindow->renderWidget,                 SIGNAL(MouseDrag(int, int)),
                        this,                                       SLOT(MouseDrag(int, int)));
    QObject::connect(   renderWindow->renderWidget,                 SIGNAL(MouseUp(int, int)),
                        this,                                       SLOT(MouseUp(int, int)));

    // signal for zoom slider
    QObject::connect(   renderWindow->zoomSlider,                   SIGNAL(valueChanged(int)),
                        this,                                       SLOT(zoomChanged(int)));

    // signal for x translate sliders
    QObject::connect(   renderWindow->xTranslateSlider,             SIGNAL(valueChanged(int)),
                        this,                                       SLOT(xTranslateChanged(int)));

    // signal for y translate slider
    QObject::connect(   renderWindow->yTranslateSlider,             SIGNAL(valueChanged(int)),
                        this,                                       SLOT(yTranslateChanged(int)));

    // signal for vertex scale slider
    QObject::connect(   renderWindow->vertexSizeSlider,	            SIGNAL(valueChanged(int)),
                        this,                                       SLOT(vertexSizeChanged(int)));

    // signal for check box for flat normals
    QObject::connect(   renderWindow->flatNormalsBox,          		SIGNAL(stateChanged(int)),
                        this,                                       SLOT(flatNormalsCheckChanged(int)));

    // signal for check box for showing vertices
    QObject::connect(   renderWindow->showVerticesBox, 		        SIGNAL(stateChanged(int)),
                        this,                                       SLOT(showVerticesCheckChanged(int)));

    QObject::connect(   renderWindow->increaseSubdivisionButton,    SIGNAL(released()),
                        this,                                       SLOT(subdivisionIncreased()));

    QObject::connect(   renderWindow->decreaseSubdivisionButton, 	SIGNAL(released()),
                        this,                                       SLOT(subdivisionDecreased()));

    QObject::connect( renderWindow->exportFileType, SIGNAL(currentIndexChanged(int)),
                      this,                         SLOT(fileTypeChanged(int)));

    QObject::connect(renderWindow->exportSurface, SIGNAL(released()),
                     this,                        SLOT(exportCheckChanged()));

    // copy the rotation matrix from the widgets to the model
    renderParameters->rotationMatrix = renderWindow->modelRotator->RotationMatrix();
    renderParameters->lightMatrix = renderWindow->lightRotator->RotationMatrix();
    } // RenderController::RenderController()

// slot for responding to arcball rotation for object
void RenderController::objectRotationChanged()
    { // RenderController::objectRotationChanged()
    // copy the rotation matrix from the widget to the model
    renderParameters->rotationMatrix = renderWindow->modelRotator->RotationMatrix();
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::objectRotationChanged()

// slot for responding to arcball rotation for light
void RenderController::lightRotationChanged()
    { // RenderController::lightRotationChanged()
    // copy the rotation matrix from the widget to the model
    renderParameters->lightMatrix = renderWindow->lightRotator->RotationMatrix();
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::lightRotationChanged()

// slot for responding to zoom slider
void RenderController::zoomChanged(int value)
    { // RenderController::zoomChanged()
    // compute the new scale
    float newZoomScale = pow(10.0, (float) value / 100.0);

    // clamp it
    if (newZoomScale < ZOOM_SCALE_MIN)
        newZoomScale = ZOOM_SCALE_MIN;
    else if (newZoomScale > ZOOM_SCALE_MAX)
        newZoomScale = ZOOM_SCALE_MAX;

    // and reset the value  
    renderParameters->zoomScale = newZoomScale;
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::zoomChanged()

// slot for responding to x translate sliders
void RenderController::xTranslateChanged(int value)
    { // RenderController::xTranslateChanged()
    // reset the model's x translation (slider ticks are 1/100 each)
    renderParameters->xTranslate = (float) value / 100.0;

    // clamp it
    if (renderParameters->xTranslate < TRANSLATE_MIN)
        renderParameters->xTranslate = TRANSLATE_MIN;
    else if (renderParameters->xTranslate > TRANSLATE_MAX)
        renderParameters->xTranslate = TRANSLATE_MAX;
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::xTranslateChanged()

// slot for responding to y translate slider
void RenderController::yTranslateChanged(int value)
    { // RenderController::tTranslateChanged()
    // reset the model's y translation (slider ticks are 1/100 each)
    renderParameters->yTranslate =  (float) value / 100.0;

    // clamp it
    if (renderParameters->yTranslate < TRANSLATE_MIN)
        renderParameters->yTranslate = TRANSLATE_MIN;
    else if (renderParameters->yTranslate > TRANSLATE_MAX)
        renderParameters->yTranslate = TRANSLATE_MAX;
    
    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::yTranslateChanged()

// slot for responding to vertex size slider
void RenderController::vertexSizeChanged(int value)
    { // RenderController::vertexSizeChanged()
    // reset the model's y translation (slider ticks are 1/512 each)
    renderParameters->vertexSize = (float) value / 512;

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::vertexSizeChanged()
    
// slot for toggling whether vertices are shown
void RenderController::showVerticesCheckChanged(int state)
    { // RenderController::showVerticesCheckChanged()
    // reset the model's flag
    renderParameters->showVertices = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::showVerticesCheckChanged()
     
// slot for toggling flat normals
void RenderController::flatNormalsCheckChanged(int state)
    { // RenderController::flatNormalsCheckChanged()
    // reset the model's flag
    renderParameters->useFlatNormals = (state == Qt::Checked); 

    // reset the interface
    renderWindow->ResetInterface();
    } // RenderController::flatNormalsCheckChanged()

//Slot for increasing subdivision level
void RenderController::subdivisionIncreased()
{
    //Increase subdivisionLevel
    renderParameters->subdivisionLevel++;

    //If the subdivision level is greater than the number of meshes in the vector, create a new one
    if (directedEdgeSurface->size() < (renderParameters->subdivisionLevel+1))
    {
        //Perform another iteration of loop subdivision and add to the array
        int finalIndex = directedEdgeSurface->size() - 1;
        directedEdgeSurface->push_back(directedEdgeSurface->at(finalIndex).loopSubdivision());
    }

    //Reset Interface
    renderWindow->ResetInterface();
}

//Slot for decreasing subdivision level
void RenderController::subdivisionDecreased()
{
    //Decrease subdivision level (must be > 0)
    if (renderParameters->subdivisionLevel != 0)
    {
        renderParameters->subdivisionLevel--;
    }

    //Reset Interface
    renderWindow->ResetInterface();
}

//Slot for changing the combo box value
void RenderController::fileTypeChanged(int value)
{
    renderParameters->chosenFileType = value;

    renderWindow->ResetInterface();
}

//Slot for changing state of export button
void RenderController::exportCheckChanged()
{
    int type;
    std::string title;
    if(renderParameters->chosenFileType == 2)
    {
        type = OBJ;
        title = "surface.obj";
    }


    else if(renderParameters->chosenFileType == 1)
    {
       type = DIREDGENORMAL;
       title = "surface.diredgenormal";
    }

    else if(renderParameters->chosenFileType == 0)
    {
        type = TRI;
        title = "surface.tri";
    }


    std::ofstream outputFile(title);
    directedEdgeSurface->at(renderParameters->subdivisionLevel).WriteObjectStream(outputFile, type);

    renderWindow->ResetInterface();
}
     
// slots for responding to mouse manipulation
// these are general purpose signals which pass the mouse moves to the controller
// after scaling to the notional unit sphere
void RenderController::MouseDown(int whichButton, int x, int y)
    { // RenderController::MouseDown()
    } // RenderController::MouseDown()
    
// note that Continue & End assume the button has already been set
void RenderController::MouseDrag(int x, int y)
    { // RenderController::MouseDrag()
    } // RenderController::MouseDrag()

void RenderController::MouseUp(int x, int y)
    { // RenderController::MouseUp()
    } // RenderController::MouseUp()
