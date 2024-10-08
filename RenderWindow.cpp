/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Window
//  -----------------------------
//  
//  The render window class is really just a container
//  for tracking the visual hierarchy.  While it is
//  entirely possible to use Qt Creator, I try to avoid
//  over-commitment to it because I need to write code in
//  multiple environments, some of which are not well-suited
//  to IDEs in general, let alone Qt Creator
//
//  Also, by providing sample code, the didactic purpose of 
//  showing how things fit together is better served.
//
/////////////////////////////////////////////////////////////////

#include "RenderWindow.h"
#include "RenderParameters.h"

// constructor
RenderWindow::RenderWindow
        (
        // the object to be rendered
        std::vector<DirectedEdgeSurface>     *newDirectedEdgeSurface,
        // the model object storing render parameters
        RenderParameters        *newRenderParameters,
        // the title for the window (with default value)
        const char              *windowName
        )
    // call the inherited constructor
    // NULL indicates that this widget has no parent
    // i.e. that it is a top-level window
    :
    // member instantiation
    QWidget(NULL),
    directedEdgeSurface(newDirectedEdgeSurface),
    renderParameters(newRenderParameters)
    { // RenderWindow::RenderWindow()
    // set the window's title
    setWindowTitle(QString(windowName));
    
    // initialise the grid layout
    windowLayout = new QGridLayout(this);
    
    // create all of the widgets, starting with the custom render widgets
    renderWidget                = new RenderWidget              (newDirectedEdgeSurface,     newRenderParameters,        this);

    // construct custom arcball Widgets
    lightRotator                = new ArcBallWidget             (                       this);
    modelRotator                = new ArcBallWidget             (                       this);

    // construct standard QT widgets
    // modelling options
	showVerticesBox				= new QCheckBox					("Show Vertices",		this);
	flatNormalsBox				= new QCheckBox					("Flat Normals",		this);
    
    // spatial sliders
    xTranslateSlider            = new QSlider                   (Qt::Horizontal,        this);
    yTranslateSlider            = new QSlider                   (Qt::Vertical,          this);
    zoomSlider                  = new QSlider                   (Qt::Vertical,          this);

	// slider for vertex size
    vertexSizeSlider            = new QSlider                   (Qt::Horizontal,        this);

    //Buttons to increase/decrease subdivision
    increaseSubdivisionButton  = new QPushButton                ("Next Subdivision Level", this);
    decreaseSubdivisionButton  = new QPushButton                 ("Previous Subdivision Level", this);
    exportSurface               = new QPushButton               ("Export File", this);

    //QComboBox for the file format selection
    exportFileType = new QComboBox(this);
    exportFileType->addItems({".tri", ".diredgenormal", ".obj"});

    // labels for sliders and arcballs
    modelRotatorLabel           = new QLabel                    ("Model",               this);
    lightRotatorLabel           = new QLabel                    ("Light",               this);
    yTranslateLabel             = new QLabel                    ("Y",                   this);
    zoomLabel                   = new QLabel                    ("Zm",                  this);
    vertexSizeLabel				= new QLabel					("Vertex Size",			this);
    
    // add all of the widgets to the grid               Row         Column      Row Span    Column Span
    
    // the top two widgets have to fit to the widgets stack between them
    int nStacked = 6;
    
    windowLayout->addWidget(renderWidget,               0,          1,          nStacked,   1           );
    windowLayout->addWidget(yTranslateSlider,           0,          2,          nStacked,   1           );
    windowLayout->addWidget(zoomSlider,                 0,          4,          nStacked,   1           );

    // the stack in the middle
    windowLayout->addWidget(lightRotator,               0,          3,          1,          1           );
    windowLayout->addWidget(lightRotatorLabel,          1,          3,          1,          1           );
    windowLayout->addWidget(modelRotator,               2,          3,          1,          1           );
    windowLayout->addWidget(modelRotatorLabel,          3,          3,          1,          1           );
    windowLayout->addWidget(flatNormalsBox,		        4,          3,          1,          1           );
    windowLayout->addWidget(showVerticesBox,	        5,          3,          1,          1           );

    // Translate Slider Row
    windowLayout->addWidget(xTranslateSlider,           nStacked,   1,          1,          1           );
    windowLayout->addWidget(yTranslateLabel,            nStacked,   2,          1,          1           );
    // nothing in column 3
    windowLayout->addWidget(zoomLabel,                  nStacked,   4,          1,          1           );
    
    // Vertex Size slider row
    windowLayout->addWidget(vertexSizeSlider,			nStacked+1,	1,			1,			1			);
    windowLayout->addWidget(vertexSizeLabel,			nStacked+1,	3,			1,			1			);

    //Add buttons underneath the vertex size sliders
    windowLayout->addWidget(increaseSubdivisionButton,  nStacked+2, 1,          1,          1           );
    windowLayout->addWidget(decreaseSubdivisionButton,  nStacked+3, 1,          1,          1           );

    windowLayout->addWidget(exportFileType, nStacked+2, 3, 1, 1);
    windowLayout->addWidget(exportSurface, nStacked + 3, 3, 1, 1);

    // now reset all of the control elements to match the render parameters passed in
    ResetInterface();
    } // RenderWindow::RenderWindow()

// routine to reset interface
// sets every visual control to match the model
// gets called by the controller after each change in the model
void RenderWindow::ResetInterface()
    { // RenderWindow::ResetInterface()
    // set check boxes
    showVerticesBox    	    ->setChecked        (renderParameters   ->  showVertices);
    flatNormalsBox    	    ->setChecked        (renderParameters   ->  useFlatNormals);

    //Set Combo Box
    exportFileType->setCurrentIndex(renderParameters->chosenFileType);
    
    // set sliders
    // x & y translate are scaled to notional unit sphere in render widgets
    // but because the slider is defined as integer, we multiply by a 100 for all sliders
    xTranslateSlider        ->setMinimum        ((int) (TRANSLATE_MIN                               * PARAMETER_SCALING));
    xTranslateSlider        ->setMaximum        ((int) (TRANSLATE_MAX                               * PARAMETER_SCALING));
    xTranslateSlider        ->setValue          ((int) (renderParameters -> xTranslate              * PARAMETER_SCALING));
    
    yTranslateSlider        ->setMinimum        ((int) (TRANSLATE_MIN                               * PARAMETER_SCALING));
    yTranslateSlider        ->setMaximum        ((int) (TRANSLATE_MAX                               * PARAMETER_SCALING));
    yTranslateSlider        ->setValue          ((int) (renderParameters -> yTranslate              * PARAMETER_SCALING));

    // zoom slider is a logarithmic scale, so we want a narrow range
    zoomSlider              ->setMinimum        ((int) (ZOOM_SCALE_LOG_MIN                          * PARAMETER_SCALING));
    zoomSlider              ->setMaximum        ((int) (ZOOM_SCALE_LOG_MAX                          * PARAMETER_SCALING));
    zoomSlider              ->setValue          ((int) (log10(renderParameters -> zoomScale)        * PARAMETER_SCALING));

	// vertex scale is a linear variable in %
	vertexSizeSlider		->setMinimum		(0);
	vertexSizeSlider		->setMaximum		(512);
	vertexSizeSlider		->setValue			(512 * renderParameters -> vertexSize);

    // now flag them all for update 
    renderWidget            ->update();
    modelRotator            ->update();
    lightRotator            ->update();
    xTranslateSlider        ->update();
    yTranslateSlider        ->update();
    zoomSlider              ->update();
    vertexSizeSlider        ->update();
    showVerticesBox		    ->update();
    flatNormalsBox		    ->update();
    exportFileType->update();
    } // RenderWindow::ResetInterface()
