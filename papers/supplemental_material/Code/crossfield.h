/*
Copyright or © or Copr. Emmanuel Iarussi, David Bommes, Adrien Bousseau 
BendFields: Regularized Curvature Fields from Rough Concept Sketches (2015)

emmanueliarussi (AT) gmail (DOT) com
bommes (AT) aices (DOT) rwth-aachen (DOT) de
adrien (DOT) bousseau (AT) inria (DOT) fr

This software is a computer program whose purpose is to compute cross
fields over sketches using the approach especified in BendFields paper.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef CROSSFIELD_H
#define CROSSFIELD_H

// Includes from the project
#include "tangentmap.h"
#include "unknownsindexer.h"
#include "vec2doperations.h"
#include "labeledmap.h"
#include "periodjumpfield.h"

// External libraries / headers (Solvers, IO, Debugging)
#include <cv.h>
#include <highgui.h>
#include <QDebug>
#include <QList>
#include <QPoint>
#include <QRgb>
#include <QColor>
#include <QDomDocument>
#include <QDomElement>
#include <QFile>
#include "math.h"

using namespace cv;
class UnknownsIndexer;
class PeriodJumpField;

// This class handles crossfield representations
// both exposed in the paper (angle and vector
// based)
class CrossField
{
    // Tangent field
    double ** field;

    // Angle (alpha,beta) representation
    double ** alpha;
    double ** beta;

    // Vector (v0,v1) representation
    Vec2d ** v0_Field;
    Vec2d ** v1_Field;

    // Maps of the important areas in the Crossfield
    // Cells with contains a curvature line information from the sketch
    // Provided by the user manually in the current implementation
    bool  ** curvatureLines;
    // Cells with contains a stroke (but not curvature. Ex: outlines)
    bool ** boundaries;
    // Cells with contains a stroke (thresholded for removing invisible strokes)
    bool  ** constrainedLines;
    // Weight on each constrained cell (depending on the stroke intensity)
    double  ** strokeWeight;
    // For keeping track of the orientation flips
    bool  ** constraintFliped;

    // Size
    int h,w;

    // Create all arrays needed to store the crossfield
    void createField(int,int);

    // No corners image
    Mat noCorners;

public:

    // Constructor
    CrossField();

    // Memory disposer
    ~CrossField();

    // Constructor (just dimensions needed)
    CrossField(int,int);

    // Init the alpha/beta field, given a single angle (called thita)
    void initialice(LabeledMap *, Mat,TangentMap * tangents);

    void initialiceThita( TangentMap * tangents );

    // Getters and Setters for both representations
    inline double getThita(int i,int j) { return this->field[i][j]; }
    inline double getThita(QPoint p) { return this->field[p.x()][p.y()]; }
    inline double getAlpha(QPoint p) { return this->alpha[p.x()][p.y()]; }
    inline double getBeta(QPoint p) {  return this->beta[p.x()][p.y()]; }
    inline double getAlpha(int i,int j) { return this->alpha[i][j]; }
    inline double getBeta(int i,int j) { return this->beta[i][j]; }

    inline Vec2d getV0(int i, int j) { return this->v0_Field[i][j]; }
    inline Vec2d getV1(int i, int j) { return this->v1_Field[i][j]; }
    inline Vec2d getV2(int i, int j) { return -(this->v0_Field[i][j]); }
    inline Vec2d getV3(int i, int j) { return -(this->v1_Field[i][j]); }

    void setAlphaBetha(QList<QPoint>,double,double);
    void setAlphaBetha(double,double, int i, int j);

    inline int height(){ return this->h; }
    inline int width(){ return this->w; }

    void setV0V1(double u0, double v0, double u1, double v1, int i, int j);
    void setV0(double u0, double v0, int i, int j);
    void setV1(double u1, double v1, int i, int j);

    void setThita(int,int,double);
    void cleanThita();

    void markAsCurvatureLine(int i, int j);
    void markAsConstraintFlipped(int i, int j) { this->constraintFliped[i][j] = true; }
    void unmarkAsBoundary(int i, int j) { this->boundaries[i][j] = false; }
    bool isCurvatureLine(int i, int j) { return this->curvatureLines[i][j]; }
    bool isConstrainedLine(int i, int j) { return this->constrainedLines[i][j]; }
    bool isConstraintFliped(int i, int j) { return this->constraintFliped[i][j]; }
    bool isBoundaryPixel(int i, int j) { return this->boundaries[i][j]; }
    double getStrokeWeight(int i, int j) { return this->strokeWeight[i][j]; }

    // Set constraints
    // It sets the constrained areas given the image and
    // A version without corners (see section 4 of the paper)
    void setConstraintsMap(Mat noCornersImage, Mat sketch);

    void updateConstraintsWithUserInput();

    // For normalizing vector representation
    void normalize();

    // Flipping cross representation for normal computation
    void rotateCrosses(PeriodJumpField * pjumpfield, Mat mask, UnknownsIndexer * index);

    // Checks counterclockwise ordering of the vectors in the cross
    bool checkPositiveAngles(UnknownsIndexer * index);

    // Generate the (v0,v1) field from (alpha,beta)
    void generate_V0V1_field();

    // Export and load imformation from files
    void saveToXML(QString path);
    void saveObj(QString path);
    void loadFromXML(QString path);
    void saveInitMap(QString);

};

#endif // CROSSFIELD_H


