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

#ifndef BENDFIELD_H
#define BENDFIELD_H

// Includes from the project
#include "crossfield.h"
#include "periodjumpfield.h"
#include "glwidget.h"
#include "unknownsindexer.h"

// External libraries / headers (Solvers, IO, Debugging)
#include <QDebug>
#include <eigen/Eigen/Eigen>
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>
#include <eigen/Eigen/SparseQR>
#include <eigen/Eigen/IterativeLinearSolvers>
#include <eigen/Eigen/Dense>
#include <vector>
#include <cmath>
#include "math.h"

using namespace cv;
using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;
typedef SparseMatrix<double, Eigen::ColMajor> SpMat;

// This class takes a initialized period jump field
// and a crossfield (with some sparse constraints defined)
// and returns an interpolated crossfield with the BendField
// energy.

class BendField
{
    // Crossfield and topology
    CrossField * crossfield;
    PeriodJumpField * pjumpfield;

    // Unknowns index
    UnknownsIndexer * index;

    // Mask of the selected region to solve for
    Mat mask;

    // Sizes
    int h,w,sizeX,numerOfEntries;

    // Build system from the bendfield energy (eq 6 in the paper)
    void buildCovariant(SpMat &A, VectorXd& b, int sizeX);

    // Builds a Harmonic energy for the new representation (u,v)
    // instead of the used in the previous step (alpha,beta)
    // It's only called once, to initialice the field
    // It relies in the computed P
    void buildHarmonic(SpMat &A, VectorXd& b, int sizeX);

    // Compute error
    void computeError(SpMat A, VectorXd b,VectorXd x,int sizeX);

    // Saves the harmonic solution of the system stored in X
    // to the crossfield representation
    void saveIntoField(VectorXd x);

    // Maps the 2D crossfield to X 1D array for u0 variables
    int getColInAfor_u0(int i,int j);

    // Maps the 2D crossfield to X 1D array for v0 variables
    int getColInAfor_v0(int i,int j);

    // Maps the 2D crossfield to X 1D array for u1 variables
    int getColInAfor_u1(int i,int j);

    // Maps the 2D crossfield to X 1D array for v1 variables
    int getColInAfor_v1(int i,int j);

    // Maps the X 1D array to 2D crossfield
    QPoint getPixelInField(int row);


public:

    // Constructor
    BendField();

    // Constructor
    // Init all the vectors and prepares the system for solving the bendfield energy
    // Get's :
    // Crossfield
    // Period Jumps (filled, no modification in the bendfield optimization step)
    // Unknown Index (indexing constrained cels in the crosfield)
    // Mask (the mask for selected region to solve for)
    BendField(CrossField * ,PeriodJumpField *, UnknownsIndexer *, Mat);

    // Iterative method for solving the BendField
    void smoothBendField(GLWidget * glwidget);

};

#endif // BENDFIELD_H
