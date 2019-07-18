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

#include "bendfield.h"

// Weight on the harmonic energy for
// initialization
#define WEIGHTHARMONICSMOOTH 1.0

// Weight on the bendfield energy for
// regularization (E_h in the paper, eq 6)
#define WEIGHTREGULARIZATION 0.01

// Weight on the bendfield energy
#define WEIGHTBENDFIELD 0.1

// Weight on the bendfield energy for the constraints
// (W_strokes in the paper, eq 6)
// Other values to try 0.04 - 0.03 - 0.025 - 0.01
#define WEIGHTCONSTRAINT 0.025

// Constructor
BendField::BendField() {}

// Constructor
// Init all the vectors and prepares the system for solving the bendfield energy
// Get's :
// Crossfield
// Period Jumps (filled, no modification in the bendfield optimization step)
// Unknown Index (indexing constrained cels in the crosfield)
// Mask (the mask for selected region to solve for)
BendField::BendField(CrossField * c,PeriodJumpField * p, UnknownsIndexer * idx, Mat m)
{
    // Point crossfield and pjumpfield
    this->crossfield = c;
    this->pjumpfield = p;

    // Mask
    this->mask = m;

    // Index
    this->index = idx;

    // Sizes
    this->h = this->crossfield->height();
    this->w = this->crossfield->width();

    // Number of unknowns and pixels
    this->numerOfEntries = idx->getNumberOfPixels();
    this->sizeX = idx->getNumberOfUknowns();
}

// Computes the residual of the system as a measure of convergence
void BendField::computeError(SpMat A, VectorXd b,VectorXd x,int sizeX)
{
    VectorXd error = (A*x)-b;
    double residual = 0;

    for(int i = 0; i < sizeX; i++)
    {
        residual += (error(i)*error(i));
    }

    qDebug() << "Residual " << residual;
}

// Maps the 2D crossfield to X 1D array for u0 variables
int BendField::getColInAfor_u0(int i,int j)
{
    return this->index->getFromIndex(i,j,0);
}

// Maps the 2D crossfield to X 1D array for v0 variables
int BendField::getColInAfor_v0(int i,int j)
{
    return this->index->getFromIndex(i,j,1);
}

// Maps the 2D crossfield to X 1D array for u1 variables
int BendField::getColInAfor_u1(int i,int j)
{
    return this->index->getFromIndex(i,j,2);
}

// Maps the 2D crossfield to X 1D array for v1 variables
int BendField::getColInAfor_v1(int i,int j)
{
    return this->index->getFromIndex(i,j,3);
}

// Maps the X 1D array to 2D crossfield
QPoint BendField::getPixelInField(int row)
{
    return index->getPixelIndexed(row);
}

// Build system from the bendfield energy (eq 6 in the paper)
void BendField::buildCovariant(SpMat &A, VectorXd& b, int sizeX)
{
    // Triplets init
    vector<T>  * tripletList = new vector<T>;
    tripletList->reserve(sizeX*6);
    tripletList->clear();

    // for values in b
    QList<QPair<int,double> > b_values;

    int lastIndexA = 0;

    // Set up the constraints
    // For every pixel
    for(int i = 2; i < this->h-2; i++)
    {
        for(int j = 2; j < this->w-2; j++)
        {
            // White pixel in the mask
            if(index->isUnknown(i,j)==true)
            {
                // Get the coordinates for every variable involved
                // Mapping 2D -> 1D
                // V0 Center
                int u0_c = getColInAfor_u0(i,j);
                int v0_c = getColInAfor_v0(i,j);

                // V0 Right Neighbourd
                int u0_r = getColInAfor_u0(i,j+1);
                int v0_r = getColInAfor_v0(i,j+1);

                // V0 Bottom Neighbourd
                int u0_b = getColInAfor_u0(i+1,j);
                int v0_b = getColInAfor_v0(i+1,j);

                // V1 Center
                int u1_c = getColInAfor_u1(i,j);
                int v1_c = getColInAfor_v1(i,j);

                // V1 Right Neighbourd
                int u1_r = getColInAfor_u1(i,j+1);
                int v1_r = getColInAfor_v1(i,j+1);

                // V1 Bottom Neighbourd
                int u1_b = getColInAfor_u1(i+1,j);
                int v1_b = getColInAfor_v1(i+1,j);

                // Get the values for every variable involved

                // V0 Center
                double val_u0_c = 0;
                double val_v0_c = 0;

                // V1 Center
                double val_u1_c = 0;
                double val_v1_c = 0;

                val_u0_c = this->crossfield->getV0(i,j)(0);
                val_v0_c = this->crossfield->getV0(i,j)(1);

                val_u1_c = this->crossfield->getV1(i,j)(0);
                val_v1_c = this->crossfield->getV1(i,j)(1);

                // Get p's (topology)
                int rightP = (this->pjumpfield->getRightP(i,j))%4;
                int bottomP = (this->pjumpfield->getBottomP(i,j))%4;

                // Handling flips (see Vector representation in paper)
                // To handle flips
                int t1 = 0;
                int t2 = 0;
                int t3 = 0;
                int t4 = 0;
                int t5 = 0;
                int t6 = 0;
                int t7 = 0;
                int t8 = 0;

                // sign to flip V0 or V1
                double t1_sign = 1.0;
                double t2_sign = 1.0;
                double t3_sign = 1.0;
                double t4_sign = 1.0;
                double t5_sign = 1.0;
                double t6_sign = 1.0;
                double t7_sign = 1.0;
                double t8_sign = 1.0;

                // RIGHT
                // V0_C --> V0_R(u0,v0)
                // V1_C --> V1_R(u1,v1)
                if(rightP==0)
                {
                   t1 = u0_r;
                   t3 = v0_r;
                   t5 = u1_r;
                   t7 = v1_r;
                }
                // V0_C --> -V1_R(u1,v1)
                // V1_C --> V0_R(u0,v0)
                else if((rightP==1)||(rightP==-3))
                {
                    t1 = u1_r;
                    t3 = v1_r;
                    t5 = u0_r;
                    t7 = v0_r;

                    t1_sign = -1.0;
                    t3_sign = -1.0;
                }
                // V0_C --> V1_R(u1,v1)
                // V1_C --> -V0_R(u0,v0)
                else if((rightP==-1)||(rightP==3))
                {
                    t1 = u1_r;
                    t3 = v1_r;
                    t5 = u0_r;
                    t7 = v0_r;

                    t5_sign = -1.0;
                    t7_sign = -1.0;
                }
                // V0_C --> -V0_R(u0,v0)
                // V1_C --> -V1_R(u1,v1)
                else if((rightP==2)||(rightP==-2))
                {
                    t1 = u0_r;
                    t3 = v0_r;
                    t5 = u1_r;
                    t7 = v1_r;

                    t1_sign = -1.0;
                    t3_sign = -1.0;
                    t5_sign = -1.0;
                    t7_sign = -1.0;
                }

                // BOTTOM
                // V0 --> V0(u0,v0)
                // V1 --> V1(u1,v1)
                if(bottomP==0)
                {
                    t2 = u0_b;
                    t4 = v0_b;
                    t6 = u1_b;
                    t8 = v1_b;
                }
                // V0 --> -V1(u1,v1)
                // V1 --> V0(u0,v0)
                else if((bottomP==1)||(bottomP==-3))
                {
                    t2 = u1_b;
                    t4 = v1_b;
                    t6 = u0_b;
                    t8 = v0_b;

                    t2_sign = -1.0;
                    t4_sign = -1.0;
                }
                // V0 --> V1(u1,v1)
                // V1 --> -V0(u0,v0)
                else if((bottomP==-1)||(bottomP==3))
                {
                    t2 = u1_b;
                    t4 = v1_b;
                    t6 = u0_b;
                    t8 = v0_b;

                    t6_sign = -1.0;
                    t8_sign = -1.0;
                }
                // V0 --> -V0(u0,v0)
                // V1 --> -V1(u1,v1)
                else if((bottomP==2)||(bottomP==-2))
                {
                    t2 = u0_b;
                    t4 = v0_b;
                    t6 = u1_b;
                    t8 = v1_b;

                    t2_sign = -1.0;
                    t4_sign = -1.0;
                    t6_sign = -1.0;
                    t8_sign = -1.0;
                }

                // ADD BENDFIELD CONSTRAINTS (when possible, checking borders)

                if(((!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                    &&!((crossfield->isBoundaryPixel(i,j+1))&&(!crossfield->isCurvatureLine(i,j+1))))
                    &&!((crossfield->isBoundaryPixel(i+1,j))&&(!crossfield->isCurvatureLine(i+1,j))))
                {

                    // Handle bottom boundary
                    if((this->index->isUnknown(i,j+1)==true)&&(this->index->isUnknown(i+1,j)==true))
                    {
                        // First term
                        // (u0r-u0c) * <u1c> +
                        tripletList->push_back(T(lastIndexA,t1,WEIGHTBENDFIELD * val_u1_c * t1_sign));
                        tripletList->push_back(T(lastIndexA,u0_c,WEIGHTBENDFIELD * val_u1_c * (-1.0) ));

                        // (u0b-u0c) * <v1c> ^2
                        tripletList->push_back(T(lastIndexA,t2,WEIGHTBENDFIELD * val_v1_c * t2_sign));
                        tripletList->push_back(T(lastIndexA,u0_c,WEIGHTBENDFIELD * val_v1_c * (-1.0)));
                        lastIndexA++;

                        // (v0r-v0c) * <u1c> +
                        tripletList->push_back(T(lastIndexA,t3,WEIGHTBENDFIELD * val_u1_c * t3_sign));
                        tripletList->push_back(T(lastIndexA,v0_c,WEIGHTBENDFIELD * val_u1_c * (-1.0)));

                        // (v0b-v0c) * <v1c> ^2
                        tripletList->push_back(T(lastIndexA,t4,WEIGHTBENDFIELD * val_v1_c * t4_sign));
                        tripletList->push_back(T(lastIndexA,v0_c,WEIGHTBENDFIELD * val_v1_c * (-1.0)));
                        lastIndexA++;

                        // Second term
                        // (u1r-u1c) * <u0c> +
                        tripletList->push_back(T(lastIndexA,t5,WEIGHTBENDFIELD * val_u0_c * t5_sign));
                        tripletList->push_back(T(lastIndexA,u1_c,WEIGHTBENDFIELD * val_u0_c * (-1.0)));

                        // (u1b-u1c) * <v0c> ^2
                        tripletList->push_back(T(lastIndexA,t6,WEIGHTBENDFIELD * val_v0_c * t6_sign));
                        tripletList->push_back(T(lastIndexA,u1_c,WEIGHTBENDFIELD * val_v0_c * (-1.0)));
                        lastIndexA++;

                        // (v1r-v1c) * <u0c> +
                        tripletList->push_back(T(lastIndexA,t7,WEIGHTBENDFIELD * val_u0_c * t7_sign));
                        tripletList->push_back(T(lastIndexA,v1_c,WEIGHTBENDFIELD * val_u0_c * (-1.0)));

                        // (v1b-v1c) * <v0c> ^2
                        tripletList->push_back(T(lastIndexA,t8,WEIGHTBENDFIELD * val_v0_c * t8_sign));
                        tripletList->push_back(T(lastIndexA,v1_c,WEIGHTBENDFIELD * val_v0_c * (-1.0)));
                        lastIndexA++;
                    }
                }

                // ADD STROKE CONSTRAINTS

                if(this->crossfield->isConstrainedLine(i,j))
                {

                    double m_weightConstraint = WEIGHTCONSTRAINT*this->crossfield->getStrokeWeight(i,j);

                    if(!this->crossfield->isConstraintFliped(i,j))
                    {
                        // u0c = u0
                        tripletList->push_back(T(lastIndexA,u0_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*cos(crossfield->getThita(i,j))));
                        lastIndexA++;

                        // v0c = v0
                        tripletList->push_back(T(lastIndexA,v0_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*sin(crossfield->getThita(i,j))));
                        lastIndexA++;
                    }
                    else
                    {
                        // u1c = u0
                        tripletList->push_back(T(lastIndexA,u1_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*cos(crossfield->getThita(i,j))* (-1.0)));
                        lastIndexA++;

                        // v1c = v0
                        tripletList->push_back(T(lastIndexA,v1_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*sin(crossfield->getThita(i,j))* (-1.0)));
                        lastIndexA++;
                    }

                }

                // ADD REGULARIZATION CONSTRAINTS

                if((!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                    &&!((crossfield->isBoundaryPixel(i,j+1))&&(!crossfield->isCurvatureLine(i,j+1))))
                {
                    // Handle bottom boundary
                    if(this->index->isUnknown(i,j+1)==true)
                    {
                        //(u0r-u0c)
                        tripletList->push_back(T(lastIndexA,t1,WEIGHTREGULARIZATION * t1_sign));
                        tripletList->push_back(T(lastIndexA,u0_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;

                        //+(v0r-v0c)
                        tripletList->push_back(T(lastIndexA,t3,WEIGHTREGULARIZATION * t3_sign));
                        tripletList->push_back(T(lastIndexA,v0_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;

                        //(u1r-u1c)
                        tripletList->push_back(T(lastIndexA,t5,WEIGHTREGULARIZATION * t5_sign));
                        tripletList->push_back(T(lastIndexA,u1_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;

                        //+(v1r-v1c)
                        tripletList->push_back(T(lastIndexA,t7,WEIGHTREGULARIZATION * t7_sign));
                        tripletList->push_back(T(lastIndexA,v1_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;
                    }
                }

                if((!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                     &&!((crossfield->isBoundaryPixel(i+1,j))&&(!crossfield->isCurvatureLine(i+1,j))))
                {
                    // Handle bottom boundary
                    if(this->index->isUnknown(i+1,j)==true)
                    {
                        //(u0b-u0c)
                        tripletList->push_back(T(lastIndexA,t2,WEIGHTREGULARIZATION * t2_sign));
                        tripletList->push_back(T(lastIndexA,u0_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;

                        //(v0b-v0c)
                        tripletList->push_back(T(lastIndexA,t4,WEIGHTREGULARIZATION * t4_sign));
                        tripletList->push_back(T(lastIndexA,v0_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;

                        //(u1b-u1c)
                        tripletList->push_back(T(lastIndexA,t6,WEIGHTREGULARIZATION * t6_sign));
                        tripletList->push_back(T(lastIndexA,u1_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;

                        //(v1b-v1c)
                        tripletList->push_back(T(lastIndexA,t8,WEIGHTREGULARIZATION * t8_sign));
                        tripletList->push_back(T(lastIndexA,v1_c,WEIGHTREGULARIZATION*(-1.0)));
                        lastIndexA++;
                    }
                }
            }
        }
    }

    // Set the system vectors, A and b
    A.resize(lastIndexA,sizeX);
    A.reserve(VectorXd::Constant(sizeX,6));
    A.setFromTriplets(tripletList->begin(),tripletList->end());
    A.makeCompressed();

    b.resize(lastIndexA);

    for(int i = 0; i < lastIndexA; i++)
    {
       b(i) = 0;
    }

    for(int i = 0; i < b_values.size();i++)
    {
       b(b_values.at(i).first) += b_values.at(i).second;
    }

    delete tripletList;
}

// Iterative method for solving the BendField
void BendField::smoothBendField(GLWidget * glwidget)
{
    // The matrix for solver
    SpMat A;
    SpMat A_transpose;
    SpMat A_final;

    VectorXd b;
    VectorXd b_final;

    VectorXd x;

    // Inits
    b.resize(sizeX);
    b_final.resize(sizeX);
    x.resize(sizeX);

    // Iterative algorithm to solve for non linear energy
    int iter = 0;

    // While not converge
    // Empirical test: 10 iteration is more than what's really needed
    while(iter < 10)
    {
        // Clean and init
        for(int i = 0; i < sizeX; i++)
        {       
            x(i) = 0;
            b(i) = 0;
            b_final(i) = 0;
        }

        // The first iteration is for building a harmonic solution
        if(iter == 0)
        {
            // Build the system (harmonic)
            qDebug() << "Building harmonics..";
            buildHarmonic(A,b,sizeX);
        }
        else
        {
            // Build the system (bendfield)
            qDebug() << "Building bendfield..";
            buildCovariant(A,b,sizeX);
        }

        qDebug() << "Built!";

        // Solve the system (normal equation)
        // Find Transpose
        A_transpose =  SpMat(A.transpose());

        // Multiplication
        A_final = A_transpose * A;
        A_final.makeCompressed();

        // The same for b
        b_final.resize(sizeX);
        b_final = A_transpose * b;

        // Solver
        // Test with 2 tolerances
        // in case it can't find a solution
        ConjugateGradient<SpMat> chol;
        chol.compute(A_final);

       if(chol.info()!=Success) {
           qDebug() << "Compute failed" ;
           return;
        }
        else
        {
           qDebug() << "Computed!" ;
        }

        chol.setTolerance(0.001);
        // Solve
        x = chol.solve(b_final);

        if(chol.info()!=Success) {
            qDebug() << "Solver failed" ;
            return;
         }
         else
         {
            qDebug() << "Solved!" ;
         }

        // Save in the corresponding vector field
        this->saveIntoField(x);

        // Show residual
        this->computeError(A,b,x,sizeX);

        iter += 1;

        qDebug() << "iter " << iter;


        // If last iteration
        // Orient consistently the crosses in the field
        if(iter==9)
        {
            this->crossfield->rotateCrosses(pjumpfield, this->mask, this->index);
        }

         glwidget->repaint();
    }

    // Normalize
    this->crossfield->normalize();

    qDebug() << "Exit.. ";
}

// Saves the harmonic solution of the system stored in X
// to the crossfield representation
void BendField::saveIntoField(VectorXd x)
{
    for(int k = 0; k < numerOfEntries; k++)
    {
       // The position in the field
       QPoint fieldPosition = this->getPixelInField(k);

       // u,v
       double u_0 = (double)x(k);
       double v_0 = (double)x((this->numerOfEntries) + k);
       double u_1 = (double)x((2*this->numerOfEntries) + k);
       double v_1 = (double)x((3*this->numerOfEntries) + k);

       // if not constrained
       this->crossfield->setV0V1(u_0,v_0,u_1,v_1,fieldPosition.x(),fieldPosition.y());

    }
}


// Builds a Harmonic energy for the new representation (u,v)
// instead of the used in the previous step (alpha,beta)
// It's only called once, to initialice the field
// It relies in the computed P
void BendField::buildHarmonic(SpMat &A, VectorXd& b, int sizeX)
{
    // Triplets init
    vector<T>  * tripletList = new vector<T>;
    tripletList->reserve(sizeX*6);
    tripletList->clear();

    Mat SmoothRight = Mat::zeros(h,w, CV_8UC3);
    Mat SmoothLeft = Mat::zeros(h,w, CV_8UC3);
    Mat Stroke = Mat::zeros(h,w, CV_8UC3);

    // for values in b
    QList<QPair<int,double> > b_values;

    int lastIndexA = 0;

    // Inside strokes, the regularization weight is smaller

    // Set up the constraints
    // For every pixel
    for(int i = 2; i < this->h-2; i++)
    {
        for(int j = 2; j < this->w-2; j++)
        {
            // White pixel in the mask
            if(this->index->isUnknown(i,j)==true)
            {
                // Get the coordinates for every variable involved
                // V0 Center
                int u0_c = getColInAfor_u0(i,j);
                int v0_c = getColInAfor_v0(i,j);

                // V0 Right Neighbourd
                int u0_r = getColInAfor_u0(i,j+1);
                int v0_r = getColInAfor_v0(i,j+1);

                // V0 Bottom Neighbourd
                int u0_b = getColInAfor_u0(i+1,j);
                int v0_b = getColInAfor_v0(i+1,j);

                // V1 Center
                int u1_c = getColInAfor_u1(i,j);
                int v1_c = getColInAfor_v1(i,j);

                // V1 Right Neighbourd
                int u1_r = getColInAfor_u1(i,j+1);
                int v1_r = getColInAfor_v1(i,j+1);

                // V1 Bottom Neighbourd
                int u1_b = getColInAfor_u1(i+1,j);
                int v1_b = getColInAfor_v1(i+1,j);

                // Get p's (topology)
                int rightP = (this->pjumpfield->getRightP(i,j))%4;
                int bottomP = (this->pjumpfield->getBottomP(i,j))%4;

                // To handle flips
                int t1 = 0;
                int t2 = 0;
                int t3 = 0;
                int t4 = 0;
                int t5 = 0;
                int t6 = 0;
                int t7 = 0;
                int t8 = 0;

                // sign to flip V0 or V1
                double t1_sign = 1.0;
                double t2_sign = 1.0;
                double t3_sign = 1.0;
                double t4_sign = 1.0;
                double t5_sign = 1.0;
                double t6_sign = 1.0;
                double t7_sign = 1.0;
                double t8_sign = 1.0;


                // RIGHT
                // V0 --> V0(u0,v0)
                // V1 --> V1(u1,v1)
                if(rightP==0)
                {
                   t1 = u0_r;
                   t3 = v0_r;
                   t5 = u1_r;
                   t7 = v1_r;
                }
                // V0 --> -V1(u1,v1)
                // V1 --> V0(u0,v0)
                else if((rightP==1)||(rightP==-3))
                {
                    t1 = u1_r;
                    t3 = v1_r;
                    t5 = u0_r;
                    t7 = v0_r;

                    t1_sign = -1.0;
                    t3_sign = -1.0;
                }
                // V0 --> V1(u1,v1)
                // V1 --> -V0(u0,v0)
                else if((rightP==-1)||(rightP==3))
                {
                    t1 = u1_r;
                    t3 = v1_r;
                    t5 = u0_r;
                    t7 = v0_r;

                    t5_sign = -1.0;
                    t7_sign = -1.0;
                }
                // V0 --> -V0(u0,v0)
                // V1 --> -V1(u1,v1)
                else if((rightP==2)||(rightP==-2))
                {
                    t1 = u0_r;
                    t3 = v0_r;
                    t5 = u1_r;
                    t7 = v1_r;

                    t1_sign = -1.0;
                    t3_sign = -1.0;
                    t5_sign = -1.0;
                    t7_sign = -1.0;
                }

                // BOTTOM
                // V0 --> V0(u0,v0)
                // V1 --> V1(u1,v1)
                if(bottomP==0)
                {
                    t2 = u0_b;
                    t4 = v0_b;
                    t6 = u1_b;
                    t8 = v1_b;
                }
                // V0 --> -V1(u1,v1)
                // V1 --> V0(u0,v0)
                else if((bottomP==1)||(bottomP==-3))
                {
                    t2 = u1_b;
                    t4 = v1_b;
                    t6 = u0_b;
                    t8 = v0_b;

                    t2_sign = -1.0;
                    t4_sign = -1.0;
                }
                // V0 --> V1(u1,v1)
                // V1 --> -V0(u0,v0)
                else if((bottomP==-1)||(bottomP==3))
                {
                    t2 = u1_b;
                    t4 = v1_b;
                    t6 = u0_b;
                    t8 = v0_b;

                    t6_sign = -1.0;
                    t8_sign = -1.0;
                }
                // V0 --> -V0(u0,v0)
                // V1 --> -V1(u1,v1)
                else if((bottomP==2)||(bottomP==-2))
                {
                    t2 = u0_b;
                    t4 = v0_b;
                    t6 = u1_b;
                    t8 = v1_b;

                    t2_sign = -1.0;
                    t4_sign = -1.0;
                    t6_sign = -1.0;
                    t8_sign = -1.0;
                }

                // ADD STROKE CONSTRAINTS
                if(this->crossfield->isConstrainedLine(i,j))
                {
                    Stroke.at<Vec3b>(i,j)[1] = 255;

                    double m_weightConstraint = this->crossfield->getStrokeWeight(i,j);

                    m_weightConstraint = 0.4 * m_weightConstraint;

                    if(!this->crossfield->isConstraintFliped(i,j))
                    {
                        // u0c = u0
                        tripletList->push_back(T(lastIndexA,u0_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*cos(crossfield->getThita(i,j))));
                        lastIndexA++;

                        // v0c = v0
                        tripletList->push_back(T(lastIndexA,v0_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*sin(crossfield->getThita(i,j))));
                        lastIndexA++;
                    }
                    else
                    {
                        // u1c = u0
                        tripletList->push_back(T(lastIndexA,u1_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*cos(crossfield->getThita(i,j))* (-1.0)));
                        lastIndexA++;

                        // v1c = v0
                        tripletList->push_back(T(lastIndexA,v1_c,m_weightConstraint));
                        b_values.push_back(QPair<int,double>(lastIndexA,m_weightConstraint*sin(crossfield->getThita(i,j))* (-1.0)));
                        lastIndexA++;
                    }
                }


                // ADD REGULARIZATION CONSTRAINTS
                if((!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                  &&!((crossfield->isBoundaryPixel(i,j+1))&&(!crossfield->isCurvatureLine(i,j+1))))
                {
                    // Handle bottom boundary
                    if(this->index->isUnknown(i,j+1)==true)
                    {
                        // Handle right border
                        SmoothRight.at<Vec3b>(i,j)[2] = 255;
                        SmoothRight.at<Vec3b>(i,j+1)[2] = 255;

                        //(u0r-u0c)
                        tripletList->push_back(T(lastIndexA,t1,WEIGHTHARMONICSMOOTH* t1_sign));
                        tripletList->push_back(T(lastIndexA,u0_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;

                        //+(v0r-v0c)
                        tripletList->push_back(T(lastIndexA,t3,WEIGHTHARMONICSMOOTH* t3_sign));
                        tripletList->push_back(T(lastIndexA,v0_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;

                        //(u1r-u1c)
                        tripletList->push_back(T(lastIndexA,t5,WEIGHTHARMONICSMOOTH* t5_sign));
                        tripletList->push_back(T(lastIndexA,u1_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;

                        //+(v1r-v1c)
                        tripletList->push_back(T(lastIndexA,t7,WEIGHTHARMONICSMOOTH* t7_sign));
                        tripletList->push_back(T(lastIndexA,v1_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;


                    }
                }

                if((!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                  &&!((crossfield->isBoundaryPixel(i+1,j))&&(!crossfield->isCurvatureLine(i+1,j))))
                {
                    // Handle bottom boundary
                    if(this->index->isUnknown(i+1,j)==true)
                    {
                        SmoothLeft.at<Vec3b>(i,j)[1] = 255;
                        SmoothLeft.at<Vec3b>(i+1,j)[1] = 255;

                        //(u0b-u0c)
                        tripletList->push_back(T(lastIndexA,t2,WEIGHTHARMONICSMOOTH* t2_sign));
                        tripletList->push_back(T(lastIndexA,u0_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;

                        //(v0b-v0c)
                        tripletList->push_back(T(lastIndexA,t4,WEIGHTHARMONICSMOOTH* t4_sign));
                        tripletList->push_back(T(lastIndexA,v0_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;

                        //(u1b-u1c)
                        tripletList->push_back(T(lastIndexA,t6,WEIGHTHARMONICSMOOTH* t6_sign));
                        tripletList->push_back(T(lastIndexA,u1_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;

                        //(v1b-v1c)
                        tripletList->push_back(T(lastIndexA,t8,WEIGHTHARMONICSMOOTH* t8_sign));
                        tripletList->push_back(T(lastIndexA,v1_c,WEIGHTHARMONICSMOOTH* (-1.0)));
                        lastIndexA++;
                    }
                }
            }
        }
    }

    // Set A and b
    A.resize(lastIndexA,sizeX);
    A.reserve(VectorXd::Constant(sizeX,6));
    A.setFromTriplets(tripletList->begin(),tripletList->end());
    A.makeCompressed();

    b.resize(lastIndexA);

    for(int i = 0; i < lastIndexA; i++)
    {
       b(i) = 0;
    }

    for(int i = 0; i < b_values.size();i++)
    {
       b(b_values.at(i).first) += b_values.at(i).second;
    }

    // ** Debugging **
    imwrite("output/stroke_cov.png",Stroke);
    imwrite("output/smooth_r_cov.png",SmoothRight);
    imwrite("output/smooth_l_cov.png",SmoothLeft);
    // ** Debugging **

    delete tripletList;
}
