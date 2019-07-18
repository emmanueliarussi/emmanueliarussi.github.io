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

#include "harmoniccrossfield.h"

// Weight of the B small constraint
// called  W_b in the paper (eq. 5)
#define BETA_W 0.000001

// Weight of the harmonic smoothing
// constraints. Not used -> 1
#define WEIGHTSMOOTH 1

// Weight of the constrained crosses
// Not used -> 1
#define WEIGHTSMOOTHSTROKE 1

// How many stitchings in one step?
#define MAXACCUMSTITCH 1

// Distinguishable value of a P-Jump
// when not initialized
#define INIT_P_VALUE -60

HarmonicCrossField::HarmonicCrossField() {}

HarmonicCrossField::HarmonicCrossField(CrossField * c,PeriodJumpField * p, UnknownsIndexer * idx, Mat m)
{
    // Point crossfield and pjumpfield
    this->crossfield = c;
    this->pjumpfield = p;
    this->mask = m;

    // Index
    this->index = idx;

    // Sizes
    this->h = this->crossfield->height();
    this->w = this->crossfield->width();

    qDebug() << "SIZE PIXELS" << index->getNumberOfPixels();
    qDebug() << "SIZE X" << index->getNumberOfUknowns();
    qDebug() << "SIZE A " <<index->getNumberOfUknowns() << "x" << index->getNumberOfUknowns();;

    // Each pixel has 4 equations (2 for each edge, right and bottom)
    this->A.resize(index->getNumberOfUknowns(),index->getNumberOfUknowns());
    this->b.resize(index->getNumberOfUknowns());
    this->x.resize(index->getNumberOfUknowns());

    this->inits();

    this->A.reserve(VectorXi::Constant(index->getNumberOfUknowns(),6));

    // Set up the system
    this->setUpSystem();

    this->A.makeCompressed();
}

// Returns period jump information for a single cross [i,j] (right and bottom P)
int * HarmonicCrossField::getP(PeriodJumpField * pjumpfield, int i, int j)
{
    int * p = new int[4];

    p[0] = pjumpfield->getRightP(i,j);
    p[1] = pjumpfield->getBottomP(i,j);

    return p;
}

// Return Alpha angle from the cross [i,j] neighbours (right and bottom)
double * HarmonicCrossField::getAlpha(CrossField * crossField, int i, int j)
{
    double * alpha  = new double[2];

    alpha[0] = crossField->getAlpha(i,j+1);
    alpha[1] = crossField->getAlpha(i+1,j);

    return alpha;
}

// Return Beta angle from the cross [i,j] neighbours (right and bottom)
double * HarmonicCrossField::getBeta(CrossField * crossField, int i, int j)
{
    double * beta  = new double[2];

    beta[0] = crossField->getBeta(i,j+1);
    beta[1] = crossField->getBeta(i+1,j);

    return beta;
}

// Saves the harmonic solution of the system stored in X
// to the crossfield representation
void HarmonicCrossField::saveIntoCrossfield()
{
    for(int i = 0; i < index->getNumberOfPixels(); i++)
    {
        // The position in the field (mapping the X 1D array to 2D crossfield)
        QPoint fieldPosition = this->getPixelInField(i);

        // Alpha and Beta
        double alpha = x(i);
        double beta =  x(index->getNumberOfPixels() + i);

        // Set it in the crossfield
        this->crossfield->setAlphaBetha(alpha,beta,fieldPosition.x(),fieldPosition.y());
    }
}

// Maps the 2D crossfield to X 1D array for Alpha variables
int HarmonicCrossField::getColInAforAlpha(int i,int j)
{
    return this->index->getFromIndex(i,j,0);
}

// Maps the 2D crossfield to X 1D array for Beta variables
int HarmonicCrossField::getColInAforBeta(int i,int j)
{
    return this->index->getFromIndex(i,j,1);
}

// Maps the X 1D array to 2D crossfield
QPoint HarmonicCrossField::getPixelInField(int row)
{    
    return index->getPixelIndexed(row);
}

// Initialization of the arrays for the Eigen solver (X,b)
void  HarmonicCrossField::inits()
{
    // Init x, B
    for(int i = 0; i < index->getNumberOfUknowns(); i++)
    {
        // = 0
        this->x(i) = 0;
        this->b(i) = 0;
    }
}

// Set up of the harmonic system  (eq. 5 in the paper)
void HarmonicCrossField::setUpSystem()
{
    // ** Debugging **
    // ** Map of the constrained pixels for debugging **
    Mat BetaSmall = Mat::zeros(h,w, CV_8UC3);
    Mat SmoothRight = Mat::zeros(h,w, CV_8UC3);
    Mat SmoothBottom = Mat::zeros(h,w, CV_8UC3);
    Mat Stroke = Mat::zeros(h,w, CV_8UC3);
    Mat edges_r = Mat::zeros(h,w, CV_8UC3);

    // Triplets for insertion in the Sparse Matrix A
    tripletList = new vector<T>;
    tripletList->reserve(index->getNumberOfUknowns()*6);

    // For every pixel (No border is considered)
    for(int i = 2; i < this->h-2; i++)
    {
        for(int j = 2; j < this->w-2; j++)
        {
            // If it's a pixel to solve for: The index considers it
            if(this->index->isUnknown(i,j)==true)
            {
                // Given the two neigbours ĵ P values
                int rightP = this->pjumpfield->getRightP(i,j);
                int bottomP = this->pjumpfield->getBottomP(i,j);

                // And the i, î (mapping from 2D crossfield to 1D in X vector)
                int col_i_alpha = getColInAforAlpha(i,j);
                int col_i_beta = getColInAforBeta(i,j);

                // Set up X with alpha an beta
                x(col_i_alpha) = this->crossfield->getAlpha(i,j);
                x(col_i_beta) = this->crossfield->getBeta(i,j);

                // The same for the other neighbour
                // The col corresponding to j
                int col_j_alpha = getColInAforAlpha(i,j+1);
                int col_j_beta = getColInAforBeta(i,j+1);               

                // Set up constraints :
                // Using Eigen triplets mechanism : http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html

                // CONSTRAINT BETA_i SMALL

                // A(î,î) += BETA_W * 2
                tripletList->push_back(T(col_i_beta,col_i_beta,BETA_W*2));

                // b(î) += 0
                b(col_i_beta) += 0;

                // ** Debugging **
                BetaSmall.at<Vec3b>(i,j)[0] = 255;
                // ** Debugging **

                // CONSTRAINTS SMOOTHNESS

                // Add values to A and b for every defined (p=0) edge              
                if(((rightP==0)&&!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                    &&!((crossfield->isBoundaryPixel(i,j+1))&&(!crossfield->isCurvatureLine(i,j+1))))
                {
                    // Handle right boundary
                    if(this->index->isUnknown(i,j+1)==true)
                    {

                        // ** Debugging **
                        SmoothRight.at<Vec3b>(i,j)[1] = 255;
                        SmoothRight.at<Vec3b>(i,j+1)[1] = 255;
                        // ** Debugging **

                        // New values in A for alpha
                        // A(i,i) += 4
                        tripletList->push_back(T(col_i_alpha,col_i_alpha,WEIGHTSMOOTH* 4));
                        // A(j,j) += 4
                        tripletList->push_back(T(col_j_alpha,col_j_alpha,WEIGHTSMOOTH* 4));
                        // A(i,j) += -4
                        tripletList->push_back(T(col_i_alpha,col_j_alpha,WEIGHTSMOOTH* -4));
                        // A(j,i) += -4
                        tripletList->push_back(T(col_j_alpha,col_i_alpha,WEIGHTSMOOTH* -4));
                        // Crossed between î and j, ĵ and i are = 0
                        // A(î,î) += 4*(1-2*b)^2 = 4*(-1)^2 = 4
                        tripletList->push_back(T(col_i_beta,col_i_beta, WEIGHTSMOOTH* 4));
                        // A(ĵ,ĵ) += 4
                        tripletList->push_back(T(col_j_beta,col_j_beta, WEIGHTSMOOTH* 4));
                        // A(î,ĵ) += 4*(2*b-1) | b=0 => 4*(2*b-1) = 4*(-1) = -4
                        tripletList->push_back(T(col_i_beta,col_j_beta, WEIGHTSMOOTH* -4));
                        // A(ĵ,î) += 4*(2*b-1) | b=0 => 4*(2*b-1) = 4*(-1) = -4
                        tripletList->push_back(T(col_j_beta,col_i_beta, WEIGHTSMOOTH* -4));


                        // New row in B for alpha
                        // b(i) += -2*PI*p = 0
                        b(col_i_alpha) += 0;
                        // b(j) += -2*PI*p = 0
                        b(col_j_alpha) += 0;
                        // b(î) += 0
                        b(col_i_beta) += 0;
                        // b(ĵ) += 0
                        b(col_j_beta) += 0;
                    }
                }
                // The col corresponding to j
                col_j_alpha = getColInAforAlpha(i+1,j);
                col_j_beta = getColInAforBeta(i+1,j);              

                if(((bottomP==0)&&!((crossfield->isBoundaryPixel(i,j))&&(!crossfield->isCurvatureLine(i,j))))
                     &&!((crossfield->isBoundaryPixel(i+1,j))&&(!crossfield->isCurvatureLine(i+1,j))))
                {
                    // Handle bottom boundary
                    if(this->index->isUnknown(i+1,j)==true)
                    {
                        // ** Debugging **
                        SmoothBottom.at<Vec3b>(i,j)[1] = 255;
                        SmoothBottom.at<Vec3b>(i+1,j)[1] = 255;
                        // ** Debugging **

                        // New values in A for alpha
                        // A(i,i) += 4
                        tripletList->push_back(T(col_i_alpha,col_i_alpha,WEIGHTSMOOTH* 4));
                        // A(j,j) += 4
                        tripletList->push_back(T(col_j_alpha,col_j_alpha,WEIGHTSMOOTH* 4));
                        // A(i,j) += -4
                        tripletList->push_back(T(col_i_alpha,col_j_alpha,WEIGHTSMOOTH* -4));
                        // A(j,i) += -4
                        tripletList->push_back(T(col_j_alpha,col_i_alpha,WEIGHTSMOOTH* -4));
                        // Crossed between î and j, ĵ and i are = 0
                        // A(î,î) += 4*(1-2*b)^2 = 4*(-1)^2 = 4
                        tripletList->push_back(T(col_i_beta,col_i_beta,WEIGHTSMOOTH* 4));
                        // A(ĵ,ĵ) += 4
                        tripletList->push_back(T(col_j_beta,col_j_beta,WEIGHTSMOOTH* 4));
                        // A(î,ĵ) += 4*(2*b-1) | b=0 => 4*(2*b-1) = 4*(-1) = -4
                        tripletList->push_back(T(col_i_beta,col_j_beta,WEIGHTSMOOTH* -4));
                        // A(ĵ,î) += 4*(2*b-1) | b=0 => 4*(2*b-1) = 4*(-1) = -4
                        tripletList->push_back(T(col_j_beta,col_i_beta,WEIGHTSMOOTH* -4));

                        // New row in B for alpha
                        // b(i) += -2*PI*p = 0
                        b(col_i_alpha) += 0;
                        // b(j) += -2*PI*p = 0
                        b(col_j_alpha) += 0;
                        // b(î) += 0
                        b(col_i_beta) += 0;
                        // b(ĵ) += 0
                        b(col_j_beta) += 0;
                    }
                }       

                // CONSTRAINT STROKE
                // If it is a constrained pixel  W* (a+b)=  W* Titha tripletList.push_back(T(,,));
                if(this->crossfield->isConstrainedLine(i,j))
                {
                    // ** Debugging **
                    Stroke.at<Vec3b>(i,j)[2] = 255;
                    // ** Debugging **

                    // keep the weighting from the stroke
                    double weightConstraint = this->crossfield->getStrokeWeight(i,j); //*0.4

                    // A(i,i) += 2
                    tripletList->push_back(T(col_i_alpha,col_i_alpha, weightConstraint* 2));
                    // A(î,î) += 2
                    tripletList->push_back(T(col_i_beta,col_i_beta, weightConstraint* 2));
                    // A(i,î) += 2
                    tripletList->push_back(T(col_i_alpha,col_i_beta, weightConstraint* 2));
                    // A(î,i) += 2
                    tripletList->push_back(T(col_i_beta,col_i_alpha,weightConstraint*2));

                    // b(i) += 2*thita
                    b(col_i_alpha) +=  weightConstraint * 2 * this->crossfield->getThita(i,j);
                    // b(î) += 2*thita
                    b(col_i_beta) += weightConstraint * 2 * this->crossfield->getThita(i,j);
                }
            }
       }
    }


    // ** Debugging **
    // ** Map of the constrained pixels for debugging **
    imwrite("output/strokeConstrained.png",Stroke);
    imwrite("output/smoothRight.png",SmoothRight);
    imwrite("output/smoothBottom.png",SmoothBottom);
    imwrite("output/betaSmall.png",BetaSmall);
    imwrite("output/edges_r.png",edges_r);
    // ** Debugging **

    // Set matrix from list
    this->A.setFromTriplets(tripletList->begin(), tripletList->end());

}

// Best P between 2 neighbour crosses computation (check out smoothness measure in the paper for details)
int HarmonicCrossField::getBestP(double alpha_i,double beta_i,double alpha_j,double beta_j, double & bestCost)
{
    double a  = (double)((alpha_j-alpha_i)/(0.5*M_PI));

    int i0 = (int)std::ceil (a);
    int i1 = (int)std::floor(a);

    double intpart;
    double fractpart;

    fractpart = modf(a, &intpart);

    if(fractpart != 0)
    {
        double cost_i0 = estimateCost(alpha_i,beta_i,alpha_j,beta_j,i0);
        double cost_i1 = estimateCost(alpha_i,beta_i,alpha_j,beta_j,i1);

        if(cost_i0 < cost_i1)
        {
            bestCost = cost_i0;
            return i0;
        }
        else
        {
            bestCost = cost_i1;
            return i1;
        }
    }
    else
    {

        double cost_i0 = estimateCost(alpha_i,beta_i,alpha_j,beta_j,(a-1));
        double cost_i1 = estimateCost(alpha_i,beta_i,alpha_j,beta_j,a);
        double cost_i2 = estimateCost(alpha_i,beta_i,alpha_j,beta_j,(a+1));

        bestCost = std::min(std::min(cost_i0,cost_i1),cost_i2);

        if(cost_i0 == bestCost)
        {
            return a-1;
        }
        else if(cost_i1 == bestCost)
        {
            return a;
        }
        else if (cost_i2 == bestCost)
        {
            return a+1;
        }
    }
}

double HarmonicCrossField::estimateCost(double alpha_i,double beta_i,double alpha_j,double beta_j,int p_ij)
{
    // P mod 2, avoiding negative values of p
    int p_mode_2 = ((p_ij%2)+2)%2;

    double weight_a = alpha_i + (p_ij * (double)(M_PI/2)) + (1-2*(p_mode_2))*beta_i - (alpha_j + beta_j);
    double weight_b = alpha_i + (p_ij * (double)(M_PI/2)) - (1-2*(p_mode_2))*beta_i - (alpha_j - beta_j);

    return (double) ((weight_a*weight_a) + (weight_b*weight_b));
}

// Implementation of Gauss Seidel method for solving
bool HarmonicCrossField::localGaussSeidel(QVector<int> residualsQueue, QSet<int> elementsInQueue)
{
    // To update x.
    // In case of convergence, x = x_aux and return true
    VectorXd x_aux(this->x);

    // Gauss-Seidel
    int iter = 0;

    while((!residualsQueue.empty())&&(iter < 25))
    {
        iter = iter + 1;

        // take the current non-zero residual item
        int current_nonZero_index = residualsQueue.first();
        residualsQueue.pop_front();

        // new residual r_i = b_i - SUM A_ik x_i
        double r_iter = this->b(current_nonZero_index) - (this->A.col(current_nonZero_index).dot(x_aux)) ;

        if(abs(r_iter) > 0.1)
        {
            // Update x
            x_aux(current_nonZero_index) = x_aux(current_nonZero_index) + (r_iter / this->A.coeff(current_nonZero_index,current_nonZero_index));

            for (SpMat::InnerIterator it(this->A,current_nonZero_index); it; ++it)
            {
                if(!elementsInQueue.contains(it.row()))
                {
                    residualsQueue.push_front(it.row());
                    elementsInQueue.insert(it.row());
                }
            }
        }
    }

    // If it converged
    if(residualsQueue.empty())
    {
        this->x = x_aux;
        return true;
    }
    else
    {
        return false;
    }
}

// In Order insertion implementation (for Gauss Seidel method)
void HarmonicCrossField::insertInOrder(QList<StitchData> & list, StitchData data)
{
    if(list.empty())
    {
        list.push_front(data);
    }
    else
    {
        int imax = list.size()-1;
        int imin = 0;


        // continue searching while [imin,imax] is not empty
        while (imax >= imin)
        {
            /* calculate the midpoint for roughly equal partition */
            int imid = (imin + imax)/2;

            // determine which subarray to search
            if (list.at(imid).cost <= data.cost)
            {
              // change min index to search upper subarray
              imin = imid + 1;
            }
            else if (list.at(imid).cost > data.cost)
            {
              // change max index to search lower subarray
              imax = imid - 1;
            }
        }

        list.insert(imin,data);
    }

}

// In Order insertion implementation (for Gauss Seidel method)
void HarmonicCrossField::insertInOrder_integer(QList<int> & list, int i)
{
    if(list.empty())
    {
        list.push_front(i);
    }
    else
    {
        int imax = list.size()-1;
        int imin = 0;

        // continue searching while [imin,imax] is not empty
        while (imax >= imin)
        {
            /* calculate the midpoint for roughly equal partition */
            int imid = (imin + imax)/2;

            // determine which subarray to search
            if (list.at(imid) > i)
            {
              // change min index to search upper subarray
              imin = imid + 1;
            }
            else if (list.at(imid) <= i)
            {
              // change max index to search lower subarray
              imax = imid - 1;
            }
        }

        list.insert(imin,i);
    }

}

// Iterative stitching method for solving the field
void HarmonicCrossField::smoothWithIterativeGreedy(GLWidget * glwidget)
{
    QTime totalTime;
    totalTime.start();

    // List of all non defined p (to be "stitched" = find a propper value for p)
    QVector<QPair<QPoint,QPoint> > * nonDefinedP = pjumpfield->getNonDefinedPList();

    // Print initial energy of the system
    this->printEnergy(crossfield,pjumpfield,0);

    // Stats for debugging
    int avgGauss = 0;
    int avgGradient = 0;
    int gaussCounter = 0;
    int gradientCounter = 0;

    // Fill X by solving once exactly:
    // performs a Cholesky factorization of A
    SimplicialCholesky<SpMat> chol(this->A);
    this->x = chol.solve(this->b);

    qDebug()  << "TO STITCH " << nonDefinedP->size();

    // While there is still some edges to stich
    while(!nonDefinedP->empty())
    {
        // List of all the stitches to do in one iteration
        QList<StitchData> toStitch;

        // For every candidate edge, find the best one to stich (less cost)
        // and add to the full system
        for(int i = 0; i < nonDefinedP->size(); i++)
        {
            // Get alphas and betas
            QPoint pixel_i(nonDefinedP->at(i).first.x(),nonDefinedP->at(i).first.y());
            QPoint pixel_j(nonDefinedP->at(i).second.x(),nonDefinedP->at(i).second.y());

            // Obtain alphas and betas from the solution vector
            double alpha_i = x(getColInAforAlpha(pixel_i.x(),pixel_i.y()));
            double beta_i = x(getColInAforBeta(pixel_i.x(),pixel_i.y()));

            double alpha_j = x(getColInAforAlpha(pixel_j.x(),pixel_j.y()));
            double beta_j = x(getColInAforBeta(pixel_j.x(),pixel_j.y()));

            double cost = INFINITY;

            // Find the best p and stimated quality comes in "cost"
            int approx_p = this->getBestP(alpha_i,beta_i,alpha_j,beta_j,cost);

            StitchData data;
            data.i = i;
            data.cost = cost;
            data.approx_p = approx_p;

            // Insert in order
            this->insertInOrder(toStitch,data);
        }

        // Residuals x
        QVector<int> residualsQueue;
        QSet<int> elementsInQueue;
        QVector<int> modifiedXIndex;

        VectorXd residuals;
        residuals = ((A * x) - b).cwiseAbs();

        double accum_cost = 0;

        QList<int> toRemove;

        // Stitch only first one
        // if(!toStitch.empty())
        // ** Add constraints to A
        while((accum_cost < MAXACCUMSTITCH)&&(!toStitch.empty()))
        {

            // Get p
            int best_p = toStitch.first().approx_p;

            // Accum cost
            accum_cost += toStitch.first().cost;

            // Stich best candidate
            // Remove it from the list of non defined
            QPair<QPoint,QPoint> bestCandidate = nonDefinedP->at(toStitch.first().i);

            // To remove it after from the
            insertInOrder_integer(toRemove,toStitch.first().i);

            // Remove edge from toStitch list
            toStitch.removeFirst();

            // Save p
            this->pjumpfield->setP(bestCandidate.first,bestCandidate.second,best_p);

            // P mod 2, avoiding negative values of p
            int best_p_mode_2 = ((best_p%2)+2)%2;

            // Add equations to the system (stitching)

            int first_i = bestCandidate.first.x();
            int first_j = bestCandidate.first.y();
            int second_i = bestCandidate.second.x();
            int second_j = bestCandidate.second.y();

            if((!((crossfield->isBoundaryPixel(first_i,first_j))&&(!crossfield->isCurvatureLine(first_i,first_j))))
              &&!((crossfield->isBoundaryPixel(second_i,second_j))&&(!crossfield->isCurvatureLine(second_i,second_j))))
            {                
                // The col corresponding to i
                int col_i_alpha = getColInAforAlpha(first_i,first_j);
                int col_i_beta = getColInAforBeta(first_i,first_j);

                // The col corresponding to j
                int col_j_alpha = getColInAforAlpha(second_i,second_j);
                int col_j_beta = getColInAforBeta(second_i,second_j);

                // To update residuals
                if(!elementsInQueue.contains(col_i_alpha))
                {
                    modifiedXIndex.push_back(col_i_alpha);
                    elementsInQueue.insert(col_i_alpha);
                }
                if(!elementsInQueue.contains(col_i_beta))
                {
                     modifiedXIndex.push_back(col_i_beta);
                     elementsInQueue.insert(col_i_beta);
                }
                if(!elementsInQueue.contains(col_j_alpha))
                {
                    modifiedXIndex.push_back(col_j_alpha);
                    elementsInQueue.insert(col_j_alpha);
                }
                if(!elementsInQueue.contains(col_j_beta))
                {
                    modifiedXIndex.push_back(col_j_beta);
                    elementsInQueue.insert(col_j_beta);
                }


                // CONSTRAINTS SMOOTHNESS

                // New values in A for alpha
                // A(i,i) += 4
                tripletList->push_back(T(col_i_alpha,col_i_alpha,WEIGHTSMOOTH* 4));
                // A(j,j) += 4
                tripletList->push_back(T(col_j_alpha,col_j_alpha,WEIGHTSMOOTH* 4));
                // A(i,j) += -4
                tripletList->push_back(T(col_i_alpha,col_j_alpha,WEIGHTSMOOTH* -4));
                // A(j,i) += -4
                tripletList->push_back(T(col_j_alpha,col_i_alpha,WEIGHTSMOOTH* -4));
                // Crossed between î and j, ĵ and i are = 0
                // A(î,î) += 4*(1-2*b)^2 = 4*(-1)^2 = 4
                tripletList->push_back(T(col_i_beta,col_i_beta, WEIGHTSMOOTH* 4 * (1-2*(best_p_mode_2)) * (1-2*(best_p_mode_2))));
                // A(ĵ,ĵ) += 4
                tripletList->push_back(T(col_j_beta,col_j_beta, WEIGHTSMOOTH* 4));
                // A(î,ĵ) += 4*(2*b-1) | b=0 => 4*(2*b-1) = 4*(-1) = -4
                tripletList->push_back(T(col_i_beta,col_j_beta, WEIGHTSMOOTH*  4*((2*(best_p_mode_2))-1)));
                // A(ĵ,î) += 4*(2*b-1) | b=0 => 4*(2*b-1) = 4*(-1) = -4
                tripletList->push_back(T(col_j_beta,col_i_beta, WEIGHTSMOOTH* 4*((2*(best_p_mode_2))-1)));

                // New row in B for alpha
                // b(i) += -2*PI*p = 0
                this->b(col_i_alpha) += WEIGHTSMOOTH* -2*M_PI*best_p;
                // b(j) += -2*PI*p = 0
                this->b(col_j_alpha) += WEIGHTSMOOTH* 2*M_PI*best_p;
                // b(î) += 0
                this->b(col_i_beta) += 0;
                // b(ĵ) += 0
                this->b(col_j_beta) += 0;
            }

        }

        // Remove from nonDefinedP (from the back to the front, to avoid altering the ordering)
        for(int i = 0; i < toRemove.size(); i++)
        {
            nonDefinedP->remove(toRemove.at(i));
        }

        // Set matrix from list
        this->A.setFromTriplets(tripletList->begin(), tripletList->end());
        this->A.makeCompressed();

        // Stack of options to solve the system
        // First, local Gauss Seidel by stimating residuals (just update x locally)

        // Residual = Old x - updated A times x - b
        residuals = residuals - ((A * x) - b).cwiseAbs();

        elementsInQueue.clear();

        // Push Residuals (Now filter by error)
        for(int k = 0; k < modifiedXIndex.size(); k++)
        {
            if(abs(residuals(modifiedXIndex.at(k))) > 0.1)
            {
                residualsQueue.push_front(modifiedXIndex.at(k));
                elementsInQueue.insert(modifiedXIndex.at(k));
            }
        }

        // Solve with the stack of solvers

        // First try local update (Fast)
        qDebug() << "TRY WITH GAUSS " << nonDefinedP->size();

        QTime t;
        t.start();
        bool gauss = localGaussSeidel(residualsQueue,elementsInQueue);
        avgGauss += t.elapsed();
        gaussCounter++;

        // Then, try Conjugate Gradient (Fast--)
        if(!gauss)
        {
            qDebug() << "TRY WITH GRADIENT " << nonDefinedP->size();

            VectorXd auxX = x;
            t.restart();
            // If it didn't converged
            // then, try global iterative
            ConjugateGradient<SparseMatrix<double> > cg;
            cg.compute(this->A);
            cg.setMaxIterations(1000);
            cg.setTolerance(0.001);
            auxX = cg.solveWithGuess(this->b,this->x);

            avgGradient += t.elapsed();
            gradientCounter++;

            // Then, try Cholesky (Slow)
            if(cg.info()==Success)
            {                
                this->x = auxX;
            }
            else
            {
                qDebug() << "TRY WITH CHOLESKY " << nonDefinedP->size();
                // If it didn't converged
                // then, try global cholesky
                SimplicialCholesky<SpMat> chol2(this->A);
                this->x = chol2.solve(this->b);
            }
        }

//        // FOR PLOTTING
//        if(toStitch.size()<400)
//        {
//            // Save into crossfield
//            saveIntoCrossfield();

//            // Generate V0-V1 field
//            this->crossfield->generate_V0V1_field();

//            // Print information
//            this->printEnergy(crossfield,pjumpfield,0);

//            // plot (this is here, and not in the MainWindow class, because sometimes I want to paint at each iteration of the algorithm, for debugging reasons)
//            glwidget->repaint();
//        }

    }

    // Run one last time exactly (Cholesky)
    // performs a Cholesky factorization of A
    SimplicialCholesky<SpMat> chol3(this->A);
    this->x = chol3.solve(this->b);
    if(chol3.info()==Success)
    {
        qDebug() << "Last exact solution: SUCCESS";
    }
    else
    {
        qDebug() << "Solve CG with a lot more presition";

        VectorXd auxX = x;
        ConjugateGradient<SparseMatrix<double> > cg;
        cg.compute(this->A);
        cg.setMaxIterations(1000000);
        cg.setTolerance(0.00001);
        auxX = cg.solveWithGuess(this->b,this->x);

        if(cg.info()==Success)
        {
            qDebug() << "SUCCESS!";
            this->x = auxX;
        }
        else
        {
            qDebug() << "WARNING: FAIL SOLVING";
        }
    }

    // Save into crossfield
    saveIntoCrossfield();

    // Generate V0-V1 field (change representation to be able to run BendField algorithm
    this->crossfield->generate_V0V1_field();

    // Print information
    this->printEnergy(crossfield,pjumpfield,0);

    // Plot on the interface
    // DEBUGGING
    glwidget->repaint();

    // Print stats
    if(gaussCounter != 0)
    {
        qDebug() << "GAUSS: " << (double) (avgGauss / gaussCounter);
    }
    if(gradientCounter!=0)
    {
        qDebug() << "GRADIENT: " << (double) (avgGradient / gradientCounter);
    }
    qDebug() << "UNKNOWNS: " << index->getNumberOfUknowns();
    qDebug() << "TOTAL TIME " << totalTime.elapsed() / 1000.;

}




// Computes and prints the total energy of the current optimization process
// Saves an error map (color coded) showing the less accurate regions
double HarmonicCrossField::printEnergy(CrossField * crossField,PeriodJumpField * pjumpfield,int iterations)
{
    // Output error image
    Mat output = Mat::zeros(crossField->height(),crossField->width(), CV_64FC3);

    // Energy accum
    double finalEnergy = 0;

    // Current status of the stitching process
    int stitchedCount = 0;

    // To keep track if the max error
    double maxValue = -INFINITY;

    // For each pair of crosses in the cross field (No border is considered)
    for(int i = 2; i < crossField->height()-2; i++)
    {
        for(int j = 2; j < crossField->width()-2; j++)
        {
            if((j!=(this->w-2)-1)&&(i!=(this->h-2)-1))
            {
                // Compute energy
                // See smoothness measure in the paper
                double energy = 0;
                int * p = this->getP(pjumpfield, i, j);
                double * alpha  = this->getAlpha(crossField,i,j);
                double * beta  = this->getBeta(crossField,i,j);

                // Only right and bottom are considered
                for(int k = 0; k < 2; k++)
                {
                    // If it's a valid edge (p has been initialized)
                    if(p[k]!=INIT_P_VALUE)
                    {
                        double energy = estimateCost(crossField->getAlpha(i,j),crossField->getBeta(i,j),alpha[k],beta[k],p[k]);
                        energy +=  energy;
                        stitchedCount++;
                    }
                }

                // Keep track of max error for normalization
                if(energy > maxValue)
                {
                    maxValue = energy;
                }

                // Save error on the red channel
                output.at<cv::Vec3d>(i,j)[2] = energy;

                // Contribute to the final energy
                finalEnergy += energy;
            }
        }
    }

    // Print the data
    qDebug() << " TOTAL ENERGY: " << finalEnergy << " - STITCHED: " << stitchedCount << " - RELATIVE " << (double) (finalEnergy / stitchedCount);

    // Save the image
    Mat realOutput = Mat::zeros(crossField->height(),crossField->width(), CV_8UC3);
    output = output / maxValue;
    realOutput = output * 255;
    imwrite("output/energy.png",realOutput);

    return finalEnergy;
}
