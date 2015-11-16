#ifndef RANDOMNORMALMATCHING_H_
#define RANDOMNORMALMATCHING_H_

#include <flann/flann.hpp>
#include "obcore/math/linalg/linalg.h"
#include "obvision/registration/Trace.h"
#include "obvision/registration/icp/PointToLineEstimator2D.h"
#include <cmath>
#include "omp.h"

namespace obvious
{

/**
 * @class RandomNormalMatching
 * @brief Matching algorithm with PCA alignment and RANSAC scheme
 * @author Stefan May
 **/
class RandomNormalMatching
{
public:
  /**
   * Constructor
   * @param trials number of trials / matching guesses
   * @param epsThresh threshold for rating good matches
   * @param phiMax maximum rotation
   * @param sizeControlSet approximate set of control set
   */

	RandomNormalMatching(
			unsigned int trials = 100,
			double epsThresh = 0.15,
			unsigned int sizeControlSet = 140,
			double zhit = 0.45,
			double zphi = 0.0,
			double zshort = 0.25,
			double zmax = 0.05,
			double zrand = 0.25,
			double percentagePointsInC = 0.9,
			double rangemax = 20,
			double sigphi = M_PI / 180.0 * 3,
			double sighit = 0.2,
			double lamshort = 0.08,
			double maxAngleDiff = 3.0,
			double maxAnglePenalty = 0.5);

  /**
   * Destructor
   */
  virtual ~RandomNormalMatching();

  /**
   * Activate internal trace writing. While the method iterate is executed, all states are traced.
   */
  void activateTrace();

  /**
   * Deactivate internal trace writing.
   */
  void deactivateTrace();

  /**
   * Matching method
   * @param M Matrix for model points. Points are accessed by rows. e.g. x = M(p, 0) y= M(p,1)
   * @param maskM Mask for matrix M. Valid points are identified with true-value. It has the same size as M.getCols()
   * @param NM Matrix with normals for model data set
   * @param S Matrix for scene points
   * @param maskS Mask for matrix S
   * @param phiMax Maximum allowed rotation as output
   * @param transMax Maximum allowed translation
   * @param resolution Angular resolution of the laser scan
   * @return 3x3 registration matrix
   */
  obvious::Matrix match(const obvious::Matrix* M,
                        const bool* maskM,
                        const obvious::Matrix* NM,
                        const obvious::Matrix* S,
                        const bool* maskS,
                        double phiMax = M_PI / 4.0,
                        const double transMax = 1.5,
                        const double resolution = 0.0);

	/**
	 * Matching method
	 * @param M Matrix for model points. Points are accessed by rows. e.g. x = M(p, 0) y= M(p,1)
	 * @param maskM Mask for matrix M. Valid points are identified with true-value. It has the same size as M.getCols()
	 * @param NM Matrix with normals for model data set
	 * @param S Matrix for scene points
	 * @param maskS Mask for matrix S
	 * @param phiMax Maximum allowed rotation as output
	 * @param transMax Maximum allowed translation
	 * @param resolution Angular resolution of the laser scan
	 * @return 3x3 registration matrix
	 */
	obvious::Matrix match2(const obvious::Matrix* M,
						  const bool* maskM,
						  const obvious::Matrix* NM,
						  const obvious::Matrix* S,
						  const bool* maskS,
						  double phiMax = M_PI / 4.0,
						  const double transMax = 1.5,
						  const double resolution = 0.0);

	/**
	* get probability of two single rangefinder scans (of two rays)
	* @param m distance of model point
	* @param s distance of scene/control point
	* @param phiDiff difference of the two angles
	*/
	double probabilityOfTwoSingleScans(double m, double s, double phiDiff);

  /**
   * Serialize assignment to trace folder
   * @param folder trace folder (must not be existent)
   */
  void serializeTrace(const char* folder);

private:

  // Calculate normals of point set
  void calcNormals(const Matrix* M, Matrix* N, const bool* maskIn, bool* maskOut);

  // Calculate angle of normals
  void calcPhi(const Matrix* N, const bool* mask, double* phi);

  // Subsample mask for better performance
  void subsampleMask(bool* mask, unsigned int size, double probability);

  // extract valid sample indices from matrix giving a validity mask
  vector<unsigned int> extractSamples(const obvious::Matrix* M, const bool* mask);

  // init kd-tree for fast NN search in model
  void initKDTree(const obvious::Matrix* M, vector<unsigned int> valid);

  // pick control set for RANSAC in-/outlier detection
  obvious::Matrix* pickControlSet(const obvious::Matrix* M, vector<unsigned int> idxValid, vector<unsigned int> &idxControl);

  // probability model variables
  double _zhit;
  double _zphi;
  double _zshort;
  double _zmax;
  double _zrand;

  double _percentagePointsInC;

  double _rangemax;
  double _sigphi;
  double _sighit;
  double _lamshort;

  double _maxAngleDiff;
  double _maxAnglePenalty;

  // squared distance threshold
  double _scaleDistance;

  // normalization weight for orientation rating
  double _scaleOrientation;

  // number of trials
  unsigned int _trials;

  // approximate control set
  unsigned int _sizeControlSet;

  // tree for accelerating NN search
  flann::Index<flann::L2<double> >* _index;
  flann::Matrix<double>* _model;

  // Trace module
  Trace* _trace;

  // Number of samples investigated for PCA in local neighborhood
  int _pcaSearchRange;

  // Min number of valid samples PCA
  unsigned int _pcaMinSamples;
};

}

#endif /* RANDOMNORMALMATCHING_H_ */
