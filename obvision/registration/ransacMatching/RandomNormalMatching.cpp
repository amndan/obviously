#include "RandomNormalMatching.h"

#include <math.h>
#include "obcore/base/System.h"
#include "obcore/base/Logger.h"
#include "obcore/math/mathbase.h"
#include <limits>

#include "obcore/base/Timer.h"

using namespace std;

namespace obvious {

#define NORMALCONSENSUS 1
#define USEKNN 1

RandomNormalMatching::RandomNormalMatching(unsigned int trials,
		double epsThresh, unsigned int sizeControlSet) {
	_scaleDistance = 1.0 / (epsThresh * epsThresh);
	_scaleOrientation = 0.33;
	_trials = trials;
	_sizeControlSet = sizeControlSet;
	_model = NULL;
	_index = NULL;
	_trace = NULL;
	_pcaSearchRange = 10;
	_pcaMinSamples = 3;

	_zhit = 0.45;
	_zshort = 0.25;
	_zmax = 0.05;
	_zrand = 0.25;

	_phit = 0;
	_pshort = 0;
	_pmax = 0;
	_prand = 0;

	_rangemax = 20;
	_sighit = 0.4;
	_lamshort = 0.08;
}

RandomNormalMatching::~RandomNormalMatching() {
	if (_model) {
		delete _model;
		_model = NULL;
		delete _index;
		_index = NULL;
	}
	if (_trace)
		delete _trace;
	_trace = NULL;
}

void RandomNormalMatching::activateTrace() {
	if (!_trace)
		_trace = new Trace(2);
}

void RandomNormalMatching::deactivateTrace() {
	if (_trace)
		delete _trace;
	_trace = NULL;
}

vector<unsigned int> RandomNormalMatching::extractSamples(
		const obvious::Matrix* M, const bool* mask) {
	vector<unsigned int> validIndices;
	for (unsigned int i = _pcaSearchRange / 2;
			i < M->getRows() - _pcaSearchRange / 2; i++) {
		if (mask[i])
			validIndices.push_back(i);
	}
	return validIndices;
}

void RandomNormalMatching::initKDTree(const obvious::Matrix* M,
		vector<unsigned int> idxValid) {
	// Build FLANN tree for fast access to nearest neighbors
	unsigned int cols = M->getCols();
	unsigned int rows = idxValid.size();
	double** mData;
	obvious::System<double>::allocate(rows, cols, mData);
	for (unsigned int r = 0; r < rows; r++) {
		mData[r][0] = (*M)(idxValid[r], 0);
		mData[r][1] = (*M)(idxValid[r], 1);
	}
	if (_model) {
		delete _model;
		_model = NULL;
		delete _index;
		_index = NULL;
	}
	_model = new flann::Matrix<double>(&mData[0][0], rows, 2);
	flann::KDTreeSingleIndexParams p;
	_index = new flann::Index<flann::L2<double> >(*_model, p);
	_index->buildIndex();
	obvious::System<double>::deallocate(mData);
}

obvious::Matrix* RandomNormalMatching::pickControlSet(const obvious::Matrix* M,
		vector<unsigned int> idxValid, vector<unsigned int> &idxControl) {
	unsigned int sizeControlSet = _sizeControlSet;
	if ((idxValid.size()) < sizeControlSet) {
		LOGMSG(DBG_DEBUG,
				"Size of scene smaller than control set ... reducing size to " << idxValid.size());
		sizeControlSet = idxValid.size();
	}
	obvious::Matrix* C = new obvious::Matrix(3, sizeControlSet);
	vector<unsigned int> idxTemp = idxValid;
	unsigned int ctr = 0;
	while (idxControl.size() < sizeControlSet) {
		unsigned int r = rand() % idxTemp.size();
		unsigned int idx = idxTemp[r];
		idxControl.push_back(idx);
		idxTemp.erase(idxTemp.begin() + r);

		(*C)(0, ctr) = (*M)(idx, 0);
		(*C)(1, ctr) = (*M)(idx, 1);
		(*C)(2, ctr++) = 1.0;
	}
	return C;
}

void RandomNormalMatching::calcNormals(const Matrix* M, Matrix* N,
		const bool* maskIn, bool* maskOut) {
	int points = M->getRows();

	// mask borders at which we cannot calculate normal vectors
	for (int i = 0; i < _pcaSearchRange / 2; i++)
		maskOut[i] = false;
	for (int i = points - _pcaSearchRange / 2; i < points; i++)
		maskOut[i] = false;

	for (int i = _pcaSearchRange / 2; i < points - _pcaSearchRange / 2; i++) {
		if (maskIn[i]) {
			unsigned int cnt = 0;

			for (int j = -_pcaSearchRange / 2; j < _pcaSearchRange / 2; j++)
				if (maskIn[i + j])
					cnt++;

			if (cnt > _pcaMinSamples) {
				Matrix A(cnt, 2);
				cnt = 0;
				for (int j = -_pcaSearchRange / 2; j < _pcaSearchRange / 2;
						j++) {
					if (maskIn[i + j]) {
						A(cnt, 0) = (*M)(i + j, 0);
						A(cnt, 1) = (*M)(i + j, 1);
						cnt++;
					}
				}
				Matrix* Axes = A.pcaAnalysis();
				// longer axis
				double xLong = (*Axes)(0, 1) - (*Axes)(0, 0);
				double yLong = (*Axes)(0, 3) - (*Axes)(0, 2);
				// shorter axis
				double xShort = (*Axes)(1, 1) - (*Axes)(1, 0);
				double yShort = (*Axes)(1, 3) - (*Axes)(1, 2);
				// rate axes lengths -> main axis needs to be twice as long as second axis
				double lenLongSqr = xLong * xLong + yLong * yLong;
				double lenShortSqr = xShort * xShort + yShort * yShort;

				if (lenShortSqr > 1e-6 && (lenLongSqr / lenShortSqr) < 4.0) {
					maskOut[i] = false;
					continue;
				}

				// shorter axis is normal
				double len = sqrt(lenShortSqr);
				if (((*M)(i, 0) * xShort + (*M)(i, 1) * yShort) < 0.0) {
					(*N)(i, 0) = xShort / len;
					(*N)(i, 1) = yShort / len;
				} else {
					(*N)(i, 0) = -xShort / len;
					(*N)(i, 1) = -yShort / len;
				}

				delete Axes;
			} else
				maskOut[i] = false;
		}
	}
}

void RandomNormalMatching::calcPhi(const Matrix* N, const bool* mask,
		double* phi) {
	if (mask == NULL) {
		for (unsigned int i = 0; i < N->getRows(); i++)
			phi[i] = atan2((*N)(i, 1), (*N)(i, 0));
	} else {
		for (unsigned int i = 0; i < N->getRows(); i++) {
			if (mask[i]) {
				phi[i] = atan2((*N)(i, 1), (*N)(i, 0));
			} else {
				phi[i] = -1e6;
			}
		}
	}
}

void RandomNormalMatching::subsampleMask(bool* mask, unsigned int size,
		double probability) {
	if (probability > 1.0)
		probability = 1.0;
	if (probability < 0.0)
		probability = 0.0;
	int probability_thresh = (int) (1000.0 - probability * 1000.0 + 0.5);
	for (unsigned int i = 0; i < size; i++) {
		if ((rand() % 1000) < probability_thresh) {
			mask[i] = 0;
		}
	}
}

obvious::Matrix RandomNormalMatching::match(const obvious::Matrix* M,
		const bool* maskM, const obvious::Matrix* NM, const obvious::Matrix* S,
		const bool* maskS, double phiMax, const double transMax,
		const double resolution) {
	obvious::Matrix TBest(3, 3);
	TBest.setIdentity();

	const int pointsInM = M->getRows();
	const int pointsInS = S->getRows();

	if (pointsInM != pointsInS) {
		LOGMSG(DBG_ERROR,
				"Model and scene need to be of same size, size of M: " << pointsInM << ", size of S: " << pointsInS);
		return TBest;
	}

	if (pointsInM < 3) {
		LOGMSG(DBG_ERROR,
				"Model and scene contain too less points, size of M: " << pointsInM << ", size of S: " << pointsInS);
		return TBest;
	}

	// ----------------- Model ------------------
	obvious::Matrix* NMpca = new Matrix(pointsInM, 2); // Normals for model
	double* phiM = new double[pointsInM];    // Orientation of model points
	bool* maskMpca = new bool[pointsInM];      // Validity mask of model points

	memcpy(maskMpca, maskM, pointsInM * sizeof(bool));

	if (NM) {
		calcPhi(NM, maskM, phiM);
	} else // if normals are not supplied
	{
		calcNormals(M, NMpca, maskM, maskMpca);
		calcPhi(NMpca, maskMpca, phiM);
	}
	vector<unsigned int> idxMValid = extractSamples(M, maskMpca);

#if USEKNN
	initKDTree(M, idxMValid);
#endif
	// -------------------------------------------

	// ----------------- Scene -------------------
	obvious::Matrix* NSpca = new Matrix(pointsInS, 2); // Normals for scene
	double* phiS = new double[pointsInS];    // Orientation of scene points
	bool* maskSpca = new bool[pointsInS];      // Validity mask of scene points
	memcpy(maskSpca, maskS, pointsInS * sizeof(bool));

	// Determine number of valid samples in local scene neighborhood
	// only from these points a valid orientation is computable
	unsigned int validPoints = 0;
	for (int i = 0; i < pointsInS; i++)
		if (maskSpca[i])
			validPoints++;

	// Probability of point masking
	double probability = 180.0 / (double) validPoints;
	if (probability < 0.99)
		subsampleMask(maskSpca, pointsInS, probability);

	calcNormals(S, NSpca, maskS, maskSpca);
	calcPhi(NSpca, maskSpca, phiS);

	vector<unsigned int> idxSValid = extractSamples(S, maskSpca);
	// -------------------------------------------

	// --------------- Control set ---------------
	vector<unsigned int> idxControl; //represents the indices of points used for Control in S.
	obvious::Matrix* Control = pickControlSet(S, idxSValid, idxControl);
	obvious::Matrix* NControl = new obvious::Matrix(idxControl.size(), 2);
	for (unsigned int i = 0; i < Control->getCols(); i++) {
		(*NControl)(i, 0) = (*NSpca)(idxControl[i], 0);
		(*NControl)(i, 1) = (*NSpca)(idxControl[i], 1);
	}
	unsigned int pointsInC = Control->getCols();
	unsigned int cntMatchThresh = pointsInC / 3; // TODO: Determine meaningful parameter
	double* phiControl = new double[pointsInC]; // Orientation of control points
	calcPhi(NControl, NULL, phiControl);
	// -------------------------------------------//

	// Determine frustum, i.e., direction of leftmost and rightmost model point
	double thetaMin = -((double) pointsInM - 1.0) / 2.0 * resolution; // theoretical bounding
	double thetaBoundMin = atan2((*M)(idxMValid.front(), 1),
			(*M)(idxMValid.front(), 0)); // real bounding
	double thetaBoundMax = atan2((*M)(idxMValid.back(), 1),
			(*M)(idxMValid.back(), 0));  // real bounding

	LOGMSG(DBG_DEBUG,
			"Valid points in scene: " << idxSValid.size() << ", valid points in model: " << idxMValid.size() << ", Control set: " << Control->getCols());
	LOGMSG(DBG_DEBUG,
			"Model phi min:: " << rad2deg(thetaBoundMin) << ", Model phi max: " << rad2deg(thetaBoundMax));

	if (idxSValid.size() < 3) {
		LOGMSG(DBG_ERROR,
				"Too less valid points in scene, matchable size: " << idxSValid.size());
		return TBest;
	}

	if (idxMValid.size() < 3) {
		LOGMSG(DBG_ERROR,
				"Too less valid points in model, matchable size: " << idxMValid.size());
		return TBest;
	}

	// Check for maximum meaningful trials
	unsigned int trials = _trials;
	if (idxMValid.size() < _trials)
		trials = idxMValid.size();

	if (_trace) {
		_trace->reset();
		_trace->setModel(M, idxMValid);
		_trace->setScene(S, idxSValid);
	}

	// Calculate search "radius", i.e., maximum difference in polar indices because of rotation
	phiMax = min(phiMax, M_PI * 0.5);
	int span;
	if (resolution > 1e-6) {
		span = floor(phiMax / resolution);
		if (span > (int) pointsInM)
			span = (int) pointsInM;
	} else {
		LOGMSG(DBG_ERROR,
				"Resolution not properly set: resolution = " << resolution);
		return TBest;
	}

	srand(time(NULL));

	double bestRatio = 0.0;
	unsigned int bestCnt = 0;
	double bestErr = 1e12;

#ifndef DEBUG
	// trace is only possible for single threaded execution
	if (_trace) {
		omp_set_num_threads(1);
		LOGMSG(DBG_WARN,
				"Configured single-threaded execution due to application of trace module");
	}
#endif

	//Timer t;
	//t.start();
	vector<unsigned int> idxTrials = idxMValid;
#pragma omp parallel
	{
		bool* maskControl = new bool[pointsInC];
		double* thetaControl = new double[pointsInC];

#pragma omp for
		for (unsigned int trial = 0; trial < trials; trial++) {

			int idx;
#pragma omp critical
			{
				const int randIdx = rand() % (idxTrials.size());
				idx = idxTrials[randIdx];

				// remove chosen element to avoid picking same index a second time
				idxTrials.erase(idxTrials.begin() + randIdx);
			}

			// leftmost scene point
			const int iMin = max(idx - span, _pcaSearchRange / 2);
			// rightmost scene point
			const int iMax = min(idx + span, pointsInS - _pcaSearchRange / 2);

			for (int i = iMin; i < iMax; i++) {
#if STRUCTAPPROACH
				if(samplesS[i].valid)
#else
				if (maskSpca[i])
#endif
				{

#if STRUCTAPPROACH
					double phi = samplesM[idx].orientation - samplesS[i].orientation;
#else
					double phi = phiM[idx] - phiS[i];
#endif
					if (phi > M_PI)
						phi -= 2.0 * M_PI;
					else if (phi < -M_PI)
						phi += 2.0 * M_PI;

					if (fabs(phi) < phiMax) {
						obvious::Matrix T =
								obvious::MatrixFactory::TransformationMatrix33(
										phi, 0, 0);

						// Calculate translation
						const double sx = (*S)(i, 0);
						const double sy = (*S)(i, 1);
						T(0, 2) = (*M)(idx, 0) - (T(0, 0) * sx + T(0, 1) * sy);
						T(1, 2) = (*M)(idx, 1) - (T(1, 0) * sx + T(1, 1) * sy);

						// Transform control set
						obvious::Matrix STemp = T * (*Control);
						unsigned int pointsInControl = STemp.getCols();

						// Determine number of control points in field of view
						unsigned int maxCntMatch = 0;
						for (unsigned int j = 0; j < pointsInControl; j++) {
							thetaControl[j] = atan2(STemp(1, j), STemp(0, j));
							if (thetaControl[j] > thetaBoundMax
									|| thetaControl[j] < thetaBoundMin) {
								maskControl[j] = false;
							} else {
								maskControl[j] = true;
								maxCntMatch++;
							}
						}

						// Determine how many nearest neighbors (model <-> scene) are close enough
						unsigned int cntMatch = 0;
						flann::Matrix<int> indices(new int[1], 1, 1);
						flann::Matrix<double> dists(new double[1], 1, 1);
						double errSum = 0;
						//double scoreSum = 0.0;

						for (unsigned int s = 0; s < pointsInControl; s++) {
							// clip points outside of model frustum
							if (maskControl[s]) {

#if USEKNN
								// find nearest neighbor of control point
								double q[2];
								q[0] = STemp(0, s);
								q[1] = STemp(1, s);
								flann::Matrix<double> query(q, 1, 2);
								flann::SearchParams p(-1, 0.0);
								_index->knnSearch(query, indices, dists, 1, p);
								const int idxQuery = idxMValid[indices[0][0]];
								double distConsensus = dists[0][0];
#else
								// speeded-up NN search through back projection
								const int idxQuery = round((thetaControl[s]-thetaMin) / resolution);

								if(!maskM[idxQuery]) continue;

								double distX = (*M)(idxQuery, 0) - STemp(0, s);
								double distY = (*M)(idxQuery, 1) - STemp(1, s);
								double distConsensus = distX*distX + distY*distY;
#endif

#if NORMALCONSENSUS
								// Experimental idea: rate matching results additionally with normal consensus
								// consensus score is in range [0, 1] -> perfect match = 0
								double normalConsensus = (1.0
										- cos(
												phiM[idxQuery] - phiControl[s]
														- phi)) / 2.0;
								// Normalized error (weight distance and normal consensus)
								double err = distConsensus * _scaleDistance
										+ normalConsensus * _scaleOrientation;
#else
								double err = distConsensus*_scaleDistance;
#endif

								errSum += err;
								if (err < 1.0)
									cntMatch++;
							}
						}

						delete[] indices.ptr();
						delete[] dists.ptr();

						if (cntMatch <= cntMatchThresh)
							continue;

						// Experimental rating
						double ratio = (double) cntMatch / (double) maxCntMatch;

#pragma omp critical
						{
							// Rating from Markus Kuehn
							double equalThres = 1e-5;
							bool rateCondition = ((ratio - bestRatio)
									> equalThres) && (cntMatch > bestCnt);
							bool similarityCondition = fabs(
									(ratio - bestRatio) < equalThres)
									&& (cntMatch == bestCnt)
									&& errSum < bestErr;
							bool goodMatch = rateCondition
									|| similarityCondition;

							if (goodMatch) {
								bestRatio = ratio;
								bestCnt = cntMatch;
								bestErr = errSum;
								TBest = T;
							}

						}

						if (_trace) {
							//trace is only possible for single threaded execution
							vector<unsigned int> idxM;
							idxM.push_back(idx);
							vector<unsigned int> idxS;
							idxS.push_back(i);
							_trace->addAssignment(M, idxM, S, idxS, &STemp,
									errSum, trial);
						}

					}                // if phiMax
				} // if maskS
			} // for i
		} // for trials

		delete[] maskControl;

	} // OMP

	//cout << "elapsed: " << t.elapsed() << endl;
	//t.reset();

	delete NMpca;
	delete NSpca;
	delete[] phiM;
	delete[] phiS;
	delete[] phiControl;
	delete[] maskMpca;
	delete[] maskSpca;

	delete Control;

	return TBest;
}

obvious::Matrix RandomNormalMatching::match2(const obvious::Matrix* M,
		const bool* maskM, const obvious::Matrix* NM, const obvious::Matrix* S,
		const bool* maskS, double phiMax, const double transMax,
		const double resolution) {
	obvious::Matrix TBest(3, 3);
	TBest.setIdentity();

	const int pointsInM = M->getRows();
	const int pointsInS = S->getRows();

	if (pointsInM != pointsInS) {
		LOGMSG(DBG_ERROR,
				"Model and scene need to be of same size, size of M: " << pointsInM << ", size of S: " << pointsInS);
		return TBest;
	}

	if (pointsInM < 3) {
		LOGMSG(DBG_ERROR,
				"Model and scene contain too less points, size of M: " << pointsInM << ", size of S: " << pointsInS);
		return TBest;
	}

	// ----------------- Model ------------------
	obvious::Matrix* NMpca = new Matrix(pointsInM, 2); // Normals for model
	double* phiM = new double[pointsInM];    // Orientation of model points
	bool* maskMpca = new bool[pointsInM];      // Validity mask of model points

	memcpy(maskMpca, maskM, pointsInM * sizeof(bool));

	if (NM) {
		calcPhi(NM, maskM, phiM);
	} else // if normals are not supplied
	{
		calcNormals(M, NMpca, maskM, maskMpca);
		calcPhi(NMpca, maskMpca, phiM);
	}
	vector<unsigned int> idxMValid = extractSamples(M, maskMpca);

#if USEKNN
	initKDTree(M, idxMValid);
#endif
	// -------------------------------------------

	// ----------------- Scene -------------------
	obvious::Matrix* NSpca = new Matrix(pointsInS, 2); // Normals for scene
	double* phiS = new double[pointsInS];    // Orientation of scene points
	bool* maskSpca = new bool[pointsInS];      // Validity mask of scene points
	memcpy(maskSpca, maskS, pointsInS * sizeof(bool));

	// Determine number of valid samples in local scene neighborhood
	// only from these points a valid orientation is computable
	unsigned int validPoints = 0;
	for (int i = 0; i < pointsInS; i++)
		if (maskSpca[i])
			validPoints++;

	// Probability of point masking
	double probability = 180.0 / (double) validPoints;
	if (probability < 0.99)
		subsampleMask(maskSpca, pointsInS, probability);

	calcNormals(S, NSpca, maskS, maskSpca);
	calcPhi(NSpca, maskSpca, phiS);

	vector<unsigned int> idxSValid = extractSamples(S, maskSpca);
	// -------------------------------------------

	// --------------- Control set ---------------
	vector<unsigned int> idxControl; //represents the indices of points used for Control in S.
	obvious::Matrix* Control = pickControlSet(S, idxSValid, idxControl);
	obvious::Matrix* NControl = new obvious::Matrix(idxControl.size(), 2);
	for (unsigned int i = 0; i < Control->getCols(); i++) {
		(*NControl)(i, 0) = (*NSpca)(idxControl[i], 0);
		(*NControl)(i, 1) = (*NSpca)(idxControl[i], 1);
	}
	unsigned int pointsInC = Control->getCols();
	unsigned int cntMatchThresh = pointsInC * 0.8; // TODO: Determine meaningful parameter
	double* phiControl = new double[pointsInC]; // Orientation of control points
	calcPhi(NControl, NULL, phiControl);
	// -------------------------------------------//

	// Determine frustum, i.e., direction of leftmost and rightmost model point
	double thetaMin = -((double) pointsInM - 1.0) / 2.0 * resolution; // theoretical bounding
	double thetaBoundMin = atan2((*M)(idxMValid.front(), 1),
			(*M)(idxMValid.front(), 0)); // real bounding
	double thetaBoundMax = atan2((*M)(idxMValid.back(), 1),
			(*M)(idxMValid.back(), 0));  // real bounding

	LOGMSG(DBG_DEBUG,
			"Valid points in scene: " << idxSValid.size() << ", valid points in model: " << idxMValid.size() << ", Control set: " << Control->getCols());
	LOGMSG(DBG_DEBUG,
			"Model phi min:: " << rad2deg(thetaBoundMin) << ", Model phi max: " << rad2deg(thetaBoundMax));

	if (idxSValid.size() < 3) {
		LOGMSG(DBG_ERROR,
				"Too less valid points in scene, matchable size: " << idxSValid.size());
		return TBest;
	}

	if (idxMValid.size() < 3) {
		LOGMSG(DBG_ERROR,
				"Too less valid points in model, matchable size: " << idxMValid.size());
		return TBest;
	}

	// Check for maximum meaningful trials
	unsigned int trials = _trials;
	if (idxMValid.size() < _trials)
		trials = idxMValid.size();

	if (_trace) {
		_trace->reset();
		_trace->setModel(M, idxMValid);
		_trace->setScene(S, idxSValid);
	}

	// Calculate search "radius", i.e., maximum difference in polar indices because of rotation
	phiMax = min(phiMax, M_PI * 0.5);
	int span;
	if (resolution > 1e-6) {
		span = floor(phiMax / resolution);
		if (span > (int) pointsInM)
			span = (int) pointsInM;
	} else {
		LOGMSG(DBG_ERROR,
				"Resolution not properly set: resolution = " << resolution);
		return TBest;
	}

	srand(time(NULL));

	double bestRatio = 0.0;
	unsigned int bestCnt = 0;
	double bestErr = 1e-200;

#ifndef DEBUG
	// trace is only possible for single threaded execution
	if (_trace) {
		omp_set_num_threads(1);
		LOGMSG(DBG_WARN,
				"Configured single-threaded execution due to application of trace module");
	}
#endif

//Timer t;
//t.start();
	vector<unsigned int> idxTrials = idxMValid;

	bool* maskControl = new bool[pointsInC];
	double* thetaControl = new double[pointsInC];

	// ToDo: create angle array from model

	std::vector<double> anglesModel;
	std::vector<double> distModel;

	for(unsigned int i = 0; i < idxMValid.size(); i++){
		//double phi;
		anglesModel.push_back(atan2((*M)(idxMValid[i], 1), (*M)(idxMValid[i], 0)));
		distModel.push_back( sqrt( pow(((*M)(idxMValid[i], 0)), 2) + pow(((*M)(idxMValid[i],1)), 2) ) );
	}


	for (unsigned int trial = 0; trial < trials; trial++) {

		int idx;
#pragma omp critical
		{
			const int randIdx = rand() % (idxTrials.size());
			idx = idxTrials[randIdx];

			// remove chosen element to avoid picking same index a second time
			idxTrials.erase(idxTrials.begin() + randIdx);
		}

		// leftmost scene point
		const int iMin = max(idx - span, _pcaSearchRange / 2);
		// rightmost scene point
		const int iMax = min(idx + span, pointsInS - _pcaSearchRange / 2);

		for (int i = iMin; i < iMax; i++) {

#if STRUCTAPPROACH
			if(samplesS[i].valid)
#else
			if (maskSpca[i])
#endif
			{

#if STRUCTAPPROACH
				double phi = samplesM[idx].orientation - samplesS[i].orientation;
#else
				double phi = phiM[idx] - phiS[i];
#endif
				if (phi > M_PI)
					phi -= 2.0 * M_PI;
				else if (phi < -M_PI)
					phi += 2.0 * M_PI;

				if (fabs(phi) < phiMax) {
					obvious::Matrix T =
							obvious::MatrixFactory::TransformationMatrix33(phi,
									0, 0);

					// Calculate translation
						// todo: dan: trying to use phi + dx, dy directly for transformation via rightsided multiplication!?
					const double sx = (*S)(i, 0);
					const double sy = (*S)(i, 1);
					T(0, 2) = (*M)(idx, 0) - (T(0, 0) * sx + T(0, 1) * sy);
					T(1, 2) = (*M)(idx, 1) - (T(1, 0) * sx + T(1, 1) * sy);

					// Transform control set
					obvious::Matrix STemp = T * (*Control);
					unsigned int pointsInControl = STemp.getCols();

					// Determine number of control points in field of view
					unsigned int maxCntMatch = 0;
					for (unsigned int j = 0; j < pointsInControl; j++) {
						thetaControl[j] = atan2(STemp(1, j), STemp(0, j));
						if (thetaControl[j] > thetaBoundMax || thetaControl[j] < thetaBoundMin) {
							maskControl[j] = false;
						} else {
							maskControl[j] = true;
							maxCntMatch++;
						}
					}

					std::vector<double> probOfSingleScans;

					if (maxCntMatch > cntMatchThresh){

						// Rating dan_tob
						for (unsigned int s = 0; s < pointsInControl; s++) { // whole control set
							// clip points outside of model frustum
							if (maskControl[s]) { // if point is in field of view

								double angle = atan2((STemp)(1, s), (STemp)(0, s));
								double distance = sqrt( ((STemp)(0, s)) * ((STemp)(0, s)) + ((STemp)(1, s)) * ((STemp)(1, s)) );

								double minAngleDiff = 2 * M_PI;
								int idxMinAngleDiff;
								double diff;

								for(int i = 0; i < anglesModel.size(); i++){ // find right model point to control point
									diff = abs(angle - anglesModel[i]);
									if ( diff < minAngleDiff ){ // find min angle
										minAngleDiff = diff;
										idxMinAngleDiff = i;
									}
								}

								//cout <<  "min angle " << minAngleDiff << endl;


								double probOfSingleScan;
								probOfSingleScan = probabilityOfTwoSingleScans(distModel[idxMinAngleDiff], distance);
								probOfSingleScans.push_back(probOfSingleScan);

								//cout << "angle model|scene: " << anglesModel[idxMinAngleDiff] * 180.0 / M_PI << " | " <<  angle * 180.0 / M_PI  <<
								//		"; dist model|scene: " << distModel[idxMinAngleDiff] << " | " << distance << " prob: "<< probOfSingleScan << endl;

							}// if point is in field of view
						} // whole control set

						double probOfActualScan = 1;

						for(int i = 0; i < probOfSingleScans.size(); i++){
							probOfActualScan = probOfActualScan * probOfSingleScans[i];
						}


						if(probOfActualScan > bestErr){
							TBest = T;
							bestErr = probOfActualScan;
							//cout << "new errSum: " << probOfActualScan << " trial: " << trial << endl;

							if (_trace) {
								//trace is only possible for single threaded execution
								vector<unsigned int> idxM;
								idxM.push_back(idx);
								vector<unsigned int> idxS;
								idxS.push_back(i);
								_trace->addAssignment(M, idxM, S, idxS, &STemp, 10e100 * probOfActualScan,
										trial);
							}
						}
					}// if phiMax
				} // if maskS
			}
		} // for i
	} // for trials

	delete[] maskControl;

	//cout << "elapsed: " << t.elapsed() << endl;
	//t.reset();

	delete NMpca;
	delete NSpca;
	delete[] phiM;
	delete[] phiS;
	delete[] phiControl;
	delete[] maskMpca;
	delete[] maskSpca;

	delete Control;

	return TBest;
}

double RandomNormalMatching::probabilityOfTwoSingleScans(double m, double s) {

	// hit
	if (s < _rangemax) {
		_phit = (1) / (sqrt(2 * M_PI * pow(_sighit, 2))) * pow(M_E, ((-0.5 * pow((m - s), 2)) / (pow(_sighit, 2)))) ;
	}

	// short
	if (s < m){
	double n = (1)/(1 - pow(M_E, (-_lamshort * m)));
	_pshort = n * _lamshort * pow(M_E, (-_lamshort * s)) ;
	}

	// max
	if (s >= _rangemax){
	_pmax = 1;
	}

	// rand
	if (s < _rangemax){
	_prand = 1 / _rangemax;
	}

	return _zhit * _phit + _zshort * _pshort + _zmax * _pmax + _zrand * _prand;
}


void RandomNormalMatching::serializeTrace(const char* folder) {
	if (_trace)
		_trace->serialize(folder);
	else
		LOGMSG(DBG_ERROR, "Trace not activated");
}

}
