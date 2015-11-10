/**
 * Sample application showing the usage of the OBVIOUS 2D range image implementation.
 * @author Stefan May
 * @date 23.09.2014
 */

#include <string.h>
#include <iostream>
#include <math.h>
#include <fstream>

#include "obcore/math/mathbase.h"
#include "obcore/math/linalg/linalg.h"
#include "obcore/base/Timer.h"

#include "obvision/registration/ransacMatching/RansacMatching.h"
#include "obvision/registration/ransacMatching/RandomNormalMatching.h"

using namespace std;
using namespace obvious;

#define DATASETSIZE 1081

int main(int argc, char** argv) {
	Timer timer;
	timer.start();

	char* sz;
	std::vector<double> x_m;
	std::vector<double> y_m;
	std::vector<double> m_m;
	std::vector<double> x_s;
	std::vector<double> y_s;
	std::vector<double> m_s;

	//load data
	string line;

	cout << "open model file" << endl;
	ifstream model("/home/amndan/Schreibtisch/maske/01/449_1420815145_10.4716_match2/rawData.dat");
	if (model.is_open()) {
		while (getline(model, line)) {
			//cout << line << '\n';
			const char * cline = line.c_str();
			x_m.push_back(std::strtod(cline, &sz));
			y_m.push_back(std::strtod(sz + 1, &sz));
			m_m.push_back((bool)std::strtod(sz + 1, &sz));
			x_s.push_back(std::strtod(sz + 1, &sz));
			y_s.push_back(std::strtod(sz + 1, &sz));
			m_s.push_back((bool)std::strtod(sz + 1, &sz));
			cout << "x: " << x_m[x_m.size()-1] << ", y: " << y_m[y_m.size()-1] << ", m: " << m_m[m_m.size()-1] << endl;
		}
		model.close();
	} else
		cout << "Unable to open file";

	cout << x_s.size() << endl;
	cout << x_m.size() << endl;


	// Model coordinates
	obvious::Matrix* M = new obvious::Matrix(x_m.size(), 2);
	obvious::Matrix* S = new obvious::Matrix(x_s.size(), 2);
	bool maskM[x_m.size()];
	bool maskS[x_s.size()];

	for (unsigned int i = 0; i < x_m.size(); i++) {

		(*M)(i, 0) = x_m[i];
		(*M)(i, 1) = y_m[i];
		(*S)(i, 0) = x_s[i];
		(*S)(i, 1) = y_s[i];

		maskM[i] = m_m[i];
		maskS[i] = m_s[i];

		cout << "masks: " << "m: " << m_m[i] << " s: " << m_s[i] << endl;

//		maskM[i] = true;
//		maskS[i] = true;
	}

	//Model Normals
//  obvious::Matrix* N = new obvious::Matrix(DATASETSIZE, 2);
//  // compute mean of components build by left and right neighbors
//  for(int i=1; i<DATASETSIZE-1; i++)
//  {
//    double xleft  = (*M)(i, 0)   - (*M)(i-1, 0);
//    double xright = (*M)(i+1, 0) - (*M)(i, 0);
//    double yleft  = (*M)(i, 1)   - (*M)(i-1, 1);
//    double yright = (*M)(i+1, 1) - (*M)(i, 1);
//
//
//    // x-component of normal
//    double xNormal = -(yright + yleft) * 0.5;
//    // y-component of normal
//    double yNormal = (xright + xleft) * 0.5;
//    //Normalize
//    double lengthNormal = sqrt(xNormal*xNormal + yNormal*yNormal);
//    (*N)(i-1, 0) = xNormal / lengthNormal;
//    (*N)(i-1, 1) = yNormal / lengthNormal;
//  }
//
//  // left bound
//  (*N)(0, 0) = -((*M)(1, 1) - (*M)(0, 1));
//  (*N)(0, 1) = (*M)(1, 0) - (*M)(0, 0);
//
//  // right bound
//  (*N)(DATASETSIZE-1, 0) = -((*M)(DATASETSIZE-1, 1) - (*M)(DATASETSIZE-2, 1));
//  (*N)(DATASETSIZE-1, 1) = (*M)(DATASETSIZE-1, 0) - (*M)(DATASETSIZE-2, 0);


	//for(unsigned int i=0; i<S.getRows(); i++)
	//  cout << S(i, 0) << endl;



	unsigned int trials = 20;
	double epsThresh = 0.02;
	unsigned int sizeControlSet = 100;



	//RansacMatching matcher(trials, epsThresh, sizeControlSet);
	 RandomNormalMatching matcher(trials, epsThresh, sizeControlSet);
	  matcher.activateTrace();

	//Matrix F = matcher.match(M, maskM, &S, maskS, deg2rad(45.0), 1.5 , deg2rad(0.25));
	Matrix F = matcher.match2(M, maskM, NULL, S, maskS, deg2rad(45.0), 1.5, deg2rad(0.25));

	matcher.serializeTrace("/tmp/trace2");

	F.invert();
	cout << endl << "Estimated transformation:" << endl;
	F.print();

	cout << "elapsed: " << timer.elapsed() << " s" << endl;
	return 0;
}
