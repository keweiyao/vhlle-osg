#pragma once
#include <cmath>
#include <string>
#include <iostream>
//class TGraph ;

class EoS
{
public :
	virtual ~EoS() { return ; } // why did I do this?
	virtual void eos(double e, double nb, double nq, double ns,
		double &T, double &mub, double &muq, double &mus, double &p) = 0 ;
	virtual double p(double e, double nb, double nq, double ns) = 0;
	double s(double e, double nb, double nq, double ns) ;
	double inline cs2(void) {return 0.333333;};
	double inline cs(void) {return 0.57735;};
	virtual double cs2(double e) = 0;
	virtual void gete(double s, double& e, double nb) = 0;
};

