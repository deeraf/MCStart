/*
 * LogNormalDistribution.h
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */

#ifndef LOGNORMALDISTRIBUTION_H_
#define LOGNORMALDISTRIBUTION_H_

#include "ParticleClassVector.h"
#include "ParticleSystem.h"


class LogNormalDistribution {
public:
	LogNormalDistribution() {
	}
	LogNormalDistribution(double sigma, double x50, double tolerant, double N, double c,double wasserDampMasse,
			double glycerinDampMasse, double stickstoffDampMasse);
	void SetParameter(double sigma, double x50, double n);
	double q(double x, double x50);
	ParticleClassVector Initialize();
private:
	double sigma;
	double totalNumberOfParticles;
	double x50;
	double maxNumberOfParticlesInOneClass;
	double tolerant;
	ParticleSystem ps;
};

#endif /* LOGNORMALDISTRIBUTION_H_ */
