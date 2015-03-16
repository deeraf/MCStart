/*
 * Calculator.h
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */

#ifndef CALCULATOR_H_
#define CALCULATOR_H_

#include <stdlib.h>
#include <vector>
#include "ParticleClassVector.h"
#include "ParticleSystem.h"
#include "Logger.h"
#include "LogNormalDistribution.h"

class Calculator {
public:
	Calculator(double simulationTime, double tolerant, unsigned int maxNumberOfParticles, double interval, double wasserDampMasse, double glycerinDampMasse, double stickstoffDampMasse);
	Calculator(LogNormalDistribution logNormalDistribution, double simulationTime, double tolerant, unsigned int maxNumberOfParticles, double interval, double wasserDampMasse, double glycerinDampMasse, double stickstoffDampMasse);
	virtual ~Calculator();
	double CalK(double x1, double x2);
	double CalVolumen(unsigned int max);
	double GetDeltaTime(double keimbildungsRate, double coagulationsRate);
	double qRandom();
	void DoSimulation();
	void DoWachstum(double deltaT, double currentT, double vol);
	void DoOnlyWachstumForOnlyWachtum(double deltaT, double currentT, double vol);
	void DoOnlyWachstum();
private:
	ParticleClassVector pcv;
	ParticleSystem ps;
	Logger logger;
	double volumen;
	double totalSimulationTime;
	double tolerant;
	unsigned int maxNumberOfParticles;
};

#endif /* CALCULATOR_H_ */
