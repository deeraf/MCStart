/*
 * LogNormalDistribution.cpp
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */

#include "LogNormalDistribution.h"
#include "Particle.h"
#include "ParticleClassVector.h"
#include "ParticleSystem.h"
#include "Logger.h"
#include "math.h"
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

LogNormalDistribution::LogNormalDistribution(double s, double x, double t, double N, double c,double wasserDampMasse,
		double glycerinDampMasse, double stickstoffDampMasse) {
	sigma = s;
	x50 = x;
	tolerant = t;
	totalNumberOfParticles = N;
	maxNumberOfParticlesInOneClass = c;
	ps = ParticleSystem(wasserDampMasse, glycerinDampMasse, stickstoffDampMasse);
}

double LogNormalDistribution::q(double x, double x50) {
	double result = 0.0;

	double a = sqrt(2 * M_PI) * log(sigma);
	double b = log(x / x50) / log(sigma);

	if (x != 0) {
		result = 1 / a * 1 / x * exp(-0.5 * pow(b, 2));
	} else {
		result = 0;
	}

	return result;
}

ParticleClassVector LogNormalDistribution::Initialize() {
	ParticleClassVector pcv;

	double deltaX = 1E-15;
	double xi = 0;
	double interval = 1 / totalNumberOfParticles;
	double currentR = 0;
	double nextR = 0;
	double *xiArray = (double*) malloc(sizeof(double) * totalNumberOfParticles);
	double diff1;
	double diff2;

	for (int i = 0; i < totalNumberOfParticles; i++) {
		nextR = (i + 0.5) * interval;

		while (currentR < nextR) {
			xi += deltaX;
			currentR += q(xi - deltaX * 0.5, x50) * deltaX;
		}

		xiArray[i] = xi;
		cout << "deltaA: " << i << " xi " << xi << endl;
	}

	do {
		int countOfParticleInClass = 0;
		int index = 0;

		pcv = ParticleClassVector(maxNumberOfParticlesInOneClass, tolerant);

		pcv.InsertNewParticleClass();

		for (int j = 0; j < totalNumberOfParticles; j++) {
			if (countOfParticleInClass < maxNumberOfParticlesInOneClass) {
				pcv.AddParticle(index, xiArray[j]);
				countOfParticleInClass++;
			} else {
				pcv.InsertNewParticleClass();
				index++;
				double durchmesser = xiArray[j];
				unsigned int numberOfWasser = ps.CalWasserNumber(durchmesser);
				unsigned int numberOfGlycerin = ps.CalGlycerinNumber(durchmesser);
				double kritischMolWasser = ps.CalkritischeMolWasser();
				double kritischMolGlycerin = ps.CalkritischeMolGlyzerin();
				double zusammensetzung = ps.CalKritischGlycerinZusammensetzung();
//				double zusammensetzung = 0.0;
				pcv.AddParticle(index, xiArray[j], zusammensetzung, kritischMolWasser, kritischMolGlycerin, numberOfWasser, numberOfGlycerin);
				pcv.AddParticle(index, xiArray[j]);
				countOfParticleInClass = 1;
			}
		}

		for (int k = 0; k < pcv.GetNumberOfParticleClasses(); k++) {
			pcv.GetParticleClassAt(k).CalAverageDiameter();
		}

		pcv.UpdateBoundary();

		diff1 = (pcv.GetParticleClassAt(0).GetMaxDiameter() - pcv.GetParticleClassAt(0).GetMinDiameter()) / pcv.GetParticleClassAt(0).GetAverageDiameter();
		diff2 = (pcv.GetParticleClassAt(pcv.GetNumberOfParticleClasses() - 1).GetMaxDiameter() - pcv.GetParticleClassAt(pcv.GetNumberOfParticleClasses() - 1).GetMinDiameter())
				/ pcv.GetParticleClassAt(pcv.GetNumberOfParticleClasses() - 1).GetAverageDiameter();

		maxNumberOfParticlesInOneClass -= 1;
	} while (diff1 > tolerant || diff2 > tolerant);

	for (int l = 0; l < pcv.GetNumberOfParticleClasses(); l++) {
		for (int m = 0; m < pcv.GetParticleClassAt(l).GetSize(); m++) {
			pcv.GetParticleClassAt(l).GetParticleAt(m).SetWasserMol(ps.CalkritischeMolWasser());
			pcv.GetParticleClassAt(l).GetParticleAt(m).SetGlycerinMol(ps.CalkritischeMolGlyzerin());
			cout<<pcv.GetParticleClassAt(l).GetParticleAt(m).GetWasserMol()<<"  "<<pcv.GetParticleClassAt(l).GetParticleAt(m).GetGlycerinMol()<<endl;
		}
	}

	//Logger log;
	//log.GetMaxNumberOfParticlesAccordingTotTolerant(++maxNumberOfParticlesInOneClass);

	return pcv;
}
