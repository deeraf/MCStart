/*
 * Logger.h
 *
 *  Created on: May 18, 2013
 *      Author: amifan
 */

#ifndef LOGGER_H_
#define LOGGER_H_

//#include <mgl.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "string.h"
#include "Particle.h"
#include "ParticleClass.h"
#include "ParticleClassVector.h"

using namespace std;

class Logger {
public:
	Logger();
	virtual ~Logger();
	void DrawPlot(double *x, double *qx, int n,int name, double time);
	void WriteToFileParticle(ParticleClassVector pcv, double volumen, int name);
	void PrintTraceCoagulation(double numberOfParticles, double numberOfCoagulation, double ithDiameter, double jthDiameter, Particle newParticle, double nextCriterion,
			double nextRandom, double S, double currentTime);
	void PrintTraceKeimbildung(double numberOfParticles, double numberOfKeimbildung, double diameter, double zusammensetzung, double keimbildungsrate, double currentTime);
	void WriteToFileCoagulation(double numberOfParticles, double numberOfCoagulation, double ithDiameter, double jthDiameter, Particle newParticle, double nextCriterion,
			double nextRandom, double S, double currentTime);
	void WriteToFileKeimbildung (double numberOfParticles, double numberOfKeimbildung, double diameter, double zusammensetzung, double keimbildungsrate, double currentTime);
	void WriteToFileWachstum(double deltaWasserDampMasse, double deltaGlycerinDampMasse, double wasserDampMasse, double glycerinDampMasse, double currentTime);
	void WriteToFileKeimbildungMasse(double deltaWasserDampMasse, double deltaGlycerinDampMasse, double wasserDampMasse, double glycerinDampMasse, double currentTime);
	void WriteToFileMittelWert(double zusammensetzung, double diameter, double numberOfAllParticles, double volumen, double currentTime);
	void WriteToFile3D(ParticleClassVector pcv, int name);
	void WriteToFileDebug(ParticleClassVector pcv, int name);
	void WriteToFileParticleCount(ParticleClassVector pcv, double volumen, int name);


private:
	char *GetCharString(double i, string extension);
	char *GetCharString(double i);
};

#endif  /*LOGGER_H_*/
