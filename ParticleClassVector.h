/*
 * ParticleVector.h
 *
 *  Created on: Apr 10, 2013
 *      Author: amifan
 */

#ifndef PARTICLEVECTOR_H_
#define PARTICLEVECTOR_H_

#include <vector>
#include "ParticleClass.h"
#include "ParticleSystem.h"

using namespace std;

class ParticleClassVector {
public:
	ParticleClassVector();
	ParticleClassVector(double interval, double tolerant);
	virtual ~ParticleClassVector();
	ParticleClass& GetParticleClassAt(int i);
	double *GetXResultAsInterval();
	double *GetXResultAsIntervalNew();
	double *GetQxResultAsInterval();
	double *GetQxResultAsIntervalNew();
	double GetNumberOfAllParticles();
	double GetTotalRate();
	double CalK(double x1, double x2);
	Particle UpdateParticleAndRate(unsigned int i, unsigned int j, double rate1, double rate2, double criterion);
	unsigned int GetNumberOfParticleClasses();
	unsigned int AddParticleAndUpdateRate(Particle p);
	unsigned int GetNumberOfParticles();
	vector<Particle> GetAllParticles();
	vector<Particle>::iterator Erase(vector<ParticleClass>::iterator i, vector<ParticleClass>::iterator j);
	void InsertNewParticleClass();
	void InitializeCumulativeRate();
	void Clear();
	void ReduceParticles();
	void DuplicateParticles();
	static bool CompareByAverageDiameter(const ParticleClass& x, const ParticleClass& y);
	static bool CompareByDiameter(const Particle& x, const Particle& y);
	static bool CompareByCumulativeRate(const ParticleClass& x, const ParticleClass& y);
	double GetMittelDurchmesser();
	double GetMittelZusammensetzung();
	void AddParticle(int i, double p);
	void AddParticle(int i, double p, double zusammensetzung, double kritischMolWasser, double kritischMolGlycerin, unsigned int numberOfWasser, unsigned int numberOfGlycerin);
	void UpdateBoundary();



protected:
	vector<ParticleClass> pcv;
	vector<ParticleClass>::iterator it;
	ParticleSystem ps;
	double *xResult;
	double *qxResult;
	double totalRate;
	double tolerant;
	double interval;
	unsigned int totalNumberOfParticles;
	double mittelDurchmesser;
	double mittelZusammensetzung;

	double rightBoundary;
	double lastClassRange;
};

#endif /* PARTICLEVECTOR_H_ */
