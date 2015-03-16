/*
 * ParticleClass.h
 *
 *  Created on: Jul 4, 2013
 *      Author: amifan
 */

#ifndef PARTICLECLASS_H_
#define PARTICLECLASS_H_

#include "Particle.h"
#include "string.h"
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

class ParticleClass {

public:
	ParticleClass();
	virtual ~ParticleClass();
	double GetCumulativeRate();
	double GetAverageDiameter();
	double GetLeftBoundary();
	double GetRightBoundary();
	double GetDiameterAt(unsigned int index);
	unsigned int GetSize();
	Particle& GetParticleAt(unsigned int index);
	static bool CompareByDiameter(const Particle& x, const Particle& y);
	void SetCumulativeRate(double cumulativeRate);
	void UpdateCumulativeRate(double deltaRate);
	void PushBackParticle(Particle particle);
	void CalAverageDiameter();
	void SetAverageDiameter(double diameter);
	void UpdateAverageDiameter(double diameter);
	void SetLeftBoundary(double boundary);
	void SetRightBoundary(double boundary);
	void EraseParticleAt(unsigned int index);
	void SortMe();
	double GetMinDiameter();
	double GetMaxDiameter();

	double averageDiameter;
	double cumulativeRate;
private:
	vector<Particle> pv;
	vector<Particle>::iterator it;
	double leftBoundary;
	double rightBoundary;
};

#endif /* PARTICLECLASS_H_ */
