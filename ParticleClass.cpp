/*
 * ParticleClass.cpp
 *
 *  Created on: Jul 4, 2013
 *      Author: amifan
 */

#include "ParticleClass.h"
#include "string.h"
#include <algorithm>

using namespace std;

ParticleClass::ParticleClass() {
	averageDiameter = 0.0;
	cumulativeRate = 0.0;
	leftBoundary = -1.0;
	rightBoundary = -1.0;
}

ParticleClass::~ParticleClass() {

}

bool ParticleClass::CompareByDiameter(const Particle& x, const Particle& y) {
	return x.diameter < y.diameter;
}
;

double ParticleClass::GetCumulativeRate() {
	return cumulativeRate;
}

void ParticleClass::UpdateCumulativeRate(double deltaRate) {
	cumulativeRate += deltaRate;
}

void ParticleClass::SetCumulativeRate(double rate) {
	cumulativeRate = rate;
}

double ParticleClass::GetAverageDiameter() {
	return averageDiameter;
}

void ParticleClass::CalAverageDiameter() {
	unsigned int count = pv.size();
	double temp = 0.0;

	for (unsigned int i = 0; i < count; i++) {
		temp += pv[i].GetDiameter();
	}

	if (count == 0) {
		averageDiameter = 0;
	} else {
		averageDiameter = temp / count;
	}
}

void ParticleClass::UpdateAverageDiameter(double newDiameter) {
	averageDiameter = (averageDiameter * (pv.size() - 1) + newDiameter) / pv.size();
}

void ParticleClass::SetAverageDiameter(double diameter) {
	averageDiameter = diameter;
}

unsigned int ParticleClass::GetSize() {
	return pv.size();
}

void ParticleClass::PushBackParticle(Particle particle) {
	pv.push_back(particle);
}

double ParticleClass::GetLeftBoundary() {
	return leftBoundary;
}

void ParticleClass::SetLeftBoundary(double boundary) {
	leftBoundary = boundary;
}

double ParticleClass::GetRightBoundary() {
	return rightBoundary;
}

void ParticleClass::SetRightBoundary(double boundary) {
	rightBoundary = boundary;
}

void ParticleClass::EraseParticleAt(unsigned int i) {

	if (i == pv.size()) {
		i--;
	}
	double deltaDiameter = pv.at(i).GetDiameter();

	if (pv.size() > 1) {
		averageDiameter = (averageDiameter * pv.size() - deltaDiameter) / (pv.size() - 1);
	} else {
		averageDiameter = 0.0;
	}
	it = pv.begin();
	pv.erase(it + i);
}

double ParticleClass::GetDiameterAt(unsigned int i) {
	return pv.at(i).GetDiameter();
}

Particle& ParticleClass::GetParticleAt(unsigned int i) {
	return pv.at(i);
}

void ParticleClass::SortMe() {
	sort(pv.begin(), pv.end(), CompareByDiameter);
}

double ParticleClass::GetMinDiameter() {
	return pv.at(0).GetDiameter();
}

double ParticleClass::GetMaxDiameter() {
	return pv.at(pv.size() - 1).GetDiameter();
}

