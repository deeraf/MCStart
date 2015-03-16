/*
 * ParticleVector.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: amifan
 */

#include "ParticleClassVector.h"
#include "ParticleClass.h"
#include "Particle.h"
#include "Logger.h"
#include "math.h"
#include <algorithm>
#include <iostream>

using namespace std;

ParticleClassVector::ParticleClassVector() {

}

ParticleClassVector::ParticleClassVector(double classInterval, double t) {
	totalNumberOfParticles = 0;
	interval = classInterval;
	ps = ParticleSystem();
}

ParticleClassVector::~ParticleClassVector() {

}

bool ParticleClassVector::CompareByAverageDiameter(const ParticleClass& x, const ParticleClass& y) {
	return x.averageDiameter < y.averageDiameter;
}
;

bool ParticleClassVector::CompareByDiameter(const Particle& x, const Particle& y) {
	return x.diameter < y.diameter;
}
;

bool ParticleClassVector::CompareByCumulativeRate(const ParticleClass& x, const ParticleClass& y) {
	return x.cumulativeRate < y.cumulativeRate;
}
;

void ParticleClassVector::InsertNewParticleClass() {
	pcv.push_back(ParticleClass());
}

void ParticleClassVector::AddParticle(int i, double p) {
	pcv[i].PushBackParticle(Particle(p));
}

void ParticleClassVector::AddParticle(int i, double p, double zusammensetzung, double kritischMolWasser, double kritischMolGlycerin, unsigned int numberOfWasser, unsigned int numberOfGlycerin){
	pcv[i].PushBackParticle(Particle(p, zusammensetzung, kritischMolWasser, kritischMolGlycerin, numberOfWasser, numberOfGlycerin));
}


unsigned int ParticleClassVector::AddParticleAndUpdateRate(Particle p) {
	bool foundClass = false;
	unsigned int index = 0;
	double oldDiameter = 0.0;
	double deltaRate = 0.0;

	for (unsigned int i = 0; i < pcv.size(); i++) {
		if (p.GetDiameter() >= pcv[i].GetLeftBoundary() && p.GetDiameter() <= pcv[i].GetRightBoundary()) {
			pcv[i].PushBackParticle(p);
			oldDiameter = pcv[i].GetAverageDiameter();
			pcv[i].UpdateAverageDiameter(p.GetDiameter());
			totalNumberOfParticles++;
			foundClass = true;
			index = i;
			break;
		}
	}

	if (!foundClass) {
		int n = p.GetDiameter() / interval;
		InsertNewParticleClass();
		index = pcv.size() - 1;
		pcv[index].SetLeftBoundary(n * interval);
		pcv[index].SetRightBoundary((n + 1) * interval);
		pcv[index].PushBackParticle(p);
		oldDiameter = pcv[index].GetAverageDiameter();
		pcv[index].UpdateAverageDiameter(p.GetDiameter());
		totalNumberOfParticles++;
	}

	double temp1 = 0.0;
	double temp2 = 0.0;

	if (index == 0) {
		temp1 = pcv[0].GetCumulativeRate();
	} else {
		temp1 = pcv[index].GetCumulativeRate() - pcv[index - 1].GetCumulativeRate();
	}

	temp2 = CalK(pcv[index].GetAverageDiameter(), pcv[index].GetAverageDiameter()) * pcv[index].GetSize()
			* (pcv[index].GetSize() - 1) * 0.5;

	for (unsigned int k = index + 1; k < pcv.size(); k++) {
		temp2 += CalK(pcv[index].GetAverageDiameter(), pcv[k].GetAverageDiameter()) * pcv[index].GetSize()
				* pcv[k].GetSize();
	}

	for (unsigned int j = 0; j < pcv.size(); j++) {
		if (j < index) {
			deltaRate += CalK(pcv[index].GetAverageDiameter(), pcv[j].GetAverageDiameter()) * pcv[index].GetSize()
					* pcv[j].GetSize()
					- CalK(oldDiameter, pcv[j].GetAverageDiameter()) * (pcv[index].GetSize() - 1) * pcv[j].GetSize();
		} else if (j == index) {
			deltaRate += temp2 - temp1;
		}

		pcv[j].UpdateCumulativeRate(deltaRate);
	}

	return index;
}

ParticleClass& ParticleClassVector::GetParticleClassAt(int i) {
	if (i == -1) {
		return pcv[pcv.size() - 1];
	} else {
		return pcv[i];
	}
}

double *ParticleClassVector::GetXResultAsInterval() {
	vector<Particle> temp = GetAllParticles();

	double minimal = temp[0].GetDiameter();
	double maximal = temp[temp.size() - 1].GetDiameter();
	int interval = (maximal - minimal) / 1E-10;

	xResult = (double*) malloc(sizeof(double) * interval);

	for (int i = 0; i < interval; i++) {
		xResult[i] = (minimal + i * 1E-10) * 1E9;
	}

	return xResult;
}

double *ParticleClassVector::GetXResultAsIntervalNew() {
	int count = GetNumberOfParticleClasses();

	xResult = (double*) malloc(sizeof(double) * count);

	sort(pcv.begin(), pcv.end(), CompareByAverageDiameter);

	xResult[0] = 0;
	for (int i = 0; i < count; i++) {
		xResult[i] = pcv[i].GetAverageDiameter() * 1E9;
	}

	return xResult;
}

double *ParticleClassVector::GetQxResultAsInterval() {
	vector<Particle> temp = GetAllParticles();

	double minimal = temp[0].GetDiameter();
	double maximal = temp[temp.size() - 1].GetDiameter();
	int interval = (maximal - minimal) / 1E-10;

	qxResult = (double*) malloc(sizeof(double) * interval);

	for (int i = 0; i < interval; i++) {
		int n = 0;
		for (unsigned int j = 0; j < temp.size(); j++) {
			if (temp[j].GetDiameter() > minimal + i * 1E-10 && temp[j].GetDiameter() <= minimal + (i + 1) * 1E-10) {
				n += 1;
			}
		}
		qxResult[i] = n;
	}

	return qxResult;
}

double *ParticleClassVector::GetQxResultAsIntervalNew() {
	int count = GetNumberOfParticleClasses();

	qxResult = (double*) malloc(sizeof(double) * count);

	sort(pcv.begin(), pcv.end(), CompareByAverageDiameter);

	for (int i = 0; i < count; i++) {
		qxResult[i] = pcv[i].GetSize();
	}

	return qxResult;
}

double ParticleClassVector::GetNumberOfAllParticles() {
	double size = 0;

	for (unsigned int i = 0; i < pcv.size(); i++) {
		size += pcv[i].GetSize();
	}

	return size;
}

unsigned int ParticleClassVector::GetNumberOfParticleClasses() {
	return pcv.size();
}

unsigned int ParticleClassVector::GetNumberOfParticles() {
	return totalNumberOfParticles;
}

double ParticleClassVector::GetTotalRate() {
	if (pcv.size() > 0) {
		return pcv[pcv.size() - 1].GetCumulativeRate();
	} else {
		return 0.0;
	}
}

void ParticleClassVector::InitializeCumulativeRate() {
	double rate = 0.0;
	double cumulativeRate = 0.0;

	for (unsigned int k = 0; k < pcv.size(); k++) {
		pcv[k].CalAverageDiameter();
	}

	for (unsigned int i = 0; i < pcv.size(); i++) {
		cumulativeRate += CalK(pcv[i].GetAverageDiameter(), pcv[i].GetAverageDiameter()) * pcv[i].GetSize()
				* (pcv[i].GetSize() - 1) * 0.5;

		for (unsigned int j = i + 1; j < pcv.size(); j++) {
			rate += CalK(pcv[i].GetAverageDiameter(), pcv[j].GetAverageDiameter()) * pcv[j].GetSize();
		}

		cumulativeRate += rate * pcv[i].GetSize();

		pcv[i].SetCumulativeRate(cumulativeRate);

		rate = 0.0;
	}
}

Particle ParticleClassVector::UpdateParticleAndRate(unsigned int i, unsigned int j, double rate1, double rate2,
		double criterion) {
	unsigned int indexI = (criterion - rate1) / (rate2 - rate1) * pcv[i].GetSize();
	unsigned int indexJ = (criterion - rate1) / (rate2 - rate1) * pcv[j].GetSize();
	unsigned int newParticleIndex = 0;

	double newGlycerinMol = pcv[i].GetParticleAt(indexI).GetGlycerinMol()
			+ pcv[j].GetParticleAt(indexJ).GetGlycerinMol();
	double newWasserMol = pcv[i].GetParticleAt(indexI).GetWasserMol() + pcv[j].GetParticleAt(indexJ).GetWasserMol();
	double newZusammensetzung = newGlycerinMol / (newGlycerinMol + newWasserMol);

	unsigned int newNumberOfWasser = pcv[i].GetParticleAt(indexI).GetNumberOfWasser()
			+ pcv[j].GetParticleAt(indexJ).GetNumberOfWasser();
	unsigned int newNumberOfGlycerin = pcv[i].GetParticleAt(indexI).GetNumberOfGlycerin()
			+ pcv[j].GetParticleAt(indexJ).GetNumberOfGlycerin();

	double newParticleDiameter = 6
			* (newNumberOfGlycerin * ps.GlycerinMolekuelmasse + newNumberOfWasser * ps.Wassermolekuelmasse) / M_PI
			/ ((1 - newZusammensetzung) * ps.WasserDicht + newZusammensetzung * ps.GlycerinDicht);
	newParticleDiameter = pow(newParticleDiameter, 0.333333);

	double ithOldAverageDiameter = pcv[i].GetAverageDiameter();
	double jthOldAverageDiameter = pcv[j].GetAverageDiameter();
	double oldAverageDiameter = 0.0;

	bool foundClass = false;

	Particle newParticle;

	double ithDeltaRate = 0.0;
	double jthDeltaRate = 0.0;
	double newDeltaRate = 0.0;

	pcv[i].EraseParticleAt(indexI);
	pcv[j].EraseParticleAt(indexJ);

	for (unsigned int a = 0; a < pcv.size(); a++) {
		if (newParticleDiameter >= pcv[a].GetLeftBoundary() && newParticleDiameter <= pcv[a].GetRightBoundary()) {
			oldAverageDiameter = pcv[a].GetAverageDiameter();
			newParticle = Particle(newParticleDiameter, newZusammensetzung, newWasserMol, newGlycerinMol,
					newNumberOfWasser, newNumberOfGlycerin);
			pcv[a].PushBackParticle(newParticle);
			pcv[a].UpdateAverageDiameter(newParticleDiameter);
			newParticleIndex = a;
			foundClass = true;
			break;
		}
	}

	if (!foundClass) {
		int n = newParticleDiameter / interval;
		InsertNewParticleClass();
		pcv[pcv.size() - 1].SetLeftBoundary(n * interval);
		pcv[pcv.size() - 1].SetRightBoundary((n + 1) * interval);
		pcv[pcv.size() - 1].PushBackParticle(
				Particle(newParticleDiameter, newZusammensetzung, newWasserMol, newGlycerinMol, newNumberOfWasser,
						newNumberOfGlycerin));
		pcv[pcv.size() - 1].CalAverageDiameter();

		Clear();

		InitializeCumulativeRate();

		newParticleIndex = pcv.size() - 1;

		totalNumberOfParticles--;

		return pcv[pcv.size() - 1].GetParticleAt(0);
	}

	if (newParticleIndex == i || newParticleIndex == j || i == j) {
		InitializeCumulativeRate();
	} else {
		double temp1 = 0.0;
		if (i == 0) {
			temp1 = pcv[0].GetCumulativeRate();

		} else {
			temp1 = pcv[i].GetCumulativeRate() - pcv[i - 1].GetCumulativeRate();
		}

		double temp2 = CalK(pcv[i].GetAverageDiameter(), pcv[i].GetAverageDiameter()) * pcv[i].GetSize()
				* (pcv[i].GetSize() - 1) * 0.5;
		for (unsigned int c = i + 1; c < pcv.size(); c++) {
			temp2 += CalK(pcv[i].GetAverageDiameter(), pcv[c].GetAverageDiameter()) * pcv[i].GetSize()
					* pcv[c].GetSize();
		}
		ithDeltaRate = temp2 - temp1;

		double temp4 = 0.0;
		if (j == 0) {
			temp4 = pcv[0].GetCumulativeRate();
		} else {
			temp4 = pcv[j].GetCumulativeRate() - pcv[j - 1].GetCumulativeRate();
		}

		double temp5 = CalK(pcv[j].GetAverageDiameter(), pcv[j].GetAverageDiameter()) * pcv[j].GetSize()
				* (pcv[j].GetSize() - 1) * 0.5;
		for (unsigned int d = j + 1; d < pcv.size(); d++) {
			temp5 += CalK(pcv[j].GetAverageDiameter(), pcv[d].GetAverageDiameter()) * pcv[j].GetSize()
					* pcv[d].GetSize();
		}
		jthDeltaRate = temp5 - temp4;

		double temp7 = 0.0;
		if (newParticleIndex == 0) {
			temp7 = pcv[0].GetCumulativeRate();
		} else {
			temp7 = pcv[newParticleIndex].GetCumulativeRate() - pcv[newParticleIndex - 1].GetCumulativeRate();

		}

		double temp8 = CalK(pcv[newParticleIndex].GetAverageDiameter(), pcv[newParticleIndex].GetAverageDiameter())
				* pcv[newParticleIndex].GetSize() * (pcv[newParticleIndex].GetSize() - 1) * 0.5;
		for (unsigned int e = newParticleIndex + 1; e < pcv.size(); e++) {
			temp8 += CalK(pcv[newParticleIndex].GetAverageDiameter(), pcv[e].GetAverageDiameter())
					* pcv[newParticleIndex].GetSize() * pcv[e].GetSize();
		}
		newDeltaRate = temp8 - temp7;

		double deltaRate = 0.0;
		for (unsigned int f = 0; f < pcv.size(); f++) {
			if (f == i) {
				deltaRate += ithDeltaRate;
			}

			if (f == j) {
				deltaRate += jthDeltaRate;
			}

			if (f == newParticleIndex) {
				deltaRate += newDeltaRate;
			}

			if (f < i) {
				deltaRate -= (CalK(pcv[f].GetAverageDiameter(), ithOldAverageDiameter) * pcv[f].GetSize()
						* (pcv[i].GetSize() + 1)
						- CalK(pcv[f].GetAverageDiameter(), pcv[i].GetAverageDiameter()) * pcv[f].GetSize()
								* pcv[i].GetSize());
			}

			if (f < j) {
				deltaRate -= (CalK(pcv[f].GetAverageDiameter(), jthOldAverageDiameter) * pcv[f].GetSize()
						* (pcv[j].GetSize() + 1)
						- CalK(pcv[f].GetAverageDiameter(), pcv[j].GetAverageDiameter()) * pcv[f].GetSize()
								* pcv[j].GetSize());
			}

			if (f < newParticleIndex) {
				deltaRate += CalK(pcv[f].GetAverageDiameter(), pcv[newParticleIndex].GetAverageDiameter())
						* pcv[f].GetSize() * pcv[newParticleIndex].GetSize()
						- CalK(pcv[f].GetAverageDiameter(), oldAverageDiameter) * pcv[f].GetSize()
								* (pcv[newParticleIndex].GetSize() - 1);
			}

			pcv[f].UpdateCumulativeRate(deltaRate);
		}
	}

	Clear();

	totalNumberOfParticles--;

	return newParticle;
}

void ParticleClassVector::ReduceParticles() {
	vector<Particle> tempPV;

	for (unsigned int i = 0; i < pcv.size(); i++) {
		for (unsigned int j = 0; j < pcv[i].GetSize(); j++) {
			tempPV.push_back(pcv[i].GetParticleAt(j));
		}
	}

	sort(tempPV.begin(), tempPV.end(), CompareByDiameter);

	pcv.clear();
	int n = tempPV[0].GetDiameter() / interval;
	InsertNewParticleClass();
	pcv[0].SetLeftBoundary(n * interval);
	pcv[0].SetRightBoundary((n + 1) * interval);
	pcv[0].PushBackParticle(tempPV[0]);
	pcv[0].CalAverageDiameter();

	unsigned int currentClassIndex = 0;
	for (unsigned int i = 2; i < tempPV.size(); i = i + 2) {
		if (tempPV[i].GetDiameter() >= pcv[currentClassIndex].GetLeftBoundary()
				&& tempPV[i].GetDiameter() <= pcv[currentClassIndex].GetRightBoundary()) {
			pcv[currentClassIndex].PushBackParticle(tempPV[i]);
			pcv[currentClassIndex].UpdateAverageDiameter(tempPV[i].GetDiameter());
		} else {
			int m = tempPV[i].GetDiameter() / interval;
			InsertNewParticleClass();
			currentClassIndex++;
			pcv[currentClassIndex].SetLeftBoundary(m * interval);
			pcv[currentClassIndex].SetRightBoundary((m + 1) * interval);
			pcv[currentClassIndex].PushBackParticle(tempPV[i]);
			pcv[currentClassIndex].CalAverageDiameter();
		}
	}

	InitializeCumulativeRate();

	totalNumberOfParticles = GetNumberOfAllParticles();
}

void ParticleClassVector::DuplicateParticles() {
	vector<Particle> tempPV;

	for (unsigned int i = 0; i < pcv.size(); i++) {
		for (unsigned int j = 0; j < pcv[i].GetSize(); j++) {
			tempPV.push_back(pcv[i].GetParticleAt(j));
		}
	}

	sort(tempPV.begin(), tempPV.end(), CompareByDiameter);

	pcv.clear();
	int n = tempPV[0].GetDiameter() / interval;
	InsertNewParticleClass();
	pcv[0].SetLeftBoundary(n * interval);
	pcv[0].SetRightBoundary((n + 1) * interval);
	pcv[0].PushBackParticle(tempPV[0]);
	pcv[0].PushBackParticle(tempPV[0]);
	pcv[0].CalAverageDiameter();

	unsigned int currentClassIndex = 0;
	for (unsigned int i = 1; i < tempPV.size(); i++) {
		if (tempPV[i].GetDiameter() >= pcv[currentClassIndex].GetLeftBoundary()
				&& tempPV[i].GetDiameter() <= pcv[currentClassIndex].GetRightBoundary()) {
			pcv[currentClassIndex].PushBackParticle(tempPV[i]);
			pcv[currentClassIndex].UpdateAverageDiameter(tempPV[i].GetDiameter());
			pcv[currentClassIndex].PushBackParticle(tempPV[i]);
			pcv[currentClassIndex].UpdateAverageDiameter(tempPV[i].GetDiameter());
		} else {
			int m = tempPV[i].GetDiameter() / interval;
			InsertNewParticleClass();
			currentClassIndex++;
			pcv[currentClassIndex].SetLeftBoundary(m * interval);
			pcv[currentClassIndex].SetRightBoundary((m + 1) * interval);
			pcv[currentClassIndex].PushBackParticle(tempPV[i]);
			pcv[currentClassIndex].PushBackParticle(tempPV[i]);
			pcv[currentClassIndex].CalAverageDiameter();
		}
	}

	InitializeCumulativeRate();

	totalNumberOfParticles = GetNumberOfAllParticles();
}

vector<Particle> ParticleClassVector::GetAllParticles() {
	vector<Particle> tempPV;

	for (unsigned int i = 0; i < pcv.size(); i++) {
		for (unsigned int j = 0; j < pcv[i].GetSize(); j++) {
			tempPV.push_back(Particle(pcv[i].GetParticleAt(j)));
		}
	}

	sort(tempPV.begin(), tempPV.end(), CompareByDiameter);

	return tempPV;
}

void ParticleClassVector::Clear() {
	int firstEmptyClassIndex = -1;
	int secondEmptyClassIndex = -1;
	bool foundFirstEmptyClass = false;

	for (unsigned int i = 0; i < pcv.size(); i++) {
		if (pcv[i].GetSize() == 0 && foundFirstEmptyClass == false) {
			firstEmptyClassIndex = i;
			foundFirstEmptyClass = true;
		} else if (pcv[i].GetSize() == 0 && foundFirstEmptyClass == true) {
			secondEmptyClassIndex = i;
		}
	}

	if (firstEmptyClassIndex >= 0 && secondEmptyClassIndex >= 0) {
		it = pcv.begin();
		pcv.erase(it + firstEmptyClassIndex);
		pcv.erase(it + secondEmptyClassIndex - 1);
	} else if (firstEmptyClassIndex >= 0) {
		it = pcv.begin();
		pcv.erase(it + firstEmptyClassIndex);
	}
}

double ParticleClassVector::GetMittelDurchmesser() {
	double durchmesser = 0.0;
	unsigned int totalNum = 0;

	for (unsigned int i = 0; i < pcv.size(); i++) {
			for (unsigned int j = 0; j < pcv[i].GetSize(); j++) {
				durchmesser += pcv[i].GetParticleAt(j).GetDiameter();
				totalNum ++;
			}
		}

	if (totalNum > 0) {
		return durchmesser / totalNum;
	} else {
		return 0.0;
	}
}

double ParticleClassVector::GetMittelZusammensetzung() {
	double zusammensetzung = 0.0;
	unsigned int totalNum = 0;

	for (unsigned int i = 0; i < pcv.size(); i++) {
			for (unsigned int j = 0; j < pcv[i].GetSize(); j++) {
				zusammensetzung += pcv[i].GetParticleAt(j).GetZusammensetzung();
				totalNum ++;
			}
		}

	if (totalNum > 0) {
		return zusammensetzung / totalNum;
	} else {
		return 0.0;
	}
}

double ParticleClassVector::CalK(double firstDiameter, double secondDiameter) {
	if (firstDiameter == 0 || secondDiameter == 0) {
		return 0;
	} else {
		double ro = 1.13E3;
		double kb = 1.3807E-23;
		double T = 298.15;
		double visko = 7E-6;
		double m1 = M_PI / 6 * ro * pow(firstDiameter, 3);
		double m2 = M_PI / 6 * ro * pow(secondDiameter, 3);
		double D1 = kb * T / (3 * visko * M_PI * firstDiameter);
		double D2 = kb * T / (3 * visko * M_PI * secondDiameter);
		double cquer1 = sqrt(8 * kb * T / M_PI / m1);
		double cquer2 = sqrt(8 * kb * T / M_PI / m2);
		double l1 = 8 * D1 / M_PI / cquer1;
		double l2 = 8 * D2 / M_PI / cquer2;
		double g1 = (pow(firstDiameter + l1, 3) - pow(pow(firstDiameter, 2) + pow(l1, 2), 3 / 2)) / 3 / firstDiameter
				/ l1 - firstDiameter;
		double g2 = (pow(secondDiameter + l2, 3) - pow(pow(secondDiameter, 2) + pow(l2, 2), 3 / 2)) / 3 / secondDiameter
				/ l2 - secondDiameter;
		double g12 = sqrt(pow(g1, 2) + pow(g2, 2));
		double c12 = sqrt(pow(cquer1, 2) + pow(cquer2, 2));
		double F = 1 - 2 * g12 / (firstDiameter + secondDiameter + 2 * g12)
				+ 8 * (D1 + D2) / c12 / (firstDiameter + secondDiameter);
		double k = 2 * M_PI * (D1 + D2) * (firstDiameter + secondDiameter) / F;

		return k;
	}
}

void ParticleClassVector::UpdateBoundary() {
	pcv[0].SortMe();

	double x1 = pcv[0].GetMinDiameter();
	double x2 = pcv[0].GetMaxDiameter();
	double x3 = x1;
	double x4 = x2;

	for (unsigned int i = 1; i < pcv.size(); i++) {
		pcv[i].SortMe();

		x3 = pcv[i].GetMinDiameter();
		x4 = pcv[i].GetMaxDiameter();
		pcv[i - 1].SetLeftBoundary(x1);
		pcv[i - 1].SetRightBoundary((x2 + x3) / 2);
		x1 = (x2 + x3) / 2;
		x2 = x4;
	}

	pcv[pcv.size() - 1].SetLeftBoundary(x1);
	pcv[pcv.size() - 1].SetRightBoundary(x2);

	rightBoundary = x2;
	lastClassRange = x2 - x1;
}

