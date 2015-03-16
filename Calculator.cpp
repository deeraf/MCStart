/*
 * Calculator.cpp
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */
#include "math.h"
#include "Logger.h"
#include "Particle.h"
#include "Calculator.h"
#include "ParticleSystem.h"
#include "LogNormalDistribution.h"

using namespace std;

Calculator::Calculator(double simulationTime, double classTolerant, unsigned int maxParticleNumber, double interval, double wasserDampMasse,
		double glycerinDampMasse, double stickstoffDampMasse) {
	srand(time(NULL));
	ps = ParticleSystem(wasserDampMasse, glycerinDampMasse, stickstoffDampMasse);
	pcv = ParticleClassVector(interval, classTolerant);
	totalSimulationTime = simulationTime;
	tolerant = classTolerant;
	maxNumberOfParticles = maxParticleNumber;
}

Calculator::Calculator(LogNormalDistribution logNormalDistribution, double simulationTime, double classTolerant, unsigned int maxParticleNumber, double interval, double wasserDampMasse,
		double glycerinDampMasse, double stickstoffDampMasse) {
	srand(time(NULL));
	ps = ParticleSystem(wasserDampMasse, glycerinDampMasse, stickstoffDampMasse);
	pcv = logNormalDistribution.Initialize();;
	totalSimulationTime = simulationTime;
	tolerant = classTolerant;
	maxNumberOfParticles = maxParticleNumber;
}

Calculator::~Calculator() {

}

double Calculator::qRandom() {
	double u1, u2, s;
	double r = 1.0;

	do {
		do {
			//generate 2 uniform random variables on [-1, 1]
			u1 = rand() / (RAND_MAX / 2.0) - 1.0;
			u2 = rand() / (RAND_MAX / 2.0) - 1.0;

			s = u1 * u1 + u2 * u2;
		} while (s == 0 || s >= 1);

		//generate 2 normal distributed variables with mean 0 and standard deviation 1, using Box-Muller in polar form
		double n1 = u1 * sqrt((-2.0 * log(s)) / s);
		//double n2 = u2 * sqrt((-2.0 * log(s)) / s);

		//generate 2 lognormal distributed variables with mean 0 and standard deviation 1
		r = exp(0 + 2 * n1);
		//	double ln2 = exp(0+ 1*n2);
	} while (r > 1);

	return r;
}

void Calculator::DoSimulation() {
	double nextRandom = 0.0;
	//double nextqRandom = 0.0;
	double nextCriterion = 0.0;
	double currentTime = 0.0;
	double currentDeltaT = 0.0;
	double coagulationRate = 0.0;
	double keimbildungRate = 0.0;
	double rate1 = 0.0;
	double rate2 = 0.0;
	double iniVolumen = 0.0;
	// Time to draw the next picture 5E-7
	double currentPaintTime = 5E-7;
	double currentFinePaintTime = 1E-4;
	double currentFinePaintInterval = 1E-4;
	unsigned int finePaintCount = 1;
	bool malTwo = true;
	Particle newParticle;

	unsigned int numberOfCoagulation = 0;
	unsigned int numberOfKeimbildung = 0;
	unsigned int numberOfAction = 0;
	unsigned int currentPictureIndex = 0;

	bool foundParticlePair = false;
	bool keimbildungGenug = false;

	volumen = CalVolumen(maxNumberOfParticles);
	iniVolumen = volumen;

	while (currentTime <= totalSimulationTime) {
		if (currentTime >= 1E-4) {
			if (currentTime >= currentFinePaintTime) {
				//logger.DrawPlot(pcv.GetXResultAsIntervalNew(), pcv.GetQxResultAsIntervalNew(), pcv.GetNumberOfParticleClasses(), currentPictureIndex, currentTime);
				logger.WriteToFileParticle(pcv, volumen, currentPictureIndex);
				logger.WriteToFileParticleCount(pcv, volumen, currentPictureIndex);
				logger.WriteToFile3D(pcv, currentPictureIndex++);

				if (finePaintCount <= 9) {
					currentFinePaintTime += currentFinePaintInterval;
					finePaintCount ++;
				}else {
					currentFinePaintInterval *= 10;
					finePaintCount = 1;
					currentFinePaintTime += currentFinePaintInterval;
				}
			}
		}else {
			if (currentTime >= currentPaintTime) {
				//logger.DrawPlot(pcv.GetXResultAsIntervalNew(), pcv.GetQxResultAsIntervalNew(), pcv.GetNumberOfParticleClasses(), currentPictureIndex, currentTime);
				logger.WriteToFileParticle(pcv, volumen, currentPictureIndex);
				logger.WriteToFileParticleCount(pcv, volumen, currentPictureIndex);
				logger.WriteToFile3D(pcv, currentPictureIndex++);
				if (malTwo) {
					currentPaintTime *= 2;
					malTwo = false;
				} else {
					currentPaintTime *= 5;
					malTwo = true;
				}
			}
		}

		if (pcv.GetNumberOfParticles() >= maxNumberOfParticles) {
			pcv.ReduceParticles();
			volumen = volumen * 0.5;
		} else if (keimbildungGenug && pcv.GetNumberOfParticles() <= maxNumberOfParticles * 0.25) {
			pcv.DuplicateParticles();
			volumen = volumen * 2;
		}

		if (numberOfAction >= 500) {
			logger.WriteToFileMittelWert(pcv.GetMittelZusammensetzung(), pcv.GetMittelDurchmesser(), pcv.GetNumberOfAllParticles(), volumen, currentTime);
			numberOfAction = 0;
		}

		nextRandom = (double) rand() / (double) RAND_MAX;
		//nextqRandom = qRandom();

		ps.UpdateKritischGlycerinZusammensetzung();

		keimbildungRate = ps.CalKeimbildungsrate();
		coagulationRate = pcv.GetTotalRate();

		if (keimbildungRate / (keimbildungRate + coagulationRate / pow(volumen,2)) < nextRandom || ps.wasserDampfMasse <= 0
				|| ps.glycerinDampfMasse <= 0) {
			nextCriterion = coagulationRate * nextRandom;
			foundParticlePair = false;

			for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
				if (foundParticlePair) {
					break;
				} else if (pcv.GetParticleClassAt(i).GetCumulativeRate() >= nextCriterion) {
					if (i > 0) {
						nextCriterion -= pcv.GetParticleClassAt(i - 1).GetCumulativeRate();
					}

					rate1 = rate2 = 0.0;
					for (unsigned int j = i; j < pcv.GetNumberOfParticleClasses(); j++) {
						rate1 = rate2;

						if (i == j) {
							rate2 += pcv.CalK(pcv.GetParticleClassAt(i).GetAverageDiameter(), pcv.GetParticleClassAt(i).GetAverageDiameter())
									* pcv.GetParticleClassAt(i).GetSize() * (pcv.GetParticleClassAt(i).GetSize() - 1) * 0.5;
						} else {
							rate2 += pcv.CalK(pcv.GetParticleClassAt(i).GetAverageDiameter(), pcv.GetParticleClassAt(j).GetAverageDiameter())
									* pcv.GetParticleClassAt(i).GetSize() * pcv.GetParticleClassAt(j).GetSize();
						}

						if (rate1 < nextCriterion && rate2 >= nextCriterion) {
							if (i == j && pcv.GetParticleClassAt(i).GetSize() < 2) {
								foundParticlePair = true;
								break;
							} else {
								double id = pcv.GetParticleClassAt(i).GetAverageDiameter();
								double jd = pcv.GetParticleClassAt(j).GetAverageDiameter();
								newParticle = pcv.UpdateParticleAndRate(i, j, rate1, rate2, nextCriterion);
								foundParticlePair = true;
								numberOfCoagulation++;
								numberOfAction++;
								currentDeltaT = GetDeltaTime(keimbildungRate, coagulationRate);
								currentTime += currentDeltaT;
								logger.PrintTraceCoagulation(pcv.GetNumberOfAllParticles(), numberOfCoagulation, id, jd, newParticle, nextCriterion, nextRandom,
								coagulationRate, currentTime);
								logger.WriteToFileCoagulation(pcv.GetNumberOfAllParticles(), numberOfCoagulation,
								 id, jd, newParticle, nextCriterion,nextRandom, coagulationRate,currentTime);
								break;
							}
						}
					}
				}
			}

			if (foundParticlePair == false) {
				pcv.InitializeCumulativeRate();
			}
		} else {
			if (ps.CheckNewBuiltParticle()) {
				Particle p = ps.GetNewBuiltParticle();
				pcv.AddParticleAndUpdateRate(p);
				double deltaKeimbildungWasserDampMasse = p.GetWasserMol() * ps.WasserMolMasse / volumen;
				double deltaKeimbildungGlycerinDampMasse = p.GetGlycerinMol() * ps.GlycerinMolMasse / volumen;

				if (ps.wasserDampfMasse - deltaKeimbildungWasserDampMasse < 0) {
					ps.wasserDampfMasse = 0;
				} else {
					ps.wasserDampfMasse -= deltaKeimbildungWasserDampMasse;
				}

				if (ps.glycerinDampfMasse - deltaKeimbildungGlycerinDampMasse < 0) {
					ps.glycerinDampfMasse = 0;
				} else {
					ps.glycerinDampfMasse -= deltaKeimbildungGlycerinDampMasse;
				}

				currentDeltaT = GetDeltaTime(keimbildungRate, coagulationRate);
				currentTime += currentDeltaT;
				numberOfKeimbildung++;
				numberOfAction++;

				if (numberOfKeimbildung >= maxNumberOfParticles * 0.25) {
					keimbildungGenug = true;
				}

				logger.PrintTraceKeimbildung(pcv.GetNumberOfAllParticles(), numberOfKeimbildung, p.GetDiameter(), p.GetZusammensetzung(),
				keimbildungRate, currentTime);
				logger.WriteToFileKeimbildung(pcv.GetNumberOfAllParticles(), numberOfKeimbildung, p.GetDiameter(), p.GetZusammensetzung(), keimbildungRate, currentTime);
				logger.WriteToFileKeimbildungMasse(deltaKeimbildungWasserDampMasse, deltaKeimbildungGlycerinDampMasse, ps.wasserDampfMasse, ps.glycerinDampfMasse, currentTime);
			}
		}

		if (pcv.GetNumberOfParticleClasses() > 0 && ps.wasserDampfMasse > 0 && ps.glycerinDampfMasse > 0) {
			double c = pow(ps.CalGlycerinPartielDruck() / ps.GlycerinUebersaettigungsdruck, ps.kritischGlycerinZusammensetzung)
					* pow(ps.CalWasserPartielDruck() / ps.WasserUebersaettigungsdruck, (1 - ps.kritischGlycerinZusammensetzung));
			double tempwasserDampfMasse = ps.wasserDampfMasse;
			double tempglycerinDampfMasse = ps.glycerinDampfMasse;

			for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
				double x = ps.CalzunehmenMolekuelGlycerin1(pcv.GetParticleClassAt(i).GetParticleAt(0));
				double y = ps.CalzunehmenMolekuelWasser1(pcv.GetParticleClassAt(i).GetParticleAt(0));
				tempglycerinDampfMasse -= x * ps.GlycerinMolekuelmasse * pcv.GetParticleClassAt(i).GetSize();
				tempwasserDampfMasse -= y * ps.Wassermolekuelmasse * pcv.GetParticleClassAt(i).GetSize();
			}

			tempwasserDampfMasse *= (currentDeltaT / volumen);
			tempglycerinDampfMasse *= (currentDeltaT / volumen);

			double tempGlycerinPartielDruck =
					(tempglycerinDampfMasse / ps.GlycerinMolMasse)
							/ (tempglycerinDampfMasse / ps.GlycerinMolMasse + tempwasserDampfMasse / ps.WasserMolMasse
									+ ps.stickStoffDampfMasse / ps.StickStoffMolMasse) * ps.GesamtDruck;
			double tempWasserPartielDruck =
					(tempwasserDampfMasse / ps.WasserMolMasse)
							/ (tempglycerinDampfMasse / ps.GlycerinMolMasse + tempwasserDampfMasse / ps.WasserMolMasse
									+ ps.stickStoffDampfMasse / ps.StickStoffMolMasse) * ps.GesamtDruck;
			double d = pow(tempGlycerinPartielDruck / ps.GlycerinUebersaettigungsdruck,
					ps.CalKritischGlycerinZusammensetzungTemp(tempwasserDampfMasse, tempglycerinDampfMasse))
					* pow(tempWasserPartielDruck / ps.WasserUebersaettigungsdruck,
							(1 - ps.CalKritischGlycerinZusammensetzungTemp(tempwasserDampfMasse, tempglycerinDampfMasse)));

			if (((c - d) / c) < ps.MaxDeltaUebersaettigungsrate) {
				Particle minParticle = pcv.GetParticleClassAt(0).GetParticleAt(0);
				double massenAenderungsrate = ps.Calmassenaenderungsrate(minParticle, currentDeltaT);

				if (massenAenderungsrate < ps.miniDeltaMasse) {
					// Do nothing
				} else if (massenAenderungsrate > ps.maxDeltaMasse) {
					DoWachstum(currentDeltaT * 0.5, currentTime, volumen);
					DoWachstum(currentDeltaT * 0.5, currentTime, volumen);
					pcv.InitializeCumulativeRate();
				} else {
					DoWachstum(currentDeltaT, currentTime, volumen);
					pcv.InitializeCumulativeRate();

				}
			} else {
				DoWachstum(currentDeltaT * 0.5, currentTime, volumen);
				DoWachstum(currentDeltaT * 0.5, currentTime, volumen);
				pcv.InitializeCumulativeRate();

			}
		}
	}

	//logger.DrawPlot(pcv.GetXResultAsIntervalNew(), pcv.GetQxResultAsIntervalNew(), pcv.GetNumberOfParticleClasses(), currentPictureIndex, currentTime);
	logger.WriteToFileParticleCount(pcv, volumen, currentPictureIndex);
	logger.WriteToFileParticle(pcv, volumen, currentPictureIndex++);
	logger.WriteToFile3D(pcv, currentPictureIndex++);

}

void Calculator::DoOnlyWachstum() {
	double currentTime = 0.0;
	double currentDeltaT = 0.0;
	double coagulationRate = 0.0;
	volumen = 1E-17;//CalVolumen(maxNumberOfParticles);
	pcv.InitializeCumulativeRate();

	//logger.DrawPlot(pcv.GetXResultAsIntervalNew(), pcv.GetQxResultAsIntervalNew(), pcv.GetNumberOfParticleClasses(), 999, currentTime);
	logger.WriteToFileParticle(pcv, volumen, 999);


	while (currentTime <= totalSimulationTime) {
		coagulationRate = pcv.GetTotalRate();
		currentDeltaT = 1E-09;//GetDeltaTime(0, coagulationRate);
		currentTime += currentDeltaT;
		DoOnlyWachstumForOnlyWachtum(currentDeltaT * 0.5, currentTime, volumen);
		cout<<"!!!!!!!"<<currentTime<<"!!!!!!!"<<endl;
	}
}

void Calculator::DoOnlyWachstumForOnlyWachtum(double dt, double ct, double vol) {
	double deltaT = dt;
	double deltaGlycerinMol = 0.0;
	double deltaWasserMol = 0.0;

	if (ps.wasserDampfMasse <= 0 || ps.glycerinDampfMasse <= 0) {

	} else {
		for (unsigned int i = 1; i < pcv.GetNumberOfParticleClasses(); i++) {
			for (unsigned int j = 0; j < pcv.GetParticleClassAt(i).GetSize(); j++) {
				deltaGlycerinMol -= pcv.GetParticleClassAt(i).GetParticleAt(j).GetGlycerinMol();
				deltaWasserMol -= pcv.GetParticleClassAt(i).GetParticleAt(j).GetWasserMol();

				ps.CalWachstumOnlyNewParticle(pcv.GetParticleClassAt(i).GetParticleAt(j), deltaT);

				deltaGlycerinMol += pcv.GetParticleClassAt(i).GetParticleAt(j).GetGlycerinMol();
				deltaWasserMol += pcv.GetParticleClassAt(i).GetParticleAt(j).GetWasserMol();
			}
		}
		double deltaWachstumWasserDampMasse = deltaWasserMol * ps.WasserMolMasse / vol;
		double deltaWachstumGlycerinDampMasse = deltaGlycerinMol * ps.GlycerinMolMasse / vol;

		if (ps.wasserDampfMasse - deltaWachstumWasserDampMasse < 0) {
			ps.wasserDampfMasse = 0;
		} else {
			ps.wasserDampfMasse -= deltaWachstumWasserDampMasse;
		}

		if (ps.glycerinDampfMasse - deltaWachstumGlycerinDampMasse < 0) {
			ps.glycerinDampfMasse = 0;
		} else {
			ps.glycerinDampfMasse -= deltaWachstumGlycerinDampMasse;
		}

		if (deltaWachstumGlycerinDampMasse > 0 && deltaWachstumWasserDampMasse > 0) {
			logger.WriteToFileWachstum(deltaWachstumWasserDampMasse, deltaWachstumGlycerinDampMasse, ps.wasserDampfMasse, ps.glycerinDampfMasse, ct);
		}

		cout<<deltaWasserMol<<"  "<<deltaGlycerinMol<<endl;


	}
}

void Calculator::DoWachstum(double dt, double ct, double vol) {
	double deltaT = dt;
	double deltaGlycerinMol = 0.0;
	double deltaWasserMol = 0.0;

	if (ps.wasserDampfMasse <= 0 || ps.glycerinDampfMasse <= 0) {

	} else {
		for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
			for (unsigned int j = 0; j < pcv.GetParticleClassAt(i).GetSize(); j++) {
				deltaGlycerinMol -= pcv.GetParticleClassAt(i).GetParticleAt(j).GetGlycerinMol();
				deltaWasserMol -= pcv.GetParticleClassAt(i).GetParticleAt(j).GetWasserMol();

			//	cout<<"1 "<< pcv.GetParticleClassAt(i).GetParticleAt(j).GetGlycerinMol()<<"  "<<pcv.GetParticleClassAt(i).GetParticleAt(j).GetWasserMol()<<endl;

				ps.CalWachstumNewParticle(pcv.GetParticleClassAt(i).GetParticleAt(j), deltaT);

				deltaGlycerinMol += pcv.GetParticleClassAt(i).GetParticleAt(j).GetGlycerinMol();
				deltaWasserMol += pcv.GetParticleClassAt(i).GetParticleAt(j).GetWasserMol();

				//cout<<"2 "<< pcv.GetParticleClassAt(i).GetParticleAt(j).GetGlycerinMol()<<"  "<<pcv.GetParticleClassAt(i).GetParticleAt(j).GetWasserMol()<<endl;

			}
		}
		double deltaWachstumWasserDampMasse = deltaWasserMol * ps.WasserMolMasse / vol;
		double deltaWachstumGlycerinDampMasse = deltaGlycerinMol * ps.GlycerinMolMasse / vol;

		if (ps.wasserDampfMasse - deltaWachstumWasserDampMasse < 0) {
			ps.wasserDampfMasse = 0;
		} else {
			ps.wasserDampfMasse -= deltaWachstumWasserDampMasse;
		}

		if (ps.glycerinDampfMasse - deltaWachstumGlycerinDampMasse < 0) {
			ps.glycerinDampfMasse = 0;
		} else {
			ps.glycerinDampfMasse -= deltaWachstumGlycerinDampMasse;
		}

		if (deltaWachstumGlycerinDampMasse > 0 && deltaWachstumWasserDampMasse > 0) {
			logger.WriteToFileWachstum(deltaWachstumWasserDampMasse, deltaWachstumGlycerinDampMasse, ps.wasserDampfMasse, ps.glycerinDampfMasse, ct);
		}

		cout<<deltaWasserMol<<"  "<<deltaGlycerinMol<<endl;


	}
}

double Calculator::GetDeltaTime(double kr, double cr) {
	return 1 / (cr / volumen + kr * volumen);
}

double Calculator::CalVolumen(unsigned int max) {
	ps.UpdateKritischGlycerinZusammensetzung();

	double diameter = ps.CalKeimGroesse();
	double beta = pcv.CalK(diameter, diameter);
	double C = 0.1;
	double Jo = ps.CalKeimbildungsrate();

	return pow(beta * C * max * max / Jo, 0.5);
}
