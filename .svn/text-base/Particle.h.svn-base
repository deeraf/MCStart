/*
 * Particle.h
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

class Particle {
public:
	Particle();
	Particle(double diameter);
	Particle(double d, double z, double wm,double gm, unsigned int wn, unsigned int gn);
	double GetDiameter();
	double GetWasserMol();
	double GetGlycerinMol();
	double GetZusammensetzung();
	unsigned int GetNumberOfWasser();
	unsigned int GetNumberOfGlycerin();
	void SetDiameter(double d);
	void SetWasserMol(double wm);
	void SetGlycerinMol(double gm);
	void SetZusammensetzung(double z);
	void SetNumberOfWasser(unsigned int num);
	void SetNumberOfGlycerin(unsigned int num);
	double diameter;
private:
	double zusammensetzung;
	double wasserMol;
	double glycerinMol;
	unsigned int numberOfWasser;
	unsigned int numberOfGlycerin;
};

#endif /* PARTICLE_H_ */
