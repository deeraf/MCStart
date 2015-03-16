/*
 * Particle.cpp
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */

#include "Particle.h"

#include <iostream>

using namespace std;

Particle::Particle() {
}

Particle::Particle(double d) {
	diameter = d;
}

Particle::Particle(double d, double z, double wm, double gm, unsigned int wn, unsigned int gn) {
	diameter = d;
	zusammensetzung = z;
	wasserMol = wm;
	glycerinMol = gm;
	numberOfWasser = wn;
	numberOfGlycerin = gn;
}

double Particle::GetDiameter() {
	return diameter;
}

void Particle::SetDiameter(double d) {
	diameter = d;
}

double Particle::GetZusammensetzung() {
	return zusammensetzung;
}

void Particle::SetZusammensetzung(double z) {
	zusammensetzung = z;
}

double Particle::GetWasserMol() {
	return wasserMol;
}

void Particle::SetWasserMol(double wm) {
	wasserMol = wm;
}

double Particle::GetGlycerinMol() {
	return glycerinMol;
}

void Particle::SetGlycerinMol(double gm) {
	glycerinMol = gm;
}

unsigned int Particle::GetNumberOfGlycerin() {
	return numberOfGlycerin;
}

void Particle::SetNumberOfGlycerin(unsigned int num) {
	numberOfGlycerin = num;
}

unsigned int Particle::GetNumberOfWasser() {
	return numberOfWasser;
}

void Particle::SetNumberOfWasser(unsigned int num) {
	numberOfWasser = num;
}

