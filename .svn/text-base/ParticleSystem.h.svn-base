/*
 * ParticleSystem.h
 *
 *  Created on: Aug 13, 2013
 *      Author: amifan
 */

#ifndef PARTICLESYSTEM_H_
#define PARTICLESYSTEM_H_

#include "Particle.h"

using namespace std;

class ParticleSystem {
public:
	ParticleSystem();
	ParticleSystem(double wasserDampMasse, double glycerinDampMasse, double stickstoffDampMasse);
	virtual ~ParticleSystem();
	static const double WasserMolMasse = 1.8E-2;
	static const double GlycerinMolMasse = 9.21E-2;
	static const double WasserDicht = 1E3;
	static const double GlycerinDicht = 1.26E3;
	static const double WasserOberflaeche = 7.3E-2;
	static const double GlycerinOberflaeche = 6.4E-2;
	static const double WasserMolekuelVol = 2.98902E-29;
	static const double GlycerinMolekuelVol = 1.21379E-28;
	static const double Wassermolekuelmasse = 2.98902E-26;
	static const double GlycerinMolekuelmasse = 1.52938E-25;
	static const double WasserUebersaettigungsdruck = 3161.5;
	static const double GlycerinUebersaettigungsdruck = 0.1;
	static const double UmgebungsTemperatur = 298.15;
	static const double GesamtDruck = 1E5;
	static const double StickStoffMolMasse = 2.9E-2;
	static const double BoltzmannConstant = 1.38066E-23;
	static const double AnzahlproMol = 6.022045E+23;
	static const double Gasconstant = 8.3144;
	static const double miniDeltaMasse = 0.08;
	static const double maxDeltaMasse = 0.5;
	static const double MaxDeltaUebersaettigungsrate = 0.05;

	double wasserPartialDruck;
	double glycerinPartialDruck;
	double wasserUebersaettigungsgrad;
	double glycerinUebersaettigungsgrad;
	double keimbildungsRate;

	double wasserDampfMasse;
	double glycerinDampfMasse;
	double stickStoffDampfMasse;
	double kritischGlycerinZusammensetzung;

// Folgende Funktionen zur Berechnung von CalKritischGlycerinZusammensetzung und CalKritischGlycerinZusammensetzungTemp
double func(double PaGly,double PsGly,double PaH2O,double PsH2O,double v1,double v2,double x,double *df);
double Krit_Gly_Newt(double PaGly,double PsGly,double PaH2O,double PsH2O,double v1,double v2,
				double xl,double xr,int jmax,double eps,int *fehler);
//Die folgende Gleichungen gehören zur Keimbildung;
	double CalWasserPartielDruck();
	double CalGlycerinPartielDruck();
	double CalGlycerinUebersaettigungsgrad ();
	double CalWasserUebersaettigungsgrad ();
	double CalKritischGlycerinZusammensetzung();
	double CalKritischGlycerinZusammensetzung3();
	double CalKritischGlycerinZusammensetzungTemp(double tempWasserDampMasse, double tempGlycerinDampMasse);
	double CalKeimGroesse();
	double CalkritischeMolWasser();
	double CalkritischeMolWasser1(Particle p);
	unsigned int CalWasserNumber(double diameter);
	double CalkritischeMolGlyzerin();
	double CalkritischeMolGlyzerin1(Particle p);
	unsigned CalGlycerinNumber(double diameter);
	double CalKeimbildungsrate ();
	void UpdateKritischGlycerinZusammensetzung();
	double GetkritischGlycerinZusammensetzung();
	Particle GetNewBuiltParticle();

// Aktualisierung nach Keimbildung;
	double CalKeimnewDampfmasseGlycerin();
	double CalKeimnewDampfmasseWasser();
	double CalkeimnewpartiellDruckGlycerin();
	double CalkeimnewpartiellDruckWasser();

// Die folgenden Gleichungen gehören zum Wachstum (Nach Keimbildung);
	double CalzunehmenMolekuelGlycerin1(Particle p);
	double CalzunehmenMolekuelWasser1(Particle p);
	double CalWachstumneuMolGlycerin1(Particle p, double t);
	double CalWachstumneuMolWasser1(Particle p,double t);
	double CalWachstumneuZusammensetzungGlycerin1(Particle p, double t);
	double CalWachstumneuGrosse1(Particle p,double t);

// Berechnung für Wachstumsüberprüfung (Nach Keimbildung);
	double CalWachstumneumasseglycerin ();
	double CalWachstumneumassewasser();
	double Calmassenaenderungsrate(Particle p, double deltaT);
	double CalWachstumneuUebersaetigungsgradGlycerin();
	double CalWachstumneuUebersaetigungsgradWasser();
	double CalUebersaettigungsgradaenderungsrate();
	void CalWachstumNewParticle(Particle &p,double deltaT);
	void CalWachstumOnlyNewParticle(Particle &p,double deltaT);

	bool CheckNewBuiltParticle();

	Particle tempParticle;
};

#endif /* PARTICLESYSTEM_H_ */
