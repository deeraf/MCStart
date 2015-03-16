/*
 * ParticleSystem.cpp
 *
 *  Created on: Aug 13, 2013
 *      Author: amifan
 */

#include "ParticleSystem.h"
#include "math.h"
#include <iostream>

using namespace std;

ParticleSystem::ParticleSystem() {
	wasserDampfMasse = 0.0;
	glycerinDampfMasse = 0.0;
	stickStoffDampfMasse = 0.0;
	kritischGlycerinZusammensetzung = 0.0;
}

ParticleSystem::ParticleSystem(double wdm, double gdm, double sdm) {
	wasserDampfMasse = wdm;
	glycerinDampfMasse = gdm;
	stickStoffDampfMasse = sdm;
	kritischGlycerinZusammensetzung = CalKritischGlycerinZusammensetzung();
}

ParticleSystem::~ParticleSystem() {

}
// Die folgenden Funktionen realisieren den Berechnung der kritischen Glycerinzusammensetzung
// als Loesung x der Gleichung v1*log(PaGly/(PsGly*x)) - v2*log(PaH2O/(PsH2O*(1-x))) = 0  ---> f(x) = 0
// Die nicht algebraisch berechenbare Loesung geschieht ueber
// kritischGlycerinZusammensetzung = Krit_Gly_Newt(PaGly,PsGly,PaH2O,PsH2O,v1,v2,xl,xr,jmax,eps,&fehler);
// als Newton-Verfahren
// Dabei wird die zu berechnende Funktion inklusive ihrer Ableitung in func(x,&df) realisiert
// df ist der Ableitungswert f'(x)

double ParticleSystem::func(double paGly,double psGly,double paH2O,double psH2O,double v1,double v2,double x,double *df) {
// Diese Funktion berechnet die Funktionswerte sowie
// die Werte der ersten Ableitung
	  *df= v2 / (x - 1) - v1 / x;
	  // Ableitung
	  return v1*(log(paGly)-log(psGly)-log(x)) - v2*(log(paH2O)-log(psH2O)-log(1-x));
	 // Interessant: In der Ableitung kommen nur noch die Parameter v1 und v2 vor
} // zu berechnende Funktion double func(....)

double ParticleSystem::Krit_Gly_Newt(double PaGly,double PsGly,double PaH2O,double PsH2O,double v1,double v2,
		double xl,double xr,int jmax,double eps,int *fehler) {
	// Diese Funktion bestimmt die reelle Nullstelle der Funktion 'func'
	// im Intervall [x1,x2] mit der relativen Genauigkeit 'xacc'.
	// Nach Beendigung der Rechnung kann die Variable 'fehler' einen der
	// folgenden Werte haben:

	//         fehler=0       Die Newton-Raphson Iteration ist ok.
	//         fehler=1       Kein Vorzeichenwechsel der Funktion 'func'
	//                        im Intervall [x1,x2].
	//         fehler=2       Keine Konvergenz waehrend 'jmax' Iter.schritten.
	//         fehler=3       Waehrend der Iteration springt der x Wert
	//                        aus dem Intervall [x1,x2].

	  int j;
	  double x,df,dx;

	  x=0.5*(xl+xr);

	  if((func(PaGly,PsGly,PaH2O,PsH2O,v1,v2,xl,&df)*func(PaGly,PsGly,PaH2O,PsH2O,v1,v2,xr,&df)) > 0.0) {
	    *fehler=1;
	    // ERROR(1): kein Vorzeichenwechsel in [x1,x2]
	    return x;
	  }

	  *fehler=2;
	  j=0;

	  do {
	    dx=func(PaGly,PsGly,PaH2O,PsH2O,v1,v2,x,&df)/df;
	    x=x-dx;

	    if(((xl-x)*(x-xr))<0.0) *fehler=3;
	    else {
	      if(fabs(dx/x) < eps) *fehler=0;
	    }
	    j++;
	  } while ((j<=jmax) && (*fehler==2));

	  // ERROR(2):  keine Konvergenz
	  // ERROR(3):  x aus [x1,x2] gesprungen

	  return x;
} // Newton-Verfahren zur Loesung von f(x) = 0

// Die folgenden Gleichungen gehören zur Keimbildung;

double ParticleSystem::CalGlycerinPartielDruck() {
	return (glycerinDampfMasse / GlycerinMolMasse)
			/ (glycerinDampfMasse / GlycerinMolMasse 
+ wasserDampfMasse / WasserMolMasse 
+ stickStoffDampfMasse / StickStoffMolMasse) * GesamtDruck;
}

double ParticleSystem::CalWasserPartielDruck() {
	return (wasserDampfMasse / WasserMolMasse)
			/ (glycerinDampfMasse / GlycerinMolMasse 
+ wasserDampfMasse / WasserMolMasse 
+ stickStoffDampfMasse / StickStoffMolMasse) * GesamtDruck;
}

double ParticleSystem::CalGlycerinUebersaettigungsgrad() {
	return CalGlycerinPartielDruck() / GlycerinUebersaettigungsdruck / GetkritischGlycerinZusammensetzung();
}

double ParticleSystem::CalWasserUebersaettigungsgrad() {
	return CalWasserPartielDruck() / WasserUebersaettigungsdruck / (1-GetkritischGlycerinZusammensetzung());
}

void ParticleSystem::UpdateKritischGlycerinZusammensetzung() {
	kritischGlycerinZusammensetzung = CalKritischGlycerinZusammensetzung();
}

double ParticleSystem::GetkritischGlycerinZusammensetzung() {
	return kritischGlycerinZusammensetzung;
}
double ParticleSystem::CalKritischGlycerinZusammensetzung() {
	double wasserDruck = 0.0;
	double glycerinDruck = 0.0;
	double wasserMolVolumen = 0.0;
	double glycerinMolVolumen = 0.0;
	double glycerinpartialdruck = CalGlycerinPartielDruck();
//	double glycerinpartialdruck = 0.0;
	double Wasserpartialdruck = CalWasserPartielDruck();
//	double Wasserpartialdruck = 0.0;

	double hvar = 0.0;
	int newtonfehler;

	// wasserDruck = CalWasserUebersaettigungsgrad();
	wasserDruck = Wasserpartialdruck/WasserUebersaettigungsdruck;
	// Wasserpartialdruck = wasserDruck * WasserUebersaettigungsdruck;

	// glycerinDruck = CalGlycerinUebersaettigungsgrad();
	glycerinDruck = glycerinpartialdruck/GlycerinUebersaettigungsdruck;
	// glycerinpartialdruck = glycerinDruck * GlycerinUebersaettigungsdruck;

	wasserMolVolumen = WasserMolMasse / WasserDicht;

	glycerinMolVolumen = GlycerinMolMasse / GlycerinDicht;

	if (wasserDruck <= 1 || glycerinDruck <= 1) {
		return 0;
	} else {
		// alte Ver. hvar = 1 / (log(glycerinDruck) * wasserVolumen / log(wasserDruck) / glycerinVolumen + 1);
		hvar = Krit_Gly_Newt(glycerinpartialdruck,GlycerinUebersaettigungsdruck,
				Wasserpartialdruck,WasserUebersaettigungsdruck,wasserMolVolumen,glycerinMolVolumen,
				0.1,0.9,32,0.000000001,&newtonfehler);
		return hvar;
	}
}

double ParticleSystem::CalKritischGlycerinZusammensetzungTemp(double tempWasserDampMasse, double tempGlycerinDampMasse) {
	double wasserDruck = 0.0;
	double glycerinDruck = 0.0;
	double wasserMolVolumen = 0.0;
	double glycerinMolVolumen = 0.0;
	double hvar = 0.0;
	int newtonfehler = 0;

	double tempWasserPartielDruck = (tempWasserDampMasse / WasserMolMasse)
			/ (tempGlycerinDampMasse / GlycerinMolMasse + tempWasserDampMasse / WasserMolMasse + stickStoffDampfMasse / StickStoffMolMasse) * GesamtDruck;

	double tempGlycerinPartielDruck = (tempGlycerinDampMasse / GlycerinMolMasse)
			/ (tempGlycerinDampMasse / GlycerinMolMasse + tempWasserDampMasse / WasserMolMasse + stickStoffDampfMasse / StickStoffMolMasse) * GesamtDruck;

	wasserDruck = tempWasserPartielDruck / WasserUebersaettigungsdruck;

	glycerinDruck = tempGlycerinPartielDruck / GlycerinUebersaettigungsdruck;

	wasserMolVolumen = WasserMolMasse / WasserDicht;

	glycerinMolVolumen = GlycerinMolMasse / GlycerinDicht;

	if (wasserDruck <= 1 || glycerinDruck <= 1) {
		return 0;
	} else {
		// alte Ver. hvar = 1 / (log(glycerinDruck) * wasserVolumen / log(wasserDruck) / glycerinVolumen + 1);
		hvar = Krit_Gly_Newt(tempGlycerinPartielDruck,GlycerinUebersaettigungsdruck,
				tempWasserPartielDruck,WasserUebersaettigungsdruck,wasserMolVolumen,glycerinMolVolumen,
				0.1,0.9,32,0.000000001,&newtonfehler);
		return hvar;
	}
}

double ParticleSystem::CalKeimGroesse() {
	double mittelOberflaeche = 0.0;
	double mittelVolumen = 0.0;
	double mitteluebersaettigungsgrad = 0.0;
	double wasseruebersaettigungsgrad = 0.0;
	double glyzerinuebersaettigungsgrad = 0.0;

	mittelOberflaeche = (1 - kritischGlycerinZusammensetzung) * WasserOberflaeche + kritischGlycerinZusammensetzung * GlycerinOberflaeche;

	mittelVolumen = (1 - kritischGlycerinZusammensetzung) * WasserMolekuelVol + kritischGlycerinZusammensetzung * GlycerinMolekuelVol;

	wasseruebersaettigungsgrad = CalWasserUebersaettigungsgrad();

	glyzerinuebersaettigungsgrad = CalGlycerinUebersaettigungsgrad();

	mitteluebersaettigungsgrad = pow(wasseruebersaettigungsgrad, (1 - kritischGlycerinZusammensetzung))
			* pow(glyzerinuebersaettigungsgrad, kritischGlycerinZusammensetzung);

	if (mitteluebersaettigungsgrad > 1) {
		return 4 * mittelOberflaeche * mittelVolumen / (BoltzmannConstant * UmgebungsTemperatur * log(mitteluebersaettigungsgrad));
	} else {
		return 0;
	}
}

double ParticleSystem::CalkritischeMolGlyzerin() {
	double mitteldichte = 0.0;
	double durchmesser = CalKeimGroesse();

	mitteldichte = (1 - kritischGlycerinZusammensetzung) * WasserDicht + kritischGlycerinZusammensetzung * GlycerinDicht;

	return M_PI * pow(durchmesser, 3) * mitteldichte * kritischGlycerinZusammensetzung / 6 / AnzahlproMol
			/ ((GlycerinMolekuelmasse - Wassermolekuelmasse) * kritischGlycerinZusammensetzung + Wassermolekuelmasse);
}

double ParticleSystem::CalkritischeMolGlyzerin1(Particle p) {
	double mitteldichte = 0.0;
	double Durchmesser = p.GetDiameter();

	mitteldichte = (1 - p.GetZusammensetzung()) * WasserDicht + p.GetZusammensetzung() * GlycerinDicht;

	return M_PI * pow(Durchmesser, 3) * mitteldichte * p.GetZusammensetzung() / 6 / AnzahlproMol
			/ ((GlycerinMolekuelmasse - Wassermolekuelmasse) * p.GetZusammensetzung() + Wassermolekuelmasse);
}

unsigned int ParticleSystem::CalGlycerinNumber(double durchmesser) {
	unsigned int result = 0;
	double mitteldichte = 0.0;

	mitteldichte = (1 - kritischGlycerinZusammensetzung) * WasserDicht + kritischGlycerinZusammensetzung * GlycerinDicht;

	result = ceil(
			M_PI * pow(durchmesser, 3) * mitteldichte * kritischGlycerinZusammensetzung / 6
					/ ((GlycerinMolekuelmasse - Wassermolekuelmasse) * kritischGlycerinZusammensetzung + Wassermolekuelmasse));

	if (result <= 1) {
		return 1;
	} else {
		return result;
	}
}

double ParticleSystem::CalkritischeMolWasser() {
	double mitteldichte = 0.0;
	double durchmesser = CalKeimGroesse();

	mitteldichte = (1 - kritischGlycerinZusammensetzung) * WasserDicht + kritischGlycerinZusammensetzung * GlycerinDicht;

	return M_PI * pow(durchmesser, 3) * mitteldichte * (1 - kritischGlycerinZusammensetzung) / 6 / AnzahlproMol
			/ ((GlycerinMolekuelmasse - Wassermolekuelmasse) * kritischGlycerinZusammensetzung + Wassermolekuelmasse);
}

double ParticleSystem::CalkritischeMolWasser1(Particle p) {
	double mitteldichte = 0.0;
	double durchmesser = p.GetDiameter();

	mitteldichte = (1 - p.GetZusammensetzung()) * WasserDicht + p.GetZusammensetzung() * GlycerinDicht;

	return M_PI * pow(durchmesser, 3) * mitteldichte * (1 - p.GetZusammensetzung()) / 6 / AnzahlproMol
			/ ((GlycerinMolekuelmasse - Wassermolekuelmasse) * p.GetZusammensetzung() + Wassermolekuelmasse);
}

unsigned int ParticleSystem::CalWasserNumber(double durchmesser) {
	unsigned int result = 0;
	double mitteldichte = 0.0;

	mitteldichte = (1 - kritischGlycerinZusammensetzung) * WasserDicht + kritischGlycerinZusammensetzung * GlycerinDicht;

	result = ceil(
			M_PI * pow(durchmesser, 3) * mitteldichte * (1 - kritischGlycerinZusammensetzung) / 6
					/ ((GlycerinMolekuelmasse - Wassermolekuelmasse) * kritischGlycerinZusammensetzung + Wassermolekuelmasse));

	if (result <= 1) {
		return 1;
	} else {
		return result;
	}
}

Particle ParticleSystem::GetNewBuiltParticle() {
	return tempParticle;
}

bool ParticleSystem::CheckNewBuiltParticle() {
	double durchmesser = CalKeimGroesse();
	unsigned int numberOfWasser = CalWasserNumber(durchmesser);
	unsigned int numberOfGlycerin = CalGlycerinNumber(durchmesser);
	double kritischMolWasser = CalkritischeMolWasser();
	double kritischMolGlycerin = CalkritischeMolGlyzerin();

	durchmesser = (numberOfGlycerin * GlycerinMolekuelmasse + numberOfWasser * Wassermolekuelmasse) * 6 / M_PI
			/ ((1 - kritischGlycerinZusammensetzung) * WasserDicht + kritischGlycerinZusammensetzung * GlycerinDicht);
	durchmesser = pow(durchmesser, 0.333333);

	if (durchmesser > 0 && numberOfWasser > 0 && numberOfGlycerin > 0 && kritischGlycerinZusammensetzung > 0 && kritischMolWasser > 0
			&& kritischMolGlycerin > 0) {
		tempParticle = Particle(durchmesser, kritischGlycerinZusammensetzung, kritischMolWasser, kritischMolGlycerin, numberOfWasser, numberOfGlycerin);
		return true;
	} else {
		return false;
	}
}

double ParticleSystem::CalKeimbildungsrate() {
	double deltaG = 0.0;
	double sin = 0.0;
	double cos = 0.0;
	double wasserDruck = CalWasserPartielDruck();
	double glycerinDruck = CalGlycerinPartielDruck();
	double wasserAnteil = 0.0;
	double glycerinAnteil = 0.0;
	double mitteloberflaeche = 0.0;
	double mittelvolumen = 0.0;

	if (kritischGlycerinZusammensetzung == 0) {
		return 0.0;
	} else {
		sin = 1 / pow((1 + pow(kritischGlycerinZusammensetzung, -2)), 0.5);
		cos = 1 / pow((1 + pow(kritischGlycerinZusammensetzung, 2)), 0.5);

		wasserAnteil = wasserDruck / pow((2 * M_PI * Wassermolekuelmasse * BoltzmannConstant * UmgebungsTemperatur), 0.5);
		glycerinAnteil = glycerinDruck / pow((2 * M_PI * GlycerinMolekuelmasse * BoltzmannConstant * UmgebungsTemperatur), 0.5);

		mitteloberflaeche = (1 - kritischGlycerinZusammensetzung) * WasserOberflaeche + kritischGlycerinZusammensetzung * GlycerinOberflaeche;
		mittelvolumen = (1 - kritischGlycerinZusammensetzung) * WasserMolekuelVol + kritischGlycerinZusammensetzung * GlycerinMolekuelVol;

		deltaG = M_PI * mitteloberflaeche * pow(CalKeimGroesse(), 2) / 3;

		return  1e-3*2 * mittelvolumen * wasserAnteil * glycerinAnteil / (wasserAnteil * pow(sin, 2) + glycerinAnteil * pow(cos, 2))
				* ((wasserDruck + glycerinDruck) / BoltzmannConstant / UmgebungsTemperatur)
				* pow((mitteloberflaeche / UmgebungsTemperatur / BoltzmannConstant), 0.5) * exp(-deltaG / BoltzmannConstant / UmgebungsTemperatur);
	}
}

// Aktualisierung nach Keimbildung;
double ParticleSystem::CalKeimnewDampfmasseGlycerin() {
	return glycerinDampfMasse - CalkritischeMolGlyzerin() * GlycerinMolMasse;
}

double ParticleSystem::CalKeimnewDampfmasseWasser() {
	return wasserDampfMasse - CalkritischeMolWasser() * WasserMolMasse;
}

double ParticleSystem::CalkeimnewpartiellDruckGlycerin() {
	double wasserAnteil = CalKeimnewDampfmasseWasser() / WasserMolMasse;
	double glycerinAnteil = CalKeimnewDampfmasseGlycerin() / GlycerinMolMasse;
	double stickStoffAnteil = stickStoffDampfMasse / StickStoffMolMasse;

	return glycerinAnteil / (glycerinAnteil + wasserAnteil + stickStoffAnteil) * GesamtDruck;
}

double ParticleSystem::CalkeimnewpartiellDruckWasser() {
	double wasserAnteil = CalKeimnewDampfmasseWasser() / WasserMolMasse;
	double glycerinAnteil = CalKeimnewDampfmasseGlycerin() / GlycerinMolMasse;
	double stickStoffAnteil = stickStoffDampfMasse / StickStoffMolMasse;

	return wasserAnteil / (glycerinAnteil + wasserAnteil + stickStoffAnteil) * GesamtDruck;
}

//Die folgenden Gleichungen gehören zum Wachstum (Nach Keimbildung);
double ParticleSystem::CalzunehmenMolekuelGlycerin1(Particle p) {
	double result = 0.0;
	double mitteloberflaeche = 0.0;
	double PdG = 0.0;
	double glycerinteil = 0.0;
	double glycerinpartialdruck = CalGlycerinPartielDruck();

	mitteloberflaeche = (1 - p.GetZusammensetzung()) * WasserOberflaeche + p.GetZusammensetzung() * GlycerinOberflaeche;

	PdG = p.GetZusammensetzung() * GlycerinUebersaettigungsdruck
			* exp(4 * mitteloberflaeche * GlycerinMolMasse / GlycerinDicht / UmgebungsTemperatur / Gasconstant / p.GetDiameter());

	glycerinteil = 2 * M_PI * GlycerinMolekuelmasse * BoltzmannConstant * UmgebungsTemperatur;
	result = M_PI * pow(p.GetDiameter(), 2) * ((glycerinpartialdruck - PdG) / pow(glycerinteil, 0.5));

	if (result >= 0) {
		return result;
	} else {
		return 0;
	}
}

double ParticleSystem::CalzunehmenMolekuelWasser1(Particle p) {
	double result = 0.0;
	double mitteloberflaeche = 0.0;
	double PwG = 0.0;
	double Wasserteil = 0.0;
	double Wasserpartialdruck = CalWasserPartielDruck();

	mitteloberflaeche = (1 - p.GetZusammensetzung()) * WasserOberflaeche + p.GetZusammensetzung() * GlycerinOberflaeche;

	PwG = (1 - p.GetZusammensetzung()) * WasserUebersaettigungsdruck
			* exp(4 * mitteloberflaeche * WasserMolMasse/WasserDicht / UmgebungsTemperatur / Gasconstant / p.GetDiameter());

	Wasserteil = 2 * M_PI * Wassermolekuelmasse * BoltzmannConstant * UmgebungsTemperatur;

	result = M_PI * pow(p.GetDiameter(), 2) * ((Wasserpartialdruck - PwG) / pow(Wasserteil, 0.5));

	if (result >= 0) {
		return result;
	} else {
		return 0;
	}
}

void ParticleSystem::CalWachstumOnlyNewParticle(Particle &p, double deltaT) {
	double zunehmenMolekuelWasser = CalzunehmenMolekuelWasser1(p);
	double zunehmenMolekuelGlycerin = CalzunehmenMolekuelGlycerin1(p);
//	cout<<"!!!!!"<<zunehmenMolekuelWasser<<" "<<zunehmenMolekuelGlycerin<<"!!!!!!!!!"<<endl;
	double oldGlycerinMolekuel = p.GetNumberOfGlycerin();
	double oldWasserMolekuel = p.GetNumberOfWasser();
//	cout<<"!!!!!"<<oldGlycerinMolekuel<<" "<<oldWasserMolekuel<<"!!!!!!!!!"<<endl;
	double newGlycerinMolekuel = oldGlycerinMolekuel
			+ zunehmenMolekuelGlycerin * deltaT;
	double newWasserMolekuel = oldWasserMolekuel
			+ zunehmenMolekuelWasser * deltaT;
	double newGlycerinMol = newGlycerinMolekuel / AnzahlproMol;
	double newWasserMol = newWasserMolekuel / AnzahlproMol;
	double newZusammensetzung = newGlycerinMol
			/ (newGlycerinMol + newWasserMol);
	double altmitteldichte = (1 - p.GetZusammensetzung()) * WasserDicht
			+ p.GetZusammensetzung() * GlycerinDicht;
	double newDicht = newZusammensetzung * GlycerinDicht
			+ (1 - newZusammensetzung) * WasserDicht;
	double newDurchmesser = ((zunehmenMolekuelGlycerin * GlycerinMolekuelmasse
			+ zunehmenMolekuelWasser * Wassermolekuelmasse) * deltaT
			+ M_PI * pow(p.GetDiameter(), 3) * altmitteldichte / 6) * 6 / M_PI
			/ newDicht;
	newDurchmesser = pow(newDurchmesser, 0.333333);

	cout << newDurchmesser << "  " << p.GetDiameter() << endl;

	p.SetGlycerinMol(newGlycerinMol);
	p.SetWasserMol(newWasserMol);
	p.SetNumberOfGlycerin(newGlycerinMolekuel);
	p.SetNumberOfWasser(newWasserMolekuel);
}

void ParticleSystem::CalWachstumNewParticle(Particle &p, double deltaT) {
	double zunehmenMolekuelWasser = 0.1*CalzunehmenMolekuelWasser1(p);
	double zunehmenMolekuelGlycerin = 0.1*CalzunehmenMolekuelGlycerin1(p);
	double oldGlycerinMolekuel = p.GetNumberOfGlycerin();
	double oldWasserMolekuel = p.GetNumberOfWasser();
	double newGlycerinMolekuel = oldGlycerinMolekuel + zunehmenMolekuelGlycerin * deltaT;
	double newWasserMolekuel = oldWasserMolekuel + zunehmenMolekuelWasser * deltaT;
	double newGlycerinMol = newGlycerinMolekuel / AnzahlproMol;
	double newWasserMol = newWasserMolekuel / AnzahlproMol;
	double newZusammensetzung = newGlycerinMol / (newGlycerinMol + newWasserMol);
	double altmitteldichte = (1 - p.GetZusammensetzung()) * WasserDicht + p.GetZusammensetzung() * GlycerinDicht;
	double newDicht = newZusammensetzung * GlycerinDicht + (1 - newZusammensetzung) * WasserDicht;
	double newDurchmesser = ((zunehmenMolekuelGlycerin * GlycerinMolekuelmasse + zunehmenMolekuelWasser * Wassermolekuelmasse) * deltaT
			+ M_PI * pow(p.GetDiameter(), 3) * altmitteldichte / 6) * 6 / M_PI / newDicht;
	newDurchmesser = pow(newDurchmesser, 0.333333);

	if (newDurchmesser > p.GetDiameter()) {
		p.SetDiameter(newDurchmesser);
		p.SetGlycerinMol(newGlycerinMol);
		p.SetWasserMol(newWasserMol);
		p.SetNumberOfGlycerin(newGlycerinMolekuel);
		p.SetNumberOfWasser(newWasserMolekuel);
		p.SetZusammensetzung(newZusammensetzung);
	//	cout<<newWasserMolekuel<<" "<<newDurchmesser<<endl;
		//cout<<"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"<<endl;
	}
}

double ParticleSystem::CalWachstumneuMolGlycerin1(Particle p, double deltaT) {
	double result = 0.0;

	result = p.GetGlycerinMol() + CalzunehmenMolekuelGlycerin1(p) * deltaT / AnzahlproMol;

	if (result >= p.GetGlycerinMol()) {
		return result;
	} else {
		return p.GetGlycerinMol();
	}
}

double ParticleSystem::CalWachstumneuMolWasser1(Particle p, double deltaT) {
	double result = 0.0;

	result = p.GetWasserMol() + CalzunehmenMolekuelWasser1(p) * deltaT / AnzahlproMol;

	if (result >= p.GetWasserMol()) {
		return result;
	} else {
		return p.GetWasserMol();
	}
}

double ParticleSystem::CalWachstumneuZusammensetzungGlycerin1(Particle p, double deltaT) {
	double result = 0.0;
	double Glycerinmol = CalWachstumneuMolGlycerin1(p, deltaT);
	double Wassermol = CalWachstumneuMolWasser1(p, deltaT);

	result = Glycerinmol / (Glycerinmol + Wassermol);

	if (result >= 0) {
		return result;
	} else {
		return p.GetZusammensetzung();
	}
}

double ParticleSystem::CalWachstumneuGrosse1(Particle p, double deltaT) {
	double result = 0.0;
	double altmitteldichte = 0.0;
	double neuzusammensetzung = CalWachstumneuZusammensetzungGlycerin1(p, deltaT);
	double neumitteldichte = 0.0;
	double zunehmendMolekuel = 0.0;

	altmitteldichte = (1 - p.GetZusammensetzung()) * WasserDicht + p.GetZusammensetzung() * GlycerinDicht;

	neumitteldichte = (1 - neuzusammensetzung) * WasserDicht + neuzusammensetzung * GlycerinDicht;

	zunehmendMolekuel = 6 * (CalzunehmenMolekuelGlycerin1(p) * GlycerinMolekuelmasse + CalzunehmenMolekuelWasser1(p) * Wassermolekuelmasse) * deltaT / M_PI;

	result = pow((pow(CalKeimGroesse(), 3) * altmitteldichte + zunehmendMolekuel) / neumitteldichte, 0.33333);

	if (result >= p.GetDiameter()) {
		return result;
	} else {
		return p.GetDiameter();
	}
}

//Berechnung für Wachstumsüberprüfung (nach Keimbildung)

double ParticleSystem::Calmassenaenderungsrate(Particle p, double deltaT) {
	double result = 0.0;
	double deltam = 0.0;
	double altmitteldichte = 0.0;

	altmitteldichte = (1 - p.GetZusammensetzung()) * WasserDicht + p.GetZusammensetzung() * GlycerinDicht;

	deltam = (CalzunehmenMolekuelGlycerin1(p) * GlycerinMolekuelmasse + CalzunehmenMolekuelWasser1(p) * Wassermolekuelmasse) * deltaT;

	result = 6 * deltam / (M_PI * pow(p.GetDiameter(), 3) * altmitteldichte);

	return result;
}
