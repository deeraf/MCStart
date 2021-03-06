/*
 * Logger.cpp
 *
 *  Created on: May 18, 2013
 *      Author: amifan
 */

#include "Logger.h"
#include "string.h"
#include <algorithm>
#include <cmath>
#include <math.h>

using namespace std;

Logger::Logger() {
}

Logger::~Logger() {
}

/*void Logger::DrawPlot(double *x, double *qx, int n, int name, double time) {
	mglGraph gr;
	mglData x_mgl;
	mglData qx_mgl;

	double ymax = 0;

	for (int i = 0; i < n; i++) {
		if (qx[i] > ymax) {
			ymax = qx[i];
		}
	}

	gr.SetRanges(mglPoint(0.01, 0), mglPoint(1000, ymax));
	gr.SetQuality(2);
	gr.SetSize(1600, 1200);
	gr.SetFunc("lg(x)", "", 0);
	gr.SetTicks('x', 0);
	gr.SetTicks('y', 0);
	gr.Box();
	x_mgl.Set(x, n);
	qx_mgl.Set(qx, n);
	gr.Plot(x_mgl, qx_mgl, "");
	gr.Axis();
	gr.Grid("xy", "g;");
	gr.Title(GetCharString(time));
	gr.WritePNG(GetCharString(name, ".png"), "", false);
}*/

char *Logger::GetCharString(double i, string ext) {
	stringstream index;
	index << i;

	string fileNameStr = index.str();

	if (i < 10) {
		fileNameStr.insert(0, 1, '0');
	}

	string extension = ext;
	int fileNameSize = 4 + fileNameStr.length();
	char *fileName = (char*) malloc(sizeof(char) * fileNameSize);
	for (unsigned int i = 0; i < fileNameStr.length(); i++) {
		fileName[i] = fileNameStr[i];
	}

	for (unsigned int j = 0; j <= extension.size(); j++) {
		fileName[j + fileNameStr.length()] = extension[j];
	}

	return fileName;
}

char *Logger::GetCharString(double i) {
	stringstream is;

	if (i == 0.0) {
		is << "Initialize";
	} else {
		is << i;
	}

	string str = is.str();

	char *charStr = (char*) malloc(sizeof(char) * str.length());
	for (unsigned int i = 0; i < str.length(); i++) {
		charStr[i] = str[i];
	}

	return charStr;
}

void Logger::PrintTraceCoagulation(double numberOfParticles, double numberOfCoagulation, double ithDiameter, double jthDiameter, Particle newParticle,
		double nextCriterion, double nextRandom, double S, double currentTime) {
	cout << "Number of all particles: " << numberOfParticles << " Number of coagulations: " << numberOfCoagulation << " Di: " << ithDiameter << " Dj: "
			<< jthDiameter << " New diameter: " << newParticle.GetDiameter() << " Water molecular: " << newParticle.GetNumberOfWasser()
			<< " Glycerin molecular: " << newParticle.GetNumberOfGlycerin() << " Criterion: " << nextCriterion << " Random: " << nextRandom << " S: " << S
			<< " Current time: " << currentTime << endl;
}

void Logger::PrintTraceKeimbildung(double numberOfParticles, double numberOfKeimbildung, double diameter, double zusammensetzung, double keimbildungsrate,
		double currentTime) {
	cout << "Number of all particles: " << numberOfParticles << " Number of Keimbildung: " << numberOfKeimbildung << " Diameter: " << diameter
			<< " Zusammensetzung: " << zusammensetzung << " Keimbildungsrate: " << keimbildungsrate << " Current time: " << currentTime << endl;
}

void Logger::WriteToFileParticle(ParticleClassVector pcv, double volumen, int name) {
	ofstream outputFile;
	ofstream outputFile2;

	outputFile.open(GetCharString(name, ".txt"), ios::app);
	outputFile2.open("0_MCOutput.txt");

	outputFile << "Current volume: " << volumen << endl;

	for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
		outputFile2 << pcv.GetParticleClassAt(i).GetSize() << endl;
		for (unsigned int j = 0; j < pcv.GetParticleClassAt(i).GetSize(); j++) {
			outputFile << pcv.GetParticleClassAt(i).GetParticleAt(j).GetDiameter() << "	" << pcv.GetParticleClassAt(i).GetParticleAt(j).GetZusammensetzung()
					<< endl;
		}
	}

	outputFile.close();
	outputFile2.close();
}

void Logger::WriteToFileParticleCount(ParticleClassVector pcv, double volumen, int name) {
	ofstream outputFile;

	outputFile.open(GetCharString(name, "_C.txt"), ios::app);

	outputFile << "Current volume: " << volumen << endl;

	vector<Particle> tempPV;

	for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
		for (unsigned int j = 0; j < pcv.GetParticleClassAt(i).GetSize(); j++) {
			tempPV.push_back(pcv.GetParticleClassAt(i).GetParticleAt(j));
		}
	}

	sort(tempPV.begin(), tempPV.end(), pcv.CompareByDiameter);

	double currentParticleDiameter = tempPV.at(0).GetDiameter();
	unsigned int currentParticleCount = 0;

	for (unsigned int k = 0; k < tempPV.size(); k++) {
		if (abs(tempPV.at(k).GetDiameter() - currentParticleDiameter) < 1E-10) {
			currentParticleCount ++;
		} else {
			outputFile<<currentParticleDiameter<<"	"<<currentParticleCount<<endl;
			if (k != tempPV.size() -1){
				currentParticleDiameter = tempPV.at(k).GetDiameter();
			}
			currentParticleCount = 1;
		}
	}

	outputFile.close();
}

void Logger::WriteToFileCoagulation(double numberOfParticles, double numberOfCoagulation, double ithDiameter, double jthDiameter, Particle newParticle,
		double nextCriterion, double nextRandom, double S, double currentTime) {
	ofstream outputFile;

	/*outputFile.open("AllDebugInformation.txt", ios::app);

	 outputFile << "Number of all particles: " << numberOfParticles << " Number of coagulations: " << numberOfCoagulation
	 << " Di: " << ithDiameter << " Dj: " << jthDiameter << " New diameter: " << newParticle.GetDiameter()
	 << " Water molecular: " << newParticle.GetNumberOfWasser() << " Glycerin molecular: "
	 << newParticle.GetNumberOfGlycerin() << " Criterion: " << nextCriterion << " Random: " << nextRandom
	 << " S: " << S << " Current time: " << currentTime << endl;

	 outputFile.close();*/

	outputFile.open("AllCoagulationInformation.txt", ios::app);

	outputFile << "Number of all particles: " << numberOfParticles << " Number of coagulations: " << numberOfCoagulation << " Di: " << ithDiameter << " Dj: "
			<< jthDiameter << " New diameter: " << newParticle.GetDiameter() << " Water molecular: " << newParticle.GetNumberOfWasser()
			<< " Glycerin molecular: " << newParticle.GetNumberOfGlycerin() << " Criterion: " << nextCriterion << " Random: " << nextRandom << " S: " << S
			<< " Current time: " << currentTime << endl;

	outputFile.close();
}

void Logger::WriteToFileKeimbildung(double numberOfParticles, double numberOfKeimbildung, double diameter, double zusammensetzung, double keimbildungsrate,
		double currentTime) {
	ofstream outputFile;

	/*outputFile.open("AllDebugInformation.txt", ios::app);

	 outputFile << "Number of all particles: " << numberOfParticles << " Number of Keimbildung: " << numberOfKeimbildung
	 << " Diameter: " << diameter << " Zusammensetzung: " << zusammensetzung << " Keimbildungsrate: "
	 << keimbildungsrate << " Current time: " << currentTime << endl;

	 outputFile.close();*/

	outputFile.open("AllKeimbildungInformation.txt", ios::app);

	outputFile << "Number of all particles: " << numberOfParticles << " Number of Keimbildung: " << numberOfKeimbildung << " Diameter: " << diameter
			<< " Zusammensetzung: " << zusammensetzung << " Keimbildungsrate: " << keimbildungsrate << " Current time: " << currentTime << endl;

	outputFile.close();
}

void Logger::WriteToFileWachstum(double deltaWasserDampMasse, double deltaGlycerinDampMasse, double wasserDampMasse, double glycerinDampMasse,
		double currentTime) {
	ofstream outputFile;

	/*outputFile.open("AllDebugInformation.txt", ios::app);

	 outputFile << " Delta WasserDampMasse: " << deltaWasserDampMasse
	 << " Delta GlycerinDampMasse: " << deltaGlycerinDampMasse << " Current WasserDampMasse: " << wasserDampMasse
	 << " Current GlycerinDampMasse: " << glycerinDampMasse << " Current time: " << currentTime << endl;

	 outputFile.close();*/

	outputFile.open("MasseAfterWachstum.txt", ios::app);

	outputFile << " Delta WasserDampMasse: " << deltaWasserDampMasse << " Delta GlycerinDampMasse: " << deltaGlycerinDampMasse << " Current WasserDampMasse: "
			<< wasserDampMasse << " Current GlycerinDampMasse: " << glycerinDampMasse << " Current time: " << currentTime << endl;

	outputFile.close();
}

void Logger::WriteToFileKeimbildungMasse(double deltaWasserDampMasse, double deltaGlycerinDampMasse, double wasserDampMasse, double glycerinDampMasse,
		double currentTime) {
	ofstream outputFile;

	outputFile.open("MasseAfterKeimbildung.txt", ios::app);

	outputFile << " Delta WasserDampMasse: " << deltaWasserDampMasse << " Delta GlycerinDampMasse: " << deltaGlycerinDampMasse << " Current WasserDampMasse: "
			<< wasserDampMasse << " Current GlycerinDampMasse: " << glycerinDampMasse << " Current time: " << currentTime << endl;

	outputFile.close();
}

void Logger::WriteToFileMittelWert(double zusammensetzung, double diameter, double numberOfAllParticles, double volumen, double currentTime) {
	ofstream outputFile;

	outputFile.open("Mittelwert.txt", ios::app);

	outputFile << zusammensetzung << "	" << diameter << "	" << numberOfAllParticles << "	" << volumen << "	" << currentTime << endl;

	outputFile.close();
}

void Logger::WriteToFileDebug(ParticleClassVector pcv, int name) {
	ofstream outputFile;

	outputFile.open("Debug.txt", ios::app);

	for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
		for (unsigned int j = 0; j < pcv.GetParticleClassAt(i).GetSize(); j++) {
			outputFile << pcv.GetParticleClassAt(i).GetParticleAt(j).GetDiameter() << "	"
					<< pcv.GetParticleClassAt(i).GetParticleAt(j).GetZusammensetzung() << endl;
		}
	}

	outputFile.close();
}

void Logger::WriteToFile3D(ParticleClassVector pcv, int name) {
	ofstream outputFile;

	vector<Particle> tempPV;

	for (unsigned int i = 0; i < pcv.GetNumberOfParticleClasses(); i++) {
		for (unsigned int j = 0; j < pcv.GetParticleClassAt(i).GetSize(); j++) {
			tempPV.push_back(pcv.GetParticleClassAt(i).GetParticleAt(j));
		}
	}

	sort(tempPV.begin(), tempPV.end(), pcv.CompareByDiameter);
	double maxDiameter = tempPV.at(tempPV.size() - 1).GetDiameter();
	unsigned int numberOfClass = maxDiameter / 1E-9 + 1;
	unsigned int numberOfClass1 = 0;
	unsigned int numberOfClass2 = 0;
	unsigned int numberOfClass3 = 0;
	unsigned int numberOfClass4 = 0;
	unsigned int numberOfClass5 = 0;
	unsigned int numberOfClass6 = 0;
	unsigned int numberOfClass7 = 0;
	unsigned int numberOfClass8 = 0;
	unsigned int numberOfClass9 = 0;
	unsigned int numberOfClass10 = 0;
	unsigned int numberOfClass11 = 0;
	unsigned int numberOfClass12 = 0;
	unsigned int numberOfClass13 = 0;
	unsigned int numberOfClass14 = 0;
	unsigned int numberOfClass15 = 0;
	unsigned int numberOfClass16 = 0;
	unsigned int numberOfClass17 = 0;
	unsigned int numberOfClass18 = 0;
	unsigned int numberOfClass19 = 0;
	unsigned int numberOfClass20 = 0;
	unsigned int numberOfClass21 = 0;
	unsigned int numberOfClass22 = 0;
	unsigned int numberOfClass23 = 0;
	unsigned int numberOfClass24 = 0;
	unsigned int numberOfClass25 = 0;
	unsigned int numberOfClass26 = 0;
	unsigned int numberOfClass27 = 0;
	unsigned int numberOfClass28 = 0;
	unsigned int numberOfClass29 = 0;
	unsigned int numberOfClass30 = 0;
	unsigned int numberOfClass31 = 0;
	unsigned int numberOfClass32 = 0;
	unsigned int numberOfClass33 = 0;
	unsigned int numberOfClass34 = 0;
	unsigned int numberOfClass35 = 0;
	unsigned int numberOfClass36 = 0;
	unsigned int numberOfClass37 = 0;
	unsigned int numberOfClass38 = 0;
	unsigned int numberOfClass39 = 0;
	unsigned int numberOfClass40 = 0;
	unsigned int numberOfClass41 = 0;
	unsigned int numberOfClass42 = 0;
	unsigned int numberOfClass43 = 0;
	unsigned int numberOfClass44 = 0;
	unsigned int numberOfClass45 = 0;
	unsigned int numberOfClass46 = 0;
	unsigned int numberOfClass47 = 0;
	unsigned int numberOfClass48 = 0;
	unsigned int numberOfClass49 = 0;
	unsigned int numberOfClass50 = 0;

	unsigned int currentIndex = 0;

	outputFile.open(GetCharString(name, "_3D.txt"), ios::app);

	outputFile << "	d\\z	" << "	0.01	" << "	0.02	" << "	0.03	" << "	0.04	" << "	0.05	" << "	0.06	" << "	0.07	" << "	0.08	" << "	0.09	"
			<< "	0.10	" << "	0.11	" << "	0.12	" << "	0.13	" << "	0.14	" << "	0.15	" << "	0.16	" << "	0.17	" << "	0.18	" << "	0.19	"
			<< "	0.20	" << "	0.21	" << "	0.22	" << "	0.23	" << "	0.24	" << "	0.25	" << "	0.26	" << "	0.27	" << "	0.28	" << "	0.29	"
			<< "	0.30	" << "	0.31	" << "	0.32	" << "	0.33	" << "	0.34	" << "	0.35	" << "	0.36	" << "	0.37	" << "	0.38	" << "	0.39	"
			<< " 	0.40	" << " 	0.41	" << " 	0.42	" << " 	0.43	" << " 	0.44	" << "	0.45	" << " 	0.46	" << " 	0.47	" << " 	0.48	" << " 	0.49	"
			<< "	rest	" << endl;

	for (unsigned int i = 0; i < numberOfClass; i++) {
		double leftBoundary = 1E-9 * i;
		double rightBoundary = 1E-9 * (i + 1);
		outputFile << leftBoundary << " - " << rightBoundary << "	";

		while (tempPV.at(currentIndex).GetDiameter() >= leftBoundary && tempPV.at(currentIndex).GetDiameter() <= rightBoundary) {
			if (tempPV.at(currentIndex).GetZusammensetzung() > 0 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.01) {
				numberOfClass1++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.01 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.02) {
				numberOfClass2++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.02 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.03) {
				numberOfClass3++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.03 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.04) {
				numberOfClass4++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.04 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.05) {
				numberOfClass5++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.05 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.06) {
				numberOfClass6++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.06 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.07) {
				numberOfClass7++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.07 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.08) {
				numberOfClass8++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.08 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.09) {
				numberOfClass9++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.09 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.10) {
				numberOfClass10++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.10 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.11){
				numberOfClass11++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.11 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.12) {
				numberOfClass12++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.12 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.13) {
				numberOfClass13++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.13 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.14) {
				numberOfClass14++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.14 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.15) {
				numberOfClass15++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.15 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.16) {
				numberOfClass16++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.16 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.17) {
				numberOfClass17++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.17 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.18) {
				numberOfClass18++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.18 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.19) {
				numberOfClass19++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.19 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.20) {
				numberOfClass20++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.20 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.21){
				numberOfClass21++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.21 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.22) {
				numberOfClass22++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.22 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.23) {
				numberOfClass23++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.23 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.24) {
				numberOfClass24++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.24 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.25) {
				numberOfClass25++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.25 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.26) {
				numberOfClass26++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.26 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.27) {
				numberOfClass27++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.27 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.28) {
				numberOfClass28++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.28 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.29) {
				numberOfClass29++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.29 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.30) {
				numberOfClass30++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.30 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.31){
				numberOfClass31++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.31 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.32) {
				numberOfClass32++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.32 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.33) {
				numberOfClass33++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.33 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.34) {
				numberOfClass34++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.34 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.35) {
				numberOfClass35++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.35 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.36) {
				numberOfClass36++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.36 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.37) {
				numberOfClass37++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.37 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.38) {
				numberOfClass38++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.38 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.39) {
				numberOfClass39++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.39 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.40) {
				numberOfClass40++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.40 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.41){
				numberOfClass41++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.41 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.42) {
				numberOfClass42++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.42 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.43) {
				numberOfClass43++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.43 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.44) {
				numberOfClass44++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.44 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.45) {
				numberOfClass45++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.45 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.46) {
				numberOfClass46++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.46 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.47) {
				numberOfClass47++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.47 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.48) {
				numberOfClass48++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.48 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.49) {
				numberOfClass49++;
			} else if (tempPV.at(currentIndex).GetZusammensetzung() > 0.49 && tempPV.at(currentIndex).GetZusammensetzung() <= 0.50) {
				numberOfClass50++;
			} else {
				numberOfClass++;
			}

			if (currentIndex < tempPV.size() - 1) {
				currentIndex++;
			} else {
				break;
			}
		}

		outputFile << numberOfClass1 << "	" << numberOfClass2 << "	" << numberOfClass3 << "	" << numberOfClass4 << "	" << numberOfClass5 << "	" << numberOfClass6 << "	" << numberOfClass7 << "	" << numberOfClass8 << "	" << numberOfClass9 << "	" << numberOfClass10
				<< numberOfClass11 << "	" << numberOfClass12 << "	" << numberOfClass13 << "	" << numberOfClass14 << "	" << numberOfClass15 << "	" << numberOfClass16 << "	" << numberOfClass17 << "	" << numberOfClass18 << "	" << numberOfClass19 << "	" << numberOfClass20
				<< numberOfClass21 << "	" << numberOfClass22 << "	" << numberOfClass23 << "	" << numberOfClass24 << "	" << numberOfClass25 << "	" << numberOfClass26 << "	" << numberOfClass27 << "	" << numberOfClass28 << "	" << numberOfClass29 << "	" << numberOfClass30
				<< numberOfClass31 << "	" << numberOfClass32 << "	" << numberOfClass33 << "	" << numberOfClass34 << "	" << numberOfClass35 << "	" << numberOfClass36 << "	" << numberOfClass37 << "	" << numberOfClass38 << "	" << numberOfClass39 << "	" << numberOfClass40
				<< numberOfClass41 << "	" << numberOfClass42 << "	" << numberOfClass43 << "	" << numberOfClass44 << "	" << numberOfClass45 << "	" << numberOfClass46 << "	" << numberOfClass47 << "	" << numberOfClass48 << "	" << numberOfClass49 << "	" << numberOfClass50
				<< endl;

		numberOfClass1 = 0;
		numberOfClass2 = 0;
		numberOfClass3 = 0;
		numberOfClass4 = 0;
		numberOfClass5 = 0;
		numberOfClass6 = 0;
		numberOfClass7 = 0;
		numberOfClass8 = 0;
		numberOfClass9 = 0;
		numberOfClass10 = 0;
		numberOfClass11 = 0;
		numberOfClass12 = 0;
		numberOfClass13 = 0;
		numberOfClass14 = 0;
		numberOfClass15 = 0;
		numberOfClass16 = 0;
		numberOfClass17 = 0;
		numberOfClass18 = 0;
		numberOfClass19 = 0;
		numberOfClass20 = 0;
		numberOfClass21 = 0;
		numberOfClass22 = 0;
		numberOfClass23 = 0;
		numberOfClass24 = 0;
		numberOfClass25 = 0;
		numberOfClass26 = 0;
		numberOfClass27 = 0;
		numberOfClass28 = 0;
		numberOfClass29 = 0;
		numberOfClass30 = 0;
		numberOfClass31 = 0;
		numberOfClass32 = 0;
		numberOfClass33 = 0;
		numberOfClass34 = 0;
		numberOfClass35 = 0;
		numberOfClass36 = 0;
		numberOfClass37 = 0;
		numberOfClass38 = 0;
		numberOfClass39 = 0;
		numberOfClass40 = 0;
		numberOfClass41 = 0;
		numberOfClass42 = 0;
		numberOfClass43 = 0;
		numberOfClass44 = 0;
		numberOfClass45 = 0;
		numberOfClass46 = 0;
		numberOfClass47 = 0;
		numberOfClass48 = 0;
		numberOfClass49 = 0;
		numberOfClass50 = 0;


	}

	outputFile.close();
}
