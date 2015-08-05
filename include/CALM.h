#ifndef _CALM_H_
#define _CALM_H_

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include "TRandom2.h"
#include "TGenPhaseSpace.h"
#include "ParticleType.h"
#include "ParticleDB.h"
#include "Particle.h"
#include "THGlobal.h"
#include "TF1.h"
//--------reggae
#include "reggae.h"
#include "specrel.h"

//#include <fstream>

using namespace std;

enum eEventType {GLOBAL, MINIJETS_GLOBAL, MINIJETS_LOCAL};

class CALM {
	public:
		CALM();
		~CALM();
		int   GenerateParticles(ParticleDB* aPartDB,int aMultBinMin, int aMultBinMax, double aEnergy, std::list<Particle>* aParticles, eEventType aEventType = GLOBAL);
		void   SetMultiplicities(ParticleDB* aDB,int aMultBinMin, int aMultBinMax);
		void   Randomize();
	private:
		TRandom2* mRandom;
      // values constant (what can be generated etc)
		int		 mNpart; // number of kinds of particles (pions, kaons, nuclides etc.)
		double* mNmean; // mean values for the particle yields dN/dy of each particle count -> taken from ALICE
		double mRapidityInterval; //interval of rapidity
		double* mXYZ; // to generate distance from primary vertex 
		int* mNpartkinds; // number of particles for each kind
		string** mNames; // names of particles to be generated
      // values for this event
      vector<string> mParticlesThisEvent; 
};

#endif
