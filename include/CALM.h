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
//#include <fstream>
//--------reggae
#include "reggae.h"
#include "specrel.h"

using namespace std;

enum eEventType {GLOBAL, MINIJETS_GLOBAL, MINIJETS_LOCAL,GLOBAL_REGGAE, MINIJETS_GLOBAL_REGGAE, MINIJETS_LOCAL_REGGAE};

class CALM {
	public:
		CALM();
		~CALM();
		int   GenerateParticles(ParticleDB* aPartDB,int aMultBinMin, int aMultBinMax, std::list<Particle>* aParticles, eEventType aEventType = GLOBAL);
		void   Randomize();
	private:
      void ReadParameters();
		TRandom2* mRandom;
      // values constant (what can be generated etc)
		int		 mNpart; // number of kinds of particles (pions, kaons, nuclides etc.)
		double* mNmean; // mean values for the particle yields dN/dy of each particle count -> taken from ALICE
		double mRapidityInterval; //interval of rapidity
		double* mXYZ; // to generate distance from primary vertex 
		int* mNpartkinds; // number of particles for each kind
		string** mNames; // names of particles to be generated
		TF1* Ptot; //Total momentum to be distributed among particles
      // values for this event
      vector<string> mParticlesThisEvent;
      int mMultMin;  //min multiplicity for CALM
      int mMultMax;  //max multiplicity for CALM
      eEventType mEventType;
};

#endif
