/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#ifndef _TH2_EVENT_H_
  #define _TH2_EVENT_H_

#include <list>
#include <vector>
#include <TRandom2.h>
#include "ParticleDB.h"
#include "Particle.h"
#include "CALM.h"
#include <sstream>
#include <string>

class Event {

   public:
      Event();
      Event(ParticleDB* aDB, CALM* aCALM);
      ~Event();

      void      Reset(int aEventIter=0);
      std::list<Particle>* GetParticleList();
      ParticleDB*     GetParticleDB() const;
      unsigned int    GetEventID() const;

      void      GeneratePrimordials(int aSeed=0);
      void      Randomize();

   private:
      void ReadParameters();

      std::list<Particle>	mParticles;
      std::vector<int>	mMultiplicities;
      unsigned int	mEventID;
      ParticleDB*		mPartDB;
      CALM*		mCALM;
      TRandom2*		mRandom;
      int mMultMin;  //min multiplicity for CALM
      int mMultMax;  //max multiplicity for CALM
      eEventType mEventType;
};

#endif
