/********************************************************************************
 *																																							*
 *                 THERMINATOR 2: THERMal heavy-IoN generATOR 2                               *
 *																																							*
 * Version:                                                                                                     *
 *			Release, 2.0.3, 1 February 2011                                                            *
 *																																							*
 * Authors:                                                                                                     *
 *			Mikolaj Chojnacki  (Mikolaj.Chojnacki@ifj.edu.pl)											*
 *			Adam Kisiel           (kisiel@if.pw.edu.pl)                                             *
 *			Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)										*
 *			Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)										*
 *																																							*
 * Project homepage:																														*
 *			http://therminator2.ifj.edu.pl/                                                            *
 *																																							*
 * For the detailed description of the program and further references                *
 * to the description of the model please refer to															*
 * http://arxiv.org/abs/1102.0273                                                                      *
 *																																							*
 * This code can be freely used and redistributed. However if you decide to       *
 * make modifications to the code, please, inform the authors.									*
 * Any publication of results obtained using this code must include the           *
 * reference to arXiv:1102.0273 and the published version of it, when                *
 * available.                                                                                                   *
 *																																							*
 ********************************************************************************/

#include <sstream>
#include <TDatime.h>
#include "Crc32.h"
#include "Configurator.h"
#include "Event.h"
#include "THGlobal.h"

extern Configurator* sMainConfig;
extern TString	sTimeStamp;
extern int	sRandomize;

using namespace std;

Event::Event()
   : mPartDB(0), mCALM(0), mRandom(0), mMultMin(10), mMultMax(20), mEventType(GLOBAL)
{
   Reset();
}

Event::Event(ParticleDB* aDB, CALM* aCALM)
: mPartDB(aDB), mCALM(aCALM), mMultMin(10), mMultMax(20)
{
   mRandom = new TRandom2();
   Reset();
   ReadParameters();
}

Event::~Event()
{
   mParticles.clear();
   delete mRandom;
}

void Event::Reset(int aEventIter)
{
   ostringstream oss;
   Crc32 tEventID;

   mParticles.clear();
   Particle::ZeroEID();

   oss << sTimeStamp.Data() << "Event: " << aEventIter;
   tEventID.Update(oss.str().data(), oss.str().length());
   tEventID.Finish();
   mEventID = tEventID.GetValue();
}

list<Particle>* Event::GetParticleList()
{
   return &mParticles;
}

ParticleDB* Event::GetParticleDB() const
{
   return mPartDB;
}

unsigned int Event::GetEventID() const
{
   return mEventID;
}

void Event::GeneratePrimordials(int aSeed)
{
   int control=99;
   do
     {
       control = mCALM->GenerateParticles(mPartDB, mMultMin, mMultMax, &mParticles, mEventType);
     }
   while( control == 99 );

}

void Event::Randomize()
{
   TDatime tDate;

   mRandom->SetSeed(tDate.Get() / 2 * 3);
}

void Event::ReadParameters()
{
   std::string tMultMin, tMultMax, tEventType;
   int tMultMinInt, tMultMaxInt;
   int tEventTypeEnum;
   std::stringstream tConvert;
   try {
      tMultMin	= sMainConfig->GetParameter("MultiplicityMin");
      tMultMax	= sMainConfig->GetParameter("MultiplicityMax");
      tEventType	= sMainConfig->GetParameter("EventType");
      tConvert<<tMultMin<<' '<<tMultMax<<' '<<' '<<tEventType;
      tConvert>>tMultMinInt>>tMultMaxInt>>tEventTypeEnum;
      if(tMultMinInt<tMultMaxInt)
      {
         mMultMin = tMultMinInt;
         mMultMax = tMultMaxInt;
      }
      mEventType = (eEventType) tEventTypeEnum;
      PRINT_DEBUG_1("<Event::ReadParameters>\tSetting multiplicity range ("<<mMultMin<<","<<mMultMax);
   }
   catch (TString tError) {
      PRINT_DEBUG_1("<Event::ReadParameters>\tError reading parameters. Using default ones.");
   }
}
