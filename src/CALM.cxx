#include "CALM.h"

CALM::CALM(): mRandom(0), mNames(0), mNmean(0)
{
   mRandom = new TRandom2(0);
   mNpart = 4; //particle types (pions, kaons, protons, lambdas)
   //double Nmean[] = {8.94, 1.1, 0.648, 0.19};
   double Nmean[] = {1.493, 0.183, 0.083, 0.048}; //charged particle yields per rapidity unit from 900 GeV data from http://arxiv.org/pdf/1504.00024v1.pdf (ALICE), lambdas from http://arxiv.org/pdf/1012.3257v2.pdf (ALICE)
   double RapidityInterval = 5; //rapidity <-2.5;2.5>
   double XYZ[] = {5.,5.,5.};
   int Npartkinds[]  = {3,4,4,2};
   string Names[] = {
      "pi0139plu","pi0139min","pi0135zer",
      "Ka0492plu","Ka0492min","Ka0492zer","Ka0492zrb",
      "pr0938plu","pr0938plb","ne0939zer","ne0939zrb",
      "Lm1115zer","Lm1115zrb" };
   int it=0;
   mNmean = new double[mNpart];
   mNpartkinds = new int[mNpart];
   mNames = new string*[mNpart];
   mRapidityInterval = RapidityInterval;
   for(int i=0;i<mNpart;i++)
   {
      mNmean[i]=Nmean[i];
      mNpartkinds[i]=Npartkinds[i];
      mNames[i] = new string[Npartkinds[i]];
      for(int j=0;j<Npartkinds[i];j++)
      {
         mNames[i][j]=Names[it];
         PRINT_DEBUG_2("name["<<it<<":"<<j<<","<<i<<"] = "<<Names[it]);
         it++;
      }
   }
   mXYZ = new double[3];
   for(int i=0;i<3;i++)
      mXYZ[i] = XYZ[i];
}
CALM::~CALM()
{
   delete mRandom;
}

int CALM::GenerateParticles(ParticleDB* aPartDB, int aMultBinMin, int aMultBinMax, double aEnergy, list<Particle>* aParticles, eEventType aEventType)
{
   int Nrand[mNpart]; // number of particles generated (for each kind) - from Poisson distribution
   int Npart[mNpart][aMultBinMax]; // particle to be generated
   int Nsum, Qsum, Bsum, Ssum;
   int tmpInt;
   ParticleType* tParticleType;
   //_______distributing the total number of particles for each kind and for the specific particles
   //_______GLOBAL CONSERVATION LAWS - or one minijet for minijets with local conservation
   do
   {
      Nsum = 0;
      // generating the number of particles in each kind
      for(int i=0; i<mNpart;++i)
      {
         Nrand[i]=mRandom->Poisson(mNmean[i]*mRapidityInterval*mNpartkinds[i]);
         Nsum+=Nrand[i];
      }
   }
   while(Nsum<aMultBinMin || Nsum>aMultBinMax || (Nrand[1]+Nrand[3])%2!=0 || (Nrand[2]+Nrand[3])%2!=0);
   do
   {
      Qsum = 0;
      Ssum = 0;
      Bsum = 0;
      // generating the number of specific particles within each kind
      // check of the charge, strangeness and baryon number
      for(int i=0;i<mNpart;++i)
      {
         for(int j=0;j<Nrand[i];++j)
         {
            Npart[i][j]=(int)mRandom->Uniform(mNpartkinds[i]);
            tParticleType = aPartDB->GetParticleType(mNames[i][Npart[i][j]].c_str() );
            if ( mNames[i][Npart[i][j]].find("plu")!=std::string::npos ) Qsum++;
            else if ( mNames[i][Npart[i][j]].find("min")!=std::string::npos || mNames[i][Npart[i][j]].find("plb")!=std::string::npos ) Qsum--;
            else if ( mNames[i][Npart[i][j]].find("zer")!=std::string::npos || mNames[i][Npart[i][j]].find("zrb")!=std::string::npos ) ;
            tmpInt = tParticleType->GetNumberQ()-tParticleType->GetNumberAQ()+tParticleType->GetNumberS()-tParticleType->GetNumberAS();
            if( tmpInt ==3 ) Bsum++;
            else if( tmpInt ==-3 ) Bsum--;
            tmpInt = tParticleType->GetNumberS()-tParticleType->GetNumberAS();
            if( tmpInt ==1 ) Ssum--; //  for quark s: S=-1
            else if( tmpInt ==-1 ) Ssum++;
         }
      }
   }
   while(Qsum!=0 || Ssum!=0 || Bsum!=0);
   //________rewriting the particles into one list
   for(int i=0;i<mNpart;++i)
   {
      for(int j=0;j<Nrand[i];++j)
      {
         mParticlesThisEvent.push_back(mNames[i][Npart[i][j]] );
      }
   }
   //________XYZ generating
   double XYZrand[Nsum][3];
   for(int j=0; j<Nsum;++j)
   {
      for(int i=0; i<3;++i)
      {
         XYZrand[j][i]=mRandom->Gaus(0,mXYZ[i]);
      }
   }
   //________Genbod part
   // generate total momentum for given energy
   double TotEnergy;
   int control=0;
      switch(aEventType)
      {
      case GLOBAL:
      default:
      {
         TLorentzVector en;
         TGenPhaseSpace event;
         Particle* tParticle;
         double masses[Nsum];
         for (int i=0;i<Nsum;i++)
            masses[i]=aPartDB->GetParticleType(mParticlesThisEvent[i].c_str())->GetMass();
         do
         {
            // generate total momentum
            TF1* Ptot = new TF1("Ptot","4.33538e-02*TMath::Landau(x,3.24886e+00,2.17010e+00)*exp(8.34570e-03*x)",0,100);
            TotEnergy = Ptot->GetRandom();
            delete Ptot;
            en.SetE(TotEnergy);
            control++;
         }
         while( !(event.SetDecay(en, Nsum, masses) || control>10) );
         if (control>10)
         {
            mParticlesThisEvent.clear();
            return 99;
         }
         double weight=0;
         TLorentzVector* tmp;
         double tmpweight;
         int controltmp=0;
         do
         {
            weight = event.Generate();
            if(weight != weight) weight=0;
            tmpweight = mRandom->Uniform(1.e-13);
            controltmp++;
            if(controltmp>1e6) break;
         }
         while(weight==0 || (tmpweight>weight ));
         if (controltmp>=1e6)
            return 99;
         // saving all the particles (their momenta)
         for(int i=0;i<Nsum;i++)
         {
            tmp = event.GetDecay(i);
            tParticle = new Particle(aPartDB->GetParticleType(mParticlesThisEvent[i].c_str()));
            tParticle->SetParticlePX(tmp->E() ,tmp->Px(),tmp->Py(), tmp->Pz(),
                                     0,XYZrand[i][0],XYZrand[i][1],XYZrand[i][2],
                                     weight, 0);
            aParticles->push_back(*tParticle);
            PRINT_DEBUG_2(mParticlesThisEvent[i]<<" , "<<endl);
            delete tParticle;
         }
         break;
      }
      case MINIJETS_GLOBAL:
      {
         TGenPhaseSpace event1, event0;
         Particle* tParticle;
         double weight0, weight1;
         int it=0;
         vector<double> masses[2];
         vector<string> names[2];
         do
         {
            if(masses[0].size() > 0 || masses[1].size()>0 )
            {
               masses[0].clear();
               masses[1].clear();
               names[0].clear();
               names[1].clear();
            }
            for(int i=0;i<Nsum;++i)
            {
               if (mRandom->Integer(2))
               {
                  masses[1].push_back( aPartDB->GetParticleType(mParticlesThisEvent[i].c_str() )->GetMass() );
                  names[1].push_back( mParticlesThisEvent[i].c_str() );
               }
               else
               {
                  masses[0].push_back( aPartDB->GetParticleType(mParticlesThisEvent[i].c_str() )->GetMass() );
                  names[0].push_back( mParticlesThisEvent[i].c_str() );
               }
            }
         }while( masses[0].size() < 4 || masses[1].size() < 4);
         double masses0 [masses[0].size()];
         double masses1 [masses[1].size()];
         for(int j=0;j<masses[0].size();++j) masses0[j] = masses[0][j];
         for(int j=0;j<masses[1].size();++j) masses1[j] = masses[1][j];
         TLorentzVector* tmp;
         TLorentzVector en;
         double divideEn[] = {1,1}; // 0: energy of particles, 1: boostenergy
         do
         {
            // generate total momentum
            TF1* Ptot = new TF1("Ptot","4.33538e-02*TMath::Landau(x,3.24886e+00,2.17010e+00)*exp(8.34570e-03*x)",0,100);
            TotEnergy = Ptot->GetRandom();
            delete Ptot;
            en.SetE(TotEnergy*(divideEn[0]/(2.*(divideEn[0]+divideEn[1]))));
            control++;
         }
         while( !( ((event0.SetDecay(en, masses[0].size(), masses0)) && (event1.SetDecay(en, masses[1].size(), masses1)) ) || control >10) );
         if (control>10)
         {
            mParticlesThisEvent.clear();
            return 99;
         }
         double tmpweight;
         int controltmp=0;
         do
         {
            weight0 = event0.Generate();
            weight1 = event1.Generate();
            if( (weight0 != weight0) || (weight1 != weight1) )
            {
               weight1=0;
               weight0=0;
            }
            tmpweight = mRandom->Uniform(1.e-13);
            controltmp++;
            if(controltmp>1e6) break;
         }
         while(weight0==0 || weight1==0 || (tmpweight>weight0*weight1) );
         if (controltmp>=1e6)
            return 99;
         // generate boost momentum
         double phi, eta, theta, p1[3], p2[3], Ejet1, Ejet2;
         phi = mRandom->Uniform(0,2*TMath::Pi());
         eta = mRandom->Uniform(-2.,2.);
         theta = 2*TMath::ATan(TMath::Exp(-eta));
         Ejet1 = TotEnergy*(divideEn[1]/(2.*(divideEn[0]+divideEn[1])))/masses[0].size();
         p1[0] = Ejet1 * TMath::Sin(theta) * TMath::Sin(phi) ;
         p1[1] = Ejet1 * TMath::Sin(theta) * TMath::Cos(phi) ;
         p1[2] = Ejet1 * TMath::Cos(theta) ;
         Ejet2 = TotEnergy*(divideEn[1]/(2.*(divideEn[0]+divideEn[1])))/masses[1].size();
         p2[0] = Ejet2 * TMath::Sin(theta) * TMath::Sin(phi) ;
         p2[1] = Ejet2 * TMath::Sin(theta) * TMath::Cos(phi) ;
         p2[2] = Ejet2 * TMath::Cos(theta) ;
         for(int i=0;i<masses[0].size();i++)
         {
            tmp = event0.GetDecay(i);
            tParticle = new Particle(aPartDB->GetParticleType( names[0][i] ));
            tParticle->SetParticlePX(tmp->E()+Ejet1 ,tmp->Px()+p1[0],tmp->Py()+p1[1], tmp->Pz()+p1[2],
                                     0,XYZrand[i][0],XYZrand[i][1],XYZrand[i][2],
                                     weight0*weight1, 0);
            aParticles->push_back(*tParticle);
            delete tParticle;
         }
         for(int i=0;i<masses[1].size();i++)
         {
            tmp = event1.GetDecay(i);
            tParticle = new Particle(aPartDB->GetParticleType( names[1][i] ));
            tParticle->SetParticlePX(tmp->E()+Ejet2 ,tmp->Px()-p2[0],tmp->Py()-p2[1], tmp->Pz()-p2[2],
                                     0,XYZrand[masses[0].size()+i][0],XYZrand[masses[0].size()+i][1],XYZrand[masses[0].size()+i][2],
                                     weight0*weight1, 0);
            aParticles->push_back(*tParticle);
            delete tParticle;
         }
         break;
      }
      case MINIJETS_LOCAL:
      {
         TGenPhaseSpace event1, event0;
         Particle* tParticle;
         double weight0, weight1;
         int it=0;
         vector<double> masses[2];
         vector<string> names[2];
         int Qjet[2],Bjet[2],Sjet[2];
         do
         {
            if(masses[0].size() > 0 || masses[1].size() > 0 )
            {
               masses[0].clear();
               masses[1].clear();
               names[0].clear();
               names[1].clear();
            }
            for(int it_clean=0;it_clean<3;it_clean++)
            {
               Qjet[it_clean] = 0;
               Sjet[it_clean] = 0;
               Bjet[it_clean] = 0;
            }
            for(int i=0;i<Nsum;++i)
            {
               if (mRandom->Integer(2))
               {
                  tParticleType = aPartDB->GetParticleType( mParticlesThisEvent[i].c_str() );
                  masses[1].push_back( tParticleType->GetMass() );
                  names[1].push_back( mParticlesThisEvent[i].c_str() );
                  if ( mParticlesThisEvent[i].find("plu")!=std::string::npos ) Qjet[1]++;
                  else if ( mParticlesThisEvent[i].find("min")!=std::string::npos || mParticlesThisEvent[i].find("plb")!=std::string::npos ) Qjet[1]--;
                  else if ( mParticlesThisEvent[i].find("zer")!=std::string::npos || mParticlesThisEvent[i].find("zrb")!=std::string::npos ) ;
                  tmpInt = tParticleType->GetNumberQ()-tParticleType->GetNumberAQ()+tParticleType->GetNumberS()-tParticleType->GetNumberAS();
                  if( tmpInt ==3 ) Bjet[1]++;
                  else if( tmpInt ==-3 ) Bjet[1]--;
                  tmpInt = tParticleType->GetNumberS()-tParticleType->GetNumberAS();
                  if( tmpInt ==1 ) Sjet[1]--; //  for quark s: S=-1
                  else if( tmpInt ==-1 ) Sjet[1]++;
               }
               else
               {
                  tParticleType = aPartDB->GetParticleType( mParticlesThisEvent[i].c_str() );
                  masses[0].push_back( tParticleType->GetMass() );
                  names[0].push_back( mParticlesThisEvent[i].c_str() );
                  if ( mParticlesThisEvent[i].find("plu")!=std::string::npos ) Qjet[0]++;
                  else if ( mParticlesThisEvent[i].find("min")!=std::string::npos || mParticlesThisEvent[i].find("plb")!=std::string::npos ) Qjet[0]--;
                  else if ( mParticlesThisEvent[i].find("zer")!=std::string::npos || mParticlesThisEvent[i].find("zrb")!=std::string::npos ) ;
                  tmpInt = tParticleType->GetNumberQ()-tParticleType->GetNumberAQ()+tParticleType->GetNumberS()-tParticleType->GetNumberAS();
                  if( tmpInt ==3 ) Bjet[0]++;
                  else if( tmpInt ==-3 ) Bjet[0]--;
                  tmpInt = tParticleType->GetNumberS()-tParticleType->GetNumberAS();
                  if( tmpInt ==1 ) Sjet[0]--; //  for quark s: S=-1
                  else if( tmpInt ==-1 ) Sjet[0]++;
               }
            }
            control++;
            if(control>100) break;
         }while( Qjet[0]!=0 || Sjet[0]!=0 || Bjet[0]!=0 || Qjet[1]!=0 || Sjet[1]!=0 || Bjet[1]!=0  );
         if (control>100)
         {
            mParticlesThisEvent.clear();
            return 99;
         }
         else
            control =0;
         double masses0 [masses[0].size()];
         double masses1 [masses[1].size()];
         for(int j=0;j<masses[0].size();++j) masses0[j] = masses[0][j];
         for(int j=0;j<masses[1].size();++j) masses1[j] = masses[1][j];
         TLorentzVector* tmp;
         TLorentzVector en;
         double divideEn[] = {1,1}; // 0: energy of particles, 1: boostenergy
         do
         {
            // generate total momentum
            TF1* Ptot = new TF1("Ptot","4.33538e-02*TMath::Landau(x,3.24886e+00,2.17010e+00)*exp(8.34570e-03*x)",0,100);
            TotEnergy = Ptot->GetRandom();
            delete Ptot;
            en.SetE(TotEnergy*(divideEn[0]/(2.*(divideEn[0]+divideEn[1]))));
            control++;
         }
         while( !( ((event0.SetDecay(en, masses[0].size(), masses0)) && (event1.SetDecay(en, masses[1].size(), masses1)) ) || control >10) );
         if (control>10)
         {
            mParticlesThisEvent.clear();
            return 99;
         }
         double tmpweight;
         int controltmp=0;
         do
         {
            weight0 = event0.Generate();
            weight1 = event1.Generate();
            if( (weight0 != weight0) || (weight1 != weight1) )
            {
               weight1=0;
               weight0=0;
            }
            tmpweight = mRandom->Uniform(1.e-13);
            controltmp++;
            if(controltmp>1e6) break;
         }
         while(weight0==0 || weight1==0 || (tmpweight>weight0*weight1) );
         if (controltmp>=1e6)
            return 99;
         // generate boost momentum
         double phi, eta, theta, p1[3], p2[3], Ejet1, Ejet2;
         phi = mRandom->Uniform(0,2*TMath::Pi());
         eta = mRandom->Uniform(-2.,2.);
         theta = 2*TMath::ATan(TMath::Exp(-eta));
         Ejet1 = TotEnergy*(divideEn[1]/(2.*(divideEn[0]+divideEn[1])))/masses[0].size();
         p1[0] = TotEnergy/4./masses[0].size() * TMath::Sin(theta) * TMath::Sin(phi) ;
         p1[1] = TotEnergy/4./masses[0].size() * TMath::Sin(theta) * TMath::Cos(phi) ;
         p1[2] = TotEnergy/4./masses[0].size() * TMath::Cos(theta) ;
         Ejet2 = TotEnergy*(divideEn[1]/(2.*(divideEn[0]+divideEn[1])))/masses[1].size();
         p2[0] = TotEnergy/4./masses[1].size() * TMath::Sin(theta) * TMath::Sin(phi) ;
         p2[1] = TotEnergy/4./masses[1].size() * TMath::Sin(theta) * TMath::Cos(phi) ;
         p2[2] = TotEnergy/4./masses[1].size() * TMath::Cos(theta) ;
         for(int i=0;i<masses[0].size();i++)
         {
            tmp = event0.GetDecay(i);
            tParticle = new Particle(aPartDB->GetParticleType( names[0][i] ));
            tParticle->SetParticlePX(tmp->E()+Ejet1 ,tmp->Px()+p1[0],tmp->Py()+p1[1], tmp->Pz()+p1[2],
                                     0,XYZrand[i][0],XYZrand[i][1],XYZrand[i][2],
                                     weight0*weight1, 0);
            aParticles->push_back(*tParticle);
            delete tParticle;
         }
         for(int i=0;i<masses[1].size();i++)
         {
            tmp = event1.GetDecay(i);
            tParticle = new Particle(aPartDB->GetParticleType( names[1][i] ));
            tParticle->SetParticlePX(tmp->E()+Ejet2 ,tmp->Px()-p2[0],tmp->Py()-p2[1], tmp->Pz()-p2[2],
                                     0,XYZrand[masses[0].size()+i][0],XYZrand[masses[0].size()+i][1],XYZrand[masses[0].size()+i][2],
                                     weight0*weight1, 0);
            aParticles->push_back(*tParticle);
            delete tParticle;
         }
         break;
      }
      }
   mParticlesThisEvent.clear();
   return 0;
}
