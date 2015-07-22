/********************************************************************************
 *                                                                              *
 *                      CALM: ConservAtion Laws Model                           *
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

#include <fstream>
#include <iostream>
#include <TString.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "EventGenerator.h"

Configurator *sMainConfig;
TString	sMainINI;
TString sHyperXML;
TString	sEventDIR;
TString	sTimeStamp;
int	sRandomize;
int	sIntegrateSample;
int	sParentPID;

void ReadParameters();
void ReadSHARE(ParticleDB* aPartDB);
void CheckSHARE(ParticleDB* aPartDB);
void MessageIntro();
void MessageHelp();
void MessageVersion();
void CopyINIFile();
void AddLogEntry(const char* aEntry);

int main(int argc, char **argv)
{
  ParticleDB*	  tPartDB; 
  EventGenerator* tEventGen;

  sMainINI   = "./events.ini";
  sParentPID = 0;
  sHyperXML  = "";

  if (argc > 1) {
    TString tDummy;
    for(int i=1; i<argc;i++) {
      tDummy = argv[i];
      if((tDummy == "-h") || (tDummy == "--help")) {
        MessageHelp();
        return 0;
      } else if((tDummy == "-v") || (tDummy == "--version")) {
        MessageVersion();
        return 0;
      } else if (tDummy.EndsWith(".xml"))
        sHyperXML = tDummy;
      else if (tDummy.EndsWith(".ini"))
	sMainINI  = tDummy;
      else if (tDummy.IsDigit())
	sParentPID = tDummy.Atoi();
    }
  }

  MessageIntro();
  sMainConfig = new Configurator;
  ReadParameters();
  
  {
    char tBuff[2*kFileNameMaxChar];
    std::sprintf(tBuff,"[input]\t%s\t%i",sMainINI.Data(),sParentPID); 
    AddLogEntry(tBuff);
  }

  tPartDB     = new ParticleDB();
  ReadSHARE(tPartDB);

  
  tEventGen   = new EventGenerator(tPartDB); 
  tEventGen->GenerateEvents();
  tEventGen->SetEventsTemp();
 
  delete tEventGen;
  delete tPartDB;
  delete sMainConfig;

  return 0;
}

// ##############################################################
// #			---===[ END OF MAIN ]===--		#
// ##############################################################

void ReadParameters()
{
  TDatime tDate;
  Parser* tParser;
  
  tDate.Set();
  sTimeStamp = tDate.AsSQLString();
  
  tParser = new Parser(sMainINI);
  tParser->ReadINI(sMainConfig);
  delete tParser;
  try {
    sRandomize = sMainConfig->GetParameter("Randomize").Atoi();
  } catch (TString tError) {
    PRINT_MESSAGE("<calm_events::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
}

void ReadSHARE(ParticleDB* aPartDB)
{
  TString tShareDir;
  Parser* tParser;
  
  try {
    tShareDir = sMainConfig->GetParameter("ShareDir"); tShareDir.Prepend("./");
  } catch (TString tError) {
    PRINT_DEBUG_1("<Parser::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find SHARE input file location.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }

  tParser = new Parser((tShareDir + "particles.data").Data());
  tParser->ReadSHAREParticles(aPartDB);
  delete tParser;
  
  tParser = new Parser((tShareDir + "decays.data").Data());
  tParser->ReadSHAREDecays(aPartDB);
  delete tParser;
}

void MessageIntro()
{
  PRINT_MESSAGE("  ***********************************************************************"	);
  PRINT_MESSAGE("  *\t\tCALM 2 version "<<_CALM_VERSION_<<"\t\t\t*"	);
  PRINT_MESSAGE("  *\t\t\t\t\t\t\t\t\t*"							);
  PRINT_MESSAGE("  * authors: M.Chojnacki, A.Kisiel, W.Florkowski, W.Broniowski\t\t*"		);
  PRINT_MESSAGE("  * cite as: arXiv:1102.0273\t\t\t\t\t\t*"					);
  PRINT_MESSAGE("  * webpage: http://therminator2.ifj.edu.pl/\t\t\t\t*"				);
  PRINT_MESSAGE("  ***********************************************************************"	);
}

void MessageHelp()
{
  PRINT_MESSAGE("Usage:");
  PRINT_MESSAGE("calm_events [EVENTS_INI] [PPID] [HYPER_XML]");
  PRINT_MESSAGE("calm2_events [OPTION]");
  PRINT_MESSAGE("  [EVENTS_INI]\t\tmain settings file;\t\tdefault: events.ini");
  PRINT_MESSAGE("  [PPID]\t\tparent's system process ID;\tdefault: 0");
  PRINT_MESSAGE("  [HYPER_XML]\tlocation of the hypersurface XML file;\tdefault:");
  PRINT_MESSAGE("  [OPTION]");
  PRINT_MESSAGE("    -h | --help\t\tthis screen");
  PRINT_MESSAGE("    -v | --version\tversion information");
}

void MessageVersion()
{
  PRINT_MESSAGE("version:\tCALM 2 version "<<_CALM_VERSION_);
  PRINT_MESSAGE("compiled with:\t"<<_CXX_VER_<<", ROOT("<<_ROOT_VER_<<")");
  std::cout <<  "  preprocessor: ";
#ifdef _DEBUG_LEVEL_
  std::cout << "DEBUG="<<_DEBUG_LEVEL_<<" ";
#endif
#ifdef _PARTICLE_DECAYER_RESCALE_CHANNELS_
  std::cout << "RESCALE_CHANNELS="<<_PARTICLE_DECAYER_RESCALE_CHANNELS_<<" ";
#endif
#ifdef _PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_
  std::cout << "DISABLE_THREE_BODY_DECAYS="<<_PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_<<" ";
#endif
#ifdef _PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_
  std::cout << "DISABLE_TWO_BODY_DECAYS="<<_PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_<<" ";
#endif
  std::cout << std::endl;
}



void AddLogEntry(const char* aEntry)
{
  TString tLogName;
  TDatime tDate;
  ofstream tFile;
  
  tDate.Set();  
  try {
    tLogName = sMainConfig->GetParameter("LogFile"); tLogName.Prepend("./");
  }
  catch (TString tError) {
    return;
  }
  
  tFile.open(tLogName, std::ios_base::app);
  if (static_cast<long>(tFile.tellp()) == 0) {
    tFile << "# CALM Log File"<<std::endl;
  }
  tFile << '['<<tDate.AsSQLString()<<"]\tcalm_events\t"<<sParentPID<<'\t';
  tFile << aEntry << std::endl;
  tFile.close();
}


