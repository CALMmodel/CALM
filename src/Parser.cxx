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

#include <fstream>
#include <sstream>
#include <math.h>
#include <TMath.h>
#include "THGlobal.h"
#include "Parser.h"

using namespace std;

const double factorials[7] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0 };

Parser::Parser()
: mFileName(""), mFile(0)
{
}

Parser::Parser(const char* aFileName)
: mFileName(aFileName)
{
  mFile.open(mFileName);
  if((!mFile) && (!mFile.is_open())) {
    PRINT_MESSAGE("<Parser::Parser>\tFile "<<mFileName<<" not found.");
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
  PRINT_DEBUG_2("<Parser::Parser>\tReading form file "<<mFileName);
}

Parser::~Parser()
{
  mFile.close();
}

void Parser::ReadINI(Configurator* aINI)
{
  Parameter  	tOpt;
  TString	tStr;
  int		tLine;
  int		i;
  char    	tBuff[200];

  PRINT_DEBUG_3("<Parser::ReadINI>");
  mFile.seekg(0, std::ios::beg);
  tLine = 0;
  while ((!mFile.eof())) {
    tLine++;
    mFile.getline(tBuff,200);
    tStr = tBuff;
    if((tStr.IsNull()) || (tStr.IsWhitespace()) || (tStr[0] == '#') || (tStr[0] == ';')) {
      PRINT_DEBUG_3("\tLine "<<tLine<<" Ignoring  : "<<tStr.Data());
      continue;
    } else if(tStr.Contains('[')) {
      PRINT_DEBUG_2("\tLine "<<tLine<<" Section   : "<<tStr.Data());
      continue;
    } else if(!tStr.Contains('=')) {
      PRINT_DEBUG_1("\tLine "<<tLine<<" WARNING   : "<<tStr.Data());
      PRINT_DEBUG_1("\tWrong format. Treating as commentary.")
      continue;
    }
    tStr.ReplaceAll(" ","");
    tStr.ReplaceAll("\t","");
    tOpt.keyword = "";
    tOpt.value = "";
    for(i=0; tStr[i] != '='; i++)
      tOpt.keyword += tStr[i];
    for(i++; i<tStr.Length(); i++)
      tOpt.value += tStr[i];
    aINI->AddParameter(&tOpt);
    PRINT_DEBUG_2("\tLine "<<tLine<<" Parameter : "<<tOpt.keyword<<" = "<<tOpt.value);
  }
}

void Parser::ReadSHAREParticles(ParticleDB* aDB)
{
  istringstream* iss;
  ParticleType*	 tPartBuf;
  char   buff[200];
  char   name[20];
  double mass, gamma, spin, I, I3, Nq, Ns, Naq, Nas, Nc, Nac, MC;
  int    number = 0;   

  while (!mFile.eof()) {
    mFile.getline(buff,200);
    if (!(*buff) || (*buff == '#'))
      continue;
    iss = new istringstream(buff);   
    (*iss) >> name >> mass >> gamma >> spin >> I >> I3 >> Nq >> Ns >> Naq >> Nas >> Nc >> Nac >> MC;
    number++;
    PRINT_DEBUG_2('\t'<<number<<" "<<name<<" "<<mass<<" "<<gamma<<" "<<spin<<" "<<I<<" "<<I3<<" "<<Nq<<" "<<Naq<<" "<<Ns<<" "<<Nas<<" "<<Nc<<" "<<Nac<<" "<<MC);
    tPartBuf = new ParticleType();
    tPartBuf->SetNumber(number);
    tPartBuf->SetName(name);
    tPartBuf->SetMass(mass);
    tPartBuf->SetGamma(gamma);
    tPartBuf->SetSpin(spin);
    tPartBuf->SetBarionN(static_cast<int> ((Nq + Ns + Nc)/3. - (Naq + Nas + Nac)/3.) );
    tPartBuf->SetI(I);
    tPartBuf->SetI3(I3);
    tPartBuf->SetStrangeN(static_cast<int> (Nas - Ns));
    tPartBuf->SetCharmN(static_cast<int> (Nc - Nac));
    tPartBuf->SetNumberQ(static_cast<int> (Nq));
    tPartBuf->SetNumberAQ(static_cast<int> (Naq));
    tPartBuf->SetNumberS(static_cast<int> (Ns));
    tPartBuf->SetNumberAS(static_cast<int> (Nas));
    tPartBuf->SetNumberC(static_cast<int> (Nc));
    tPartBuf->SetNumberAC(static_cast<int> (Nac));
    tPartBuf->SetPDGCode(static_cast<int> (MC));
    aDB->AddParticleType(tPartBuf);
    delete iss;
  }
}


double Parser::SHAREClebschGordan(double aJot,  double aEm,
				  double aJot1, double aEm1,
				  double aJot2, double aEm2)
{
  int mint, maxt;
  double cgc = 0.0;
  int tIter;
  double coef;

  maxt = lrint(aJot1 + aJot2 - aJot);
  mint = 0;
  if (lrint(aJot1 - aEm1)	< maxt) maxt = lrint(aJot1 - aEm1);
  if (lrint(aJot2 + aEm2)	< maxt) maxt = lrint(aJot2 + aEm2);
  if (lrint(-(aJot-aJot2+aEm1))	> mint) mint = lrint(-(aJot-aJot2+aEm1));
  if (lrint(-(aJot-aJot1-aEm2))	> mint) mint = lrint(-(aJot-aJot1-aEm2));
  
  PRINT_DEBUG_3("\t\t\tmint " << mint << " j1 " <<  aJot1 << " m1 " << aEm1);
  PRINT_DEBUG_3("\t\t\tmaxt " << maxt << " j2 " <<  aJot2 << " m2 " << aEm2);

  for (tIter = mint; tIter<=maxt; tIter ++) {
    coef = TMath::Power(-1, tIter);
    PRINT_DEBUG_3("\t\t\tcoef1 " << coef);
    coef *= TMath::Sqrt( (2*aJot+1)
            * factorials[lrint(aJot1+aEm1)]
	    * factorials[lrint(aJot1-aEm1)]
	    * factorials[lrint(aJot2+aEm2)]
	    * factorials[lrint(aJot2-aEm2)]
	    * factorials[lrint(aJot+aEm)]
	    * factorials[lrint(aJot-aEm)]
	    );
    PRINT_DEBUG_3("\t\t\tcoef2 " << coef);
    coef /= factorials[tIter]
            * factorials[lrint(aJot1+aJot2-aJot-tIter)] 
	    * factorials[lrint(aJot1-aEm1-tIter)] 
	    * factorials[lrint(aJot2+aEm2-tIter)] 
	    * factorials[lrint(aJot-aJot2+aEm1+tIter)] 
	    * factorials[lrint(aJot-aJot1-aEm2+tIter)];
    PRINT_DEBUG_3("\t\t\tcoef3 " << coef);
    cgc += coef;
  }
  cgc *= SHAREDeltaJ(aJot1, aJot2, aJot);

  return cgc;
}

double Parser::SHAREDeltaJ(double aJot1, double aJot2, double aJot)
{
  double res = TMath::Sqrt( 1.0
               * factorials[lrint( aJot1 + aJot2 - aJot)]
               * factorials[lrint( aJot1 - aJot2 + aJot)]
               * factorials[lrint(-aJot1 + aJot2 + aJot)]
               / factorials[lrint( aJot1 + aJot2 + aJot + 1)]
	       );
  return res;
}

