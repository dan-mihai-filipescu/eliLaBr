#include <TROOT.h>
#include <TChain.h>
#include <TSelector.h>

#include "TApplication.h"
#include <TClassTable.h>
#include "TClass.h"

#include "vectors.h"
#include "PileItUp.h"

#include <stdio.h>
#include <string.h>

TChain *ELIChain;
PileItUp pileitup;

using namespace std;

FILE *inNames;
const char* FileInName;
char FileIn[100];
stringstream myStringStream;
std::string inputFileName;

long int nentries, firstentry;

#ifndef __CINT__
void StandaloneApplication(int argc, char** argv)
{
   if (!TClassTable::GetDict("TSTLvector"))
     {
     gSystem->Load("libTvectors.so");
     }
   FileInName = "files";
   myStringStream << FileInName;
   inputFileName = myStringStream.str();
   cout<<inputFileName<<endl;
   inNames = fopen(inputFileName.c_str(),"r");
   
   ELIChain = new TChain("eliLaBr");

   fscanf(inNames,"%li %li",&firstentry,&nentries);
   printf("FirstEntry = %li, Nentries = %li\n",firstentry,nentries);
   while(!feof(inNames))
     {
     fscanf(inNames,"%s",FileIn);
     if(feof(inNames)) break;
     printf("INPUT File: %s\n",FileIn);
     ELIChain->Add(FileIn);
     }
   fclose(inNames);

   ELIChain->Process(&pileitup,"",nentries,firstentry);
}

// This is the standard main of C++ starting a ROOT application
int main(int argc, char** argv)
{
   /*if(argc != 5)
     {
     printf(" Usage : %s <ADC_config_file> <input_data_file> <output_data_file> <run#>\n", argv[0]);
     exit(0);
     }*/

   TApplication app("Root Application", &argc, argv);
   //gROOT->SetBatch(kTRUE);
   StandaloneApplication(app.Argc(), app.Argv());
   //StandaloneApplication();

   app.Run();

   return 0;
}
#endif
