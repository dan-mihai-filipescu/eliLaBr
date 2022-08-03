#include <string>
#include <iostream>
using namespace std;

char *vector_char(int imin, int imax)
{
	 char *p;
	 p = (char *) malloc((size_t) ((imax-imin+1)*sizeof(char)));
	 if (!p) {
	            printf("Vector: eroare de alocare !\n");
	            exit(1);
	         }
	 return p-imin;
}

void read2()
{
  int nchan = 5000;
  stringstream ESingle;

  TFile* in_file=new TFile("eliGammaSource.root"); 
  
  char* rootName;   rootName = vector_char(0,200);  
  char* aplName;    aplName  = vector_char(0,200);
  int  vacuum;
  int  N_sim;
  ifstream in;
  in.open("incident_fileName",ios::in);  
  in>>rootName;   printf("\n file name is : %s",rootName);
  in>>vacuum;     printf("\n vacuum boolean is %i",vacuum);
  in>>N_sim;      printf("\n simulated for %i events",N_sim);
  in.close();        

  sscanf(rootName,"%11s",aplName);  
  if (vacuum==1) ESingle<<aplName<<"_v_apl";
       else ESingle<<aplName<<"_apl";
  cout<<"\n new file name is: "<<ESingle.str().c_str()<<"\n"; 
  printf("\n");
  
  TH1F *l1 = (TH1F*)gDirectory->Get("H1trackLength");
  ofstream ofile(ESingle.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++)
    {
     double y=l1->GetBinContent(i);
     ofile<<y/N_sim<<endl;
    }
  ofile.close(); 
  
}



/*
stringstream ESingle,ESingleN,ESingleNorm,
               ESingleInc,ESingleIncN,ESingleIncNorm,
	       AnglePol,AnglePolN;
double y, norm = 1.e+6, total;
TFile* in_file=new TFile("eliGammaSource.root"); 
//in_file->cd("eliLaBr");

  TH1F *h1 = (TH1F*)gDirectory->Get("H1Sngl");
  //h1->Draw();
  ESingle<<"ESingle"<<nfile;
  ofstream ofile(ESingle.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=h1->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

  TH1F *h2 = (TH1F*)gDirectory->Get("H1SnglN");
  //h2->Draw();
  //h2->Scale(norm*100./total);
  total = h2->Integral();
  ESingleN<<"ESingleN"<<nfile;
  ofstream ofile(ESingleN.str().c_str(), ios::out);
  ESingleNorm<<"ESingleNorm"<<nfile;
  ofstream ofile1(ESingleNorm.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=h2->GetBinContent(i);
     ofile<<y<<endl;
     ofile1<<y*norm/total<<endl;}
  ofile.close();
  ofile1.close();

  TH1F *p1 = (TH1F*)gDirectory->Get("H1Inc");
  //p1->Draw();
  ESingleInc<<"ESingleInc"<<nfile;
  ofstream ofile(ESingleInc.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=p1->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

  TH1F *p2 = (TH1F*)gDirectory->Get("H1IncN");
  //p2->Draw()
  // p2->Scale(norm*100./total);
  total = p2->Integral();
  ESingleIncN<<"ESingleIncN"<<nfile;
  ofstream ofile(ESingleIncN.str().c_str(), ios::out);
  ESingleIncNorm<<"ESingleIncNorm"<<nfile;
  ofstream ofile1(ESingleIncNorm.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=p2->GetBinContent(i);
     ofile<<y<<endl;
     ofile1<<y*norm/total<<endl;}
  ofile.close();
  ofile1.close();

  TH1F *ang = (TH1F*)gDirectory->Get("Angle");
  //ang->Draw();
  AnglePol<<"PolAngle"<<nfile;
  ofstream ofile(AnglePol.str().c_str(), ios::out);
  for (Int_t i = 0; i < 3600; i++) 
    {y=ang->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

  TH1F *ang1 = (TH1F*)gDirectory->Get("AngleN");
  //ang1->Draw();
  AnglePolN<<"PolAngleN"<<nfile;
  ofstream ofile(AnglePolN.str().c_str(), ios::out);
  for (Int_t i = 0; i < 3600; i++) 
    {y=ang1->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();
*/
