#include <string>
#include <iostream>
using namespace std;

void read2_sum()
{
const int n_sim = 10;
int nfile = 100;
int nchan = 2000;
float Emin=0., Emax=20.;
int nangle = 3600.;
float anglemin=0., anglemax=360.;
stringstream ESingle,ESingleN,ESingleNorm,
               ESingleInc,ESingleIncN,ESingleIncNorm,
	       AnglePol,AnglePolN;
stringstream ss;

double y, norm = 1.e+6, total;
TFile* in_file=new TFile(); 
//in_file->cd("eliLaBr");

TH1F *h1_sum = new TH1F("h1_sum","Energy deposit",                nchan,  Emin,     Emax);
TH1F *h2_sum = new TH1F("h2_sum","Normalized Energy deposit",     nchan,  Emin,     Emax);
TH1F *p1_sum = new TH1F("p1_sum","Incident Energy",               nchan,  Emin,     Emax);
TH1F *p2_sum = new TH1F("p2_sum","Normalized Incident Energy",    nchan,  Emin,     Emax);
TH1F *a1_sum = new TH1F("a1_sum","Polarization angle",            nangle, anglemin, anglemax);
TH1F *a2_sum = new TH1F("a2_sum","Normalized Polarization Angle", nangle, anglemin, anglemax);

TH1F *h1[n_sim];
TH1F *h2[n_sim];
TH1F *p1[n_sim];
TH1F *p2[n_sim];
TH1F *ang[n_sim];
TH1F *ang1[n_sim];

int Icnt;

    for (Icnt=1; Icnt <= n_sim; Icnt++)
    {
     ss<<"./run"<<Icnt<<"/eliGammaSource.root";
     in_file->Open(ss.str().c_str());
     
    h1[Icnt] = (TH1F*)gDirectory->Get("H1Sngl");
    h1_sum->Add(h1[Icnt],1.);
    h2[Icnt] = (TH1F*)gDirectory->Get("H1SnglN");
    h2_sum->Add(h2[Icnt],1.);
    p1[Icnt] = (TH1F*)gDirectory->Get("H1Inc");
    p1_sum->Add(p1[Icnt],1.);
    p2[Icnt] = (TH1F*)gDirectory->Get("H1IncN");
    p2_sum->Add(p2[Icnt],1.);
    ang[Icnt] = (TH1F*)gDirectory->Get("Angle");
    a1_sum->Add(ang[Icnt],1.);
    ang1[Icnt] = (TH1F*)gDirectory->Get("AngleN");
    a2_sum->Add(ang1[Icnt],1.);
    ss.str("");
    in_file->Close();
    }

  //TH1F *h1 = (TH1F*)gDirectory->Get("H1Sngl");
  //h1->Draw();
  ESingle<<"ESingle"<<nfile;
  ofstream ofile(ESingle.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=h1_sum->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

  //TH1F *h2 = (TH1F*)gDirectory->Get("H1SnglN");
  //h2->Draw();
  //h2->Scale(norm*100./total);
  total = h2_sum->Integral();
  ESingleN<<"ESingleN"<<nfile;
  ofstream ofile(ESingleN.str().c_str(), ios::out);
  ESingleNorm<<"ESingleNorm"<<nfile;
  ofstream ofile1(ESingleNorm.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=h2_sum->GetBinContent(i);
     ofile<<y<<endl;
     ofile1<<y*norm/total<<endl;}
  ofile.close();
  ofile1.close();

  //TH1F *p1 = (TH1F*)gDirectory->Get("H1Inc");
  //p1->Draw();
  ESingleInc<<"ESingleInc"<<nfile;
  ofstream ofile(ESingleInc.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=p1_sum->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

  //TH1F *p2 = (TH1F*)gDirectory->Get("H1IncN");
  //p2->Draw()
  // p2->Scale(norm*100./total);
  total = p2_sum->Integral();
  ESingleIncN<<"ESingleIncN"<<nfile;
  ofstream ofile(ESingleIncN.str().c_str(), ios::out);
  ESingleIncNorm<<"ESingleIncNorm"<<nfile;
  ofstream ofile1(ESingleIncNorm.str().c_str(), ios::out);
  for (Int_t i = 0; i < nchan; i++) 
    {y=p2_sum->GetBinContent(i);
     ofile<<y<<endl;
     ofile1<<y*norm/total<<endl;}
  ofile.close();
  ofile1.close();

  //TH1F *ang = (TH1F*)gDirectory->Get("Angle");
  //ang->Draw();
  AnglePol<<"PolAngle"<<nfile;
  ofstream ofile(AnglePol.str().c_str(), ios::out);
  for (Int_t i = 0; i < nangle; i++) 
    {y=a1_sum->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

  //TH1F *ang1 = (TH1F*)gDirectory->Get("AngleN");
  //ang1->Draw();
  AnglePolN<<"PolAngleN"<<nfile;
  ofstream ofile(AnglePolN.str().c_str(), ios::out);
  for (Int_t i = 0; i < nangle; i++) 
    {y=a2_sum->GetBinContent(i);
     ofile<<y<<endl;}
  ofile.close();

}
