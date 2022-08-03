#include <string>
#include <iostream>
using namespace std;

void read2()
{
stringstream ESingleNorm, ESingleInc;

double y, norm = 1.e+6, total;
TFile* in_file=new TFile("eliGammaSource.root"); 

  TH1F *h1 = (TH1F*)gDirectory->Get("H1SnglN");
  total = h1->Integral();
  ESingleNorm<<"SingleN";
  ofstream ofile1(ESingleNorm.str().c_str(), ios::out);
  for (Int_t i = 0; i < 5000; i++)
    {
     y=h1->GetBinContent(i);
     ofile1<<y*norm/total<<endl;
    }
  ofile1.close();

  TH1F *p1 = (TH1F*)gDirectory->Get("H1IncN");
  ESingleInc<<"Inc";
  ofstream ofile(ESingleInc.str().c_str(), ios::out);
  for (Int_t i = 0; i < 5000; i++)
    {
     y=p1->GetBinContent(i);
     ofile<<y<<endl;
    }
  ofile.close();
}
