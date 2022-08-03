#include <string>
#include <iostream>
using namespace std;

void read_BeamSpot()
{
stringstream ProjeX, ProjeY;
double x,y;

   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");

const int nbin = 300;
float xmin=-2.5, xmax=2.5, Emin=4., Emax=20.;

    TCanvas* c=new TCanvas("Canvas","Canvas",800,800);
    TLegend legend(0.08,0.85,0.2,0.9,"Z = +0.0 m");
    char* label=new char;


TFile* in_file=new TFile("eliGammaSource.root");

TH2F *Hgamma_pos = (TH2F*)gDirectory->Get("GammaBeam-profile");
TH2F *En_radius =  (TH2F*)gDirectory->Get("En-radius");
TH2F *En_X =       (TH2F*)gDirectory->Get("En-X");
TH2F *En_Y =       (TH2F*)gDirectory->Get("En-Y");

    Hgamma_pos->SetMinimum(0.00);
    Hgamma_pos->SetTitle("Beam profile");
    Hgamma_pos->GetXaxis()->SetTitle("Horizontal (mm)");
    Hgamma_pos->GetXaxis()->SetLabelSize(0.04);
    Hgamma_pos->GetXaxis()->CenterTitle();
    //Hgamma_pos->GetXaxis()->SetRangeUser(xmin,xmax);
    Hgamma_pos->GetYaxis()->SetTitle("Vertical (mm)");
    Hgamma_pos->GetYaxis()->SetLabelSize(0.04);
    Hgamma_pos->GetYaxis()->CenterTitle();
    //Hgamma_pos->GetYaxis()->SetRangeUser(xmin,xmax);

    En_radius->SetMinimum(0.00);
    En_radius->SetTitle("Beam profile");
    En_radius->GetXaxis()->SetTitle("Beam radius (mm)");
    En_radius->GetXaxis()->SetLabelSize(0.04);
    En_radius->GetXaxis()->CenterTitle();
    //En_radius->GetXaxis()->SetRangeUser(0.,xmax);
    En_radius->GetYaxis()->SetTitle("Energy (MeV)");
    En_radius->GetYaxis()->SetLabelSize(0.04);
    En_radius->GetYaxis()->CenterTitle();
    //En_radius->GetYaxis()->SetRangeUser(Emin,Emax);

    En_X->SetMinimum(0.00);
    En_X->SetTitle("Beam profile");
    En_X->GetXaxis()->SetTitle("Beam X (mm)");
    En_X->GetXaxis()->SetLabelSize(0.04);
    En_X->GetXaxis()->CenterTitle();
    //En_X->GetXaxis()->SetRangeUser(xmin,xmax);
    En_X->GetYaxis()->SetTitle("Energy (MeV)");
    En_X->GetYaxis()->SetLabelSize(0.04);
    En_X->GetYaxis()->CenterTitle();
    //En_X->GetYaxis()->SetRangeUser(Emin,Emax);

    En_Y->SetMinimum(0.00);
    En_Y->SetTitle("Beam profile");
    En_Y->GetXaxis()->SetTitle("Beam Y (mm)");
    En_Y->GetXaxis()->SetLabelSize(0.04);
    En_Y->GetXaxis()->CenterTitle();
    //En_Y->GetXaxis()->SetRangeUser(xmin,xmax);
    En_Y->GetYaxis()->SetTitle("Energy (MeV)");
    En_Y->GetYaxis()->SetLabelSize(0.04);
    En_Y->GetYaxis()->CenterTitle();
    //En_Y->GetYaxis()->SetRangeUser(Emin,Emax);

c->Clear();
c->Divide(2,2);
legend.DrawClone();

c->cd(1);Hgamma_pos->DrawClone("Colz");
c->cd(2);En_radius->DrawClone("Colz");
c->cd(3);En_X->DrawClone("Colz");
c->cd(4);En_Y->DrawClone("Colz");
//delete in_file;

TH1F* averageEnergy = new TH1F ( "Average Energy", "Gamma beam energy distribution on target surface - Average Energy",nbin, 0., xmax );
TH1F* projection_on_r = new TH1F ( "Projection", "Gamma beam energy distribution on target surface -  Projection",nbin, 0., xmax );

averageEnergy->Add(En_radius->ProfileX(),1.);
averageEnergy->GetXaxis()->SetTitle("Beam radius (mm)");
averageEnergy->GetXaxis()->SetLabelSize(0.04);
averageEnergy->GetXaxis()->CenterTitle();
averageEnergy->GetYaxis()->SetRangeUser(Emin,Emax);
projection_on_r->Add(En_radius->ProjectionX(),1.);
projection_on_r->GetXaxis()->SetTitle("Beam radius (mm)");
projection_on_r->GetXaxis()->SetLabelSize(0.04);
projection_on_r->GetXaxis()->CenterTitle();

TCanvas* projections=new TCanvas("Canvas1","Canvas1",800,400);
projections->Clear();
projections->Divide(2);
projections->cd(1);averageEnergy->DrawClone();
projections->cd(2);projection_on_r->DrawClone();

TH1F* projection_on_X = new TH1F ( "ProjectionX", "Gamma beam distribution on target surface -  ProjectionX",nbin, -5., 5. );
TH1F* projection_on_Y = new TH1F ( "ProjectionY", "Gamma beam distribution on target surface -  ProjectionY",nbin, -5., 5. );

projection_on_X->Add(Hgamma_pos->ProjectionX(),1.);
projection_on_X->GetXaxis()->SetTitle("Beam profile X (mm)");
projection_on_X->GetXaxis()->SetLabelSize(0.04);
projection_on_X->GetXaxis()->CenterTitle();

  ProjeX<<"ProjeX.txt";
  ofstream ofile1(ProjeX.str().c_str(), ios::out);
  for (Int_t i = 0; i < nbin; i++)
    {
     x=projection_on_X->GetBinCenter(i);
     y=projection_on_X->GetBinContent(i);
     ofile1<<x<<'\t'<<y<<endl;
    }
  ofile1.close();

projection_on_Y->Add(Hgamma_pos->ProjectionY(),1.);
projection_on_Y->GetXaxis()->SetTitle("Beam profile Y (mm)");
projection_on_Y->GetXaxis()->SetLabelSize(0.04);
projection_on_Y->GetXaxis()->CenterTitle();

  ProjeY<<"ProjeY.txt";
  ofstream ofile2(ProjeY.str().c_str(), ios::out);
  for (Int_t i = 0; i < nbin; i++)
    {
     x=projection_on_X->GetBinCenter(i);
     y=projection_on_Y->GetBinContent(i);
     ofile2<<x<<'\t'<<y<<endl;
    }
  ofile2.close();


TCanvas* projections1=new TCanvas("Canvas2","Canvas2",800,400);
projections1->Clear();
projections1->Divide(2);
projections1->cd(1);projection_on_X->DrawClone();
projections1->cd(2);projection_on_Y->DrawClone();
}
