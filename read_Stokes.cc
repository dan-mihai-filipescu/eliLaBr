void read_Stokes()
{
   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");

    //TCanvas* c=new TCanvas("Canvas","Canvas",800,800);
    //TLegend legend(0.08,0.85,0.2,0.9,"Z = +3.0 m");
    //char* label=new char;

const int nbin = 1000;
float xmin=0., xmax=20., ymin=-1.01, ymax=1.01;

TFile* in_file=new TFile("eliGammaSource.root");


TCanvas* cc=new TCanvas("Canvas1","Canvas1",1200,800);

TH2F *Stokes1 = (TH2F*)gDirectory->Get("stokes1");

    Stokes1->SetMaximum(600.);
    Stokes1->SetTitle("Stokes par. #1ex");
    Stokes1->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes1->GetXaxis()->SetLabelSize(0.04);
    Stokes1->GetXaxis()->CenterTitle();
    Stokes1->GetXaxis()->SetRangeUser(xmin,xmax);
    Stokes1->GetYaxis()->SetTitle("Stokes value");
    Stokes1->GetYaxis()->SetLabelSize(0.04);
    Stokes1->GetYaxis()->CenterTitle();
    Stokes1->GetYaxis()->SetRangeUser(ymin,ymax);
//TH1F *Stokes1_proj = *Stokes1->ProfileX();

TH2F *Stokes2 = (TH2F*)gDirectory->Get("stokes2");

    Stokes2->SetMaximum(600.);
    Stokes2->SetTitle("Stokes par. #2ex");
    Stokes2->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes2->GetXaxis()->SetLabelSize(0.04);
    Stokes2->GetXaxis()->CenterTitle();
    Stokes2->GetXaxis()->SetRangeUser(xmin,xmax);
    Stokes2->GetYaxis()->SetTitle("Stokes value");
    Stokes2->GetYaxis()->SetLabelSize(0.04);
    Stokes2->GetYaxis()->CenterTitle();
    Stokes2->GetYaxis()->SetRangeUser(ymin,ymax);
//TH1F *Stokes2_proj = *Stokes2->ProfileX();

TH2F *Stokes3 = (TH2F*)gDirectory->Get("stokes3");

    Stokes3->SetMaximum(600.);
    Stokes3->SetTitle("Stokes par. #3ex");
    Stokes3->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes3->GetXaxis()->SetLabelSize(0.04);
    Stokes3->GetXaxis()->CenterTitle();
    Stokes3->GetXaxis()->SetRangeUser(xmin,xmax);
    Stokes3->GetYaxis()->SetTitle("Stokes value");
    Stokes3->GetYaxis()->SetLabelSize(0.04);
    Stokes3->GetYaxis()->CenterTitle();
    Stokes3->GetYaxis()->SetRangeUser(ymin,ymax);
//TH1F *Stokes3_proj = *Stokes3->ProfileX();

TH2F *PolAngle = (TH2F*)gDirectory->Get("PolAngle");

    PolAngle->SetMaximum(1200.);
    PolAngle->SetTitle("Pol. Angle distrib.");
    PolAngle->GetXaxis()->SetTitle("Energy (MeV)");
    PolAngle->GetXaxis()->SetLabelSize(0.04);
    PolAngle->GetXaxis()->CenterTitle();
    PolAngle->GetXaxis()->SetRangeUser(xmin,xmax);
    PolAngle->GetYaxis()->SetTitle("Angle");
    PolAngle->GetYaxis()->SetLabelSize(0.04);
    PolAngle->GetYaxis()->CenterTitle();
    PolAngle->GetYaxis()->SetRangeUser(0.,360.);
//TH1F *Stokes3_proj = *Stokes3->ProfileX();

cc->Clear();
cc->Divide(3,2);
cc->cd(1);Stokes1->DrawClone("Colx");
cc->cd(2);Stokes2->DrawClone("Colx");
cc->cd(3);Stokes3->DrawClone("Colx");
cc->cd(4);Stokes1->ProfileX()->DrawClone();
cc->cd(5);Stokes2->ProfileX()->DrawClone();
cc->cd(6);Stokes3->ProfileX()->DrawClone();


//Hgamma_pos->Draw("Colx");

//TH1F *Polarization = (TH1F*)Hgamma_pos->ProjectionY();
//TCanvas* cp=new TCanvas("CanvasP","CanvasP",800,800);
//Hgamma_pos->ProfileX()->DrawClone();

TCanvas* cc1=new TCanvas("Canvas2","Canvas2",800,800);
PolAngle->DrawClone("Colz");
}
