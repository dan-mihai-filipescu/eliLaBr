void read3()
{
   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");

    TCanvas* c=new TCanvas("Canvas","Canvas",800,800);
    TLegend legend(0.08,0.85,0.2,0.9,"Z = +3.0 m");
    char* label=new char;


TFile* in_file=new TFile("eliGammaSource.root");

TH2F *Hele_pos = (TH2F*)gDirectory->Get("ele_pos");
TH2F *Hlaser_pos = (TH2F*)gDirectory->Get("laser_pos");
TH2F *Hhoriz_xxp = (TH2F*)gDirectory->Get("horiz_xxp");
TH2F *Hvert_yyp = (TH2F*)gDirectory->Get("vert_yyp");
TH2F *Hlaser_imp = (TH2F*)gDirectory->Get("laser_imp");

    Hele_pos->SetMinimum(-0.001);
    //Hele_pos->SetMaximum(3000.);
    Hele_pos->SetTitle("Electron Beam profile");
    Hele_pos->GetXaxis()->SetTitle("Horizontal (mm)");
    Hele_pos->GetXaxis()->SetLabelSize(0.04);
    Hele_pos->GetXaxis()->CenterTitle();
    Hele_pos->GetXaxis()->SetRangeUser(-2.,2.);
    Hele_pos->GetYaxis()->SetTitle("Vertical(mm)");
    Hele_pos->GetYaxis()->SetLabelSize(0.04);
    Hele_pos->GetYaxis()->CenterTitle();
    Hele_pos->GetYaxis()->SetRangeUser(-2.,2.);

    Hlaser_pos->SetMinimum(-0.001);
    //Hlaser_pos->SetMaximum(3000.);
    //Hlaser_pos->SetTitle("Electron Beam profile");
    Hlaser_pos->GetXaxis()->SetTitle("Horizontal (mm)");
    Hlaser_pos->GetXaxis()->SetLabelSize(0.04);
    Hlaser_pos->GetXaxis()->CenterTitle();
    Hlaser_pos->GetXaxis()->SetRangeUser(-2.,2.);
    Hlaser_pos->GetYaxis()->SetTitle("Vertical(mm)");
    Hlaser_pos->GetYaxis()->SetLabelSize(0.04);
    Hlaser_pos->GetYaxis()->CenterTitle();
    Hlaser_pos->GetYaxis()->SetRangeUser(-2.,2.);

    Hhoriz_xxp->SetMinimum(-0.001);
    //Hhoriz_xxp->SetMaximum(3000.);
    Hhoriz_xxp->SetTitle("Horizontal spacephase");
    Hhoriz_xxp->GetXaxis()->SetTitle("X position (mm)");
    Hhoriz_xxp->GetXaxis()->SetLabelSize(0.04);
    Hhoriz_xxp->GetXaxis()->CenterTitle();
    Hhoriz_xxp->GetXaxis()->SetRangeUser(-2.,2.);
    Hhoriz_xxp->GetYaxis()->SetTitle("X angle (mrad)");
    Hhoriz_xxp->GetYaxis()->SetLabelSize(0.04);
    Hhoriz_xxp->GetYaxis()->CenterTitle();
    Hhoriz_xxp->GetYaxis()->SetRangeUser(-0.0006,0.0006);

    Hvert_yyp->SetMinimum(-0.001);
    //Hvert_yyp->SetMaximum(3000.);
    Hvert_yyp->SetTitle("Vertical spacephase");
    Hvert_yyp->GetXaxis()->SetTitle("Y position (mm)");
    Hvert_yyp->GetXaxis()->SetLabelSize(0.04);
    Hvert_yyp->GetXaxis()->CenterTitle();
    Hvert_yyp->GetXaxis()->SetRangeUser(-2.,2.);
    Hvert_yyp->GetYaxis()->SetTitle("Y angle (mrad)");
    Hvert_yyp->GetYaxis()->SetLabelSize(0.04);
    Hvert_yyp->GetYaxis()->CenterTitle();
    Hvert_yyp->GetYaxis()->SetRangeUser(-0.0006,0.0006);

    Hlaser_imp->SetMinimum(-0.001);
    //Hlaser_imp->SetMaximum(3000.);
    Hlaser_imp->SetTitle("Laser impulse");
    Hlaser_imp->GetXaxis()->SetTitle("Horizontal (mm)");
    Hlaser_imp->GetXaxis()->SetLabelSize(0.04);
    Hlaser_imp->GetXaxis()->CenterTitle();
    Hlaser_imp->GetXaxis()->SetNdivisions(3,2,2,kFALSE);
    Hlaser_imp->GetYaxis()->SetTitle("Vertical(mm)");
    Hlaser_imp->GetYaxis()->SetLabelSize(0.04);
    Hlaser_imp->GetYaxis()->CenterTitle();
    Hlaser_imp->GetYaxis()->SetNdivisions(3,2,2,kFALSE);

//    Hhoriz_xxp->SetMinimum(-0.001);

c->Clear();
c->Divide(2,2);
legend.DrawClone();

c->cd(1);Hele_pos->DrawClone("Contz");
c->cd(2);Hlaser_pos->DrawClone("Contz");
c->cd(3);Hhoriz_xxp->DrawClone("Contz");
c->cd(4);Hvert_yyp->DrawClone("Contz");
//c->cd(4);Hlaser_imp->DrawClone("Contz");


TCanvas* cc=new TCanvas("Canvas1","Canvas1",800,800);

TH2F *Stokes1 = (TH2F*)gDirectory->Get("stokes1");

    Stokes1->SetMaximum(150.);
    Stokes1->SetTitle("Stokes par. #1ex");
    Stokes1->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes1->GetXaxis()->SetLabelSize(0.04);
    Stokes1->GetXaxis()->CenterTitle();
    Stokes1->GetXaxis()->SetRangeUser(0.,20.);
    Stokes1->GetYaxis()->SetTitle("Stokes value");
    Stokes1->GetYaxis()->SetLabelSize(0.04);
    Stokes1->GetYaxis()->CenterTitle();
    Stokes1->GetYaxis()->SetRangeUser(-1.01,1.01);
//TH1F *Stokes1_proj = *Stokes1->ProfileX();

TH2F *Stokes2 = (TH2F*)gDirectory->Get("stokes2");

    Stokes2->SetMaximum(150.);
    Stokes2->SetTitle("Stokes par. #2ex");
    Stokes2->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes2->GetXaxis()->SetLabelSize(0.04);
    Stokes2->GetXaxis()->CenterTitle();
    Stokes2->GetXaxis()->SetRangeUser(0.,20.);
    Stokes2->GetYaxis()->SetTitle("Stokes value");
    Stokes2->GetYaxis()->SetLabelSize(0.04);
    Stokes2->GetYaxis()->CenterTitle();
    Stokes2->GetYaxis()->SetRangeUser(-1.01,1.01);
//TH1F *Stokes2_proj = *Stokes2->ProfileX();

cc->Clear();
cc->Divide(2,2);
cc->cd(1);Stokes1->DrawClone("Colx");
cc->cd(2);Stokes2->DrawClone("Colx");
cc->cd(3);Stokes1->ProfileX()->DrawClone();
cc->cd(4);Stokes2->ProfileX()->DrawClone();


//Hgamma_pos->Draw("Colx");

//TH1F *Polarization = (TH1F*)Hgamma_pos->ProjectionY();
//TCanvas* cp=new TCanvas("CanvasP","CanvasP",800,800);
//Hgamma_pos->ProfileX()->DrawClone();

TH2F *Hgamma_pos = (TH2F*)gDirectory->Get("hist2");
TCanvas* cc1=new TCanvas("Canvas2","Canvas2",800,800);

    //Hgamma_pos->SetMaximum(150.);
    Hgamma_pos->SetTitle("Beam profile");
    Hgamma_pos->GetXaxis()->SetTitle("Horizontal (mm)");
    Hgamma_pos->GetXaxis()->SetLabelSize(0.04);
    Hgamma_pos->GetXaxis()->CenterTitle();
    Hgamma_pos->GetXaxis()->SetRangeUser(-2.,2.);
    Hgamma_pos->GetYaxis()->SetTitle("Vertical (mm)");
    Hgamma_pos->GetYaxis()->SetLabelSize(0.04);
    Hgamma_pos->GetYaxis()->CenterTitle();
    Hgamma_pos->GetYaxis()->SetRangeUser(-2.,2.);
Hgamma_pos->Draw("Colz");

//TCanvas* cp=new TCanvas("CanvasP","CanvasP",800,800);
//Hgamma_pos->ProfileX()->DrawClone();

TH2F *Htarget_ET = (TH2F*)gDirectory->Get("target-ET");
TCanvas* cc2=new TCanvas("Canvas3","Canvas3",800,800);

    //Hgamma_pos->SetMaximum(150.);
    Htarget_ET->SetTitle("Target-ET");
    Htarget_ET->GetXaxis()->SetTitle("Energy (MeV)");
    Htarget_ET->GetXaxis()->SetLabelSize(0.04);
    Htarget_ET->GetXaxis()->CenterTitle();
    Htarget_ET->GetXaxis()->SetRangeUser(0.,20.);
    Htarget_ET->GetYaxis()->SetTitle("Time (ns)");
    Htarget_ET->GetYaxis()->SetLabelSize(0.04);
    Htarget_ET->GetYaxis()->CenterTitle();
    Htarget_ET->GetYaxis()->SetRangeUser(0.,500.);
Htarget_ET->Draw("Colz");

//TH1F *Hangle = (TH1F*)gDirectory->Get("Angle");
//Hangle->GetXaxis()->SetRangeUser(0.,360.);
//Hangle->Draw();
}
