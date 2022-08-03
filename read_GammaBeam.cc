void read_GammaBeam()
{
   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");
   
float xmin=-2., xmax=2., ymin=-0.0006., ymax=0.0006;

    TCanvas* cg=new TCanvas("Canvas","Canvas",800,800);
    TLegend legendg(0.08,0.85,0.2,0.9,"Z = +3.0 m");
    char* labelg=new char;


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
    Hele_pos->GetXaxis()->SetRangeUser(xmin,xmax);
    Hele_pos->GetYaxis()->SetTitle("Vertical(mm)");
    Hele_pos->GetYaxis()->SetLabelSize(0.04);
    Hele_pos->GetYaxis()->CenterTitle();
    Hele_pos->GetYaxis()->SetRangeUser(xmin,xmax);

    Hlaser_pos->SetMinimum(-0.001);
    //Hlaser_pos->SetMaximum(3000.);
    //Hlaser_pos->SetTitle("Electron Beam profile");
    Hlaser_pos->GetXaxis()->SetTitle("Horizontal (mm)");
    Hlaser_pos->GetXaxis()->SetLabelSize(0.04);
    Hlaser_pos->GetXaxis()->CenterTitle();
    Hlaser_pos->GetXaxis()->SetRangeUser(xmin,xmax);
    Hlaser_pos->GetYaxis()->SetTitle("Vertical(mm)");
    Hlaser_pos->GetYaxis()->SetLabelSize(0.04);
    Hlaser_pos->GetYaxis()->CenterTitle();
    Hlaser_pos->GetYaxis()->SetRangeUser(xmin,xmax);

    Hhoriz_xxp->SetMinimum(-0.001);
    //Hhoriz_xxp->SetMaximum(3000.);
    Hhoriz_xxp->SetTitle("Horizontal spacephase");
    Hhoriz_xxp->GetXaxis()->SetTitle("X position (mm)");
    Hhoriz_xxp->GetXaxis()->SetLabelSize(0.04);
    Hhoriz_xxp->GetXaxis()->CenterTitle();
    Hhoriz_xxp->GetXaxis()->SetRangeUser(xmin,xmax);
    Hhoriz_xxp->GetYaxis()->SetTitle("X angle (mrad)");
    Hhoriz_xxp->GetYaxis()->SetLabelSize(0.04);
    Hhoriz_xxp->GetYaxis()->CenterTitle();
    Hhoriz_xxp->GetYaxis()->SetRangeUser(ymin,ymax);

    Hvert_yyp->SetMinimum(-0.001);
    //Hvert_yyp->SetMaximum(3000.);
    Hvert_yyp->SetTitle("Vertical spacephase");
    Hvert_yyp->GetXaxis()->SetTitle("Y position (mm)");
    Hvert_yyp->GetXaxis()->SetLabelSize(0.04);
    Hvert_yyp->GetXaxis()->CenterTitle();
    Hvert_yyp->GetXaxis()->SetRangeUser(xmin,xmax);
    Hvert_yyp->GetYaxis()->SetTitle("Y angle (mrad)");
    Hvert_yyp->GetYaxis()->SetLabelSize(0.04);
    Hvert_yyp->GetYaxis()->CenterTitle();
    Hvert_yyp->GetYaxis()->SetRangeUser(ymin,ymax);

    //Hlaser_imp->SetMinimum(-0.001);
    //Hlaser_imp->SetMaximum(3000.);
    Hlaser_imp->SetTitle("Laser impulse");
    Hlaser_imp->GetXaxis()->SetTitle("Horizontal imp");
    Hlaser_imp->GetXaxis()->SetLabelSize(0.04);
    Hlaser_imp->GetXaxis()->CenterTitle();
    Hlaser_imp->GetXaxis()->SetNdivisions(3,2,2,kFALSE);
    Hlaser_imp->GetXaxis()->SetRangeUser(-0.3,0.3);
    Hlaser_imp->GetYaxis()->SetTitle("Verticalimp");
    Hlaser_imp->GetYaxis()->SetLabelSize(0.04);
    Hlaser_imp->GetYaxis()->CenterTitle();
    Hlaser_imp->GetYaxis()->SetNdivisions(3,2,2,kFALSE);
    Hlaser_imp->GetYaxis()->SetRangeUser(-0.3,0.3);

//    Hhoriz_xxp->SetMinimum(-0.001);

cg->Clear();
cg->Divide(2,2);
legendg.DrawClone();

cg->cd(1);Hele_pos->DrawClone("Contz");
cg->cd(2);Hlaser_pos->DrawClone("Contz");
cg->cd(3);Hhoriz_xxp->DrawClone("Contz");
cg->cd(4);Hvert_yyp->DrawClone("Contz");
//cg->cd(4);Hlaser_imp->DrawClone("Contz");


TCanvas* cg1=new TCanvas("Canvas2","Canvas2",800,800);
Hlaser_imp->DrawClone("Colz");
}
