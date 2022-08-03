void read_GammaBeam_sum()
{
    gStyle->SetOptFit(kTRUE);
    gStyle->SetOptFit(0000);
    gStyle->SetStatBorderSize(0);
    gStyle->SetStatX(.89);
    gStyle->SetStatY(.89);
    gStyle->SetPalette(1);
/*
   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");
*/

const int n_sim = 10;
const int nbin = 300;   
float xmin=-2., xmax=2., ymin=-0.0006., ymax=0.0006;

    TCanvas* cg=new TCanvas("Canvas","Canvas",800,800);
    TLegend legendg(0.08,0.85,0.2,0.9,"Z = +3.0 m");
    char* labelg=new char;
    stringstream ss;


TFile* in_file=new TFile("eliGammaSource.root");

TH2F *Hele_pos_sum   = new TH2F("ele_pos_sum",  "Electron Beam profile", nbin, xmin, xmax, nbin, xmin, xmax);
TH2F *Hlaser_pos_sum = new TH2F("laser_pos_sum","Laser Beam profile",    nbin, xmin, xmax, nbin, xmin, xmax);
TH2F *Hhoriz_xxp_sum = new TH2F("horiz_xxp_sum","Horizontal spacephase", nbin, xmin, xmax, nbin, ymin, ymax);
TH2F *Hvert_yyp_sum  = new TH2F("vert_yyp_sum", "Vertical spacephase",   nbin, xmin, xmax, nbin, ymin, ymax);
TH2F *Hlaser_imp_sum = new TH2F("laser_imp_sum","Laser impulse",         nbin, -0.3,  0.3, nbin, -0.3,  0.3);

TH2F *Hele_pos[n_sim];
TH2F *Hlaser_pos[n_sim];
TH2F *Hhoriz_xxp[n_sim];
TH2F *Hvert_yyp[n_sim];
TH2F *Hlaser_imp[n_sim];

int Icnt;

    for (Icnt=1; Icnt <= n_sim; Icnt++)
    {
     ss<<"./run"<<Icnt<<"/eliGammaSource.root";
     in_file->Open(ss.str().c_str());
     
    Hele_pos[Icnt] = (TH2F*)gDirectory->Get("ele_pos");
    Hele_pos_sum->Add(Hele_pos[Icnt],1.);
    Hlaser_pos[Icnt] = (TH2F*)gDirectory->Get("laser_pos");
    Hlaser_pos_sum->Add(Hlaser_pos[Icnt],1.);
    Hhoriz_xxp[Icnt] = (TH2F*)gDirectory->Get("horiz_xxp");
    Hhoriz_xxp_sum->Add(Hhoriz_xxp[Icnt],1.);
    Hvert_yyp[Icnt] = (TH2F*)gDirectory->Get("vert_yyp");
    Hvert_yyp_sum->Add(Hvert_yyp[Icnt],1.);
    Hlaser_imp[Icnt] = (TH2F*)gDirectory->Get("laser_imp");
    Hlaser_imp_sum->Add(Hlaser_imp[Icnt],1.);
    ss.str("");
    in_file->Close();
    }

    Hele_pos_sum->SetMinimum(-0.001);
    //Hele_pos_sum->SetMaximum(3000.);
    Hele_pos_sum->SetTitle("Electron Beam profile");
    Hele_pos_sum->GetXaxis()->SetTitle("Horizontal (mm)");
    Hele_pos_sum->GetXaxis()->SetLabelSize(0.04);
    Hele_pos_sum->GetXaxis()->CenterTitle();
    Hele_pos_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Hele_pos_sum->GetYaxis()->SetTitle("Vertical(mm)");
    Hele_pos_sum->GetYaxis()->SetLabelSize(0.04);
    Hele_pos_sum->GetYaxis()->CenterTitle();
    Hele_pos_sum->GetYaxis()->SetRangeUser(xmin,xmax);

    Hlaser_pos_sum->SetMinimum(-0.001);
    //Hlaser_pos_sum->SetMaximum(3000.);
    //Hlaser_pos_sum->SetTitle("Electron Beam profile");
    Hlaser_pos_sum->GetXaxis()->SetTitle("Horizontal (mm)");
    Hlaser_pos_sum->GetXaxis()->SetLabelSize(0.04);
    Hlaser_pos_sum->GetXaxis()->CenterTitle();
    Hlaser_pos_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Hlaser_pos_sum->GetYaxis()->SetTitle("Vertical(mm)");
    Hlaser_pos_sum->GetYaxis()->SetLabelSize(0.04);
    Hlaser_pos_sum->GetYaxis()->CenterTitle();
    Hlaser_pos_sum->GetYaxis()->SetRangeUser(xmin,xmax);

    Hhoriz_xxp_sum->SetMinimum(-0.001);
    //Hhoriz_xxp_sum->SetMaximum(3000.);
    Hhoriz_xxp_sum->SetTitle("Horizontal spacephase");
    Hhoriz_xxp_sum->GetXaxis()->SetTitle("X position (mm)");
    Hhoriz_xxp_sum->GetXaxis()->SetLabelSize(0.04);
    Hhoriz_xxp_sum->GetXaxis()->CenterTitle();
    Hhoriz_xxp_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Hhoriz_xxp_sum->GetYaxis()->SetTitle("X angle (mrad)");
    Hhoriz_xxp_sum->GetYaxis()->SetLabelSize(0.04);
    Hhoriz_xxp_sum->GetYaxis()->CenterTitle();
    Hhoriz_xxp_sum->GetYaxis()->SetRangeUser(ymin,ymax);

    Hvert_yyp_sum->SetMinimum(-0.001);
    //Hvert_yyp_sum->SetMaximum(3000.);
    Hvert_yyp_sum->SetTitle("Vertical spacephase");
    Hvert_yyp_sum->GetXaxis()->SetTitle("Y position (mm)");
    Hvert_yyp_sum->GetXaxis()->SetLabelSize(0.04);
    Hvert_yyp_sum->GetXaxis()->CenterTitle();
    Hvert_yyp_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Hvert_yyp_sum->GetYaxis()->SetTitle("Y angle (mrad)");
    Hvert_yyp_sum->GetYaxis()->SetLabelSize(0.04);
    Hvert_yyp_sum->GetYaxis()->CenterTitle();
    Hvert_yyp_sum->GetYaxis()->SetRangeUser(ymin,ymax);

    //Hlaser_imp_sum->SetMinimum(-0.001);
    //Hlaser_imp_sum->SetMaximum(3000.);
    Hlaser_imp_sum->SetTitle("Laser impulse");
    Hlaser_imp_sum->GetXaxis()->SetTitle("Horizontal imp");
    Hlaser_imp_sum->GetXaxis()->SetLabelSize(0.04);
    Hlaser_imp_sum->GetXaxis()->CenterTitle();
    Hlaser_imp_sum->GetXaxis()->SetNdivisions(3,2,2,kFALSE);
    Hlaser_imp_sum->GetXaxis()->SetRangeUser(-0.3,0.3);
    Hlaser_imp_sum->GetYaxis()->SetTitle("Verticalimp");
    Hlaser_imp_sum->GetYaxis()->SetLabelSize(0.04);
    Hlaser_imp_sum->GetYaxis()->CenterTitle();
    Hlaser_imp_sum->GetYaxis()->SetNdivisions(3,2,2,kFALSE);
    Hlaser_imp_sum->GetYaxis()->SetRangeUser(-0.3,0.3);

//    Hhoriz_xxp->SetMinimum(-0.001);

cg->Clear();
cg->Divide(2,2);
legendg.DrawClone();

cg->cd(1);Hele_pos_sum->DrawClone("Contz");
cg->cd(2);Hlaser_pos_sum->DrawClone("Contz");
cg->cd(3);Hhoriz_xxp_sum->DrawClone("Contz");
cg->cd(4);Hvert_yyp_sum->DrawClone("Contz");
//cg->cd(4);Hlaser_imp->DrawClone("Contz");


TCanvas* cg1=new TCanvas("Canvas2","Canvas2",800,800);
Hlaser_imp_sum->DrawClone("Colz");
}
