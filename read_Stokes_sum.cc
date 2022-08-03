void read_Stokes_sum()
{
   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");

    //TCanvas* c=new TCanvas("Canvas","Canvas",800,800);
    //TLegend legend(0.08,0.85,0.2,0.9,"Z = +3.0 m");
    //char* label=new char;

const int n_sim = 10;
const int nbin = 1000;
float xmin=0., xmax=20., ymin=-1.01, ymax=1.01;

TFile* in_file=new TFile();


TCanvas* cc=new TCanvas("Canvas1","Canvas1",1200,800);
stringstream ss;

TH2F *Stokes1_sum  = new TH2F ("stokes1_sum",  "Stokes1",  nbin, xmin, xmax, nbin, ymin, ymax);
TH2F *Stokes2_sum  = new TH2F ("stokes2_sum",  "Stokes2",  nbin, xmin, xmax, nbin, ymin, ymax);
TH2F *Stokes3_sum  = new TH2F ("stokes3_sum",  "Stokes3",  nbin, xmin, xmax, nbin, ymin, ymax);
TH2F *PolAngle_sum = new TH2F ("PolAngle_sum", "PolAngle", nbin, xmin, xmax,  300,   0., 180.);

TH2F *Stokes1[n_sim]; 
TH2F *Stokes2[n_sim]; 
TH2F *Stokes3[n_sim]; 
TH2F *PolAngle[n_sim];

int Icnt;


    for (Icnt=1; Icnt <= n_sim; Icnt++)
    {
     ss<<"./run"<<Icnt<<"/eliGammaSource.root";
     in_file->Open(ss.str().c_str());

    Stokes1[Icnt] = (TH2F*)gDirectory->Get("stokes1");
    Stokes1_sum->Add(Stokes1[Icnt],1.);
    Stokes2[Icnt] = (TH2F*)gDirectory->Get("stokes2");
    Stokes2_sum->Add(Stokes2[Icnt],1.);
    Stokes3[Icnt] = (TH2F*)gDirectory->Get("stokes3");
    Stokes3_sum->Add(Stokes3[Icnt],1.);
    PolAngle[Icnt] = (TH2F*)gDirectory->Get("PolAngle");
    PolAngle_sum->Add(PolAngle[Icnt],1.);
    ss.str("");
    in_file->Close();
    }

//    Stokes1_sum->SetMaximum(600.);
    Stokes1_sum->SetTitle("Stokes par. #1ex");
    Stokes1_sum->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes1_sum->GetXaxis()->SetLabelSize(0.04);
    Stokes1_sum->GetXaxis()->CenterTitle();
    Stokes1_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Stokes1_sum->GetYaxis()->SetTitle("Stokes value");
    Stokes1_sum->GetYaxis()->SetLabelSize(0.04);
    Stokes1_sum->GetYaxis()->CenterTitle();
    Stokes1_sum->GetYaxis()->SetRangeUser(ymin,ymax);
//TH1F *Stokes1_proj = *Stokes1_sum->ProfileX();

//    Stokes2_sum->SetMaximum(600.);
    Stokes2_sum->SetTitle("Stokes par. #2ex");
    Stokes2_sum->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes2_sum->GetXaxis()->SetLabelSize(0.04);
    Stokes2_sum->GetXaxis()->CenterTitle();
    Stokes2_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Stokes2_sum->GetYaxis()->SetTitle("Stokes value");
    Stokes2_sum->GetYaxis()->SetLabelSize(0.04);
    Stokes2_sum->GetYaxis()->CenterTitle();
    Stokes2_sum->GetYaxis()->SetRangeUser(ymin,ymax);
//TH1F *Stokes2_proj = *Stokes2_sum->ProfileX();

//    Stokes3_sum->SetMaximum(600.);
    Stokes3_sum->SetTitle("Stokes par. #3ex");
    Stokes3_sum->GetXaxis()->SetTitle("Energy (MeV)");
    Stokes3_sum->GetXaxis()->SetLabelSize(0.04);
    Stokes3_sum->GetXaxis()->CenterTitle();
    Stokes3_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Stokes3_sum->GetYaxis()->SetTitle("Stokes value");
    Stokes3_sum->GetYaxis()->SetLabelSize(0.04);
    Stokes3_sum->GetYaxis()->CenterTitle();
    Stokes3_sum->GetYaxis()->SetRangeUser(ymin,ymax);
//TH1F *Stokes3_proj = *Stokes3_sum->ProfileX();

//    PolAngle_sum->SetMaximum(1200.);
    PolAngle_sum->SetTitle("Pol. Angle distrib.");
    PolAngle_sum->GetXaxis()->SetTitle("Energy (MeV)");
    PolAngle_sum->GetXaxis()->SetLabelSize(0.04);
    PolAngle_sum->GetXaxis()->CenterTitle();
    PolAngle_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    PolAngle_sum->GetYaxis()->SetTitle("Angle");
    PolAngle_sum->GetYaxis()->SetLabelSize(0.04);
    PolAngle_sum->GetYaxis()->CenterTitle();
    PolAngle_sum->GetYaxis()->SetRangeUser(0.,360.);
//TH1F *Stokes3_proj = *Stokes3_sum->ProfileX();

cc->Clear();
cc->Divide(3,2);
cc->cd(1);Stokes1_sum->DrawClone("Colx");
cc->cd(2);Stokes2_sum->DrawClone("Colx");
cc->cd(3);Stokes3_sum->DrawClone("Colx");
cc->cd(4);Stokes1_sum->ProfileX()->DrawClone();
cc->cd(5);Stokes2_sum->ProfileX()->DrawClone();
cc->cd(6);Stokes3_sum->ProfileX()->DrawClone();


//Hgamma_pos->Draw("Colx");

//TH1F *Polarization = (TH1F*)Hgamma_pos->ProjectionY();
//TCanvas* cp=new TCanvas("CanvasP","CanvasP",800,800);
//Hgamma_pos->ProfileX()->DrawClone();

TCanvas* cc1=new TCanvas("Canvas2","Canvas2",800,800);
PolAngle_sum->DrawClone("Colz");
}
