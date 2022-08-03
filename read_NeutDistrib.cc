void read_NeutDistrib()
{
   gStyle->SetCanvasPreferGL(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemri");

    TCanvas* cd=new TCanvas("Canvas","Canvas",800,800);
    //TLegend legend(0.08,0.85,0.2,0.9,"Z = +3.0 m");
    //char* label=new char;


TFile* in_file=new TFile("eliGammaSource.root");

TH2F* N_distrib = (TH2F*)gDirectory->Get("source_aniso");

    N_distrib->SetMaximum(150.);
    N_distrib->SetTitle("Neutron emission distrib.");
    N_distrib->GetXaxis()->SetTitle("Theta (deg)");
    N_distrib->GetXaxis()->SetLabelSize(0.04);
    N_distrib->GetXaxis()->CenterTitle();
    N_distrib->GetXaxis()->SetRangeUser(0.,180.);
    N_distrib->GetYaxis()->SetTitle("Phy (deg)");
    N_distrib->GetYaxis()->SetLabelSize(0.04);
    N_distrib->GetYaxis()->CenterTitle();
    N_distrib->GetYaxis()->SetRangeUser(0.,360.);

N_distrib->DrawClone("Colz");
}
