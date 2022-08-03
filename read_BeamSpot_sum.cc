
void read_BeamSpot_sum()
{
stringstream ProjeX, ProjeY;
double x,y;

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

const int n_sim = 12;
const int nbin = 300;
float xmin=-5., xmax=5., Emin=4., Emax=20.;

    TCanvas* c=new TCanvas("Canvas","Canvas",800,800);
    TLegend legend(0.08,0.85,0.2,0.9,"Z = +0.0 m");
    char* label=new char;
    stringstream ss;


TFile* in_file=new TFile();

TH2F *Hgamma_pos_sum = new TH2F("GammaBeam-profile_sum",    "Beam profile", nbin, xmin, xmax, nbin, xmin, xmax);
TH2F *En_radius_sum =  new TH2F("En-radius_sum","Beam profile", nbin,   0., 2.5, nbin, Emin, Emax);
TH2F *En_X_sum =       new TH2F("En-X_sum",     "Beam profile", nbin, xmin, xmax, nbin, Emin, Emax);
TH2F *En_Y_sum =       new TH2F("En-Y_sum",     "Beam profile", nbin, xmin, xmax, nbin, Emin, Emax);

TH2F *Hgamma_pos[n_sim+1];
TH2F *En_radius[n_sim+1];
TH2F *En_X[n_sim+1];
TH2F *En_Y[n_sim+1];

int Icnt;


    for (Icnt=1; Icnt <= n_sim; Icnt++)
    {
     ss<<"./run"<<Icnt<<"/eliGammaSource.root";
     in_file->Open(ss.str().c_str());
     
    Hgamma_pos[Icnt] = (TH2F*)gDirectory->Get("GammaBeam-profile");
    Hgamma_pos_sum->Add(Hgamma_pos[Icnt],1.);
    cout<<"Gata prima!"<<endl;
    En_radius[Icnt] = (TH2F*)gDirectory->Get("En-radius");
    En_radius_sum->Add(En_radius[Icnt],1.);
    cout<<"Gata a doua!"<<En_radius[Icnt]->GetXaxis()->GetXmax()<< " "<<En_radius[Icnt]->GetYaxis()->GetXmin()<<endl;
    En_X[Icnt] = (TH2F*)gDirectory->Get("En-X");
    En_X_sum->Add(En_X[Icnt],1.);
    cout<<"Gata a treia!"<<endl;
    En_Y[Icnt] = (TH2F*)gDirectory->Get("En-Y");
    En_Y_sum->Add(En_Y[Icnt],1.);
    cout<<"Gata a patra!"<<endl;
    ss.str("");
    in_file->Close();
    }

    Hgamma_pos_sum->SetMinimum(0.00);
    //Hgamma_pos_sum->SetMaximum(1400.);
    Hgamma_pos_sum->SetTitle("Beam profile");
    Hgamma_pos_sum->GetXaxis()->SetTitle("Horizontal (mm)");
    Hgamma_pos_sum->GetXaxis()->SetLabelSize(0.04);
    Hgamma_pos_sum->GetXaxis()->CenterTitle();
    //Hgamma_pos_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    Hgamma_pos_sum->GetYaxis()->SetTitle("Vertical (mm)");
    Hgamma_pos_sum->GetYaxis()->SetLabelSize(0.04);
    Hgamma_pos_sum->GetYaxis()->CenterTitle();
    //Hgamma_pos_sum->GetYaxis()->SetRangeUser(xmin,xmax);

    En_radius_sum->SetMinimum(0.0);
    En_radius_sum->SetTitle("Beam profile");
    En_radius_sum->GetXaxis()->SetTitle("Beam radius (mm)");
    En_radius_sum->GetXaxis()->SetLabelSize(0.04);
    En_radius_sum->GetXaxis()->CenterTitle();
    //En_radius_sum->GetXaxis()->SetRangeUser(0.,xmax);
    En_radius_sum->GetYaxis()->SetTitle("Energy (MeV)");
    En_radius_sum->GetYaxis()->SetLabelSize(0.04);
    En_radius_sum->GetYaxis()->CenterTitle();
    //En_radius_sum->GetYaxis()->SetRangeUser(Emin,Emax);

    En_X_sum->SetMinimum(0.0);
    En_X_sum->SetTitle("Beam profile");
    En_X_sum->GetXaxis()->SetTitle("Beam X (mm)");
    En_X_sum->GetXaxis()->SetLabelSize(0.04);
    En_X_sum->GetXaxis()->CenterTitle();
    //En_X_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    En_X_sum->GetYaxis()->SetTitle("Energy (MeV)");
    En_X_sum->GetYaxis()->SetLabelSize(0.04);
    En_X_sum->GetYaxis()->CenterTitle();
    //En_X_sum->GetYaxis()->SetRangeUser(Emin,Emax);

    En_Y_sum->SetMinimum(0.0);
    En_Y_sum->SetTitle("Beam profile");
    En_Y_sum->GetXaxis()->SetTitle("Beam Y (mm)");
    En_Y_sum->GetXaxis()->SetLabelSize(0.04);
    En_Y_sum->GetXaxis()->CenterTitle();
    //En_Y_sum->GetXaxis()->SetRangeUser(xmin,xmax);
    En_Y_sum->GetYaxis()->SetTitle("Energy (MeV)");
    En_Y_sum->GetYaxis()->SetLabelSize(0.04);
    En_Y_sum->GetYaxis()->CenterTitle();
    //En_Y_sum->GetYaxis()->SetRangeUser(Emin,Emax);

c->Clear();
c->Divide(2,2);
legend.DrawClone();

c->cd(1);Hgamma_pos_sum->DrawClone("Colz");
c->cd(2);En_radius_sum->DrawClone("Colz");
c->cd(3);En_X_sum->DrawClone("Colz");
c->cd(4);En_Y_sum->DrawClone("Colz");
//delete in_file;

TH1F* averageEnergy = new TH1F ( "Average Energy", "Gamma beam energy distribution on target surface - Average Energy",nbin, 0., 2.5);
TH1F* projection_on_r = new TH1F ( "Projection", "Gamma beam energy distribution on target surface -  Projection",nbin, 0., 2.5);

averageEnergy->Add(En_radius_sum->ProfileX(),1.);
averageEnergy->GetXaxis()->SetTitle("Beam radius (mm)");
averageEnergy->GetXaxis()->SetLabelSize(0.04);
averageEnergy->GetXaxis()->CenterTitle();
averageEnergy->GetYaxis()->SetRangeUser(Emin,Emax);
projection_on_r->Add(En_radius_sum->ProjectionX(),1.);
projection_on_r->GetXaxis()->SetTitle("Beam radius (mm)");
projection_on_r->GetXaxis()->SetLabelSize(0.04);
projection_on_r->GetXaxis()->CenterTitle();

TCanvas* projections=new TCanvas("Canvas1","Canvas1",800,400);
projections->Clear();
projections->Divide(2);
projections->cd(1);averageEnergy->DrawClone();
projections->cd(2);projection_on_r->DrawClone();

TH1F* projection_on_X = new TH1F ( "ProjectionX", "Gamma beam distribution on target surface -  ProjectionX",nbin, xmin, xmax );
TH1F* projection_on_Y = new TH1F ( "ProjectionY", "Gamma beam distribution on target surface -  ProjectionY",nbin, xmin, xmax );

projection_on_X->Add(Hgamma_pos_sum->ProjectionX(),1.);
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

projection_on_Y->Add(Hgamma_pos_sum->ProjectionY(),1.);
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
