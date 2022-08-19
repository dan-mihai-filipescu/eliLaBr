
#include "G4VUserEventInformation.hh"
#include "eliLaBr_EventAction.hh"
#include "eliLaBr_EventInformation.hh"


using namespace std;

 TFile		 *eliLaBr_EventAction::fRootFile		= NULL;
 TTree		 *eliLaBr_EventAction::fTree			= NULL;


eliLaBr_EventAction::eliLaBr_EventAction():fFileName("nai_9_MeV.root")
{
    NMaxHit = NMAXPEAKS;
    NHistMax = 2000;
    //THistMax = 0.5*ms;
    THistMax = 1.0*ms;
    nrbin = NRBIN;
//    hist_writed = false;
    for(int i=0; i<Ndet; i++)
      {
      ss.str("");
      ss<<"EDepinTime"<<i;
      vEDepinTime[i] = new TH1D(ss.str().c_str(),ss.str().c_str(),NHistMax,0,THistMax);
      }
    ss.str("");
    ss<<"EDepinTimeTarget";
    vEDepinTimeTarget = new TH1D(ss.str().c_str(),ss.str().c_str(),NHistMax,0,THistMax);
    //vEDepinTimeTarget = new TH1D(ss.str().c_str(),ss.str().c_str(),NHistMax,115.06,115.08);
    s = new TSpectrum(NMaxHit);
    VectorOfPeaks = new Tvectors();
    VectorOfTarget = new Tvectors();
    PeaksInCounters = new TSTLvector();

    fEdep.Set(NMAXPOINTS);
    fKinEne.Set(NMAXPOINTS);
    fGlobalTime.Set(NMAXPOINTS);
//    fKinEne.Set(NMAXPOINTS);
    partType.Set(NMAXPOINTS);
//	det_nb_step.Set(NMAXPOINTS);
//	fPosX.Set(NMAXPOINTS);
//	fPosY.Set(NMAXPOINTS);
//	fPosZ.Set(NMAXPOINTS);


	fRootFile = new TFile(fFileName,"recreate");
cout<<"\n\t**************   ROOT File opened: " << fFileName << "   **************\n";
	fRootFile->cd();
	fTree = new TTree("eliLaBr","eliLaBr"); 

	fTree -> Branch ("evNo"		,&fEvNo		   ,"evNo/I"		); 
	fTree -> Branch ("nPoints"	,&fNPoints	   ,"nPoints/I"		); 
		fTree -> Branch ("Eini"     	,&Eini        		,"Eini/D"          	);
//		fTree -> Branch ("CosXini"    	,&CosXini       	,"CosXini/D"       	);
//		fTree -> Branch ("CosYini"    	,&CosYini       	,"CosYini/D"       	);
//		fTree -> Branch ("CosZini"    	,&CosZini       	,"CosZini/D"       	); 
//		fTree -> Branch ("Nmini"     	,&part        		,"Nmini/I"         	);

    fTree -> Branch ("Eini_target"	,&Eini_target	   ,"Eini_target/D"	);
    fTree -> Branch ("E_target"	,&E_target	   ,"E_target/D"	);
    fTree -> Branch ("time_ini"	,&time_ini	   ,"time_ini/D"	);
    fTree -> Branch ("tg_time"	,&tg_time	   ,"tg_time/D"	);
//    fTree -> Branch ("tg_time_min"	,&tg_time_min	   ,"tg_time_min/D"	);
//    fTree -> Branch ("tg_time_max"	,&tg_time_max	   ,"tg_time_max/D"	);

    fTree -> Branch ("mult"		,&mult		   ,"mult/I"		);
    fTree -> Branch ("mult_hit"		,&mult_hit		   ,"mult_hit/I"		);
    fTree -> Branch ("event_weight"	,&ev_weight	   ,"ev_weight/D"	);
//	fTree -> Branch ("init_ang"	,&init_ang	   ,"init_ang/D"	); 
    fTree -> Branch ("det_nb"	,det_nb		   ,"det_nb[mult]/I"	);
    fTree -> Branch ("det_hit"	,det_hit		   ,"det_hit[mult_hit]/I"	);
//	fTree -> Branch ("angle"	,angle		   ,"angle[mult]/D"	); 
    fTree -> Branch ("E_br"		,E_br		   ,"E_br[mult]/D"	);
    fTree -> Branch ("E_hit"		,E_hit		   ,"E_hit[mult_hit]/D"	);
    fTree -> Branch ("G_time"    ,G_time       ,"G_time[mult]/D");
//    fTree -> Branch ("partType"     ,partType.GetArray()   	,"partType[nPoints]/I"  );
    fTree -> Branch ("part_hit"     ,part_hit   	,"part_hit[mult_hit]/I"  );
//	fTree -> Branch ("mult_tgt"	,&mult_tgt	   ,"mult_tgt/I"	); 
//	fTree -> Branch ("tgt_nb"	,tgt_nb		   ,"tgt_nb[mult_tgt]/I"); 
//	fTree -> Branch ("E_tgt"	,&E_tgt		   ,"E_tgt/D"		); 

//	fTree -> Branch ("In_Kin_En"	,&In_Kin_En	   ,"In_Kin_En/D"	); 
//    fTree -> Branch ("eDep"		,fEdep.GetArray()  ,"eDep[nPoints]/D"	);
//    fTree -> Branch ("eKin"		,fKinEne.GetArray()  ,"eKin[nPoints]/D"	);
//	fTree -> Branch ("det_nb_step"	,det_nb_step.GetArray()  ,"det_nb_step[nPoints]/I"	); 
//	fTree -> Branch ("posX"		,fPosX.GetArray()  ,"posX[nPoints]/D"	); 
//	fTree -> Branch ("posY"		,fPosY.GetArray()  ,"posY[nPoints]/D"	); 
//	fTree -> Branch ("posZ"		,fPosZ.GetArray()  ,"posZ[nPoints]/D"	); 

//	fTree -> Branch ("mult_br"	,&mult_br	   ,"mult_br/I"		); 
    fTree -> Branch ("PeaksInCountersBranch", "PeaksInCounters"	,PeaksInCounters  , 32000, 3	);
//    fTree -> Branch ("VectorOfPeaksBranch", "VectorOfPeaks"	,VectorOfPeaks  , 32000, 3);
    fTree -> Branch ("VectorOfTargetBranch", "VectorOfTarget"	,VectorOfTarget  , 32000, 3	);


}


eliLaBr_EventAction::~eliLaBr_EventAction()
{
    fRootFile->cd();
    fTree->Write();	fRootFile->Close();
    cout<<"\n\t\t***   ROOT File - " << fFileName << " - closed   ***\t\t by EventAction \n";
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if(G4AnalysisManager::Instance())
    {
    G4cout<<"Norm = "<<norm<<G4endl;
    G4cout<<"Total electrons: "<<evNb<<G4endl;
    if(evNb > 0.)
       {
       norm = norm/evNb;
       //norm = 1.;
       analysisManager->ScaleH1(2,norm);
       analysisManager->ScaleH1(4,norm);
       analysisManager->ScaleH1(6,norm);
       analysisManager->ScaleH1(8,norm);
       G4cout<<"Incident spectra info:"<<G4endl;
       G4cout<< " MEAN = " << G4BestUnit(analysisManager->GetH1(4)->mean(), "Energy")<<G4endl;
       G4cout<< " RMS  = " << G4BestUnit(analysisManager->GetH1(4)->rms(), "Energy")<<G4endl;
       //
       // complete cleanup
       analysisManager->Write();
       analysisManager->CloseFile();
       }
    delete G4AnalysisManager::Instance();
    }
    for(int i=0; i<Ndet; i++)
      delete vEDepinTime[i];
    delete vEDepinTimeTarget;
    delete VectorOfPeaks;
    delete VectorOfTarget;
    delete PeaksInCounters;
}


void eliLaBr_EventAction::BeginOfEventAction(const G4Event* anEvent)
{
// G4cout << flush << "\n\n\n BEGIN______EVENT!!!!!!!!!!\n\n\n"  << G4endl << flush;
	fEvNo = anEvent->GetEventID();
	if ( fEvNo%100000 == 0 ) G4cout << "\nbegin event.... "<<fEvNo << G4endl << flush;
    for(int i=0; i<Ndet; i++)
      for(int j=1; j<=NHistMax; j++)
        {
	//printf("i = %d; j = %d;\n",i,j);
        vEDepinTime[i]->SetBinContent(j,0.);
	}
    for(int j=1; j<=NHistMax; j++)
      {
      //printf("i = %d; j = %d;\n",i,j);
      vEDepinTimeTarget->SetBinContent(j,0.);
      }
    VectorOfPeaks->Clear();
    VectorOfTarget->Clear();
    PeaksInCounters->Clear();
}


void eliLaBr_EventAction::EndOfEventAction(const G4Event* anEvent)
{
    G4Event* otherEvent = (G4Event*)anEvent;
    eliLaBr_EventInformation* anInfo = (eliLaBr_EventInformation*)(otherEvent->GetUserInformation());
    weight = anInfo->Get_weight();
    norm = anInfo->Get_norm();
    polAngle = anInfo->Get_angle();
    n_type = anInfo->Get_src_n_type();
    N_Det = anInfo->Get_N_Det();
    //G4cout<<"TEST: "<<anInfo->Get_norm()<<G4endl;

    evNb = anEvent->GetEventID();
//    G4cout << "Event # " << evNb << G4endl;

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4HCofThisEvent* HCE = anEvent -> GetHCofThisEvent();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

mult = 0;
mult_hit=0;
Det_touch = false;
//pt fiecare detector initializez sarcina colectata si enDep cu zero
E_target = 0.;
time_ini = tg_time_max = tg_time_min = tg_time = 0.;
tg_hit_flag = false;
for ( int b=0; b<N_Det; b++ )	// b - detector number
    {
      br_q[b]  = 0.;
      E_br[b] = 0.;
      br_time_min[b] = 0.;
      br_time_max[b] = 0.;
      G_time[b] = 0.;
      br_hit_flag[b] = false;  }

//incarc informatiile despre track
eliLaBr_TrackerHitsCollection* SHC = 0;

if(HCE)
    {
      SHC = (eliLaBr_TrackerHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("trackerCollection")));
    }	
fNPoints = SHC->entries();	

	Eini = anEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
/*	CosXini = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum().x();
	CosYini = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum().y();
	CosZini = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum().z();
*/
// init_ang = acos( anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum().z() ) * 180. / pi; 

//pentru fiecare punct de interactie generat in trackul unei particule initiale
for(int i=0; i<fNPoints; i++)
{
//In ce tip de volum a interactionat
//Cata energie a depozitat?
//In care detector?
//Ce tip de particula?
    fEdep[i] = fKinEne[i] = 0.; e_x=0; partType[i]=0;
    Vname = (*SHC)[i]->GetVname();
    fEdep[i] = (*SHC)[i]->GetEdep();
    fKinEne[i] = (*SHC)[i]->GetKinEn();
//    fKinEne[i] = (*SHC)[i]->GetKinEn();
//        detN = (*SHC)[i]->Get0RepNb();
    detN = (*SHC)[i]->GetRepNb();
    fGlobalTime[i] = (*SHC)[i]->GetGlobalTime();
    //G4cout<<"Volume = "<<Vname<<"; Nr = "<<detN <<G4endl;
//	G4ThreeVector detpos = (*SHC)[i]->GetVolPos();
//	fPosX[i] = detpos[0]; 
//	fPosY[i] = detpos[1]; 
//	fPosZ[i] = detpos[2]; 

//pos = anEvent->GetPrimaryVertex()->GetPosition();
pos = (*SHC)[i]->GetPos();
//analysisManager->FillH2(1, pos[0], pos[1]);


//	G4ThreeVector pos = (*SHC)[i]->GetPos(); 
//	fPosX[i] = pos[0]; 
//	fPosY[i] = pos[1]; 
//	fPosZ[i] = pos[2]; 
    partName = (*SHC)[i] -> GetPartName();
//    G4cout<<"partName="<< partName << "   Edep="<<fEdep[i]<<G4endl<<flush;
    if 	(partName == "neutron")		    partType[i] = 1;
	else if (partName == "gamma")		partType[i] = 2; 
    else if (partName == "e-")		    partType[i] = 3;
    else if (partName == "proton")		partType[i] = 4;
    else if (partName == "e+")          partType[i] = 5;
    else if (partName == "deuteron")    partType[i] = 6;
    else if (partName == "alpha")       partType[i] = 7;
    else if (partName == "C12[0.0]" || partName == "C13[0.0]" ) partType[i] = 8;
    else 					            partType[i] = 9;
//  G4cout << "volume " << Vname << "\tdetector " << detN << flush;
//  G4cout << "\tposition\t" << detpos[0] <<"\t"<< detpos[1] <<"\t"<< detpos[2] << flush; 
//  G4cout << "\tangle\t" << detpos.angle(G4ThreeVector(0,0,1)) * 180. / pi << G4endl << flush; 
//  G4cout << "Step status: "<<(*SHC)[i]->GetStepStat() << G4endl;
//     G4cout <<"partName="<< partName << "   Edep="<<fEdep[i]<< "  KinEn="<< (*SHC)[i]->GetKinEn()<< "\tdetector " << detN << " Step status: "<<(*SHC)[i]->GetStepStat() << G4endl;

//Daca a interactionat in material LaBr 
//adun energia depozitata la pasul asta
//la energia depozitata in total in detector 
//la trackul asta
    if 	( Vname == "he3_sensitive" )
	  {	
        br_q[detN] += fEdep[i];
        if (fEdep[i] >= 0.) {
            if(br_hit_flag[detN] == false)  {br_time_min[detN] = br_time_max[detN] = fGlobalTime[i]; br_hit_flag[detN]=true;}
            else {
                if(br_time_min[detN] > fGlobalTime[i]) br_time_min[detN] = fGlobalTime[i];
                if(br_time_max[detN] < fGlobalTime[i]) br_time_max[detN] = fGlobalTime[i];}
            vEDepinTime[detN]->Fill(fGlobalTime[i], fEdep[i]);
	    }
        //det_nb_step[i] = detN;
        }
    if 	( Vname == "target" )
      {
        E_target += fEdep[i];
        if (fEdep[i] >= 0.) {
            if(tg_hit_flag == true) {
                if(tg_time_min > fGlobalTime[i]) tg_time_min = fGlobalTime[i];
                if(tg_time_max < fGlobalTime[i]) tg_time_max = fGlobalTime[i];
		}
	    vEDepinTimeTarget->Fill(fGlobalTime[i], fEdep[i]);
        } }
//    G4cout<<"Point:"<< i << " volume " << Vname <<  "\tdetector " << detN << flush;
//    G4cout<<"partName="<< partName << "   Edep="<<fEdep[i]<<G4endl<<flush;

//Daca inainte de a intra in detector am trecut granita
if((*SHC)[i]->GetVolInStat()==1)
    {
    if 	( Vname == "he3_sensitive" ) {
    Det_touch = true;
    det_hit[mult_hit] = detN;
    E_hit[mult_hit] = (*SHC)[i]->GetKinEn();
    part_hit[mult_hit] = partType[i];
    mult_hit++;
    }
    if ( (Vname == "target") && (partType[i] == 2) ) {
        Eini_target = (*SHC)[i]->GetKinEn();
        analysisManager->FillH1(3, Eini_target);
        analysisManager->FillH1(4, Eini_target,weight);
        analysisManager->FillH1(5, polAngle);
        analysisManager->FillH1(6, polAngle, weight);
        analysisManager->FillH2(1, pos[0], pos[1],weight);
        analysisManager->FillH2(20, pos[0], pos[1],(anInfo->Get_PrimaryPartPolarization()).getX());
        analysisManager->FillH2(21, pos[0], pos[1],(anInfo->Get_PrimaryPartPolarization()).getY());
        analysisManager->FillH2(22, pos[0], pos[1],(anInfo->Get_PrimaryPartPolarization()).getZ());
        analysisManager->FillH2(2, sqrt(pos[0]*pos[0] + pos[1]*pos[1]), Eini_target, weight/sqrt(pos[0]*pos[0] + pos[1]*pos[1]));
        analysisManager->FillH2(3, pos[0], Eini_target, weight);
        analysisManager->FillH2(4, pos[1], Eini_target, weight);
        if(tg_hit_flag == false) {tg_hit_flag=true; time_ini = tg_time_min = tg_time_max = fGlobalTime[i];}}
    }

//Am terminat aici cu analiza pas cu pas a traiectoriei pentru eventul asta
}

if (E_target>0.)
  for(int i=0; i<fNPoints; i++)
    {
    if ((*SHC)[i]->GetVname() == "target")
      {
      pos = (*SHC)[i]->GetPos();
      if( ((*SHC)[i]->GetVolInStat()==1) && (partType[i] == 2) ) analysisManager->FillH2(17, pos[0], pos[1],weight);
      if (fEdep[i] >= 0.) analysisManager->FillH2(18, pos[0], pos[1],weight*fEdep[i]/E_target);
      }
    }

/*for (int i_hit=0; i_hit<mult_hit; i_hit++)
    {
    G4cout<< "Mult_hit = "<<i_hit<<"  Det_hit = "<<det_hit[i_hit]<<"  E_hit = "<<E_hit[i_hit]<<"  part_hit = "<<part_hit[i_hit]<<G4endl;
    }*/
//G4cout << "\n\n\n\tdet_1_6 "<< br_q[1][6] << G4endl << flush; 

//MULT imi arata cati detectori au sarcina acumulata nenula
if ( fNPoints > 0 ) 
	{
    if (E_target>0.)
       {//E_target = gRandom->Gaus( E_target , 0.0280*sqrt(E_target) );
        nfound = s->Search(vEDepinTimeTarget,1,"nobackground noMarkov goff nodraw",0.0001);
	xpeaks = s->GetPositionX();
	if(nfound>0) QuickSort(xpeaks, 0, nfound-1);
	imin = imax = 0;
         for (Int_t p=0;p<nfound;p++)
           {
	   //printf("Det = %d; Peak = %d;\n",d,p);
	   Float_t xp = xpeaks[p];
	   Int_t bin = vEDepinTimeTarget->GetXaxis()->FindBin(xp);
           imin = bin - nrbin;
           if(p!=0)
             {if (imin <= imax) imin = imax+1;}
              else
              {if (imin < 0) imin = 0;}
           imax = bin + nrbin;
           if(imax>(NHistMax-1)) imax = NHistMax -1; 
	   Float_t xcen = vEDepinTimeTarget->GetXaxis()->GetBinCenter(bin);
	   //if (n_type == 0) xcen = xcen - time_ini;
	   //Float_t yp = vEDepinTimeTarget->GetBinContent(bin);
	   Float_t yp = vEDepinTimeTarget->Integral(imin,imax);
           VectorOfTarget->Add(xcen,yp);
	   //printf("At index %d we have Pos %f and Val %f\n",VectorOfPeaks.GetNpeaks(),VectorOfPeaks.GetPeakPos(VectorOfPeaks.GetNpeaks()-1),VectorOfPeaks.GetPeakVal(VectorOfPeaks.GetNpeaks()-1));
	   }
        //VectorOfTarget->Clear();
	
/*       if(!hist_writed)
         if(E_target > (0.8*MeV))
	 if(E_target > (0.95*MeV))
	 if(E_target != vEDepinTimeTarget->Integral())
	 if(nfound>=2)
	 if(nfound>=0)
	 {
	 hist_writed = true;
	 fRootFile->cd();
	 vEDepinTimeTarget->Write();
	 printf("A histogram was written! :)\n");
	 //printf("Detector: %d\n",d);
	 printf("Energy: %f MeV\n",E_target/MeV);
	 printf("Energy = %f; Integral = %f\n",E_target,vEDepinTimeTarget->Integral());
	 for (Int_t p=0;p<nfound;p++) 
	   {
           Float_t xp = xpeaks[p];
           Int_t bin = vEDepinTimeTarget->GetXaxis()->FindBin(xp);
	   Float_t xcen = vEDepinTimeTarget->GetXaxis()->GetBinCenter(bin);
           Float_t yp = vEDepinTimeTarget->GetBinContent(bin);
	   printf("Bin = %f; Value = %f;\n",xcen,yp);
	   }
	 for(int i=0; i<fNPoints; i++)
	   {
	   Vname = (*SHC)[i]->GetVname();
	   if 	( Vname == "target" )
	     printf("I = %d; Edep = %f; Time = %f\n",i,(*SHC)[i]->GetEdep(),(*SHC)[i]->GetGlobalTime());
	   }
	 }*/

	
        tg_time = (tg_time_min + tg_time_max)/2.;
        analysisManager->FillH1(1, E_target);
        analysisManager->FillH1(2, E_target,weight);
        //analysisManager->FillH2(16, pos[0], pos[1],weight);
        //if(tg_hit_flag == true) analysisManager->FillH2(9, Eini_target, tg_time, weight);
    }
    /* ====== We decided NOT to Contract any more Vectors within GEANT simulations! ====== */
    /* ================ We are doing this only during Data Analysis phase ================ */
    /*VectorOfTarget->AvgTimeContract = AVGTIMECONTRACT;
    VectorOfTarget->Contract(PACKCHAN);*/
    if(tg_hit_flag == true)
        {
        //analysisManager->FillH2(2, Eini_target,sqrt(pos[0]*pos[0] + pos[1]*pos[1]),weight);
        analysisManager->FillH2(16, Eini_target, time_ini, weight);
        }

     for ( int d=0; d<N_Det; d++ )	// detector...
        {
         if ( br_q[d]>0 )
		{ // G4cout << "we're IN!!! " << G4endl << flush; 
		 det_nb[mult] = d; 
//           E_br[mult]  = gRandom->Gaus( br_q[d] , 0.0280*sqrt(br_q[d]) );  // ~2% @ 10MeV (fwhm) LaBr3
//         E_br[mult]  = gRandom->Gaus( br_q[d] , 0.0426*sqrt(br_q[d]) );  // ~5% @ 1MeV (fwhm) LaBr3
//         E_br[mult]  = gRandom->Gaus( br_q[d] , 0.0148*sqrt(br_q[d]) );  // ~3.5% @ 1MeV (fwhm) LaBr3
       
//         E_br[mult]  = gRandom->Gaus( br_q[d] , 0.0213*sqrt(br_q[d]) );  // ~5% @ 1MeV (fwhm) LaBr3
//         E_br[mult]  = gRandom->Gaus( br_q[d] , 0.0813*sqrt(br_q[d]) );  // ~?% @ 1MeV (fwhm) LaBr3
//		 E_br[mult]  = gRandom->Gaus( br_q[d] , sqrt(0.016*br_q[d]+0.0035)/2.35 );   // ~14% @ 1 MeV (fwhm) NaI
       nfound = s->Search(vEDepinTime[d],1,"nobackground noMarkov goff nodraw",0.0001);
       //Int_t nfound = s->Search(vEDepinTime[d],1,"nobackground noMarkov",0.0001);
       xpeaks = s->GetPositionX();
       if(nfound>0) QuickSort(xpeaks, 0, nfound-1);
       //VectorOfPeaks.npeaks = nfound;
       imin = imax = 0;
       for (Int_t p=0;p<nfound;p++)
         {
	 //printf("Det = %d; Peak = %d;\n",d,p);
	 Float_t xp = xpeaks[p];
	 Int_t bin = vEDepinTime[d]->GetXaxis()->FindBin(xp);
         imin = bin - nrbin;
         if(p!=0)
           {if (imin <= imax) imin = imax+1;}
           else
           {if (imin < 0) imin = 0;}
         imax = bin + nrbin;
         if(imax>(NHistMax-1)) imax = NHistMax -1; 
	 Float_t xcen = vEDepinTime[d]->GetXaxis()->GetBinCenter(bin);
	 if (n_type == 0) xcen = xcen - time_ini;
	 Float_t yp = vEDepinTime[d]->GetBinContent(bin);
         VectorOfPeaks->Add(xcen,yp);
	 //printf("At index %d we have Pos %f and Val %f\n",VectorOfPeaks.GetNpeaks(),VectorOfPeaks.GetPeakPos(VectorOfPeaks.GetNpeaks()-1),VectorOfPeaks.GetPeakVal(VectorOfPeaks.GetNpeaks()-1));
	 }
       PeaksInCounters->Add(*VectorOfPeaks);
       VectorOfPeaks->Clear();
	 
       /*if(!hist_writed)
         if(br_q[d] > (1.0*MeV))
	 if(nfound>=2)
	 if(nfound==0)
	 {
	 hist_writed = true;
	 fRootFile->cd();
	 vEDepinTime[d]->Write();
	 printf("A histogram was written! :)\n");
	 printf("Detector: %d\n",d);
	 printf("Energy: %f MeV",br_q[d]/MeV);
	 for (Int_t p=0;p<nfound;p++) 
	   {
           Float_t xp = xpeaks[p];
           Int_t bin = vEDepinTime[d]->GetXaxis()->FindBin(xp);
	   Float_t xcen = vEDepinTime[d]->GetXaxis()->GetBinCenter(bin);
           Float_t yp = vEDepinTime[d]->GetBinContent(bin);
	   printf("Bin = %f; Value = %f;\n",xcen,yp);
	   }
	 }*/
       E_br[mult]  = br_q[d];
         if (n_type == 0) G_time[mult] = (br_time_min[d] + br_time_max[d])/2. - time_ini;
           else G_time[mult] = (br_time_min[d] + br_time_max[d])/2.;
		 mult++;	
		}
	    }
    /* ====== We decided NOT to Contract any more Vectors within GEANT simulations! ====== */
    /* ================ We are doing this only during Data Analysis phase ================ */
	/*PeaksInCounters->AvgTimeContract = AVGTIMECONTRACT;
	PeaksInCounters->Contract(PACKCHAN);*/
	// fTree->Fill(); 
	}
/*if ( mult > 0 )  fTree->Fill();	*/
ev_weight = weight;
//if(Det_touch) fTree->Fill();
if (((n_type == 0)&&(tg_hit_flag=true))||(n_type == 1)) fTree->Fill();
PeaksInCounters->Clear();
VectorOfTarget->Clear();
}


void eliLaBr_EventAction::SetSaveThreshold(G4int save)
{
}

#ifdef OLD_ROOT
void eliLaBr_EventAction::QuickSort(Float_t *x, int left, int right)
#else
void eliLaBr_EventAction::QuickSort(Double_t *x, int left, int right)
#endif
{
int i = left, j = right;
#ifdef OLD_ROOT
Double_t ax;
Double_t pivot = x[(left + right) / 2];
#else
Float_t ax;
Float_t pivot = x[(left + right) / 2];
#endif

/* partition */
while (i <= j)
  {
  while (x[i] < pivot)
    i++;
  while (x[j] > pivot)
    j--;
  if (i <= j)
    {
    ax = x[i];
    x[i] = x[j];
    x[j] = ax;
    i++;
    j--;
    }
  }

/* recursion */
if (left < j)
  QuickSort(x, left, j);
if (i < right)
  QuickSort(x, i, right);
}
