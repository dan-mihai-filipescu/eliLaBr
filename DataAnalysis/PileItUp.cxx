#define PileItUp_cxx
// The class definition in PileItUp.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("PileItUp.C")
// Root > T->Process("PileItUp.C","some options")
// Root > T->Process("PileItUp.C+")
//

#define PileItUp_cxx

#include "PileItUp.h"

using namespace std;

void PileItUp::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   FILE *inputFile;
   inputFile = fopen("PileItUp.in","r");
   if(!fscanf(inputFile,"%f",&PileUpInp))
     {
     printf(" ERROR - invalid Pile-up number\n");
     exit(0);
     }
   printf("PileUpInp = %f\n",PileUpInp);
   if(PileUpInp<0)
     if((int)roundf(PileUpInp)==0)
       {
       printf(" ERROR - invalid Pile-up number\n");
       exit(0);
       }
   if(PileUpInp<0) ev_index = abs((int)roundf(PileUpInp));
     else
     {
     PoissonAvg = PileUpInp;
     do
      {ev_index = gRandom->Poisson( PoissonAvg );}
      while(ev_index==0);
     }

   if(!fscanf(inputFile,"%i",&runnum)) {printf(" ERROR - invalid GASP run number\n"); exit(0);}
   if(fscanf(inputFile,"%i %i %i",&NdetA, &NdetB, &NdetC) != 3) {printf(" ERROR - invalid numbers of detectors on each ring\n"); exit(0);}
   if(!fscanf(inputFile,"%i",&PartType)) {printf(" ERROR - invalid PartType number\n"); exit(0);}
   if(!fscanf(inputFile,"%i",&Option_Time)) {printf(" ERROR - invalid Option_Time number\n"); exit(0);}
   if(!fscanf(inputFile,"%i",&Option_Record)) {printf(" ERROR - invalid Option_Record number\n"); exit(0);}
   if(!fscanf(inputFile,"%i",&Option_OutFiles)) {printf(" ERROR - invalid Option_OutFiles number\n"); exit(0);}
   /* ====================================== OPTIONS ====================================== */
   //PartType = 1;  // 1 - neutron
                    // 2 - gamma
                    // 3 - e-
                    // 4 - proton
                    // 5 - e+
                    // 6 - deuteron
                    // 7 - alpha
                    // 8 - C12[0.0] & C13[0.0]
                    // 9 - other
   //Option_Time = 1; // 0 - Record Time spectra only if Selected particle crosses the border
                      // 1 - Record Time spectra for All the particles
   //Option_Record = 1;   // 0 - Record Incident  energy into Particle files
                          // 1 - Record Deposited energy into Particle files
   //Option_OutFiles = 0;   // 0 - Do NOT Output individual detector spectra files
                            // 1 - Output individual detector spectra files
   /* ===================================================================================== */
   if(!fscanf(inputFile,"%i",&NHistMax)) {printf(" ERROR - invalid NHistMax number\n"); exit(0);}
   if(!fscanf(inputFile,"%lf",&EneDivFact)) {printf(" ERROR - invalid EneDivFact\n"); exit(0);}
   if(!fscanf(inputFile,"%lf",&TimeMultFact)) {printf(" ERROR - invalid TimeMultFact\n"); exit(0);}
   if(!fscanf(inputFile,"%f",&PackChan)) {printf(" ERROR - invalid PackChan\n"); exit(0);}
   if(!fscanf(inputFile,"%d",&AvgTimeContract)) {printf(" ERROR - invalid AvgTimeContract\n"); exit(0);}
   if(!fscanf(inputFile,"%d",&HistEneTime)) {printf(" ERROR - invalid HistEneTime\n"); exit(0);}
   fclose(inputFile);

   EHistMax = double(NHistMax)/EneDivFact;
   THistMax = double(NHistMax)*TimeMultFact;
   iRecord = 1;
   iRecord_indv = 1;
   Ndet = NdetA+NdetB+NdetC;

   grh = (GaspRecHeader*) calloc (1, sizeof(GaspRecHeader));
   G = (GaspEvent*) calloc(BUFFER, sizeof(GaspEvent));
   for(z=0; z<BUFFER; z++) G[z].elem = (short *) calloc(Ndet*(NParam+1)+NDetTypes+2,sizeof(short));
   zero = (short*) calloc(BUFFER, sizeof(short));
   sprintf(printfile, "gasp.out");
   fp=fopen(printfile, "wt");
   grh[0].head[0]=32;
   //grh[0].head[1]=1;
   grh[0].head[1]=iRecord;
   grh[0].head[2]=runnum;
   grh[0].head[3]=0x4748;  // HG = 0x4748, XG = 0x4758
   grh[0].head[4]=16;

   recSize=16;
   z = 0;
   iGasp = 0;
   fwrite (&grh[0], 1, grh[0].head[4]*sizeof(short), fp);
   for(ii=0; ii<BUFFER; ii++) zero[ii]=0;
   fwrite (zero, 1, (grh[0].head[0]*1024/sizeof(short)-recSize)*sizeof(short), fp);
   
   iRecord++;
   grh[0].head[0]=16;
   grh[0].head[1]=iRecord;
   grh[0].head[3]=0x4758;  // HG = 0x4748, XG = 0x4758
   fwrite (&grh[0], 1, grh[0].head[4]*sizeof(short), fp);

   grh_indv = (GaspRecHeader*) calloc (1, sizeof(GaspRecHeader));
   G_indv = (GaspEvent*) calloc(BUFFER, sizeof(GaspEvent));
   for(z_indv=0; z_indv<BUFFER; z_indv++) G_indv[z_indv].elem = (short *) calloc(Ndet*(NParam+1)+NDetTypes_indv+2,sizeof(short));
   sprintf(printfile, "gasp_indv.out");
   fp_indv=fopen(printfile, "wt");
   grh_indv[0].head[0]=32;
   //grh[0].head[1]=1;
   grh_indv[0].head[1]=iRecord_indv;
   grh_indv[0].head[2]=runnum;
   grh_indv[0].head[3]=0x4748;  // HG = 0x4748, XG = 0x4758
   grh_indv[0].head[4]=16;

   recSize_indv=16;
   z_indv = 0;
   fwrite (&grh_indv[0], 1, grh_indv[0].head[4]*sizeof(short), fp_indv);
   for(ii=0; ii<BUFFER; ii++) zero[ii]=0;
   fwrite (zero, 1, (grh_indv[0].head[0]*1024/sizeof(short)-recSize_indv)*sizeof(short), fp_indv);
   
   iRecord_indv++;
   grh_indv[0].head[0]=16;
   grh_indv[0].head[1]=iRecord_indv;
   grh_indv[0].head[3]=0x4758;  // HG = 0x4748, XG = 0x4758
   fwrite (&grh_indv[0], 1, grh_indv[0].head[4]*sizeof(short), fp_indv);

        if (PartType == 1) ps<<"_n";
   else if (PartType == 2) ps<<"_g";
   else if (PartType == 3) ps<<"_em";
   else if (PartType == 4) ps<<"_p";
   else if (PartType == 5) ps<<"_ep";
   else if (PartType == 6) ps<<"_d";
   else if (PartType == 7) ps<<"_a";
   else if (PartType == 8) ps<<"_C";

   for (Icnt=0; Icnt < Ndet; Icnt++)
      {
      ss<<"ESingle"<<Icnt;
      vESingle[Icnt] = new TH1F(ss.str().c_str(),ss.str().c_str(),NHistMax,0,EHistMax);
      ss.str("");
      ss<<"EPileUp"<<Icnt;
      vEPileUp[Icnt] = new TH1F(ss.str().c_str(),ss.str().c_str(),NHistMax,0,EHistMax);
      ss.str("");
      ss<<"ESingleInc"<<Icnt;
      vESingleInc[Icnt] = new TH1F(ss.str().c_str(),ss.str().c_str(),NHistMax,0,EHistMax);
      ss.str("");
      ss<<"Time"<<Icnt;
      vTime[Icnt] = new TH1F(ss.str().c_str(),ss.str().c_str(),NHistMax,0,THistMax);
      ss.str("");
      }
   //for (Icnt=0; Icnt<Ndet; Icnt++) {cout<<"Icnt = "<<Icnt<<endl;vESingle[Icnt].Fill(50);cout<<vESingle[Icnt].GetBinContent(50)<<endl;}
   HistS = new TH1F("HistS","HistS",NHistMax,0,EHistMax);
   HistSp = new TH1F("HistSp","HistSp",NHistMax,0,EHistMax);
   HistP = new TH1F("HistP","HistP",NHistMax,0,EHistMax);
   HistT = new TH1F("HistT","HistT",NHistMax,0,THistMax);
   CHistS = new TH1F("CHistS","CHistS",NHistMax,0,EHistMax);
   CHistSp = new TH1F("CHistSp","CHistSp",NHistMax,0,EHistMax);
   CHistP = new TH1F("CHistP","CHistP",NHistMax,0,EHistMax);
   CHistT = new TH1F("CHistT","CHistT",NHistMax,0,THistMax);
   BHistS = new TH1F("BHistS","BHistS",NHistMax,0,EHistMax);
   BHistSp = new TH1F("BHistSp","BHistSp",NHistMax,0,EHistMax);
   BHistP = new TH1F("BHistP","BHistP",NHistMax,0,EHistMax);
   BHistT = new TH1F("BHistT","BHistT",NHistMax,0,THistMax);
   AHistS = new TH1F("AHistS","AHistS",NHistMax,0,EHistMax);
   AHistSp = new TH1F("AHistSp","AHistSp",NHistMax,0,EHistMax);
   AHistP = new TH1F("AHistP","AHistP",NHistMax,0,EHistMax);
   AHistT = new TH1F("AHistT","AHistT",NHistMax,0,THistMax);
   TargetP = new TH1F("TargetP","TargetP",NHistMax,0,EHistMax*100.);
   ET_TotalProj = new TH2F("ET_TotalProj","ET_TotalProj",NHistMax,0,EHistMax,NHistMax,0,THistMax);

   prev_ev = 0; 
   for (Icnt=0; Icnt < Ndet; Icnt++) {sum[Icnt]  = 0.;}
   counter = 0; 
   ev = 0;
   mult0 = 0;

   VectorOfPeaks_Tot = new Tvectors();
   VectorOfPeaks_R1 = new Tvectors();
   VectorOfPeaks_R2 = new Tvectors();
   VectorOfPeaks_R3 = new Tvectors();
   VectorOfPeaks = new Tvectors();
   PileUpVector = new TSTLvector();
   EmptyVector = new Tvectors();
   PileUpTarget =  new Tvectors();

   for (Icnt = 0; Icnt < Ndet; Icnt++)
     {
     PileUpVector->Add(*EmptyVector);
     }

   TString option = GetOption();

}

void PileItUp::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t PileItUp::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either PileItUp::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   b_E_target->GetEntry(entry);
   b_tg_time->GetEntry(entry);
   b_mult->GetEntry(entry); 
   b_mult_hit->GetEntry(entry);
   b_ev_weight->GetEntry(entry);
   b_E_br->GetEntry(entry); 
   b_E_hit->GetEntry(entry);
   b_G_time->GetEntry(entry);
   b_det_nb->GetEntry(entry);
   b_det_hit->GetEntry(entry);
   b_part_hit->GetEntry(entry);
   b_evNo->GetEntry(entry);
   b_PeaksInCountersBranch->GetEntry(entry);
   b_VectorOfTargetBranch->GetEntry(entry);
   //m.GetEntry(i);
   //Show(entry);
   //PeaksInCounters = (m.PeaksInCountersBranch);
   //printf("evNo = %d\n",evNo);
   //printf("fNelements = %d\n",PeaksInCountersBranch->GetNelements());

   /* ====== We decided NOT to Contract any more Vectors within GEANT simulations! ====== */
   /* ================ We are doing this only during Data Analysis phase ================ */
   PeaksInCountersBranch->AvgTimeContract = AvgTimeContract;
   PeaksInCountersBranch->Contract(PackChan);
   VectorOfTargetBranch->AvgTimeContract = AvgTimeContract;
   VectorOfTargetBranch->Contract(PackChan);

   VectorOfPeaks_Tot->Clear();
   VectorOfPeaks_R1->Clear();
   VectorOfPeaks_R2->Clear();
   VectorOfPeaks_R3->Clear();

   //PileUpTarget->Add(*VectorOfTargetBranch);
   if(HistEneTime)
     {
     for (Icnt=0; Icnt < VectorOfTargetBranch->GetNpeaks(); Icnt++)
       {
       //Energy  = VectorOfTargetBranch->GetPeakVal(Icnt);
       Energy  = gRandom->Gaus( VectorOfTargetBranch->GetPeakVal(Icnt) , sqrt(0.016*VectorOfTargetBranch->GetPeakVal(Icnt)+0.0035)/2.35 );   // ~14% @ 1 MeV (fwhm) NaI
       PileUpTarget->Add(VectorOfTargetBranch->GetPeakPos(Icnt),Energy);
       }
     }
     else
     {
     Energy  = gRandom->Gaus( E_target , sqrt(0.016*E_target+0.0035)/2.35 );   // ~14% @ 1 MeV (fwhm) NaI
     PileUpTarget->Add(tg_time,Energy);
     }

   for (Icnt=0; Icnt < Ndet; Icnt++) { Ene[Icnt] = 0.; det_status[Icnt] = false;}

   mult_indv = 0;
   if ( mult > 0 ) 
     {
     for (Icnt = 0; Icnt < mult; Icnt++)
       {
       VectorOfPeaks->Clear();
       *VectorOfPeaks = PeaksInCountersBranch->Get(Icnt);
       if(VectorOfPeaks->GetNpeaks()>0) mult_indv++;
       //printf("Icnt = %d\n",Icnt);
       /*if(VectorOfPeaks->GetNpeaks()==0)
         {
         printf("Event = %i, Number of peaks = %i\n",evNo, PeaksInCountersBranch->Get(Icnt).GetNpeaks());
         printf("Icnt = %i; Mult = %i\n",Icnt,mult);
         }*/
       VectorOfPeaks_Tot->Add(*VectorOfPeaks);
       if                                  (det_nb[Icnt]< NdetA              ) VectorOfPeaks_R1->Add(*VectorOfPeaks);
       if( (NdetA       <=det_nb[Icnt]) && (det_nb[Icnt]<(NdetA+NdetB))      ) VectorOfPeaks_R2->Add(*VectorOfPeaks);
       if(((NdetA+NdetB)<=det_nb[Icnt]) && (det_nb[Icnt]<(NdetA+NdetB+NdetC))) VectorOfPeaks_R3->Add(*VectorOfPeaks);
       PileUpVector->Set(det_nb[Icnt],*VectorOfPeaks);
       if(HistEneTime)
         {
         for (Jcnt = 0; Jcnt < VectorOfPeaks->GetNpeaks(); Jcnt++)
	   {
           Ene[det_nb[Icnt]] = VectorOfPeaks->GetPeakVal(Jcnt);
           vESingle[det_nb[Icnt]]->Fill(Ene[det_nb[Icnt]],ev_weight);
           sum[det_nb[Icnt]] += Ene[det_nb[Icnt]];
           //if(0<=det_nb[Icnt] && det_nb[Icnt]<=3)
           ET_TotalProj->Fill(Ene[det_nb[Icnt]],VectorOfPeaks->GetPeakPos(Jcnt),ev_weight);
           if(Option_Time==1) vTime[det_nb[Icnt]]->Fill(VectorOfPeaks->GetPeakPos(Jcnt),ev_weight);
           //cout<<"New EVENT!"<<endl;
           //cout<<"Icnt = "<<Icnt<<" DetNb "<<det_nb[Icnt]<<" DepEne "<<E_br[Icnt]<<endl;
           }
	 }
	 else
	 {
	 Ene[det_nb[Icnt]] = E_br[Icnt];
	 vESingle[det_nb[Icnt]]->Fill(Ene[det_nb[Icnt]],ev_weight);
	 sum[det_nb[Icnt]] += Ene[det_nb[Icnt]];
	 ET_TotalProj->Fill(Ene[det_nb[Icnt]],G_time[Icnt],ev_weight);
	 if(Option_Time==1) vTime[det_nb[Icnt]]->Fill(G_time[Icnt],ev_weight);
	 }
       }
     //VectorOfPeaks_Tot->Contract(PackChan);
     //printf("R1 npeaks = %d\n",VectorOfPeaks_R1->GetNpeaks());
     //printf("R2 npeaks = %d\n",VectorOfPeaks_R2->GetNpeaks());
     //printf("R3 npeaks = %d\n",VectorOfPeaks_R3->GetNpeaks());
     VectorOfPeaks_Tot->QuickSort();
     //if(VectorOfPeaks_R1->GetNpeaks()>0) 
     VectorOfPeaks_R1->QuickSort();
     //if(VectorOfPeaks_R2->GetNpeaks()>0) 
     VectorOfPeaks_R2->QuickSort();
     //if(VectorOfPeaks_R3->GetNpeaks()>0) 
     VectorOfPeaks_R3->QuickSort();
     /*for (Icnt = 0; Icnt < VectorOfPeaks_Tot->GetNpeaks(); Icnt++)
       ET_TotalProj->Fill(VectorOfPeaks_Tot->GetPeakVal(Icnt),VectorOfPeaks_Tot->GetPeakPos(Icnt),ev_weight);*/
     }

   /************************* Write GASP file *******************************************/
   iGasp=1+2+NDetTypes; //includes separator, header and (dettype) detector types
//   for(j=0; j<NDetTypes; j++) 
//     {
//     if (j==0)
   if(VectorOfTargetBranch->GetNpeaks() > 0)
     {
     G[z].elem[1] = VectorOfTargetBranch->GetPeakVal(0)*EneDivFact/10.;
     G[z].elem[2] = VectorOfTargetBranch->GetPeakPos(0)/TimeMultFact;
     }
     else
     {
     G[z].elem[1] = 0;
     G[z].elem[2] = 0;
     }
   //G[z].elem[2] = (unsigned int)mult;
   G[z].elem[3] = (unsigned int)VectorOfPeaks_Tot->GetNpeaks();
   G[z].elem[4] = (unsigned int)VectorOfPeaks_R1->GetNpeaks();
   G[z].elem[5] = (unsigned int)VectorOfPeaks_R2->GetNpeaks();
   G[z].elem[6] = (unsigned int)VectorOfPeaks_R3->GetNpeaks();
   //for (Icnt=0; Icnt<mult; Icnt++)
   for (Icnt=0; Icnt<VectorOfPeaks_Tot->GetNpeaks(); Icnt++)
     {
     //G[z].elem[iGasp++] = (unsigned int)det_nb[Icnt];
     G[z].elem[iGasp++] = (unsigned int)Icnt;
     //G[z].elem[iGasp++] = (unsigned int)(E_br[Icnt]*EneDivFact);
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_Tot->GetPeakVal(Icnt)*EneDivFact);
     //G[z].elem[Gaspi++] = (unsigned int)(G_time[Icnt]/TimeMultFact);
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_Tot->GetPeakPos(Icnt)/TimeMultFact);
     }
   for (Icnt=0; Icnt<VectorOfPeaks_R1->GetNpeaks(); Icnt++)
     {
     G[z].elem[iGasp++] = (unsigned int)Icnt;
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_R1->GetPeakVal(Icnt)*EneDivFact);
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_R1->GetPeakPos(Icnt)/TimeMultFact);
     }
   for (Icnt=0; Icnt<VectorOfPeaks_R2->GetNpeaks(); Icnt++)
     {
     G[z].elem[iGasp++] = (unsigned int)Icnt;
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_R2->GetPeakVal(Icnt)*EneDivFact);
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_R2->GetPeakPos(Icnt)/TimeMultFact);
     }
   for (Icnt=0; Icnt<VectorOfPeaks_R3->GetNpeaks(); Icnt++)
     {
     G[z].elem[iGasp++] = (unsigned int)Icnt;
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_R3->GetPeakVal(Icnt)*EneDivFact);
     G[z].elem[iGasp++] = (unsigned int)(VectorOfPeaks_R3->GetPeakPos(Icnt)/TimeMultFact);
     }
   G[z].evSize = iGasp-1; // does not include the separator
   G[z].elem[0] = 0xf000 + G[z].evSize; // separator GaspEvent
   
//   for(j=0;j<z;j++) printf("z = %i; j = %i; G = %i\n",z,j,G[z].elem[j]);
   
   if (recSize < (grh[0].head[0]*1024/(short)sizeof(short) - (G[z].evSize+1)) )
     {
     //z++;
     recSize += (G[z].evSize+1);
     }
     else
     {
     for (j=0; j<z; j++) fwrite (G[j].elem, 1, (G[j].evSize+1)*sizeof(short), fp);
     fwrite (zero, 1, (grh[0].head[0]*1024/sizeof(short)-recSize)*sizeof(short), fp);
     G[0] = G[z];
     z = 0;
     iRecord++;
     grh[0].head[1] =iRecord; 
     fwrite (&grh[0], 1, grh[0].head[4]*sizeof(short), fp);
     recSize = 16;
     recSize += (G[0].evSize+1);
     }
   z++;
   /*------------------------------------------------------------------*/
   iGasp=1+1+NDetTypes_indv; //includes separator, header and (dettype) detector types
//   for(j=0; j<NDetTypes; j++) 
//     {
//     if (j==0)
   G_indv[z_indv].elem[1] = 1.0;
   if(HistEneTime)
     G_indv[z_indv].elem[2] = (unsigned int)mult_indv;
     else
     G_indv[z_indv].elem[2] = (unsigned int)mult;
   
//   if ( mult_indv > 0 )
   for (Icnt=0; Icnt<mult; Icnt++)
     {
     if(HistEneTime)
       {
       VectorOfPeaks->Clear();
       *VectorOfPeaks = PeaksInCountersBranch->Get(Icnt);
       if(VectorOfPeaks->GetNpeaks()>0)
         {
         G_indv[z_indv].elem[iGasp++] = (unsigned int)det_nb[Icnt];
         G_indv[z_indv].elem[iGasp++] = (unsigned int)(VectorOfPeaks->GetPeakVal(0)*EneDivFact);
         G_indv[z_indv].elem[iGasp++] = (unsigned int)(VectorOfPeaks->GetPeakPos(0)/TimeMultFact);
         }
       }
       else
       {
       G_indv[z_indv].elem[iGasp++] = (unsigned int)det_nb[Icnt];
       G_indv[z_indv].elem[iGasp++] = (unsigned int)(E_br[Icnt]*EneDivFact);
       G_indv[z_indv].elem[iGasp++] = (unsigned int)(G_time[Icnt]/TimeMultFact);
       }
     }
   G_indv[z_indv].evSize = iGasp-1; // does not include the separator
   G_indv[z_indv].elem[0] = 0xf000 + G_indv[z_indv].evSize; // separator GaspEvent
   
//   for(j=0;j<z;j++) printf("z = %i; j = %i; G = %i\n",z,j,G[z].elem[j]);
   
   if (recSize_indv < (grh_indv[0].head[0]*1024/(short)sizeof(short) - (G_indv[z_indv].evSize+1)) )
     {
     //z++;
     recSize_indv += (G_indv[z_indv].evSize+1);
     }
     else
     {
     for (j=0; j<z_indv; j++) fwrite (G_indv[j].elem, 1, (G_indv[j].evSize+1)*sizeof(short), fp_indv);
     fwrite (zero, 1, (grh_indv[0].head[0]*1024/sizeof(short)-recSize_indv)*sizeof(short), fp_indv);
     G_indv[0] = G_indv[z];
     z_indv = 0;
     iRecord_indv++;
     grh_indv[0].head[1] =iRecord_indv; 
     fwrite (&grh_indv[0], 1, grh_indv[0].head[4]*sizeof(short), fp_indv);
     recSize_indv = 16;
     recSize_indv += (G_indv[0].evSize+1);
     }
   z_indv++;
   /********************************************************************/

   if (mult_hit > 0) for (Icnt=0; Icnt < mult_hit; Icnt++)
     {
     //cout<<"MultHit = "<<Icnt<<"  particle = "<<part_hit[Icnt] <<endl;
     if(part_hit[Icnt]==PartType)
       {
       if(Option_Record == 0) vESingleInc[det_hit[Icnt]]->Fill(E_hit[Icnt],ev_weight);
       if (mult > 0)
         for (ii=0; ii<mult; ii++)
	   if(det_nb[ii]==det_hit[Icnt] && det_status[det_hit[Icnt]]==false)
	     {
	     if(Option_Record == 1) vESingleInc[det_hit[Icnt]]->Fill(E_br[ii],ev_weight);
             if(Option_Time==0) vTime[det_hit[Icnt]]->Fill(G_time[ii],ev_weight);
             det_status[det_hit[Icnt]]=true;
             //cout<<"Icnt1 = "<<Icnt<<" ii = "<<ii<<" DetNb "<<det_hit[Icnt]<<" DepEne "<<E_br[ii]<<endl;
             }
       }
     }

   // this_ev = evNo%5; 
   if ( counter == ev_index -1 )
     {
     PileUpVector->AvgTimeContract = AvgTimeContract;
     PileUpVector->Contract(PackChan);
     PileUpTarget->AvgTimeContract = AvgTimeContract;
     PileUpTarget->Contract(PackChan);
     //for (Icnt=0; Icnt < Ndet; Icnt++) {if (sum[Icnt]>0.) {vEPileUp[Icnt]->Fill(sum[Icnt],ev_weight);};}
     for (Icnt=0; Icnt < Ndet; Icnt++)
       {
       VectorOfPeaks->Clear();
       *VectorOfPeaks = PileUpVector->Get(Icnt);
       for (Jcnt = 0; Jcnt < VectorOfPeaks->GetNpeaks(); Jcnt++)
         {
         vEPileUp[Icnt]->Fill(VectorOfPeaks->GetPeakVal(Jcnt),ev_weight);
         }
       }
     for (Jcnt = 0; Jcnt < PileUpTarget->GetNpeaks(); Jcnt++)
       {
       TargetP->Fill(PileUpTarget->GetPeakVal(Jcnt),ev_weight);
       //printf("Jcnt = %d; PileUp val = %f; PileUp pos = %f;\n",Jcnt,PileUpTarget->GetPeakVal(Jcnt),PileUpTarget->GetPeakPos(Jcnt));
       }
     for (Icnt=0; Icnt < Ndet; Icnt++) {sum[Icnt]  = 0.;}
     PileUpVector->ClearElements();
     PileUpTarget->Clear();
     counter = 0;
     if(PileUpInp>0)
       do
       {ev_index = gRandom->Poisson( PoissonAvg );}
       while(ev_index==0);
     //printf("EvIndex = %d;\n",ev_index);
     }
     else
     {
     //if ( mult > 0 ) {for (Icnt=0; Icnt < Ndet; Icnt++) if (Ene[Icnt]>0.){{sum[Icnt] += Ene[Icnt]; };};}
     counter++;
     }

   if(mult==0) mult0++;
   if ( ev%1000000 == 0 )  cout << "\n processing ev. # " << ev ; 
   //cout << "\n processing ev. # " << ev <<endl;
   ev++;

   return kTRUE;
}

void PileItUp::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void PileItUp::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   //TCanvas* cp=new TCanvas("Canvas","Canvas",800,800);
   new TCanvas("Canvas","Canvas",800,800);
   ET_TotalProj->GetXaxis()->SetTitle("Energy (MeV)");
   ET_TotalProj->GetXaxis()->CenterTitle();
   ET_TotalProj->GetYaxis()->SetTitle("Time (ns)");
   ET_TotalProj->GetYaxis()->CenterTitle();
   ET_TotalProj->Draw("Colx");

   double y;
   cout << "\n Mult0 "<< mult0<<endl;
   for (Icnt=0; Icnt < Ndet; Icnt++)
     {
     if(Option_OutFiles == 1)
       {
       ss<<"ESingle"<<Icnt;
       ofstream ofiles(ss.str().c_str(), ios::out);
       ss.str("");
       ss<<"EPileUp"<<Icnt;
       ofstream ofile(ss.str().c_str(), ios::out);
       ss.str("");
       if(Option_Record == 0) ss<<"ESingleInc"<<ps.str().c_str()<<Icnt;
         else ss<<"ESingle"<<ps.str().c_str()<<Icnt;
       ofstream ofilesi(ss.str().c_str(), ios::out);
       ss.str("");
       if(Option_Time==0) ss<<"Time"<<ps.str().c_str()<<Icnt;
         else ss<<"Time"<<Icnt;
       ofstream ofilet(ss.str().c_str(), ios::out);
       ss.str("");

       for (iChan = 0; iChan < NHistMax; iChan++)
         {
         y=vEPileUp[Icnt]->GetBinContent(iChan);
         ofile<<y<<endl;
         y=vESingle[Icnt]->GetBinContent(iChan);
         ofiles<<y<<endl;
         y=vESingleInc[Icnt]->GetBinContent(iChan);
         ofilesi<<y<<endl;
         y=vTime[Icnt]->GetBinContent(iChan);
         ofilet<<y<<endl;
         }
       ofile.close();  
       ofiles.close();  
       ofilesi.close();
       ofilet.close();
       }
     HistS->Add(vESingle[Icnt],1.);
     HistSp->Add(vESingleInc[Icnt],1.);
     HistP->Add(vEPileUp[Icnt],1.);
     HistT->Add(vTime[Icnt],1.);
     if(0<=Icnt && Icnt<NdetA)
       {
       AHistS->Add(vESingle[Icnt],1.);
       AHistSp->Add(vESingleInc[Icnt],1.);
       AHistP->Add(vEPileUp[Icnt],1.);
       AHistT->Add(vTime[Icnt],1.);
       }
     if(NdetA<=Icnt && Icnt<(NdetA+NdetB))
       {
       BHistS->Add(vESingle[Icnt],1.);
       BHistSp->Add(vESingleInc[Icnt],1.);
       BHistP->Add(vEPileUp[Icnt],1.);
       BHistT->Add(vTime[Icnt],1.);
       }
     if((NdetA+NdetB)<=Icnt && Icnt<(NdetA+NdetB+NdetC))
       {
       CHistS->Add(vESingle[Icnt],1.);
       CHistSp->Add(vESingleInc[Icnt],1.);
       CHistP->Add(vEPileUp[Icnt],1.);
       CHistT->Add(vTime[Icnt],1.);
       }
     }

   ofstream ofileHistS("HistS", ios::out);
   if(Option_Record == 0) ss<<"HistSInc"<<ps.str().c_str();
     else ss<<"HistS"<<ps.str().c_str();
   ofstream ofileHistSp(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileHistP("HistP", ios::out);
   if(Option_Time == 0) ss<<"HistT"<<ps.str().c_str();
     else ss<<"HistT";
   ofstream ofileHistT(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileCHistS("CHistS", ios::out);
   if(Option_Record == 0) ss<<"CHistSInc"<<ps.str().c_str();
     else ss<<"CHistS"<<ps.str().c_str();
   ofstream ofileCHistSp(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileCHistP("CHistP", ios::out);
   if(Option_Time == 0) ss<<"CHistT"<<ps.str().c_str();
     else ss<<"CHistT";
   ofstream ofileCHistT(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileBHistS("BHistS", ios::out);
   if(Option_Record == 0) ss<<"BHistSInc"<<ps.str().c_str();
     else ss<<"BHistS"<<ps.str().c_str();
   ofstream ofileBHistSp(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileBHistP("BHistP", ios::out);
   if(Option_Time == 0) ss<<"BHistT"<<ps.str().c_str();
     else ss<<"BHistT";
   ofstream ofileBHistT(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileAHistS("AHistS", ios::out);
   if(Option_Record == 0) ss<<"AHistSInc"<<ps.str().c_str();
     else ss<<"AHistS"<<ps.str().c_str();
   ofstream ofileAHistSp(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileAHistP("AHistP", ios::out);
   if(Option_Time == 0) ss<<"AHistT"<<ps.str().c_str();
     else ss<<"AHistT";
   ofstream ofileAHistT(ss.str().c_str(), ios::out);
   ss.str("");
   ofstream ofileTargetP("TargetP", ios::out);
   ss.str("");

   for (iChan = 0; iChan < NHistMax; iChan++)
     {
     y=HistS->GetBinContent(iChan);   ofileHistS  <<y<<endl;
     y=HistSp->GetBinContent(iChan);  ofileHistSp <<y<<endl;
     y=HistP->GetBinContent(iChan);   ofileHistP  <<y<<endl;
     y=HistT->GetBinContent(iChan);   ofileHistT  <<y<<endl;
     y=CHistS->GetBinContent(iChan);  ofileCHistS <<y<<endl;
     y=CHistSp->GetBinContent(iChan); ofileCHistSp<<y<<endl;
     y=CHistP->GetBinContent(iChan);  ofileCHistP <<y<<endl;
     y=CHistT->GetBinContent(iChan);  ofileCHistT <<y<<endl;
     y=BHistS->GetBinContent(iChan);  ofileBHistS <<y<<endl;
     y=BHistSp->GetBinContent(iChan); ofileBHistSp<<y<<endl;
     y=BHistP->GetBinContent(iChan);  ofileBHistP <<y<<endl;
     y=BHistT->GetBinContent(iChan);  ofileBHistT <<y<<endl;
     y=AHistS->GetBinContent(iChan);  ofileAHistS <<y<<endl;
     y=AHistSp->GetBinContent(iChan); ofileAHistSp<<y<<endl;
     y=AHistP->GetBinContent(iChan);  ofileAHistP <<y<<endl;
     y=AHistT->GetBinContent(iChan);  ofileAHistT <<y<<endl;
     y=TargetP->GetBinContent(iChan); ofileTargetP<<y<<endl;
     }

   ofileHistS.close();
   ofileHistSp.close();
   ofileHistP.close();
   ofileHistT.close();
   ofileCHistS.close();
   ofileCHistSp.close();
   ofileCHistP.close();
   ofileCHistT.close();
   ofileBHistS.close();
   ofileBHistSp.close();
   ofileBHistP.close();
   ofileBHistT.close();
   ofileAHistS.close();
   ofileAHistSp.close();
   ofileAHistP.close();
   ofileAHistT.close();
   ofileTargetP.close();

   for (iGasp=0; iGasp<z; iGasp++) fwrite (G[iGasp].elem, 1, (G[iGasp].evSize+1)*sizeof(short), fp);
   fwrite (zero, 1, (grh[0].head[0]*1024/sizeof(short)-recSize)*sizeof(short), fp);
   fclose(fp);

   for (iGasp=0; iGasp<z_indv; iGasp++) fwrite (G_indv[iGasp].elem, 1, (G_indv[iGasp].evSize+1)*sizeof(short), fp_indv);
   fwrite (zero, 1, (grh_indv[0].head[0]*1024/sizeof(short)-recSize_indv)*sizeof(short), fp_indv);
   fclose(fp_indv);

   delete PeaksInCountersBranch;
   delete VectorOfPeaks_Tot;
   delete VectorOfPeaks_R1;
   delete VectorOfPeaks_R2;
   delete VectorOfPeaks_R3;
   delete VectorOfPeaks;
   printf("PeaksInCounters deleted!\n");
   delete PileUpVector;
   printf("PileUpVector deleted!\n");
   delete EmptyVector;
   delete PileUpTarget;
}

