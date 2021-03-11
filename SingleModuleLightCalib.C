//#include <unistd>
#include <iostream>
#include "TClassTable.h"

#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <cmath>
#include <vector>


//#include <stdio>
//#include <stdint.h>
//#include <string>
//#include <stdlib>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"
#include "TROOT.h"

#define CN 64
#define nch_read 15


TFile *f_light = NULL;
TTree *ts = NULL;

//g492_542
Float_t p[nch_read]={0,0,-12.3,-7.27,-10.10,-9.72,-3.45,-5.47,0,-3.94,-7.89,0.0,-6.29,-6.53,-1.04};
Float_t g[nch_read]={1e9,1e9,111.18,117.36,109.12,115.87,107.18,118.84,1e9,119.12,120.19,113.50,106.27,115.92,121.04};

 void Calib(const char* file_light){
	 /**
	 Update light data file with branches containing number of p.e.
	 
	 INPUT
	 	- file_light: Light data file (.root) created by ADCViewer
	 */
	 
	gROOT->Reset();
	f_light = new TFile(file_light,"UPDATE");
	ts = (TTree*)f_light->Get("rlog");
	
	//Declaration of leaves types
    UInt_t    arc_evnr;
    ULong64_t arc_evtime;
    char      arc_runtag[64];
    UInt_t    arc_sn;
    UInt_t    arc_cycle;
    Float_t   amplitude[CN];
    Float_t   integral[CN];
    Float_t   l_time[CN];
    ULong64_t wfts[CN];

	//Declare new variables
    Float_t TrueACL1_tot; //Total collected photons at the ACL detector, ADC units
    Float_t TrueLCM1_tot; //Total collected photons at the LCM detector, ADC units
    Float_t TrueACL[nch_read]; //Array of collected photons at each SiPM, ADC units
	

   	TBranch *bacl1;
   	TBranch *bacl2;
   	TBranch *bacl3;

   

   //add branches for true values
   bacl1 = ts->Branch("tacl1_tot",&TrueACL1_tot,"tacl1_tot/F");
   bacl3 = ts->Branch("tlcm1_tot",&TrueLCM1_tot,"tlcm1_tot/F");
   std::string str = "ch00/F";
   for(int i=1;i < nch_read;i++){
	   char buf[6];
	   sprintf(buf,":ch%02d", i);
	   str+=buf;
   }
   bacl2 = ts->Branch("tacl",&TrueACL,str.c_str());


   // Set branch addresses.
   ts->SetBranchAddress("event",&arc_evnr);
   ts->SetBranchAddress("utime",&arc_evtime);
   ts->SetBranchAddress("sn",&arc_sn);
   ts->SetBranchAddress("cycle",&arc_cycle);
   ts->SetBranchAddress("runtag",&arc_runtag);
   ts->SetBranchAddress("amplitude",&amplitude);
   ts->SetBranchAddress("integral",&integral);
   ts->SetBranchAddress("time",&l_time);
   ts->SetBranchAddress("wfts",&wfts);

   for(int en=0; en<ts->GetEntries(); en++)
   {
     ts->GetEntry(en);
	 TrueACL[0]=0;
	 TrueACL[1]=0;
	 for(int i=2;i<nch_read;i++){
		 if((i-1)%7==0) TrueACL[i]=0;
		 else TrueACL[i]=(integral[i]+p[i])/g[i];
	 }
	 TrueACL1_tot=TrueACL[2]+TrueACL[3]+TrueACL[4]+TrueACL[5]+TrueACL[6]+TrueACL[7];
   	 TrueLCM1_tot=TrueACL[9]+TrueACL[10]+TrueACL[11]+TrueACL[12]+TrueACL[13]+TrueACL[14];
     printf("ACL1  photons true:%3.1f \n",TrueACL1_tot);
         
     bacl1->Fill();
     bacl2->Fill();
     bacl3->Fill();


   } 
  ts->Write("rlog",TObject::kOverwrite);
  f_light->Close();

 }
