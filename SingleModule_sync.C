#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <limits>


#define CN 64
#define nch_read 15
#define max_diff 50000
#define MAXHITS 5000

TFile *f_light = NULL;
TTree *ts = NULL;
TFile *f_tracks = NULL;
TTree *tr = NULL;

int it_count_offs=0;
int it_count_sync=0;

TH1I *h1 = new TH1I("h1","",1000,1355000,1360000); //Histo to find initial timestamp offset between the two systems
TH1I *h2 = new TH1I("#Delta t","",201,-max_diff,max_diff); //Histo for distribution of time differences

//Light event variables
UInt_t    l_evnr;
ULong64_t l_sync_ts;
ULong64_t l_utime;
ULong64_t l_taisec;
ULong64_t l_tainsec;
char      l_runtag[64];
UInt_t    l_sn;
UInt_t    l_cycle;
Float_t   l_amplitude[CN];
Float_t   l_integral[CN];
Float_t   l_time[CN];
ULong64_t l_wfts[CN];
Float_t   l_tacl1_tot; //collected photons at the detector, ADC units
Float_t   l_tlcm1_tot; //collected photons at the detector, ADC units
Float_t   l_tacl[nch_read]; //collected photons at the detector, ADC units

//Track event variables
Int_t         t_eventID;
Int_t	      t_sync_ts;
Int_t         t_event_start_t; //clock count 10MHz clock
Int_t         t_event_end_t;
Int_t         t_event_duration;
Int_t         t_event_unix_ts;
Int_t         t_event_nhits;
Float_t       t_event_q;
Float_t       t_event_q_raw;
Int_t         t_event_ntracks;
Int_t         t_event_n_ext_trigs;
Float_t		  t_event_hits_x[MAXHITS];       
Float_t		  t_event_hits_y[MAXHITS];
Float_t		  t_event_hits_z[MAXHITS];
Int_t		  t_event_hits_ts[MAXHITS];
Float_t		  t_event_hits_q[MAXHITS];
Int_t		  t_trackID;
Float_t       t_track_start_pos_x;
Float_t       t_track_start_pos_y;
Float_t       t_track_start_pos_z;
Float_t       t_track_start_pos_t;
Float_t       t_track_end_pos_x;
Float_t       t_track_end_pos_y;
Float_t       t_track_end_pos_z;
Float_t       t_track_end_pos_t;
Float_t   	  t_track_length;
Int_t		  t_track_nhits;
Float_t		  t_track_q;
Float_t		  t_track_q_raw;
Float_t		  t_track_theta;
Float_t		  t_track_phi;
Float_t		  t_track_residual_x;
Float_t		  t_track_residual_y;
Float_t		  t_track_residual_z;
Float_t		  t_track_hits_x[MAXHITS];       
Float_t		  t_track_hits_y[MAXHITS];
Float_t		  t_track_hits_z[MAXHITS];
Int_t		  t_track_hits_ts[MAXHITS];
Float_t		  t_track_hits_q[MAXHITS];
Int_t 		  t_trigID;
Int_t		  t_trig_type;
Float_t       t_racl1_tot; //reconstructed photons
Float_t       t_rlcm1_tot; //reconstructed photons
Float_t       t_cgx;
Float_t       t_cgy;
Float_t       t_cgz;

//Further variables
Int_t         l_nentries;
Int_t         t_nentries;
Int_t 	      t_offset;
ULong64_t 	  l_offset;
Int_t 	      t_time_max;
ULong64_t 	  l_time_max;

void sync(Long64_t t_shift,int n_iter=2, bool verbose=false);

void open(const char* file_light,const char* file_tracks,bool allvars=true){
	/**
	* Open light and tracks data file and load to TTrees ts and tr.
	*
	* Input:
	*	- file_light: Light data file with calibration branches (.root) (check Calibration.C by Livio Calivers)
	*	- file_tracks: Charge tracks data file (.root) (check script by Roman Berner)
	*   - allvars: Load all possible variables (bool)
	*/
	
	//Read Light file
	f_light = new TFile(file_light,"READ");
	ts = (TTree*)f_light->Get("rlog");
	l_nentries = ts->GetEntries();
	
    ts->SetBranchAddress("event",&l_evnr);
    ts->SetBranchAddress("utime",&l_utime);
    ts->SetBranchAddress("taisec",&l_taisec);
	if(allvars){
	    ts->SetBranchAddress("tainsec",&l_tainsec);
	    ts->SetBranchAddress("runtag",&l_runtag);
	    ts->SetBranchAddress("sn",&l_sn);
	    ts->SetBranchAddress("cycle",&l_cycle);
	    ts->SetBranchAddress("amplitude",&l_amplitude);
	    ts->SetBranchAddress("integral",&l_integral);
	    ts->SetBranchAddress("time",&l_time);
	    ts->SetBranchAddress("wfts",&l_wfts);
		ts->SetBranchAddress("tacl1_tot",&l_tacl1_tot);
        ts->SetBranchAddress("tlcm1_tot",&l_tlcm1_tot);
		ts->SetBranchAddress("tacl",&l_tacl);
	}	
	//Read Track file
	f_tracks = new TFile(file_tracks,"READ");
	tr = (TTree*)f_tracks->Get("tracks");
	t_nentries = tr->GetEntries();
	
    tr->SetBranchAddress("eventID",&t_eventID);
    tr->SetBranchAddress("event_start_t",&t_event_start_t);
	if(allvars){
	    tr->SetBranchAddress("event_end_t",&t_event_end_t);
	    tr->SetBranchAddress("event_duration",&t_event_duration);
	    tr->SetBranchAddress("event_unix_ts",&t_event_unix_ts);
	    tr->SetBranchAddress("event_nhits",&t_event_nhits);
	    tr->SetBranchAddress("event_q",&t_event_q);
	    tr->SetBranchAddress("event_q_raw",&t_event_q_raw);
	    tr->SetBranchAddress("event_ntracks",&t_event_ntracks);
	    tr->SetBranchAddress("event_n_ext_trigs",&t_event_n_ext_trigs);
		tr->SetBranchAddress("event_hits_x",&t_event_hits_x);
		tr->SetBranchAddress("event_hits_y",&t_event_hits_y);
		tr->SetBranchAddress("event_hits_z",&t_event_hits_z);
		tr->SetBranchAddress("event_hits_ts",&t_event_hits_ts);
		tr->SetBranchAddress("event_hits_q",&t_event_hits_q);
	    tr->SetBranchAddress("trackID",&t_trackID);
	    tr->SetBranchAddress("track_start_pos_x",&t_track_start_pos_x);
	    tr->SetBranchAddress("track_start_pos_y",&t_track_start_pos_y);
	    tr->SetBranchAddress("track_start_pos_z",&t_track_start_pos_z);
	    tr->SetBranchAddress("track_start_pos_t",&t_track_start_pos_t);
	    tr->SetBranchAddress("track_end_pos_x",&t_track_end_pos_x);
	    tr->SetBranchAddress("track_end_pos_y",&t_track_end_pos_y);
	    tr->SetBranchAddress("track_end_pos_z",&t_track_end_pos_z);
	    tr->SetBranchAddress("track_end_pos_t",&t_track_end_pos_t);
	    tr->SetBranchAddress("track_length",&t_track_length);
	    tr->SetBranchAddress("track_nhits",&t_track_nhits);
	    tr->SetBranchAddress("track_q",&t_track_q);
	    tr->SetBranchAddress("track_q_raw",&t_track_q_raw);
	    tr->SetBranchAddress("track_theta",&t_track_theta);
	    tr->SetBranchAddress("track_phi",&t_track_phi);
	    tr->SetBranchAddress("track_residual_x",&t_track_residual_x);
	    tr->SetBranchAddress("track_residual_y",&t_track_residual_y);
	    tr->SetBranchAddress("track_residual_z",&t_track_residual_z);
		tr->SetBranchAddress("track_hits_x",&t_track_hits_x);
		tr->SetBranchAddress("track_hits_y",&t_track_hits_y);
		tr->SetBranchAddress("track_hits_z",&t_track_hits_z);
		tr->SetBranchAddress("track_hits_ts",&t_track_hits_ts);
		tr->SetBranchAddress("track_hits_q",&t_track_hits_q);
	    tr->SetBranchAddress("trigID",&t_trigID);
	    tr->SetBranchAddress("trig_type",&t_trig_type);
	    tr->SetBranchAddress("acl1",&t_racl1_tot);
	    tr->SetBranchAddress("lcm1",&t_rlcm1_tot);
	    tr->SetBranchAddress("cgx",&t_cgx);
	    tr->SetBranchAddress("cgy",&t_cgy);
	    tr->SetBranchAddress("cgz",&t_cgz);
	}
	ts->GetEntry(0);
	tr->GetEntry(0);
	t_offset=t_event_start_t; //offset of first event in us
	l_offset=l_utime; //offset of first event in us
	ts->GetEntry(l_nentries-1);
	tr->GetEntry(t_nentries-1);
	t_time_max = t_event_start_t;
	l_time_max = l_utime;
}


void find_offset(long lower=-100000000, long upper=100000000,int n_iter=0, bool verbose=false){
	/**
	* Plot time stamp difference of each charge event to any light event within a sliding time window and call sync function (if n_iter!=0)
	* If 
	* !!! open(...) has to be run first to load trees !!!
	
	* Input:
	*	- lower: lowest timestamp offset considered (10MHz timestamp, 1e7 = 1second)
	*	- upper: highest timestamp offset considered (10MHz timestamp, 1e7 = 1second)
	*	- n_iter: Number of times the scan is repeated around the detected peak, if n_iter=0 the scan is run a single time and the sync function is not called
	*	- verbose: Enable verbose mode (bool)
	*
	* Output:
	* Plot of ts difference distribution in us
	*/
	
	it_count_offs++;
	h1= new TH1I("h1","",1000,lower,upper);
	//open(file_light,file_tracks,0);
	ts->GetEntry(0);
	tr->GetEntry(0);
	t_offset = t_event_start_t; //offset of first event in us
	l_offset = l_utime; //offset of first event in us
	
	
	Int_t t_evcount=0;
	Int_t evnr_temp=-1;
	
	Long_t current_diff=0; //in us
	
	Long_t first_diff_max=lower;
	Int_t first_counter=0;
	
	Long_t last_diff_max=upper;
	Int_t last_counter=0;
	
	
	for(int i=0;i<400;i++){
		tr->GetEntry(i);
		//printf("fist: %d	",t_eventID);
		if(t_eventID!=evnr_temp){
			evnr_temp = t_eventID;
			//printf("last: %d\n",t_eventID);
			t_evcount++;
			ts->GetEntry(first_counter);
			
      if(verbose==1) printf("i: %d \n",i);

			if(verbose==1) printf("utime: %d \n",(Int_t)((l_utime-l_offset)*625/100000));
			if(verbose==1) printf("t_event_start: %d \n",(t_event_start_t-t_offset));
			
			current_diff=(Long_t)((l_utime-l_offset)*625/100000)-(t_event_start_t-t_offset);
			if(verbose==1) printf("Init: %ld \n",current_diff);
			
			while(current_diff<first_diff_max){
				first_counter++;
				ts->GetEntry(first_counter);
				current_diff=((l_utime-l_offset)*625/100000)-(t_event_start_t-t_offset);
			}
			if(verbose==1) printf("First: %d (%ld) \n",first_counter,current_diff);
			
			last_counter=first_counter;
			while(current_diff<last_diff_max&&last_counter<l_nentries){
				h1->Fill(current_diff);
				last_counter++;
				ts->GetEntry(last_counter);
				current_diff=((l_utime-l_offset)*625/100000)-(t_event_start_t-t_offset);
			}
			if(verbose==1) printf("Last: %d (%ld) \n",last_counter,current_diff);
			if(verbose==1) printf("%d \n",last_counter-first_counter);
			/*
			for(int j=0;j<l_nentries;j++){
				ts->GetEntry(j);
				current_diff=(t_event_start_t-t_offset)/10-((l_utime-l_offset)/1000000*625);
				h1->Fill(current_diff);
			}
			*/
		}
	}
	gStyle->SetOptStat(111111);
	printf("Total charge events: %d\n",t_evcount);
	TCanvas* c1 = new TCanvas("c1","c1");
	c1->SetRightMargin(0.09);
	c1->SetLeftMargin(0.15);
	c1->SetBottomMargin(0.15);
	h1->Draw();
	
	TSpectrum *s = new TSpectrum(1);
	  Int_t nfound = s->Search(h1,2,"",0.50);
	  Double_t *xpeaks;
	  xpeaks = s->GetPositionX();
	  printf("peak: %.0f\n",xpeaks[0]);
	  Double_t range = 0.1*xpeaks[0];
	  if(range<100000) range=100000;
	  //printf("%lld",(Long64_t)(xpeaks[0]));
	  if(it_count_offs<n_iter) find_offset(xpeaks[0]-range,xpeaks[0]+range,n_iter);
	  if(it_count_offs==n_iter) {
		  sync((Long64_t)(xpeaks[0]));
		  it_count_offs=0;
	  }
	  if(n_iter==0) it_count_offs=0; 
}

void sync(Long64_t t_shift_init,int n_iter=0, bool verbose=false){
	/**
	* Produce a root file with synchronized light and track events.
	* !!! open(...) has to be run first to load trees !!!
	* 
	*
	* Input:
	*	- t_shift_init: Time shift between tracks and light timestamp as determined by find_offset(...) (10Mhz timestamp).
	*	- n_iter: Number of iterations the sync is repeated to correct for deviation of t_shift_init to true offset.
	*/
	Long64_t t_shift = t_shift_init;
	it_count_sync++;
	h2->Reset();
	
	Int_t overflow_const = (ULong64_t)std::numeric_limits<Int_t>::max();

	ts->BuildIndex("utime");
	
	ts->GetEntry(0);
    string out_file = "sync_wLightReco_" + string(l_runtag) + ".root";
    TFile* f_out = new TFile(out_file.c_str(),"RECREATE"); //create output file
    TTree * t_out = new TTree("t_out","t_out");
	
    std::string str1 = "ch00/F";
    for(int i=1;i < nch_read;i++){
 	   char buf1[6];
 	   sprintf(buf1,":ch%02d", i);
 	   str1+=buf1;
    }
	
    std::string str2 = "ch00/F";
    for(int i=1;i < CN;i++){
 	   char buf2[6];
 	   sprintf(buf2,":ch%02d", i);
 	   str2+=buf2;
    }
	
    t_out->Branch("l_event",&l_evnr);
    t_out->Branch("l_sync_ts",&l_sync_ts);
    t_out->Branch("l_utime",&l_utime);
    t_out->Branch("l_taisec",&l_taisec);
    t_out->Branch("l_tainsec",&l_tainsec);
    t_out->Branch("l_runtag",&l_runtag,"l_runtag/C");
    t_out->Branch("l_sn",&l_sn);
    t_out->Branch("l_cycle",&l_cycle);
    t_out->Branch("l_amplitude",&l_amplitude,str2.c_str());
    t_out->Branch("l_integral",&l_integral,str2.c_str());
    t_out->Branch("l_time",&l_time,str2.c_str());
    t_out->Branch("l_wfts",&l_wfts,str2.c_str());
	t_out->Branch("l_tacl1_tot",&l_tacl1_tot);
   t_out->Branch("l_tlcm1_tot",&l_tlcm1_tot);
	t_out->Branch("l_tphotons",&l_tacl,str1.c_str());
	
    t_out->Branch("t_eventID",&t_eventID);
    t_out->Branch("t_sync_ts",&t_sync_ts);
    t_out->Branch("t_event_start_t",&t_event_start_t);
    t_out->Branch("t_event_end_t",&t_event_end_t);
    t_out->Branch("t_event_duration",&t_event_duration);
    t_out->Branch("t_event_unix_ts",&t_event_unix_ts);
    t_out->Branch("t_event_nhits",&t_event_nhits);
    t_out->Branch("t_event_q",&t_event_q);
    t_out->Branch("t_event_q_raw",&t_event_q_raw);
    t_out->Branch("t_event_ntracks",&t_event_ntracks);
    t_out->Branch("t_event_n_ext_trigs",&t_event_n_ext_trigs);
	t_out->Branch("t_event_hits_x",&t_event_hits_x,"event_hits_x[t_event_nhits]/F");
	t_out->Branch("t_event_hits_y",&t_event_hits_y,"event_hits_y[t_event_nhits]/F");
	t_out->Branch("t_event_hits_z",&t_event_hits_z,"event_hits_z[t_event_nhits]/F");
	t_out->Branch("t_event_hits_ts",&t_event_hits_ts,"event_hits_ts[t_event_nhits]/F");
	t_out->Branch("t_event_hits_q",&t_event_hits_q,"event_hits_q[t_event_nhits]/F");
    t_out->Branch("t_trackID",&t_trackID);
    t_out->Branch("t_track_start_pos_x",&t_track_start_pos_x);
    t_out->Branch("t_track_start_pos_y",&t_track_start_pos_y);
    t_out->Branch("t_track_start_pos_z",&t_track_start_pos_z);
    t_out->Branch("t_track_start_pos_t",&t_track_start_pos_t);
    t_out->Branch("t_track_end_pos_x",&t_track_end_pos_x);
    t_out->Branch("t_track_end_pos_y",&t_track_end_pos_y);
    t_out->Branch("t_track_end_pos_z",&t_track_end_pos_z);
    t_out->Branch("t_track_end_pos_t",&t_track_end_pos_t);
    t_out->Branch("t_track_length",&t_track_length);
    t_out->Branch("t_track_nhits",&t_track_nhits);
    t_out->Branch("t_track_q",&t_track_q);
    t_out->Branch("t_track_q_raw",&t_track_q_raw);
    t_out->Branch("t_track_theta",&t_track_theta);
    t_out->Branch("t_track_phi",&t_track_phi);
    t_out->Branch("t_track_residual_x",&t_track_residual_x);
    t_out->Branch("t_track_residual_y",&t_track_residual_y);
    t_out->Branch("t_track_residual_z",&t_track_residual_z);
	t_out->Branch("t_track_hits_x",&t_track_hits_x,"track_hits_x[t_track_nhits]/F");
	t_out->Branch("t_track_hits_y",&t_track_hits_y,"track_hits_y[t_track_nhits]/F");
	t_out->Branch("t_track_hits_z",&t_track_hits_z,"track_hits_z[t_track_nhits]/F");
	t_out->Branch("t_track_hits_ts",&t_track_hits_ts,"track_hits_ts[t_track_nhits]/F");
	t_out->Branch("t_track_hits_q",&t_track_hits_q,"track_hits_q[t_track_nhits]/F");
    t_out->Branch("t_trigID",&t_trigID);
    t_out->Branch("t_trig_type",&t_trig_type);
    t_out->Branch("t_acl1",&t_racl1_tot);
    t_out->Branch("t_lcm1",&t_rlcm1_tot);
    t_out->Branch("t_cgx",&t_cgx);
    t_out->Branch("t_cgy",&t_cgy);
    t_out->Branch("t_cgz",&t_cgz);
	
	Int_t        evnr_temp=0;  		//eventID of previous event
	ULong64_t    l_ts_exp;			//expected timestamp for light data
	tr->GetEntry(0);
	Long64_t     t_sync_ts_old=(t_event_start_t-t_offset)/10;;	//charge ts in us from runstart previous event
	Int_t        j=0;				//running variable for light events
	Long64_t	 diff_best=0;
	Long64_t	 diff_temp=0;
	bool		 l_skipped_flag=1;
	Int_t 		 t_evcount=0;			//count # of track events
	Int_t 		 t_neve_early=0;		//count # track events before first light event
	Int_t		 last_overflow_ts=0;	//save last tracks ts reset time to avoid corrupt event to trigger reset
	
	
	
	//printf("l_offset: %llu\n",l_offset);
	//printf("t_offset: %d\n",t_offset);
	//printf("t_shift: %lld\n",t_shift);
	
	for(Int_t i=0;i<t_nentries;i++){
		tr->GetEntry(i);
		if(t_eventID != evnr_temp) {
			t_evcount++;
		    evnr_temp = t_eventID;
			//printf("eventID: %d\n",t_eventID);
			
			//check if charge ts had overflow
			t_sync_ts=(t_event_start_t-t_offset);
			//printf("t_sync_ts: %lld %lld\n",t_sync_ts,t_sync_ts_old);
			
			//if(t_sync_ts-t_sync_ts_old>1000000) continue;
			if(t_sync_ts<t_sync_ts_old) {								//check if tracks ts had overflow
				if((t_event_unix_ts-last_overflow_ts)>100){
					t_shift += overflow_const;
					last_overflow_ts=t_event_unix_ts;
					printf("Event %d : Charge timestamp overflow->Reset shift\n",t_eventID);
					printf("last_overflow_ts: %d\n",last_overflow_ts);
				}else{
					printf("Event %d : Corrupt timestamp->Skip\n",t_eventID);
					t_sync_ts_old=t_sync_ts;
					continue;
				}
			}
			t_sync_ts_old=t_sync_ts;
			
			
			//calculated expected light ts
			if(t_sync_ts + t_shift>0){ 				//check if light run started after charge run
				l_ts_exp = (t_sync_ts + t_shift)/625*100000+l_offset;
			}
			else{ 
				t_neve_early++;
				continue;
			}
			
			//FIND MATCHING LIGHT EVENT
			//printf("l_ts_exp: %lld\n",l_ts_exp);
			j=ts->GetEntryNumberWithBestIndex(l_ts_exp);
			if(verbose) printf("i: %d, eventid: %d, j: %d\n",i,t_eventID,j);
			ts->GetEntry(j);
			//printf("Found ts: %llu\n",l_utime);
			//if(j==l_nentries-1){
			//	printf("Last light event->Break\n");
			//	 break; //check if light data finished
			//}
			
			l_sync_ts = (l_utime-l_offset)*625/100000;
			//printf("l_sync_ts: %lld\n",l_sync_ts);
			diff_temp = ((Long64_t)l_ts_exp-(Long64_t)l_utime)*625/100000;
			diff_best = diff_temp;
			//printf("Diff to exp: %lld \n",diff_temp);
			
			ts->GetEntry(j+1);
			l_sync_ts = (l_utime-l_offset)*625/100000;
			//printf("l_sync_ts: %lld\n",l_sync_ts);
			diff_temp = ((Long64_t)l_ts_exp-(Long64_t)l_utime)*625/100000;
			//printf("Diff to exp: %lld \n",diff_temp);
			if(TMath::Abs(diff_temp)<TMath::Abs(diff_best)){  				//check if matches better
				diff_best=diff_temp;
			} else ts->GetEntry(j);
			l_sync_ts = (l_utime-l_offset)*625/100000;
			
			
			if(TMath::Abs(diff_best)>max_diff) {							//check tolerance
				printf("Event %d : No matching light event found -> skip\n",t_eventID);
				l_skipped_flag=1;
				continue;
			}
			else {
				t_out->Fill();				//write event to output file
				h2->Fill(diff_best);
				l_skipped_flag=0;
			}
		}
		else {
		    evnr_temp=t_eventID;
			if(!l_skipped_flag) t_out->Fill(); //Since multiple tracks per event->Write all tracks belonging to matched event.
			continue;
		}
		
	}	
	gStyle->SetOptStat(111111);
	TCanvas* c2 = new TCanvas("c2","c2");
	c2->SetRightMargin(0.09);
	c2->SetLeftMargin(0.15);
	c2->SetBottomMargin(0.15);
	h2->Draw();
	
	printf("Number of track events: 		%d (%d before first light event)\n",t_evcount,t_neve_early);
	printf("Number of light events: 		%d (%d after last track event)\n",l_nentries,j-l_nentries);
	printf("Number of synchronised events:		%.0f \n",h2->GetEntries());
	
	
    t_out->Write("t_out",TObject::kOverwrite);
    f_out->Close();
	Long64_t shift = (Long64_t)(h2->GetMean());
	printf("%llu\n",shift);
	if(it_count_sync<n_iter) sync(t_shift_init-shift,n_iter);
}


