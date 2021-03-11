#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TEveViewer.h"
#include "TSystem.h"
#include "TGLViewer.h"
#include "TMath.h"
 
#include "TEveViewer.h"
#include "TEvePointSet.h"

#define GRIDNODES 10 // number of nodes for grid in each direction

TCanvas *c=0;
TTree *tr;
char str[128];
TFile* _file0=0;
TEveLine *track;
TEveLine *atrack;
TGeoNavigator *navig=0;
TGeoManager *geom=0;
TGeoNode *node=0;
TGeoNode *tcnt1_node;
TGeoNode *tcnt2_node;
Int_t Nentries;
int gpileup;
Int_t curev=0;
TEveRGBAPalette* pal;
//TGeoManager *geom;
//TEvePointSet* ps=0;
TEveBoxSet* q=0;
TGLAnnotation* ann;
Int_t curnode=0;
Int_t prevnode=0;


//DEFINE LIGHT R/O PARAMETERS

//ORIENTATION: pixelplane: x-y plane
//			   arclight: parallel to y-z plane
// 			   drift: -z direction
//		       origin: middel of pixel_plane

Float_t TPC_size_x = 600;
Float_t TPC_size_y = 1200;
Float_t TPC_size_z = 320;

Float_t TPC_shift_x = 150;
Float_t TPC_shift_y = -150;
Float_t TPC_shift_z = TPC_size_z/2.;

//LCM parameter
Float_t lcm1_size_x = 10;
Float_t lcm1_size_y = 300;
Float_t lcm1_size_z = 280;

Float_t lcm1_shift_x = -150;
Float_t lcm1_shift_y = 0;
Float_t lcm1_shift_z = 140;

//ArCLight parameter
Float_t acl1_size_x = 10;
Float_t acl1_size_y = 300;
Float_t acl1_size_z = 280;

Float_t acl1_shift_x = -150;
Float_t acl1_shift_y = 300;
Float_t acl1_shift_z = 140;


//INITILISE CHARGE EVENT VARIABLES
Int_t           eventID;
Int_t           event_start_t;
Int_t           event_end_t;
Int_t           event_duration;
Int_t           event_unix_ts;
Int_t           event_nhits;
Float_t         event_q;
Float_t         event_q_raw;
Int_t           event_ntracks;
Int_t           event_n_ext_trigs;
Int_t		    trackID;
Float_t         track_start_pos_x;
Float_t         track_start_pos_y;
Float_t         track_start_pos_z;
Float_t         track_start_pos_t;
Float_t         track_end_pos_x;
Float_t         track_end_pos_y;
Float_t         track_end_pos_z;
Float_t         track_end_pos_t;
Float_t   	    track_length;
Int_t		    track_nhits;
Float_t		    track_q;
Float_t		    track_q_raw;
Float_t		    track_theta;
Float_t		    track_phi;
Float_t		    track_residual_x;
Float_t		    track_residual_y;
Float_t		    track_residual_z;
Int_t 		    trigID;
Int_t		    trig_type;

//New variables
Float_t 		lcm1_photons; //collected photons at the detector
Float_t 		acl1_photons; //collected photons at the detector
Float_t 		cgx;
Float_t 		cgy;
Float_t 		cgz;

TBranch 		*blcm1;
TBranch 		*bacl1;
TBranch 		*bcgx;
TBranch 		*bcgy;
TBranch 		*bcgz;



int open(char *fname)
{
	/**
	* Open track data file and load to TTree tr.
	*
	* Input:
	*	- fname: Charge tracks data file (.root) (check script by Roman Berner)
	*/
 	_file0 = new TFile(fname,"UPDATE");
 	if(_file0) { 
	    printf("%s opened.\n",str); 
	    tr=(TTree*)(_file0->Get("tracks"));
	   	//add branches
	   	blcm1 = tr->Branch("lcm1",&lcm1_photons,"lcm1/F");
	   	bacl1 = tr->Branch("acl1",&acl1_photons,"acl1/F");
	   	bcgx = tr->Branch("cgx",&cgx,"cgx/F");
	   	bcgy = tr->Branch("cgy",&cgy,"cgy/F");
	   	bcgz = tr->Branch("cgz",&cgz,"cgz/F");

		// Set branch addresses.
		tr->SetBranchAddress("eventID",&eventID);
		tr->SetBranchAddress("event_start_t",&event_start_t);
		tr->SetBranchAddress("event_end_t",&event_end_t);
		tr->SetBranchAddress("event_duration",&event_duration);
		tr->SetBranchAddress("event_unix_ts",&event_unix_ts);
		tr->SetBranchAddress("event_nhits",&event_nhits);
		tr->SetBranchAddress("event_q",&event_q);
		tr->SetBranchAddress("event_q_raw",&event_q_raw);
		tr->SetBranchAddress("event_ntracks",&event_ntracks);
		tr->SetBranchAddress("event_n_ext_trigs",&event_n_ext_trigs);
		tr->SetBranchAddress("trackID",&trackID);
		tr->SetBranchAddress("track_start_pos_x",&track_start_pos_x);
		tr->SetBranchAddress("track_start_pos_y",&track_start_pos_y);
		tr->SetBranchAddress("track_start_pos_z",&track_start_pos_z);
		tr->SetBranchAddress("track_start_pos_t",&track_start_pos_t);
		tr->SetBranchAddress("track_end_pos_x",&track_end_pos_x);
		tr->SetBranchAddress("track_end_pos_y",&track_end_pos_y);
		tr->SetBranchAddress("track_end_pos_z",&track_end_pos_z);
		tr->SetBranchAddress("track_end_pos_t",&track_end_pos_t);
		tr->SetBranchAddress("track_length",&track_length);
		tr->SetBranchAddress("track_nhits",&track_nhits);
		tr->SetBranchAddress("track_q",&track_q);
		tr->SetBranchAddress("track_q_raw",&track_q_raw);
		tr->SetBranchAddress("track_theta",&track_theta);
		tr->SetBranchAddress("track_phi",&track_phi);
		tr->SetBranchAddress("track_residual_x",&track_residual_x);
		tr->SetBranchAddress("track_residual_y",&track_residual_y);
		tr->SetBranchAddress("track_residual_z",&track_residual_z);
		tr->SetBranchAddress("trigID",&trigID);
		tr->SetBranchAddress("trig_type",&trig_type);
 
	    return tr->GetEntries();
 	}
 	printf("Error: %s can't be opened.\n",str); return 0; 
}

void SingleModule3DLightReco(char *fname, int pileup=0)
{
	/**
	* Initialise Light reco, build geometry, load first track
	
	* Input:
		- fname: Charge tracks data file (.root) (check script by Roman Berner)
		- pileupe: Pile up tracks in viewer (=1) 
	*/
    gpileup=pileup;
    Nentries=open(fname);

    track=new TEveLine(2);
    track->SetNextPoint(0,0,0);
    track->SetNextPoint(10,10,10);

	//--- Definition of a simple geometry
   	gSystem->Load("libGeom");
   	geom = new TGeoManager("SingleCube",      "SingleCube");
   	Int_t i;
   	Float_t z;
   	//--- define some materials
   	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
   	TGeoMaterial *matAr = new TGeoMaterial("Ar", 26.98,26,1.4);
    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Aluminium",2, matAl);
    TGeoMedium *Ar = new TGeoMedium("Argon",3, matAr);
    //TGeoMedium *LAr = new TGeoMedium("Liquid Argon",146,146,0,1,8,90,-1,-1,0.5,-1);

	//--- make the top container volume
	TGeoVolume *top = geom->MakeBox("TOP", Vacuum, 1000., 1000., 1000.);
	geom->SetTopVolume(top);

	// Make the elementary assembly of the whole structure
	TGeoVolume *tpc = new TGeoVolumeAssembly("TPC");
	i=2;

	TGeoVolume *argon = geom->MakeBox("ARGON", Ar,TPC_size_x/2. ,TPC_size_y/2. ,TPC_size_z/2.);
	argon->SetTransparency(70);
	argon->SetVisibility(kTRUE);
	argon->SetLineColor(kBlue);

	tpc->AddNode(argon,i,new TGeoTranslation(TPC_shift_x,TPC_shift_y,TPC_shift_z));  i++;

	TGeoVolume *pixel_board=geom->MakeBox("PixelBoard", Al, 150,150,0.1);
	pixel_board->SetTransparency(70);
	pixel_board->SetVisibility(kTRUE);
	pixel_board->SetLineColor(kBlue);

	tpc->AddNode(pixel_board,i,new TGeoTranslation(0,0,0)); i++;


	TGeoVolume *cath = geom->MakeBox("CATHODE", Al, TPC_size_x/2. ,TPC_size_y/2. ,0.1);
	cath->SetTransparency(90);
	cath->SetLineColor(kYellow);
	tpc->AddNode(cath, i, new TGeoTranslation(TPC_shift_x,TPC_shift_y,TPC_size_z)); i++;

	TGeoVolume *acl1=geom->MakeBox("ACL1 module", Al, acl1_size_x/2., acl1_size_y/2., acl1_size_z/2.);
	acl1->SetTransparency(50);
	acl1->SetVisibility(kTRUE);
	acl1->SetLineColor(kGreen);

	tpc->AddNode(acl1,i,new TGeoTranslation(acl1_shift_x-acl1_size_x/2.,acl1_shift_y,acl1_shift_z));
	i++;

	TGeoVolume *lcm1=geom->MakeBox("LCM1 module", Al, lcm1_size_x/2., lcm1_size_y/2., lcm1_size_z/2.);
	lcm1->SetTransparency(50);
	lcm1->SetVisibility(kTRUE);
	lcm1->SetLineColor(kGreen);

	tpc->AddNode(lcm1,i,new TGeoTranslation(lcm1_shift_x-lcm1_size_x/2.,lcm1_shift_y,lcm1_shift_z));
	i++;


	top->AddNode(tpc, 0, new TGeoTranslation(0,0,0));

	/*
	TGeoVolume *lcm=geom->MakeBox("LCM module", Al, lcm_size_x/2., lcm_size_y/2., lcm_size_z/2.);
	lcm->SetTransparency(50);
	lcm->SetVisibility(kTRUE);
	lcm->SetLineColor(kGreen);
	*/

	//light->AddNode(lcm,i,new TGeoTranslation(lcm_pitch,lcm_y_offset-0.5-lcm_size_y, lcm_z_offset));   //third LCM module not in use


	//--- close the geometry
	geom->CloseGeometry();

	TEveManager::Create();

	TGeoNode* node = gGeoManager->GetTopNode();
	TEveGeoTopNode* en = new TEveGeoTopNode(gGeoManager, node);
	en->SetVisLevel(4);
	en->GetNode()->GetVolume()->SetVisibility(kFALSE);

	gEve->AddGlobalElement(en);

    Double_t cc[3]={0,0,0};
    gEve->GetDefaultGLViewer()->SetPerspectiveCamera(TGLViewer::kCameraPerspXOY, 1, 0,cc,0,900);
    gEve->GetDefaultGLViewer()->SetClearColor(kWhite);
 
 	gEve->Redraw3D(kTRUE);

   	en->ExpandIntoListTreesRecursively();
   	TGLViewer* v = gEve->GetDefaultGLViewer();
   	ann = new TGLAnnotation(v, "Event", 0.1, 0.9);
   	ann->SetTextSize(0.03);// % of window diagonal

   	gEve->AddElement(track);
/*
	fViewer0=gEve->GetDefaultGLViewer();
	TGLOverlayButton *but1 = new TGLOverlayButton(fViewer0, "fw", 10.0, 10.0, 65.0, 26.0); 
	but1->Connect("Clicked()", 0, 0, "PrevEvent()");
	TGLOverlayButton *but2 = new TGLOverlayButton(fViewer0, "bw", 85.0, 10.0, 115.0, 26.0); 
	but2->Connect("Clicked()", 0, 0, "NextEvent()");
	TEveText* t = new TEveText("DADA");
	t->PtrMainTrans()->RotateLF(1, 3, TMath::PiOver2());
	t->SetMainColor(kOrange-2);
	t->SetFontSize(64);
	t->SetFontMode(TGLFont::kExtrude);
	t->SetLighting(kTRUE);
	gEve->AddElement(t);
*/


	TObjArray *navigators = geom->GetListOfNavigators(); 
	navig=geom->GetCurrentNavigator(); 

	curev=0;
	lcm1_photons=0;
	acl1_photons=0;
	AddTrack(curev);
}

void NextEvent()
{ 
	/**
	* Load next event
	*/
	
	curev++;
 	AddTrack(curev);
}
void PrevEvent()
{ 
	/**
	* Load previous event
	*/
 	if(curev==0) return;
 	curev--;
 	AddTrack(curev);
}

void fullrun(char *fname){
	SingleModule3DLightReco(fname,1);
	for (int i=1;i<Nentries;i++){
		NextEvent();
	}
	tr->Write(0,TObject::kOverwrite);
	_file0->Close();
}


void AddTrack(int evt)
{
	printf("Processing entry %d...",evt);
	tr->GetEntry(evt);
	printf("Event: %d, Track: %d Event time %llu\n",eventID,trackID,event_unix_ts);

 
    if(track>0 && gpileup==0)
    {
    gEve->PreDeleteElement(track);
    delete track;
    track=0;
    }  



	node=0;

	lcm1_photons=0;
	acl1_photons=0;

	//if(1){
	if(track_nhits>10 && track_length>10){ 						//Apply cut
		Float_t cgx_temp=(track_end_pos_x+track_start_pos_x)/2.;
		Float_t cgy_temp=(track_end_pos_y+track_start_pos_y)/2.;
		Float_t cgz_temp=(track_end_pos_z+track_start_pos_z)/2.;

		Float_t vx=(track_end_pos_x-track_start_pos_x)/track_length;
		Float_t vy=(track_end_pos_y-track_start_pos_y)/track_length;
		Float_t vz=(track_end_pos_z-track_start_pos_z)/track_length;

		track=new TEveLine(2);
		track->Reset(2);
		//Extended Tracks
		//track->SetNextPoint(cgx-(cgz-150)*vx,cgy-(cgz-150)*vy,cgz-(cgz-150)*vz);
		//track->SetNextPoint(cgx-(cgz+150)*vx,cgy-(cgz+150)*vy,cgz-(cgz+150)*vz);
		//Recorded Tracks
		track->SetNextPoint(track_start_pos_x,track_start_pos_y,track_start_pos_z);
		track->SetNextPoint(track_end_pos_x,track_end_pos_y,track_end_pos_z);
		track->SetLineColor(kRed);
		gEve->AddElement(track);


		//now go step by step and highlight track elements in argon

		Float_t curx,cury,curz;
		Int_t istep=0;
		Float_t step=1.0; //1mm step
		curx=cgx_temp; cury=cgy_temp; curz=cgz_temp;
		geom->FindNode(curx,cury,curz); 
		node=geom->GetCurrentNode();
		curnode=node->GetNumber();
		prevnode=curnode;
		//printf("Initial node: %s\n",node->GetName());
		atrack=new TEveLine(2);
		atrack->SetPoint(0,curx,cury,curz);
		//first up from the cg (z to positive)
		while(curz>0 && curz < TPC_size_z && sqrt((curx-TPC_shift_x)*(curx-TPC_shift_x))< TPC_size_x/2. && sqrt((cury-TPC_shift_y)*(cury-TPC_shift_y))< TPC_size_y/2. ) //Go along track as long as within TPC boundary
		{
			istep++;
			curx=cgx_temp+vx*istep*step;
			cury=cgy_temp+vy*istep*step;
			curz=cgz_temp+vz*istep*step;
			geom->FindNode(curx,cury,curz); 
			node=geom->GetCurrentNode();
			curnode=node->GetNumber();
			//   printf("Step %d curnode: %s, curnode number%d, curxyz: %f %f %f\n",istep,node->GetName(), curnode,curx,cury,curz);
			//   if(curnode==prevnode && curnode==1) {printf("Adding point..\n"); atrack->SetPoint(atrack->GetN(),curx,cury,curz);}
			if(curnode!=prevnode && curnode==2) { atrack=new TEveLine(2); atrack->SetPoint(0,curx,cury,curz);}
			if(curnode==prevnode && curnode==2) TracePhotons(curx,cury,curz);
			if(curnode!=prevnode && prevnode==2) {atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);} //last point in a track
			prevnode=curnode;
		}
		if(prevnode==2) {
			printf("END: Adding end point, adding track to gEve..\n"); 
			atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);
		} //last point in a track

		cgx=curx;
		cgy=cury;
		cgy=cury;

		//now down from the cg (z to negative)
		curx=cgx_temp; cury=cgy_temp; curz=cgz_temp;
		geom->FindNode(curx,cury,curz); 
		node=geom->GetCurrentNode();
		curnode=node->GetNumber();
		prevnode=curnode;
		//printf("Initial node: %s\n",node->GetName());
		atrack=new TEveLine(2);
		atrack->SetPoint(0,curx,cury,curz);
		istep=0; 
		while(curz>0 && curz < TPC_size_z && sqrt(curx*curx)< TPC_size_x/2. && sqrt(cury*cury)< TPC_size_y/2. )
		{
			istep++;
			curx=cgx_temp-vx*istep*step;
			cury=cgy_temp-vy*istep*step;
			curz=cgz_temp-vz*istep*step;
			geom->FindNode(curx,cury,curz); 
			node=geom->GetCurrentNode();
			curnode=node->GetNumber();
			// printf("Step %d curnode: %s, curnode number%d, curxyz: %f %f %f\n",istep,node->GetName(), curnode,curx,cury,curz);
			//   if(curnode==prevnode && curnode==1) {printf("Adding point..\n"); atrack->SetPoint(atrack->GetN(),curx,cury,curz);}
			if(curnode!=prevnode && curnode==2) { atrack=new TEveLine(2); atrack->SetPoint(0,curx,cury,curz);}
			if(curnode==prevnode && curnode==2) TracePhotons(curx,cury,curz);
			if(curnode!=prevnode && prevnode==2) { atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);} //last point in a track
			prevnode=curnode;
		}
		if(prevnode==2) {
			printf("END: Adding end point, adding track to gEve..\n"); 
			atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);
		} //last point in a track

		cgx=(cgx+curx)/2.;
		cgy=(cgy+cury)/2.;
		cgz=(cgz+curz)/2.;


   }
   
   printf("lcm1 detected %f virtual photons in this event\n",lcm1_photons);


   blcm1->Fill();
   bacl1->Fill();
   bcgx->Fill();
   bcgy->Fill();
   bcgz->Fill();

   sprintf(str,"Event %d",eventID);
   ann->SetText(str);
   gEve->Redraw3D(kFALSE);

}

TracePhotons(Float_t px, Float_t py, Float_t pz, Int_t rays=0)
{
	Float_t ldx,ldy,ldz; // coordinates of grid nodes at the light detectors
	Float_t ldd; //distance from emitting point to grid node
	Float_t  domega; //solid angle
	Double_t v0[3];
	Double_t v1[3];
	Double_t rv[3]; //ray vector
	Double_t *nv; //normal vector
	Double_t cosine=1; //cosine of the ray-to-normal angle
	Double_t cosiney=1; //cosine of the ray-to-normal angle
	Double_t cosinez=1; //cosine of the ray-to-normal angle
	TEveLine * rtrace;
	//TGeoRotation *rotls = new TGeoRotation("rotls",0,0,0);
	//printf("Photons emitted at %f %f %f\n",px,py,pz);

	for(Float_t iy=0; iy<GRIDNODES; iy++) for(Float_t iz=0; iz<GRIDNODES; iz++){
		// lcm1
		//calculate coordinates and distance
		//coordinates in a light detector plane 
		ldx= lcm1_shift_x;
		ldy= lcm1_shift_y-lcm1_size_y/2.0 + (iy+0.5)/GRIDNODES*lcm1_size_y;
		ldz= lcm1_shift_z-lcm1_size_z/2.0 + (iz+0.5)/GRIDNODES*lcm1_size_z;

		if(rays==1){
			rtrace=new TEveLine(2);
			rtrace->SetPoint(0,px,py,pz);
			rtrace->SetPoint(1,ldx,ldy,ldz);
			gEve->AddElement(rtrace);
			}
		rv[0]=px-ldx; rv[1]=py-ldy; rv[2]=pz-ldz;
		ldd=sqrt( rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2] );
		//check if the emission point is actially seen from the detector
			geom->SetCurrentDirection(-rv[0],-rv[1],-rv[2]); 
			geom->FindNode(px,py,pz);
			node=geom->FindNextBoundary();
			nv=geom->FindNormal(kTRUE);
			
			//Calculate solid angle
			cosine=-(nv[0]*rv[0]+nv[1]*rv[1]+nv[2]*rv[2])/(ldd*1.0); //normal vector has unit length
			//cosiney=1./sqrt(rv[1]*rv[1]/(rv[0]*rv[0])+1.);
			//cosinez=1./sqrt(rv[2]*rv[2]/(rv[0]*rv[0])+1.);    
			domega=cosine*4*TMath::ASin((lcm1_size_y/GRIDNODES)*(lcm1_size_z/GRIDNODES)/sqrt((4.0*ldd*ldd+(lcm1_size_y/GRIDNODES)*(lcm1_size_y/GRIDNODES))*(4.0*ldd*ldd+(lcm1_size_z/GRIDNODES)*(lcm1_size_z/GRIDNODES)))); // * fabs(ldy-py) /ldd;
			//domega=4*TMath::ASin(cosiney*cosinez*(lcm1_size_y/GRIDNODES)*(lcm1_size_z/GRIDNODES)/sqrt((4.0*ldd*ldd+cosiney*cosiney*(lcm1_size_y/GRIDNODES)*(lcm1_size_y/GRIDNODES))*(4.0*ldd*ldd+cosinez*cosinez*(lcm1_size_z/GRIDNODES)*(lcm1_size_z/GRIDNODES))));

			lcm1_photons=lcm1_photons+1.*domega/(4*3.14152);
			if(rays==2){
				rtrace=new TEveLine(2);
				rtrace->SetPoint(0,px,py,pz);
				rtrace->SetPoint(1,ldx,ldy,ldz);
				gEve->AddElement(rtrace);
				printf("domega %f cosine %f\n",domega,cosine);
			}
	}

	for(Float_t iy=0; iy<GRIDNODES; iy++) for(Float_t iz=0; iz<GRIDNODES; iz++){
		// ACL1
		//calculate coordinates and distance
		//coordinates in a light detector plane 
		ldx= acl1_shift_x;
		ldy= acl1_shift_y-acl1_size_y/2.0 + (iy+0.5)/GRIDNODES*acl1_size_y;
		ldz= acl1_shift_z-acl1_size_z/2.0 + (iz+0.5)/GRIDNODES*acl1_size_z;

		if(rays==1){
			rtrace=new TEveLine(2);
			rtrace->SetPoint(0,px,py,pz);
			rtrace->SetPoint(1,ldx,ldy,ldz);
			gEve->AddElement(rtrace);
		}
		rv[0]=px-ldx; rv[1]=py-ldy; rv[2]=pz-ldz;
		ldd=sqrt( rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2] );
		//check if the emission point is actially seen from the detector
		geom->SetCurrentDirection(-rv[0],-rv[1],-rv[2]);
		geom->FindNode(px,py,pz);
		node=geom->FindNextBoundary();
		nv=geom->FindNormal(kTRUE);
		
		//Calculate solid angle
		cosine=-(nv[0]*rv[0]+nv[1]*rv[1]+nv[2]*rv[2])/(ldd*1.0); //normal vector has unit length
		//cosiney=1./sqrt(rv[1]*rv[1]/(rv[0]*rv[0])+1.);
		//cosinez=1./sqrt(rv[2]*rv[2]/(rv[0]*rv[0])+1.);
		domega=cosine*4*TMath::ASin((acl1_size_y/GRIDNODES)*(acl1_size_z/GRIDNODES)/sqrt((4.0*ldd*ldd+(acl1_size_y/GRIDNODES)*(acl1_size_y/GRIDNODES))*(4.0*ldd*ldd+(acl1_size_z/GRIDNODES)*(acl1_size_z/GRIDNODES))));
		//domega=4*TMath::ASin(cosiney*cosinez*(acl1_size_y/GRIDNODES)*(acl1_size_z/GRIDNODES)/sqrt((4.0*ldd*ldd+cosiney*cosiney*(acl1_size_y/GRIDNODES)*(acl1_size_y/GRIDNODES))*(4.0*ldd*ldd+cosinez*cosinez*(acl1_size_z/GRIDNODES)*(acl1_size_z/GRIDNODES))));
		acl1_photons=acl1_photons+1.*domega/(4*3.14152);
		if(rays==2){
			rtrace=new TEveLine(2);
				rtrace->SetPoint(0,px,py,pz);
				rtrace->SetPoint(1,ldx,ldy,ldz);
				gEve->AddElement(rtrace);
				printf("domega %f cosine %f\n",domega,cosine);
		}
	}
}

 
