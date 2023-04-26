#ifndef ROOT_TEveManager
#define ROOT_TEveManager

#include <TEveTrack.h>
#include <TEveVSDStructs.h>
#include <TEveManager.h>
#include <TEveViewer.h>
#include <TSystem.h>
#include <TGLViewer.h>
#include <TMath.h>
 
#include <TEvePointSet.h>

#include <vector>
#include <cmath>
#include <string.h>
#include <time.h>


#define nsipm 64 // number of light detection modules for each type

TCanvas *c=0;
TTree *tr;
TTree *evtr;
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

int cut_count=0;
clock_t t0, t1, t2, t3;


//DEFINE LIGHT R/O PARAMETERS

//ORIENTATION: pixelplane: x-y plane
//			   arclight: parallel to y-z plane
// 			   drift: -z direction
//		       origin: middel of pixel_plane

Float_t TPC_size_x = 300;
Float_t TPC_size_y = 300;
Float_t TPC_size_z = 300;

Float_t TPC_shift_x = 0;
Float_t TPC_shift_y = 0;
Float_t TPC_shift_z = TPC_size_z/2.;

//LCM parameter
Float_t sipm_size_x = 6;
Float_t sipm_size_y = 6;
Float_t sipm_size_z = 0.1;

Float_t sipm_shift_x = 0;
Float_t sipm_shift_y = 0;
Float_t sipm_shift_z = 0;

Int_t nsipm_x = 8;
Int_t nsipm_y = 8;

Float_t sipm_space_x = TPC_size_x / nsipm_x;
Float_t sipm_space_y = TPC_size_y / nsipm_y;

Float_t lyield = 22538; // ph/MeV at 0.5 kV/cm




//INITILISE CHARGE EVENT VARIABLES
Int_t		eventID;
//Int_t		event_ntracks;
//Int_t		trackID;
//Int_t		track_eventID;
Float_t         track_start_pos_x;
Float_t         track_start_pos_y;
Float_t         track_start_pos_z;
//Float_t         track_start_pos_t;
Float_t         track_end_pos_x;
Float_t         track_end_pos_y;
Float_t         track_end_pos_z;
//Float_t         track_end_pos_t;
Float_t			track_dx;
Float_t			track_dy;
Float_t			track_dz;
Float_t   	    track_length;
//Int_t		    track_nhits;
Float_t		    track_q;
//Float_t		    track_q_raw;
//Float_t		    track_theta;
//Float_t		    track_phi;
//Int_t 		    trigID;
//Int_t		    trig_type;

//New variables
Float_t 		sipm_vis[nsipm]={}; //collected visibility at the detector
Float_t 		sipm_npe[nsipm]={}; //collected visibility at the detector
Float_t 		cgx;
Float_t 		cgy;
Float_t 		cgz;
Int_t			lightreco_flag;
Float_t			dqdx;

TBranch 		*bsipm;
TBranch			*bnpe;
TBranch 		*bcgx;
TBranch 		*bcgy;
TBranch 		*bcgz;
TBranch 		*blrf;



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
	    printf("opened.\n",str); 
		//evtr=(TTree*)(_file0->Get("events"));
	        tr=(TTree*)(_file0->Get("tracks"));
	   	//add branches
	   	bsipm = tr->Branch("sipm_vis",&sipm_vis,"sipm_vis[64]/F");
		bnpe = tr->Branch("sipm_npe",&sipm_npe,"sipm_npe[64]/F");
	   	bcgx = tr->Branch("cgx",&cgx,"cgx/F");
	   	bcgy = tr->Branch("cgy",&cgy,"cgy/F");
	   	bcgz = tr->Branch("cgz",&cgz,"cgz/F");
	   	blrf = tr->Branch("lightreco_flag",&lightreco_flag,"lightreco_flag/I");		

		// Set branch addresses.
		tr->SetBranchAddress("eventID",&eventID);
//		evtr->SetBranchAddress("event_ntracks",&event_ntracks);
//		tr->SetBranchAddress("trackID",&trackID);
//		tr->SetBranchAddress("track_eventID",&track_eventID);
		tr->SetBranchAddress("track_start_pos_x",&track_start_pos_x);
		tr->SetBranchAddress("track_start_pos_y",&track_start_pos_y);
		tr->SetBranchAddress("track_start_pos_z",&track_start_pos_z);
		tr->SetBranchAddress("track_end_pos_x",&track_end_pos_x);
		tr->SetBranchAddress("track_end_pos_y",&track_end_pos_y);
		tr->SetBranchAddress("track_end_pos_z",&track_end_pos_z);
//		tr->SetBranchAddress("track_length",&track_length);
//		tr->SetBranchAddress("track_dx",&track_dx);
//		tr->SetBranchAddress("track_dy",&track_dy);
//		tr->SetBranchAddress("track_dz",&track_dz);
		tr->SetBranchAddress("track_q",&track_q);
//		tr->SetBranchAddress("track_nhits",&track_nhits);
//		tr->SetBranchAddress("track_theta",&track_theta);
//		tr->SetBranchAddress("track_phi",&track_phi);

		//evtr->BuildIndex("eventID");
 
	    return tr->GetEntries();
 	}
 	printf("Error: %s can't be opened.\n",str); return 0; 
}

void SingleModule3DLightReco(char *fname, int pileup=1)
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
	//TGeoVolume *top = geom->MakeBox("TOP", Vacuum, 1000., 1000., 1000.);
	//geom->SetTopVolume(top);

	// Make the elementary assembly of the whole structure
	TGeoVolume *tpc = new TGeoVolumeAssembly("TPC");
	i=2;

	TGeoVolume *argon = geom->MakeBox("ARGON", Ar,TPC_size_x/2. ,TPC_size_y/2. ,TPC_size_z/2.);
	argon->SetTransparency(70);
	argon->SetVisibility(kTRUE);
	argon->SetLineColor(kBlue);

	//tpc->AddNode(argon,i,new TGeoTranslation(TPC_shift_x,TPC_shift_y,TPC_shift_z));  i++;

	TGeoVolume *pixel_board=geom->MakeBox("PixelBoard", Al, 150,150,1);
	pixel_board->SetTransparency(70);
	pixel_board->SetVisibility(kTRUE);
	pixel_board->SetLineColor(kBlue);


	tpc->AddNode(pixel_board,i,new TGeoTranslation(0,0,-0.5)); i++;
	
	TGeoVolume *cath = geom->MakeBox("CATHODE", Al, TPC_size_x/2. ,TPC_size_y/2. ,0.1);
	cath->SetTransparency(90);
	cath->SetLineColor(kYellow);
	tpc->AddNode(cath, i, new TGeoTranslation(TPC_shift_x,TPC_shift_y,TPC_size_z)); i++;

	TGeoVolume *sipm=geom->MakeBox("SIPM module", Al, sipm_size_x/2., sipm_size_y/2., sipm_size_z/2.);
	sipm->SetTransparency(50);
	sipm->SetVisibility(kTRUE);
	sipm->SetLineColor(kGreen);

	for(int j=0;j<nsipm_x;j++){
        for(int l=0;l<nsipm_y;l++){
		    tpc->AddNode(sipm,i,new TGeoTranslation(sipm_shift_x-TPC_size_x/2.+(j+0.5)*sipm_space_x,sipm_shift_y-TPC_size_y/2.+(l+0.5)*sipm_space_y ,sipm_shift_z));
		    i++;
        }
	}


	//top->AddNode(tpc, 0, new TGeoTranslation(0,0,0));
	geom->SetTopVolume(tpc);
	//geom->AddNode(tpc, 0, new TGeoTranslation(0,0,0));
	

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
	memset(sipm_vis, 0, sizeof(sipm_vis));
	memset(sipm_npe, 0, sizeof(sipm_npe));
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
	printf("%d Events found\n",Nentries);
	//Nentries = 1;
	for (int i=0;i<Nentries;i++){
		//printf("Processing Entry %d of %d\r",i,Nentries);
		NextEvent();
	}
	tr->Write(0,TObject::kOverwrite);
	//evtr->Write(0,TObject::kOverwrite);
	_file0->Close();
}


void AddTrack(int evt)
{
	
	tr->GetEntry(evt);
	//evtr->GetEntryWithIndex(track_eventID);
	//printf("Event: %d, Track: %d Event time %llu\n",eventID,trackID,event_unix_ts);
	//while(track_dy<1000){
	//	evt++;
	//	tr->GetEntry(evt);
	//}
	lightreco_flag=0;
 
    if(track>0 && gpileup==0)
    {
    gEve->PreDeleteElement(track);
    delete track;
    track=0;
    }  



	node=0;

	memset(sipm_vis, 0, sizeof(sipm_vis));
	memset(sipm_npe, 0, sizeof(sipm_npe));


	if(1){
	//if(track_dy>1000 && event_ntracks<5){ 
	//	printf("tr: %d ev: %d\n",track_eventID,eventID);
	//	printf("%d\n",++cut_count);
	//}					//Apply cut
	//if(0){
		printf("Processing entry %d of %d\n",evt,Nentries);
		lightreco_flag=1;				
		cgx = (track_end_pos_x+track_start_pos_x)/2.;
		cgy = (track_end_pos_y+track_start_pos_y)/2.;
		cgz = (track_end_pos_z+track_start_pos_z)/2.;

		track_dx=track_end_pos_x-track_start_pos_x;
		track_dy=track_end_pos_y-track_start_pos_y;
		track_dz=track_end_pos_z-track_start_pos_z;

		track_length = sqrt(track_dx*track_dx+track_dy*track_dy+track_dz*track_dz);
		dqdx = track_q/track_length;
		//printf("tr_len: %f tr_dqdx: %f",track_length,dqdx);

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
		Float_t istep=0.5;
		Float_t step=1.0; //1mm step
		curx=track_start_pos_x+vx*istep*step;
		cury=track_start_pos_y+vy*istep*step;
		curz=track_start_pos_z+vz*istep*step;
		geom->FindNode(curx,cury,curz); 
		node=geom->GetCurrentNode();
		curnode=node->GetNumber();
		prevnode=curnode;
		//printf("Initial node: %s\n",node->GetName());
		atrack=new TEveLine(2);
		atrack->SetPoint(0,curx,cury,curz);
		//first up from the cg (z to positive)
		do{
			geom->FindNode(curx,cury,curz); 
			node=geom->GetCurrentNode();
			curnode=node->GetNumber();
			//   printf("Step %d curnode: %s, curnode number%d, curxyz: %f %f %f\n",istep,node->GetName(), curnode,curx,cury,curz);
			//   if(curnode==prevnode && curnode==1) {printf("Adding point..\n"); atrack->SetPoint(atrack->GetN(),curx,cury,curz);}
			if(curnode!=prevnode && curnode==2) { atrack=new TEveLine(2); atrack->SetPoint(0,curx,cury,curz);}
			if(istep==1) {TracePhotons(curx,cury,curz,0); printf("?");}
			if(istep!=1) TracePhotons(curx,cury,curz,0);
			//if(curnode==prevnode && curnode==2) TracePhotons(curx,cury,curz,0);
			if(curnode!=prevnode && prevnode==2) {atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);} //last point in a track
			prevnode=curnode;
			istep++;
			curx=track_start_pos_x+vx*istep*step;
			cury=track_start_pos_y+vy*istep*step;
			curz=track_start_pos_z+vz*istep*step;
		}while(curx < track_end_pos_x && cury < track_end_pos_y && curz < track_end_pos_z);
		if(prevnode==2) {
			printf("END: Adding end point, adding track to gEve..\n"); 
			atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);
			
		} //last point in a track

		/**
		//now down from the cg (z to negative)
		curx=cgx; cury=cgy; curz=cgz;
		geom->FindNode(curx,cury,curz); 
		node=geom->GetCurrentNode();
		curnode=node->GetNumber();
		prevnode=curnode;
		//printf("Initial node: %s\n",node->GetName());
		atrack=new TEveLine(2);
		atrack->SetPoint(0,curx,cury,curz);
		istep=0; 
		while(curx > track_start_pos_x && cury > track_start_pos_y && curz > track_start_pos_z)
		{
			istep++;
			curx=cgx-vx*istep*step;
			cury=cgy-vy*istep*step;
			curz=cgz-vz*istep*step;
			geom->FindNode(curx,cury,curz); 
			node=geom->GetCurrentNode();
			curnode=node->GetNumber();
			// printf("Step %d curnode: %s, curnode number%d, curxyz: %f %f %f\n",istep,node->GetName(), curnode,curx,cury,curz);
			//   if(curnode==prevnode && curnode==1) {printf("Adding point..\n"); atrack->SetPoint(atrack->GetN(),curx,cury,curz);}
			if(curnode!=prevnode && curnode==2) { atrack=new TEveLine(2); atrack->SetPoint(0,curx,cury,curz);}
			if(istep==1) {TracePhotons(curx,cury,curz,2); printf("?");}
			if(istep!=1) TracePhotons(curx,cury,curz,2);
			//if(curnode==prevnode && curnode==2&&istep==1) TracePhotons(curx,cury,curz,1); //debugging case
			//if(curnode==prevnode && curnode==2&&istep!=1) TracePhotons(curx,cury,curz,0);
			if(curnode!=prevnode && prevnode==2) { atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);} //last point in a track
			prevnode=curnode;
		}
		if(prevnode==2) {
			printf("END: Adding end point, adding track to gEve..\n"); 
			atrack->SetPoint(1,curx,cury,curz), gEve->AddElement(atrack);
		} //last point in a track

		**/

   }
   

   bsipm->Fill();
   bcgx->Fill();
   bcgy->Fill();
   bcgz->Fill();
   blrf->Fill();

   //sprintf(str,"Event %d",eventID);
   //ann->SetText(str);
   gEve->Redraw3D(kFALSE);

}

TracePhotons(Float_t px, Float_t py, Float_t pz, Int_t rays=2)
{
	Float_t ldx,ldy,ldz; // coordinates of grid nodes at the light detectors
	Float_t ldd; //distance from emitting point to grid node
	Float_t  domega; //solid angle
	Double_t v0[3];
	Double_t v1[3];
	Double_t rv[3]; //ray vector
	Double_t nv[3]; //normal vector
	Double_t cosine=1; //cosine of the ray-to-normal angle
	Double_t cosiney=1; //cosine of the ray-to-normal angle
	Double_t cosinez=1; //cosine of the ray-to-normal angle
	TEveLine * rtrace;
	//TGeoRotation *rotls = new TGeoRotation("rotls",0,0,0);
	//printf("Photons emitted at %f %f %f\n",px,py,pz);

	for(int j=0;j<nsipm_x;j++){
	for(int l=0;l<nsipm_y;l++){
			t0=clock();
		//for(Float_t iy=0; iy<GRIDNODES; iy++) for(Float_t ix=0; iz<GRIDNODES; ix++){
			// ACL1
			//calculate coordinates and distance
			//coordinates in a light detector plane
			//t0 = clock(); 
			ldx= sipm_shift_x-TPC_size_x/2.+(j+0.5)*sipm_space_x;
			ldy= sipm_shift_y-TPC_size_y/2.+(l+0.5)*sipm_space_y;
			ldz= sipm_shift_z;
			

			//if(rays==1&&j==4){
			//	rtrace=new TEveLine(2);
			//	rtrace->SetPoint(0,px,py,pz);
			//	rtrace->SetPoint(1,ldx,ldy,ldz);
			//	gEve->AddElement(rtrace);
			//}
			rv[0]=px-ldx; rv[1]=py-ldy; rv[2]=pz-ldz;
			ldd=sqrt( rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2] );
			//check if the emission point is actially seen from the detector
			//geom->SetCurrentDirection(-rv[0],-rv[1],-rv[2]);
			//geom->FindNode(px,py,pz);
			//node=geom->FindNextBoundary();
			//nv=geom->FindNormal(kTRUE);
			nv[0]=0;
			nv[1]=0;
			nv[2]=1;
			//t1=clock();
			//Calculate solid angle
			//printf("ldx: %f ldy: %f ldz: %f",rv[0],rv[1],rv[2]);
			cosine=(nv[0]*rv[0]+nv[1]*rv[1]+nv[2]*rv[2])/(ldd*1.0); //normal vector has unit length
			//t2=clock();
			//cosiney=1./sqrt(rv[1]*rv[1]/(rv[0]*rv[0])+1.);
			//cosinez=1./sqrt(rv[2]*rv[2]/(rv[0]*rv[0])+1.);
			//domega = cosine*(sipm_size_x)*(sipm_size_y)/(ldd*ldd);
			domega=cosine*4*TMath::ASin((sipm_size_x)*(sipm_size_y)/sqrt((4.0*ldd*ldd+(sipm_size_x)*(sipm_size_x))*(4.0*ldd*ldd+(sipm_size_y)*(sipm_size_y))));
			//domega=4*TMath::ASin(cosiney*cosinez*(acl1_size_y/GRIDNODES)*(acl1_size_z/GRIDNODES)/sqrt((4.0*ldd*ldd+cosiney*cosiney*(acl1_size_y/GRIDNODES)*(acl1_size_y/GRIDNODES))*(4.0*ldd*ldd+cosinez*cosinez*(acl1_size_z/GRIDNODES)*(acl1_size_z/GRIDNODES))));
			//if(strcmp(node->GetName(), "CATHODE_18")!=0) sipm_vis[j*nsipm_x+l]=sipm_vis[j*nsipm_x+l]+1.*domega/(4*3.14152);
			sipm_vis[j*nsipm_x+l]+=1.*domega/(4*3.14152);
			//t3=clock();
			//printf("ini: %f, cos: %f, domega: %f\n", (double)(t1-t0), (double)(t2-t1), (double)(t3-t2));
			if(rays==2){
				rtrace=new TEveLine(2);
				trace->SetPoint(0,px,py,pz);
				rtrace->SetPoint(1,ldx,ldy,ldz);
				gEve->AddElement(rtrace);
				printf("domega %.10f cosine %f\n",domega,cosine);
			}
		t1=clock();
		sipm_npe[j*nsipm_x+l] = sipm_vis[j*nsipm_x+l] * dqdx * lyield;
		t2=clock();
		//printf("vis: %f, npe: %f\n", (double)(t1-t0), (double)(t2-t1));
		//printf("npe: %f\n", sipm_npe[j*nsipm_x+l]);
		//}
	}
	}
}

 