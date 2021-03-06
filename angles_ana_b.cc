#include "fwdet_res.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"

#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"

#include "hloop.h"
#include "hcategory.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphErrors.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

bool primary_vertex=false;//should I take into consideration onlu primary particles?
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & beamVector, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
//do przeniesienia do osobnego pliku**************************************************
    double trackDistance(HParticleCand* track1, HParticleCand*  track2)
    {
      double dist;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
      return dist;
    }

    double trackDistance(HParticleCand* track1, HFwDetCand*  track2)
    {
      double dist;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);

      base_2.setX(track2->getPointX());
      base_2.setY(track2->getPointY());
      base_2.setZ(track2->getPointZ());
      dir_2.setX(track2->getDirTx());
      dir_2.setY(track2->getDirTy());
      dir_2.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
  
      dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
      return dist;
    }

    HGeomVector trackVertex(HParticleCand* track1, HParticleCand*  track2)
    {
      HGeomVector ver;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
      return ver;
    }

    HGeomVector trackVertex(HParticleCand* track1, HFwDetCand*  track2)
    {
      HGeomVector ver;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);

      base_2.setX(track2->getPointX());
      base_2.setY(track2->getPointY());
      base_2.setZ(track2->getPointZ());
      dir_2.setX(track2->getDirTx());
      dir_2.setY(track2->getDirTy());
      dir_2.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
  
      ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
      return ver;
    }

    Int_t getMotherIndex(HGeantKine* particle)
    {
      Int_t trackID=particle->getTrack();
      HGeantKine* particleParent=particle->getParent(trackID);
      Int_t parentID=0;
      if(particleParent!=0)
	parentID=particleParent->getID();
      return parentID;
      
      //return particle->getGeneratorInfo1();
    }


    bool isLepton(HParticleCand* particle)
    {
      double delta=0.05;
      double mquality=particle->getRichMatchingQuality();

      double dphi=particle->getDeltaPhi();
      double dtheta=particle->getDeltaTheta();
      
          
      if(particle->isFlagBit(kIsUsed))
	{
	  if(mquality==-1 || mquality>5)
	    return false;
	  if(particle->getBeta()<(1-delta) || particle->getBeta()>(1+delta))
	    return false;
	  // if(dtheta<-0.4 || dtheta>0.4)
	  // return false;
	  // if(dphi*TMath::Sin(particle->getTheta())>0.4 || dphi*TMath::Sin(particle->getTheta())<-0.4)
	  // return false;
	}
      else
	return false;
     
      return true;
    }
    //***********************************************************************************

 
 	


Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
    if (!loop->setInput(""))
    {                                                    // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    //////////////////////////////////////////////////////////////////////////////
    //      Fast tree builder for creating of ntuples                            //
    //////////////////////////////////////////////////////////////////////////////

    loop->printCategories();    // print all categories found in input + status

    HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fFwDetStrawCal = nullptr;
    fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");
    if (!fFwDetStrawCal)
    {
        cout << "No catFwDetStrawCal!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fwDetCatSim = nullptr;
    fwDetCatSim = HCategoryManager::getCategory(catFwDetCand, kTRUE, "catFwDetCand");
    if (!fwDetCatSim)
    {
        cout << "No catFwDetCand!" << endl;
	//exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * particleCatSim= nullptr;
    particleCatSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if(!particleCatSim)
      {
	cout<< "No catParticleCandSim!"<<endl;
      }

    Int_t entries = loop->getEntries();
    //     //setting numbers of events regarding the input number of events by the user
    if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;

    //     // specify output file
    TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();
    //
    cout << "NEW ROOT TREE , vertex_ana" << endl;
    //
    
    //---------------------------HISTOGRAMS-----------------------------------------------------

    TCanvas* cLambdaFW=new TCanvas("cLambdaFW","Pion-FwDet specta for FwDet");
    TH1F* hLdistance=new TH1F("hLdistance","Distance between pion in HADES and any particle in FwDet",1000,0,300);
    TH1F* hLmass=new TH1F("hLmass","Mass for pairs pion-FwDet",1000,800,2200);
    TH1F* hLmassDist=new TH1F("hLmassDist","Mass for pairs pion-FwDet after distance cut",1000,800,2200);
    TH1F* hLmassDistMass=new TH1F("hLmassDistMass","Mass for pairs pion-FwDet after distance and mass cut",40,1100,1140);
    TH1F* hLFwDetCh2=new TH1F("hLFwDetCh2","Ch2 for FwDet tracks",1800,0,900);
    TH2F* h2Lvertex=new TH2F("h2Lvertex","Lambda 1115 vertex;z[mm];r[mm]",300,-100,200,200,0,200);

    TCanvas* cLambdaH=new TCanvas("cLambdaH","Pion-FwDet specta for HADES");
    TH1F* hLHdistance=new TH1F("hLdistance","Distance between pion in HADES and any particle in FwDet",1000,0,300);
    TH1F* hLHmass=new TH1F("hLmass","Mass for pairs pion-FwDet",1000,800,2200);
    TH1F* hLHmassDist=new TH1F("hLmassDist","Mass for pairs pion-FwDet after distance cut",1000,800,2200);
    TH1F* hLHmassDistMass=new TH1F("hLmassDistMass","Mass for pairs pion-FwDet after distance and mass cut",40,1100,1140);
    //H1F* hLHFwDetCh2=new TH1F("hLFwDetCh2","Ch2 for FwDet tracks",1800,0,900);
    TH2F* h2LHvertex=new TH2F("h2Lvertex","Lambda 1115 vertex;z[mm];r[mm]",300,-100,200,200,0,200);
    
    TCanvas* cLeptonId=new TCanvas("cLeptonId","Lepton identyfication");
    TH2F* h2LIdPhidTheta=new TH2F("h2LIdPhidTheta","d Phi vs sin(Th) dTh; dTH; sin(Th)*dPhi",4000,-40,40,2000,-40,40);
    TH1F* h1LImatchqu=new TH1F("h1LImatchqu","Match quality between mdc and rich",320,-2,30);
    TH2F* h2LIbethamom=new TH2F("h2LIbethamom","Betha vs p",1000,-1000,1000,1500,-0.2,1.3);
    TH1I* hLIparticleid=new TH1I("hLIparticleid","PID for particles in geant kine",30,0,30);
    TH1F* hLIleptonPID=new TH1F("hLIleptonPID","Nuber of leptons identyfied by PID",700,0,1400);
    TH1F* hLIleptonIdef=new TH1F("hLIleptonIdef","Nuber of leptons identyfied by lepton identyfication",700,0,1400);
    TH1F* hLImatchqLepton=new TH1F("hLImatchqLepton","Match quality between mdc and rich for leptons",320,-2,30);
    
    TCanvas* cIventInfo=new TCanvas("cIventInfo","info about evet");
    TH1I* hIIparticlecand=new TH1I("hIIparticlecand","number of particles in particle cand",20,0,20);
    TH1I* hIIparticleFwDet=new TH1I("hIIparticlecand","number of particles in FwDetCand",20,0,20);
    TH2F* h2IIleptonsfromPV=new TH2F("h2IIleptonsfromPV","momentum vs theta nangle for leptons from pv ;mom;theta",1000,0,3600,180,0,180);
    TH2F* h2IIprotons=new TH2F("h2IIprotonsfromPV","momentum vs theta angle for protons different then pv ;mom;theta",1000,0,5000,180,0,180);
    TH2F* h2IIpions=new TH2F("h2IIpions","momentum vs theta angle for pions from Lambda 1520 and bg ;mom;theta",1000,0,5000,180,0,180);
    TH2F* h2IIleptonsInAcceptance=new TH2F("h2IIleptonsInAcceptance","momentum vs theta nangle for leptons from pv in Hades acceptance ;mom;theta",1000,0,3600,180,0,180);

    TCanvas* cParticlesRelation=new TCanvas("cParticlesRelation","Relations between particles when some of them were detected");
    TH2F* h2PRleptons=new TH2F("h2PRleptons","momentum vs theta nangle for leptons when Lambda was detected ;mom;theta",200,0,3600,60,0,180);
    TH2F* h2PRprotons=new TH2F("h2PRprotons","momentum vs theta angle for protons when diLepton was detected ;mom;theta",200,0,3600,60,0,180);
    TH2F* h2PRpions=new TH2F("h2PRpions","momentum vs theta angle for pions when diLepton was detected ;mom;theta",200,0,3600,60,0,180);
    TH1I* hPRmotherindex=new TH1I("hPRmotherindex","mother index for all particles",30,0,30);
    TH1I* hPRLambdaDoughter=new TH1I("hPRLambdaDoughter","particles from Lambda",30,0,30);

    TCanvas* cDiLepton=new TCanvas("cDiLepton","Di-lepton spectrum");
    TH1F* hDLdistance=new TH1F("hDLdistance","Distance between leptons",1000,0,300);
    TH1F* hDLmass=new TH1F("hDLmass","Di lepton spectrum for all pairs",400,0,800);
    TH1F* hDLmassDist=new TH1F("hDLmassDist","Di lepton spectrum for tracks closer then x",400,0,800);
    TH1F* hDLmassDistOA=new TH1F("hDLmassDistOA","Di lepton spectrum for tracks closer then x and opening angle more then y",400,0,800);
    TH1F* hDLopeningangle=new TH1F("hDLopeningangle","Di lepton opening angle spectrum",720,0,180);

    TCanvas* cLambda1520=new TCanvas("cLambda1520","Pion-FwDet specta");
    TH1F* hL1520distance=new TH1F("hL1520distance","Distance between di-lepton vector and lambda vector",150,0,300);
    TH1F* hL1520mass=new TH1F("hL1520mass","Mass for all pairs di-Lepton, lambda",600,1,-1);
    TH1F* hL1520massDist=new TH1F("hL1520massDist","Mass for pairs di-Lepton Lambda after distance cut",600,1,-1);
    TH2F* h2L1520vertex=new TH2F("h2L1520vertex","Reconstructed vertex for Lambda 1520 decay",50,-200,50,10,0,50);

    TCanvas* cEff=new TCanvas("cEff","Detection efficiency for different particles");
    TH1F* hEprotons4Pi=new TH1F("hEprotons4Pi","Protons from kine; theta",180,0,180);
    TH1F* hEprotonsdet=new TH1F("hEprotonsdet","Protons registered in detector; theta",180,0,180);
    TH1F* hEprotonsEff=new TH1F("hEprotonsEff","Detection efficiency; theta;efficiency",180,0,180);
    TH1F* hEpions4Pi=new TH1F("hEpions4Pi","Pions from kine; theta",180,0,180);
    TH1F* hEpionsdet=new TH1F("hEpionsdet","Pions registered in detector; theta",180,0,180);
    TH1F* hEpionsEff=new TH1F("hEpionsEff","Detection efficiency; theta;efficiency",180,0,180);
    TH1F* hEleptons4Pi=new TH1F("hEleptons4Pi","Leptons from kine; theta",180,0,180);
    TH1F* hEleptonsdet=new TH1F("hEleptonsdet","Leptons registered in detector; theta",180,0,180);
    TH1F* hEleptonsEff=new TH1F("hEleptonsEff","Detection efficiency; theta;efficiency",180,0,180);
    TH2F* h2Eproton4Pi=new TH2F("h2Eproton4Pi","Proton emitted in 4Pi; phi; theta",100,0,360,50,0,180);
    TH2F* h2Eprotondet=new TH2F("h2Eprotondet","Proton detected in 4Pi; phi; theta",100,0,360,50,0,180);
    TH2F* h2EprotonEff=new TH2F("h2EprotonEff","Proton efficiency in 4Pi; phi; theta",100,0,360,50,0,180);
    TH2F* h2Epion4Pi=new TH2F("h2Epion4Pi","Pion emitted in 4Pi; phi; theta",100,0,360,50,0,180);
    TH2F* h2Epiondet=new TH2F("h2Epiondet","Pion detected in 4Pi; phi; theta",100,0,360,50,0,180);
    TH2F* h2EpionEff=new TH2F("h2EpionEff","Pion efficiency in 4Pi; phi; theta",100,0,360,50,0,180);

    TCanvas* cEffFromLambda=new TCanvas("cEffFromLambda","Detection efficiency for particles from Lambda(1115 decay)");
    TH1F* hEFLprotons4Pi=new TH1F("hEFLprotons4Pi","Protons from kine; theta",180,0,180);
    TH1F* hEFLprotonsdet=new TH1F("hEFLprotonsdet","Protons registered in detector; theta",180,0,180);
    TH1F* hEFLprotonsEff=new TH1F("hEFLprotonsEff","Detection efficiency; theta;efficiency",180,0,180);
    
    TCanvas* cFwDet=new TCanvas("cFWDet","Foward tracker calibration");
    TH1F* hFDdistanceLambda=new TH1F("hFDdistanceLambda","Distance between pion in HADES and proton from Lambda(1115)",1000,0,300);
    TH1F* hFDdistanceBg=new TH1F("hFDdistanceBg","Distance between pion in HADES and proton from background",1000,0,300);
    TH2F* h2FDdetectedProtons=new TH2F("h2FDdetectedProtons","all protons detected inFwDet; mom; theta",200,0,4000,40,0,40);
    TH2F* h2FDsimProtons=new TH2F("h2FDsimProtons","all protons simulated inFwDet; mom; theta",200,0,4000,40,0,40);

    int zmin=-50;
    int zmax=200;
    int zbin=25;
    TCanvas* cEffInZ=new TCanvas("cEffInZ","Detection Efficiency in function of Z");
    TH1F* hEIZpionEff=new TH1F("hEIZpionEff","pion efficiency in function of Z vertex coordinate",zbin,zmin,zmax);
    TH1F* hEIZpionDet=new TH1F("hEIZpionDet","pions detected in function of Z vertex coordinate",zbin,zmin,zmax);
    TH1F* hEIZpionSim=new TH1F("hEIZpionSim","pions simulated in function of Z vertex coordinate",zbin,zmin,zmax);
    TH1F* hEIZprotonEff=new TH1F("hEIZprotonEff","proton efficiency in function of Z vertex coordinate",zbin,zmin,zmax);
    TH1F* hEIZprotonDet=new TH1F("hEIZprotonDet","protons detected in function of Z vertex coordinate",zbin,zmin,zmax);
    TH1F* hEIZprotonSim=new TH1F("hEIZprotonSim","protons simulated in function of Z vertex coordinate",zbin,zmin,zmax);

    TCanvas* cLambdaReal=new TCanvas("cLambdaReal","Values for pairs from Lambda(1115) decay");
    TH1F* hLRmass=new TH1F("hLRmass","Mass for pion-proton from real Lambda(1115) decay",1000,1,-1);
TH1F* hLRdist=new TH1F("hLRdist","Distance for pion-proton from real Lambda(1115) decay",100,1,-1);

 TCanvas* cDetectedParticles=new TCanvas("cDetectedParticles","particles detected");
 TH1F* hDPfromL1115=new TH1F("hDPfromL1115","Particles detected form L(1115) dacay",30,0,30);

 TCanvas* cDiLeptonExcl=new TCanvas("cDiLeptonExcl","Di lepton Lambda, exclusive");
 TH1F* hDLEoa=new TH1F("hDLEoa","Di lepton spedtrum",150,0,450);
 TH1F* hDLEwoa=new TH1F("hDLEwoa","Di lepton spedtrum, for opening angle > 4 deg",150,0,450);
 //TH1F* hVRIZ [7];
    //TF1* fVRIZgaus[7];
    //TGraphErrors* gVRIZ=new TGraphErrors(7);
    /*for(int i=0; i<7;i++)
      {
	ostringstream hname;
	hname << "resolution_in_h_" << i;
	ostringstream gname;
	gname << "gaus_no_" << i;
	hVRIZ[i]=new TH1F(hname.str().c_str(),hname.str().c_str(),100,-300,300);
	fVRIZgaus[i]=new TF1(gname.str().c_str(),"gaus",-300,300);
	}*/
    //---------------------------------------------------------------------------------------
    
    //event loop *************************************************
    //*********************************************************
    for (Int_t i = 0; i < entries; i++)                   
    {
        loop->nextEvent(i);         // get next event. categories will be cleared before
        if(i%5000==0)
	  cout<<"event no. "<<i<<endl;
        	// Geant pairs ana
	int hnum=particleCatSim->getEntries();
	int fnum=fwDetCatSim->getEntries();
	int knum=fCatGeantKine->getEntries();
	//HParticleCandSim* lep0=nullptr;
	HParticleCandSim* lep1=nullptr;
	HParticleCandSim* lep2=nullptr;
	HParticleCandSim* pion=nullptr;
	HParticleCandSim* lep=nullptr;
	HParticleCandSim* protonH=nullptr;
	HFwDetCandSim*  proton=nullptr;
	HGeantKine* kine=nullptr;

	HParticleTool tool;
	
	HGeomVector vertexL;
	HGeomVector vertexDL;
	HGeomVector vertexL1520;
	HGeomVector dirL;
	HGeomVector dirDL;

	double min_dist_dl=10;
	double min_dist_l=20;
	double min_angle=4;

	hIIparticlecand->Fill(hnum);
	hIIparticleFwDet->Fill(fnum);
	
	//lepton identification*********************************
	for(int n=0;n<hnum;n++)
	  {
	    lep=HCategoryManager::getObject(lep, particleCatSim,n);
	    if(lep->getRichMatchingQuality()!=-1 && lep->isFlagBit(kIsUsed))
	      {
		double dphi=lep->getDeltaPhi();
		double dtheta=lep->getDeltaTheta();
		double mom=lep->getMomentum();
		h1LImatchqu->Fill(lep->getRichMatchingQuality());
		h2LIdPhidTheta->Fill(dtheta,dphi);
		h2LIbethamom->Fill(lep->getMomentum()*lep->getCharge(),lep->getBeta());

		if(lep->getGeantPID()==2 || lep->getGeantPID()==3)
		  {
		    hLIleptonPID->Fill(mom);
		  }
		if(!(lep->getGeantPID()==2 || lep->getGeantPID()==3))	    
		  {
		    hLImatchqLepton->Fill(lep->getRichMatchingQuality());
		  }
		if(isLepton(lep))
		  {
		    hLIleptonIdef->Fill(mom);
		  }
	      }
	  }
	//end of lepton ident***********************************

	//find any lepton
	for(int j=0;j<hnum;j++) //first particle(electron)
	  {
	    lep1=HCategoryManager::getObject(lep1, particleCatSim,j);
	    Int_t lep1ID=lep1->getGeantPID();

	    if((lep1ID==2 || lep1ID==3) && lep1->getRichMatchingQuality()!=-1 && lep1->getGeantParentTrackNum()==0 && lep1->isFlagBit(kIsUsed))
	      {
		lep1->calc4vectorProperties(HPhysicsConstants::mass(lep1->getGeantPID()));
		h2IIleptonsInAcceptance->Fill(lep1->getMomentum(),lep1->getTheta());
		hEleptonsdet->Fill(lep1->getTheta());
	      }
	    else
	      continue; //kick off bad tracks

	    //DiLepton part*************************************
	    if(lep1ID==2)//electron
	      for(int k=0;k<hnum;k++)//second particle (positon)
		{
		  lep2=HCategoryManager::getObject(lep2, particleCatSim,k);
		  Int_t lep2ID=lep2->getGeantPID();
		  
		  if(lep2ID==3 && lep2->getRichMatchingQuality()!=-1 && lep2->getGeantParentTrackNum()==0 && lep2->isFlagBit(kIsUsed))
		    {
		      //ideal energy******************************
		      /*lep1->SetPx(lep1->getGeantxMom());
		      lep1->SetPy(lep1->getGeantyMom());
		      lep1->SetPz(lep1->getGeantzMom());

		      lep2->SetPx(lep2->getGeantxMom());
		      lep2->SetPy(lep2->getGeantyMom());
		      lep2->SetPz(lep2->getGeantzMom());

				  
		      lep1->SetE(TMath::Sqrt(lep1->getGeantTotalMom()*lep1->getGeantTotalMom()+HPhysicsConstants::mass(lep1->getGeantPID())*HPhysicsConstants::mass(lep1->getGeantPID())));
		      lep2->SetE(TMath::Sqrt(lep2->getGeantTotalMom()*lep2->getGeantTotalMom()+HPhysicsConstants::mass(lep2->getGeantPID())*HPhysicsConstants::mass(lep2->getGeantPID())));
		      */ //end ideal energy**************************
				  
				 
		      double diLeptonM=(*lep1+*lep2).M();
		      vertexDL=trackVertex(lep1,lep2);
		      dirDL.setXYZ((*lep1+*lep2).X(),(*lep1+*lep2).Y(),(*lep1+*lep2).Z());
			  
		      hDLmass->Fill(diLeptonM);
		      hDLdistance->Fill(trackDistance(lep1,lep2));
		      hDLopeningangle->Fill(tool.getOpeningAngle(lep1,lep2));
		      		      	 
		      hDLmassDist->Fill(diLeptonM);
		      if(tool.getOpeningAngle(lep1,lep2)>min_angle)
			{
			  hDLmassDistOA->Fill(diLeptonM);
			}
		    }
		}//end of second particle	      
	  }//end of first particle

	//Lambda(1115)
	for(int l=0;l<hnum;l++)//third particle (pion)
	  {
	    pion=HCategoryManager::getObject(pion, particleCatSim,l);
	    if(pion->getGeantPID()==9 && pion->isFlagBit(kIsUsed) && pion->getGeantParentPID()==18)//pion from lambda decay
	      {
		hEpionsdet->Fill(pion->getTheta());
		h2Epiondet->Fill(pion->getPhi(),pion->getTheta());
		hEIZpionDet->Fill(pion->getGeantzVertex());
	      }
	    if(pion->getGeantPID()==14 && pion->isFlagBit(kIsUsed) && (pion->getGeantParentPID()==18 || pion->getGeantCreationMechanism()==0))//proton from lambda decay and primary production
	      {
		hEprotonsdet->Fill(pion->getTheta());
		h2Eprotondet->Fill(pion->getPhi(),pion->getTheta());
	      }
	    if(pion->getGeantPID()==14 && pion->isFlagBit(kIsUsed) && pion->getGeantParentPID()==18)//proton from Lambda(1115) decay
	      {
		hEFLprotonsdet->Fill(pion->getTheta());
	      }
	    if(pion->isFlagBit(kIsUsed) && pion->getGeantParentPID()==18)
	      hDPfromL1115->Fill(pion->getGeantPID());
	    for(int m=0;m<fnum;m++)//fourth particle (any in fwDet)
	      {		
		proton=HCategoryManager::getObject(proton,fwDetCatSim,m);

		hLFwDetCh2->Fill(proton->getChi2());
		if(pion->getGeantPID()==9 && pion->isFlagBit(kIsUsed) /*&& proton->getChi2()<10*/)
		  {
		    pion->calc4vectorProperties(HPhysicsConstants::mass(pion->getGeantPID()));
		    proton->calc4vectorProperties(HPhysicsConstants::mass(14));

		    //ideal energy***********************************************
		    /*pion->SetPx(pion->getGeantxMom());
		    pion->SetPy(pion->getGeantyMom());
		    pion->SetPz(pion->getGeantzMom());

		    pion->SetE(TMath::Sqrt(pion->getGeantTotalMom()*pion->getGeantTotalMom()+HPhysicsConstants::mass(pion->getGeantPID())*HPhysicsConstants::mass(pion->getGeantPID())));
		    *///ideal energy**********************************************
		 		      
		    double lambdaM=(*pion+*proton).M();
		    double lambdaD=trackDistance(pion,proton);

		    vertexL=trackVertex(pion,proton);
		    dirL.setXYZ((*pion+*proton).X(),(*pion+*proton).Y(),(*pion+*proton).Z());
			      
		    hLmass->Fill(lambdaM);
		    hLdistance->Fill(lambdaD);
		    //recognize proton orygin----------------
		    int tn=proton->getTrack();
		    int kineID=-1;
		    int kineparentID=-1;
		    for(int p=0;p<knum;p++)
		      {
			kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
			int kineT=kine->getTrack();
			if(kineT==tn)
			  {
			    kineID=kine->getID();
			    kineparentID=getMotherIndex(kine);
			    break;
			  }
		      }
		    if(kineID==14 && kineparentID==18)
		      hFDdistanceLambda->Fill(lambdaD);

		    if(kineID==14 && kineparentID!=18)
		      hFDdistanceBg->Fill(lambdaD);
		    //---------------------------------------
		    if(kineparentID==18 && pion->getGeantParentPID()==18)//proton and pion from L(1115)
		      {
			hLRmass->Fill(lambdaM);
			hLRdist->Fill(lambdaD);
		      }
		    
		    if(lambdaD<min_dist_l /*&& proton->Theta()>1*TMath::Pi()/180*/)
		      {
			hLmassDist->Fill(lambdaM);
			h2Lvertex->Fill(vertexL.Z(),TMath::Sqrt(vertexL.X()*vertexL.X()+vertexL.Y()*vertexL.Y()));
		
			if(lambdaM<1120 && lambdaM>1108)
			  {
			    // isLambda=1;
			    hLmassDistMass->Fill(lambdaM);

			    //Lambda(1520)*********************
			    for(int j=0;j<hnum;j++)
			      for(int k=0;k<hnum;k++)
				{
				  lep1=HCategoryManager::getObject(lep1, particleCatSim,j);
				  Int_t lep1ID=lep1->getGeantPID();
				  lep2=HCategoryManager::getObject(lep2, particleCatSim,k);
				  Int_t lep2ID=lep2->getGeantPID();

				  if(lep2ID==3 && lep2->getRichMatchingQuality()!=-1 && lep2->getGeantParentTrackNum()==0 && lep2->isFlagBit(kIsUsed))//Take only accepted positon from PV 
				    if(lep1ID==2 && lep1->getRichMatchingQuality()!=-1 && lep1->getGeantParentTrackNum()==0 && lep1->isFlagBit(kIsUsed))//Take only accepted electron from PV
				      {
					TLorentzVector lvLambda=*pion+*proton;
					TLorentzVector lvDiLepton=*lep1+*lep2;
					double mass_diLepton=lvDiLepton.M();
					double mass_1520=(lvLambda+lvDiLepton).M();
					double distance_1520=tool.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
					hDLEoa->Fill(mass_diLepton);
					if(tool.getOpeningAngle(lep1,lep2)>min_angle)//opening angle condition
					  {
					    vertexL1520=tool.calcVertexAnalytical(vertexL,dirL,vertexDL,dirDL);
					    hDLEwoa->Fill(mass_diLepton);
					    hL1520mass->Fill(mass_1520);
					    hL1520distance->Fill(distance_1520);
					    if(distance_1520<min_dist_l)
					      {
						hL1520massDist->Fill(mass_1520);
						h2L1520vertex->Fill(vertexL1520.Z(),TMath::Sqrt(vertexL1520.X()*vertexL1520.X()+vertexL1520.Y()*vertexL1520.Y()));
					      }
					  }
				      }
				}
			    //Lambda(1520)-the end************
			  }
		      }
		  }

	      }//end of proton in FwDet
	    for(int z=0; z<hnum; z++)//proton in HADES
	      {
		pion=HCategoryManager::getObject(pion, particleCatSim,l);
		protonH=HCategoryManager::getObject(protonH,particleCatSim,z);

		if(pion->getGeantPID()==9 && pion->isFlagBit(kIsUsed) && protonH->getGeantPID()==14 && protonH->isFlagBit(kIsUsed))
		  {
		    pion->calc4vectorProperties(HPhysicsConstants::mass(pion->getGeantPID()));
		    protonH->calc4vectorProperties(HPhysicsConstants::mass(14));

		    //ideal energy***********************************************
		    /*pion->SetPx(pion->getGeantxMom());
		    pion->SetPy(pion->getGeantyMom());
		    pion->SetPz(pion->getGeantzMom());

		    pion->SetE(TMath::Sqrt(pion->getGeantTotalMom()*pion->getGeantTotalMom()+HPhysicsConstants::mass(pion->getGeantPID())*HPhysicsConstants::mass(pion->getGeantPID())));

		    protonH->SetPx(protonH->getGeantxMom());
		    protonH->SetPy(protonH->getGeantyMom());
		    protonH->SetPz(protonH->getGeantzMom());

		    protonH->SetE(TMath::Sqrt(protonH->getGeantTotalMom()*protonH->getGeantTotalMom()+HPhysicsConstants::mass(protonH->getGeantPID())*HPhysicsConstants::mass(protonH->getGeantPID())));
		    */ //ideal energy**********************************************
		 		      
		    double lambdaM=(*pion+*protonH).M();
		    double lambdaD=trackDistance(pion,protonH);

		    vertexL=trackVertex(pion,protonH);
		    dirL.setXYZ((*pion+*protonH).X(),(*pion+*protonH).Y(),(*pion+*protonH).Z());
			      
		    hLHmass->Fill(lambdaM);
		    hLHdistance->Fill(lambdaD);

		    if(pion->getGeantParentPID()==18 && protonH->getGeantParentPID()==18)
		      {
			hLRmass->Fill(lambdaM);
			hLRdist->Fill(lambdaD);
		      }
			
		    if(lambdaD<min_dist_l)
		      {
			hLHmassDist->Fill(lambdaM);
			h2LHvertex->Fill(vertexL.Z(),TMath::Sqrt(vertexL.X()*vertexL.X()+vertexL.Y()*vertexL.Y()));
			if(lambdaM<1120 && lambdaM>1108)//dilepton maching
			  {
			    hLHmassDistMass->Fill(lambdaM);

			    //Lambda(1520)*********************
			    for(int j=0;j<hnum;j++)
			      for(int k=0;k<hnum;k++)
				{
				  lep1=HCategoryManager::getObject(lep1, particleCatSim,j);
				  Int_t lep1ID=lep1->getGeantPID();
				  lep2=HCategoryManager::getObject(lep2, particleCatSim,k);
				  Int_t lep2ID=lep2->getGeantPID();
				  
				  if(lep2ID==3 && lep2->getRichMatchingQuality()!=-1 && lep2->getGeantParentTrackNum()==0 && lep2->isFlagBit(kIsUsed))//Take only accepted positon from PV 
				    if(lep1ID==2 && lep1->getRichMatchingQuality()!=-1 && lep1->getGeantParentTrackNum()==0 && lep1->isFlagBit(kIsUsed))//Take only accepted electron from PV
				      {
					TLorentzVector lvLambda=*pion+*protonH;
					    TLorentzVector lvDiLepton=*lep1+*lep2;
					    double mass_dilepton=lvDiLepton.M();
					    double mass_1520=(lvLambda+lvDiLepton).M();
					    double distance_1520=tool.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
					    vertexL1520=tool.calcVertexAnalytical(vertexL,dirL,vertexDL,dirDL);
					    hDLEoa->Fill(mass_dilepton);
					if(tool.getOpeningAngle(lep1,lep2)>min_angle)//opening angle condition
					  {
					    hDLEwoa->Fill(mass_dilepton);					  
					    hL1520mass->Fill(mass_1520);
					    hL1520distance->Fill(distance_1520);
					    if(distance_1520<min_dist_l)
					      {
						hL1520massDist->Fill(mass_1520);
						h2L1520vertex->Fill(vertexL1520.Z(),TMath::Sqrt(vertexL1520.X()*vertexL1520.X()+vertexL1520.Y()*vertexL1520.Y()));
					      }
					  }
				      }
				}
			    //Lambda(1520)-the end************
			  }//end of dilepton maching
		      }
		  }
	      }//end of proton in HADES
	    }//end of third particle (pion)	

	//Only FwDet-----------------------------------------
	for(int m=0;m<fnum;m++)
	  {		
	    proton=HCategoryManager::getObject(proton,fwDetCatSim,m);
   	    int tn=proton->getTrack();
	    int kineID=-1;
	    int kineparentID=-1;
	    
	    for(int p=0;p<knum;p++)
	      {
		kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
		int kineT=kine->getTrack();
		if(kineT==tn)
		  {
		    kineID=kine->getID();
		    kineparentID=getMotherIndex(kine);
		    break;
		  }
	      }
	    if(kineID==14)
	      {
		proton->calc4vectorProperties(HPhysicsConstants::mass(14));
		hEprotonsdet->Fill(proton->Theta()*180/TMath::Pi());
		h2FDdetectedProtons->Fill(proton->Rho(),proton->Theta()*180/TMath::Pi());
		h2Eprotondet->Fill(proton->Phi()*180/TMath::Pi(),proton->Theta()*180/TMath::Pi());
	      }
	    if(kineID==14 && kineparentID==18)
	      {
		hEFLprotonsdet->Fill(proton->Theta()*180/TMath::Pi());
		hDPfromL1115->Fill(14);
	      }
	  }
	//FwDet, the End--------------------------------------

	
	//kine analysis*****************************************
	for(int p=0;p<knum;p++)
	  {
	    kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
	    int kineID=kine->getID();
	    int mech=kine->getMechanism();
	    int kineparentID=getMotherIndex(kine);
	    HGeomVector lambdaVertex;

	    hPRmotherindex->Fill(kineparentID);

	    if(kineparentID==18)
	      hPRLambdaDoughter->Fill(kineID);

	    if(mech==0 && (kineID==2 || kineID==3))
	      {
		h2IIleptonsfromPV->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
		hEleptons4Pi->Fill(kine->getThetaDeg());
	      }
	    if(kineID==14 && kineparentID==18)//proton from lambda
	      {
		h2IIprotons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
		hEprotons4Pi->Fill(kine->getThetaDeg());
		h2Eproton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
		hEFLprotons4Pi->Fill(kine->getThetaDeg());
	      }
	    if(kineID==14 && mech==0)//proton from primary vertex
	      {
		//h2IIprotons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
		hEprotons4Pi->Fill(kine->getThetaDeg());
		h2Eproton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
	      }
	    if(kineID==9 && kineparentID==18)//Pi- from Lambda
	      {
		kine->getVertex(lambdaVertex);
		h2IIpions->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
		hEpions4Pi->Fill(kine->getThetaDeg());
		h2Epion4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
		hEIZpionSim->Fill(lambdaVertex.getZ());
	      }
	    if(kineID==14 && kine->getThetaDeg()<6.5)
	      h2FDsimProtons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	  }
	//******************************************************end kine

    } // end eventloop
    //**********************************************


    
    //draw alll
    cLambdaFW->Divide(2,3);
    cLambdaFW->cd(1);
    hLdistance->Draw();
    cLambdaFW->cd(2);
    hLmass->Draw();
    cLambdaFW->cd(3);
    hLmassDist->Draw();
    cLambdaFW->cd(4);
    hLmassDistMass->Draw();
    cLambdaFW->cd(5);
    hLFwDetCh2->Draw();
    cLambdaFW->cd(6);
    h2Lvertex->Draw("COLZ");

    cLambdaH->Divide(2,3);
    cLambdaH->cd(1);
    hLHdistance->Draw();
    cLambdaH->cd(2);
    hLHmass->Draw();
    cLambdaH->cd(3);
    hLHmassDist->Draw();
    cLambdaH->cd(4);
    hLHmassDistMass->Draw();
    cLambdaH->cd(5);
    // hLHFwDetCh2->Draw();
    cLambdaH->cd(6);
    h2LHvertex->Draw("COLZ");

    
    cLeptonId->Divide(2,3);
    cLeptonId->cd(1);
    h2LIbethamom->Draw("COLZ");
    cLeptonId->cd(2);
    h2LIdPhidTheta->Draw("COLZ");
    cLeptonId->cd(3);
    h1LImatchqu->Draw();
    cLeptonId->cd(5);
    hLImatchqLepton->Draw();
    cLeptonId->cd(4);
    hLIleptonPID->Draw();
    cLeptonId->cd(6);
    hLIleptonIdef->Draw();

    cIventInfo->Divide(2,3);
    cIventInfo->cd(1);
    hIIparticleFwDet->Draw();
    cIventInfo->cd(2);
    hIIparticlecand->Draw();
    cIventInfo->cd(3);
    h2IIleptonsfromPV->Draw("COLZ");
    cIventInfo->cd(4);
    h2IIpions->Draw("COLZ");
    cIventInfo->cd(5);
    h2IIprotons->Draw("COLZ");
    cIventInfo->cd(6);
    h2IIleptonsInAcceptance->Draw("COLZ");

    cParticlesRelation->Divide(2,3);
    cParticlesRelation->cd(1);
    h2PRprotons->Draw("COLZ");
    cParticlesRelation->cd(2);
    h2PRpions->Draw("COLZ");
    cParticlesRelation->cd(3);
    h2PRleptons->Draw("COLZ");
    cParticlesRelation->cd(4);
    hPRmotherindex->Draw();
    cParticlesRelation->cd(5);
    hPRLambdaDoughter->Draw();

    cDiLepton->Divide(2,3);
    cDiLepton->cd(1);
    hDLmass->Draw();
    cDiLepton->cd(2);
    hDLmassDist->Draw();
    cDiLepton->cd(3);
    hDLmassDistOA->Draw();
    cDiLepton->cd(5);
    hDLdistance->Draw();
    cDiLepton->cd(6);
    hDLopeningangle->Draw();

    cLambda1520->Divide(2,2);
    cLambda1520->cd(1);
    hL1520mass->Draw();
    cLambda1520->cd(2);
    hL1520distance->Draw();
    cLambda1520->cd(3);
    hL1520massDist->Draw();
    cLambda1520->cd(4);
    h2L1520vertex->Draw("COLZ");

    cEff->Divide(4,3);
    cEff->cd(1);
    hEprotons4Pi->Draw();
    hEprotonsdet->SetLineColor(kRed);
    hEprotonsdet->Draw("SAME");
    cEff->cd(2);
    hEprotonsdet->Draw();
    cEff->cd(3);
    hEprotonsEff->Divide(hEprotonsdet,hEprotons4Pi);
    hEprotonsEff->Draw();
    cEff->cd(4);
    h2EprotonEff->Divide(h2Eprotondet,h2Eproton4Pi);
    h2EprotonEff->Draw("COLZ");
    cEff->cd(5);
    hEpions4Pi->Draw();
    hEpionsdet->SetLineColor(kRed);
    hEpionsdet->Draw("SAME");
    cEff->cd(6);
    hEpionsdet->Draw();
    cEff->cd(7);
    hEpionsEff->Divide(hEpionsdet,hEpions4Pi);
    hEpionsEff->Draw();
    cEff->cd(8);
    h2EpionEff->Divide(h2Epiondet,h2Epion4Pi);
    h2EpionEff->Draw("COLZ");
    cEff->cd(9);
    hEleptons4Pi->Draw();
    hEleptonsdet->SetLineColor(kRed);
    hEleptonsdet->Draw("SAME");
    cEff->cd(10);
    hEleptonsdet->Draw();
    cEff->cd(11);
    hEleptonsEff->Divide(hEleptonsdet,hEleptons4Pi);
    hEleptonsEff->Draw();

    cFwDet->Divide(2,2);
    cFwDet->cd(1);
    hFDdistanceBg->Draw();
    cFwDet->cd(2);
    hFDdistanceLambda->Draw();
    cFwDet->cd(1);
    hFDdistanceLambda->SetLineColor(kRed);
    hFDdistanceLambda->Draw("SAME");
    cFwDet->cd(3);
    h2FDdetectedProtons->Draw("COLZ");
    cFwDet->cd(4);
    h2FDsimProtons->Draw("COLZ");
    
    cEffInZ->Divide(3,2);
    cEffInZ->cd(1);
    hEIZpionSim->Draw();
    cEffInZ->cd(2);
    hEIZpionDet->Draw();
    cEffInZ->cd(3);
    hEIZpionEff->Divide(hEIZpionDet,hEIZpionSim);
    hEIZpionEff->Draw();
    cEffInZ->cd(4);
    hEIZprotonSim->Draw();
    cEffInZ->cd(5);
    hEIZprotonDet->Draw();
    cEffInZ->cd(6);
    hEIZprotonEff->Draw();

    cEffFromLambda->Divide(3);
    cEffFromLambda->cd(1);
    hEFLprotons4Pi->Draw();
    cEffFromLambda->cd(2);
    hEFLprotonsdet->Draw();
    cEffFromLambda->cd(3);
    hEFLprotonsEff->Divide(hEFLprotonsdet,hEFLprotons4Pi);
    hEFLprotonsEff->Draw();

    cLambdaReal->Divide(2);
    cLambdaReal->cd(1);
    hLRmass->Draw();
    cLambdaReal->cd(2);
    hLRdist->Draw();

    cDetectedParticles->cd();
    hDPfromL1115->Draw();

    cDiLeptonExcl->Divide(2);
    cDiLeptonExcl->cd(1);
    hDLEoa->Draw();
    cDiLeptonExcl->cd(2);
    hDLEwoa->Draw();
    //save all
    hEprotons4Pi->Write();
    hEprotonsEff->Write();
    hEpionsEff->Write();
    hEleptonsEff->Write();

    
    cLeptonId->Write();
    cIventInfo->Write();
    cParticlesRelation->Write();
    cDiLepton->Write();
    cLambdaFW->Write();
    cLambdaH->Write();
    cLambda1520->Write();
    cEff->Write();
    cFwDet->Write();
    cEffInZ->Write();
    cEffFromLambda->Write();
    cLambdaReal->Write();
    cDetectedParticles->Write();
    cDiLeptonExcl->Write();


    hLHmassDist->Write();
    h2LHvertex->Write("COLZ");
    hLmassDist->Write();
    h2Lvertex->Write("COLZ");
    
    output_file->Close();
    
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();
    
    return 0;
}










