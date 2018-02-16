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
    TCanvas* cDiLepton=new TCanvas("cDiLepton","Di-lepton spectrum");
    TH1F* hDLdistance=new TH1F("hDLdistance","Distance between leptons",1000,0,300);
    TH1F* hDLmass=new TH1F("hDLmass","Di lepton spectrum for all pairs",400,0,800);
    TH1F* hDLmassDist=new TH1F("hDLmassDist","Di lepton spectrum for tracks closer then x",400,0,800);
    TH1F* hDLmassDistOA=new TH1F("hDLmassDistOA","Di lepton spectrum for tracks closer then x and opening angle more then y",400,0,800);
    TH1F* hDLopeningangle=new TH1F("hDLopeningangle","Di lepton opening angle spectrum",720,0,180);

    TCanvas* cLambda=new TCanvas("cLambda","Pion-FwDet specta");
    TH1F* hLdistance=new TH1F("hLdistance","Distance between pion in HADES and any particle in FwDet",1000,0,300);
    TH1F* hLmass=new TH1F("hLmass","Mass for pairs pion-FwDet",1000,800,2200);
    TH1F* hLmassDist=new TH1F("hLmassDist","Mass for pairs pion-FwDet after distance cut",1000,800,2200);
    TH1F* hLmassDistMass=new TH1F("hLmassDistMass","Mass for pairs pion-FwDet after distance and mass cut",40,1100,1140);
    TH1F* hLFwDetCh2=new TH1F("hLFwDetCh2","Ch2 for FwDet tracks",1800,0,900);

    TCanvas* cLambda1520=new TCanvas("cLambda1520","Pion-FwDet specta");
    TH1F* hL1520distance=new TH1F("hL1520distance","Distance between di-lepton vector and lambda vector",150,0,300);
    TH1F* hL1520mass=new TH1F("hL1520mass","Mass for all pairs di-Lepton, lambda",600,1,-1);
    TH1F* hL1520massDist=new TH1F("hL1520massDist","Mass for pairs di-Lepton Lambda after distance cut",600,1,-1);
    TH2F* h2L1520vertex=new TH2F("h2L1520vertex","Reconstructed vertex for Lambda 1520 decay",50,-200,50,10,0,50);

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
      
    int nLambdas=0;

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
	HParticleCandSim* lep0=nullptr;
	HParticleCandSim* lep1=nullptr;
	HParticleCandSim* lep2=nullptr;
	HParticleCandSim* pion=nullptr;
	HParticleCandSim* lep=nullptr;
	HParticleCandSim* protonH=nullptr;
	HFwDetCandSim*  proton=nullptr;
	HGeantKine* kine=nullptr;
	Int_t isDilepton=0;
	Int_t isLambda=0;
	Int_t isLambda1520=0;
	Int_t leptonparticle1=-1;
	Int_t leptonparticle2=-1;
	Int_t lambdaparticle1=-1;
	Int_t lambdaparticle2=-1;
	Int_t leptonparticle=-1;
	HParticleTool tool;
	
	HGeomVector vertexL;
	HGeomVector vertexDL;
	HGeomVector vertexL1520;
	HGeomVector dirL;
	HGeomVector dirDL;

	
	double min_dist_dl=10;
	double min_dist_l=50;
	double min_angle=4;

	hIIparticlecand->Fill(hnum);
	hIIparticleFwDet->Fill(fnum);
	
	//lepton identification*********************************
	for(int n=0;n<hnum;n++)
	  {
	    lep=HCategoryManager::getObject(lep, particleCatSim,n);
	    if(lep->getRichMatchingQuality()!=-1)
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

	//di-lepton and lambda reconstruction
	if(hnum>2)
	  {
	    leptonparticle=-1;
	    for(int j=0; j<hnum; j++)
	      {
		for(int k=j+1; k<hnum; k++)
		  {
		    for(int l=0; l<hnum; l++)
		      {
			for(int m=0; m<fnum; m++)
			  {			  
			    lep1=HCategoryManager::getObject(lep1, particleCatSim,j);
			    lep2=HCategoryManager::getObject(lep2, particleCatSim,k);
			    pion=HCategoryManager::getObject(pion, particleCatSim,l);
			    proton=HCategoryManager::getObject(proton,fwDetCatSim,m);
			    lep0=HCategoryManager::getObject(lep0, particleCatSim,j);

			    //lepton
			    if(leptonparticle!=j)//prevent double couting
			      if(lep0->getRichMatchingQuality()!=-1)
				if(lep0->getGeantPID()==2 || lep0->getGeantPID()==3)
				  {
				    leptonparticle=j;
				    lep0->calc4vectorProperties(HPhysicsConstants::mass(lep0->getGeantPID()));
				    h2IIleptonsInAcceptance->Fill(lep0->getMomentum(),lep0->getTheta());
				  }
			  
			    //DiLepton
			    if(j!=leptonparticle1 && k!=leptonparticle2)//prevent double counting
			      {
				if((lep1->getGeantPID()==2 && lep2->getGeantPID()==3)||(lep1->getGeantPID()==3 && lep2->getGeantPID()==2))
				  {
				    if((lep1->getGeantParentTrackNum()==0 && lep2->getGeantParentTrackNum()==0) && (lep1->isFlagBit(kIsUsed) && lep2->isFlagBit(kIsUsed)))//only primary vertex
				      {
					leptonparticle1=j;
					leptonparticle2=k;

					lep1->calc4vectorProperties(HPhysicsConstants::mass(lep1->getGeantPID()));
					lep2->calc4vectorProperties(HPhysicsConstants::mass(lep2->getGeantPID()));
				  
					//ideal energy******************************
					/*lep1->SetPx(lep1->getGeantxMom());
					lep1->SetPy(lep1->getGeantyMom());
					lep1->SetPz(lep1->getGeantzMom());

					lep2->SetPx(lep2->getGeantxMom());
					lep2->SetPy(lep2->getGeantyMom());
					lep2->SetPz(lep2->getGeantzMom());

				  
					lep1->SetE(TMath::Sqrt(lep1->getGeantTotalMom()*lep1->getGeantTotalMom()+HPhysicsConstants::mass(lep1->getGeantPID())*HPhysicsConstants::mass(lep1->getGeantPID())));
					lep2->SetE(TMath::Sqrt(lep2->getGeantTotalMom()*lep2->getGeantTotalMom()+HPhysicsConstants::mass(lep2->getGeantPID())*HPhysicsConstants::mass(lep2->getGeantPID())));*/
					//end ideal energy**************************
				  
				 
					double diLeptonM=(*lep1+*lep2).M();
					vertexDL=trackVertex(lep1,lep2);
					dirDL.setXYZ((*lep1+*lep2).X(),(*lep1+*lep2).Y(),(*lep1+*lep2).Z());
			  
					hDLmass->Fill(diLeptonM);
					hDLdistance->Fill(trackDistance(lep1,lep2));
					hDLopeningangle->Fill(tool.getOpeningAngle(lep1,lep2));
		      
					// if(trackDistance(lep1,lep2)<min_dist_dl)// only opening angle cut
			 
					hDLmassDist->Fill(diLeptonM);
					if(tool.getOpeningAngle(lep1,lep2)>min_angle)
					  {
					    hDLmassDistOA->Fill(diLeptonM);
					    isDilepton=1;
					    //cout<<"isDilepton"<<endl;
					  }
					else
					  isDilepton=0;
				      }
				    else
				      isDilepton=0;
				  }
				else
				  isDilepton=0;
			      }
		      
			    //Lambda part
			    if(lambdaparticle1!=l && lambdaparticle2!=m)//prevent double couting
			      {
				hLFwDetCh2->Fill(proton->getChi2());
				if(pion->getGeantPID()==9 &&/* pion->isFlagBit(kIsUsed) &&*/ proton->getChi2()<10)
				  {
				    lambdaparticle1=l;
				    lambdaparticle2=m;
				    pion->calc4vectorProperties(HPhysicsConstants::mass(pion->getGeantPID()));
				    proton->calc4vectorProperties(HPhysicsConstants::mass(14));

				    //ideal energy***********************************************
				    /*pion->SetPx(pion->getGeantxMom());
				    pion->SetPy(pion->getGeantyMom());
				    pion->SetPz(pion->getGeantzMom());

				    pion->SetE(TMath::Sqrt(pion->getGeantTotalMom()*pion->getGeantTotalMom()+HPhysicsConstants::mass(pion->getGeantPID())*HPhysicsConstants::mass(pion->getGeantPID())));
				    */
				    //ideal energy**********************************************
		       
			      
			      
				    double lambdaM=(*pion+*proton).M();
				    double lambdaD=trackDistance(pion,proton);

				    vertexL=trackVertex(pion,proton);
				    dirL.setXYZ((*pion+*proton).X(),(*pion+*proton).Y(),(*pion+*proton).Z());

			      
				    hLmass->Fill(lambdaM);
				    hLdistance->Fill(lambdaD);
			
				    if(lambdaD<min_dist_l)
				      {
					hLmassDist->Fill(lambdaM);
					if(lambdaM<1120 && lambdaM>1108)
					  {
					    isLambda=1;
					    hLmassDistMass->Fill(lambdaM);
					    //cout<<"isLambda"<<endl;
					  }
					else
					  isLambda=0;
				      }
				    else
				      isLambda=0;
				  }
				else
				  isLambda=0;
			      }

			    //Lambda 1520 part
		    		    
			    if(isLambda==1 && isDilepton==1)
			      {
				//cout<<"isLambda1520"<<endl<<endl;
				nLambdas++;
				TLorentzVector lvLambda=*pion+*proton;
				TLorentzVector lvDiLepton=*lep1+*lep2;

				double mass_1520=(lvLambda+lvDiLepton).M();
				double distance_1520=tool.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
				vertexL1520=tool.calcVertexAnalytical(vertexL,dirL,vertexDL,dirDL);

				hL1520mass->Fill(mass_1520);
				hL1520distance->Fill(distance_1520);
				if(distance_1520<min_dist_l)
				  {
				    hL1520massDist->Fill(mass_1520);
				    h2L1520vertex->Fill(vertexL1520.Z(),TMath::Sqrt(vertexL1520.X()*vertexL1520.X()+vertexL1520.Y()*vertexL1520.Y()));
				    isLambda1520=1;
				  }
				else
				  isLambda1520=0;
			      }
			    else
			      isLambda1520=0;
			  }
		      }
		  }
	      }
	  }
	//kine analysis*****************************************
	for(int p=0;p<knum;p++)
	  {
	    kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
	    int kineID=kine->getID();
	    int mech=kine->getMechanism();
	    int kineparentID=getMotherIndex(kine);

	    hPRmotherindex->Fill(kineparentID);

	    if(kineparentID==18)
	      hPRLambdaDoughter->Fill(kineID);

	    if(mech==0 && (kineID==2 || kineID==3))
	      {
		h2IIleptonsfromPV->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
		//if(kine->isInAcceptance())
		//h2IIleptonsInAcceptance->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	      }
	    if(kineID==14 && kineparentID==18)//proton
	      h2IIprotons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	    if(kineID==9 && kineparentID==18)//Pi-
	      h2IIpions->Fill(kine->getTotalMomentum(),kine->getThetaDeg());

	    /* no sens in new version--------------------------------------------
	       if(isDilepton==1)
	       {
	       if(kineID==9 && kineparentID==18)
	       h2PRpions->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	       if(kineID==14 && kineparentID==18)
	       h2PRprotons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	       }
	       if(isLambda==1)
	       {
	       if(kineID==2 || kineID==3)
	       h2PRleptons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	       }
	       --------------------------------------------------------------------*/
	  }
	//******************************************************

    } // end eventloop
    //**********************************************


    cout<<"number of Lambdas(1520)'s: "<<nLambdas<<endl;
    //draw alll
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

    cLambda->Divide(2,3);
    cLambda->cd(1);
    hLdistance->Draw();
    cLambda->cd(2);
    hLmass->Draw();
    cLambda->cd(3);
    hLmassDist->Draw();
    cLambda->cd(4);
    hLmassDistMass->Draw();
    cLambda->cd(5);
    hLFwDetCh2->Draw();

    cLambda1520->Divide(2,2);
    cLambda1520->cd(1);
    hL1520mass->Draw();
    cLambda1520->cd(2);
    hL1520distance->Draw();
    cLambda1520->cd(3);
    hL1520massDist->Draw();
    cLambda1520->cd(4);
    h2L1520vertex->Draw("COLZ");

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
    
    //save all
    cDiLepton->Write();
    cLambda->Write();
    cLambda1520->Write();
    cLeptonId->Write();
    cIventInfo->Write();
    cParticlesRelation->Write();
    
    output_file->Close();
    
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();
    
    return 0;
}
