#include "fwdet_res.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hvectorcand.h"
#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"

#include <TCanvas.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;


HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & beamVector, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
    if (!loop->setInput(""))
    {                                                    // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }

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

    HCategory * fCatVectorCand = nullptr;
    fCatVectorCand = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");
    if (!fCatVectorCand)
    {
        cout << "No catVectorCand!" << endl;
	//exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fCatParticleCandSim= nullptr;
    fCatParticleCandSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if(!fCatParticleCandSim)
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
    cout << "NEW ROOT TREE " << endl;
    //
    //crete histograms
    TCanvas* cMomentumDistr = new TCanvas("cMomentumDistr","Momentum distribution");
    TH1F* hMomPhiFW= new TH1F("hMomPhiFW","Phi-coordinate for momentum recorded in FW",100,1,-1);
    TH1F* hMomThetaFW= new TH1F("hMomThetaFW","Theta-coordinate for momentum recorded in FW",100,1,-1);
    TH1F* hMomPhiH= new TH1F("hMomPhiH","Phi-coordinate for momentum recorded in HADES",100,1,-1);
    TH1F* hMomThetaH= new TH1F("hMomThetaH","Theta-coordinate for momentum recorded in HADES",100,1,-1);

    TCanvas* cDistance= new TCanvas("cDistance","Distance between tracks from simulation");
    TH1F* hDistanceAll= new TH1F("hDistanceAll","Distance between all tracks",500,1,-1);
    TH1F* hDistanceCut= new TH1F("hDistanceCut","Distance between tracks after cut",40,1,-1);
    TH1F* hDistanceMassCut= new TH1F("hDistanceMassCut","Distance between tracks after cut for lambda mass",20,1,-1);
    
    TCanvas* cVertex=new TCanvas("cVertex","Vertex coordinates");
    TH1F* hVerZ=new TH1F("hVerZ","Z-coordinate of vetex",100,1,-1);
    TH1F* hVerZmassCut=new TH1F("hVerZmassCut","Z-coordinate of vetex after mass cut",40,1,-1);

    TCanvas* cMass=new TCanvas("cMass","invariant mass");
    TH1F* hMasSum=new TH1F("hMasSum","Invariant mass spektrum",500,700,2000);
    
    //event loop *************************************************
    //*********************************************************
    for (Int_t i = 0; i < entries; i++)                   
    {
        loop->nextEvent(i);         // get next event. categories will be cleared before
        if(i%5000==0)
	  cout<<"event no. "<<i<<endl;
        HParticleCandSim* particlecand =nullptr;
       	HVectorCand* fwdetstrawvec = nullptr;
	HParticleTool particle_tool;
 	//vector candidate reconstraction
	Int_t vcnt=0;
	Int_t pcnt=0;
	if (fCatVectorCand)
	  {
	    vcnt = fCatVectorCand->getEntries();
	    for (int j = 0; j < vcnt; ++j)
	      {
                fwdetstrawvec = HCategoryManager::getObject(fwdetstrawvec, fCatVectorCand, j);
		hMomPhiFW->Fill(fwdetstrawvec->getHadesPhi());
		hMomThetaFW->Fill(fwdetstrawvec->getHadesTheta());
	      }
	  }
	if(fCatParticleCandSim)
	  {
	    pcnt=fCatParticleCandSim->getEntries();
	    for(int i=0;i<pcnt;i++)
	      {
		particlecand = HCategoryManager::getObject(particlecand, fCatParticleCandSim,i);
		hMomPhiH->Fill(particlecand->getPhi()/180*TMath::Pi());
		hMomThetaH->Fill(particlecand->getTheta()/180*TMath::Pi());
	      }
	  }

	//all possible tracks' combinations from HADES and FW


	for(int j=0;j<vcnt;j++)
	  for(int i=0;i<pcnt;i++)
	    {
	      fwdetstrawvec=HCategoryManager::getObject(fwdetstrawvec, fCatVectorCand, j);
	      particlecand = HCategoryManager::getObject(particlecand, fCatParticleCandSim,i);      
	      HGeomVector base_FW;
	      base_FW.setX(fwdetstrawvec->getX());
	      base_FW.setY(fwdetstrawvec->getY());
	      base_FW.setZ(fwdetstrawvec->getZ());
	      HGeomVector dir_FW;
	      dir_FW.setX(fwdetstrawvec->getTx());
	      dir_FW.setY(fwdetstrawvec->getTy());
	      dir_FW.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
	      
	      HGeomVector base_H;
	      HGeomVector dir_H;
	      particle_tool.calcSegVector(particlecand->getZ(),particlecand->getR(),TMath::DegToRad()*particlecand->getPhi(),TMath::DegToRad()*particlecand->getTheta(),base_H,dir_H);
	      double distance=particle_tool.calculateMinimumDistance(base_FW,dir_FW,base_H,dir_H);
	      hDistanceAll->Fill(distance);
	      //distance cut
	      if(distance<50)
		{
		  hDistanceCut->Fill(distance);
		  HGeomVector vertex;
		  vertex=particle_tool.calcVertexAnalytical(base_FW,dir_FW,base_H,dir_H);
		  hVerZ->Fill(vertex.getZ());

		  fwdetstrawvec->calc4vectorProperties(938);
		  particlecand->calc4vectorProperties(140);

		  TLorentzVector sum_mass = *fwdetstrawvec + *particlecand;
		  //  sum_mass.SetPxPyPzE(fwdetstrawvec->Px()+particlecand->Px(),fwdetstrawvec->Py()+particlecand->Py(),fwdetstrawvec->Pz()+particlecand->Pz(),fwdetstrawvec->E()+particlecand->E());
		  hMasSum->Fill(sum_mass.M());
		  if(sum_mass.M()<1200 && sum_mass.M()>1070)
		    {
		      hDistanceMassCut->Fill(distance);
		      hVerZmassCut->Fill(vertex.getZ());
		    }
		}
	    }
	
    } // end eventloop
	//***********************************************************************************
	
    	




    //draw all
    cMomentumDistr->Divide(2,2);
    cMomentumDistr->cd(1);
    hMomPhiFW->Draw();
    cMomentumDistr->cd(2);
    hMomThetaFW->Draw();
    cMomentumDistr->cd(3);
    hMomPhiH->Draw();
    cMomentumDistr->cd(4);
    hMomThetaH->Draw();

    cDistance->Divide(2,2);
    cDistance->cd(1);
    hDistanceAll->Draw();
    cDistance->cd(2);
    hDistanceCut->Draw();
    cDistance->cd(3);
    hDistanceMassCut->Draw();
    
    cVertex->Divide(2);
    cVertex->cd(1);
    hVerZ->Draw();
    cVertex->cd(2);
    hVerZmassCut->Draw();

    cMass->cd();
    hMasSum->Draw();
    
    //save all
    cMomentumDistr->Write();
    cDistance->Write();
    cVertex->Write();
    cMass->Write();


    output_file->Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}
