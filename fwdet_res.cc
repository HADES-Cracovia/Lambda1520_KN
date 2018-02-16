
#include "fwdet_res.h"
#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"
#include <TCanvas.h>
#include <TStyle.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;



HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & beamVector, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
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
    //     loop->printChain();            // print all files in the chain
    //     loop->Print();

    //     HEventHeader * header = loop->getEventHeader();
    HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");

    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    // HCategory * fFwDetStrawCal = nullptr;
    // fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");

    // if (!fFwDetStrawCal)
    // {
    //  cout << "No catFwDetStrawCal!" << endl;
    //  exit(EXIT_FAILURE);  // do you want a brute force exit ?
    // }

    HCategory * fCatVectorCand = nullptr;
    fCatVectorCand = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");

    // HCategory * fCatVectorCandSim = nullptr;
    // fCatVectorCandSim = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");

    if (!fCatVectorCand)
    {
        cout << "No catVectorCand!" << endl;
	//exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    //
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
    TCanvas* cAnglesKine=new TCanvas("cAnglesKine","Angle distribution for simulation");
    TH1F* hThetaKine=new TH1F("hThetaKine","theta for simulation",1000,1,-1);
    TH1F* hPhiKine=new TH1F("hPhiKine","phi for simulation",1000,1,-1);

    TCanvas* cAnglesRec=new TCanvas("cAnglesRec","Angle distribution for reconstruction");
    TH1F* hThetaRec=new TH1F("hThetaRec","theta for reconstruction",800,1,-1);
    TH1F* hPhiRec=new TH1F("hPhiRec","phi for reconstruction",500,1,-1);
    TH1F* hThetaRecToSim=new TH1F("hThetaRecToSim","Difference of theta between theta in vertex and at first plane of detector",1000,-1,1);
    TH1F* hPhiRecToSim=new TH1F("hPhiRecToSim","Difference of phi between vertex and at first plane of detector",1000,-1,1);

    TCanvas* cRes=new TCanvas("cRes","Resolution");
    TH1F* hThetaRes=new TH1F("hThetaRes","Resolution in theta",600,-0.1,0.1);
    TH1F* hPhiRes=new TH1F("hPhiRes","Resolution in phi",200,-0.1,0.1);
    TH1F* hXRes= new TH1F("hXRes","Resolution in X direction",200,-8,8);
    TH1F* hYRes= new TH1F("hYRes","Resolution in Y direction",200,-8,8);

    TCanvas* cResB=new TCanvas("cResB","Resolution calculated only for the same tracks");
    TH1F* hThetaResB=new TH1F("hThetaResB","Resolution in theta, calculated only for the same tracks;#theta [rad]",1000,-0.01,0.01);
    TH1F* hPhiResB=new TH1F("hPhiResB","Resolution in phi, calculated only for the same tracks",1000,-0.1,0.1);
    TH1F* hXResB= new TH1F("hXResB","Resolution in X direction, calculated only for the same tracks",400,-20,20);
    TH1F* hYResB= new TH1F("hYResB","Resolution in Y direction, calculated only for the same tracks",400,-20,20);

    TCanvas* cAngles2D=new TCanvas("cAngles2D","Simulated vs reconstracted angles");
    TH2F* h2Phi= new TH2F("h2Phi","Phi reconstructed vs. simulated; sim; rec",1000,1,-1,1000,1,-1);
    TH2F* h2Theta= new TH2F("h2Theta","Theta reconstructed vs. simulated; sim; rec",1000,0,0.5,1000,0,0.19);
    TH2F* h2PhiG= new TH2F("h2PhiG","Phi reconstructed vs. simulated in geant; recG; recDet",1000,1,-1,1000,1,-1);
    TH2F* h2ThetaG= new TH2F("h2ThetaG","Theta reconstructed vs. simulated in geant; recG; recDet",1000,0,0.5,1000,0,0.19);

    TH2F* h2XG= new TH2F("h2XG","X reconstructed vs. simulated in geant; recG; recDet",1000,1,-1,1000,1,-1);
    TH2F* h2YG= new TH2F("h2YG","Y reconstructed vs. simulated in geant; recG; recDet",1000,1,-1,1000,1,-1);
    
    TCanvas* cResGeant=new TCanvas("ResGeant", "Resolution of decector, taking into account geant simulation value in STS1");
    TH1F* hPhiG=new TH1F("hPhiG","Phi resolution reconstructed at firstv layer",400,-0.05,0.05);
    TH1F* hThetaG=new TH1F("hThetaG","Theta resolution reconstructed at first layer",400,-0.01,0.01);
    TH1F* hXG=new TH1F("hXG","X- coordinaet resolution at first layer accorging to Geant simulation",400,-2,2);
    TH1F* hYG=new TH1F("hYG","Y- coordinaet resolution at first layer accorging to Geant simulation",400,-2,2);

    TCanvas* cParameters=new TCanvas("Parameters","Some parameters");
    TH1F* hZ0=new TH1F("hZ0", "Z0 value histogram",400,3100,3200);
    TH1F* hZ0G=new TH1F("hZ0G","Z0 from geant simulation",400,3100,3200);

    TH2F* h2Xpar= new TH2F("h2Xpar","X vectorcandsim vs. simulated in Kine; geant@first layer; hKine",1000,-400,400,1000,-400,400);
    TH2F* h2Ypar= new TH2F("h2Ypar","Y vectorcandsim vs. simulated in Kine; geant@first layer; hKine",1000,-400,400,1000,-400,400);

    TCanvas* cXYcheck= new TCanvas("cXYcheck","Diffeence between x and y coordinate measured and calculated");
    TH1F* hXcheck=new TH1F("hXcheck","X- coordinate difference between calculated and taken from simulations",400,-5,5);
    TH1F* hYcheck=new TH1F("hYcheck","Y- coordinate difference between calculated and taken from simulations",400,-5,5);

    TCanvas* cVertex=new TCanvas("cVertex","Vertex coordinate");
    TH1F* hVx=new TH1F("hVx","x-coordinate",800,1,-1);
    TH1F* hVy=new TH1F("hVy","y-coordinate",800,1,-1);
    TH1F* hVz=new TH1F("hVz","z-coordinate",800,1,-1);
    
    //event loop *************************************************
    //*********************************************************
    for (Int_t i = 0; i < entries; i++)                   
    {
      if(i % 1409 ==0)
	cout<<"Event No. "<<i<<endl;
        /*Int_t nbytes =*/  loop->nextEvent(i);         // get next event. categories will be cleared before
        
        int geant_hKine_cnt = fCatGeantKine->getEntries();

        HGeantKine * gKine = nullptr;
	HVectorCandSim* fwdetstrawvec = nullptr;
	//HVectorCand* fwdetstrawvecSim= new HVectorCand;
	
        // kine
	for (int j = 0; j < geant_hKine_cnt; j++)
        {
            gKine = HCategoryManager::getObject(gKine, fCatGeantKine, j);
	    Float_t theta;
	    Float_t phi;

	    phi=gKine->getPhiDeg()/180*TMath::Pi();
	    theta=gKine->getThetaDeg()/180*TMath::Pi();
	    hThetaKine->Fill(theta);
	    hPhiKine->Fill(phi);
	 }

	//vector candidate reconstraction
	Int_t vcnt=0;
	if (fCatVectorCand)
        {
	  vcnt = fCatVectorCand->getEntries();
	    for (int j = 0; j < vcnt; ++j)
            {
	      fwdetstrawvec = (HVectorCandSim *) HCategoryManager::getObject(fwdetstrawvec, fCatVectorCand, j);
		//TVector3 mom;
		//mom=fwdetstrawvec->Vect();
	      hPhiRec->Fill(fwdetstrawvec->getHadesPhi());//getPhiH());
	      hThetaRec->Fill(fwdetstrawvec->getHadesTheta());//getThetaH());
              
            }
        }
	//resolution reconstruction
	for(int i=0; i<geant_hKine_cnt;i++)
	  for(int j=0; j<vcnt;j++)
	    {
	      double z0;//distance to STS1
	      fwdetstrawvec = (HVectorCandSim *) HCategoryManager::getObject(fwdetstrawvec, fCatVectorCand, j);
              z0=fwdetstrawvec->getZ();
	      hZ0->Fill(z0);
	      //variables reconstructed from detector
	      Float_t ty =fwdetstrawvec->getTy();
	      Float_t tx =fwdetstrawvec->getTx();
	      Float_t rr=TMath::Sqrt(tx*tx+ty*ty);
	      Float_t theta_rec=TMath::ATan(TMath::Sqrt(ty*ty+tx*tx));
	      Float_t phi_rec=TMath::ACos(tx/rr);//TMath::ATan(ty/tx);;
	      if(ty<0)
		phi_rec=2.0*TMath::Pi()-phi_rec;
	      //variables from geant simulation at first layer of STS1
	      TVector3 vectG(fwdetstrawvec->getPx1(),fwdetstrawvec->getPy1(),fwdetstrawvec->getPz1());
	      Float_t theta_recG=vectG.Theta();
	      Float_t phi_recG=vectG.Phi();//from -pi to pi
	      if(phi_recG<0)
		phi_recG=2*TMath::Pi()+phi_recG;
	      //if(tyG<0)
		//phi_recG=2.0*TMath::Pi()-phi_recG;
	      
	      Int_t track_ID_rec=fwdetstrawvec->getTrack();
		
	      gKine = HCategoryManager::getObject(gKine, fCatGeantKine,i);
	      Float_t theta_sim;
	      Float_t phi_sim;
	      Int_t track_ID_sim=gKine->getTrack();
	      phi_sim=gKine->getPhiDeg()/180*TMath::Pi();
	      theta_sim=gKine->getThetaDeg()/180*TMath::Pi();

	     
	      // fwdeststrawvecSim= HCategoryManager::getObject(fwdetstrawvecSim, fCatVectorCandSim);
	      
	 	      
	      //Theta and Phi uncertentity				
	      Float_t dif_theta;
	      Float_t dif_phi;
	      //if(theta_sim-theta_rec>TMath::Pi())
	      //	dif_theta=theta_sim-theta_rec-2*TMath::Pi();
	      // else
	      dif_theta=theta_sim-theta_rec;
	      // if(phi_sim-phi_rec>TMath::Pi())
	      //	dif_phi=phi_sim-phi_rec-2*TMath::Pi();
	      //else
	      dif_phi=phi_sim-phi_rec;
	      hThetaRes->Fill(dif_theta);		//value in rad
	      hPhiRes->Fill(dif_phi);			//value in rad
	      //X and Y uncertentity
	      if(theta_sim<TMath::Pi()/18)//only particles with theta <10 deg
		{
		  double xsim=z0*TMath::Tan(theta_sim)*TMath::Cos(phi_sim);
		  double ysim=z0*TMath::Tan(theta_sim)*TMath::Sin(phi_sim);
		  hXRes->Fill(xsim-(fwdetstrawvec->getX()));
		  hYRes->Fill(ysim-(fwdetstrawvec->getY()));
		}
	      if(track_ID_rec==track_ID_sim)
		{
		 
		  Float_t vertexX;
		  Float_t vertexY;
		  Float_t vertexZ;

		  gKine->getVertex(vertexX,vertexY,vertexZ);
		  hVx->Fill(vertexX);
		  hVy->Fill(vertexY);
		  hVz->Fill(vertexZ);

		  double xsim=(z0-vertexZ)*TMath::Tan(theta_sim)*TMath::Cos(phi_sim)+vertexX;
		  double ysim=(z0-vertexZ)*TMath::Tan(theta_sim)*TMath::Sin(phi_sim)+vertexY;
		  
		  hThetaResB->Fill(dif_theta);
		  hPhiResB->Fill(dif_phi);
		  hXResB->Fill(xsim-(fwdetstrawvec->getX()));
		  hYResB->Fill(ysim-(fwdetstrawvec->getY()));
		  hThetaRecToSim->Fill(theta_recG-theta_sim);
		  hPhiRecToSim->Fill(phi_recG-phi_sim);
		  h2Phi->Fill(phi_sim,phi_rec);
		  h2Theta->Fill(theta_sim,theta_rec);

		  //comparison between geant simulation at STS1 and signal from detector
		  hXG->Fill(fwdetstrawvec->getX()-fwdetstrawvec->getX1());
		  hYG->Fill(fwdetstrawvec->getY()-fwdetstrawvec->getY1());
		  hZ0G->Fill(fwdetstrawvec->getZ1());
		  hPhiG->Fill(phi_rec-phi_recG);
		  hThetaG->Fill(theta_rec-theta_recG);
		  h2PhiG->Fill(phi_recG,phi_rec);
		  h2ThetaG->Fill(theta_recG,theta_rec);
		  h2XG->Fill(fwdetstrawvec->getX1(),fwdetstrawvec->getX());
		  h2YG->Fill(fwdetstrawvec->getY1(),fwdetstrawvec->getY());

		  h2Xpar->Fill(fwdetstrawvec->getX1(),(fwdetstrawvec->getZ1()-vertexZ)*TMath::Tan(theta_sim)*TMath::Cos(phi_sim)+vertexX);
		  h2Ypar->Fill(fwdetstrawvec->getY1(),(fwdetstrawvec->getZ1()-vertexZ)*TMath::Tan(theta_sim)*TMath::Sin(phi_sim)+vertexY);
		 
		  hXcheck->Fill(fwdetstrawvec->getX1()-(fwdetstrawvec->getZ1()-vertexZ)*TMath::Tan(theta_sim)*TMath::Cos(phi_sim)+vertexX);
		  hYcheck->Fill(fwdetstrawvec->getY1()-(fwdetstrawvec->getZ1()-vertexZ)*TMath::Tan(theta_sim)*TMath::Sin(phi_sim)+vertexY);
		}
	    }
    } // end eventloop
	//***********************************************************************************
	
    //Drawing

    
    cAnglesKine->Divide(2);
    cAnglesKine->cd(1);
    hThetaKine->Draw();
    cAnglesKine->cd(2);
    hPhiKine->Draw();

    cAnglesRec->Divide(2,2);
    cAnglesRec->cd(1);
    hThetaRec->Draw();
    cAnglesRec->cd(2);
    hPhiRec->Draw();
    cAnglesRec->cd(3);
    hThetaRecToSim->Draw();
    cAnglesRec->cd(4);
    hPhiRecToSim->Draw();
	
    cRes->Divide(2,2);
    cRes->cd(1);
    hThetaRes->SetAxisRange(-0.15,0.15);
    hThetaRes->Draw();
    cRes->cd(2);
    hPhiRes->SetAxisRange(-2.5,2.5);
    hPhiRes->Draw();
    cRes->cd(3);
    hXRes->SetAxisRange(-200,200);
    hXRes->Draw();
    hXRes->SetAxisRange(-5,5,"X");
    cRes->cd(4);
    hYRes->SetAxisRange(-200,200);
    hYRes->SetAxisRange(-5,5,"X");
    hYRes->Draw();

    cResB->Divide(2,2);
    cResB->cd(1);
    hThetaResB->SetAxisRange(-0.15,0.15);
    hThetaResB->Draw();
    cResB->cd(2);
    hPhiResB->SetAxisRange(-2.5,2.5);
    hPhiResB->Draw();
    cResB->cd(3);
    hXResB->SetAxisRange(-200,200);
    hXResB->Draw();
    cResB->cd(4);
    hYResB->SetAxisRange(-200,200);
    hYResB->Draw();

    cAngles2D->Divide(2,2);
    cAngles2D->cd(1);
    h2Theta->Draw();
    cAngles2D->cd(2);
    h2Phi->Draw();
    cAngles2D->cd(3);
    h2PhiG->Draw();
    cAngles2D->cd(4);
    h2ThetaG->Draw();

    cResGeant->Divide(2,3);
    cResGeant->cd(1);
    hPhiG->Draw();
    cResGeant->cd(2);
    hThetaG->Draw();
    cResGeant->cd(3);
    hXG->Draw();
    cResGeant->cd(4);
    hYG->Draw();
    cResGeant->cd(5);
    h2XG->Draw();
    cResGeant->cd(6);
    h2YG->Draw();

    cParameters->Divide(2,2);
    cParameters->cd(1);
    hZ0->Draw();
    cParameters->cd(2);
    hZ0G->Draw();
    cParameters->cd(3);
    h2Xpar->Draw();
    h2Xpar->SetAxisRange(-400,400,"Y");
    h2Xpar->SetAxisRange(-400,400,"X");
    cParameters->cd(4);
    h2Ypar->Draw();
    h2Ypar->SetAxisRange(-400,400,"Y");
    h2Ypar->SetAxisRange(-400,400,"X");

    cXYcheck->Divide(2);
    cXYcheck->cd(1);
    hXcheck->Draw();
    cXYcheck->cd(2);
    hYcheck->Draw();

    cVertex->Divide(3);
    cVertex->cd(1);
    hVx->Draw();
    cVertex->cd(2);
    hVy->Draw();
    cVertex->cd(3);
    hVz->Draw();
    //End drawing
    
    //fitting
    TF1* fTheta=new TF1("fTheta","gaus(0)+pol5(3)",-1,1);
    TF1* fPhi=new TF1("fPhi","gaus(0)+pol2(3)",-1,1);
    TF1* fX=new TF1("fX","gaus(0)+pol2(3)",-20,20);
    TF1* fY=new TF1("fY","gaus(0)+pol2(3)",-20,20);

 
    fTheta->SetParameter(0,0.6*hThetaRes->GetMaximum());
    fTheta->SetParameter(1,0);
    fTheta->SetParameter(2,0.001);
    fTheta->SetRange(-0.04,0.04);
 
    fPhi->SetParameter(0,0.6*hPhiRes->GetMaximum());
    fPhi->SetParameter(1,0);
    fPhi->SetParameter(2,0.01);
    fPhi->SetRange(-0.04,0.04);

    fX->SetParameter(0,0.6*hXRes->GetMaximum());
    fX->SetParameter(1,0);
    fX->SetParameter(2,10);

    fY->SetParameter(0,0.6*hYRes->GetMaximum());
    fY->SetParameter(1,0);
    fY->SetParameter(2,10);
   
    //save canvases
    cRes->cd(1);
    hThetaRes->Fit(fTheta,"","",-0.1,0.1);
    cRes->cd(2);
    hPhiRes->Fit(fPhi,"","",-0.8,0.8);
    cRes->cd(3);
    hYRes->Fit(fY,"","",-200,200);
    cRes->cd(4);
    hXRes->Fit(fX,"","",-200,200);

    cRes->Write();
    cResB->Write();
    cAnglesRec->Write();
    cAnglesKine->Write();
    cAngles2D->Write();
    cResGeant->Write();
    cParameters->Write();
    cXYcheck->Write();
    cVertex->Write();
	
    output_file->Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}
