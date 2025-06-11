#include "HSCPAnalyzer_pMSSM.h"
#include "RazorHelper.h"
#include "HSCPTree_METStudies.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"
#include "TH3F.h"
// #include "TEfficiency.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TProfile3D.h"
#include "TLorentzVector.h"
#include "TMath.h"


using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;


void HSCPAnalyzer_pMSSM::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  // Usage: 
  // ./RazorRun lists/short.txt  HSCPAnalyzer_pMSSM -f=PPStau_M-557_new2.root -d=no
 
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1


  int option = options%10;

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }

  // const float ELE_MASS = 0.000511;
  // const float MU_MASS  = 0.105658;
  // const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "test";

  }

  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "HSCP_tree.root";
  TFile *outFile;
  outFile = new TFile(outfilename.c_str(), "RECREATE");

  //======================
  //  define histograms
  //=====================

  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  //===================================
  // various other histograms desired
  //===================================

  TH1F* h_genHSCPPt = new TH1F("h_genHSCPPt","Gen HSCP Pt;Pt [GeV];",100,0,4000);
  TH1F* h_genHSCPEta = new TH1F("h_genHSCPEta","Gen HSCP Eta;Eta;",100,-6,6);
  // TH1F* h_genHSCPPhi = new TH1F("h_genHSCPPhi","Gen HSCP Phi;Phi;",100,-3.5,3.5);
  TH1F* h_genHSCPBeta = new TH1F("h_genHSCPBeta","Gen HSCP Beta;Beta;",100,0,1);
  // TH1I* h_genHSCPCharge = new TH1I("h_genHSCPCharge","Gen HSCP Charge;Charge [e];",7,-3.5,3.5);

  Int_t n_eta_bins = 100, n_pt_bins = 100, n_beta_bins = 100;
  Double_t eta_edges[101], pt_edges[101], beta_edges[101];
  
  // Create bin edges
  for(Int_t i = 0; i <= n_eta_bins; i++) eta_edges[i] = -2.5 + i * 5.0/n_eta_bins;
  for(Int_t i = 0; i <= n_pt_bins; i++) pt_edges[i] = i * 4000.0/n_pt_bins;
  for(Int_t i = 0; i <= n_beta_bins; i++) beta_edges[i] = i * 1.0/n_beta_bins;
  
  // // Create TEfficiency object
  // TEfficiency* eff_3d = new TEfficiency("eff_3d", "3D Efficiency;#eta;p_{T};#beta",
  //                                       n_eta_bins, eta_edges,
  //                                       n_pt_bins, pt_edges,
  //                                       n_beta_bins, beta_edges);

  TH3F* h_pass = new TH3F("h_pass", "Passing Events;#eta;p_{T};#beta", 
                        n_eta_bins, eta_edges, 
                        n_pt_bins, pt_edges, 
                        n_beta_bins, beta_edges);

  TH3F* h_total = new TH3F("h_total", "Total Events;#eta;p_{T};#beta", 
                          n_eta_bins, eta_edges, 
                          n_pt_bins, pt_edges, 
                          n_beta_bins, beta_edges);


  TProfile* eff_vs_pt = new TProfile("eff_vs_pt", "Efficiency vs Pt;Pt (GeV);Efficiency", 100, 0, 4000);
  TProfile* eff_vs_eta = new TProfile("eff_vs_eta", "Efficiency vs Eta;Eta;Efficiency", 50, -2.5, 2.5);
  TProfile* eff_vs_beta = new TProfile("eff_vs_beta", "Efficiency vs Beta;Beta;Efficiency", 50, 0, 1);


  TProfile3D* eff_profile = new TProfile3D("eff_profile", "Efficiency Profile;#eta;p_{T};#beta", 
                                         n_eta_bins, eta_edges, 
                                         n_pt_bins, pt_edges, 
                                         n_beta_bins, beta_edges);


  //histograms needed for trig eff plots

  // const Int_t NBINS = 20;
  // Double_t edges[NBINS + 1] = {0.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,700.0,800.0,1000.0,1200.0,1400.0,1600.0,1800.0,2000.0,2500.0,3000.0};

  
  // const Int_t NBINS2 = 22;
  // Double_t edges2[NBINS2 + 1] = {0.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,700.0,800.0,1000.0,1200.0,1400.0,1600.0,1800.0,2000.0,2500.0,3000.0,15000.0,30000.0};

 
  // const Int_t NBINS3 = 15;
  // Double_t edges3[NBINS3 + 1] = {0.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,700.0,800.0,1000.0,2500.0,5000.0};

  //====================================================
  // histograms needed for numb of events in categories
  //====================================================

  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    //begin event
    if(jentry % 10000 == 0) {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;


    auto preselection_pass = false;
    float hscp_eta = -9999.0;
    float hscp_pt = -9999.0;
    float hscp_beta = -9999.0;

    for(unsigned int i = 0; i < Pt->size(); i++) {
      if (
      Pt->at(i) > 55 
      && (eta->at(i) < 1.0 && eta->at(i) > -1.0)
      && NOPH->at(i) >= 2 && FOVH->at(i) > 0.8 
      && NOM->at(i) >= 10 && isHighPurity->at(i) == true 
      && Chi2->at(i)/Ndof->at(i) < 5 
      && dZ->at(i) < 0.1 && dXY->at(i) < 0.02 
      && PFMiniIso_relative->at(i) < 0.02 
      && EoverP->at(i) < 0.3 
      && track_genTrackMiniIsoSumPt->at(i) < 15 
      && PtErr->at(i)/(Pt->at(i)*Pt->at(i)) < 0.0008 
      && (1 - ProbQ_noL1->at(i)) > 0.3
      // Ias_StripOnly->at(i);
      ) {
        preselection_pass = true;
      }      
    }

    if (preselection_pass) {
      NEvents->Fill(1, GeneratorWeight);
      h_pass->Fill(hscp_eta, hscp_pt, hscp_beta);
    }
    h_total->Fill(hscp_eta, hscp_pt, hscp_beta);

    // gParticles
    int total_NFinalState1 = 0;
    // int HSCP_charge_ = 999;
    int index_FinalState1[10000] = {-1}; //next time vector of ints
    for(unsigned int i = 0; i < gParticlePt->size(); i++) {
      if(gParticleStatus->at(i) == 1) {
        if(total_NFinalState1 >= 10000) {
          std::cout<<"more than 10000 particles with status 1!"<<std::endl;
        } else {
          index_FinalState1[total_NFinalState1] = i;
          total_NFinalState1++;
        }
      }
    }

    // First loop: find the maximum pT HSCP particle
    float max_pt = -1;
    int max_pt_index = -1;

    for(int i = 0; i < total_NFinalState1; i++){
      //require the final state particle to have ID in the millions, i.e. require the particle to be the R-hadron
      if(fabs(gParticleId->at(index_FinalState1[i])) > 999999){  
        float HSCP_pt_ = gParticlePt->at(index_FinalState1[i]);
        
        if(HSCP_pt_ > max_pt){
            max_pt = HSCP_pt_;
            max_pt_index = i;
        }
      }	//end of if statement to check for HSCP gParticle_ID greater than a million
    } // end of loop to find max pT HSCP

    // Second part: fill histograms only for the max pT HSCP particle
    if(max_pt_index >= 0){
      hscp_pt = gParticlePt->at(index_FinalState1[max_pt_index]);
      hscp_eta = gParticleEta->at(index_FinalState1[max_pt_index]);
      hscp_beta = gParticleBeta->at(index_FinalState1[max_pt_index]);
      
      h_genHSCPPt->Fill(hscp_pt);
      h_genHSCPEta->Fill(hscp_eta);
      h_genHSCPBeta->Fill(hscp_beta);
      eff_vs_pt->Fill(hscp_pt, preselection_pass);
      eff_vs_eta->Fill(hscp_eta, preselection_pass);
      eff_vs_beta->Fill(hscp_beta, preselection_pass);
      // eff_3d->Fill(preselection_pass, hscp_eta, hscp_pt, hscp_beta);
      eff_profile->Fill(hscp_eta, hscp_pt, hscp_beta, preselection_pass);
    }
  }// end of jentry loop

  // TH3F* h_efficiency = (TH3F*)h_pass->Clone("h_efficiency");
  // h_efficiency->Divide(h_total);

  if (!isData) {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    NEvents->Write();
    h_genHSCPPt->Write();
    h_genHSCPEta->Write();
    h_genHSCPBeta->Write();  
    // eff_3d->Write();  
    eff_vs_pt->Write();
    eff_vs_eta->Write();
    eff_vs_beta->Write();
    eff_profile->Write();
    outFile->Close();
    
  } 
}
