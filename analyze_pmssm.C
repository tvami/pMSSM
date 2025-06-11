#include <iostream>
#include <string>
#include <vector>

void analyze_pmssm(const char* filename) {
// Usage:
// root -l
// .L analyze_pmssm.C
// analyze_pmssm("NanoGen.root", 36)
// or just
// root -l -q 'analyze_pmssm.C("PMSSM_set_1_LL_TuneCP2_13TeV_FS_NANOv9.root")'

  


  TFile* inputFile = TFile::Open("input.root", "READ");
  if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Error: Could not open input file!" << std::endl;
      return;
  }

  // Retrieve the TEfficiency object
  TProfile3D* eff_profile = (TProfile3D*)inputFile->Get("eff_profile");
  if (!eff_profile) {
      std::cerr << "Error: Could not find TProfile3D object!" << std::endl;
      inputFile->Close();
      return 1;
  }



  gROOT->SetBatch(kTRUE);
    // Open the ROOT file
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    cout << "Error: Cannot open file " << filename << endl;
    return;
  }
  
    // Get the tree
  TTree* tree = (TTree*)file->Get("Events");
  if (!tree) {
    cout << "Error: Cannot find tree in file" << endl;
    file->Close();
    return;
  }
  
    // Print tree structure
    //  cout << "Tree structure:" << endl;
    //  tree->Print();
    //  return;
  
    // Create histograms with dynamic titles
  TString title_base = Form("SUSY particles");
  TH1F* h_particle_pdgId = new TH1F("h_particle_pdgId", Form("PDG ID of %s;PDG ID;Count", title_base.Data()), 55, -15, 40);
  TH1F* h_particle_pt = new TH1F("h_particle_pt", Form("p_{T} of %s;p_{T} [GeV];Count", title_base.Data()), 100, 0, 2000);
  TH1F* h_particle_eta = new TH1F("h_particle_eta", Form("#eta of %s;#eta;Count", title_base.Data()), 100, -5, 5);
  TH1F* h_particle_beta = new TH1F("h_particle_beta", Form("#beta of %s;#beta;Count", title_base.Data()), 100, 0, 1);
  TH1F* h_particle_final_pt = new TH1F("h_particle_final_pt", Form("p_{T} of %s;p_{T} [GeV];Count", title_base.Data()), 100, 0, 2000);
  TH1F* h_particle_final_eta = new TH1F("h_particle_final_eta", Form("#eta of %s;#eta;Count", title_base.Data()), 100, -5, 5);
  TH1F* h_particle_final_beta = new TH1F("h_particle_final_beta", Form("#beta of %s;#beta;Count", title_base.Data()), 100, 0, 1);

    // Counters
  int total_events = 0;
  int events_with_mother = 0;
  int total_daughters = 0;
  int total_pairs = 0;
  
    // Get number of entries
  Long64_t nEntries = tree->GetEntries();
  cout << "Processing " << nEntries << " events..." << endl;
  int numberAllCharginos = 0;

  tree->SetScanField(0);

  // Print out the branches in the tree
  // TObjArray* branches = tree->GetListOfBranches();
  // cout << "Branches in the tree:" << endl;
  // for (int i = 0; i < branches->GetEntries(); i++) {
  //   TBranch* branch = (TBranch*)branches->At(i);
  //   cout << "  " << branch->GetName() << " (Type: " << branch->GetClassName() << ")" << endl;
  // }
  

    // Now let's manually loop through events
  // nEntries = 100;
  for (Long64_t iEvent = 0; iEvent < nEntries; iEvent++) {
    tree->GetEntry(iEvent);
    total_events++;
    if (total_events % 10000 == 0) std::cout << "We are at event " << total_events << std::endl;

    
      // Get branch values using GetLeaf
    TLeaf* leaf_nGenPart = tree->GetLeaf("nGenPart");
    TLeaf* leaf_pdgId = tree->GetLeaf("GenPart_pdgId");
    TLeaf* leaf_mother = tree->GetLeaf("GenPart_genPartIdxMother");
    TLeaf* leaf_pt = tree->GetLeaf("GenPart_pt");
    TLeaf* leaf_eta = tree->GetLeaf("GenPart_eta");
    TLeaf* leaf_phi = tree->GetLeaf("GenPart_phi");
    TLeaf* leaf_mass = tree->GetLeaf("GenPart_mass");
    TLeaf* leaf_status = tree->GetLeaf("GenPart_status");
    TLeaf* leaf_statusFlags = tree->GetLeaf("GenPart_statusFlags");

    std::vector<std::string> branchNames = {
      "GenModel_pMSSM_MCMC_519_92470.slha",
      "GenModel_pMSSM_MCMC_145_121736.slha",
      "GenModel_pMSSM_MCMC_92_80210.slha",
      "GenModel_pMSSM_MCMC_391_23728.slha",
      "GenModel_pMSSM_MCMC_225_103376.slha",
      "GenModel_pMSSM_MCMC_592_67140.slha",
      "GenModel_pMSSM_MCMC_420_13716.slha",
      "GenModel_pMSSM_MCMC_180_8666.slha",
      "GenModel_pMSSM_MCMC_187_7018.slha",
      "GenModel_pMSSM_MCMC_144_74108.slha",
      "GenModel_pMSSM_MCMC_396_36518.slha",
      "GenModel_pMSSM_MCMC_415_84583.slha",
      "GenModel_pMSSM_MCMC_227_76920.slha",
      "GenModel_pMSSM_MCMC_141_46210.slha",
      "GenModel_pMSSM_MCMC_408_46291.slha",
      "GenModel_pMSSM_MCMC_172_66754.slha",
      "GenModel_pMSSM_MCMC_286_6622.slha",
      "GenModel_pMSSM_MCMC_99_40581.slha",
      "GenModel_pMSSM_MCMC_7_109751.slha",
      "GenModel_pMSSM_MCMC_337_78455.slha",
      "GenModel_pMSSM_MCMC_454_103548.slha",
      "GenModel_pMSSM_MCMC_348_103344.slha",
      "GenModel_pMSSM_MCMC_559_20332.slha",
      "GenModel_pMSSM_MCMC_130_62470.slha",
      "GenModel_pMSSM_MCMC_409_108221.slha",
      "GenModel_pMSSM_MCMC_219_41160.slha",
      "GenModel_pMSSM_MCMC_333_56703.slha",
      "GenModel_pMSSM_MCMC_350_42125.slha",
      "GenModel_pMSSM_MCMC_496_55914.slha",
      "GenModel_pMSSM_MCMC_555_61336.slha",
      "GenModel_pMSSM_MCMC_313_111220.slha",
      "GenModel_pMSSM_MCMC_548_78340.slha",
      "GenModel_pMSSM_MCMC_595_105415.slha",
      "GenModel_pMSSM_MCMC_521_34214.slha",
      "GenModel_pMSSM_MCMC_429_32850.slha",
      "GenModel_pMSSM_MCMC_337_79479.slha"
    };

    for (const auto& branchName : branchNames) {
        TLeaf* leaf = tree->GetLeaf(branchName.c_str());
        if (leaf) {
            // Get the value of the leaf
            double value = leaf->GetValue();
            if (value != 0) { // Check if the leaf is "true" (non-zero)
                // Extract the numbers after "MCMC_"
                size_t pos = branchName.find("MCMC_");
                if (pos != std::string::npos) {
                    std::string numberStr = branchName.substr(pos + 5); // Skip "MCMC_"
                    size_t underscorePos = numberStr.find('_');
                    if (underscorePos != std::string::npos) {
                        std::string firstNumberStr = numberStr.substr(0, underscorePos);
                        std::string secondNumberStr = numberStr.substr(underscorePos + 1);
                        size_t dotPos = secondNumberStr.find('.'); // Remove ".slha" if present
                        if (dotPos != std::string::npos) {
                            secondNumberStr = secondNumberStr.substr(0, dotPos);
                        }
                        int firstNumber = std::stoi(firstNumberStr);
                        int secondNumber = std::stoi(secondNumberStr);
                        // Output both numbers
                        // std::cout << "Leaf is valid and true. First number: " << firstNumber
                                  // << ", Second number: " << secondNumber << std::endl;
                    }
                }
                break; // Exit the loop after finding the first valid leaf
            }
        } else {
            std::cerr << "Error: Leaf not found for branch: " << branchName << std::endl;
        }
    }

    if (!leaf_nGenPart || !leaf_pdgId || !leaf_mother) {
      cout << "Error: Cannot find required leaves" << endl;
      continue;
    }
    
    int nGenPart = (int)leaf_nGenPart->GetValue();
    int numChargino = 0;
    
      // Find mother particles
    for (int j = 0; j < nGenPart; j++) {
      int pdgId = (int)leaf_pdgId->GetValue(j);    
      int particle_pdg = (int)leaf_pdgId->GetValue(j);
      float particle_pt = leaf_pt ? leaf_pt->GetValue(j) : 0;
      float particle_eta = leaf_eta ? leaf_eta->GetValue(j) : 0;
      float particle_phi = leaf_phi ? leaf_phi->GetValue(j) : 0;
      float particle_mass = leaf_mass ? leaf_mass->GetValue(j) : 0;
      int particle_status = leaf_status ? (int)leaf_status->GetValue(j) : 0;
      UShort_t particle_statusFlags = leaf_statusFlags ? (UShort_t)leaf_statusFlags->GetValue(j) : 0;



      // if (abs(particle_pdg) != 2000015 || abs(particle_pdg) != 1000015) continue;
      // if (abs(particle_pdg) < 1000000) continue;
      // keep both charginos
      // if ((abs(particle_pdg) != 1000024) && (abs(particle_pdg) != 1000037)) continue;
      // keep only the first chargino
      if ((abs(particle_pdg) != 1000024)) continue;
      // if (!(particle_statusFlags & (1 << 0))) continue;
      if (!(particle_statusFlags & (1 << 13))) continue;
      if (particle_eta > 1.0 || particle_eta < -1.0) continue;

      // lets calculate the total momentum of the particle
      float particle_px = particle_pt * cos(particle_phi);
      float particle_py = particle_pt * sin(particle_phi);
      float particle_pz = particle_pt * sinh(particle_eta);
      float particle_p = sqrt(particle_px * particle_px + particle_py * particle_py + particle_pz * particle_pz);
      // let's calculate the energy of the particle
      float particle_energy = sqrt(particle_p * particle_p + particle_mass * particle_mass);
      float particle_beta = particle_p / particle_energy;

      h_particle_pt->Fill(particle_pt);
      h_particle_eta->Fill(particle_eta);
      h_particle_beta->Fill(particle_beta);


      // Get the bin corresponding to the given eta, pt, and beta
      // int binX = eff_3d->GetPassedHistogram()->GetXaxis()->FindBin(particle_eta);
      // int binY = eff_3d->GetPassedHistogram()->GetYaxis()->FindBin(particle_pt);
      // int binZ = eff_3d->GetPassedHistogram()->GetZaxis()->FindBin(particle_beta);
      // std::cout << " binX = " << binX << ", binY = " << binY << ", binZ = " << binZ << std::endl;
      // Access the efficiency for the given bin
      // double efficiency = eff_3d->GetEfficiency(eff_3d->GetGlobalBin(binX, binY, binZ));
      int bin = eff_profile->FindBin(particle_eta, particle_pt, particle_beta);
      double efficiency = eff_profile->GetBinContent(bin);

      // if (efficiency < 0.5 ) continue;

      // Print the efficiency
      std::cout << "Efficiency for eta=" << particle_eta << ", pt=" << particle_pt << ", beta=" << particle_beta
                << " is: " << efficiency << std::endl;



      numChargino++;
      numberAllCharginos++;
      int diffID = (particle_pdg > 0) ? particle_pdg  - 1000000 : particle_pdg + 1000000;
      h_particle_pdgId->Fill(diffID);
      h_particle_final_pt->Fill(particle_pt);
      h_particle_final_eta->Fill(particle_eta);
      h_particle_final_beta->Fill(particle_beta);



      
      // cout << "  Particle " << j << ": PDG=" << particle_pdg
      // << ", pT=" << particle_pt
      // << ", eta=" << particle_eta
      // << ", phi=" << particle_phi
      // << ", mass=" << particle_mass
      // << ", status=" << particle_status
      // << ", statusFlag=" << particle_statusFlags
      // << endl;
        
      } // close loop on genparticles
      std::cout << "Found " << numChargino << " charginos in event " << iEvent << std::endl;
    } // close loop on events
  
    // Print summary
  cout << "\n=== ANALYSIS SUMMARY ===" << endl;
  cout << "Total events processed: " << total_events << endl;
  cout << "Total charginos found: " << numberAllCharginos << " / " <<  nEntries << endl;
  
    // Create canvas with 5 plots
  TCanvas* c1 = new TCanvas("c1", Form("PDG Daughter Analysis"), 1600, 1000);
  c1->Divide(3, 2);
  
  c1->cd(1);
  // h_particle_pdgId->Draw();
  h_particle_eta->Draw();
  
  c1->cd(2);
  h_particle_pt->Draw();
  
  c1->cd(3);
  h_particle_beta->Draw();

  c1->cd(4);
  h_particle_final_eta->Draw();

  c1->cd(5);
  h_particle_final_pt->Draw();

  c1->cd(6);
  h_particle_final_beta->Draw();
  
    // Save canvas
  c1->SaveAs(Form("pmssm_analysis.png"));
  
    // Clean up
  file->Close();
  delete file;
  
  cout << "\nAnalysis complete!" << endl;
}