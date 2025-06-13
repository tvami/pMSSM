#include <iostream>
#include <string>
#include <vector>
#include <THnSparse.h>

void analyze_pmssm(const char* filename) {
// Usage:
// root -l
// .L analyze_pmssm.C
// analyze_pmssm("NanoGen.root", 36)
// or selectedCharginoIndexust
// root -l -q 'analyze_pmssm.C("PMSSM_set_1_LL_TuneCP2_13TeV_FS_NANOv9.root")'

  gROOT->SetBatch(kTRUE);

  // Open the input file containing the efficiency map
  // TFile* inputFile = TFile::Open("triggerAndPres.root", "READ");
  // TFile* inputFile = TFile::Open("triggerAndPresAndSelection.root", "READ");
  TFile* inputFile = TFile::Open("triggerAndPresAndSelectionTight.root", "READ");
  if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Error: Could not open input file!" << std::endl;
      return;
  }

  // Get efficiency map from the template
  TProfile3D* eff_profile = (TProfile3D*)inputFile->Get("eff_profile");
  if (!eff_profile) {
      std::cerr << "Error: Could not find TProfile3D obselectedCharginoIndexect!" << std::endl;
      inputFile->Close();
      return 1;
  }

  // Open the input ROOT file
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    cout << "Error: Cannot open file " << filename << endl;
    return;
  }

  // Make the output root file
  TString inputFileNameStr(filename);
  TString baseName = inputFileNameStr(inputFileNameStr.Last('/') + 1, inputFileNameStr.Length());
  baseName.ReplaceAll(".root", ""); // Remove .root extension

  // Create the output file name
  TString outputFileName = baseName + "_output.root";

  // Create the output ROOT file
  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile->IsOpen()) {
      std::cerr << "Error: Could not open output file " << outputFileName << std::endl;
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
  // TH2F* h_models = new TH2F("h_models", "Models;Model 1;Model 2", 600, 0., 600., 125000, 5000.,130000.);
  TH1F* h_models = new TH1F("h_models", "Models;ID", 60, 0., 60.);
  TH1F* h_models_final = new TH1F("h_models_final", "Models;ID", 60, 0., 60.);

  TH1F* h_models_mother = new TH1F("h_models_mother", "Models;ID", 60, 0., 600.);
  TH1F* h_models_mother_final = new TH1F("h_models_mother_final", "Models;ID", 60, 0., 600.);
  TH1F* h_models_daughter = new TH1F("h_models_daughter", "Models;ID", 1448, 0., 144855.);
  TH1F* h_models_daughter_final = new TH1F("h_models_daughter_final", "Models;ID", 1448, 0., 144855.);

// Define binning for pMSSMmotherID
  const int pMSSMmotherID_nbins = 600;
  const double pMSSMmotherID_low = 0.5;
  const double pMSSMmotherID_up = 600.5;

  // Define binning for pMSSMdaughterID
  const int pMSSMdaughterID_nbins = 144855;
  const double pMSSMdaughterID_low = 0.5;
  const double pMSSMdaughterID_up = 144855.5;

  // Define binning for signal region yields (with maxEfficiency)
  const int sr_nbins = 2;
  const double sr_low = -0.5;
  const double sr_up = 1.5;

  // Create THnSparse for signal region yields
  int bins[3] = {pMSSMmotherID_nbins, pMSSMdaughterID_nbins, sr_nbins};
  double lowedges[3] = {pMSSMmotherID_low, pMSSMdaughterID_low, sr_low};
  double upedges[3] = {pMSSMmotherID_up, pMSSMdaughterID_up, sr_up};

  THnSparseD* thnSparseYields = new THnSparseD(
      "thnSparseYields",
      "",
      3,
      bins,
      lowedges,
      upedges
  );
  
    // Counters
  int total_events = 0;
  int events_with_mother = 0;
  int total_daughters = 0;
  int total_pairs = 0;
  
    // Get number of entries
  Long64_t nEntries = tree->GetEntries();
  cout << "Processing " << nEntries << " events..." << endl;
  int numberAllCharginos = 0;
  int numCharginoWithNotCorrectStatus = 0;
  int weightedEvents = 0;

  tree->SetScanField(0);

  // Print out the branches in the tree
  // TObselectedCharginoIndexArray* branches = tree->GetListOfBranches();
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

    
    const std::string prefix = "GenModel_pMSSM_MCMC_";
    std::vector<std::string> modelIDNames;

    // Get the list of branches in the tree
    TObjArray* branches = tree->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: No branches found in the tree." << std::endl;
        return;
    }

    // Iterate through all branches
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch) {
            std::string branchName = branch->GetName();
            // Check if the branch name contains the prefix
            if (branchName.find(prefix) != std::string::npos) {
                modelIDNames.push_back(branchName);
            }
        }
    }

    int pMSSMmotherID = 0;
    int pMSSMdaughterID = 0;

    for (const auto& modelIDName : modelIDNames) {
      TLeaf* leaf = tree->GetLeaf(modelIDName.c_str());
      if (!leaf) { continue; }
      else {
        // Get the value of the leaf
        double value = leaf->GetValue();
        if (value != 0) { // Check if the leaf is "true" (non-zero)
          // Extract the numbers after "MCMC_"
          size_t pos = modelIDName.find("MCMC_");
          if (pos != std::string::npos) {
            std::string numberStr = modelIDName.substr(pos + 5); // Skip "MCMC_"
            size_t underscorePos = numberStr.find('_');
            if (underscorePos != std::string::npos) {
              std::string pMSSMmotherIDStr = numberStr.substr(0, underscorePos);
              std::string pMSSMdaughterIDStr = numberStr.substr(underscorePos + 1);
              size_t dotPos = pMSSMdaughterIDStr.find('.'); // Remove ".slha" if present
              if (dotPos != std::string::npos) {
                  pMSSMdaughterIDStr = pMSSMdaughterIDStr.substr(0, dotPos);
              }
              pMSSMmotherID = std::stoi(pMSSMmotherIDStr);
              pMSSMdaughterID = std::stoi(pMSSMdaughterIDStr);
              h_models_mother->Fill(pMSSMmotherID);
              h_models_daughter->Fill(pMSSMdaughterID);
              double valuesTotalYields[3] = {static_cast<double>(pMSSMmotherID), static_cast<double>(pMSSMdaughterID), 0.};
              thnSparseYields->Fill(valuesTotalYields);
              // Setting each pair from modelIDNames into a bin in a histogram
              if (pMSSMmotherID == 519 && pMSSMdaughterID == 92470) {
                  h_models->Fill(1);
              } else if (pMSSMmotherID == 145; pMSSMdaughterID == 121736) {
                  h_models->Fill(2);
              } else if (pMSSMmotherID == 92 && pMSSMdaughterID == 80210) {
                  h_models->Fill(3);
              } else if (pMSSMmotherID == 391 && pMSSMdaughterID == 23728) {
                  h_models->Fill(4);
              } else if (pMSSMmotherID == 225 && pMSSMdaughterID == 103376) {
                  h_models->Fill(5);
              } else if (pMSSMmotherID == 592 && pMSSMdaughterID == 67140) {
                  h_models->Fill(6);
              } else if (pMSSMmotherID == 420 && pMSSMdaughterID == 13716) {
                  h_models->Fill(7);
              } else if (pMSSMmotherID == 180 && pMSSMdaughterID == 8666) {
                  h_models->Fill(8);
              } else if (pMSSMmotherID == 187 && pMSSMdaughterID == 7018) {
                  h_models->Fill(9);
              } else if (pMSSMmotherID == 144 && pMSSMdaughterID == 74108) {
                  h_models->Fill(10);
              } else if (pMSSMmotherID == 396 && pMSSMdaughterID == 36518) {
                  h_models->Fill(11);
              } else if (pMSSMmotherID == 415 && pMSSMdaughterID == 84583) {
                  h_models->Fill(12);
              } else if (pMSSMmotherID == 227 && pMSSMdaughterID == 76920) {
                  h_models->Fill(13);
              } else if (pMSSMmotherID == 141 && pMSSMdaughterID == 46210) {
                  h_models->Fill(14);
              } else if (pMSSMmotherID == 408 && pMSSMdaughterID == 46291) {
                  h_models->Fill(15);
              } else if (pMSSMmotherID == 172 && pMSSMdaughterID == 66754) {
                  h_models->Fill(16);
              } else if (pMSSMmotherID == 286 && pMSSMdaughterID == 6622) {
                  h_models->Fill(17);
              } else if (pMSSMmotherID == 99 && pMSSMdaughterID == 40581) {
                  h_models->Fill(18);
              } else if (pMSSMmotherID == 7 && pMSSMdaughterID == 109751) {
                  h_models->Fill(19);
              } else if (pMSSMmotherID == 337 && pMSSMdaughterID == 78455) {
                  h_models->Fill(20);
              } else if (pMSSMmotherID == 454 && pMSSMdaughterID == 103548) {
                  h_models->Fill(21);
              } else if (pMSSMmotherID == 348 && pMSSMdaughterID == 103344) {
                  h_models->Fill(22);
              } else if (pMSSMmotherID == 559 && pMSSMdaughterID == 20332) {
                  h_models->Fill(23);
              } else if (pMSSMmotherID == 130 && pMSSMdaughterID == 62470) {
                  h_models->Fill(24);
              } else if (pMSSMmotherID == 409 && pMSSMdaughterID == 108221) {
                  h_models->Fill(25);
              } else if (pMSSMmotherID == 219 && pMSSMdaughterID == 41160) {
                  h_models->Fill(26);
              } else if (pMSSMmotherID == 333 && pMSSMdaughterID == 56703) {
                  h_models->Fill(27);
              } else if (pMSSMmotherID == 350 && pMSSMdaughterID == 42125) {
                  h_models->Fill(28);
              } else if (pMSSMmotherID == 496 && pMSSMdaughterID == 55914) {
                  h_models->Fill(29);
              } else if (pMSSMmotherID == 555 && pMSSMdaughterID == 61336) {
                  h_models->Fill(30);
              } else if (pMSSMmotherID == 313 && pMSSMdaughterID == 111220) {
                  h_models->Fill(31);
              } else if (pMSSMmotherID == 548 && pMSSMdaughterID == 78340) {
                  h_models->Fill(32);
              } else if (pMSSMmotherID == 595 && pMSSMdaughterID == 105415) {
                  h_models->Fill(33);
              } else if (pMSSMmotherID == 521 && pMSSMdaughterID == 34214) {
                  h_models->Fill(34);
              } else if (pMSSMmotherID == 429 && pMSSMdaughterID == 32850) {
                  h_models->Fill(35);
              } else if (pMSSMmotherID == 337 && pMSSMdaughterID == 79479) {
                  h_models->Fill(36);
              } else if (pMSSMmotherID == 533 && pMSSMdaughterID == 7236) {
                  h_models->Fill(37);
              } else if (pMSSMmotherID == 175 && pMSSMdaughterID == 88241) {
                  h_models->Fill(38);
              } else if (pMSSMmotherID == 578 && pMSSMdaughterID == 17518) {
                  h_models->Fill(39);
              } else if (pMSSMmotherID == 329 && pMSSMdaughterID == 96582) {
                  h_models->Fill(40);
              } else if (pMSSMmotherID == 71 && pMSSMdaughterID == 123080) {
                  h_models->Fill(41);
              } else if (pMSSMmotherID == 128 && pMSSMdaughterID == 137591) {
                  h_models->Fill(42);
              } else if (pMSSMmotherID == 581 && pMSSMdaughterID == 66623) {
                  h_models->Fill(43);
              } else if (pMSSMmotherID == 284 && pMSSMdaughterID == 17369) {
                  h_models->Fill(44);
              } else if (pMSSMmotherID == 245 && pMSSMdaughterID == 110378) {
                  h_models->Fill(45);
              } else if (pMSSMmotherID == 434 && pMSSMdaughterID == 106293) {
                  h_models->Fill(46);
              } else if (pMSSMmotherID == 305 && pMSSMdaughterID == 49514) {
                  h_models->Fill(47);
              } else if (pMSSMmotherID == 520 && pMSSMdaughterID == 69892) {
                  h_models->Fill(48);
              } else if (pMSSMmotherID == 343 && pMSSMdaughterID == 56182) {
                  h_models->Fill(49);
              } else if (pMSSMmotherID == 310 && pMSSMdaughterID == 46904) {
                  h_models->Fill(50);
              } else if (pMSSMmotherID == 243 && pMSSMdaughterID == 118386) {
                  h_models->Fill(51);
              } else if (pMSSMmotherID == 547 && pMSSMdaughterID == 74010) {
                  h_models->Fill(52);
              } else if (pMSSMmotherID == 94 && pMSSMdaughterID == 82573) {
                  h_models->Fill(53);
              } else if (pMSSMmotherID == 94 && pMSSMdaughterID == 49337) {
                  h_models->Fill(54);
              } else if (pMSSMmotherID == 322 && pMSSMdaughterID == 125244) {
                  h_models->Fill(55);
              } else if (pMSSMmotherID == 93 && pMSSMdaughterID == 76578) {
                  h_models->Fill(56);
              } else {
                // std::cerr << "Warning: Unrecognized model ID: " << modelIDName << std::endl;
              }

              // Output both numbers
              // std::cout << "Leaf is valid and true. First number: " << pMSSMmotherID
                        // << ", Second number: " << pMSSMdaughterID << std::endl;
            }
          }
          break; // Exit the loop after finding the first valid leaf
        }
      } // end if leaf is valid and true
    } // end loop on model ID

    if (!leaf_nGenPart || !leaf_pdgId || !leaf_mother) {
      cout << "Error: Cannot find required leaves" << endl;
      continue;
    }
    
    int nGenPart = (int)leaf_nGenPart->GetValue();
    int numChargino = 0;
    int selectedCharginoIndex = -1;
    
      // Loop thru all gen particles particles
    float highestPtInEtaRange = -1; // Track the highest pT for charginos with eta in [-1, 1]
    int selectedCharginoIndexInEtaRange = -1; // Index of the highest pT chargino in eta range
    float highestPtOverall = -1; // Track the highest pT for all charginos
    int selectedCharginoIndexOverall = -1;

    for (int j = 0; j < nGenPart; j++) {
      int j_particle_pdg = leaf_pdgId->GetValue(j);
      float j_particle_pt = leaf_pt ? leaf_pt->GetValue(j) : 0;
      float j_particle_eta = leaf_eta ? leaf_eta->GetValue(j) : 0;
      int j_particle_statusFlags = leaf_statusFlags ? leaf_statusFlags->GetValue(j) : 0;

      // Keep only the first chargino
      if (abs(j_particle_pdg) != 1000024) continue;
      if (!(j_particle_statusFlags & (1 << 13))) continue;

      // Check if this chargino has the highest pT overall
      if (j_particle_pt > highestPtOverall) {
        highestPtOverall = j_particle_pt;
        selectedCharginoIndexOverall = j;
      }

      // Check if this chargino has the highest pT within eta range [-1, 1]
      if (j_particle_eta > -1 && j_particle_eta < 1) {
        if (j_particle_pt > highestPtInEtaRange) {
          highestPtInEtaRange = j_particle_pt;
          selectedCharginoIndexInEtaRange = j;
        }
      }
    } // end loop on gen particles

    if (selectedCharginoIndexInEtaRange != -1) {
      selectedCharginoIndex = selectedCharginoIndexInEtaRange;
    } else if (selectedCharginoIndexOverall != -1) {
      selectedCharginoIndex = selectedCharginoIndexOverall;
    } else {
      numCharginoWithNotCorrectStatus++;
      continue; // Skip to the next event if no chargino is found
    }

    // Loop through all gen particles again to find the selected chargino
    if (selectedCharginoIndex > 0) {
      numberAllCharginos++;
      int particle_pdg = leaf_pdgId->GetValue(selectedCharginoIndex);
      float particle_phi = leaf_phi ? leaf_phi->GetValue(selectedCharginoIndex) : 0;
      float particle_mass = leaf_mass ? leaf_mass->GetValue(selectedCharginoIndex) : 0;
      int particle_status = leaf_status ? (int)leaf_status->GetValue(selectedCharginoIndex) : 0;
      float particle_pt = leaf_pt ? leaf_pt->GetValue(selectedCharginoIndex) : 0;
      float particle_eta = leaf_eta ? leaf_eta->GetValue(selectedCharginoIndex) : 0;

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
      // binning: eta is 5 / 100 = 0.05, pt is 4000/100 = 40, beta is 1/100 = 0.01
      int binPlus = eff_profile->FindBin(particle_eta+0.0026, particle_pt+21, particle_beta+0.0051);
      int binMinus = eff_profile->FindBin(particle_eta-0.0026, particle_pt-21, particle_beta-0.0051);
      int binPlusPlus = eff_profile->FindBin(particle_eta+0.0052, particle_pt+42, particle_beta+0.0102);
      int binMinusMinus = eff_profile->FindBin(particle_eta-0.0052, particle_pt-42, particle_beta-0.0102);
      // std::cout << "Bin for eta=" << particle_eta << ", pt=" << particle_pt << ", beta=" << particle_beta
      //           << " is: " << bin << " (plus: " << binPlus << ", minus: " << binMinus << ")"
      // << std::endl;
      double efficiency = eff_profile->GetBinContent(bin);
      double efficiencyPlus1 = eff_profile->GetBinContent(binPlus);
      double efficiencyMinus1 = eff_profile->GetBinContent(binMinus);
      double efficiencyPlus2 = eff_profile->GetBinContent(binPlusPlus);
      double efficiencyMinus2 = eff_profile->GetBinContent(binMinusMinus);

      double maxEfficiency = 0;
      maxEfficiency = std::max({efficiency, efficiencyMinus1, efficiencyPlus1, efficiencyMinus2, efficiencyPlus2});
      weightedEvents += maxEfficiency;

      // if (maxEfficiency < 0.5 ) continue;

      // Print the efficiency
      // std::cout << "Efficiency for eta=" << particle_eta << ", pt=" << particle_pt << ", beta=" << particle_beta
      //           << " is: " << maxEfficiency << " (plus1: " << efficiencyPlus1 << ", minus1: " << efficiencyMinus1 << ")"
      //           << " (plus2: " << efficiencyPlus2 << ", minus2: " << efficiencyMinus2 << ")"
      //           << std::endl;


      /* 
      16065 events in total
      I loose 
      475  (3%) for not finding the chargino with status =1 
      1154 (8.4%) for not finding the chargino with eta in [-1, 1], but some of this is recovered by the next step
      9926 (68.7%) passes with some non-zero efficiency
      rest (21%) we loose bc of the efficiency cut
      */


      int diffID = (particle_pdg > 0) ? particle_pdg  - 1000000 : particle_pdg + 1000000;
      h_particle_pdgId->Fill(diffID, maxEfficiency);
      h_particle_final_pt->Fill(particle_pt, maxEfficiency);
      h_particle_final_eta->Fill(particle_eta, maxEfficiency);
      h_particle_final_beta->Fill(particle_beta, maxEfficiency);

      h_models_mother_final->Fill(pMSSMmotherID, maxEfficiency);
      h_models_daughter_final->Fill(pMSSMdaughterID, maxEfficiency);
      double valuesYields[3] = {static_cast<double>(pMSSMmotherID), static_cast<double>(pMSSMdaughterID), 1};
      thnSparseYields->Fill(valuesYields, maxEfficiency);

      if (pMSSMmotherID == 519 && pMSSMdaughterID == 92470) {
        h_models_final->Fill(1, maxEfficiency);
      } else if (pMSSMmotherID == 145; pMSSMdaughterID == 121736) {
        h_models_final->Fill(2, maxEfficiency);
      } else if (pMSSMmotherID == 92 && pMSSMdaughterID == 80210) {
        h_models_final->Fill(3, maxEfficiency);
      } else if (pMSSMmotherID == 391 && pMSSMdaughterID == 23728) {
        h_models_final->Fill(4, maxEfficiency);
      } else if (pMSSMmotherID == 225 && pMSSMdaughterID == 103376) {
        h_models_final->Fill(5, maxEfficiency);
      } else if (pMSSMmotherID == 592 && pMSSMdaughterID == 67140) {
        h_models_final->Fill(6, maxEfficiency);
      } else if (pMSSMmotherID == 420 && pMSSMdaughterID == 13716) {
        h_models_final->Fill(7, maxEfficiency);
      } else if (pMSSMmotherID == 180 && pMSSMdaughterID == 8666) {
        h_models_final->Fill(8, maxEfficiency);
      } else if (pMSSMmotherID == 187 && pMSSMdaughterID == 7018) {
        h_models_final->Fill(9, maxEfficiency);
      } else if (pMSSMmotherID == 144 && pMSSMdaughterID == 74108) {
        h_models_final->Fill(10, maxEfficiency);
      } else if (pMSSMmotherID == 396 && pMSSMdaughterID == 36518) {
        h_models_final->Fill(11, maxEfficiency);
      } else if (pMSSMmotherID == 415 && pMSSMdaughterID == 84583) {
        h_models_final->Fill(12, maxEfficiency);
      } else if (pMSSMmotherID == 227 && pMSSMdaughterID == 76920) {
        h_models_final->Fill(13, maxEfficiency);
      } else if (pMSSMmotherID == 141 && pMSSMdaughterID == 46210) {
        h_models_final->Fill(14, maxEfficiency);
      } else if (pMSSMmotherID == 408 && pMSSMdaughterID == 46291) {
        h_models_final->Fill(15, maxEfficiency);
      } else if (pMSSMmotherID == 172 && pMSSMdaughterID == 66754) {
        h_models_final->Fill(16, maxEfficiency);
      } else if (pMSSMmotherID == 286 && pMSSMdaughterID == 6622) {
        h_models_final->Fill(17, maxEfficiency);
      } else if (pMSSMmotherID == 99 && pMSSMdaughterID == 40581) {
        h_models_final->Fill(18, maxEfficiency);
      } else if (pMSSMmotherID == 7 && pMSSMdaughterID == 109751) {
        h_models_final->Fill(19, maxEfficiency);
      } else if (pMSSMmotherID == 337 && pMSSMdaughterID == 78455) {
        h_models_final->Fill(20, maxEfficiency);
      } else if (pMSSMmotherID == 454 && pMSSMdaughterID == 103548) {
        h_models_final->Fill(21, maxEfficiency);
      } else if (pMSSMmotherID == 348 && pMSSMdaughterID == 103344) {
        h_models_final->Fill(22, maxEfficiency);
      } else if (pMSSMmotherID == 559 && pMSSMdaughterID == 20332) {
        h_models_final->Fill(23, maxEfficiency);
      } else if (pMSSMmotherID == 130 && pMSSMdaughterID == 62470) {
        h_models_final->Fill(24, maxEfficiency);
      } else if (pMSSMmotherID == 409 && pMSSMdaughterID == 108221) {
        h_models_final->Fill(25, maxEfficiency);
      } else if (pMSSMmotherID == 219 && pMSSMdaughterID == 41160) {
        h_models_final->Fill(26, maxEfficiency);
      } else if (pMSSMmotherID == 333 && pMSSMdaughterID == 56703) {
        h_models_final->Fill(27, maxEfficiency);
      } else if (pMSSMmotherID == 350 && pMSSMdaughterID == 42125) {
        h_models_final->Fill(28, maxEfficiency);
      } else if (pMSSMmotherID == 496 && pMSSMdaughterID == 55914) {
        h_models_final->Fill(29, maxEfficiency);
      } else if (pMSSMmotherID == 555 && pMSSMdaughterID == 61336) {
        h_models_final->Fill(30, maxEfficiency);
      } else if (pMSSMmotherID == 313 && pMSSMdaughterID == 111220) {
        h_models_final->Fill(31, maxEfficiency);
      } else if (pMSSMmotherID == 548 && pMSSMdaughterID == 78340) {
        h_models_final->Fill(32, maxEfficiency);
      } else if (pMSSMmotherID == 595 && pMSSMdaughterID == 105415) {
        h_models_final->Fill(33, maxEfficiency);
      } else if (pMSSMmotherID == 521 && pMSSMdaughterID == 34214) {
        h_models_final->Fill(34, maxEfficiency);
      } else if (pMSSMmotherID == 429 && pMSSMdaughterID == 32850) {
        h_models_final->Fill(35, maxEfficiency);
      } else if (pMSSMmotherID == 337 && pMSSMdaughterID == 79479) {
        h_models_final->Fill(36, maxEfficiency);
      } else if (pMSSMmotherID == 533 && pMSSMdaughterID == 7236) {
        h_models_final->Fill(37, maxEfficiency);
      } else if (pMSSMmotherID == 175 && pMSSMdaughterID == 88241) {
        h_models_final->Fill(38, maxEfficiency);
      } else if (pMSSMmotherID == 578 && pMSSMdaughterID == 17518) {
        h_models_final->Fill(39, maxEfficiency);
      } else if (pMSSMmotherID == 329 && pMSSMdaughterID == 96582) {
        h_models_final->Fill(40, maxEfficiency);
      } else if (pMSSMmotherID == 71 && pMSSMdaughterID == 123080) {
        h_models_final->Fill(41, maxEfficiency);
      } else if (pMSSMmotherID == 128 && pMSSMdaughterID == 137591) {
        h_models_final->Fill(42, maxEfficiency);
      } else if (pMSSMmotherID == 581 && pMSSMdaughterID == 66623) {
        h_models_final->Fill(43, maxEfficiency);
      } else if (pMSSMmotherID == 284 && pMSSMdaughterID == 17369) {
        h_models_final->Fill(44, maxEfficiency);
      } else if (pMSSMmotherID == 245 && pMSSMdaughterID == 110378) {
        h_models_final->Fill(45, maxEfficiency);
      } else if (pMSSMmotherID == 434 && pMSSMdaughterID == 106293) {
        h_models_final->Fill(46, maxEfficiency);
      } else if (pMSSMmotherID == 305 && pMSSMdaughterID == 49514) {
        h_models_final->Fill(47, maxEfficiency);
      } else if (pMSSMmotherID == 520 && pMSSMdaughterID == 69892) {
        h_models_final->Fill(48, maxEfficiency);
      } else if (pMSSMmotherID == 343 && pMSSMdaughterID == 56182) {
        h_models_final->Fill(49, maxEfficiency);
      } else if (pMSSMmotherID == 310 && pMSSMdaughterID == 46904) {
        h_models_final->Fill(50, maxEfficiency);
      } else if (pMSSMmotherID == 243 && pMSSMdaughterID == 118386) {
        h_models_final->Fill(51, maxEfficiency);
      } else if (pMSSMmotherID == 547 && pMSSMdaughterID == 74010) {
        h_models_final->Fill(52, maxEfficiency);
      } else if (pMSSMmotherID == 94 && pMSSMdaughterID == 82573) {
        h_models_final->Fill(53, maxEfficiency);
      } else if (pMSSMmotherID == 94 && pMSSMdaughterID == 49337) {
        h_models_final->Fill(54, maxEfficiency);
      } else if (pMSSMmotherID == 322 && pMSSMdaughterID == 125244) {
        h_models_final->Fill(55, maxEfficiency);
      } else if (pMSSMmotherID == 93 && pMSSMdaughterID == 76578) {
        h_models_final->Fill(56, maxEfficiency);
      }



      // cout << "  Particle " << selectedCharginoIndex << ": PDG=" << particle_pdg
      // << ", pT=" << particle_pt
      // << ", eta=" << particle_eta
      // << ", phi=" << particle_phi
      // << ", mass=" << particle_mass
      // << ", status=" << particle_status
      // << ", statusFlag=" << particle_statusFlags
      // << endl;
        
    } // close condition on selectedChargion found 

      // std::cout << "Found " << numChargino << " charginos in event " << iEvent << std::endl;
      // if (numChargino == 0) {
      //   // Let's print out what happens in these events
      //   for (int selectedCharginoIndex = 0; selectedCharginoIndex < nGenPart; selectedCharginoIndex++) {   
      //     int particle_pdg = (int)leaf_pdgId->GetValue(selectedCharginoIndex);
      //     float particle_pt = leaf_pt ? leaf_pt->GetValue(selectedCharginoIndex) : 0;
      //     float particle_eta = leaf_eta ? leaf_eta->GetValue(selectedCharginoIndex) : 0;
      //     float particle_phi = leaf_phi ? leaf_phi->GetValue(selectedCharginoIndex) : 0;
      //     float particle_mass = leaf_mass ? leaf_mass->GetValue(selectedCharginoIndex) : 0;
      //     int particle_status = leaf_status ? (int)leaf_status->GetValue(selectedCharginoIndex) : 0;
      //     UShort_t particle_statusFlags = leaf_statusFlags ? (UShort_t)leaf_statusFlags->GetValue(selectedCharginoIndex) : 0;

      //     if ((abs(particle_pdg) != 1000024)) continue;
      //     std::cout << "Found " << numChargino << " charginos in event " << iEvent << std::endl;
      //     // if (!(particle_statusFlags & (1 << 13))) continue;
      //     if (particle_eta > 1.0 || particle_eta < -1.0) 
      //     {
      //       // std::cout << "Skipping particle with eta out of range: " << particle_eta << endl;
      //       continue;
      //     }
      //     numCharginoWithNotCorrectStatus++;
          
      //     cout << "  Particle " << selectedCharginoIndex << ": PDG=" << particle_pdg
      //       << ", pT=" << particle_pt
      //       << ", eta=" << particle_eta
      //       << ", phi=" << particle_phi
      //       << ", mass=" << particle_mass
      //       << ", status=" << particle_status
      //       << ", statusFlags=" << particle_statusFlags
      //       << " ["
      //       << (particle_statusFlags & (1 << 0) ? "isPrompt " : "")
      //       << (particle_statusFlags & (1 << 1) ? "isDecayedLeptonHadron " : "")
      //       << (particle_statusFlags & (1 << 2) ? "isTauDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 3) ? "isPromptTauDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 4) ? "isDirectTauDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 5) ? "isDirectPromptTauDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 6) ? "isDirectHadronDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 7) ? "isHardProcess " : "")
      //       << (particle_statusFlags & (1 << 8) ? "fromHardProcess " : "")
      //       << (particle_statusFlags & (1 << 9) ? "isHardProcessTauDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 10) ? "isDirectHardProcessTauDecayProduct " : "")
      //       << (particle_statusFlags & (1 << 11) ? "fromHardProcessBeforeFSR " : "")
      //       << (particle_statusFlags & (1 << 12) ? "isFirstCopy " : "")
      //       << (particle_statusFlags & (1 << 13) ? "isLastCopy " : "")
      //       << (particle_statusFlags & (1 << 14) ? "isLastCopyBeforeFSR " : "")
      //       << "]"
      //       << endl;
      //   }

      // } // close on if numCHargino == 0
  } // close loop on events
  
    // Print summary
  cout << "\n=== ANALYSIS SUMMARY ===" << endl;
  cout << "Total events processed: " << total_events << endl;
  cout << "Total charginos found: " << numberAllCharginos << " / " <<  nEntries << endl;
  cout << "Weighted events: " << weightedEvents << endl;
  cout << "Total charginos with incorrect status: " << numCharginoWithNotCorrectStatus << endl;
  

  
    // Create canvas with 5 plots
  TCanvas* c1 = new TCanvas("c1", Form("PDG Daughter Analysis"), 1600, 1000);
  c1->Divide(3, 2);
  
  c1->cd(1);
  h_particle_eta->SetStats(false); // Disable stats box for this histogram
  h_particle_final_eta->SetStats(false);
  h_particle_eta->Draw();
  h_particle_eta->SetLineColor(kBlue);
  h_particle_final_eta->SetLineColor(kRed);
  h_particle_final_eta->Draw("SAME");

  // Create a legend for the first plot
  TLegend* legend1 = new TLegend(0.1, 0.8, 0.9, 0.9); // Position: (x1, y1, x2, y2)
  legend1->SetBorderSize(0); // No border
  legend1->SetFillStyle(0);  // Transparent background
  legend1->SetTextSize(0.04); 
  legend1->AddEntry(h_particle_eta, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_eta->GetMean(), h_particle_eta->GetStdDev(), h_particle_eta->Integral()), "l");
  legend1->AddEntry(h_particle_final_eta, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_final_eta->GetMean(), h_particle_final_eta->GetStdDev(), h_particle_final_eta->Integral()), "l");
  legend1->Draw();

  c1->cd(2);
  h_particle_pt->SetStats(false); 
  h_particle_final_pt->SetStats(false);
  h_particle_pt->Draw();
  h_particle_final_pt->SetLineColor(kRed);
  h_particle_final_pt->Draw("SAME");

  // Create a legend for the second plot
  TLegend* legend2 = new TLegend(0.1, 0.8, 0.9, 0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->SetTextSize(0.04); 
  legend2->AddEntry(h_particle_pt, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_pt->GetMean(), h_particle_pt->GetStdDev(), h_particle_pt->Integral()), "l");
  legend2->AddEntry(h_particle_final_pt, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_final_pt->GetMean(), h_particle_final_pt->GetStdDev(), h_particle_final_pt->Integral()), "l");
  legend2->Draw();

  c1->cd(3);
  h_particle_beta->Draw();
  h_particle_beta->SetStats(false);
  h_particle_final_beta->SetStats(false);
  h_particle_final_beta->SetLineColor(kRed);
  h_particle_final_beta->Draw("SAME");

  // Create a legend for the third plot
  TLegend* legend3 = new TLegend(0.1, 0.8, 0.9, 0.9);
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->SetTextSize(0.04); 
  legend3->AddEntry(h_particle_beta, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_beta->GetMean(), h_particle_beta->GetStdDev(), h_particle_beta->Integral()), "l");
  legend3->AddEntry(h_particle_final_beta, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_final_beta->GetMean(), h_particle_final_beta->GetStdDev(), h_particle_final_beta->Integral()), "l");
  legend3->Draw();

  c1->cd(4);
  h_models_mother->Draw();
  h_models_mother->SetStats(false);
  h_models_mother_final->SetStats(false);
  h_models_mother_final->SetLineColor(kRed);
  h_models_mother_final->Draw("SAME");

  // // Add the legend here too
  // TLegend* legend4 = new TLegend(0.1, 0.8, 0.9, 0.9);
  // legend4->SetBorderSize(0);
  // legend4->SetFillStyle(0);
  // legend4->SetTextSize(0.04);
  // legend4->AddEntry(h_models_mother, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_models->GetMean(), h_models->GetStdDev(), h_models->Integral()), "l");
  // legend4->AddEntry(h_models_mother_final, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_models_final->GetMean(), h_models_final->GetStdDev(), h_models_final->Integral()), "l");
  // legend4->Draw();

  // Let's repeat what we have in cd(4) but now normalize the histograms to unit aread
  // c1->cd(5);
  // TH1F* h_models_norm = (TH1F*)h_models_mother->Clone("h_models_mother_norm");
  // TH1F* h_models_final_norm = (TH1F*)h_models_mother_final->Clone("h_models_mother_final_norm");
  // h_models_norm->Scale(1.0 / h_models_norm->Integral());
  // h_models_final_norm->Scale(1.0 / h_models_final_norm->Integral());
  // h_models_norm->Draw();
  // h_models_norm->SetStats(false);
  // h_models_final_norm->SetStats(false);
  // h_models_final_norm->SetLineColor(kRed);
  // h_models_final_norm->Draw("SAME");

  c1->cd(6);
  h_models_daughter->Draw();
  h_models_daughter->SetStats(false);
  h_models_daughter_final->SetStats(false);
  h_models_daughter_final->SetLineColor(kRed);
  h_models_daughter_final->Draw("SAME");
  TLegend* legend5 = new TLegend(0.1, 0.8, 0.9, 0.9);
  legend5->SetBorderSize(0);
  legend5->SetFillStyle(0);
  legend5->SetTextSize(0.04);
  legend5->AddEntry(h_models_daughter, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_models_daughter->GetMean(), h_models_daughter->GetStdDev(), h_models_daughter->Integral()), "l");
  legend5->AddEntry(h_models_daughter_final, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_models_daughter_final->GetMean(), h_models_daughter_final->GetStdDev(), h_models_daughter_final->Integral()), "l");
  legend5->Draw();


  // Save canvas
  c1->SaveAs(Form("pmssm_analysis.png"));
  thnSparseYields->Write();
  h_particle_pdgId->Write();
  h_particle_final_pt->Write();
  h_particle_final_eta->Write();
  h_particle_final_beta->Write();
  h_particle_pt->Write();
  h_particle_eta->Write();
  h_particle_beta->Write();
  h_models_mother->Write();
  h_models_mother_final->Write();
  h_models_daughter->Write();
  h_models_daughter_final->Write();
  h_models->Write();
  h_models_final->Write();

  
    // Clean up
  file->Close();
  delete file;

  outputFile->Close();
  delete outputFile;
  
  cout << "\nAnalysis complete!" << endl;
}