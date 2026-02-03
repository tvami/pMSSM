#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <vector>

void plot_merged_results() {
  gROOT->SetBatch(kTRUE);

  // Process specific file
  TString fileName = "merged_files/output_PMSSM_LL_2018_full_stat.root";
  std::cout << "Processing: " << fileName << std::endl;
    
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
      std::cerr << "Error opening file: " << fileName << std::endl;
      return;
    }
    
    // Get histograms from file
    TH1F* h_particle_eta = (TH1F*)file->Get("h_particle_eta");
    TH1F* h_particle_final_eta = (TH1F*)file->Get("h_particle_final_eta");
    TH1F* h_particle_pt = (TH1F*)file->Get("h_particle_pt");
    TH1F* h_particle_final_pt = (TH1F*)file->Get("h_particle_final_pt");
    TH1F* h_particle_beta = (TH1F*)file->Get("h_particle_beta");
    TH1F* h_particle_final_beta = (TH1F*)file->Get("h_particle_final_beta");
    TH1F* h_models_mother = (TH1F*)file->Get("h_models_mother");
    TH1F* h_models_mother_final = (TH1F*)file->Get("h_models_mother_final");
    TH1F* h_models_daughter = (TH1F*)file->Get("h_models_daughter");
    TH1F* h_models_daughter_final = (TH1F*)file->Get("h_models_daughter_final");
    
    // Check if histograms exist
    if (!h_particle_eta || !h_particle_final_eta || !h_particle_pt || !h_particle_final_pt ||
        !h_particle_beta || !h_particle_final_beta || !h_models_mother || !h_models_mother_final ||
        !h_models_daughter || !h_models_daughter_final) {
      std::cerr << "Error: Could not find all required histograms in " << fileName << std::endl;
      file->Close();
      return;
    }
    
    // Create canvas with 6 plots
    TCanvas* c1 = new TCanvas("c1", Form("PDG Daughter Analysis"), 1600, 1000);
    c1->Divide(3, 2);
    
    c1->cd(1);
    h_particle_eta->SetStats(false);
    h_particle_final_eta->SetStats(false);
    h_particle_eta->SetLineColor(kBlue);
    h_particle_final_eta->SetLineColor(kRed);
    double maxVal_eta = TMath::Max(h_particle_eta->GetMaximum(), h_particle_final_eta->GetMaximum());
    h_particle_eta->SetMaximum(maxVal_eta * 1.35);
    h_particle_eta->Draw();
    h_particle_final_eta->Draw("SAME");
    
    TLegend* legend1 = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetTextSize(0.04);
    legend1->AddEntry(h_particle_eta, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_eta->GetMean(), h_particle_eta->GetStdDev(), h_particle_eta->Integral()), "l");
    legend1->AddEntry(h_particle_final_eta, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_final_eta->GetMean(), h_particle_final_eta->GetStdDev(), h_particle_final_eta->Integral()), "l");
    legend1->Draw();
    
    c1->cd(2);
    h_particle_pt->SetStats(false);
    h_particle_final_pt->SetStats(false);
    h_particle_pt->SetLineColor(kBlue);
    h_particle_final_pt->SetLineColor(kRed);
    double maxVal_pt = TMath::Max(h_particle_pt->GetMaximum(), h_particle_final_pt->GetMaximum());
    h_particle_pt->SetMaximum(maxVal_pt * 1.35);
    h_particle_pt->Draw();
    h_particle_final_pt->Draw("SAME");
    
    TLegend* legend2 = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.04);
    legend2->AddEntry(h_particle_pt, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_pt->GetMean(), h_particle_pt->GetStdDev(), h_particle_pt->Integral()), "l");
    legend2->AddEntry(h_particle_final_pt, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_final_pt->GetMean(), h_particle_final_pt->GetStdDev(), h_particle_final_pt->Integral()), "l");
    legend2->Draw();
    
    c1->cd(3);
    h_particle_beta->SetStats(false);
    h_particle_final_beta->SetStats(false);
    h_particle_beta->SetLineColor(kBlue);
    h_particle_final_beta->SetLineColor(kRed);
    double maxVal_beta = TMath::Max(h_particle_beta->GetMaximum(), h_particle_final_beta->GetMaximum());
    h_particle_beta->SetMaximum(maxVal_beta * 1.35);
    h_particle_beta->Draw();
    h_particle_final_beta->Draw("SAME");
    
    TLegend* legend3 = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->SetTextSize(0.04);
    legend3->AddEntry(h_particle_beta, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_beta->GetMean(), h_particle_beta->GetStdDev(), h_particle_beta->Integral()), "l");
    legend3->AddEntry(h_particle_final_beta, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_particle_final_beta->GetMean(), h_particle_final_beta->GetStdDev(), h_particle_final_beta->Integral()), "l");
    legend3->Draw();
    
    c1->cd(4);
    h_models_mother->SetStats(false);
    h_models_mother->SetLineColor(kBlue);
    h_models_mother->SetMarkerColor(kBlue);
    h_models_mother_final->SetLineColor(kRed);
    h_models_mother_final->SetStats(false);
    h_models_mother_final->SetMarkerColor(kRed);
    
    // Set Y-axis range to accommodate both histograms
    double minVal = TMath::Min(h_models_mother->GetMinimum(), h_models_mother_final->GetMinimum());
    double maxVal = TMath::Max(h_models_mother->GetMaximum(), h_models_mother_final->GetMaximum());
    h_models_mother->SetMinimum(minVal * 0.9);
    h_models_mother->SetMaximum(maxVal * 1.35);
    
    h_models_mother->Draw("");
    h_models_mother_final->Draw("SAME");
    
    TLegend* legend4 = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->SetTextSize(0.04);
    legend4->AddEntry(h_models_mother, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_models_mother->GetMean(), h_models_mother->GetStdDev(), h_models_mother->Integral()), "l");
    legend4->AddEntry(h_models_mother_final, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_models_mother_final->GetMean(), h_models_mother_final->GetStdDev(), h_models_mother_final->Integral()), "l");
    legend4->Draw();
    
    c1->cd(6);
    h_models_daughter->SetStats(false);
    h_models_daughter_final->SetStats(false);
    h_models_daughter->SetLineColor(kBlue);
    h_models_daughter_final->SetLineColor(kRed);
    double maxVal_daughter = TMath::Max(h_models_daughter->GetMaximum(), h_models_daughter_final->GetMaximum());
    h_models_daughter->SetMaximum(maxVal_daughter * 1.35);
    h_models_daughter->Draw();
    h_models_daughter_final->Draw("SAME");
    
    TLegend* legend5 = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend5->SetBorderSize(0);
    legend5->SetFillStyle(0);
    legend5->SetTextSize(0.04);
    legend5->AddEntry(h_models_daughter, Form("Orig: Avg=%.2f, Std=%.2f, N=%.2f", h_models_daughter->GetMean(), h_models_daughter->GetStdDev(), h_models_daughter->Integral()), "l");
    legend5->AddEntry(h_models_daughter_final, Form("Final: Avg=%.2f, Std=%.2f, N=%.2f", h_models_daughter_final->GetMean(), h_models_daughter_final->GetStdDev(), h_models_daughter_final->Integral()), "l");
    legend5->Draw();
    
    // Create output filename
    TString baseName = fileName;
    baseName.ReplaceAll(".root", "");
    TString outputName = "pmssm_analysis_" + baseName + ".png";
    
    // Save canvas
    c1->SaveAs(outputName);
    std::cout << "Saved plot: " << outputName << std::endl;
    
    // Create second canvas with normalized histograms
    TCanvas* c2 = new TCanvas("c2", Form("PDG Daughter Analysis (Normalized)"), 1600, 1000);
    c2->Divide(3, 2);
    
    // Clone and normalize histograms
    TH1F* h_particle_eta_norm = (TH1F*)h_particle_eta->Clone("h_particle_eta_norm");
    TH1F* h_particle_final_eta_norm = (TH1F*)h_particle_final_eta->Clone("h_particle_final_eta_norm");
    TH1F* h_particle_pt_norm = (TH1F*)h_particle_pt->Clone("h_particle_pt_norm");
    TH1F* h_particle_final_pt_norm = (TH1F*)h_particle_final_pt->Clone("h_particle_final_pt_norm");
    TH1F* h_particle_beta_norm = (TH1F*)h_particle_beta->Clone("h_particle_beta_norm");
    TH1F* h_particle_final_beta_norm = (TH1F*)h_particle_final_beta->Clone("h_particle_final_beta_norm");
    TH1F* h_models_mother_norm = (TH1F*)h_models_mother->Clone("h_models_mother_norm");
    TH1F* h_models_mother_final_norm = (TH1F*)h_models_mother_final->Clone("h_models_mother_final_norm");
    TH1F* h_models_daughter_norm = (TH1F*)h_models_daughter->Clone("h_models_daughter_norm");
    TH1F* h_models_daughter_final_norm = (TH1F*)h_models_daughter_final->Clone("h_models_daughter_final_norm");
    
    // Normalize to unit area
    if (h_particle_eta_norm->Integral() > 0) h_particle_eta_norm->Scale(1.0 / h_particle_eta_norm->Integral());
    if (h_particle_final_eta_norm->Integral() > 0) h_particle_final_eta_norm->Scale(1.0 / h_particle_final_eta_norm->Integral());
    if (h_particle_pt_norm->Integral() > 0) h_particle_pt_norm->Scale(1.0 / h_particle_pt_norm->Integral());
    if (h_particle_final_pt_norm->Integral() > 0) h_particle_final_pt_norm->Scale(1.0 / h_particle_final_pt_norm->Integral());
    if (h_particle_beta_norm->Integral() > 0) h_particle_beta_norm->Scale(1.0 / h_particle_beta_norm->Integral());
    if (h_particle_final_beta_norm->Integral() > 0) h_particle_final_beta_norm->Scale(1.0 / h_particle_final_beta_norm->Integral());
    if (h_models_mother_norm->Integral() > 0) h_models_mother_norm->Scale(1.0 / h_models_mother_norm->Integral());
    if (h_models_mother_final_norm->Integral() > 0) h_models_mother_final_norm->Scale(1.0 / h_models_mother_final_norm->Integral());
    if (h_models_daughter_norm->Integral() > 0) h_models_daughter_norm->Scale(1.0 / h_models_daughter_norm->Integral());
    if (h_models_daughter_final_norm->Integral() > 0) h_models_daughter_final_norm->Scale(1.0 / h_models_daughter_final_norm->Integral());
    
    c2->cd(1);
    h_particle_eta_norm->SetStats(false);
    h_particle_final_eta_norm->SetStats(false);
    h_particle_eta_norm->SetLineColor(kBlue);
    h_particle_final_eta_norm->SetLineColor(kRed);
    h_particle_eta_norm->GetYaxis()->SetTitle("Normalized to Unit Area");
    double maxVal_eta_norm = TMath::Max(h_particle_eta_norm->GetMaximum(), h_particle_final_eta_norm->GetMaximum());
    h_particle_eta_norm->SetMaximum(maxVal_eta_norm * 1.35);
    h_particle_eta_norm->Draw();
    h_particle_final_eta_norm->Draw("SAME");
    
    TLegend* legend1_norm = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend1_norm->SetBorderSize(0);
    legend1_norm->SetFillStyle(0);
    legend1_norm->SetTextSize(0.04);
    legend1_norm->AddEntry(h_particle_eta_norm, Form("Orig: Avg=%.2f, Std=%.2f", h_particle_eta_norm->GetMean(), h_particle_eta_norm->GetStdDev()), "l");
    legend1_norm->AddEntry(h_particle_final_eta_norm, Form("Final: Avg=%.2f, Std=%.2f", h_particle_final_eta_norm->GetMean(), h_particle_final_eta_norm->GetStdDev()), "l");
    legend1_norm->Draw();
    
    c2->cd(2);
    h_particle_pt_norm->SetStats(false);
    h_particle_final_pt_norm->SetStats(false);
    h_particle_pt_norm->SetLineColor(kBlue);
    h_particle_final_pt_norm->SetLineColor(kRed);
    h_particle_pt_norm->GetYaxis()->SetTitle("Normalized to Unit Area");
    double maxVal_pt_norm = TMath::Max(h_particle_pt_norm->GetMaximum(), h_particle_final_pt_norm->GetMaximum());
    h_particle_pt_norm->SetMaximum(maxVal_pt_norm * 1.35);
    h_particle_pt_norm->Draw();
    h_particle_final_pt_norm->Draw("SAME");
    
    TLegend* legend2_norm = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend2_norm->SetBorderSize(0);
    legend2_norm->SetFillStyle(0);
    legend2_norm->SetTextSize(0.04);
    legend2_norm->AddEntry(h_particle_pt_norm, Form("Orig: Avg=%.2f, Std=%.2f", h_particle_pt_norm->GetMean(), h_particle_pt_norm->GetStdDev()), "l");
    legend2_norm->AddEntry(h_particle_final_pt_norm, Form("Final: Avg=%.2f, Std=%.2f", h_particle_final_pt_norm->GetMean(), h_particle_final_pt_norm->GetStdDev()), "l");
    legend2_norm->Draw();
    
    c2->cd(3);
    h_particle_beta_norm->SetStats(false);
    h_particle_final_beta_norm->SetStats(false);
    h_particle_beta_norm->SetLineColor(kBlue);
    h_particle_final_beta_norm->SetLineColor(kRed);
    h_particle_beta_norm->GetYaxis()->SetTitle("Normalized to Unit Area");
    double maxVal_beta_norm = TMath::Max(h_particle_beta_norm->GetMaximum(), h_particle_final_beta_norm->GetMaximum());
    h_particle_beta_norm->SetMaximum(maxVal_beta_norm * 1.35);
    h_particle_beta_norm->Draw();
    h_particle_final_beta_norm->Draw("SAME");
    
    TLegend* legend3_norm = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend3_norm->SetBorderSize(0);
    legend3_norm->SetFillStyle(0);
    legend3_norm->SetTextSize(0.04);
    legend3_norm->AddEntry(h_particle_beta_norm, Form("Orig: Avg=%.2f, Std=%.2f", h_particle_beta_norm->GetMean(), h_particle_beta_norm->GetStdDev()), "l");
    legend3_norm->AddEntry(h_particle_final_beta_norm, Form("Final: Avg=%.2f, Std=%.2f", h_particle_final_beta_norm->GetMean(), h_particle_final_beta_norm->GetStdDev()), "l");
    legend3_norm->Draw();
    
    c2->cd(4);
    h_models_mother_norm->SetStats(false);
    h_models_mother_norm->SetLineColor(kBlue);
    h_models_mother_norm->SetMarkerColor(kBlue);
    h_models_mother_final_norm->SetLineColor(kRed);
    h_models_mother_final_norm->SetStats(false);
    h_models_mother_final_norm->SetMarkerColor(kRed);
    h_models_mother_norm->GetYaxis()->SetTitle("Normalized to Unit Area");
    
    double minVal_norm = TMath::Min(h_models_mother_norm->GetMinimum(), h_models_mother_final_norm->GetMinimum());
    double maxVal_norm = TMath::Max(h_models_mother_norm->GetMaximum(), h_models_mother_final_norm->GetMaximum());
    h_models_mother_norm->SetMinimum(minVal_norm * 0.9);
    h_models_mother_norm->SetMaximum(maxVal_norm * 1.35);
    
    h_models_mother_norm->Draw("");
    h_models_mother_final_norm->Draw("SAME");
    
    TLegend* legend4_norm = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend4_norm->SetBorderSize(0);
    legend4_norm->SetFillStyle(0);
    legend4_norm->SetTextSize(0.04);
    legend4_norm->AddEntry(h_models_mother_norm, Form("Orig: Avg=%.2f, Std=%.2f", h_models_mother_norm->GetMean(), h_models_mother_norm->GetStdDev()), "l");
    legend4_norm->AddEntry(h_models_mother_final_norm, Form("Final: Avg=%.2f, Std=%.2f", h_models_mother_final_norm->GetMean(), h_models_mother_final_norm->GetStdDev()), "l");
    legend4_norm->Draw();
    
    c2->cd(6);
    h_models_daughter_norm->SetStats(false);
    h_models_daughter_final_norm->SetStats(false);
    h_models_daughter_norm->SetLineColor(kBlue);
    h_models_daughter_final_norm->SetLineColor(kRed);
    h_models_daughter_norm->GetYaxis()->SetTitle("Normalized to Unit Area");
    double maxVal_daughter_norm = TMath::Max(h_models_daughter_norm->GetMaximum(), h_models_daughter_final_norm->GetMaximum());
    h_models_daughter_norm->SetMaximum(maxVal_daughter_norm * 1.35);
    h_models_daughter_norm->Draw();
    h_models_daughter_final_norm->Draw("SAME");
    
    TLegend* legend5_norm = new TLegend(0.1, 0.8, 0.9, 0.9);
    legend5_norm->SetBorderSize(0);
    legend5_norm->SetFillStyle(0);
    legend5_norm->SetTextSize(0.04);
    legend5_norm->AddEntry(h_models_daughter_norm, Form("Orig: Avg=%.2f, Std=%.2f", h_models_daughter_norm->GetMean(), h_models_daughter_norm->GetStdDev()), "l");
    legend5_norm->AddEntry(h_models_daughter_final_norm, Form("Final: Avg=%.2f, Std=%.2f", h_models_daughter_final_norm->GetMean(), h_models_daughter_final_norm->GetStdDev()), "l");
    legend5_norm->Draw();
    
    // Save normalized canvas
    TString outputName_norm = outputName;
    outputName_norm.ReplaceAll(".png", "_normalized.png");
    c2->SaveAs(outputName_norm);
    std::cout << "Saved normalized plot: " << outputName_norm << std::endl;
    
    // Clean up
    delete c2;
    delete c1;
    file->Close();
    delete file;
  
  std::cout << "Processing complete!" << std::endl;
}
