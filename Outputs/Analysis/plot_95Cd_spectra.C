

void plot_95Cd_spectra(){
  gStyle->SetOptStat(0);
  TLatex* latex = new TLatex();
  TFile* f332 = new TFile("outputs_95Cd_332_realistic.root");
  TLine* line  = new TLine();
  TCanvas* cc332 = new TCanvas("cc332", "cc332", 1200, 800);
  //  cc332->cd()->SetLogy(1);
  TH1D* he_idaten332 = (TH1D*)f332->Get("he_idaten");
  TH1D* he_hpge332 = (TH1D*)f332->Get("he_hpge");
  he_idaten332->SetLineColor(1);
  he_idaten332->GetXaxis()->SetTitle("Energy (keV)");
  he_idaten332->GetYaxis()->SetTitle("Counts / 10 keV (IDATEN), 5 keV (HPGe)");
  he_hpge332->SetLineColor(2);
  he_idaten332->Draw("hist");
  he_idaten332->Rebin(10);
  //  he_idaten->GetXaxis()->SetRangeUser(0, 1600);
  he_hpge332->Rebin(5);
  he_hpge332->Draw("Same");
  TLegend* legend = new TLegend(0.7, .7, 0.9, 0.9);
  legend->AddEntry(he_idaten332, "IDATEN", "l");
  legend->AddEntry(he_hpge332, "HPGe", "l");
  legend->Draw();

  latex->DrawLatex(350, 17, "^{95}Cd (33/2^{+}), T_{1/2} = 40 #mus");
  latex->DrawLatex(350, 15, "[20-#mus correlation window]");
  
  latex->DrawLatex(175, 16, "159");
  latex->SetTextAlign(22);
  latex->DrawLatex(875, 11, "875");
  latex->DrawLatex(1294, 8, "1294");
  cc332->SaveAs("95Cd_332_isomer_spectra.png");
  
  TFile* f232 = new TFile("outputs_95Cd_232.root");
  TCanvas* cc232 = new TCanvas("cc232", "cc232", 1200, 800);
  cc232->cd()->SetLogy(1);
  TH1D* he_idaten232 = (TH1D*)f232->Get("he_idaten");
  TH1D* he_hpge232 = (TH1D*)f232->Get("he_hpge");
  he_idaten232->SetLineColor(1);
  he_idaten232->GetXaxis()->SetTitle("Energy (keV)");
  he_idaten232->GetYaxis()->SetTitle("Counts / 10 keV (IDATEN), 5 keV (HPGe)");
  he_hpge232->SetLineColor(2);
  he_idaten232->Draw("hist");
  he_idaten232->Rebin(10);
  //  he_idaten->GetXaxis()->SetRangeUser(0, 1600);
  he_hpge232->Rebin(5);
  he_hpge232->Draw("Same");
  legend->Draw();

  // latex->DrawLatex(600, 600, "^{95}Cd (23/2^{+}), T_{1/2} = 2.1 ms");
  // latex->SetTextAlign(22);
  // latex->DrawLatex(164, 700, "164");
  // latex->DrawLatex(428, 250, "428");
  // latex->SetTextAngle(90);
  // latex->DrawLatex(823, 50, "823");
  // latex->DrawLatex(936, 60, "936");
  // latex->DrawLatex(1003, 55, "1003");
  // latex->DrawLatex(1117, 45, "1117");
  // cc232->SaveAs("95Cd_232_isomer_spectra.png");

  latex->DrawLatex(600, 180, "^{95}Cd (23/2^{+}), T_{1/2} = 2.1 ms");
  latex->DrawLatex(650, 100, "[IC events of 164 or 428]");
  
  latex->SetTextAlign(22);
  latex->DrawLatex(164, 200, "164");
  latex->DrawLatex(428, 60, "428");
  latex->SetTextAngle(90);
  latex->DrawLatex(823, 14, "823");
  latex->DrawLatex(936, 15, "936");
  latex->DrawLatex(1003, 10, "1003");
  latex->DrawLatex(1117, 10, "1117");
  cc232->SaveAs("95Cd_232_isomer_spectra.png");

  TFile* f12 = new TFile("outputs_95Cd_12.root");
  TCanvas* cc12 = new TCanvas("cc12", "cc12", 1200, 800);
  //  cc12->cd()->SetLogy(1);
  TH1D* he_idaten12 = (TH1D*)f12->Get("he_idaten");
  TH1D* he_hpge12 = (TH1D*)f12->Get("he_hpge");
  he_idaten12->SetLineColor(1);
  he_idaten12->GetXaxis()->SetTitle("Energy (keV)");
  he_idaten12->GetYaxis()->SetTitle("Counts / 5 keV (IDATEN), 5 keV (HPGe)");
  he_hpge12->SetLineColor(2);
  he_idaten12->Draw("hist");
  he_idaten12->Rebin(5);
  he_idaten12->GetYaxis()->SetRangeUser(0, 18);
  he_idaten12->GetXaxis()->SetRangeUser(0, 400);
  he_hpge12->Rebin(5);
  he_hpge12->Draw("Same");
  legend->Draw();

  latex->SetTextAngle(0);
  latex->DrawLatex(150, 16, "^{95}Cd (1/2^{-}), T_{1/2} = 90 ms");
  latex->DrawLatex(150, 14, "[IC events of 77 or 267]");
  
  latex->DrawLatex(77, 5, "77");
  latex->DrawLatex(267, 14, "267");
  cc12->SaveAs("95Cd_12_isomer_spectra.png");
}
