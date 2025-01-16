#include "idaten_comm_sim_sort.C"
void gen_94Ru_comm_spectra(){
  gStyle->SetOptStat(0);
  for(int i = 0; i < 10; i++){
    teu->Add(Form("/media/jcpark/t7_shield/idaten/idaten_comm_94Ru_%d.root", 1001+i));
  }
  TCanvas* cc = new TCanvas("c94Ru", "c94Ru", 1200, 600);
  cc->Divide(2, 1);
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);
  // he_idaten->Draw();
  // // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 788, 30);

  TLatex* latex = new TLatex();
  latex->SetTextAlign(22);
  latex->SetTextSize(0.04);
  
  TH1D* ht_idaten_146_311 = new TH1D("ht_idaten_146_311", "", 400, 20, 420);  
  GetStartStopHist(teu, ht_idaten_146_311, 146, 12, 311, 15);
  cc->cd(1)->SetLogy(1);
  ht_idaten_146_311->Draw();
  TF1* f6plus = new TF1("f6plus", "[0]*0.5**(x/[1])", 50, 420);
  f6plus->SetParameter(0, ht_idaten_146_311->GetMaximum());
  f6plus->SetParameter(1, 65);
  ht_idaten_146_311->Fit(f6plus, "RN0");
  f6plus->Draw("same");
  latex->DrawLatex(200, ht_idaten_146_311->GetMaximum()*1.1, Form("T_{1/2}(6^{+}) = %.1f #pm %.1f ns (input: 65 ns)", f6plus->GetParameter(1), f6plus->GetParError(1)));
   
  TH1D* ht_idaten_311_756 = new TH1D("ht_idaten_311_756", "", 200, -2, 2);  
  GetStartStopHist(teu, ht_idaten_311_756, 311, 15, 756, 20);
  TH1D* ht_idaten_311_1431 = new TH1D("ht_idaten_311_1431", "", 200, -2, 2);  
  GetStartStopHist(teu, ht_idaten_311_1431, 311, 15, 1431, 30);
  ht_idaten_311_756->Add(ht_idaten_311_1431);
  cc->cd(2)->SetLogy(1);
  ht_idaten_311_756->Draw();
  TF1* f4plus = new TF1("f4plus", "gaus", -0.5, 0.5);
  f4plus->SetParameter(0, ht_idaten_311_756->GetMaximum());
  f4plus->SetParameter(1, 0.022);
  f4plus->SetParameter(2, 0.1);
 
  ht_idaten_311_756->Fit(f4plus, "RN0");
  f4plus->Draw("same");
  latex->DrawLatex(0, ht_idaten_311_756->GetMaximum()*1.3, Form("T_{1/2}(4^{+}) = %.1f #pm %.1f ps (input: 22 ps)", f4plus->GetParameter(1)*1000*log(2), f4plus->GetParError(1)*1000*log(2)));

  cc->SaveAs("94Ru_comm_sim_results.png");
}
