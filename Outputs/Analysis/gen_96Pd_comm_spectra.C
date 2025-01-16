#include "idaten_comm_sim_sort.C"
void gen_96Pd_comm_spectra(){
  gStyle->SetOptStat(0);
  for(int i = 0; i <= 56; i++){
    teu->Add(Form("/media/jcpark/t7_shield/idaten/idaten_comm_96Pd_%d.root", 1001+i));
  }
  TCanvas* cc = new TCanvas("c96Pd", "c96Pd", 1800, 600);
  cc->Divide(3, 1);
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
  TH1D* ht_idaten_106 = new TH1D("ht_idaten_106_325", "", 1500, 0, 1.5e4);  
  GetLaBr3TimeSpectrum(teu, ht_idaten_106, 106, 10);
  cc->cd(1)->SetLogy(1);
  ht_idaten_106->Draw();
  TF1* f8plus = new TF1("f8plus", "[0]*0.5**(x/[1])", 0, 1.5e4);
  f8plus->SetParameter(0, ht_idaten_106->GetMaximum());
  f8plus->SetParameter(1, 1800);
  ht_idaten_106->Fit(f8plus, "RN0");
  f8plus->Draw("same");
  latex->DrawLatex(7.5e3, ht_idaten_106->GetMaximum()*0.8, Form("T_{1/2}(8^{+}) = %.0f #pm %.0f ns (input: 1800 ns)", f8plus->GetParameter(1), f8plus->GetParError(1)));
	     
  TH1D* ht_idaten_106_325 = new TH1D("ht_idaten_106_325", "", 500, -5, 45);  
  GetStartStopHist(teu, ht_idaten_106_325, 106, 10, 325, 20);
  cc->cd(2)->SetLogy(1);
  ht_idaten_106_325->Draw();
  TF1* f6plus = new TF1("f6plus", "[0]*0.5**(x/[1])", 4, 45);
  f6plus->SetParameter(0, ht_idaten_106_325->GetMaximum());
  f6plus->SetParameter(1, 7.5);
  ht_idaten_106_325->Fit(f6plus, "RN0");
  f6plus->Draw("same");
  latex->DrawLatex(20, ht_idaten_106_325->GetMaximum()*1.1, Form("T_{1/2}(6^{+}) = %.2f #pm %.2f ns (input: 7.5 ns)", f6plus->GetParameter(1), f6plus->GetParError(1)));
   
  TH1D* ht_idaten_325_684 = new TH1D("ht_idaten_325_684", "", 150, -2, 13);  
  GetStartStopHist(teu, ht_idaten_325_684, 325, 20, 684, 25);
  cc->cd(3)->SetLogy(1);
  ht_idaten_325_684->Draw();
  TF1* f4plus = new TF1("f4plus", "[0]*0.5**(x/[1])", 1, 13);
  f4plus->SetParameter(0, ht_idaten_325_684->GetMaximum());
  f4plus->SetParameter(1, 1.0);
  ht_idaten_325_684->Fit(f4plus, "RN0");
  f4plus->Draw("same");
  latex->DrawLatex(6, ht_idaten_325_684->GetMaximum()*1.1, Form("T_{1/2}(4^{+}) = %.3f #pm %.3f ns (input: 1.0 ns)", f4plus->GetParameter(1), f4plus->GetParError(1)));

  cc->SaveAs("96Pd_comm_sim_results.png");
}
