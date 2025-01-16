#include "idaten_comm_sim_sort.C"
void gen_94Pd_comm_spectra(){
  gStyle->SetOptStat(0);
  for(int i = 0; i < 11; i++){
    teu->Add(Form("/media/jcpark/t7_shield/idaten/idaten_comm_94Pd_14plus_%d.root", 1001+i));
  }
  teu->Add("/media/jcpark/t7_shield/idaten/idaten_comm_94Pd_19minus.root");
  TCanvas* cc = new TCanvas("c94Pd", "c94Pd", 1600, 900);
   cc->Divide(2, 2);
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  
  // he_idaten->Draw();
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);
  // he_idaten->Draw();
  // // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 788, 30);

  TLatex* latex = new TLatex();
  latex->SetTextAlign(22);
  TH1D* ht_idaten_1651 = new TH1D("ht_idaten_1651", "", 100, 0, 2000);  
  GetLaBr3TimeSpectrum(teu, ht_idaten_1651, 1651, 30);
  cc->cd(1)->SetLogy(1);
  ht_idaten_1651->Draw();
  TF1* f19minus = new TF1("f19minus", "[0]*0.5**(x/[1])+[2]", 20, 2000);
  f19minus->SetParameter(0, ht_idaten_1651->GetMaximum());
  f19minus->SetParameter(1, 200);
  f19minus->SetParameter(2, 2);
  
  ht_idaten_1651->Fit(f19minus, "RN0");
  f19minus->Draw("same");
  latex->DrawLatex(1000, ht_idaten_1651->GetMaximum()*0.8, Form("T_{1/2}(19^{#minus}) = %.0f #pm %.0f ns (input: 200 ns)", f19minus->GetParameter(1), f19minus->GetParError(1)));
  
  TH1D* ht_idaten_96 = new TH1D("ht_idaten_96", "", 5000, 0, 5000);  
  GetLaBr3TimeSpectrum(teu, ht_idaten_96, 96, 10);
  cc->cd(2)->SetLogy(1);
  ht_idaten_96->Draw();
  TF1* f14plus = new TF1("f14plus", "[0]*0.5**(x/[1])+[2]", 1000, 5000);
  f14plus->SetParameter(0, ht_idaten_96->GetMaximum());
  f14plus->SetParameter(1, 500);
  f14plus->SetParameter(2, 2);
  
  ht_idaten_96->Fit(f14plus, "RN0");
  f14plus->Draw("same");
  latex->DrawLatex(2500, ht_idaten_96->GetMaximum()*0.8, Form("T_{1/2}(14^{+}) = %.0f #pm %.0f ns (input: 500 ns)", f14plus->GetParameter(1), f14plus->GetParError(1)));

  TH1D* ht_idaten_1092_324 = new TH1D("ht_idaten_1092_324", "", 700, -1, 6);  
  GetStartStopHist(teu, ht_idaten_1092_324, 1092, 30, 324, 15);
  cc->cd(3)->SetLogy(1);
  ht_idaten_1092_324->Draw();
  TF1* f8plus = new TF1("f8plus", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))+[3]", -1, 6);
  f8plus->SetParameter(0, 0.755);
  f8plus->SetParameter(1, ht_idaten_1092_324->GetMaximum());
  f8plus->SetParameter(2, 0.1);
  f8plus->FixParameter(3, 0.);
  ht_idaten_1092_324->Fit(f8plus, "RN0");
  f8plus->Draw("same");
  latex->DrawLatex(3, ht_idaten_1092_324->GetMaximum()*0.8, Form("T_{1/2}(8^{+}) = %.3f #pm %.3f ns (input: 0.755 ns)", f8plus->GetParameter(0), f8plus->GetParError(0)));
	     
 
   
  TH1D* ht_idaten_324_660 = new TH1D("ht_idaten_324_660", "", 400, -2, 2);  
  GetStartStopHist(teu, ht_idaten_324_660, 324, 15, 660, 15);
  TH1D* ht_idaten_324_905 = new TH1D("ht_idaten_324_905", "", 400, -2, 2);  
  GetStartStopHist(teu, ht_idaten_324_905, 324, 15, 905, 30);
  TH1D* ht_idaten_324_814 = new TH1D("ht_idaten_324_814", "", 400, -2, 2);  
  GetStartStopHist(teu, ht_idaten_324_814, 324, 15, 814, 25);

  ht_idaten_324_660->Add(ht_idaten_324_905);
  ht_idaten_324_660->Add(ht_idaten_324_814);
  
  cc->cd(4)->SetLogy(1);
  ht_idaten_324_660->Draw();
  
  TF1* f6plus = new TF1("f6plus", "gaus", -0.3, 0.6);
  //  f6plus->SetParameter(0, 0.05);
  f6plus->SetParameter(0, ht_idaten_324_660->GetMaximum());
  f6plus->SetParameter(1, 0.1);
  f6plus->FixParameter(2, 0.2);
  ht_idaten_324_660->Fit(f6plus, "RN0");
  f6plus->Draw("same");
  latex->DrawLatex(0, ht_idaten_324_660->GetMaximum()*1.5, Form("T_{1/2}(6^{+}) = %.1f #pm %.1f ps (input: 50 ps)", f6plus->GetParameter(1)*1000.*log(2), f6plus->GetParError(1)*1000.*log(2)));

  cc->SaveAs("94Pd_comm_sim_results.png");
}
