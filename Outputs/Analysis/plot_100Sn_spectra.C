const double TMIN = -1.0;
const double TMAX = 1.0;
const int REBIN_IDATEN = 2;
void plot_100Sn_spectra(){
  gStyle->SetOptStat(0);
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  TFile* f = new TFile("outputs_100Sn.root");
  TFile* fd = new TFile("outputs_100Sn_doubleTau.root");
  TLine* line  = new TLine();
  TCanvas* cc = new TCanvas("cc", "cc", 1400, 600);
  cc->Divide(2,1);
  cc->cd(1);
  TH1D* he_idaten = (TH1D*)f->Get("he_idaten");
  TH1D* he_hpge = (TH1D*)f->Get("he_hpge");
  he_idaten->SetLineColor(1);
  he_idaten->GetXaxis()->SetTitle("Energy (keV)");
  he_idaten->GetYaxis()->SetTitle("Counts / 2 keV");
  he_idaten->GetXaxis()->SetRangeUser(0, 480);
  he_hpge->SetLineColor(2);
  he_idaten->Draw("hist");
  he_hpge->Rebin(2);
  he_hpge->Draw("Same");
  TLegend* legend = new TLegend(0.7, .7, 0.9, 0.9);
  legend->AddEntry(he_idaten, "IDATEN", "l");
  legend->AddEntry(he_hpge, "HPGe", "l");
  legend->Draw();

  latex->DrawLatex(250, 40, "^{100}Sn #rightarrow ^{100}In");
  
  cc->cd(2);
  TH1D* ht_idaten_plastic_141 = (TH1D*)f->Get("ht_idaten_plastic_141");
  ht_idaten_plastic_141->SetLineColor(1);
  ht_idaten_plastic_141->Sumw2();
  ht_idaten_plastic_141->Rebin(REBIN_IDATEN);
  ht_idaten_plastic_141->Draw("hist");
  ht_idaten_plastic_141->GetXaxis()->SetTitle("T_{IDATEN} - T_{plastic} (ns)");
  ht_idaten_plastic_141->GetYaxis()->SetTitle("Counts / 0.1 ns");
  ht_idaten_plastic_141->GetXaxis()->SetRangeUser(-1.0, 1.0);
  ht_idaten_plastic_141->GetYaxis()->SetRangeUser(0, 90);
  
  TF1* f141 = new TF1("f141", "gaus", TMIN, TMAX);
  f141->SetLineColor(1);
  ht_idaten_plastic_141->Fit(f141, "QREMN0");
  f141->Draw("same");
  TH1D* ht_idaten_plastic_95 = (TH1D*)f->Get("ht_idaten_plastic_95");
  ht_idaten_plastic_95->SetLineColor(2);
  ht_idaten_plastic_95->Sumw2();
  ht_idaten_plastic_95->Rebin(REBIN_IDATEN);
  ht_idaten_plastic_95->Draw("same hist");
  TF1* f95 = new TF1("f95", "gaus", TMIN, TMAX);
  f95->SetLineColor(2);
  ht_idaten_plastic_95->Fit(f95, "QREMN0");
  f95->Draw("same");
  TH1D* ht_idaten_plastic_40 = (TH1D*)f->Get("ht_idaten_plastic_40");
  ht_idaten_plastic_40->SetLineColor(4);
  ht_idaten_plastic_40->Sumw2();
  ht_idaten_plastic_40->Rebin(REBIN_IDATEN*2);
  ht_idaten_plastic_40->Draw("same hist");
  TF1* f40 = new TF1("f40", "gaus", TMIN, TMAX);
  f40->SetLineColor(4);
  ht_idaten_plastic_40->Fit(f40, "QREMN0");
  f40->Draw("same");

  TH1D* ht_idaten_plastic_d_141 = (TH1D*)fd->Get("ht_idaten_plastic_141");
  ht_idaten_plastic_d_141->SetLineColor(1);
  ht_idaten_plastic_d_141->Sumw2();
  ht_idaten_plastic_d_141->Rebin(REBIN_IDATEN);
  ht_idaten_plastic_d_141->SetLineStyle(2);
  ht_idaten_plastic_d_141->Draw("hist same");
  ht_idaten_plastic_d_141->GetXaxis()->SetTitle("T_{IDATEN} - T_{plastic_d} (ns)");
  ht_idaten_plastic_d_141->GetYaxis()->SetTitle("Counts / 0.1 ns");
  
  TF1* fd141 = new TF1("f141", "gaus", TMIN, TMAX);
  fd141->SetLineColor(1);
  ht_idaten_plastic_d_141->Fit(fd141, "QREMN0");
  fd141->Draw("same");
  fd141->SetLineStyle(2);
  TH1D* ht_idaten_plastic_d_95 = (TH1D*)fd->Get("ht_idaten_plastic_95");
  ht_idaten_plastic_d_95->SetLineColor(2);
  ht_idaten_plastic_d_95->Sumw2();
  ht_idaten_plastic_d_95->Rebin(REBIN_IDATEN);
  ht_idaten_plastic_d_95->SetLineStyle(2);
  ht_idaten_plastic_d_95->Draw("same hist");
  TF1* fd95 = new TF1("f95", "gaus", TMIN, TMAX);
  fd95->SetLineColor(2);
  fd95->SetLineStyle(2);
  ht_idaten_plastic_d_95->Fit(fd95, "QREMN0");
  fd95->Draw("same");
  TH1D* ht_idaten_plastic_d_40 = (TH1D*)fd->Get("ht_idaten_plastic_40");
  ht_idaten_plastic_d_40->SetLineColor(4);
  ht_idaten_plastic_d_40->Sumw2();
  ht_idaten_plastic_d_40->Rebin(REBIN_IDATEN*2);
  ht_idaten_plastic_d_40->SetLineStyle(2);
  ht_idaten_plastic_d_40->Draw("same hist");
  TF1* fd40 = new TF1("f40", "gaus", TMIN, TMAX);
  fd40->SetLineColor(4);
  fd40->SetLineStyle(2);
  ht_idaten_plastic_d_40->Fit(fd40, "QREMN0");
  fd40->Draw("same");

  TLegend* legendt = new TLegend(0.7, .7, 0.9, 0.9);
  legendt->AddEntry(ht_idaten_plastic_141, "141-keV", "l");
  legendt->AddEntry(ht_idaten_plastic_95, "95-keV", "l");
  legendt->AddEntry(ht_idaten_plastic_40, "40-keV", "l");
  legendt->Draw();
  line->DrawLine(0, 0, 0, 50);

  //  latex->DrawLatex(TMIN+0.1, 70, Form("#tau_{sim} = 7 ps, #mu = %.0f #pm %.0f ps", f141->GetParameter(1)*1000, f141->GetParError(1)*1000 ));
  latex->DrawLatex(TMIN+0.2, 85, "B(M1)_{theo} #leftrightarrow #tau_{theo}");
  latex->SetTextColor(2);
  latex->DrawLatex(TMIN+0.1, 79, Form("#tau_{sim} = 19 ps, #mu = %.0f #pm %.0f ps", f95->GetParameter(1)*1000, f95->GetParError(1)*1000 ));
  latex->SetTextColor(4);
  latex->DrawLatex(TMIN+0.1, 73, Form("#tau_{sim} = 58 ps, #mu = %.0f #pm %.0f ps", f40->GetParameter(1)*1000, f40->GetParError(1)*1000 ));
  latex->SetTextColor(1);
  latex->DrawLatex(TMIN+0.2, 67, "B(M1)_{theo}/2 #leftrightarrow 2#tau_{theo}");
  //  latex->DrawLatex(TMIN+0.1, 55, Form("#tau_{sim} = 14 ps, #mu = %.0f #pm %.0f ps", fd141->GetParameter(1)*1000, fd141->GetParError(1)*1000 ));
  latex->SetTextColor(2);

  latex->DrawLatex(TMIN+0.1, 61, Form("#tau_{sim} = 38 ps, #mu = %.0f #pm %.0f ps", fd95->GetParameter(1)*1000, fd95->GetParError(1)*1000 ));
  latex->SetTextColor(4);
  latex->DrawLatex(TMIN+0.1, 55, Form("#tau_{sim} = 115 ps, #mu = %.0f #pm %.0f ps", fd40->GetParameter(1)*1000, fd40->GetParError(1)*1000 ));

  
  cc->SaveAs("100Sn_decay_spectra_doubleTau.png");
}
