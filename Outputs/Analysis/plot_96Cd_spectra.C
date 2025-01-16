void plot_96Cd_spectra(){
  gStyle->SetOptStat(0);
  TLatex* latex = new TLatex();
  TFile* f = new TFile("outputs_96Cd.root");
  TLine* line  = new TLine();
  TCanvas* cc = new TCanvas("cc", "cc", 1400, 600);
  cc->Divide(2,1);
  cc->cd(1);
  TH1D* he_idaten = (TH1D*)f->Get("he_idaten");
  TH1D* he_hpge = (TH1D*)f->Get("he_hpge");
  he_idaten->SetLineColor(1);
  he_idaten->GetXaxis()->SetTitle("Energy (keV)");
  he_idaten->GetYaxis()->SetTitle("Counts / 2 keV");
  he_hpge->SetLineColor(2);
  he_idaten->Draw("hist");
  he_idaten->Rebin(2);
  he_idaten->GetXaxis()->SetRangeUser(0, 1600);
  he_hpge->Rebin(2);
  he_hpge->Draw("Same");
  TLegend* legend = new TLegend(0.7, .7, 0.9, 0.9);
  legend->AddEntry(he_idaten, "IDATEN", "l");
  legend->AddEntry(he_hpge, "HPGe", "l");
  legend->Draw();

  latex->DrawLatex(1200, 70, "^{96}Cd");
  latex->DrawLatex(320, 120, "307");
  latex->DrawLatex(400, 85, "419/441/457");
  latex->DrawLatex(780, 40, "811");
  latex->DrawLatex(1000, 40, "1026/1040");
  latex->DrawLatex(1100, 16, "1104");
  latex->DrawLatex(1430, 10, "1561");

  cc->cd(2);
  TH1D* ht_idaten_307_419 = (TH1D*)f->Get("ht_idaten_307_419");
  ht_idaten_307_419->SetLineColor(1);
  ht_idaten_307_419->Sumw2();
  ht_idaten_307_419->Rebin(10);
  ht_idaten_307_419->Draw("hist ");
  //ht_idaten_307_419->Draw("E1");
  ht_idaten_307_419->GetXaxis()->SetTitle("T_{stop} - T_{start} (ns)");
  ht_idaten_307_419->GetYaxis()->SetTitle("Counts / ns");
  //  ht_idaten_307_419->GetXaxis()->SetRangeUser(-1.5, 1.5);
  ht_idaten_307_419->GetYaxis()->SetRangeUser(0, 60);

  TF1* fhalf_419 = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))+[3]", -3, 20);
  fhalf_419->SetParameter(0, 0.5);
  fhalf_419->SetParameter(1, 25);
  fhalf_419->SetParameter(2, 0.1);
  fhalf_419->FixParameter(3, 0.);
  
  ht_idaten_307_419->Fit(fhalf_419, "RN0");
  fhalf_419->SetLineColor(1);
  fhalf_419->Draw("same");

  
  // TF1* f141 = new TF1("f141", "gaus", -1.5, 1.5);
  // f141->SetLineColor(1);
  // ht_idaten_307_419->Fit(f141, "REN0");
  // f141->Draw("same");
  TH1D* ht_idaten_419_441 = (TH1D*)f->Get("ht_idaten_419_441");
  ht_idaten_419_441->SetLineColor(2);
  ht_idaten_419_441->Sumw2();
  ht_idaten_419_441->Rebin(5);
  //  ht_idaten_419_441->Draw("same hist");
  // TF1* f95 = new TF1("f95", "gaus", -1.5, 1.5);
  // f95->SetLineColor(2);
  // ht_idaten_419_441->Fit(f95, "REN0");
  // f95->Draw("same");
  TH1D* ht_idaten_1104_441 = (TH1D*)f->Get("ht_idaten_1104_441");
  ht_idaten_1104_441->SetLineColor(4);
  ht_idaten_1104_441->Sumw2();
  ht_idaten_1104_441->Rebin(5);

  TH1D* ht_idaten_1104_1026_1040 = (TH1D*)f->Get("ht_idaten_1104_1026_1040");
  ht_idaten_1104_1026_1040->SetLineColor(2);
  ht_idaten_1104_1026_1040->Sumw2();
  ht_idaten_1104_1026_1040->Rebin(5);

  TH1D* ht_idaten_1104_811 = (TH1D*)f->Get("ht_idaten_1104_811");
  ht_idaten_1104_811->SetLineColor(2);
  ht_idaten_1104_811->Sumw2();
  ht_idaten_1104_811->Rebin(5);

  TH1D* ht_idaten_1561_441 = (TH1D*)f->Get("ht_idaten_1561_441");
  ht_idaten_1561_441->SetLineColor(4);
  ht_idaten_1561_441->Sumw2();
  ht_idaten_1561_441->Rebin(5);

  TH1D* ht_idaten_1561_1026_1040 = (TH1D*)f->Get("ht_idaten_1561_1026_1040");
  ht_idaten_1561_1026_1040->SetLineColor(2);
  ht_idaten_1561_1026_1040->Sumw2();
  ht_idaten_1561_1026_1040->Rebin(5);

  TH1D* ht_idaten_1561_811 = (TH1D*)f->Get("ht_idaten_1561_811");
  ht_idaten_1561_811->SetLineColor(2);
  ht_idaten_1561_811->Sumw2();
  ht_idaten_1561_811->Rebin(5);

    TH1D* ht_idaten_457_441 = (TH1D*)f->Get("ht_idaten_457_441");
  ht_idaten_457_441->SetLineColor(4);
  ht_idaten_457_441->Sumw2();
  ht_idaten_457_441->Rebin(5);

  TH1D* ht_idaten_457_1026_1040 = (TH1D*)f->Get("ht_idaten_457_1026_1040");
  ht_idaten_457_1026_1040->SetLineColor(2);
  ht_idaten_457_1026_1040->Sumw2();
  ht_idaten_457_1026_1040->Rebin(5);

  TH1D* ht_idaten_457_811 = (TH1D*)f->Get("ht_idaten_457_811");
  ht_idaten_457_811->SetLineColor(2);
  ht_idaten_457_811->Sumw2();
  ht_idaten_457_811->Rebin(5);

  
  TH1D* ht_idaten_419_1026_1040 = (TH1D*)f->Get("ht_idaten_419_1026_1040");
  ht_idaten_419_1026_1040->SetLineColor(2);
  ht_idaten_419_1026_1040->Sumw2();
  ht_idaten_419_1026_1040->Rebin(5);

  TH1D* ht_idaten_419_811 = (TH1D*)f->Get("ht_idaten_419_811");
  ht_idaten_419_811->SetLineColor(2);
  ht_idaten_419_811->Sumw2();
  ht_idaten_419_811->Rebin(5);

//  ht_idaten_1104_441->Draw("same hist");
  // TF1* f40 = new TF1("f40", "gaus", -1.5, 1.5);
  // f40->SetLineColor(4);
  // ht_idaten_1104_441->Fit(f40, "REN0");
  // f40->Draw("same");
  
  TH1D* ht_idaten_441 = (TH1D*)ht_idaten_419_441->Clone("ht_idaten_441");
  //  TH1D* ht_idaten_441 = (TH1D*)ht_idaten_419_1026_1040->Clone("ht_idaten_441");
  //TH1D* ht_idaten_441 = (TH1D*)ht_idaten_1104_1026_1040->Clone("ht_idaten_441");
  
  ht_idaten_441->Add(ht_idaten_419_1026_1040);
  ht_idaten_441->Add(ht_idaten_419_811);

  
  ht_idaten_441->Add(ht_idaten_457_441);
  ht_idaten_441->Add(ht_idaten_457_1026_1040);
  ht_idaten_441->Add(ht_idaten_457_811);
  
  ht_idaten_441->Add(ht_idaten_1104_441);
  ht_idaten_441->Add(ht_idaten_1104_1026_1040);
  ht_idaten_441->Add(ht_idaten_1104_811);

  ht_idaten_441->Add(ht_idaten_1561_441);
  ht_idaten_441->Add(ht_idaten_1561_1026_1040);
  ht_idaten_441->Add(ht_idaten_1561_811);

  
  
  ht_idaten_441->Draw("hist same");
  //  ht_idaten_441->Draw("e1 same");
  
  TF1* fhalf_441 = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))", -1, 10);
  fhalf_441->SetParameter(0, 0.5);
  fhalf_441->SetParameter(1, 25);
  fhalf_441->SetParameter(2, 0.2);
  ht_idaten_441->Fit(fhalf_441, "RN0");
  fhalf_441->SetLineColor(2);
  fhalf_441->Draw("same");
  
  TLegend* legendt = new TLegend(0.1, .7, 0.9, 0.9);
  legendt->AddEntry(ht_idaten_307_419, "307(start)-419(stop)", "l");
  // legendt->AddEntry(ht_idaten_419_441, "419-441 keV", "l");
  // legendt->AddEntry(ht_idaten_1104441, "1104-441 keV", "l");
  legendt->AddEntry(ht_idaten_441, "419/457/1104/1561(start)-441/1026/1040/811(stop)", "l");

 
  legendt->Draw();
  //  line->DrawLine(0, 0, 0, 50);
  latex->DrawLatex(2, 25, Form("T_{1/2} = %.0f #pm %.0f ps", fhalf_419->GetParameter(0)*1000., fhalf_419->GetParError(0)*1000.));
  latex->SetTextColor(2);
  latex->DrawLatex(2, 20, Form("T_{1/2} = %.0f #pm %.0f ps", fhalf_441->GetParameter(0)*1000., fhalf_441->GetParError(0)*1000.));

  cc->SaveAs("96Cd_decay_spectra.png");

  TCanvas* czoom = new TCanvas("czoom", "czoom", 1200, 800);
  TH1D* hidaten_clone = (TH1D*)he_idaten->Clone("hidaten_clone");
  hidaten_clone->Draw("hist");
  hidaten_clone->GetXaxis()->SetRangeUser(380, 500);
  TF1* fall = new TF1("fall", "gaus(0)+gaus(3)+gaus(6)+pol1(9)", 380, 500);
  fall->FixParameter(1, 419);
  fall->FixParameter(4, 441);
  fall->FixParameter(7, 457);
  fall->SetParameter(0, 30);
  fall->SetParameter(3, 30);
  fall->SetParameter(6, 15);
  fall->SetParameter(2, 8);
  fall->SetParameter(5, 8);
  fall->SetParameter(8, 8);
  fall->SetParameter(9, 3);
  hidaten_clone->Fit(fall, "QRN0");
  fall->Draw("same");
  TF1* f419 = (TF1*)fall->Clone("f419");
  f419->SetLineStyle(2);
  f419->SetParameter(3, 0);
  f419->SetParameter(6, 0);
  f419->Draw("same");
  TF1* f441 = (TF1*)fall->Clone("f441");
  f441->SetLineStyle(2);
  f441->SetParameter(0, 0);
  f441->SetParameter(6, 0);
  f441->Draw("same");
  TF1* f457 = (TF1*)fall->Clone("f457");
  f457->SetLineStyle(2);
  f457->SetParameter(3, 0);
  f457->SetParameter(0, 0);
  f457->Draw("same");
  latex->SetTextColor(1);
  latex->DrawLatex(419-5, 40, "419");
  latex->DrawLatex(441-5, 40, "441");
  latex->DrawLatex(460-5, 30, "457");
  czoom->SaveAs("96Cd_419_441_457_zoomed.png");
}
