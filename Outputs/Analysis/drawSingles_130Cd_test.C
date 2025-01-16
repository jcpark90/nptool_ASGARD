TChain* tsim;

// time 1D histograms
TH1D* htfatima;
TH1D* htkhala;
TH1D* htidaten;

// energy 1D histograms
TH1D* hefatima;
TH1D* hekhala;
TH1D* heidaten;
TH1D* hege;

// projected energy 1D histograms
TH1D* hefatima_ge_gated;
TH1D* hekhala_ge_gated;
TH1D* heidaten_ge_gated;

// for HPGe data selection
const double EGATE_MIN = 1290;
const double EGATE_MAX = 1360;


void drawSingles_130Cd_test(){
  //  gStyle->SetOptStat(0);
  tsim = new TChain("SimulatedTree"); 
  tsim->Add("idaten_130Cd.root"); // specify input file name here

  TCanvas* cc = new TCanvas("cc", "cc", 1600, 800);
  
  cc->Divide(2,1);
  htfatima = new TH1D("htfatima", "", 1000, -5, 5); // time in ns
  htfatima->SetLineColor(1);
  htkhala = (TH1D*)htfatima->Clone("htkhala");

  hefatima = new TH1D("hefatima", "", 440, 0, 2200); // energy in keV
  hefatima->SetLineColor(1);
  hekhala = (TH1D*)hefatima->Clone("hekhala");
  
  hefatima_ge_gated = (TH1D*)hefatima->Clone("hefatima_ge_gated");
  hefatima_ge_gated->SetLineColor(4);
  hekhala_ge_gated = (TH1D*)hefatima_ge_gated->Clone("hekhala_ge_gated");

  
  hege = new TH1D("hege", "", 2500, 0, 2500); // HPGe
  hege->SetLineColor(2);

  tsim->Project("htfatima", "Fatima.fFATIMA_LaBr3_T_Time", Form("Fatima.fFATIMA_LaBr3_E_Energy > %f && Fatima.fFATIMA_LaBr3_E_Energy < %f", EGATE_MIN, EGATE_MAX));
  tsim->Project("htkhala", "Khala.fKHALA_LaBr3_T_Time", Form("Khala.fKHALA_LaBr3_E_Energy > %f && Khala.fKHALA_LaBr3_E_Energy < %f", EGATE_MIN, EGATE_MAX));

  tsim->Project("hefatima", "Fatima.fFATIMA_LaBr3_E_Energy");
  tsim->Project("hekhala", "Khala.fKHALA_LaBr3_E_Energy");
  tsim->Project("hege", "Tigress.fTIG_Ge_Energy");

  tsim->Project("hefatima_ge_gated", "Fatima.fFATIMA_LaBr3_E_Energy", Form("Tigress.fTIG_Ge_Energy > %f && Tigress.fTIG_Ge_Energy < %f", EGATE_MIN, EGATE_MAX));
  tsim->Project("hekhala_ge_gated", "Khala.fKHALA_LaBr3_E_Energy", Form("Tigress.fTIG_Ge_Energy > %f && Tigress.fTIG_Ge_Energy < %f", EGATE_MIN, EGATE_MAX));
  heidaten_ge_gated = (TH1D*)hefatima_ge_gated->Clone("heidaten_ge_gated");
  heidaten_ge_gated->Add(hekhala_ge_gated);
  
  htidaten = (TH1D*)htfatima->Clone("htidaten");
  htidaten->Add(htkhala);
  htidaten->GetXaxis()->SetTitle("Time (ns)");
  htidaten->GetYaxis()->SetTitle("Counts / 50 ns");

  heidaten = (TH1D*)hefatima->Clone("heidaten");
  heidaten->Add(hekhala);
  heidaten->GetXaxis()->SetTitle("Energy (keV)");
  heidaten->GetYaxis()->SetTitle("Counts / 5 keV");

  //  cc->cd(1)->SetLogy(1); // energy
  cc->cd(1)->SetLeftMargin(0.12);
  
  heidaten->Draw("hist");
  //  heidaten->GetYaxis()->SetRangeUser(0.5, 3e4);
  hege->Draw("same");
  //  heidaten_ge_gated->Rebin(2);
  //  heidaten_ge_gated->Draw("same");

  TLegend* legend = new TLegend(0.3, 0.7, 0.9, 0.9);
  legend->AddEntry(heidaten, "IDATEN singles", "l");
  legend->AddEntry(hege, "HPGe singles", "l");
  //  legend->AddEntry(heidaten_ge_gated, "IDATEN, gated on HPGe 141-keV #gamma", "l");
  legend->Draw();
  
  cc->cd(2)->SetLogy(1); // time
  htidaten->Draw("hist");
  TF1* fhalf = new TF1("fhalf", "[0]*0.5**(x/[1])+[2]", 0, 30000);
  fhalf->SetParameters(0, 1000);
  fhalf->SetParameters(1, 2100);
  fhalf->SetParameters(2, 40);
  htidaten->Fit(fhalf, "QRN0");
  fhalf->Draw("same");


  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->DrawLatex(6000, 1000, "^{96}Pd");
  latex->DrawLatex(6000, 600, Form("T_{1/2}(8^{+}) = %.2f #pm %.2f #mus (Sim.)", fhalf->GetParameter(1)/1000,  fhalf->GetParError(1)/1000));
  latex->DrawLatex(6000, 400, Form("T_{1/2}(8^{+}) = %.1f #mus (Geant4 ENSDF)", 2.1));
  
  cc->SaveAs("idaten_hpge_100Sn_energy_time_hists.png");
 

}
