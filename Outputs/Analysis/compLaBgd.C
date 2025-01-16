TChain* tfatima[5];
TChain* tkhala[5];
TChain* thpge[5];


TH1D* hfatima_fatima[5];
TH1D* hfatima_khala[5];
TH1D* hkhala_fatima[5];
TH1D* hkhala_khala[5];

TH1D* hfatima_sum[5];
TH1D* hkhala_sum[5];
TH1D* hhpge_sum[5];
void compLaBgd(){
  gStyle->SetOptStat(0);
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  cc->cd()->SetLogy(1);
  TLegend* legend = new TLegend(0.2, 0.7, 0.9, 0.9);
  legend->SetNColumns(2);
  for(int i = 0; i < 5; i++){
    tfatima[i] = new TChain("SimulatedTree");
    tfatima[i]->Add(Form("idaten_138La_bgd_fatima_%d.root", 40+20*i));
    hfatima_fatima[i] = new TH1D(Form("hfatima_fatima%d", 40+20*i), "", 320, 0, 1600);
    hfatima_fatima[i]->SetLineColor(1);
    hfatima_fatima[i]->SetLineStyle(i+1); 
    tfatima[i]->Project(Form("hfatima_fatima%d", 40+20*i), "Fatima.fFATIMA_LaBr3_E_Energy*1000");
    if (i==0){
      hfatima_fatima[i]->Draw();
      hfatima_fatima[i]->GetXaxis()->SetTitle("Energy (keV)");
      hfatima_fatima[i]->GetYaxis()->SetTitle("Counts / 5 keV");

    }
    else{
      hfatima_fatima[i]->Draw("same");
    }
    hfatima_khala[i] = new TH1D(Form("hfatima_khala%d", 40+20*i), "", 320, 0, 1600);
    hfatima_khala[i]->SetLineColor(2);
    hfatima_khala[i]->SetLineStyle(i+1); 
    tfatima[i]->Project(Form("hfatima_khala%d", 40+20*i), "Khala.fKHALA_LaBr3_E_Energy*1000");
    hfatima_khala[i]->Draw("same");
    hfatima_sum[i] = (TH1D*)hfatima_fatima[i]->Clone(Form("hfatima_sum%d", 40+20*i));
    hfatima_sum[i]->Add(hfatima_khala[i]);
    
    tkhala[i] = new TChain("SimulatedTree");
    tkhala[i]->Add(Form("idaten_138La_bgd_khala_%d.root", 40+20*i));
    hkhala_fatima[i] = new TH1D(Form("hkhala_fatima%d", 40+20*i), "", 320, 0, 1600);
    hkhala_fatima[i]->SetLineColor(4);
    hkhala_fatima[i]->SetLineStyle(i+1); 
    tkhala[i]->Project(Form("hkhala_fatima%d", 40+20*i), "Fatima.fFATIMA_LaBr3_E_Energy*1000");
    hkhala_fatima[i]->Draw("same");
  
    hkhala_khala[i] = new TH1D(Form("hkhala_khala%d", 40+20*i), "", 320, 0, 1600);
    hkhala_khala[i]->SetLineColor(6);
    hkhala_khala[i]->SetLineStyle(i+1); 
    tkhala[i]->Project(Form("hkhala_khala%d", 40+20*i), "Khala.fKHALA_LaBr3_E_Energy*1000");
    hkhala_khala[i]->Draw("same");
    hkhala_sum[i] = (TH1D*)hkhala_fatima[i]->Clone(Form("hkhala_sum%d", 40+20*i));
    hkhala_sum[i]->Add(hkhala_khala[i]);
    
 
    hhpge_sum[i] = new TH1D(Form("hhpge_sum%d", 40+20*i), "", 1500, 0, 1500);
    if (i==0)
      tkhala[i]->Project(Form("hhpge_sum%d", 40+20*i), "Tigress.fTIG_Ge_Energy*1000");

    
  }
  legend->AddEntry(hfatima_fatima[0], "F decay, F detectors", "l");
  legend->AddEntry(hfatima_khala[0], "F decay, K detectors", "l");
  legend->AddEntry(hkhala_fatima[0], "K decay, F detectors", "l");
  legend->AddEntry(hkhala_khala[0], "K decay, K detectors", "l");
  legend->Draw();
		 
  TCanvas* cs = new TCanvas("cs", "cs", 1200, 800);
  cs->cd()->SetLogy(1);
  TLegend* legend_sum = new TLegend(0.2, 0.7, 0.9, 0.9);
  legend_sum->SetNColumns(2);
  for(int i = 0; i < 5; i++){
    if (i==0)
      hfatima_sum[i]->Draw();
    else
      hfatima_sum[i]->Draw("same");
    hkhala_sum[i]->Draw("same");
  }
  hfatima_sum[0]->GetYaxis()->SetRangeUser(50, 1e5);
  legend_sum->AddEntry(hfatima_sum[0], "Decay from FATIMA");
  legend_sum->AddEntry(hkhala_sum[0], "Decay from KHALA");
  legend_sum->Draw();
  TArrow* arrow = new TArrow();
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  arrow->SetFillColor(1);
  arrow->DrawArrow(32, 500, 32, 6000, 0.01, "|>");
  arrow->DrawArrow(789, 500, 789, 2000, 0.01, "|>");
  arrow->DrawArrow(1438, 500, 1438, 3000, 0.01, "|>");
  latex->DrawLatex(32, 400, "^{138}Ba X ray");
  latex->DrawLatex(789, 400, "^{138}Ce (+ E_{#beta})");
  latex->DrawLatex(1338, 400, "^{138}Ba (+ E_{X ray})");
  
  cs->SaveAs("idaten_138La_bgd_response.png");
  // for(int i = 0; i < 100; i++)
  //   cout<<hkhala_sum[4]->GetRandom()<<endl;
  TFile* fout = new TFile("la138_bgd.root", "recreate");
  hkhala_sum[0]->SetName("hbgd_138La");
  hkhala_sum[0]->Write();
  //  hhpge_sum[0]->SetName("hhpge_138La_bgd");
  fout->Close();

}
