TChain* tfatima;
TChain* tkhala;

const int ENERGIES[24] = {100, 200, 300, 400, 500,
			  600, 700, 800, 900, 1000,
			  1100, 1200, 1300, 1400, 1500,
			  1600, 1700, 1800, 1900, 2000,
			  2500, 3000, 3500, 4000};
const int dENERGIES[24] = {20, 20, 25, 30, 35,
			   40, 45, 50, 50, 50,
			   50, 50, 50, 50, 50,
			   50, 50, 50, 50, 50,
			   50, 50, 50, 50};

const int REBIN[24] = {1, 1, 2, 2, 2,
		       2, 2, 2, 2, 5,
		       5, 5, 5, 5, 5,
		       5, 5, 5, 5, 5,
		       10, 10, 10, 10};

TH1D* hfatima[24];
TH1D* hkhala[24];
TH1D* hidaten[24];

TF1* ffatima[24];
TF1* fkhala[24];
TF1* fidaten[24];

TGraphErrors* gfatima;
TGraphErrors* gkhala;
TGraphErrors* gidaten;


void calib_singles(){
  gStyle->SetOptStat(0);
  TCanvas* cc = new TCanvas("cc", "cc", 1800, 600);
  cc->Divide(3,1);
  tfatima = new TChain("SimulatedTree");
  tkhala = new TChain("SimulatedTree");
  gfatima = new TGraphErrors();
  gfatima->SetMarkerSize(1.3);
  gfatima->SetMarkerStyle(20);
  gfatima->SetMarkerColor(2);
  gfatima->SetLineColor(2);
  gkhala = (TGraphErrors*)gfatima->Clone("gkhala");
  gkhala->SetMarkerStyle(22);
  gkhala->SetMarkerColor(4);
  gkhala->SetLineColor(4);
  gidaten = (TGraphErrors*)gfatima->Clone("gidaten");
  gidaten->SetMarkerStyle(21);
  gidaten->SetMarkerColor(1);
  gidaten->SetLineColor(1);
  
  for(int i = 0; i < 24; i++){
    tfatima->Reset();
    tfatima->Add(Form("idaten_%dkeV.root", ENERGIES[i]));
    tkhala->Reset();
    tkhala->Add(Form("idaten_%dkeV.root", ENERGIES[i]));
    
    hfatima[i] = new TH1D(Form("hfatima%d", ENERGIES[i]), "", 4500, 0, 4500);
    hkhala[i] = new TH1D(Form("hkhala%d", ENERGIES[i]), "", 4500, 0, 4500);
    
    hfatima[i]->Rebin(REBIN[i]); 
    hkhala[i]->Rebin(REBIN[i]);
    tfatima->Project(Form("hfatima%d", ENERGIES[i]), "Fatima.fFATIMA_LaBr3_E_Energy*1000");
    tkhala->Project(Form("hkhala%d", ENERGIES[i]), "Khala.fKHALA_LaBr3_E_Energy*1000");
    hidaten[i] = (TH1D*)hfatima[i]->Clone(Form("hidaten%d", ENERGIES[i]));
    hidaten[i]->Add(hkhala[i]);
    ffatima[i] = new TF1(Form("ffatima%d", ENERGIES[i]), "gaus", ENERGIES[i]-dENERGIES[i]*2/3, ENERGIES[i]+dENERGIES[i]);
    ffatima[i]->FixParameter(1, ENERGIES[i]);
    ffatima[i]->SetParameter(0, hfatima[i]->GetMaximum());
    ffatima[i]->SetParameter(2, dENERGIES[i]);
    fkhala[i] = (TF1*)ffatima[i]->Clone(Form("fkhala%d", ENERGIES[i]));
    fidaten[i] = (TF1*)ffatima[i]->Clone(Form("fidaten%d", ENERGIES[i]));
    
    hfatima[i]->Fit(ffatima[i], "QN0R");
    hkhala[i]->Fit(fkhala[i], "QN0R");
    hidaten[i]->Fit(fidaten[i], "QN0R");
    gfatima->SetPoint(i, ENERGIES[i], ffatima[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    gkhala->SetPoint(i, ENERGIES[i], fkhala[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    gidaten->SetPoint(i, ENERGIES[i], fidaten[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    
    if (i==0){
      cc->cd(1)->SetLogy(1);
      hfatima[i]->Draw();
      hfatima[i]->GetXaxis()->SetTitle("Energy (keV)");
      hfatima[i]->GetYaxis()->SetTitle("Counts / 5 keV");
      ffatima[i]->Draw("same");
      cc->cd(2)->SetLogy(1);
      hkhala[i]->Draw();
      hkhala[i]->GetXaxis()->SetTitle("Energy (keV)");
      hkhala[i]->GetYaxis()->SetTitle("Counts / 5 keV");
      fkhala[i]->Draw("same");
      cc->cd(3)->SetLogy(1);
      hidaten[i]->Draw();
      hidaten[i]->GetXaxis()->SetTitle("Energy (keV)");
      hidaten[i]->GetYaxis()->SetTitle("Counts / 5 keV");
      fidaten[i]->Draw("same");
      
    }
    else{
      cc->cd(1);
      hfatima[i]->Draw("same");
      ffatima[i]->Draw("same");
      cc->cd(2);
      hkhala[i]->Draw("same");
      fkhala[i]->Draw("same");
      cc->cd(3);
      hidaten[i]->Draw("same");
      fidaten[i]->Draw("same");
      
    }
    
  }
		 
  TCanvas* cg = new TCanvas("cg", "cg", 1200, 800);
  cg->cd()->SetGridx(1);
  cg->cd()->SetGridy(1);
  gidaten->Draw("alp");
  gidaten->GetYaxis()->SetRangeUser(0, 0.21);
  gidaten->GetXaxis()->SetLimits(0, 4100);
  gidaten->GetXaxis()->SetTitle("Energy (keV)");
  gidaten->GetYaxis()->SetTitle("Efficiency");
  gkhala->Draw("lp");
  gfatima->Draw("lp");
  gidaten->GetXaxis()->CenterTitle();
  gidaten->GetYaxis()->CenterTitle();

  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(gidaten, "IDATEN (84 ch)", "lp");
  legend->AddEntry(gfatima, "FATIMA (36 ch)", "lp");
  legend->AddEntry(gkhala, "KHALA (48 ch)", "lp");
  legend->Draw();
  //  cs->SaveAs("idaten_138La_bgd_response.png");
}
