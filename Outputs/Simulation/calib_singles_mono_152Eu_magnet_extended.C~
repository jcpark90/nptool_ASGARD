TChain* tfatima;
TChain* tkhala;
TChain* ttigress;


TChain* teu; // tree for 152Eu data
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TFatimaData.h"
#include "TTigressData.h"

const int ENERGIES[29] = {30, 40, 50, 60, 70, 80, 90,
			  100, 150, 200, 250, 300, 400, 500,
			  600, 700, 800, 900, 1000,
			  1100, 1200, 1300, 1400, 1500,
			  1600, 1700, 1800, 1900, 2000};

const int dENERGIES[29] = {10, 10, 10, 15, 15, 15, 15,
			   20, 20, 20, 25, 25, 30, 35,
			   40, 45, 45, 50, 50,
			   50, 50, 50, 50, 50,
			   50, 50, 50, 50, 50};

const int REBIN[29] = {1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 2, 2, 2, 2,
		       2, 2, 2, 2, 5,
		       5, 5, 5, 5, 5,
		       5, 5, 5, 5, 5};

const double E152EU_IDATEN[7] = {244, 344,  964,  1408, 778, 122, 40};
const double dE152EU_IDATEN[7] = {50, 50, 50, 50, 50, 20, 20};
const double I152EU_IDATEN[7] = {7.5, 26.5, 14.5, 21, 12.942,   28.58, 21.1+38.3+0.252+0.454+3.74+7.24+0.107+2.39+0.91+0.17};
double n152EU_IDATEN[7] = {260, 260,  100,  80,   10, 300, 500}; // when rebinned to 5
double b152EU_IDATEN[7] = {60,  50,  50, 20,  0, 50, 200}; // when rebinned to 5

TF1* f152eu_idaten[7];
TF1* f152eu_khala[7];
TF1* f152eu_fatima[7];

const double E152EU_HPGE[11] = {244,  867, 964, 1112, 1408, 411, 778, 1299, 122, 40, 344};
const double I152EU_HPGE[11] = {7.5,  4.2, 14.5, 13.6, 21, 2.234, 12.942,  1.623, 28.58, 59.4, 26.5};
double n152EU_HPGE[11] = {260, 40, 100, 80, 80, 80, 150, 10, 500, 1000, 1000}; // when rebinned to 5
double b152EU_HPGE[11] = {10,  10, 10, 10, 10, 10, 10, 10, 10, 50, 50}; // when rebinned to 5
TF1* f152eu_hpge[11];


TH1D* hfatima[29];
TH1D* hkhala[29];
TH1D* hidaten[29];
TH1D* hhpge[29];

TH1D* h152eu;
TH1D* h152eu_hpge;
TH1D* h152eu_khala;
TH1D* h152eu_fatima;


TF1* ffatima[29];
TF1* fkhala[29];
TF1* fidaten[29];
TF1* fhpge[29];

TGraphErrors* gfatima;
TGraphErrors* gkhala;
TGraphErrors* gidaten;
TGraphErrors* ghpge;
TGraphErrors* g152eu_idaten;
TGraphErrors* g152eu_fatima;
TGraphErrors* g152eu_hpge;
TGraphErrors* g152eu_khala;

TTreeReader fReader;



const double EGATE1 = 0.122;
const double EWIND1 = 0.020; // 20 keV

const double EGATE2 = 0.344;
const double EWIND2 = 0.030; // 30 keV

double GetCountsHpge(TH1D* hhpge, double energy, double rebin){
  double integral = 0;
  for(int i = 0; i < hhpge->GetNbinsX();i++){
    if (fabs(hhpge->GetBinCenter(i)-energy) < rebin*2)
      integral+=hhpge->GetBinContent(i);
  }
  return integral;
}

void calib_singles_mono_152Eu_magnet(){
  gStyle->SetOptStat(0);

  teu = new TChain("SimulatedTree");
  teu->Add("idaten_hpge_magnet_152Eu.root");
  
  fReader.SetTree(teu);
  
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  
  TCanvas* cc = new TCanvas("cc", "cc", 1800, 600);
  cc->Divide(2,2);
  tfatima = new TChain("SimulatedTree");
  tkhala = new TChain("SimulatedTree");
  ttigress = new TChain("SimulatedTree");
  
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
  ghpge = (TGraphErrors*)gfatima->Clone("ghpge");
  ghpge->SetMarkerStyle(29);
  ghpge->SetMarkerColor(6);
  ghpge->SetLineColor(6);
  
  h152eu = new  TH1D("h152Eu", "", 1500, 0, 1500);
  h152eu_hpge = new  TH1D("h152Eu_hpge", "", 3000, 0, 1500);
  h152eu_hpge->SetLineColor(6);
  h152eu_khala = new  TH1D("h152Eu_khala", "", 1500, 0, 1500);
  h152eu_fatima = new  TH1D("h152Eu_fatima", "", 1500, 0, 1500);

  
  teu->Project("h152Eu_hpge", "Tigress.fTIG_Ge_Energy");
  
  
  for(int i = 0; i < 29; i++){
    cout<<i<<endl;
    tfatima->Reset();
    tfatima->Add(Form("idaten_hpge_magnet_%dkeV.root", ENERGIES[i]));
    tkhala->Reset();
    tkhala->Add(Form("idaten_hpge_magnet_%dkeV.root", ENERGIES[i]));
    ttigress->Reset();
    ttigress->Add(Form("idaten_hpge_magnet_%dkeV.root", ENERGIES[i]));
    
    hfatima[i] = new TH1D(Form("hfatima%d", ENERGIES[i]), "", 4500, 0, 4500);
    hkhala[i] = new TH1D(Form("hkhala%d", ENERGIES[i]), "", 4500, 0, 4500);
    hhpge[i] = new TH1D(Form("hhpge%d", ENERGIES[i]), "", 9000, 0, 4500);
    
    hfatima[i]->Rebin(REBIN[i]); 
    hkhala[i]->Rebin(REBIN[i]);
    hhpge[i]->Rebin(REBIN[i]);
    
    tfatima->Project(Form("hfatima%d", ENERGIES[i]), "Fatima.fFATIMA_LaBr3_E_Energy");
    tkhala->Project(Form("hkhala%d", ENERGIES[i]), "Khala.fKHALA_LaBr3_E_Energy");
    ttigress->Project(Form("hhpge%d", ENERGIES[i]), "Tigress.fTIG_Ge_Energy");
    
    hidaten[i] = (TH1D*)hfatima[i]->Clone(Form("hidaten%d", ENERGIES[i]));
    hidaten[i]->Add(hkhala[i]);
    ffatima[i] = new TF1(Form("ffatima%d", ENERGIES[i]), "gaus", ENERGIES[i]-dENERGIES[i], ENERGIES[i]+dENERGIES[i]);
    ffatima[i]->FixParameter(1, ENERGIES[i]);
    ffatima[i]->SetParameter(0, hfatima[i]->GetMaximum());
    ffatima[i]->SetParameter(2, dENERGIES[i]);
    fkhala[i] = (TF1*)ffatima[i]->Clone(Form("fkhala%d", ENERGIES[i]));
    fidaten[i] = (TF1*)ffatima[i]->Clone(Form("fidaten%d", ENERGIES[i]));
    // fhpge[i] = (TF1*)ffatima[i]->Clone(Form("fhpge%d", ENERGIES[i]));
    // fhpge[i]->SetParameter(2, dENERGIES[i]/10.);
    
    hfatima[i]->Fit(ffatima[i], "QN0R");
    hkhala[i]->Fit(fkhala[i], "QN0R");
    hidaten[i]->Fit(fidaten[i], "QN0R");
    //    hhpge[i]->Fit(fhpge[i], "QN0R");
    
    gfatima->SetPoint(i, ENERGIES[i], ffatima[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    gkhala->SetPoint(i, ENERGIES[i], fkhala[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    gidaten->SetPoint(i, ENERGIES[i], fidaten[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    ghpge->SetPoint(i, ENERGIES[i], GetCountsHpge(hhpge[i], ENERGIES[i], REBIN[i])/100000);

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
      cc->cd(4)->SetLogy(1);
      hhpge[i]->Draw();
      hhpge[i]->GetXaxis()->SetTitle("Energy (keV)");
      hhpge[i]->GetYaxis()->SetTitle("Counts / 5 keV");
      //      fhpge[i]->Draw("same");
      
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
      cc->cd(4);
      hhpge[i]->Draw("same");
      //      fhpge[i]->Draw("same");
    }
    
  }
		 
  TCanvas* cg = new TCanvas("cg", "cg", 1200, 800);
  cg->cd()->SetGridx(1);
  cg->cd()->SetGridy(1);
  gidaten->Draw("alp");
  gidaten->GetYaxis()->SetRangeUser(0, 0.21);
  gidaten->GetXaxis()->SetLimits(0, 4100);
  gidaten->GetXaxis()->SetTitle("Energy (keV)");
  gidaten->GetYaxis()->SetTitle("Absolute efficiency");
  gkhala->Draw("lp");
  gfatima->Draw("lp");
  ghpge->Draw("lp");
  
  gidaten->GetXaxis()->CenterTitle();
  gidaten->GetYaxis()->CenterTitle();

  
  int i =0;
  //  fReader.SetLocalEntry(entry);
  
  while(fReader.Next()){            // new addition for gamma-gamma from Rob 
    i++;
    if((i % 100000)==0)cout << i << endl;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      h152eu->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      h152eu_fatima->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      
    }
    for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      h152eu->Fill(Khala->GetKhalaLaBr3EEnergy(n));
      h152eu_khala->Fill(Khala->GetKhalaLaBr3EEnergy(n));

    }

  }
  TCanvas* c152eu = new TCanvas("c152eu", "c152eu", 1200, 800);
  c152eu->cd()->SetLogy(1);
  h152eu->SetLineColor(1);
  h152eu->Draw("hist");
  
  h152eu_fatima->SetLineColor(2);
  //  h152eu_fatima->Rebin(5);
  h152eu_fatima->Draw("hist same");
  h152eu_khala->SetLineColor(4);
  //h152eu_khala->Rebin(10);
  h152eu_khala->Draw("hist same");
  h152eu_hpge->Draw("hist same");
  
  h152eu->GetXaxis()->SetTitle("Energy (keV)");
  h152eu->GetXaxis()->CenterTitle();

  TLegend* leg152 = new TLegend(0.55, 0.6, 0.9, 0.9);
  leg152->AddEntry(h152eu, "^{152}Eu singles IDATEN (1 keV/bin)", "l");
  leg152->AddEntry(h152eu_fatima, "^{152}Eu singles FATIMA (1 keV/bin)", "l");
  leg152->AddEntry(h152eu_khala, "^{152}Eu singles KHALA (1 keV/bin)", "l");
  leg152->AddEntry(h152eu_hpge, "^{152}Eu HPGe clovers (0.5 keV/bin)", "l");
  
  leg152->Draw();
  c152eu->SaveAs("idaten_hpge_magnet_152Eu_spectrum.png");
  
  
  leg152->Draw();
  
  g152eu_idaten = new TGraphErrors();
  g152eu_fatima = new TGraphErrors();
  g152eu_khala = new TGraphErrors();
  g152eu_hpge = new TGraphErrors();

  for(int i = 0; i < 6; i++){
    f152eu_idaten[i] = new TF1(Form("f152eu_idaten%d", i), "gaus(0)+pol1(3)", E152EU_IDATEN[i]-dE152EU_IDATEN[i], E152EU_IDATEN[i]+dE152EU_IDATEN[i]);
    f152eu_idaten[i]->SetParameter(0, n152EU_IDATEN[i]);
    f152eu_idaten[i]->SetParameter(1, E152EU_IDATEN[i]);
    f152eu_idaten[i]->SetParameter(2, 10);
    f152eu_idaten[i]->SetParameter(3, b152EU_IDATEN[i]);
    
    h152eu->Fit(f152eu_idaten[i], "QRN0");
    f152eu_idaten[i]->SetLineColor(1);
    f152eu_idaten[i]->Draw("same");
    g152eu_idaten->SetPoint(i, E152EU_IDATEN[i], fabs(sqrt(2*M_PI)*f152eu_idaten[i]->GetParameter(2)*f152eu_idaten[i]->GetParameter(0)/I152EU_IDATEN[i]/1e4));
    g152eu_idaten->SetPointError(i, 0,  g152eu_idaten->GetY()[i]*sqrt(pow(f152eu_idaten[i]->GetParError(0)/f152eu_idaten[i]->GetParameter(0), 2.)+pow(f152eu_idaten[i]->GetParError(2)/f152eu_idaten[i]->GetParameter(2), 2.)));

    f152eu_fatima[i] = new TF1(Form("f152eu_fatima%d", i), "gaus(0)+pol1(3)", E152EU_IDATEN[i]-dE152EU_IDATEN[i], E152EU_IDATEN[i]+dE152EU_IDATEN[i]);
    f152eu_fatima[i]->SetParameter(0, n152EU_IDATEN[i]);
    f152eu_fatima[i]->SetParameter(1, E152EU_IDATEN[i]);
    f152eu_fatima[i]->SetParameter(2, 10);
    f152eu_fatima[i]->SetParameter(3, b152EU_IDATEN[i]);
    h152eu_fatima->Fit(f152eu_fatima[i], "QRN0");
    f152eu_fatima[i]->SetLineColor(2);
    f152eu_fatima[i]->Draw("same");
    g152eu_fatima->SetPoint(i, E152EU_IDATEN[i], fabs(sqrt(2*M_PI)*f152eu_fatima[i]->GetParameter(2)*f152eu_fatima[i]->GetParameter(0)/I152EU_IDATEN[i]/1e4));
    g152eu_fatima->SetPointError(i, 0,  g152eu_fatima->GetY()[i]*sqrt(pow(f152eu_fatima[i]->GetParError(0)/f152eu_fatima[i]->GetParameter(0), 2.)+pow(f152eu_fatima[i]->GetParError(2)/f152eu_fatima[i]->GetParameter(2), 2.)));

    f152eu_khala[i] = new TF1(Form("f152eu_khala%d", i), "gaus(0)+pol1(3)", E152EU_IDATEN[i]-dE152EU_IDATEN[i], E152EU_IDATEN[i]+dE152EU_IDATEN[i]);
    f152eu_khala[i]->SetParameter(0, n152EU_IDATEN[i]);
    f152eu_khala[i]->SetParameter(1, E152EU_IDATEN[i]);
    f152eu_khala[i]->SetParameter(2, 10);
    f152eu_khala[i]->SetParameter(3, b152EU_IDATEN[i]);
    h152eu_khala->Fit(f152eu_khala[i], "QRN0");
    f152eu_khala[i]->SetLineColor(4);
    f152eu_khala[i]->Draw("same");
    g152eu_khala->SetPoint(i, E152EU_IDATEN[i], fabs(sqrt(2*M_PI)*f152eu_khala[i]->GetParameter(2)*f152eu_khala[i]->GetParameter(0)/I152EU_IDATEN[i]/1e4));
    g152eu_khala->SetPointError(i, 0,  g152eu_khala->GetY()[i]*sqrt(pow(f152eu_khala[i]->GetParError(0)/f152eu_khala[i]->GetParameter(0), 2.)+pow(f152eu_khala[i]->GetParError(2)/f152eu_khala[i]->GetParameter(2), 2.)));

  }
  double n1290 = 0;
  double dn1290 = 0;

  for(int i = 0; i < 11; i++){
    f152eu_hpge[i] = new TF1(Form("f152eu_hpge%d", i), "gaus(0)+pol1(3)", E152EU_HPGE[i]-15, E152EU_HPGE[i]+15);
    f152eu_hpge[i]->SetLineColor(6);
    f152eu_hpge[i]->SetParameter(0, n152EU_HPGE[i]);
    f152eu_hpge[i]->SetParameter(1, E152EU_HPGE[i]);
    f152eu_hpge[i]->SetParameter(2, 2);
    f152eu_hpge[i]->SetParameter(3, b152EU_HPGE[i]);
    h152eu_hpge->Fit(f152eu_hpge[i], "QRN0");
    f152eu_hpge[i]->Draw("same");
    g152eu_hpge->SetPoint(i, E152EU_HPGE[i], fabs(sqrt(2*M_PI)*f152eu_hpge[i]->GetParameter(2)*f152eu_hpge[i]->GetParameter(0)/I152EU_HPGE[i]/1e4*2));
    g152eu_hpge->SetPointError(i, 0,  g152eu_hpge->GetY()[i]*f152eu_hpge[i]->GetParError(0)/f152eu_hpge[i]->GetParameter(0));
  }
  
  cg->cd();
  g152eu_idaten->SetMarkerStyle(25);
  g152eu_idaten->SetMarkerColor(1);
  g152eu_idaten->SetMarkerSize(1.3);
  g152eu_idaten->Draw("p");
  
  g152eu_fatima->SetMarkerStyle(25);
  g152eu_fatima->SetMarkerColor(2);
  g152eu_fatima->SetLineColor(2);
  g152eu_fatima->SetMarkerSize(1.3);
  g152eu_fatima->Draw("p");
  g152eu_khala->SetMarkerStyle(26);
  g152eu_khala->SetMarkerColor(4);
  g152eu_khala->SetLineColor(4);
  g152eu_khala->SetMarkerSize(1.3);
  g152eu_khala->Draw("p");
  
  g152eu_hpge->SetMarkerStyle(30);
  g152eu_hpge->SetMarkerColor(6);
  g152eu_hpge->SetLineColor(6);
  
  g152eu_hpge->SetMarkerSize(1.3);
  g152eu_hpge->Draw("p");
  
  TLegend* legend = new TLegend(0.55, 0.55, 0.9, 0.9);
  legend->AddEntry(gidaten, "IDATEN one #gamma (filled), ^{152}Eu (open)", "lp");
  legend->AddEntry(gfatima, "FATIMA one #gamma (filled), ^{152}Eu (open)", "lp");
  legend->AddEntry(gkhala, "KHALA one #gamma (filled), ^{152}Eu (open)", "lp");
  legend->AddEntry(ghpge, "HPGe one #gamma (filled), ^{152}Eu (open) ", "lp");
    
  legend->Draw();
  gidaten->GetXaxis()->SetLimits(0, 2000);
  cg->SaveAs("idaten_hpge_magnet_152Eu_singles_efficiency_2000keV.png");
  
  //  gidaten->GetXaxis()->SetLimits(0, 4500);
  cg->cd()->SetLogx(1);
  gidaten->GetYaxis()->SetNdivisions(11, 5, 0);
  gidaten->GetYaxis()->SetRangeUser(0, 0.21);
  cg->SaveAs("idaten_hpge_magnet_152Eu_singles_efficiency_logx.png");
  
  //  cs->SaveAs("idaten_138La_bgd_response.png");
}
