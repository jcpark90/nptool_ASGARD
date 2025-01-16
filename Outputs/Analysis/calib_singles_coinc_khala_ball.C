TChain* tkhala;
TChain* ttigress;


TChain* teu;
TChain* tco;
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

//#include "TKhalaData.h"
#include "TTigressData.h"

const int ENERGIES[25] = {30, 40, 50, 60, 70, 80, 90,
			  100, 150, 200, 250, 300, 400, 500,
			  600, 700, 1000,
			  1100, 1200, 1300, 
			  1600, 1700, 1800, 1900, 2000};

const int dENERGIES[25] = {10, 10, 10, 15, 15, 15, 15,
			   20, 20, 20, 25, 25, 30, 35,
			   40, 45, 50,
			   50, 50, 50, 
			   50, 50, 50, 50, 50};

const int REBIN[25] = {1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 2, 2, 2, 2,
		       2, 2, 5,
		       5, 5, 5,
		       5, 5, 5, 5, 5};

const double E152EU_122[5] = {244,  867, 964, 1112, 1408};
const double I152EU_122[5] = {7.5,  0.0161+4.258+0.163+0.0669+0.084, 14.5, 13.6, 21};
double n152EU_122[5] = {260, 40, 100, 80, 80}; // when rebinned to 5
double b152EU_122[5] = {60,  40, 100, 80, 80}; // when rebinned to 5
TF1* f152eu_122[5];

const double E152EU_344[2] = {411, 778};
const double dE152EU_344[2] = {30, 40};

const double I152EU_344[2] = {2.234, 12.99+0.0265+0.191};

const double E152EU_HPGE[10] = {244,  867, 964, 1112, 1408, 411, 778, 1299, 122, 40};
const double I152EU_HPGE[10] = {7.5,  0.0161+4.258+0.163+0.0669+0.084, 14.5, 13.6+1.735+0.186, 21, 2.234, 12.942,  1.623, 28.58, 59.4};
double n152EU_HPGE[10] = {260, 40, 100, 80, 80, 80, 150, 10, 500, 1000}; // when rebinned to 5
double b152EU_HPGE[10] = {10,  10, 10, 10, 10, 10, 10, 10, 10, 50}; // when rebinned to 5

double n152EU_344[2] = {80, 150}; // when rebinned to 5
double b152EU_344[2] = {20, 5}; // when rebinned to 5
TF1* f152eu_344[2];

TF1* f152eu_hpge[8];

TH1D* hkhala[25];
TH1D* hhpge[25];

TH1D* h152eu;
TH1D* h60co;
TH1D* h152eu_hpge;

TH1D* h152eu_gg_122;
TH1D* h152eu_gg_344;
TH1D* h60co_gg_1173; // gate
TH1D* h60co_gg_1332; // gate


TF1* fkhala[25];
TF1* fhpge[25];

TGraphErrors* gkhala;
TGraphErrors* ghpge;
TGraphErrors* g152eu_122;
TGraphErrors* g152eu_hpge;
TGraphErrors* g152eu_344;
TGraphErrors* g60co_gg;

TTreeReader fReader152eu;
TTreeReader fReader60co;



const double EGATE1 = 122;
const double EWIND1 = 20; // 20 keV

const double EGATE2 = 344;
const double EWIND2 = 30; // 30 keV

const double EGATE3 = 1173;
const double EWIND3 = 40; // 20 keV

const double EGATE4 = 1332;
const double EWIND4 = 45; // 30 keV

double GetCountsHpge(TH1D* hhpge, double energy, double rebin){
  double integral = 0;
  for(int i = 0; i < hhpge->GetNbinsX();i++){
    if (fabs(hhpge->GetBinCenter(i)-energy) < rebin*2)
      integral+=hhpge->GetBinContent(i);
  }
  return integral;
}

void calib_singles_coinc_khala_ball(){
  gStyle->SetOptStat(0);

  teu = new TChain("SimulatedTree");
  teu->Add("../Simulation/khala_ball_MS2024_152Eu_2mm.root");
  fReader152eu.SetTree(teu);
  
  tco = new TChain("SimulatedTree");
  tco->Add("../Simulation/khala_ball_MS2024_60Co_2mm.root");
  fReader60co.SetTree(tco);
  
  TTreeReaderValue<TKhalaData> Khala_eu = {fReader152eu, "Khala"};
  TTreeReaderValue<TKhalaData> Khala_co = {fReader60co, "Khala"};
  //  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  
  TCanvas* cc = new TCanvas("cc", "cc", 1800, 600);
  cc->Divide(2,2);
  tkhala = new TChain("SimulatedTree");
  //  ttigress = new TChain("SimulatedTree");
  
  gkhala = new TGraphErrors();
  gkhala->SetMarkerSize(1.3);
  gkhala->SetMarkerStyle(20);
  gkhala->SetMarkerColor(2);
  gkhala->SetLineColor(2);
  // ghpge = (TGraphErrors*)gfatima->Clone("ghpge");
  // ghpge->SetMarkerStyle(23);
  // ghpge->SetMarkerColor(6);
  // ghpge->SetLineColor(6);
  
  h152eu = new  TH1D("h152Eu", "", 1500, 0, 1500);
  h152eu_hpge = new  TH1D("h152Eu_hpge", "", 3000, 0, 1500);
  h152eu_hpge->SetLineColor(6);
  
  h152eu_gg_122 = new  TH1D("h152Eu_gg_122", "", 1500, 0, 1500);
  h152eu_gg_344 = new  TH1D("h152Eu_gg_344", "", 1500, 0, 1500);

  h60co = new  TH1D("h60Co", "", 1500, 0, 1500);
  h60co_gg_1173 = new  TH1D("h60co_gg_1173", "", 1500, 0, 1500);
  h60co_gg_1332 = new  TH1D("h60co_gg_1332", "", 1500, 0, 1500);

  //  teu->Project("h152Eu_hpge", "Tigress.fTIG_Ge_Energy");
  
  for(int i = 0; i < 25; i++){
    cout<<i<<endl;
    tkhala->Reset();
    tkhala->Add(Form("../Simulation/khala_ball_ms2024_%dkeV.root", ENERGIES[i]));
    // ttigress->Reset();
    // ttigress->Add(Form("idaten_hpge_%dkeV.root", ENERGIES[i]));
    
    hkhala[i] = new TH1D(Form("hkhala%d", ENERGIES[i]), "", 4500, 0, 4500);
    //    hhpge[i] = new TH1D(Form("hhpge%d", ENERGIES[i]), "", 9000, 0, 4500);
    
    hkhala[i]->Rebin(REBIN[i]);
    //    hhpge[i]->Rebin(REBIN[i]);
    
    tkhala->Project(Form("hkhala%d", ENERGIES[i]), "Khala.fKHALA_LaBr3_E_Energy");
    //    ttigress->Project(Form("hhpge%d", ENERGIES[i]), "Tigress.fTIG_Ge_Energy");
    
    fkhala[i] = new TF1(Form("fkhala%d", ENERGIES[i]), "gaus", ENERGIES[i]-dENERGIES[i]*2/3, ENERGIES[i]+dENERGIES[i]);
    fkhala[i]->FixParameter(1, ENERGIES[i]);
    fkhala[i]->SetParameter(0, hkhala[i]->GetMaximum());
    fkhala[i]->SetParameter(2, dENERGIES[i]);
    // fhpge[i] = (TF1*)ffatima[i]->Clone(Form("fhpge%d", ENERGIES[i]));
    // fhpge[i]->SetParameter(2, dENERGIES[i]/10.);
    
    hkhala[i]->Fit(fkhala[i], "QN0R");
    //    hhpge[i]->Fit(fhpge[i], "QN0R");
    
    gkhala->SetPoint(i, ENERGIES[i], fkhala[i]->Integral(ENERGIES[i]*0.5, ENERGIES[i]*1.5)/REBIN[i]/100000);
    //    ghpge->SetPoint(i, ENERGIES[i], GetCountsHpge(hhpge[i], ENERGIES[i], REBIN[i])/100000);

    if (i==0){
      // cc->cd(1)->SetLogy(1);
      cc->cd(2)->SetLogy(1);
      hkhala[i]->Draw();
      hkhala[i]->GetXaxis()->SetTitle("Energy (keV)");
      hkhala[i]->GetYaxis()->SetTitle("Counts / 5 keV");
      fkhala[i]->Draw("same");
      //cc->cd(4)->SetLogy(1);
      
      // hhpge[i]->Draw();
      // hhpge[i]->GetXaxis()->SetTitle("Energy (keV)");
      // hhpge[i]->GetYaxis()->SetTitle("Counts / 5 keV");
      //      fhpge[i]->Draw("same");
      
    }
    else{
      cc->cd(2);
      hkhala[i]->Draw("same");
      fkhala[i]->Draw("same");
      // cc->cd(4);
      // hhpge[i]->Draw("same");
      //      fhpge[i]->Draw("same");
    }
    
  }
		 
  TCanvas* cg = new TCanvas("cg", "cg", 1200, 800);
  cg->cd()->SetGridx(1);
  cg->cd()->SetGridy(1);
  gkhala->Draw("alp");
  gkhala->GetYaxis()->SetRangeUser(0, 0.21);
  gkhala->GetXaxis()->SetLimits(0, 4100);
  gkhala->GetXaxis()->SetTitle("Energy (keV)");
  gkhala->GetYaxis()->SetTitle("Absolute efficiency");
  gkhala->Draw("lp");
  //  ghpge->Draw("lp");
  
  gkhala->GetXaxis()->CenterTitle();
  gkhala->GetYaxis()->CenterTitle();

  
  int i =0;
  //  fReader.SetLocalEntry(entry);
  
  while(fReader152eu.Next()){            // new addition for gamma-gamma from Rob 
    i++;
    if((i % 100000)==0)cout << i << endl;
    for(int n = 0; n < Khala_eu->GetKhalaLaBr3EMult();n++){
      h152eu->Fill(Khala_eu->GetKhalaLaBr3EEnergy(n));
      for(int m = n+1; m < Khala_eu->GetKhalaLaBr3EMult(); m++){  
      	if (Khala_eu->GetKhalaLaBr3EEnergy(m)>EGATE1-EWIND1 && Khala_eu->GetKhalaLaBr3EEnergy(m)<EGATE1+EWIND1)
      	  h152eu_gg_122->Fill(Khala_eu->GetKhalaLaBr3EEnergy(n));
      	else if (Khala_eu->GetKhalaLaBr3EEnergy(m)>EGATE2-EWIND2 && Khala_eu->GetKhalaLaBr3EEnergy(m)<EGATE2+EWIND2)
      	  h152eu_gg_344->Fill(Khala_eu->GetKhalaLaBr3EEnergy(n));
	if (Khala_eu->GetKhalaLaBr3EEnergy(n)>EGATE1-EWIND1 && Khala_eu->GetKhalaLaBr3EEnergy(n)<EGATE1+EWIND1)
      	  h152eu_gg_122->Fill(Khala_eu->GetKhalaLaBr3EEnergy(m));
      	else if (Khala_eu->GetKhalaLaBr3EEnergy(n)>EGATE2-EWIND2 && Khala_eu->GetKhalaLaBr3EEnergy(n)<EGATE2+EWIND2)
      	  h152eu_gg_344->Fill(Khala_eu->GetKhalaLaBr3EEnergy(m));      	
      }
    }
  }
  i = 0;
  while(fReader60co.Next()){            // new addition for gamma-gamma from Rob 
    i++;
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Khala_co->GetKhalaLaBr3EMult();n++){
      h60co->Fill(Khala_co->GetKhalaLaBr3EEnergy(n));
      for(int m = 0; m < Khala_co->GetKhalaLaBr3EMult(); m++){
	if (n==m) continue;
	if (Khala_co->GetKhalaLaBr3EEnergy(n)>EGATE3-EWIND3 && Khala_co->GetKhalaLaBr3EEnergy(n)<EGATE3+EWIND3)
      	  h60co_gg_1173->Fill(Khala_co->GetKhalaLaBr3EEnergy(m));
      	else if (Khala_co->GetKhalaLaBr3EEnergy(n)>EGATE4-EWIND4 && Khala_co->GetKhalaLaBr3EEnergy(n)<EGATE4+EWIND4)
      	  h60co_gg_1332->Fill(Khala_co->GetKhalaLaBr3EEnergy(m));
      }
    }
  }
  
  TCanvas* c152eu = new TCanvas("c152eu", "c152eu", 1200, 800);
  c152eu->cd()->SetLogy(1);
  h152eu->SetLineColor(1);
  h152eu->Draw("hist");
  h152eu_gg_122->SetLineColor(2);
  h152eu_gg_122->Rebin(5);
  h152eu_gg_122->Draw("hist same");
  h152eu_gg_344->SetLineColor(4);
  h152eu_gg_344->Rebin(10);
  h152eu_gg_344->Draw("hist same");

  TCanvas* c60co = new TCanvas("c60co", "c60co", 1200, 800);
  c60co->cd()->SetLogy(1);
  h60co->SetLineColor(1);
  h60co->Draw("hist");
    
  h60co_gg_1173->SetLineColor(2);
  h60co_gg_1173->Rebin(5);
  h60co_gg_1173->Draw("hist same");
  h60co_gg_1332->SetLineColor(4);
  h60co_gg_1332->Rebin(5);
  h60co_gg_1332->Draw("hist same");
  
  //  h152eu_hpge->Draw("hist same");
  TF1* f122_singles = new TF1("f122_singles", "gaus(0)+pol1(3)", 122-15, 122+12);
  f122_singles->SetLineColor(1);
  f122_singles->SetParameter(0, 3000);
  f122_singles->FixParameter(1, 122);
  f122_singles->SetParameter(2, 10);
  f122_singles->SetParameter(3, 100);  
  h152eu->Fit(f122_singles, "QRN0");
  f122_singles->Draw("same");
  TF1* f344_singles = new TF1("f344_singles", "gaus(0)+pol1(3)", 344-20, 344+20);
  f344_singles->SetLineColor(1);
  f344_singles->SetParameter(0, 3000);
  f344_singles->FixParameter(1, 344);
  f344_singles->SetParameter(2, 10);
  f344_singles->SetParameter(3, 100);  
  h152eu->Fit(f344_singles, "QRN0");
  f344_singles->Draw("same");
  
  double n122_singles = sqrt(2*M_PI)*f122_singles->GetParameter(2)*f122_singles->GetParameter(0);
  double e122_singles = n122_singles/0.284/1e6;
  
  double n344_singles = sqrt(2*M_PI)*f344_singles->GetParameter(2)*f344_singles->GetParameter(0);
  double e344_singles = n344_singles/0.265/1e6;

  TF1* f1173_singles = new TF1("f1173_singles", "gaus(0)+pol1(3)", EGATE3-EWIND3, EGATE3+EWIND3);
  f1173_singles->SetLineColor(1);
  f1173_singles->SetParameter(0, 3000);
  f1173_singles->FixParameter(1, EGATE3);
  f1173_singles->SetParameter(2, EWIND3);
  f1173_singles->SetParLimits(2, 0, 500);
  
  f1173_singles->SetParameter(3, 100);  
  h60co->Fit(f1173_singles, "QRN0");
  f1173_singles->SetLineColor(9);
  f1173_singles->Draw("same");
  TF1* f1332_singles = new TF1("f1332_singles", "gaus(0)+pol1(3)", EGATE4-EWIND4, EGATE4+EWIND4);
  f1332_singles->SetLineColor(1);
  f1332_singles->SetParameter(0, 3000);
  f1332_singles->FixParameter(1, EGATE4);
  f1332_singles->SetParLimits(2, 0, 500);
  f1332_singles->SetParameter(2, 10);
  f1332_singles->SetParameter(3, 100);  
  h60co->Fit(f1332_singles, "QRN0");
  f1332_singles->SetLineColor(9);
  f1332_singles->Draw("same");
  
  double n1173_singles = sqrt(2*M_PI)*f1173_singles->GetParameter(2)*f1173_singles->GetParameter(0);
  double e1173_singles = n1173_singles/1e6;
  
  double n1332_singles = sqrt(2*M_PI)*f1332_singles->GetParameter(2)*f1332_singles->GetParameter(0);
  double e1332_singles = n1332_singles/1e6;

  
  TLegend* leg152 = new TLegend(0.55, 0.6, 0.9, 0.9);
  leg152->AddEntry(h152eu, "^{152}Eu singles KHALA (1 keV/bin)", "l");
  leg152->AddEntry(h152eu_gg_122, "Gated on 122-keV #gamma (5 keV/bin)", "l");
  leg152->AddEntry(h152eu_gg_344, "Gated on 344-keV #gamma (10 keV/bin)", "l");
  //  leg152->AddEntry(h152eu_hpge, "^{152}Eu HPGe clovers (0.5 keV/bin)", "l");
  
  leg152->Draw();
  c152eu->SaveAs("khala_ball_152Eu_spectrum.png");
  
  g152eu_122 = new TGraphErrors();
  g152eu_344 = new TGraphErrors();
  g60co_gg = new TGraphErrors();
  
  g152eu_hpge = new TGraphErrors();
  
  for(int i = 0; i < 5; i++){
    f152eu_122[i] = new TF1(Form("f152eu_122%d", i), "gaus(0)+pol1(3)", E152EU_122[i]-40, E152EU_122[i]+40);
    f152eu_122[i]->SetParameter(0, n152EU_122[i]);
    f152eu_122[i]->FixParameter(1, E152EU_122[i]);
    f152eu_122[i]->SetParameter(2, 10);
    f152eu_122[i]->SetParameter(3, b152EU_122[i]);
    
    
    h152eu_gg_122->Fit(f152eu_122[i], "QRN0");
    f152eu_122[i]->Draw("same");
    g152eu_122->SetPoint(i, E152EU_122[i], sqrt(2*M_PI)*f152eu_122[i]->GetParameter(2)*f152eu_122[i]->GetParameter(0)/5/n122_singles/(I152EU_122[i]/(21+14.5+13.6+7.5)));
    g152eu_122->SetPointError(i, 0,  g152eu_122->GetY()[i]*f152eu_122[i]->GetParError(0)/f152eu_122[i]->GetParameter(0));
    
  }
  double n1290 = 0;
  double dn1290 = 0;
  for(int i = 0; i < 2; i++){
    f152eu_344[i] = new TF1(Form("f152eu_344%d", i), "gaus(0)+pol1(3)", E152EU_344[i]-dE152EU_344[i], E152EU_344[i]+dE152EU_344[i]);
    f152eu_344[i]->SetParameter(0, n152EU_344[i]);
    f152eu_344[i]->FixParameter(1, E152EU_344[i]);
    f152eu_344[i]->SetParameter(2, 10);
    f152eu_344[i]->FixParameter(3, b152EU_344[i]);
    if (i > 1)
      f152eu_344[i]->FixParameter(4, 0);
       
    h152eu_gg_344->Fit(f152eu_344[i], "QRN0");
    f152eu_344[i]->SetLineColor(4);
    f152eu_344[i]->Draw("same");
    if (i==2){
      for(int j = 0; j < h152eu_gg_344->GetNbinsX(); j++){
	if (h152eu_gg_344->GetBinCenter(j) > 1250 && h152eu_gg_344->GetBinCenter(j) < 1350)
	  n1290+=h152eu_gg_344->GetBinContent(j);
      }
      dn1290 = sqrt(n1290);
      g152eu_344->SetPoint(i, E152EU_344[i], n1290/n344_singles/(I152EU_344[i]/(2.234+12.942+1.727+1.623)*0.045/0.06));
      g152eu_344->SetPointError(i, 0,  g152eu_344->GetY()[i]*dn1290/n1290);
    }
    else{
      g152eu_344->SetPoint(i, E152EU_344[i], sqrt(2*M_PI)*f152eu_344[i]->GetParameter(2)*f152eu_344[i]->GetParameter(0)/10/n344_singles/(I152EU_344[i]/(2.234+12.942+1.727+1.623)*0.045/0.06));
      g152eu_344->SetPointError(i, 0,  g152eu_344->GetY()[i]*f152eu_344[i]->GetParError(0)/f152eu_344[i]->GetParameter(0));
    }
      
  }
  TF1* f60co_1173 = new TF1("f60co_1173", "gaus(0)+pol1(3)", EGATE4-EWIND4, EGATE4+EWIND4);
  f60co_1173->SetParameter(0, 5000);
  f60co_1173->FixParameter(1, EGATE4);
  f60co_1173->SetParameter(2, 10);
  f60co_1173->SetParameter(3, 30);
  h60co_gg_1173->Fit(f60co_1173, "QRN0");
  f60co_1173->SetLineColor(6);
  f60co_1173->Draw("same");
  
  g60co_gg->SetPoint(0, EGATE4, sqrt(2*M_PI)*f60co_1173->GetParameter(2)*f60co_1173->GetParameter(0)/5/n1173_singles);
  g60co_gg->SetPointError(0, 0,  g60co_gg->GetY()[0]*sqrt(pow(f60co_1173->GetParError(0)/f60co_1173->GetParameter(0), 2.)+pow(f60co_1173->GetParError(2)/f60co_1173->GetParameter(2), 2.)));
  cout<<"eff 1332: "<<sqrt(2*M_PI)*f60co_1173->GetParameter(2)*f60co_1173->GetParameter(0)/5/n1173_singles<<endl;
  TF1* f60co_1332 = new TF1("f60co_1332", "gaus(0)+pol1(3)", EGATE3-EWIND3, EGATE3+EWIND3);
  f60co_1332->SetParameter(0, 5000);
  f60co_1332->FixParameter(1, EGATE3);
  f60co_1332->SetParameter(2, 10);
  f60co_1332->SetParameter(3, 30);
  h60co_gg_1332->Fit(f60co_1332, "QRN0");
  f60co_1332->SetLineColor(kOrange);
  f60co_1332->Draw("same");
  
  g60co_gg->SetPoint(1, EGATE3, sqrt(2*M_PI)*f60co_1332->GetParameter(2)*f60co_1332->GetParameter(0)/5/n1332_singles);
  g60co_gg->SetPointError(1, 0,  g60co_gg->GetY()[1]*sqrt(pow(f60co_1332->GetParError(0)/f60co_1332->GetParameter(0), 2.)+pow(f60co_1332->GetParError(2)/f60co_1332->GetParameter(2), 2.)));
  cout<<"eff 1173: "<<sqrt(2*M_PI)*f60co_1332->GetParameter(2)*f60co_1332->GetParameter(0)/5/n1332_singles<<endl;

  
  // for(int i = 0; i < 10; i++){
  //   f152eu_hpge[i] = new TF1(Form("f152eu_hpge%d", i), "gaus(0)+pol1(3)", E152EU_HPGE[i]-15, E152EU_HPGE[i]+15);
  //   f152eu_hpge[i]->SetLineColor(6);
  //   f152eu_hpge[i]->SetParameter(0, n152EU_HPGE[i]);
  //   f152eu_hpge[i]->SetParameter(1, E152EU_HPGE[i]);
  //   f152eu_hpge[i]->SetParameter(2, 2);
  //   f152eu_hpge[i]->SetParameter(3, b152EU_HPGE[i]);
  //   h152eu_hpge->Fit(f152eu_hpge[i], "QRN0");
  //   f152eu_hpge[i]->Draw("same");
  //   g152eu_hpge->SetPoint(i, E152EU_HPGE[i], fabs(sqrt(2*M_PI)*f152eu_hpge[i]->GetParameter(2)*f152eu_hpge[i]->GetParameter(0)/I152EU_HPGE[i]/1e4*2));
  //   g152eu_hpge->SetPointError(i, 0,  g152eu_hpge->GetY()[i]*f152eu_hpge[i]->GetParError(0)/f152eu_hpge[i]->GetParameter(0));
  // }
  
  cg->cd();
  g152eu_122->SetMarkerStyle(28);
  g152eu_122->SetMarkerColor(1);
  g152eu_122->SetMarkerSize(1.3);
  g152eu_122->Draw("p");
  g152eu_344->SetMarkerStyle(27);
  g152eu_344->SetMarkerColor(1);
  g152eu_344->SetMarkerSize(1.3);
  g152eu_344->Draw("p");
  g60co_gg->SetMarkerStyle(26);
  g60co_gg->SetMarkerColor(1);
  g60co_gg->SetMarkerSize(1.3);
  g60co_gg->Draw("p");
  
  // g152eu_hpge->SetMarkerStyle(26);
  // g152eu_hpge->SetMarkerColor(6);
  // g152eu_hpge->SetLineColor(6);
  
  // g152eu_hpge->SetMarkerSize(1.3);
  // g152eu_hpge->Draw("p");
  
  TLegend* legend = new TLegend(0.6, 0.55, 0.9, 0.9);
  legend->AddEntry(gkhala, "KHALA (46 ch) mono", "lp");
  //  legend->AddEntry(ghpge, "HPGe (2 clovers) mono", "lp");
  
  legend->AddEntry(g152eu_122, "^{152}Eu #gamma#gamma (122-keV gate)", "lp");
  legend->AddEntry(g152eu_344, "^{152}Eu #gamma#gamma (344-keV gate)", "lp");
  legend->AddEntry(g60co_gg, "^{60}Co #gamma#gamma", "lp");

  //  legend->AddEntry(g152eu_hpge, "^{152}Eu singles #gamma HPGe", "lp");
  
  legend->Draw();
  cg->SaveAs("khala_ball_152Eu_efficiency_4000keV.png");
  gkhala->GetXaxis()->SetLimits(0, 2000);
  cg->SaveAs("khala_ball_152Eu_efficiency_2000keV.png");
  
  //  gkhala->GetXaxis()->SetLimits(0, 4500);
  //  cg->cd()->SetLogx(1);
  cg->cd();
  gkhala->GetYaxis()->SetNdivisions(11, 5, 0);
  gkhala->GetYaxis()->SetRangeUser(0, 0.35);
  cg->SaveAs("khala_ball_152Eu_efficiency_logx.png");
  
}
