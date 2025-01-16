TChain* tfatima;
TChain* tkhala;
TChain* ttigress;


TChain* teu;
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TFatimaData.h"
#include "TTigressData.h"
#include "TWAS3ABiData.h"
#include "TPlasticData.h"


TH1D* hdt_idaten;
double EGATE_IDATEN;
double EWIND_IDATEN;


TTreeReader fReader;

void start_stop_analysis_plastic_idaten_ES(){


  gStyle->SetOptStat(0);

  // set 1 for 10 ps 9/2- state in 77Cu
  EGATE_IDATEN = 946; // 946 keV to 5/2- state in 77Cu
  EWIND_IDATEN = 40; // 40 keV
  teu = new TChain("SimulatedTree");
  teu->Add("idaten77cu_30ps_81ps.root");
  fReader.SetTree(teu);


  int nbin = 200;
  double bmin = -2000;
  double bmax = 2000;
  double bin = (bmax - bmin)/nbin;

  hdt_idaten = new TH1D("hdt_idaten", "", nbin, bmin, bmax); // start-stop time difference histogram

  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  TTreeReaderValue<TWAS3ABiData> Was3abi = {fReader, "WAS3ABi"};
  TTreeReaderValue<TPlasticData> Plastic = {fReader, "Plastic"};



  int i = 0;

  while(fReader.Next()){
    i++;
    if((i % 100000)==0) cout << i << endl;
    double tpl_temp = 1.e13;
    for(int n = 0; n < Plastic->GetMult();n++){
      if (Plastic->GetTime(n) < tpl_temp)
	tpl_temp = Plastic->GetTime(n); // get the fastest time hit on the plastic detector
    }
    if (tpl_temp < 1.e13){ // check that we have a plastic hit in this event
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-EGATE_IDATEN) < EWIND_IDATEN)
	  hdt_idaten->Fill((Fatima->GetFatimaLaBr3TTime(n)-tpl_temp)*1000);
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-EGATE_IDATEN) < EWIND_IDATEN)
	  hdt_idaten->Fill((Khala->GetKhalaLaBr3TTime(n)-tpl_temp)*1000);
      }
    }
  }
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
//cc->SetLogy();
  hdt_idaten->SetLineColor(1);
  hdt_idaten->GetXaxis()->SetTitle("T_{labr3}-T_{plastic} (ps)");
hdt_idaten->GetYaxis()->SetTitle( Form("Counts / %.0f ps", bin) );

  hdt_idaten->Draw("hist");
  hdt_idaten->GetYaxis()->SetRangeUser(0.5, 1050);
  // simple exponential fit with background
  //  TF1* fhalf = new TF1("fhalf", "[0]*0.5**(x/[1])+[2]", 0, 50);
  // fhalf->SetParameter(0, 25);
  // fhalf->SetParameter(1, 7.5);
  // fhalf->SetParameter(2, 1);

  // incorporate time jitter into the T1/2 fit, which results in an exponentially modified gaussian

  TF1* fhalf;
  fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*(2*[3]+(log(2)/[0])*[2]**2-2*x))*TMath::Erfc(([3]+(log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))", -920, 1100);
  fhalf->SetParameter(0, 30);
  fhalf->SetParameter(1, 200);
  fhalf->SetParameter(2, 1);
  fhalf->SetParameter(3, 0);

  hdt_idaten->Fit(fhalf, "QRN0");
  fhalf->Draw("same");

  TF1* fhalf_clone = (TF1*)fhalf->Clone("fhalf_clone");
  fhalf_clone->SetLineColor(1);
  fhalf_clone->FixParameter(3, 0); // assume no shift
  hdt_idaten->Fit(fhalf_clone, "QRN0");				    
  fhalf_clone->Draw("same");
  

  cout<<fhalf->GetParameter(3)<<endl;

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.03);
  latex->DrawLatex(-1500, 900, "^{77}Cu from ^{77}Ni #beta- decay: beta(start), 946-keV gamma(stop)");
  latex->SetTextColor(2);
  latex->DrawLatex(-1500, 800, Form("T_{1/2}(9/2^{-}) = %.3f #pm %.4f ps (Sim. with zero offset, %.3f ps)", fhalf->GetParameter(0),  fhalf->GetParError(0), fhalf->GetParameter(3)));
  latex->SetTextColor(1);
  latex->DrawLatex(-1500, 700, Form("T_{1/2}(9/2^{-}) = %.3f #pm %.4f ps (Sim. without offset)", fhalf_clone->GetParameter(0),  fhalf_clone->GetParError(0)));
  latex->DrawLatex(-1500, 600, Form("T_{1/2}(9/2^{-}) = %.2f ps (Geant4 ENSDF)", 30.));
  cc->SaveAs("start_stop_analysis_77Cu_9halfminus.png");
}
