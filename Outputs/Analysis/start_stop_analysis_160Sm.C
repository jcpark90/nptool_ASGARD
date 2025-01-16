TChain* tfatima;
TChain* tkhala;
TChain* ttigress;


TChain* teu;
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TFatimaData.h"
#include "TTigressData.h"

TH1D* hdt_idaten;
double EGATE_START;
double EWIND_START;
double EGATE_STOP;
double EWIND_STOP;


TTreeReader fReader;

const int GATE_SET = 2; // 2 settings for examining the half-life of interest in 160Sm
const double TIME_CUTOFF = 10000; // 10 us, beyond which the spectrum may be mostly 138La background

void start_stop_analysis_160Sm(){
  gStyle->SetOptStat(0);

  if (GATE_SET==1){
    // set 1 for 7.5-ns 6+ state in 160Sm
    EGATE_START = 106; // 106 keV to 6+ state in 160Sm
    EWIND_START = 10; // 10 keV
    EGATE_STOP = 325; // 325 keV from 6+ state in 160Sm
    EWIND_STOP = 20; // 20 keV
  }
  else if (GATE_SET==2){
    // set 2 for 1.0-ns 4+ state in 160Sm
     EGATE_START = 162; // 106 keV to 6+ state in 160Sm
     EWIND_START = 10; // 20 keV
     EGATE_STOP = 71; // 325 keV from 6+ state in 160Sm
     EWIND_STOP = 8; // 30 keV
  }
  teu = new TChain("SimulatedTree");
  teu->Add("160Sm_isomer.root");
  fReader.SetTree(teu);
  
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};

  if (GATE_SET==1)
    hdt_idaten = new TH1D("hdt_idaten", "", 150, -10, 50); // start-stop time difference histogram for 7.5-ns state
  else if (GATE_SET==2)
    hdt_idaten = new TH1D("hdt_idaten", "", 170, -2, 15); // start-stop time difference histogram for 1.0-ns state
  
  
  int i = 0;
  
  while(fReader.Next()){           
    i++;
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      if (Fatima->GetFatimaLaBr3TTime(n) < TIME_CUTOFF){
	for(int m = n+1; m < Fatima->GetFatimaLaBr3EMult(); m++){
	  if (Fatima->GetFatimaLaBr3TTime(m) < TIME_CUTOFF){
	    if (fabs(Fatima->GetFatimaLaBr3EEnergy(m)-EGATE_START) < EWIND_START // start gamma in index m
		&& fabs(Fatima->GetFatimaLaBr3EEnergy(n)-EGATE_STOP) < EWIND_STOP) // stop gamma in index n
	      hdt_idaten->Fill(Fatima->GetFatimaLaBr3TTime(n)-Fatima->GetFatimaLaBr3TTime(m));
	    else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-EGATE_START) < EWIND_START   // start gamma in index n
		     && fabs(Fatima->GetFatimaLaBr3EEnergy(m)-EGATE_STOP) < EWIND_STOP) // stop gamma in index m 
	      hdt_idaten->Fill(Fatima->GetFatimaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
	  }
	}
      }
    }
    for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      if (Khala->GetKhalaLaBr3TTime(n) < TIME_CUTOFF){
	for(int m = n+1; m < Khala->GetKhalaLaBr3EMult(); m++){
	  if (Khala->GetKhalaLaBr3TTime(m) < TIME_CUTOFF){
	    if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-EGATE_START) < EWIND_START // start gamma in index m
		&& fabs(Khala->GetKhalaLaBr3EEnergy(n)-EGATE_STOP) < EWIND_STOP) // stop gamma in index n
	      hdt_idaten->Fill(Khala->GetKhalaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
	    else if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-EGATE_START) < EWIND_START   // start gamma in index n
		     && fabs(Khala->GetKhalaLaBr3EEnergy(m)-EGATE_STOP) < EWIND_STOP) // stop gamma in index m
	      hdt_idaten->Fill(Khala->GetKhalaLaBr3TTime(m)-Khala->GetKhalaLaBr3TTime(n));
	  }
	}
      }
    }
    //    loop through coincidence data for mixed detector types
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      if (Fatima->GetFatimaLaBr3TTime(n) < TIME_CUTOFF){
	for(int m = 0; m < Khala->GetKhalaLaBr3EMult(); m++){
	  if (Khala->GetKhalaLaBr3TTime(m) < TIME_CUTOFF){
	    if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-EGATE_START) < EWIND_START // start gamma in index m for khala
		&& fabs(Fatima->GetFatimaLaBr3EEnergy(n)-EGATE_STOP) < EWIND_STOP) // stop gamma in index n for fatima
	      hdt_idaten->Fill(Fatima->GetFatimaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
	    else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-EGATE_START) < EWIND_START   // start gamma in index n for fatima
		     && fabs(Khala->GetKhalaLaBr3EEnergy(m)-EGATE_STOP) < EWIND_STOP) // stop gamma in index m for khala
	      hdt_idaten->Fill(Khala->GetKhalaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
	  }
	}
      }
    }
    
  }
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  
  hdt_idaten->SetLineColor(1);
  hdt_idaten->GetXaxis()->SetTitle("T_{stop}-T_{start} (ns)");
  if (GATE_SET==1)
    hdt_idaten->GetYaxis()->SetTitle("Counts / 0.4 ns");
  else if (GATE_SET==2)
    hdt_idaten->GetYaxis()->SetTitle("Counts / 0.05 ns");
  
  hdt_idaten->Draw("hist");
  TF1* fhalf;

  // simple exponential fit with background
  fhalf = new TF1("fhalf", "[1]*0.5**(x/[0])+[2]", 1, 15);
  fhalf->SetParameter(0, 25);
  fhalf->SetParameter(1, 3.0);
  fhalf->SetParameter(2, 1);

  // incorporate time jitter into the T1/2 fit, which results in an exponentially modified gaussian
  
  // if (GATE_SET==1){
  //   fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))", 0, 50);
  //   fhalf->SetParameter(0, 7.5);
  //   fhalf->SetParameter(1, 25);
  //   fhalf->SetParameter(2, 1);
  // }
  // if (GATE_SET==2){
  //   fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))", -2, 15);
  //   fhalf->SetParameter(0, 1.0);
  //   fhalf->SetParameter(1, 3);
  //   fhalf->SetParameter(2, 1);
  //   fhalf->SetParameter(3, 0.5);
    
  // }
  
  
  hdt_idaten->Fit(fhalf, "QRN0");
  fhalf->Draw("same");

   TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  if (GATE_SET==1){
    latex->DrawLatex(20, 30, "^{160}Sm: 106(start)-325(stop)");
    latex->DrawLatex(20, 25, Form("T_{1/2}(6^{+}) = %.2f #pm %.2f ns (Sim.)", fhalf->GetParameter(0),  fhalf->GetParError(0)));
    latex->DrawLatex(20, 20, Form("T_{1/2}(6^{+}) = %.1f ns (Geant4 ENSDF)", 7.5));
    cc->SaveAs("start_stop_analysis_160Sm_4plus.png");
  }
  else if (GATE_SET==2){
    latex->DrawLatex(4, 80, "^{160}Sm: 162(start)-71(stop)");
    latex->DrawLatex(4, 65, Form("T_{1/2}(2^{+}) = %.2f #pm %.2f ns (Sim.)", fhalf->GetParameter(0),  fhalf->GetParError(0)));
    latex->DrawLatex(4, 50, Form("T_{1/2}(2^{+}) = %.1f ns (Geant4, Lit.)", 3.0));
    cc->SaveAs("start_stop_analysis_160Sm_2plus.png");
  }
}
