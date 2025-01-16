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
TTreeReader fReader;
void GetGammaProjectionLaBr3WithPlastic(TTree* tr, TH1D* hout, double egate, double ewind){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  TTreeReaderValue<TPlasticData> Plastic = {fReader, "Plastic"};

  long long i = 0;

  while(fReader.Next()){           
    i++;
    bool coinc = 0;
    int index_fatima = -1;
    int index_khala = -1;
    
    if((i % 100000)==0) cout << i << endl;
    if (Plastic->GetMult() > 0){
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index m
	  coinc = 1;
	  index_fatima = n;
	  break;
	}
      }
      if (!coinc){
	for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	  if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index n
	    coinc = 1;
	    index_khala = n;
	    break;
	  }
	}
      }
      else{
	for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	  if (n==index_fatima) continue;
	  else
	    hout->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
	}
	for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	  if (n==index_khala) continue;
	  else
	    hout->Fill(Khala->GetKhalaLaBr3EEnergy(n));
	}
	
      }
    }
    
  }
  return hout;
}
void GetGammaProjectionLaBr3(TTree* tr, TH1D* hout, double egate, double ewind){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  long long i = 0;

  while(fReader.Next()){           
    i++;
    bool coinc = 0;
    int index_fatima = -1;
    int index_khala = -1;
    
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index m
	coinc = 1;
	index_fatima = n;
	break;
      }
    }
    if (!coinc){
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index n
	  coinc = 1;
	  index_khala = n;
	  break;
	}
      }
    }
    else{
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (n==index_fatima) continue;
	else
	  hout->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	if (n==index_khala) continue;
	else
	  hout->Fill(Khala->GetKhalaLaBr3EEnergy(n));
      }
      
    }
    
  }
  return hout;
}
void GetGammaProjectionHPGe(TTree* tr, TH1D* hout, double egate, double ewind){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  long long i = 0;

  while(fReader.Next()){           
    i++;
    bool coinc = 0;
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Tigress->GetMultiplicityGe();n++){
      //      cout<<Tigress->GetGeEnergy(n)<<endl;
      if (fabs(Tigress->GetGeEnergy(n)-egate)<ewind){
	coinc = 1;
	break;
      }
    }
   
    if(coinc){
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	hout->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	hout->Fill(Khala->GetKhalaLaBr3EEnergy(n));
      }
      
    }
    
  }
  return hout;
}

void GetStartStopHist(TTree* tr, TH1D* hout, double estart, double ewind_start, double estop, double ewind_stop){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  long long entries = tr->GetEntries();
  while(fReader.Next() && i < entries/3.7){           
    i++;
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){      
      for(int m = n+1; m < Fatima->GetFatimaLaBr3EMult(); m++){  
      	if (fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	    && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
      	  hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Fatima->GetFatimaLaBr3TTime(m));
      	else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		 && fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m 
      	  hout->Fill(Fatima->GetFatimaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
      }
    }
    for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      for(int m = n+1; m < Khala->GetKhalaLaBr3EMult(); m++){  
       	if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	    && fabs(Khala->GetKhalaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
      	  hout->Fill(Khala->GetKhalaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
      	else if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		 && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m
      	  hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Khala->GetKhalaLaBr3TTime(n));
      	
      }
    }
    //    loop through coincidence data for mixed detector types
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      for(int m = 0; m < Khala->GetKhalaLaBr3EMult(); m++){
     	if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m for khala
	    && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n for fatima
      	  hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
      	else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n for fatima
		 && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m for khala
      	  hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
      }
    }

  }
  return hout;
}
void GetStartStopHistWithGamma(TTree* tr, TH1D* hout, double estart, double ewind_start, double estop, double ewind_stop, double egam, double ewind_gam){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  while(fReader.Next()){           
    i++;
    bool gamma_gated = 0;
    //   cout<<Tigress->GetMultiplicityGe()<<endl;
    for(int n = 0; n < Tigress->GetMultiplicityGe();n++){
      //      cout<<Tigress->GetGeEnergy(n)<<endl;
      if (fabs(Tigress->GetGeEnergy(n)-egam)<ewind_gam){
	gamma_gated = 1;
	break;
      }
    }
    if (gamma_gated){
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){      
	for(int m = n+1; m < Fatima->GetFatimaLaBr3EMult(); m++){  
	  if (fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	      && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
	    hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Fatima->GetFatimaLaBr3TTime(m));
	  else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		   && fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m 
	    hout->Fill(Fatima->GetFatimaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
	}
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	for(int m = n+1; m < Khala->GetKhalaLaBr3EMult(); m++){  
	  if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	      && fabs(Khala->GetKhalaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
	    hout->Fill(Khala->GetKhalaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
	  else if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		   && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m
	    hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Khala->GetKhalaLaBr3TTime(n));
	  
	}
      }
      //    loop through coincidence data for mixed detector types
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	for(int m = 0; m < Khala->GetKhalaLaBr3EMult(); m++){
	  if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m for khala
	      && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n for fatima
	    hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
	  else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n for fatima
		   && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m for khala
	    hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
	}
      }

    }
  }
  return hout;
}
void GetStartLaBr3StopGamma(TTree* tr, TH1D* hout, double estart, double ewind_start, double estop, double ewind_stop){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  bool gamma_gated = 0;
  while(fReader.Next()){           
    i++;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult(); n++){      
      for(int m = 0; m < Tigress->GetMultiplicityGe(); m++){  
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start // start gamma in index m
	    && fabs(Tigress->GetGeEnergy(m)-estop) < ewind_stop) // stop gamma in index n
	  hout->Fill(Tigress->GetGeTimeLED(m)-Fatima->GetFatimaLaBr3TTime(n));
      }
    }
    for(int n = 0; n < Khala->GetKhalaLaBr3EMult(); n++){      
      for(int m = 0; m < Tigress->GetMultiplicityGe(); m++){  
	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estart) < ewind_start // start gamma in index m
	    && fabs(Tigress->GetGeEnergy(m)-estop) < ewind_stop) // stop gamma in index n
	  hout->Fill(Tigress->GetGeTimeLED(m)-Khala->GetKhalaLaBr3TTime(n));
      }
    }
    
  }
  return hout;
}
void GetStartStopHistWithPlastic(TTree* tr, TH1D* hout, double estop, double ewind_stop){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  TTreeReaderValue<TWAS3ABiData> Was3abi = {fReader, "WAS3ABi"};
  TTreeReaderValue<TPlasticData> Plastic = {fReader, "Plastic"};

  Long64_t i = 0;
  Long64_t entries = tr->GetEntries();
  while(fReader.Next() && i < entries/2*60/7080){ // condition for reduced beam time          
    i++;
    if((i % 100000)==0) cout << i << endl;
    double tpl_temp = 1.e12;
    //    double tpl_temp = 0;
    //    double ene_temp = 30;
    for(int n = 0; n < Plastic->GetMult();n++){
      if (Plastic->GetTime(n) < tpl_temp)
	//      if (Plastic->GetTime(n) > tpl_temp)
	tpl_temp = Plastic->GetTime(n); // get the fastest time hit on the plastic detector
    }
    if (tpl_temp < 1.e12){ // check that we have a plastic hit in this event
      //if (tpl_temp > 0){ // check that we have a plastic hit in this event
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) 
	  hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-tpl_temp);
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estop) < ewind_stop) 
	  hout->Fill(Khala->GetKhalaLaBr3TTime(n)-tpl_temp);
      }
    }    
  }
  //  return hout;
}
void GetLaBr3EnergySpectrum(TTree* tr, TH1D* hlabr3){
  TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy");
  tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy");
  hlabr3->Add(hfatima);
  hlabr3->Add(hkhala);
}
void GetLaBr3EnergySpectrumWithPlastic(TTree* tr, TH1D* hlabr3){
  // TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  // TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  // tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy", "Plastic.fPlastic_Time > 0");
  // tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy", "Plastic.fPlastic_Time > 0");
  // hlabr3->Add(hfatima);
  // hlabr3->Add(hkhala);

  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  TTreeReaderValue<TWAS3ABiData> Was3abi = {fReader, "WAS3ABi"};
  TTreeReaderValue<TPlasticData> Plastic = {fReader, "Plastic"};

  int i = 0;
  
  while(fReader.Next()){           
    i++;
    if((i % 100000)==0) cout << i << endl;
    long long tpl_temp = 1.e12;
    for(int n = 0; n < Plastic->GetMult();n++){
      if (Plastic->GetTime(n) < tpl_temp)
	tpl_temp = Plastic->GetTime(n); // get the fastest time hit on the plastic detector
    }
    if (tpl_temp < 1.e12){ // check that we have a plastic hit in this event
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3TTime(n)-tpl_temp)<10)
	hlabr3->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	if(fabs(Khala->GetKhalaLaBr3TTime(n)-tpl_temp)<10)
	hlabr3->Fill(Khala->GetKhalaLaBr3EEnergy(n)); 
      }
    }    
  }
}

void GetHPGeEnergySpectrum(TTree* tr, TH1D* hHPGe){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy");
}
void GetHPGeEnergySpectrumWithPlastic(TTree* tr, TH1D* hHPGe){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy", "Plastic.fPlastic_Time > 0");
}
void start_stop_analysis_101Sn(){
  // gStyle->SetOptStat(0);

  teu = new TChain("SimulatedTree");

  // 105Te
  //   teu->Add("105Te_decay_dr0_500ps.root");
  // // teu->Add("105Te_decay_dr3_500ps.root");
  // // teu->Add("105Te_decay_dr0_1000ps.root");
  // // teu->Add("105Te_decay_dr3_1000ps.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 200, 0, 200);
  // GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 200, 0, 200);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* ht_idaten_plastic_172 = new TH1D("ht_idaten_plastic_172", "", 120, -2, 10);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_172, 172, 15);

  
  // TF1* fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))", -2, 10);
  // fhalf->SetParameter(0, 0.5);
  // fhalf->SetParameter(1, 25);
  // fhalf->SetParameter(2, 1);
  
  // ht_idaten_plastic_172->Fit(fhalf, "RN0");
  
  // TFile* f105Te = new TFile("outputs_105Te_dr0_500ps_reduced.root", "recreate");
  // //TFile* f105Te = new TFile("outputs_105Te_dr3_500ps_reduced.root", "recreate");
  // //  TFile* f105Te = new TFile("outputs_105Te_dr0_1000ps_reduced.root", "recreate");
  // //  TFile* f105Te = new TFile("outputs_105Te_dr3_1000ps_reduced.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // ht_idaten_plastic_172->Write();
  // fhalf->Write();
  
  // f105Te->Close();
  // TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  // ht_idaten_plastic_172->Draw();
  // ht_idaten_plastic_172->GetXaxis()->SetTitle("Time (ns)");
  // ht_idaten_plastic_172->GetYaxis()->SetTitle("Counts / 0.1 ns");
  // fhalf->Draw("same");

  // TLatex* latex = new TLatex();
  // latex->DrawLatex(4, 30, Form("#tau = %.3lf #pm %.3lf ns", fhalf->GetParameter(0)/log(2), fhalf->GetParError(0)/log(2)));
  // cc->SaveAs("105Te_dr0_500ps_reduced.png");
  // cc->SaveAs("105Te_dr0_500ps_reduced.pdf");
  // cc->SaveAs("105Te_dr0_500ps_reduced.eps");
  
  // 107Te
  // teu->Add("107Te_decay_dr0_500ps.root");
  // //teu->Add("107Te_decay_dr3_500ps.root");
  // //  teu->Add("107Te_decay_dr0_1000ps.root");
  // //teu->Add("107Te_decay_dr3_1000ps.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 200, 0, 200);
  // GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 200, 0, 200);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* ht_idaten_plastic_168 = new TH1D("ht_idaten_plastic_168", "", 120, -2, 10);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_168, 168, 15);
  
  // TF1* fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*((log(2)/[0])*[2]**2-2*x))*TMath::Erfc(((log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))", -2, 10);
  // fhalf->SetParameter(0, 0.5);
  // fhalf->SetParameter(1, 15);
  // fhalf->SetParameter(2, 1);
  // ht_idaten_plastic_168->Rebin(2);
  // ht_idaten_plastic_168->GetYaxis()->SetTitle("Counts / 0.2 ns");
  
  // ht_idaten_plastic_168->Fit(fhalf, "RN0");
  
  // TFile* f107Te = new TFile("outputs_107Te_dr0_500ps.root", "recreate");
  // // TFile* f107Te = new TFile("outputs_107Te_dr3_500ps.root", "recreate");
  // //TFile* f107Te = new TFile("outputs_107Te_dr0_1000ps.root", "recreate");
  // //TFile* f107Te = new TFile("outputs_107Te_dr3_1000ps.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // ht_idaten_plastic_168->Write();
  // fhalf->Write();
  
  // f107Te->Close();
  // TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  
  // ht_idaten_plastic_168->Draw();
  // ht_idaten_plastic_168->GetXaxis()->SetTitle("Time (ns)");
  // ht_idaten_plastic_168->GetYaxis()->SetTitle("Counts / 0.2 ns");
  // fhalf->Draw("same");

  // TLatex* latex = new TLatex();
  // latex->DrawLatex(4, ht_idaten_plastic_168->GetMaximum()*0.5, Form("#tau = %.3lf #pm %.3lf ns", fhalf->GetParameter(0)/log(2), fhalf->GetParError(0)/log(2)));
  // cc->SaveAs("107Te_dr0_500ps_reduced.png");
  // cc->SaveAs("107Te_dr0_500ps_reduced.pdf");
  // cc->SaveAs("107Te_dr0_500ps_reduced.eps");

  // // 102Sn
  // teu->Add("102Sn_isomer_2cm.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* ht_idaten_88_497 = new TH1D("ht_idaten_88_497", "", 400, -2, 2);  
  // GetStartStopHist(teu, ht_idaten_88_497, 88, 10, 497, 20);
  
  // TFile* f102Sn = new TFile("outputs_102Sn_2cm.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // // he_idaten_labr3_gate->Write();
  // ht_idaten_88_497->Write();
  // f102Sn->Close();

  // TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  // ht_idaten_88_497->Draw();
  // ht_idaten_88_497->Rebin(5);
  // ht_idaten_88_497->GetXaxis()->SetTitle("Time (ns)");
  // ht_idaten_88_497->GetYaxis()->SetTitle("Counts / 0.05 ns");
  // TF1* flifetime = new TF1("flifetime", "gaus", -2, 2);
  // ht_idaten_88_497->Fit(flifetime, "RN0");
  // flifetime->Draw("same");

  // TLatex* latex = new TLatex();
  // latex->DrawLatex(-1.8, ht_idaten_88_497->GetMaximum()*0.5, Form("#tau = %.3lf #pm %.3lf ns", flifetime->GetParameter(1), flifetime->GetParError(1)));
  // cc->SaveAs("102Sn_100ps_reduced.png");
  // cc->SaveAs("102Sn_100ps_reduced.pdf");
  // cc->SaveAs("102Sn_100ps_reduced.eps");
  
  //  104Sn
  teu->Add("104Sn_isomer.root");
  TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  GetLaBr3EnergySpectrum(teu, he_idaten);
  TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  GetHPGeEnergySpectrum(teu, he_hpge);  
  TH1D* ht_idaten_314_683 = new TH1D("ht_idaten_314_683", "", 400, -2, 2);  
  GetStartStopHist(teu, ht_idaten_314_683, 314, 15, 683, 25);
  
  TFile* f104Sn = new TFile("outputs_104Sn.root", "recreate");
  he_idaten->Write();
  he_hpge->Write();
  // he_idaten_labr3_gate->Write();
  ht_idaten_314_683->Write();
  f104Sn->Close();
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  ht_idaten_314_683->Draw();
  ht_idaten_314_683->Rebin(5);
  ht_idaten_314_683->GetXaxis()->SetTitle("Time (ns)");
  ht_idaten_314_683->GetYaxis()->SetTitle("Counts / 0.05 ns");
  TF1* flifetime = new TF1("flifetime", "gaus", -2, 2);
  ht_idaten_314_683->Fit(flifetime, "RN0");
  flifetime->Draw("same");
   TLatex* latex = new TLatex();
  latex->DrawLatex(-1.8, ht_idaten_314_683->GetMaximum()*0.5, Form("#tau = %.3lf #pm %.3lf ns", flifetime->GetParameter(1), flifetime->GetParError(1)));
  cc->SaveAs("104Sn_29ps_reduced.png");
  cc->SaveAs("104Sn_29ps_reduced.pdf");
  cc->SaveAs("104Sn_29ps_reduced.eps");
}
