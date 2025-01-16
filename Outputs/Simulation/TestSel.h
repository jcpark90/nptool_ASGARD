//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 17 11:13:25 2020 by ROOT version 6.17/01
// from TTree SimulatedTree/Data created / analysed with the NPTool package
// found on file: myResult1.root
//////////////////////////////////////////////////////////

#ifndef TestSel_h
#define TestSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TH2.h>

// Headers needed by this particular selector
#include "TInteractionCoordinates.h"

#include "TFatimaData.h"

#include "TInitialConditions.h"



class TestSel : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   // TTreeReaderValue<TInteractionCoordinates> InteractionCoordinates = {fReader, "InteractionCoordinates"};
   TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
   TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  
   // TTreeReaderValue<TInitialConditions> InitialConditions = {fReader, "InitialConditions"};
   TTreeReaderValue<Int_t> Run = {fReader, "Run"};


   TestSel(TTree * /*tree*/ =0) { }
   virtual ~TestSel() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();


  TH1F* h0;
   TH1F* hIdatenS;
   TH1F* hKhalaS;
   TH1F* hFatimaG;
   TH1F* hKhalaG;
   TH1F* hIdatenG;
   TH2F* mIdaten;
  TH1F* ht;
  TH1F* hFK; // fatima-khala coincidence projection
   ClassDef(TestSel,0);

};

#endif

#ifdef TestSel_cxx
void TestSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t TestSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef TestSel_cxx
