//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 11 21:41:12 2015 by ROOT version 5.34/28
// from TTree T/ev1 Tree
// found on file: pytree.root
//////////////////////////////////////////////////////////

#ifndef pyAna_h
#define pyAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxentry = 3649;
   const Int_t kMaxjunction = 2;

class pyAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //Pythia8::Event  *event;
   Int_t           startColTag;
   Int_t           entry_;
   Int_t           entry_idSave[kMaxentry];   //[entry_]
   Int_t           entry_statusSave[kMaxentry];   //[entry_]
   Int_t           entry_mother1Save[kMaxentry];   //[entry_]
   Int_t           entry_mother2Save[kMaxentry];   //[entry_]
   Int_t           entry_daughter1Save[kMaxentry];   //[entry_]
   Int_t           entry_daughter2Save[kMaxentry];   //[entry_]
   Int_t           entry_colSave[kMaxentry];   //[entry_]
   Int_t           entry_acolSave[kMaxentry];   //[entry_]
   Double_t        entry_pSave_xx[kMaxentry];   //[entry_]
   Double_t        entry_pSave_yy[kMaxentry];   //[entry_]
   Double_t        entry_pSave_zz[kMaxentry];   //[entry_]
   Double_t        entry_pSave_tt[kMaxentry];   //[entry_]
   Double_t        entry_mSave[kMaxentry];   //[entry_]
   Double_t        entry_scaleSave[kMaxentry];   //[entry_]
   Double_t        entry_polSave[kMaxentry];   //[entry_]
   Bool_t          entry_hasVertexSave[kMaxentry];   //[entry_]
   Double_t        entry_vProdSave_xx[kMaxentry];   //[entry_]
   Double_t        entry_vProdSave_yy[kMaxentry];   //[entry_]
   Double_t        entry_vProdSave_zz[kMaxentry];   //[entry_]
   Double_t        entry_vProdSave_tt[kMaxentry];   //[entry_]
   Double_t        entry_tauSave[kMaxentry];   //[entry_]
   Int_t           junction_;
   Bool_t          junction_remainsSave[kMaxjunction];   //[junction_]
   Int_t           junction_kindSave[kMaxjunction];   //[junction_]
   Int_t           junction_colSave[kMaxjunction][3];   //[junction_]
   Int_t           junction_endColSave[kMaxjunction][3];   //[junction_]
   Int_t           junction_statusSave[kMaxjunction][3];   //[junction_]
   Int_t           maxColTag;
   Int_t           savedSize;
   Int_t           savedJunctionSize;
   Double_t        scaleSave;
   Double_t        scaleSecondSave;
   string          headerList;

   // List of branches
   TBranch        *b_event_startColTag;   //!
   TBranch        *b_event_entry_;   //!
   TBranch        *b_entry_idSave;   //!
   TBranch        *b_entry_statusSave;   //!
   TBranch        *b_entry_mother1Save;   //!
   TBranch        *b_entry_mother2Save;   //!
   TBranch        *b_entry_daughter1Save;   //!
   TBranch        *b_entry_daughter2Save;   //!
   TBranch        *b_entry_colSave;   //!
   TBranch        *b_entry_acolSave;   //!
   TBranch        *b_entry_pSave_xx;   //!
   TBranch        *b_entry_pSave_yy;   //!
   TBranch        *b_entry_pSave_zz;   //!
   TBranch        *b_entry_pSave_tt;   //!
   TBranch        *b_entry_mSave;   //!
   TBranch        *b_entry_scaleSave;   //!
   TBranch        *b_entry_polSave;   //!
   TBranch        *b_entry_hasVertexSave;   //!
   TBranch        *b_entry_vProdSave_xx;   //!
   TBranch        *b_entry_vProdSave_yy;   //!
   TBranch        *b_entry_vProdSave_zz;   //!
   TBranch        *b_entry_vProdSave_tt;   //!
   TBranch        *b_entry_tauSave;   //!
   TBranch        *b_event_junction_;   //!
   TBranch        *b_junction_remainsSave;   //!
   TBranch        *b_junction_kindSave;   //!
   TBranch        *b_junction_colSave;   //!
   TBranch        *b_junction_endColSave;   //!
   TBranch        *b_junction_statusSave;   //!
   TBranch        *b_event_maxColTag;   //!
   TBranch        *b_event_savedSize;   //!
   TBranch        *b_event_savedJunctionSize;   //!
   TBranch        *b_event_scaleSave;   //!
   TBranch        *b_event_scaleSecondSave;   //!
   TBranch        *b_event_headerList;   //!

   pyAna(TTree *tree=0);
   virtual ~pyAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pyAna_cxx
pyAna::pyAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pytree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pytree.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

pyAna::~pyAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pyAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pyAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pyAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("startColTag", &startColTag, &b_event_startColTag);
   fChain->SetBranchAddress("entry", &entry_, &b_event_entry_);
   fChain->SetBranchAddress("entry.idSave", entry_idSave, &b_entry_idSave);
   fChain->SetBranchAddress("entry.statusSave", entry_statusSave, &b_entry_statusSave);
   fChain->SetBranchAddress("entry.mother1Save", entry_mother1Save, &b_entry_mother1Save);
   fChain->SetBranchAddress("entry.mother2Save", entry_mother2Save, &b_entry_mother2Save);
   fChain->SetBranchAddress("entry.daughter1Save", entry_daughter1Save, &b_entry_daughter1Save);
   fChain->SetBranchAddress("entry.daughter2Save", entry_daughter2Save, &b_entry_daughter2Save);
   fChain->SetBranchAddress("entry.colSave", entry_colSave, &b_entry_colSave);
   fChain->SetBranchAddress("entry.acolSave", entry_acolSave, &b_entry_acolSave);
   fChain->SetBranchAddress("entry.pSave.xx", entry_pSave_xx, &b_entry_pSave_xx);
   fChain->SetBranchAddress("entry.pSave.yy", entry_pSave_yy, &b_entry_pSave_yy);
   fChain->SetBranchAddress("entry.pSave.zz", entry_pSave_zz, &b_entry_pSave_zz);
   fChain->SetBranchAddress("entry.pSave.tt", entry_pSave_tt, &b_entry_pSave_tt);
   fChain->SetBranchAddress("entry.mSave", entry_mSave, &b_entry_mSave);
   fChain->SetBranchAddress("entry.scaleSave", entry_scaleSave, &b_entry_scaleSave);
   fChain->SetBranchAddress("entry.polSave", entry_polSave, &b_entry_polSave);
   fChain->SetBranchAddress("entry.hasVertexSave", entry_hasVertexSave, &b_entry_hasVertexSave);
   fChain->SetBranchAddress("entry.vProdSave.xx", entry_vProdSave_xx, &b_entry_vProdSave_xx);
   fChain->SetBranchAddress("entry.vProdSave.yy", entry_vProdSave_yy, &b_entry_vProdSave_yy);
   fChain->SetBranchAddress("entry.vProdSave.zz", entry_vProdSave_zz, &b_entry_vProdSave_zz);
   fChain->SetBranchAddress("entry.vProdSave.tt", entry_vProdSave_tt, &b_entry_vProdSave_tt);
   fChain->SetBranchAddress("entry.tauSave", entry_tauSave, &b_entry_tauSave);
   fChain->SetBranchAddress("junction", &junction_, &b_event_junction_);
   fChain->SetBranchAddress("junction.remainsSave", junction_remainsSave, &b_junction_remainsSave);
   fChain->SetBranchAddress("junction.kindSave", junction_kindSave, &b_junction_kindSave);
   fChain->SetBranchAddress("junction.colSave[3]", junction_colSave, &b_junction_colSave);
   fChain->SetBranchAddress("junction.endColSave[3]", junction_endColSave, &b_junction_endColSave);
   fChain->SetBranchAddress("junction.statusSave[3]", junction_statusSave, &b_junction_statusSave);
   fChain->SetBranchAddress("maxColTag", &maxColTag, &b_event_maxColTag);
   fChain->SetBranchAddress("savedSize", &savedSize, &b_event_savedSize);
   fChain->SetBranchAddress("savedJunctionSize", &savedJunctionSize, &b_event_savedJunctionSize);
   fChain->SetBranchAddress("scaleSave", &scaleSave, &b_event_scaleSave);
   fChain->SetBranchAddress("scaleSecondSave", &scaleSecondSave, &b_event_scaleSecondSave);
   fChain->SetBranchAddress("headerList", &headerList, &b_event_headerList);
   Notify();
}

Bool_t pyAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pyAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pyAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pyAna_cxx
