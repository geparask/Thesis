#define testlor_cxx
#include "testlor.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void testlor::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L testlor.C
//      Root > testlor t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

//Set-up Code
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	int pions=0,kaons=0,ksTOpi=0,i=0,j=0;
	int ipart1,ipart2,a,b;
	Double_t ppion[kMaxentry];
	Double_t npion[kMaxentry];
	Double_t p[kMaxentry];
	double prop=0;
	Double_t pp,E,KSmass;
   	TLorentzVector q[kMaxentry];//pinakas tupou tlorentz
	TLorentzVector w,v; 
	TVector3 l;
	
	TH1D *pzhist=new TH1D("pz","pz histogram",120,-12,12);
	TH1D *pthist=new TH1D("pt","pt histogram",50,0,1.5);
	TH1D *kshist=new TH1D("KSmass","KSmasss",120,0.2,0.8);
	TH1D *kshist2=new TH1D("KSmass","KSmasss",120,0.,5);
	TH1D *pseudohist=new TH1D("Pseudorapidity","Pseudorapidity",120,-15,15);
	
	cout << "Number of events = " << nentries << endl;
	
   //event loop
	for (Long64_t jentry=0; jentry<nentries;jentry++) {//nentries
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		cout << " event # " << jentry << ",  #particles " <<  entry_ << endl;
		// if (Cut(ientry) < 0) continue;

		//particle loop
		for (int ipart = 0; ipart < entry_ ; ipart++) {
			w.SetPx(entry_pSave_xx[ipart]);
			w.SetPy(entry_pSave_yy[ipart]);
			w.SetPz(entry_pSave_zz[ipart]);
			w.SetE(entry_pSave_tt[ipart]);
			q[ipart]=w;
			l=q[ipart].Vect();
			//~ cout<<"pseudorapidity="<<l.PseudoRapidity()<<endl;
			p[ipart]=l.Mag();
			cout<<"p="<<l.Mag()<<endl;//P!
			//~ cout<<"m="<<q[ipart].M()<<endl;//m?not right
			//~ cout<<"pt="<<q[ipart].Perp()<<endl;//PT!
			//~ cout <<"i="<<ipart<< "\tpdgID = " << entry_idSave[ipart]<<",\t Daughter1="<<entry_idSave[entry_daughter1Save[ipart]]<<",\tDaughter2="<<entry_idSave[entry_daughter2Save[ipart]]<<"\tEnergy="<<entry_pSave_tt[ipart]<<endl;
		
			//pions vs kaons
			if(((TMath::Abs(entry_idSave[ipart]))==(211))||(entry_idSave[ipart]==(111))) {
				pions++;
			} 	
			else if((entry_idSave[ipart]==130)||
				((TMath::Abs(entry_idSave[ipart]))==(321))||(entry_idSave[ipart]==310)
				||((TMath::Abs(entry_idSave[ipart]))==(311))||((TMath::Abs(entry_idSave[ipart]))==(323))
				||((TMath::Abs(entry_idSave[ipart]))==(313)) ) {
				kaons++; 
			}
			
			
			//Ks->pions
			if((entry_idSave[ipart]==310)&&(((entry_idSave[entry_daughter1Save[ipart]]==(211))&&(entry_idSave[entry_daughter2Save[ipart]]==(-1)*(211)))
			||((entry_idSave[entry_daughter1Save[ipart]]==(-1)*(211))&&(entry_idSave[entry_daughter2Save[ipart]]==(211))) )){
				
				ipart1[ksTOpi]=entry_daughter1Save[ipart];
		        ipart2[ksTOpi]=entry_daughter2Save[ipart];
		        cout<<"ipart1="<<ipart1<<"\t ipart2="<<ipart2<<endl;
				ksTOpi++;
				cout<<"ksTOpi="<<ksTOpi<<endl;
				
			}//end of ks->pions if
				
			if(entry_idSave[ipart]==211){
				ppion[i]=ipart;
				i++;}
				
			if(entry_idSave[ipart]==(-1)*(211)){
				npion[j]=ipart;
				j++;}
			
			
			//Fill Histograms
			pzhist->Fill(entry_pSave_zz[ipart]);
			pthist->Fill(q[ipart].Perp());
			pseudohist->Fill(l.PseudoRapidity()); 		
		
		}//end of particle loop
      
		
		//Ksmass 
		for(int n=0;n<i;n++){
			a=ppion[n];
			for(int m=0;m<j;m++){
				b=npion[m];
				v=q[a]+q[b];
				//~ TVector3 y=v.Vect();
				//~ pp=y.Mag();
				//~ E=v.E();
				KSmass=v.M();
				kshist2->Fill(KSmass);
					
			}//end of npion loop
		}//end of ppion loop
				
		//KStoPI
		for(int count=0;count<ksTOpi;count++){
			v=q[ipart1[ksTOpi]]+q[ipart2[ksTOpi]];
			//~ l=v.Vect();
			//~ pp=l.Mag();
			//~ E=v.E();
			KSmass=v.M();
			kshist->Fill(KSmass);
		}//end of ksTOpi loop
			
		i=0;
		j=0;
		ksTOpi=0;
		
	if(kaons!=0){
		prop=prop+((double)kaons)/((double)pions);
		}
		
	
	}//end of event loop
   
	//Wrap-up code
	TCanvas* k1 = new TCanvas("c1","Pythia8Ana",800,800);
	k1->Divide(1, 3);

	k1->cd(1);
	kshist->Draw();
	kshist->GetXaxis()->SetTitle("mass [GeV/c^2]");
	kshist->GetYaxis()->SetTitle("Number of events");

	k1->cd(2);
	kshist2->Draw();
	kshist2->GetXaxis()->SetTitle("mass [GeV/c^2]");
	kshist2->GetYaxis()->SetTitle("Number of events");



	//~ 
	//~ k1->cd(2);
	//~ pseudohist->Draw();
	//~ pseudohist->GetXaxis()->SetTitle("pseudorapidity");
	//~ pseudohist->GetYaxis()->SetTitle("Number of events");



	k1->cd(3);
	pthist->Draw();
	pthist->GetXaxis()->SetTitle("pt [GeV/c]");
	pthist->GetYaxis()->SetTitle("Number of events");

	cout<<"Ebeam1="<<entry_pSave_tt[1]<<"[GeV], Ebeam2="<<entry_pSave_tt[2]<<"[GeV]"<<endl;
	cout<<"Mass1= "<<entry_mSave[1]<<"  Mass2="<<entry_mSave[2]<<endl;
	cout<<"Kaons="<<prop/nentries<<"Pions"<<endl;
	cout<<"kaons="<<kaons<<", Pions="<<pions<<endl;
	cout<<"event found, ipart1="<<ipart1<<"\t ipart2="<<ipart2<<endl;
	cout<<"ipart1="<<ipart1<<endl;
	cout<<"p["<<ipart1<<"]="<<p[ipart1]<<endl;

}//End of Loop method
 
