#define pyAna_cxx
#include "pyAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void pyAna::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L pyAna.C
//      Root > pyAna t
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
	if (fChain == 0) return;
	//Set-up code
	Long64_t nentries = fChain->GetEntriesFast();
	int pions=0,kaons=0,ksTOpi=0,i=0,j=0;
	int ipart1,ipart2,a,b;
	Double_t p[kMaxentry];
	Double_t ppion[kMaxentry];
	Double_t npion[kMaxentry];
	Double_t pt[kMaxentry];
	double prop=0;
	Double_t pp,E,KSmass;
	//~ TLorentzVector q[kMaxentry];//pinakas tupou tlorentz
	
	
	

	TH1D *pzhist=new TH1D("pz","pz histogram",120,-12,12);
	TH1D *pxhist=new TH1D("px","px histogram-errors",120,-3,3);
	TH1D *pthist=new TH1D("pt","pt histogram",50,0,1.5);
	TH1D *kshist=new TH1D("KSmass","KSmasss",120,0.2,1);
	TH1D *kshist2=new TH1D("KSmass","KSmasss",120,0.,5);
	
	
	
	cout << "Number of events = " << nentries << endl;
	

	Long64_t nbytes = 0, nb = 0;
	
	
	for (Long64_t jentry=0; jentry<nentries;jentry++) {//<nentries
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		cout<<endl;
		cout << " event # " << jentry << ",  #particles " <<  entry_ << endl;
		
		// loop over all particles in this "event"
		
		
		//Loop code
		for (int ipart = 0; ipart < entry_ ; ipart++) {
			cout<<"ipart="<<ipart<<endl;
			//~ TLorentzVector w(entry_pSave_xx[ipart],entry_pSave_yy[ipart],entry_pSave_zz[ipart],entry_pSave_tt[ipart]);
			//~ q[ipart]=w;
			pt[ipart]=TMath::Sqrt((entry_pSave_xx[ipart]*entry_pSave_xx[ipart]) + (entry_pSave_yy[ipart]*entry_pSave_yy[ipart]));
			p[ipart]=TMath::Sqrt((entry_pSave_xx[ipart]*entry_pSave_xx[ipart]) + (entry_pSave_yy[ipart]*entry_pSave_yy[ipart])+(entry_pSave_zz[ipart]*entry_pSave_zz[ipart]));
			//~ cout<<"p="<<p[ipart]<<"\t px="<<entry_pSave_xx[ipart]<<"\t py="<<entry_pSave_yy[ipart]<<"\t pz="<<entry_pSave_zz[ipart]<<endl;
			//~ cout <<"i="<<ipart<< "\tpdgID = " << entry_idSave[ipart]<<",\t Daughter1="<<entry_idSave[entry_daughter1Save[ipart]]<<",\tDaughter2="<<entry_idSave[entry_daughter2Save[ipart]]<<"\tEnergy="<<entry_pSave_tt[ipart]<<endl;
			
			
			//~ cout<<"edw to xw"<<q[ipart].Px()<<endl;
			
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
				
				ipart1=entry_daughter1Save[ipart];
		        ipart2=entry_daughter2Save[ipart];
		        pp=TMath::Sqrt(((entry_pSave_xx[ipart1]+entry_pSave_xx[ipart2])*(entry_pSave_xx[ipart1]+entry_pSave_xx[ipart2]))+
				((entry_pSave_yy[ipart1]+entry_pSave_yy[ipart2])*(entry_pSave_yy[ipart1]+entry_pSave_yy[ipart2]))+
				((entry_pSave_zz[ipart1]+entry_pSave_zz[ipart2])*(entry_pSave_zz[ipart1]+entry_pSave_zz[ipart2])));
				E=entry_pSave_tt[ipart1]+entry_pSave_tt[ipart2];
				KSmass=TMath::Sqrt((E*E)-((pp*pp)));
				kshist->Fill(KSmass);
				
				//~ cout<<"mother-ipart="<<ipart<<"\t Emoth="<<entry_pSave_tt[ipart]<<"\t Pmoth="<<p[ipart]<<endl;
				//~ cout<<"ptotal="<<pp<<endl;				
				//~ cout<<"E1="<<entry_pSave_tt[ipart1]<<"\tE2="<<entry_pSave_tt[ipart2]<<"\tEtotal="<<E[i]<<endl;
				//~ cout<<"p1x="<<entry_pSave_xx[ipart1]<<"\tpy1="<<entry_pSave_yy[ipart1]<<"\tpz1="<<entry_pSave_zz[ipart1]<<endl;
				//~ cout<<"p2x="<<entry_pSave_xx[ipart2]<<"\tpy2="<<entry_pSave_yy[ipart2]<<"\tpz2="<<entry_pSave_zz[ipart2]<<endl;
				//~ cout<<"ksmass="<<KSmass<<endl;	
				ksTOpi++;}//end of ks->pions if
				
			if(entry_idSave[ipart]==211){
				ppion[i]=ipart;
				//~ cout<<"ppion found,ppion["<<i<<"]="<<ppion[i]<<endl;
				i++;}
				
			if(entry_idSave[ipart]==(-1)*(211)){
				npion[j]=ipart;
				j++;
				//~ cout<<"npion found,ipart="<<ipart<<endl;
				}
			
			
			//Fill Histograms
			pzhist->Fill(entry_pSave_zz[ipart]);
			pxhist->Fill(entry_pSave_xx[ipart]);
			pthist->Fill(pt[ipart]);
			
			
			
		} // for each particle in this envent
		
		//Ksmass 
			for(int n=0;n<i;n++){
				a=ppion[n];
				for(int m=0;m<j;m++){
					b=npion[m];
					
					pp=TMath::Sqrt(((entry_pSave_xx[a]+entry_pSave_xx[b])*(entry_pSave_xx[a]+entry_pSave_xx[b]))+
					((entry_pSave_yy[a]+entry_pSave_yy[b])*(entry_pSave_yy[a]+entry_pSave_yy[b]))+
					((entry_pSave_zz[a]+entry_pSave_zz[b])*(entry_pSave_zz[a]+entry_pSave_zz[b])));
					E=entry_pSave_tt[a]+entry_pSave_tt[b];
					KSmass=TMath::Sqrt((E*E)-((pp*pp)));
					kshist2->Fill(KSmass);
					
					
				}//end of npion loop
			}//end of ppion loop
			
			i=0;
			j=0;
		
	if(kaons!=0){
		prop=prop+((double)kaons)/((double)pions);
		}
		
   }// for each event
   
   
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

//~ k1->cd(2);
//~ pxhist->Draw();
//~ pxhist->GetXaxis()->SetTitle("mass [GeV/c^2]");
//~ pxhist->GetYaxis()->SetTitle("Number of events");



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
cout<<"Ksto pi="<<ksTOpi<<endl;
cout<<q.Energy()<<endl;



}// end of Loop() method
