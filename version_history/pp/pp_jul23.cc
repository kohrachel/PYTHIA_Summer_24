// This file contains code modified from PYTHIA examples main142.cc and main143.cc

// This program simulates pp collisions at 5TeV with the option to add pileup.
// The goal of this simulation is to study the effects of background and
// subtraction methods in the 2 and 3-point correlators.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1F.h"
#include "TH2F.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TStyle.h"

// ROOT, for math
#include "TMath.h"
#include "TRandom3.h"

// FastJet
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//==========================================================================

//! Hard-coded settings
//! Note: Almost all customization settings can/should be done here

//! Jet and hadron pT thresholds.
//! Will only show particles with pT > pTmin and |y| < yMax.
double pTmin_jet = 30;
double pTmin_hadron = 0.2;
double yMax = 3;
int power = -1; //anti-kT algorithm

//! Histogram jet energy range
int pT_lowbound = 120;
int pT_highbound = 140;

//! Eta cut for jet.eta()
double eta_cut = 1.9;

//! Number of events to generate
int num_events = 1000;

//! Number of jets
int num_jets_debug = 0;

//! Radius of jet
double radius = 0.5;

// Amount of pileup. Average number of inelastic pp collisions per event
// (=bunch-crossing). Set to zero to turn off pileup.
double mu = 200;
// int n_inel = gRandom->Poisson(mu);

//! Switch for activating mixed cones (for substraction)
bool mixed_cones = true;

//! Resolution cut
double res_cut = 0.01;

//! EPSILON value for r_L cuts in shape dep graph
double EPSILON = 0.005;

//! boolean for controlling whether to divide by jetpT in weightage
bool divJet = false;

//! boolean for signal/background property
int origin = -2;

using namespace Pythia8;

class MyUserInfo : public fastjet::PseudoJet::UserInfoBase{
public:
  // default ctor
  //  - pdg_id        the PDG id of the particle
  //  - origin     whether it is signal (TRUE) or background (FALSE)
  MyUserInfo(const int & pdg_id_in, const bool & origin_in) :
    _pdg_id(pdg_id_in), _origin(origin_in){}
  
  /// access to the PDG id
  int pdg_id() const { return _pdg_id;}
   
  /// access to the origin (0 = signal; -1 = bg in signal cone; 1,2,3 = mixed)
  int origin() const { return _origin;}
   
protected:
  int _pdg_id;      // the associated pdg id
  bool _origin;     // the associated property
};

//==========================================================================

int main(int argc, char* argv[]) {

  //! Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  //! Create Pythia instance and set it up to generatem hard QCD processes
  //! above pTHat = 20 GeV for pp collisions at 5 TeV.
  Pythia pythia;
  pythia.readString("HardQCD:all = on"); 
  pythia.readString("PhaseSpace:pTHatMin = 100.");
  pythia.readString("PhaseSpace:pTHatMax = -1.");
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");

  //! Pileup switches
  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");

  //! Hadronization switch
  pythia.readString("HadronLevel:Hadronize = on");

  //! Set random seed
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  
  //! Turn on all pion decay
  pythia.readString("421:mayDecay = on");
  pythia.readString("411:mayDecay = on");
  pythia.readString("111:mayDecay = on");
  pythia.readString("211:mayDecay = on");

  //! If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  //! Pileup particles
  Pythia pythiaPU;
  pythiaPU.readString("Beams:eCM = 5020.");
  pythiaPU.readString("SoftQCD:inelastic = on");
  // if (n_inel > 0) pythiaPU.init();
  if (mu > 0) pythiaPU.init();

  //! Setup fastjet. 
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, radius, power);

  //! Create file on which histogram(s) can be saved. IMPORTANT UPDATE PU AMOUNT
  //! For SLURM:
  TFile* outFile = new TFile(Form("pp_%ito%iGeV_R05_N%dK_jet_charged_PU%i_xdivJetPT_breakdown_%s.root",
				  pT_lowbound,pT_highbound,num_events/1000,200,argv[1]), 
			          "RECREATE");
  //! For interactive:
  // TFile* outFile = new TFile(Form("pp_%ito%iGeV_R05_N%dK_jet_charged_PU%i_xdivJetPT_poisson.root",
  // 				  pT_lowbound,pT_highbound,num_events/1000,100), 
  // 			          "RECREATE"); // NOTE THIS IS FOR NON SLURM


  //! EEC histogram setup ======================================================
  int bins = 100;
  int low = 0;
  int high = 1;
  TString x_axis = "log10(#Deltar)";
  TString y_axis = "1/N_{pairs}";

  //! Title
  TString t1 = "EEC_pt1 (";
  TString t2 = " GeV jets), both track p_{T} > ";
  TString t3 = " GeV";  
  TString title = t1 + pT_lowbound + "-" + pT_highbound + t2;

  //! Energy of both particles > 1 GeV 
  TH1F *EEC_pt1 = new TH1F("EEC_pt1", title + 1 + t3 + " all comb", bins, low, high);
  TH1F *EEC_pt1_SS = new TH1F("EEC_pt1_SS", title + 1 + t3 + " S+S", bins, low, high);
  TH1F *EEC_pt1_SB = new TH1F("EEC_pt1_SB", title + 1 + t3 + " S+B", bins, low, high);
  TH1F *EEC_pt1_BB = new TH1F("EEC_pt1_BB", title + 1 + t3 + " B+B", bins, low, high);

  // TH1F *EEC_pt3 = new TH1F("EEC_pt3", title + 3 + t3 + " all comb", bins, low, high);
  // TH1F *EEC_pt3_SS = new TH1F("EEC_pt3_SS", title + 3 + t3 + " S+S", bins, low, high);
  // TH1F *EEC_pt3_SB = new TH1F("EEC_pt3_SB", title + 3 + t3 + " S+B", bins, low, high);
  // TH1F *EEC_pt3_BB = new TH1F("EEC_pt3_BB", title + 3 + t3 + " B+B", bins, low, high);

  TH1F *EECs_pt1[] = {EEC_pt1, EEC_pt1_SS, EEC_pt1_SB, EEC_pt1_BB};
  // TH1F *EECs_pt3[] = {EEC_pt3, EEC_pt3_SS, EEC_pt3_SB, EEC_pt3_BB};
  
  // //! Update title and book histogram, cut on 1 particle only
  // t2 = " GeV jets), 1 track p_{T} >";
  // title = t1 + pT_lowbound + "-" + pT_highbound + t2;
  // TH1F *EEC_pt1_halfCut = new TH1F("EEC_pt1_halfCut", title + 1 + t3, bins, -3, 0);

  //! Set axis labels
  for (auto hist:EECs_pt1) {
    hist->GetXaxis()->SetTitle(x_axis);
    hist->GetYaxis()->SetTitle(y_axis);
  }

  //! E3C histogram setup ======================================================
  t1 = "E3C (";
  t2 = " GeV jets), track p_{T} >";
  t3 = " GeV";  
  title = t1 + pT_lowbound + "-" + pT_highbound + t2;
  
  //! Energy of all particles must be above the threshold
  TH1F *E3C_gen = new TH1F("E3C_gen", title + 1 + t3 + " all comb", bins, -3, 0);
  TH1F *E3C_sm1 = new TH1F("E3C_sm1", title + 1 + t3 + " SSB+SBB+BBB", bins, -3, 0);
  TH1F *E3C_sm12 = new TH1F("E3C_sm12", title + 1 + t3 + " SBB+BBB", bins, -3, 0);
  TH1F *E3C_m123 = new TH1F("E3C_m123", title + 1 + t3 + " BBB", bins, -3, 0);

  TH1F *E3Cs[] = {E3C_gen, E3C_sm1, E3C_sm12, E3C_m123};

  //! Set axis labels
  for (auto hist:E3Cs) {
    hist->GetXaxis()->SetTitle(x_axis);
    hist->GetYaxis()->SetTitle(y_axis);
  }

  //! Dummy histograms =========================================================

  //! Number of particles per event
  TH1F *num_ptcls_event = new TH1F("num_ptcls_event", "Number of Particles/Event",
			     1000, 0.0, 20000.0);

  //! Number of particles per jet
  TH1F *num_ptcls_jet = new TH1F("num_ptcls_jet", "Number of Particles/Jet",
			     150, 0.0, 150.0);
  
  //! Dummy histogram for per jet normalization
  TH1I *num_jets = new TH1I("num_jets", "Number of Jets", 1, 0, 1);

  //! Event ====================================================================

  //! Begin event loop. Generate event; skip if generation aborted.
  auto &event = pythia.event;
  for (int iEvent = 0; iEvent < num_events; ++iEvent) {
    if (!pythia.next()) continue;

    //! Identify particles. Jets are built from all stable particles after hadronization (particle-level jets).
    //! ch: charged (S+B); sig: signal, S; pu: pileup, B
    std::vector<fastjet::PseudoJet> ch_particles, m1_ptcls, m2_ptcls, m3_ptcls;
    int iptcls_ptgt1 = 0;
    
    //! Loop over all particles in the event 
    for (int i = 0; i < event.size(); ++i) {
      auto &p = event[i];
      if (p.isFinal()){

	//! label this as a signal particle
	origin = 0;
	auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	cur_ptcl.set_user_info(new MyUserInfo(p.id(), origin));

	if (p.isCharged()) {
	  ch_particles.push_back(cur_ptcl);
	}//! charged particles
      }//! final particles
    }//! all event particles loop

    //! Add in pileup
    int n_inel = gRandom->Poisson(mu);
    // printf("Overlaying particles from %i pileup interactions!\n", n_inel);
    for (int i_pu= 0; i_pu<n_inel; ++i_pu) {
      if (!pythiaPU.next()) continue;
      for (int i = 0; i < pythiaPU.event.size(); ++i) {
        auto &p = pythiaPU.event[i];
        if (p.isFinal()) {

	  //! Label this as not a signal particle
	  origin = -1;
	  auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	  cur_ptcl.set_user_info(new MyUserInfo(p.id(), origin));
	  
	  if (p.isCharged()) {
	    ch_particles.push_back(cur_ptcl);
	  }//! charged pileup particles
	}//! final pileup particles
      }//! pileup event loop
    }//! main pileup loop
    
    //! Debugging pileup:
    cout << "All ch ptcls: " << ch_particles.size() << endl;
    num_ptcls_event->Fill(ch_particles.size());

    //! Making all cuts on jets here
    vector <fastjet::PseudoJet> inclusiveJets, selectedJets;
    fastjet::ClusterSequence clustSeq(ch_particles, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(pTmin_jet); //only cluster jets above pTmin_jet
    fastjet::Selector select_pt_min = fastjet::SelectorPtMin(pT_lowbound); //set min jet pT
    fastjet::Selector select_pt_max = fastjet::SelectorPtMax(pT_highbound);//set max jet pT
    fastjet::Selector select_rap = fastjet::SelectorAbsRapMax(eta_cut);
    fastjet::Selector select_all = select_pt_min && select_rap && select_pt_max;//order of operation doesn't matter in how its defined here
    selectedJets = select_all(inclusiveJets);

    //! If there are no jets, skip this event
    if(selectedJets.size() == 0) {
      continue;
    }

    //! Update the number of jets
    num_jets_debug += selectedJets.size();
    
    //! For each jet:
    for(int i=0; i < int(selectedJets.size()); i++) {

      fastjet::PseudoJet jet = selectedJets[i];
      fastjet::PseudoJet dummy_loc(jet.eta(), jet.phi(), jet.pt(), jet.e());
      
      vector <fastjet::PseudoJet> constit = jet.constituents();

      //! Debugging pileup
      num_ptcls_jet->Fill(constit.size());
      num_jets->Fill(0);

      //! Mixed cones
      if (mixed_cones) {
	for (int im1 = 0; im1 < n_inel; im1++) {
	  if (!pythiaPU.next()) continue;
	  for (int i = 0; i < pythiaPU.event.size(); i++) {
	    auto &p = pythiaPU.event[i];

	    origin = 1;
	    auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	    cur_ptcl.set_user_info(new MyUserInfo(p.id(), origin));

	    if (p.isFinal()) {
	      if (p.isCharged()) {
		if (cur_ptcl.delta_R(jet) < radius) {
		  m1_ptcls.push_back(cur_ptcl);
		}//!delta_R
	      }//!charged m1
	    }//!final m1
	  }//!m1 event loop
	}//!m1 loop

	for (int im2 = 0; im2 < n_inel; im2++) {
	  if (!pythiaPU.next()) continue;
	  for (int i = 0; i < pythiaPU.event.size(); i++) {
	    auto &p = pythiaPU.event[i];

	    origin = 2;
	    auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	    cur_ptcl.set_user_info(new MyUserInfo(p.id(), origin));

	    if (p.isFinal()) {
	      if (p.isCharged()) {
		if (cur_ptcl.delta_R(jet) < radius) {
		  m2_ptcls.push_back(cur_ptcl);
		}//!delta_R
	      }//!charged m2
	    }//!final m2
	  }//!m2 event loop
	}//!m2 loop

	for (int im3 = 0; im3 < n_inel; im3++) {
	  if (!pythiaPU.next()) continue;
	  for (int i = 0; i < pythiaPU.event.size(); i++) {
	    auto &p = pythiaPU.event[i];

	    origin = 3;
	    auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	    cur_ptcl.set_user_info(new MyUserInfo(p.id(), origin));

	    if (p.isFinal()) {
	      if (p.isCharged()) {
		if (cur_ptcl.delta_R(jet) < radius) {
		  m3_ptcls.push_back(cur_ptcl);
		}//!delta_R
	      }//!charged m3
	    }//!final m3
	  }//!m3 event loop
	}//!m3 loop
      }//!mixed cones end

      //! Creating the signal cone =============================================

      double r_L;
      //! For each particle within the radius:
      for (int x = 0; x < int(constit.size()); x++) {
	auto jet_1 = constit[x];

	//! Jet + M1
	for (int im1 = 0; im1 < int(m1_ptcls.size()) - 1; im1++) {
	  auto bg_1 = m1_ptcls[im1];
	  auto bg_2 = m1_ptcls[im1+1];

	  if (jet_1.pt() > 1 && bg_1.pt() > 1 && bg_2.pt() > 1) {
	    r_L = TMath::Max(TMath::Max(jet_1.delta_R(bg_1),
					jet_1.delta_R(bg_2)),
			                bg_1.delta_R(bg_2));

	    if (r_L < res_cut) continue;

	    E3Cs[1]->Fill(r_L, jet_1.E() * bg_1.E() * bg_2.E());
	  }//!pt cuts
	}//!jet + m1

	//! Jet + M1 + M2
	for (int im1 = 0; im1 < int(m1_ptcls.size()); im1++) {
	  auto bg_1 = m1_ptcls[im1];

	  for (int im2 = 0; im2 < int(m2_ptcls.size()); im2++) {
	    auto bg_2 = m2_ptcls[im2];

	    if (jet_1.pt() > 1 && bg_1.pt() > 1 && bg_2.pt() > 1) {
	      r_L = TMath::Max(TMath::Max(jet_1.delta_R(bg_1),
					  jet_1.delta_R(bg_2)),
			       bg_1.delta_R(bg_2));
	      
	      if (r_L < res_cut) continue;

	      E3Cs[2]->Fill(r_L, jet_1.E() * bg_1.E() * bg_2.E());
	    }//! track pT cut
	  }//! m2
	}//! m1

	//! 2ND SIGNAL PTCL
	for (int y = x+1; y < int(constit.size()); y++) {
	  auto jet_2 = constit[y];

	  //! Get r_L
	  r_L= jet_1.delta_R(jet_2);

	  //! EEC_pt1 >1 GeV	      
	  if (jet_1.pt() > 1 && jet_2.pt() > 1) {
	    iptcls_ptgt1++;

	    //! EEC ============================================================
	    //! All combinations
	    EECs_pt1[0]->Fill(r_L, jet_1.E() * jet_2.E());
	      
	    //! breakdown of signal/breakdown contribution
	    if (jet_1.user_info<MyUserInfo>().origin() == 0) {	      

	      if (jet_2.user_info<MyUserInfo>().origin() == 0) { 
		//S+S
		EECs_pt1[1]->Fill(r_L, jet_1.E() * jet_2.E());
	      } else {
		//S+B
		EECs_pt1[2]->Fill(r_L, jet_1.E() * jet_2.E());
	      } 
	    } else { // c is not a signal ptcl
	      if (jet_2.user_info<MyUserInfo>().origin() == 0) {
		//B+S
		EECs_pt1[2]->Fill(r_L, jet_1.E() * jet_2.E());
	      } else {
		//B+B
		EECs_pt1[3]->Fill(r_L, jet_1.E() * jet_2.E());
	      }
	    }
	    
	    //! E3C SIGNAL + M1 ================================================
	    for (int im1 = 0; im1 < int(m1_ptcls.size()); im1++) {
	      auto bg_1 = m1_ptcls[im1];

	      if (bg_1.pt() > 1) { // jet_1.pt() > 1 && jet_2.pt() > 1 given 
		r_L = TMath::Max(TMath::Max(jet_1.delta_R(jet_2),
					    jet_1.delta_R(bg_1)),
				 jet_2.delta_R(bg_1));
		
		if (r_L < res_cut) continue;

		E3Cs[1]->Fill(r_L, jet_1.E() * jet_2.E() * bg_1.E());
	      }//! track pT cut 
	    }//! m1 loop
	    
	  }//! track pT cut on EEC

	  //! E3C; find the third particle (does not have to be unique)
	  for (int z = y+1; z < int(constit.size()); z++) {
	    auto jet_3 = constit[z];

	    // Get the maximum delta R between the 3 particles 
	    r_L = TMath::Max(TMath::Max(jet_1.delta_R(jet_2),
					jet_1.delta_R(jet_3)),
			                jet_2.delta_R(jet_3));
	    
	    if (r_L < res_cut) continue;

	    //! Fill the signal cone (SSS)
	    if (jet_1.pt() > 1 && jet_2.pt() > 1 && jet_3.pt() > 1) {
	      E3Cs[0]->Fill(r_L, jet_1.E() * jet_2.E() * jet_3.E());
	    }
	  }//! third particle
	}//! second particle
      }//! first particle

      //! End jet cone =========================================================

      //! M1 + M2 + M3
      for (int im1 = 0; im1 < int(m1_ptcls.size()); im1++) {
	auto bg_1 = m1_ptcls[im1];

	for (int im2 = 0; im1 < int(m2_ptcls.size()); im2++) {
	  auto bg_2 = m1_ptcls[im2];

	  for (int im3 = 0; im1 < int(m3_ptcls.size()); im3++) {
	    auto bg_3 = m3_ptcls[im3];

	    if (bg_1.pt() > 1 && bg_2.pt() > 1 && bg_3.pt() > 1) {
	      r_L = TMath::Max(TMath::Max(bg_1.delta_R(bg_2),
					  bg_1.delta_R(bg_3)),
			       bg_3.delta_R(bg_2));
	      
	      if (r_L < res_cut) continue;

	      E3Cs[3]->Fill(r_L, bg_1.E() * bg_2.E() * bg_3.E());
	    }//! track pT cut
	  }//! m3
	} //! m2
      }//! m1
      
      
    }//! jet loop
    cout << "All ptcls ptgt1: " << iptcls_ptgt1 << endl;
  }//! event loop

  //! Number of jets
  // cout << "Ensuring proper behavior of number of jets dummy histogram" << endl;
  // cout << "Number of jets (" << pT_lowbound << "-" << pT_highbound << " GeV): " << num_jets_debug << endl;
  // cout << "Number of histogram entries (" << pT_lowbound << "-" << pT_highbound << " GeV): " << int(num_jets->GetEntries()) << endl;
  // cout << "Debug info: numbers above should match exactly" << endl;
   
  //! Normalize the EEC_pt1 histograms
  // EEC_pt1->Scale(1./EEC_pt1->Integral());
  
  //! Normalize the E3C histograms
  // E3C->Scale(1./E3C->Integral());

  // //! Fill ratio histogram
  // TH1F *E3C_div_EEC_pt1 = (TH1F*) E3C->Clone("E3C_div_EEC_pt1");
  // E3C_div_EEC_pt1->Divide(EEC_pt1);
 
  //! Dummy histograms
  num_ptcls_event->Write();
  num_ptcls_jet->Write();
  num_jets->Write();
  
  //! Save the EEC_pt1(s) to file
  for (auto hist:EECs_pt1) {
    hist->Write();
  }

  //! Save the E3Cs to file
  for (auto hist:E3Cs) {
    hist->Write();
  }
  
  // for (auto hist:EECs_pt3) {
  //   hist->Write();
  // }
  
  // E3C_div_EEC_pt1->Write();
  
  //! Close file
  delete outFile;

  //! Done.
  return 0;
}
