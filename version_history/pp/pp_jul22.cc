// This file contains code modified from PYTHIA examples main142.cc and main143.cc. 
// To ensure proper behavior please add an argument when running
// (e.g. ./new143_pileup 1) see lines 124-6 for usage

// This program simulates pp collisions at 5TeV with the option to add pileup.
// The goal of this simulation is to study the effects of background and
// subtraction methods in the 2 and 3-point correlators.

// The output .root file contains the following histograms:
// -- 1 EEC_pt1 histogram
// -- 1 E3C histogram
// -- 6 E3C Shape Dependence histograms
// -- 1 number of jets histogram (dummy)

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// Angantyr
#include "Pythia8/HeavyIons.h"

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
int events = 1000;

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
bool signal_ptcl = false;

using namespace Pythia8;

class MyUserInfo : public fastjet::PseudoJet::UserInfoBase{
public:
  // default ctor
  //  - pdg_id        the PDG id of the particle
  //  - is_signal     whether it is signal (TRUE) or background (FALSE)
  MyUserInfo(const int & pdg_id_in, const bool & is_signal_in) :
    _pdg_id(pdg_id_in), _is_signal(is_signal_in){}
  
  /// access to the PDG id
  int pdg_id() const { return _pdg_id;}
   
  /// access to is signal
  bool is_signal() const { return _is_signal;}
   
protected:
  int _pdg_id;         // the associated pdg id
  bool _is_signal;     // the associated property
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
				  pT_lowbound,pT_highbound,events/1000,200,argv[1]), 
			          "RECREATE");
  //! For interactive:
  // TFile* outFile = new TFile(Form("pp_%ito%iGeV_R05_N%dK_jet_charged_PU%i_xdivJetPT_poisson.root",
  // 				  pT_lowbound,pT_highbound,events/1000,100), 
  // 			          "RECREATE"); // NOTE THIS IS FOR NON SLURM


  //! EEC_pt1 histogram setup ======================================================
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

  TH1F *EEC_pt3 = new TH1F("EEC_pt3", title + 3 + t3 + " all comb", bins, low, high);
  TH1F *EEC_pt3_SS = new TH1F("EEC_pt3_SS", title + 3 + t3 + " S+S", bins, low, high);
  TH1F *EEC_pt3_SB = new TH1F("EEC_pt3_SB", title + 3 + t3 + " S+B", bins, low, high);
  TH1F *EEC_pt3_BB = new TH1F("EEC_pt3_BB", title + 3 + t3 + " B+B", bins, low, high);

  TH1F *EECs_pt1[] = {EEC_pt1, EEC_pt1_SS, EEC_pt1_SB, EEC_pt1_BB};
  TH1F *EECs_pt3[] = {EEC_pt3, EEC_pt3_SS, EEC_pt3_SB, EEC_pt3_BB};
  
  // //! Update title and book histogram, cut on 1 particle only
  // t2 = " GeV jets), 1 track p_{T} >";
  // title = t1 + pT_lowbound + "-" + pT_highbound + t2;
  // TH1F *EEC_pt1_halfCut = new TH1F("EEC_pt1_halfCut", title + 1 + t3, bins, -3, 0);

  //! Set axis labels
  EEC_pt1->GetXaxis()->SetTitle(x_axis);
  EEC_pt1->GetYaxis()->SetTitle(y_axis);

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
  for (int iEvent = 0; iEvent < events; ++iEvent) {
    if (!pythia.next()) continue;

    //! Identify particles. Jets are built from all stable particles after hadronization (particle-level jets).
    //! ch: charged (S+B); sig: signal, S; pu: pileup, B
    std::vector<fastjet::PseudoJet> ch_particles, all_particles, sig_particles, pu_particles;
    int iptcls_ptgt1 = 0;
    
    //! Loop over all particles in the event 
    for (int i = 0; i < event.size(); ++i) {
      auto &p = event[i];
      if (p.isFinal()){

	//! label this as a signal particle
	signal_ptcl = true;
	auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	cur_ptcl.set_user_info(new MyUserInfo(p.id(), signal_ptcl));

	all_particles.push_back(cur_ptcl);
	if (p.isCharged()) {
	  ch_particles.push_back(cur_ptcl);
	  sig_particles.push_back(cur_ptcl);
	}//! charged particles
      }//! final particles
    }

    //! Add in pileup
    int n_inel = gRandom->Poisson(mu);
    // printf("Overlaying particles from %i pileup interactions!\n", n_inel);
    for (int i_pu= 0; i_pu<n_inel; ++i_pu) {
      if (!pythiaPU.next()) continue;
      for (int i = 0; i < pythiaPU.event.size(); ++i) {
        auto &p = pythiaPU.event[i];
        if (p.isFinal()) {

	  //! Label this as not a signal particle
	  signal_ptcl = false;
	  auto cur_ptcl = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e());
	  cur_ptcl.set_user_info(new MyUserInfo(p.id(), signal_ptcl));
	  
	  all_particles.push_back(cur_ptcl);
	  pu_particles.push_back(cur_ptcl);

	  if (p.isCharged()) {
	    ch_particles.push_back(cur_ptcl);
	  }//! charged pileup particles
	}//! final pileup particles
      }//! pileup event loop
    }//! main pileup loop

    if (mixed_cones) {
      
    }
    
    //! Debugging pileup:
    cout << "All ch ptcls: " << ch_particles.size() << endl;
    num_ptcls_event->Fill(ch_particles.size());
    cout << "All pu ptcls: " << pu_particles.size() << endl;

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
      
      vector <fastjet::PseudoJet> constit = jet.constituents();

      //! Debugging pileup
      num_ptcls_jet->Fill(constit.size());
      num_jets->Fill(0);

      //! For each particle within the radius:
      for (int x = 0; x < int(constit.size()); x++) {
	auto c = constit[x];

	//! EEC_pt1; find the second particle
	for (int y = x+1; y < int(constit.size()); y++) {
	  auto d = constit[y];

	  //! Get r_L
	  double r_L= c.delta_R(d);

	  //! EEC_pt1 >1 GeV	      
	  if (c.pt() > 1 && d.pt() > 1) {
	    iptcls_ptgt1++;
	    //! Fill the general EEC_pt1 histogram
	    EECs_pt1[0]->Fill(r_L, c.E()*d.E());
	      
	    //! Fill the EEC_pt1 by pair type
	    if (c.user_info<MyUserInfo>().is_signal()) {	      

	      if (d.user_info<MyUserInfo>().is_signal()) { 
		//S+S
		EECs_pt1[1]->Fill(r_L, c.E()*d.E());
	      } else {
		//S+B
		EECs_pt1[2]->Fill(r_L, c.E()*d.E());
	      } 
	    } else { // c is not a signal ptcl
	      if (d.user_info<MyUserInfo>().is_signal()) {
		//B+S
		EECs_pt1[2]->Fill(r_L, c.E()*d.E());
	      } else {
		//B+B
		EECs_pt1[3]->Fill(r_L, c.E()*d.E());
	      }
	    }
	  }//! track pT cut

	  //! EEC_pt1 >3 GeV
	  if (c.pt() > 3 && d.pt() > 3) {
	    
	    //! Fill the general EEC_pt1 histogram
	    EECs_pt3[0]->Fill(r_L, c.E()*d.E());
	      
	    //! Fill the EEC_pt1 by pair type
	    if (c.user_info<MyUserInfo>().is_signal()) {
	      if (d.user_info<MyUserInfo>().is_signal()) {
		//S+S
		EECs_pt3[1]->Fill(r_L, c.E()*d.E());
	      } else {
		//S+B
		EECs_pt3[2]->Fill(r_L, c.E()*d.E());
	      } 
	    } else { // c is not a signal ptcl
	      if (d.user_info<MyUserInfo>().is_signal()) {
		//B+S
		EECs_pt3[2]->Fill(r_L, c.E()*d.E());
	      } else {
		//B+B
		EECs_pt3[3]->Fill(r_L, c.E()*d.E());
	      }
	    }
	  }
	  
	}//! second particle
	// cout << "End second particle loop" << endl;
      }//! first particle
      // cout << "End first particle loop" << endl;
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
  for (auto hist:EECs_pt3) {
    hist->Write();
  }
  
  // E3C_div_EEC_pt1->Write();
  
  //! Close file
  delete outFile;

  //! Done.
  return 0;
}
