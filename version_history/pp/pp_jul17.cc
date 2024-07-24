// This file contains code modified from PYTHIA examples main142.cc and main143.cc. 
// To ensure proper behavior please add an argument when running
// (e.g. ./new143_pileup 1) see lines 124-6 for usage

// This program simulates pp collisions at 5TeV with the option to add pileup.
// The goal of this simulation is to study the effects of background and
// subtraction methods in the 2 and 3-point correlators.

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
int events = 100000;

//! Number of jets
int num_jets_debug = 0;

//! Radius of jet
double radius = 0.5;

// Amount of pileup. Average number of inelastic pp collisions per event
// (=bunch-crossing). Set to zero to turn off pileup.
double mu = 0;
// int n_inel = gRandom->Poisson(mu);

//! Resolution cut
double res_cut = 0.01;

//! EPSILON value for r_L cuts in shape dep graph
double EPSILON = 0.005;

bool divJet = false;

using namespace Pythia8;

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

  //! Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile(Form("pp_%ito%iGeV_R05_N%dK_jet_charged_PU%i_xdivJetPT_poisson_%s.root",
				  pT_lowbound,pT_highbound,events/1000,100,argv[1]), 
			          "RECREATE");
  // TFile* outFile = new TFile(Form("pp_%ito%iGeV_R05_N%dK_jet_charged_PU%i_xdivJetPT_poisson.root",
  // 				  pT_lowbound,pT_highbound,events/1000,100), 
  // 			          "RECREATE"); // NOTE THIS IS FOR NON SLURM


  //! EEC histogram setup ======================================================
  int bins = 100;
  TString x_axis = "log10(#Deltar)";
  TString y_axis = "1/N_{pairs}";

  //! Title
  TString t1 = "EEC (";
  TString t2 = " GeV jets), both track p_{T} >";
  TString t3 = " GeV";  
  TString title = t1 + pT_lowbound + "-" + pT_highbound + t2;

  //! Energy of both particles > 1 GeV 
  TH1F *EEC = new TH1F("EEC", title + 1 + t3, bins, -3, 0);
 
  //! Update title and book histogram, cut on 1 particle only
  t2 = " GeV jets), 1 track p_{T} >";
  title = t1 + pT_lowbound + "-" + pT_highbound + t2;
  TH1F *EEC_halfCut = new TH1F("EEC_halfCut", title + 1 + t3, bins, -3, 0);

  //! Update title and book histogram, no track pT cut
  t2 = " GeV jets), no track p_{T} cut";
  title = t1 + pT_lowbound + "-" + pT_highbound + t2;
  TH1F *EEC_noCut = new TH1F("EEC_noCut", title + 1 + t3, bins, -3, 0);

  //! Set axis labels
  EEC->GetXaxis()->SetTitle(x_axis);
  EEC->GetYaxis()->SetTitle(y_axis);

  //! E3C histogram setup ======================================================
  t1 = "E3C (";
  t2 = " GeV jets), track p_{T} >";
  t3 = " GeV";  
  title = t1 + pT_lowbound + "-" + pT_highbound + t2;
  
  //! Energy of all particles must be above the threshold
  TH1F *E3C = new TH1F("E3C", title + 1 + t3, bins, -3, 0);

  //! Set axis labels
  E3C->GetXaxis()->SetTitle(x_axis);
  E3C->GetYaxis()->SetTitle(y_axis);

  //! Shape Dependence Histograms ==============================================
  //! Analyzing the correlation between particle energies and the shape of
  //! the triangle between the three particles in the E3C
  //! For more information see Fig 3 in:
  //! https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.130.051901

  // //! r_L ~ 0.15 (QGP region)
  // TH2F *shape_dep_15 = new TH2F("shape_dep_15",
  // 				"Shape Dependence of 3-pt Correlator (r_L ~ 0.15)",
  // 				20, 0.0, 1.0, 20, 0.0, TMath::PiOver2());
  // shape_dep_15->GetXaxis()->SetTitle("#xi");
  // shape_dep_15->GetYaxis()->SetTitle("#phi");

  // //! r_L ~ 0.002 (QGP region)
  // TH2F *shape_dep_20 = new TH2F("shape_dep_20",
  // 				"Shape Dependence of 3-pt Correlator (r_L ~ 0.20)",
  // 			        20, 0.0, 1.0, 20, 0.0, TMath::PiOver2());
  // shape_dep_20->GetXaxis()->SetTitle("#xi");
  // shape_dep_20->GetYaxis()->SetTitle("#phi");

  //! r_L ~ 0.25 (QGP region)
  TH2F *shape_dep_25 = new TH2F("shape_dep_25",
				"Shape Dependence of 3-pt Correlator (r_L ~ 0.25)",
			        20, 0.0, 1.0, 20, 0.0, TMath::PiOver2());
  shape_dep_25->GetXaxis()->SetTitle("#xi");
  shape_dep_25->GetYaxis()->SetTitle("#phi");

  // //! r_L ~ 0.30 (QGP region)
  // TH2F *shape_dep_30 = new TH2F("shape_dep_30",
  // 				"Shape Dependence of 3-pt Correlator (r_L ~ 0.30)",
  // 			        20, 0.0, 1.0, 20, 0.0, TMath::PiOver2());
  // shape_dep_30->GetXaxis()->SetTitle("#xi");
  // shape_dep_30->GetYaxis()->SetTitle("#phi");  

  // //! r_L ~ 0.35 (QGP region)
  // TH2F *shape_dep_35 = new TH2F("shape_dep_35",
  // 				"Shape Dependence of 3-pt Correlator (r_L ~ 0.35)",
  // 			        20, 0.0, 1.0, 20, 0.0, TMath::PiOver2());
  // shape_dep_35->GetXaxis()->SetTitle("#xi");
  // shape_dep_35->GetYaxis()->SetTitle("#phi");

  // //! r_L ~ 0.42 (QGP region, almost at jet radius cutoff)
  // TH2F *shape_dep_42 = new TH2F("shape_dep_42",
  // 				"Shape Dependence of 3-pt Correlator (r_L ~ 0.42)",
  // 			        20, 0.0, 1.0, 20, 0.0, TMath::PiOver2());
  // shape_dep_42->GetXaxis()->SetTitle("#xi");
  // shape_dep_42->GetYaxis()->SetTitle("#phi");

  //! Dummy histograms =====================================================

  //! Number of particles per event
  TH1F *num_ptcls_event = new TH1F("num_ptcls_event", "Number of Particles/Event",
			     1000, 0.0, 20000.0);

  //! Number of particles per jet
  TH1F *num_ptcls_jet = new TH1F("num_ptcls_jet", "Number of Particles/Jet",
			     70, 0.0, 70.0);
  
  //! Dummy histogram for per jet normalization
  TH1I *num_jets = new TH1I("num_jets", "Number of Jets", 1, 0, 1);

  //! Begin event loop. Generate event; skip if generation aborted.
  auto &event = pythia.event;
  for (int iEvent = 0; iEvent < events; ++iEvent) {
    if (!pythia.next()) continue;

    //! Identify particles. Jets are built from all stable particles after hadronization (particle-level jets). 
    std::vector<fastjet::PseudoJet> ch_particles, all_particles, pu_particles;
    int iptcls_ptgt1 = 0;

    
    //! Loop over all particles in the event 
    for (int i = 0; i < event.size(); ++i) {
      
      auto &p = event[i];
      if (p.isFinal()){
	all_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
	if (p.isCharged()) {
	  ch_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
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
        if (not p.isFinal()) continue;
        all_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        pu_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));

	if (p.isCharged()) {
	  ch_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
	}//! charged pileup particles
      }//! pileup event loop
    } //! main pileup loop

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

	//! EEC; find the second particle
	for (int y = x+1; y < int(constit.size()); y++) {
	  auto d = constit[y];

	  //! Get r_L
	  double r_L= c.delta_R(d);

	  //! Fill the EEC without any cut
	  EEC_noCut->Fill(TMath::Log10(r_L), c.E()*d.E() / jet.pt()*jet.pt());

	  if ((c.pt() > 1 || d.pt() > 1) && !(c.pt() > 1 && d.pt() > 1)) {
	    EEC_halfCut->Fill(TMath::Log10(r_L), c.E()*d.E() / jet.pt()*jet.pt());
	  }

	  //! Fill the EEC & E3C histograms by energy of both particles	      
	  if (c.pt() > 1 && d.pt() > 1) {
	    iptcls_ptgt1++;
	    if (divJet) {
	      EEC->Fill(TMath::Log10(r_L), c.E()*d.E() / jet.pt()*jet.pt());
	      
	      E3C->Fill(TMath::Log10(r_L), 3*c.E()*c.E()*d.E() / jet.pt()*jet.pt()*jet.pt());
	      E3C->Fill(TMath::Log10(r_L), 3*c.E()*d.E()*d.E() / jet.pt()*jet.pt()*jet.pt());
	    } else {
	      EEC->Fill(TMath::Log10(r_L), c.E()*d.E());
	      
	      E3C->Fill(TMath::Log10(r_L), 3*c.E()*c.E()*d.E());
	      E3C->Fill(TMath::Log10(r_L), 3*c.E()*d.E()*d.E());
	    }
	  }

	  //! E3C; find the third particle (does not have to be unique)
	  for (int z = y+1; z < int(constit.size()); z++) {
	    auto e = constit[z];

	    // Get the maximum delta R between the 3 particles 
	    r_L = TMath::Max(TMath::Max(c.delta_R(d), c.delta_R(e)), d.delta_R(e));
		
	    //! Fill the E3C histograms (energy of all particles)
	    if (c.pt() > 1 && d.pt() > 1 && e.pt() > 1) {
	      if (divJet) {
		E3C->Fill(TMath::Log10(r_L), 6*c.E()*d.E()*e.E() / jet.pt()*jet.pt()*jet.pt());
	      } else {
		E3C->Fill(TMath::Log10(r_L), 6*c.E()*d.E()*e.E());
	      }
	      
	    }
	    //! shape dependence of 3pt correlator =============================   		
	    //! Get the minimum delta R between the 3 particles (e allowed to be c or d)
	    double r_S = TMath::Min(TMath::Min(c.delta_R(d), c.delta_R(e)), d.delta_R(e)); 

	    //! Get the intermediate delta R between the 3 particles
	    double r_M;
	    if (r_S == c.delta_R(d)) { // cd is the shortest side
	      r_M = TMath::Min(c.delta_R(e), d.delta_R(e)); 
	    }
	    else if (r_S == c.delta_R(e)) { // ce is the shortest side
	      r_M = TMath::Min(c.delta_R(d), d.delta_R(e));
	    }
	    else { // de is the shortest side
	      r_M = TMath::Min(c.delta_R(d), c.delta_R(e));
	    }

	    //! Note: The 3 particles have to be unique from here on
	    //! Set resolution cut
	    if (r_S < res_cut || r_M < res_cut || r_L < res_cut) continue;

	    //! Ensure no dividing by zero (i.e. unique particles)
	    //! This should not be needed
	    if (r_M == 0 || r_S == 0) continue;
	    
	    //! equations defined at page 4 of
	    //! https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.130.051901
	    double xi = r_S/r_M;
	    double phi = TMath::ASin(TMath::Sqrt(1 - (TMath::Power((r_L - r_M), 2.0) / TMath::Power(r_S, 2.0))));

	    //! Calculate weight
	    double weight;
	    if (divJet) {
	      weight = (c.pt() * d.pt() * e.pt() / jet.pt()*jet.pt()*jet.pt());
	    } else {
	      weight = (c.pt() * d.pt() * e.pt()); 
	    }

	    //! set r_L cut, fill shape dependence histogram
	    if (TMath::Abs(r_L - 0.25) < (EPSILON)) {
	      shape_dep_25->Fill(xi, phi, weight);
	    }
	    // if (TMath::Abs(r_L - 0.15) < (EPSILON)) {
	    //   shape_dep_15->Fill(xi, phi, weight);
	    // }
	    // else if (TMath::Abs(r_L - 0.20) < (EPSILON)) {
	    //   shape_dep_20->Fill(xi, phi, weight);
	    // }
	    // else if (TMath::Abs(r_L - 0.25) < (EPSILON)) {
	    //   shape_dep_25->Fill(xi, phi, weight);
	    // }
	    // else if (TMath::Abs(r_L - 0.30) < (EPSILON)) {
	    //   shape_dep_30->Fill(xi, phi, weight);
	    // }
	    // else if (TMath::Abs(r_L - 0.35) < (EPSILON)) {
	    //   shape_dep_35->Fill(xi, phi, weight);
	    // }
	    // else if (TMath::Abs(r_L - 0.42) < (EPSILON)) {
	    //   shape_dep_42->Fill(xi, phi, weight);
	    // } 

	  }//! third particle
	  cout << "End third particle loop" << endl;
	}//! second particle
	cout << "End second particle loop" << endl;
      }//! first particle
      cout << "End first particle loop" << endl;
    }//! jet loop
    cout << "All ptcls ptgt1: " << iptcls_ptgt1 << endl;
  }//! event loop

  //! Number of jets
  // cout << "Ensuring proper behavior of number of jets dummy histogram" << endl;
  // cout << "Number of jets (" << pT_lowbound << "-" << pT_highbound << " GeV): " << num_jets_debug << endl;
  // cout << "Number of histogram entries (" << pT_lowbound << "-" << pT_highbound << " GeV): " << int(num_jets->GetEntries()) << endl;
  // cout << "Debug info: numbers above should match exactly" << endl;
   
  //! Normalize the EEC histograms
  // EEC->Scale(1./EEC->Integral());
  
  //! Normalize the E3C histograms
  // E3C->Scale(1./E3C->Integral());

  //! Fill ratio histogram
  TH1F *E3C_div_EEC = (TH1F*) E3C->Clone("E3C_div_EEC");
  E3C_div_EEC->Divide(EEC);

  //! Normalize the shape dependence histogram and save it to file
  // shape_dependence->Scale(1./shape_dependence->Integral()); 
  // shape_dependence->Scale(1./num_jets_debug); 
  // shape_dep_15->Write();
  // shape_dep_20->Write();  
  shape_dep_25->Write();
  // shape_dep_30->Write();
  // shape_dep_35->Write();
  // shape_dep_42->Write();

  //! Dummy histograms
  num_ptcls_event->Write();
  num_ptcls_jet->Write();
  num_jets->Write();
  
  //! Save the EEC(s) to file
  EEC->Write();
  EEC_halfCut->Write();
  EEC_noCut->Write();
  
  //! Save the E3C(s) to file
  E3C->Write();
  
  E3C_div_EEC->Write();
  
  //! Close file
  delete outFile;

  //! Done.
  return 0;
}
