////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v2_11_03)
// File:        MyAnalysis_module.cc
//
// Generated at Tue Dec 11 09:43:00 2018 by Julia Tena Vidal using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"

#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TLorentzVector.h"

namespace TrackID {
  class MyAnalysis;
}


class TrackID::MyAnalysis : public art::EDAnalyzer {
public:
  explicit MyAnalysis(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyAnalysis(MyAnalysis const &) = delete;
  MyAnalysis(MyAnalysis &&) = delete;
  MyAnalysis & operator = (MyAnalysis const &) = delete;
  MyAnalysis & operator = (MyAnalysis &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  //                    ->  Added functions (*)
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // Tree members
  TTree * event_tree, * mcparticle_tree, * recotrack_tree ; 
  int event_id ; 

  // Truth information
  std::string fTruthLabel ; 
  int fPDG_Code, fTrack_ID, fNumDaughters, fFirstDaughter, fDaughter; 
  float fTrueParticleEnergy, fMass;
  float fpx, fpy, fpz, fpt, fp; // momentum variables
  simb::MCTrajectory True_trajectory ;
};


TrackID::MyAnalysis::MyAnalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  // Get the .fcl parameters
  fTruthLabel = p.get<std::string>("TruthLabel");
  this->reconfigure(p);
}

void TrackID::MyAnalysis::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  event_id = e.id().event();

  if( !e.isRealData()){
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()){
      // Loop over primary particle
      for( unsigned int t = 0; t < mcParticles->size(); ++t ){
	const simb::MCParticle trueParticle = mcParticles->at(t) ;
	if(trueParticle.Process() == "primary"){
	  fTrack_ID = trueParticle.TrackId() ;
	  fPDG_Code = trueParticle.PdgCode() ;
	  fTrueParticleEnergy = trueParticle.E() ;
	  fMass = trueParticle.Mass() ;
	  fpx = trueParticle.Px() ;
	  fpy = trueParticle.Py() ;
	  fpz = trueParticle.Pz() ;
	  fpt = trueParticle.Pt() ;
	  fp  = trueParticle.P() ;
	  fNumDaughters = trueParticle.NumberDaughters() ;
	  fFirstDaughter = trueParticle.FirstDaughter() ;
	  //	  fDaughter = trueParticle.Daughter() ;
	  True_trajectory = trueParticle.Trajectory() ;
	}
      }
    }
  }
  std::cout<<fNumDaughters<<std::endl;
  std::cout<<fFirstDaughter<<std::endl;
  std::cout<<fDaughter<<"\n"<<std::endl;

  // FILL TREES
  event_tree      -> Fill();
  mcparticle_tree -> Fill();
  recotrack_tree  -> Fill();

}

void TrackID::MyAnalysis::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of required member function here.
  // Here you add an external fcl file to change configuration
}

void TrackID::MyAnalysis::beginJob( )
{
  // Implementation of required member function here.
  // Define default for parameters and create variables and trees

  // Event ID
  event_id = 0 ;

  // Truth Information 
  fTrack_ID = -999 ;
  fTrueParticleEnergy = -999. ; 
  fPDG_Code = -999 ;
  fMass = -1  ;
  fpx = -999. ;
  fpy = -999. ;
  fpz = -999. ;
  fpt = -999. ;
  fp  = -999. ;
  fNumDaughters = -1 ;
  fFirstDaughter = -1 ;
  fDaughter = -1 ;
  // Declare trees and branches
  event_tree      = new TTree( "event_tree",           "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree( "mcparticle_tree",      "MC tree:    True Particle track information");
  recotrack_tree  = new TTree( "recoparticle_tree",    "Reco tree: reconstructed information of the tracks, hit level included");

  /*
   Event Tree contains general information of the event regardless of its nature
   [ ] event id
   [ ] time now
   [ ] mother particle pdg code (particle beam , sanity check )
   [ ] starting point
   [ ] End point 
   [ ] number of particles (including daughter particles) per type 
   [ ] number of tracks
   [ ] number of showers 
  */

  event_tree      -> Branch( "event_id",          &event_id, "event_id/I");

  /*
    MC Particle tree contains information at the truth level
    Must have:
    [X] ID
    [X] PDG code
    [X] Mass
    [X] Momentum 
    [ ] direction -> Can have a TLorenzVector
    [X] Energy 
    [ ] Process
    [ ] Keep hiearchy information : possible at this level?
       [X] NumberDaughters() {int type}
       [ ] Daughter( ) {int type}
    [?] True hit level information -> detector effects, smearing (future?)
       [ ] simb::MCTrajectory Trajectory 
   */
  mcparticle_tree -> Branch( "event_id",               &event_id,            "event_id/I");
  mcparticle_tree -> Branch( "Track_ID",               &fTrack_ID,           "Track_id/I");
  mcparticle_tree -> Branch( "trueEnergy",             &fTrueParticleEnergy, "TrueParticleEnergy/F");
  mcparticle_tree -> Branch( "PDG_Code",               &fPDG_Code,           "PDG_Code/I");
  mcparticle_tree -> Branch( "Mass",                   &fMass,               "Mass/F");
  mcparticle_tree -> Branch( "Px",                     &fpx,                 "px/F");
  mcparticle_tree -> Branch( "Py",                     &fpy,                 "py/F");
  mcparticle_tree -> Branch( "Pz",                     &fpz,                 "pz/F");
  mcparticle_tree -> Branch( "Pt",                     &fpt,                 "pt/F");
  mcparticle_tree -> Branch( "P",                      &fp,                  "p/F");
  mcparticle_tree -> Branch( "Num_Daughters",          &fNumDaughters,       "num_d/I");
  mcparticle_tree -> Branch( "First_Daughter",         &fFirstDaughter,      "First_d/I");
  //  mcparticle_tree -> Branch( "Daughter",               &fDaughter,           "d/I");

  /*
    Reco track tree contains information of the reconstructed track, including hit level information
    [ ] Chi2 hipotesis for: muons, pions, kaons and protons (the last two should not be seen but just in case)
    [ ] dE/dx (must be able to split it if needed)
    [ ] Residual range
    [ ] Lenght 
    [ ] Kinetic energy and missing energy 
    [ ] Hit information: number of hits, hit_tpc, hit_plane
    [ ] Precise hit information: 3D track reconstructed ( position, time and calorimetry )
    [ ] Keep hiearchy information
   */

  recotrack_tree  -> Branch( "event_id",          &event_id, "event_id/I");

  // Set directories
  event_tree->SetDirectory(0);
  mcparticle_tree->SetDirectory(0);
  recotrack_tree->SetDirectory(0);

}

void TrackID::MyAnalysis::endJob( )
{
  // Implementation of required member function here.
  // Write and close files. Other comments also fit here 
  TFile file("output_trees_trackID.root" , "RECREATE");
  event_tree      -> Write();
  mcparticle_tree -> Write();
  recotrack_tree  -> Write();

  file.Write();
  file.Close();

  delete event_tree ; 
  delete mcparticle_tree ; 
  delete recotrack_tree ; 
}

DEFINE_ART_MODULE(TrackID::MyAnalysis)
