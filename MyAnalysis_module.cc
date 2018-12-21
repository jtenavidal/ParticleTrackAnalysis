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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraph.h"
#include "Math/PositionVector3D.h"

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
  std::string fTruthLabel, m_particleLabel, m_hitfinderLabel, m_recotrackLabel;
  int fPDG_Code, fTrack_ID, fNumDaughters, fFirstDaughter, fDaughter; 
  float fTrueParticleEnergy, fMass;
  float fpx, fpy, fpz, fpt, fp; // momentum variables
  simb::MCTrajectory True_trajectory ;
  TLorentzVector MC_Track_Position ;
  double fMCLength ;
  double fTrack_position_x, fTrack_position_y, fTrack_position_z, fTrack_position_T ;


  // Reco information
  int r_pdg_primary, r_nu_daughters ;
  int r_mu_daughters, r_pi_daughters, r_e_daughters, r_p_daughters, r_other_daughters;
  recob::TrackTrajectory primary_trajectory ;
  double rVertex_x, rVertex_y, rVertex_z, rEnd_x, rEnd_y, rEnd_z, rLength, rMomentum ;
};


TrackID::MyAnalysis::MyAnalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  // Get the .fcl parameters
  fTruthLabel = p.get<std::string>("TruthLabel");
  m_particleLabel = p.get<std::string>("PFParticleModule","pandora");
  m_hitfinderLabel = p.get<std::string>("HitFinderModule","linecluster");
  m_recotrackLabel = p.get<std::string>("RecoTrackLabel","pandoraTrack");

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
	  fMCLength = True_trajectory.TotalLength() ;
	  MC_Track_Position = True_trajectory.Position( True_trajectory.size() ) ;
	  
	}
      }
      /*      for ( unsigned int j = 0; j<True_trajectory.size(); ++j){
	fTrack_position_x = True_trajectory.X( j ) ;                                  // Need to fill the tree in the right way so that every hit is stored
	fTrack_position_y = True_trajectory.Y( j ) ;
	fTrack_position_z = True_trajectory.Z( j ) ;
	fTrack_position_T = True_trajectory.T( j ) ;
	}*/

      // This should be a function called Print : need to understand what it is printing exactly. Obviously a track, but is it the primary? 

      TH3D *h_track = new TH3D("h_track", " Particle Track ", True_trajectory.size(), True_trajectory.X(0), True_trajectory.X(True_trajectory.size()), True_trajectory.size(), True_trajectory.Y(0), True_trajectory.Y(True_trajectory.size()), True_trajectory.size(), True_trajectory.Z(0), True_trajectory.Z(True_trajectory.size()));

      for( unsigned int i = 0; i < True_trajectory.size(); ++i ){
	h_track-> Fill(True_trajectory.X(i), True_trajectory.Y(i), True_trajectory.Z(i), True_trajectory.E(i));
      }
      gStyle->SetPalette(55);
      gStyle->SetNumberContours(250);

      TCanvas *c = new TCanvas();
      h_track->SetLineColor(2);
      h_track->GetXaxis()->SetTitle("X");
      h_track->GetYaxis()->SetTitle("Y");
      h_track->GetXaxis()->SetTitle("Z");
      h_track->Draw("hist");
      h_track->Draw("BOX2Z");
      //  c->SaveAs((path+"_track.root").c_str());
      c->SaveAs("particle_track.root");
      //c->Clear();
    }
  }

/**************************************************************************************************
 *  RECO INFORMATION
 *************************************************************************************************/
  // Get PFParticle Handle
  art::Handle< std::vector< recob::PFParticle > > pfParticleHandle ;
  e.getByLabel(m_particleLabel, pfParticleHandle ) ;

  //Find the hits associated to the reconstructed PFParticle
  art::Handle< std::vector< recob::Hit > > hitListHandle ;
  e.getByLabel(m_hitfinderLabel, hitListHandle);

  //Find the reco tracks
  art::Handle< std::vector< recob::Track > > trackHandle ;
  e.getByLabel(m_recotrackLabel, trackHandle ) ;
  
  // Get track associations with PFParticles from Pandora. Find all possible tracks associated to an event
  art::FindManyP< recob::Track > findTracks( pfParticleHandle, e, m_recotrackLabel );

  r_mu_daughters = 0;
  r_pi_daughters = 0;
  r_e_daughters = 0;
  r_p_daughters = 0;
  r_other_daughters = 0;
  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid()){

    for( unsigned int i = 0 ; i < pfParticleHandle->size(); ++i ){
      art::Ptr< recob::PFParticle > pfparticle( pfParticleHandle, i ) ; // Point to particle i 
    
      if( pfparticle->IsPrimary() == 1 ){
	r_pdg_primary = pfparticle->PdgCode() ;
	r_nu_daughters = pfparticle->NumDaughters();
      } else {
      	if( pfparticle->PdgCode() == 13 ) {
	  r_mu_daughters ++ ;
	} else if ( pfparticle->PdgCode() == 211 ) {
	  r_pi_daughters ++ ;
	} else if ( pfparticle->PdgCode() == 11 ) {
	  r_e_daughters ++ ;
	} else {
	  r_other_daughters ++ ;
	  }
      }
	
	if ( findTracks.at(i).size()!=0 ){
	  std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(i);

	  for( unsigned int j = 0 ; j < track_f.size() ; ++j ){
	    //	    std::cout<<"x= "<<track_f[j]->TrajectoryPoint( 0 ).position.X()<<std::endl;
	    //	    std::cout<<"chi2= "<<track_f[j]->Chi2()<<std::endl;
	    rVertex_x = track_f[j]->Vertex( ).X() ;
	    rVertex_y = track_f[j]->Vertex( ).Y() ;
	    rVertex_z = track_f[j]->Vertex( ).Z() ;
	    rEnd_x    = track_f[j]->End( ).X() ;
	    rEnd_y    = track_f[j]->End( ).Y() ;
	    rEnd_z    = track_f[j]->End( ).Z() ;
	    rLength   = track_f[j]->Length() ;
	    rMomentum = track_f[j]->MomentumAtPoint( 0 ) ; // need to clarify which momentum is it.
	  }
	} 
    }
  }


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
  fNumDaughters = -999 ;
  fFirstDaughter = -999 ;
  fDaughter = -999 ;
  fMCLength = -999 ;
  fTrack_position_x = -999. ;
  fTrack_position_y = -999. ;
  fTrack_position_z = -999. ;
  fTrack_position_T = -999.;

  // reco information
  r_pdg_primary   = 0 ;
  r_nu_daughters  = 0 ;
  r_mu_daughters = 0 ;
  r_pi_daughters = 0 ;
  r_e_daughters = 0 ;
  r_p_daughters = 0 ;
  r_other_daughters = 0 ;
  rVertex_x = -999. ;
  rVertex_y = -999. ;
  rVertex_z = -999. ;
  rEnd_x = -999. ;
  rEnd_y = -999. ;
  rEnd_z = -999. ;
  rLength = -999. ;
  rMomentum = -999. ;


  // Declare trees and branches
  event_tree      = new TTree( "event_tree",           "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree( "mcparticle_tree",      "MC tree:    True Particle track information");
  recotrack_tree  = new TTree( "recoparticle_tree",    "Reco tree: reconstructed information of the tracks, hit level included");

  event_tree      -> Branch( "event_id",          &event_id, "event_id/I");

  /**
     MC PARTICLE TREE BRANCHES :
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
  mcparticle_tree -> Branch( "MC_Length",              &fMCLength,      "Length/D");
  /*  mcparticle_tree -> Branch( "MC_track_position_x",              &fTrack_position_x,      "Position_x/D");
  mcparticle_tree -> Branch( "MC_track_position_y",              &fTrack_position_y,      "Position_y/D");
  mcparticle_tree -> Branch( "MC_track_position_z",              &fTrack_position_z,      "Position_z/D");
  mcparticle_tree -> Branch( "MC_track_position_T",              &fTrack_position_T,      "Position_T/D");
  */

  /**
     RECONSTRUCTED PARTICLE TREE BRANCHES :
   */
  recotrack_tree  -> Branch( "event_id",          &event_id,          "event_id/I");
  recotrack_tree  -> Branch( "pdg_primary",       &r_pdg_primary,     "pdg_primary/I");
  recotrack_tree  -> Branch( "nu_daughters",      &r_nu_daughters,    "nu_daughters/I");
  recotrack_tree  -> Branch( "mu_daughters",      &r_mu_daughters,    "mu_daughters/I");
  recotrack_tree  -> Branch( "pi_daughters",      &r_pi_daughters,    "pi_daughters/I");
  recotrack_tree  -> Branch( "e_daughters",       &r_e_daughters,     "e_daughters/I");
  recotrack_tree  -> Branch( "p_daughters",       &r_p_daughters,     "p_daughters/I");
  recotrack_tree  -> Branch( "other_daughters",   &r_other_daughters, "others_daughters/I");
  recotrack_tree  -> Branch( "Vertex_x",          &rVertex_x,         "rVertex_x/D");
  recotrack_tree  -> Branch( "Vertex_y",          &rVertex_y,         "rVertex_y/D");
  recotrack_tree  -> Branch( "Vertex_z",          &rVertex_z,         "rVertex_z/D");
  recotrack_tree  -> Branch( "End_x",             &rEnd_x,            "rEnd_x/D");
  recotrack_tree  -> Branch( "End_y",             &rEnd_y,            "rEnd_y/D");
  recotrack_tree  -> Branch( "End_z",             &rEnd_z,            "rEnd_z/D");
  recotrack_tree  -> Branch( "Length",            &rLength,           "rLength/D");
  recotrack_tree  -> Branch( "Momentum",          &rMomentum,         "rMomentum/D");

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
