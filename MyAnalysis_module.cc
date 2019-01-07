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
#include <sstream>
#include <string> 
#include "TVector3.h"
#include "math.h"
//#include "/nashome/j/jtenavid/Desktop/Exercises/PiMuStudy/MyWorkingDir/build_slf6.x86_64/sbndcode/sbndcode/ParticleTrackAnalysis/TrackFitterJulia.h" // Including my class 

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

  // Labels
  std::string fTruthLabel, m_particleLabel, m_hitfinderLabel, m_recotrackLabel, m_recoPIDLabel, m_recoCaloLabel;

  // Tree members
  TTree * event_tree, * mcparticle_tree, * recotrack_tree, * recoTrackInfo_tree ; 
  int event_id ; 

  // Truth information
  int fPDG_Code, fTrack_ID, fNumDaughters, fDaughter_mu, fDaughter_pi, fDaughter_e, fDaughter_p, fDaughter_n, fDaughter_photon, fDaughter_other ;
  float fTrueParticleEnergy, fMass;
  float fpx, fpy, fpz, fpt, fp; // momentum variables
  simb::MCTrajectory True_trajectory ;
  TLorentzVector MC_Track_Position ;
  double fMCLength, fTrack_position_x, fTrack_position_y, fTrack_position_z, fTrack_position_T ;

  // Reco information
  int r_pdg_primary, r_nu_daughters ;
  int r_mu_daughters, r_pi_daughters, r_e_daughters, r_p_daughters, r_n_daughters, r_photon_daughters, r_other_daughters;
  recob::TrackTrajectory primary_trajectory ;
  double rVertex_x, rVertex_y, rVertex_z, rEnd_x, rEnd_y, rEnd_z;
  int rnu_hits ;
  double r_chi2_mu, r_chi2_pi, r_chi2_p, r_PIDA, r_missenergy, r_KineticEnergy, r_Range ;
  float rLength, rMomentum ;
  std::vector< float > dEdx ;
  std::vector< std::vector< double >  > mu_f, dev_f;

  //Track example information to test future Track class. This is the information needed for the Track class. 
  // Saving it into a root file to work offline...
  double tr_x, tr_y, tr_z, tr_dEdx, tr_dQdx ; //info per hit

  // Functions
  void SaveMCTrack( simb::MCTrajectory const & mc_track, std::string const & path ) const ;
  void SaveRecoTrack( art::Ptr<recob::Track> const & reco_track, std::string const & path ) const ;
  std::vector< std::vector<double> > Straight(  art::Ptr<recob::Track> const & reco_track ) ;  
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
  m_recoCaloLabel = p.get<std::string>("RecoCaloLabel","pandoraCalo");
  m_recoPIDLabel = p.get<std::string>("RecoPIDLabel","pandoraPid");

  this->reconfigure(p);
}

void TrackID::MyAnalysis::analyze(art::Event const & e)
{
  event_id = e.id().event();
  std::stringstream t_path, r_path ;
  t_path << "Histograms/Truth/eid_"<<event_id<<"_truth_track" ;
  r_path << "Histograms/Reco/eid_"<<event_id<<"_reco_track" ;
  std::string truth_path = t_path.str();
  std::string reco_path = r_path.str();

  if( !e.isRealData()){
    /**************************************************************************************************
     *  MC INFORMATION
     *************************************************************************************************/
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);

    if(mcParticles.isValid()){
      // Loop over truth info
      fDaughter_mu = 0 ;
      fDaughter_pi = 0 ;
      fDaughter_e = 0 ;
      fDaughter_p = 0 ;
      fDaughter_n = 0 ;
      fDaughter_photon = 0 ;
      fDaughter_other = 0 ;
      
      for( unsigned int t = 0; t < mcParticles->size(); ++t ){
	
	const simb::MCParticle trueParticle = mcParticles->at(t) ;

	if( trueParticle.PdgCode() >= 1000018038 ) continue ; // Cut on PDG codes which refer to elements (Argon30 and above)
 
	if(trueParticle.Process() == "primary" ){
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
	  True_trajectory = trueParticle.Trajectory() ;
	  fMCLength = True_trajectory.TotalLength() ;
	  
	} else { // secondary particle information 
	  if      ( trueParticle.PdgCode() == 13   ) { ++fDaughter_mu ; }
	  else if ( trueParticle.PdgCode() == 211  ) { ++fDaughter_pi ; } 
	  else if ( trueParticle.PdgCode() == 11   ) { ++fDaughter_e ; } 
	  else if ( trueParticle.PdgCode() == 2212 ) { ++fDaughter_p ; }
	  else if ( trueParticle.PdgCode() == 2112 ) { ++fDaughter_n ; }
	  else if ( trueParticle.PdgCode() == 22  )  { ++fDaughter_photon ; }
	  else                                       { ++fDaughter_other ; }
	}
      }
      SaveMCTrack( True_trajectory , truth_path ) ;
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
  r_n_daughters = 0;
  r_photon_daughters = 0;
  r_other_daughters = 0;
	   
  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid()){

    for( unsigned int i = 0 ; i < pfParticleHandle->size(); ++i ){
      art::Ptr< recob::PFParticle > pfparticle( pfParticleHandle, i ) ; // Point to particle i 
    
      if( pfparticle->IsPrimary() == 1 ){
	r_pdg_primary = pfparticle->PdgCode() ;
	r_nu_daughters = pfparticle->NumDaughters();
      } else {
      	if        ( pfparticle->PdgCode() == 13   ) { ++ r_mu_daughters ;
	} else if ( pfparticle->PdgCode() == 211  ) { ++ r_pi_daughters ;
	} else if ( pfparticle->PdgCode() == 11   ) { ++ r_e_daughters ;
	} else if ( pfparticle->PdgCode() == 2212 ) { ++ r_p_daughters ;
	} else if ( pfparticle->PdgCode() == 2112 ) { ++ r_n_daughters ;
	} else if ( pfparticle->PdgCode() == 22   ) { ++ r_photon_daughters ;
	} else                                      { ++ r_other_daughters ;
	  }
      }
	
	if ( findTracks.at(i).size()!=0 ){
	  std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(i);
	  art::FindManyP< recob::Hit > findHits (  trackHandle, e, m_recotrackLabel ) ;
	  art::FindManyP< anab::Calorimetry > findCalorimetry ( trackHandle, e, m_recoCaloLabel );
	  art::FindManyP< anab::ParticleID > findPID ( trackHandle, e, m_recoPIDLabel );
	  // Loop over tracks per event
	  for( unsigned int j = 0 ; j < track_f.size() ; ++j ){
	    rVertex_x = track_f[j]->Vertex( ).X() ;
	    rVertex_y = track_f[j]->Vertex( ).Y() ;
	    rVertex_z = track_f[j]->Vertex( ).Z() ;
	    rEnd_x    = track_f[j]->End( ).X() ;
	    rEnd_y    = track_f[j]->End( ).Y() ;
	    rEnd_z    = track_f[j]->End( ).Z() ;
	    rLength   = track_f[j]->Length() ;
	    rMomentum = track_f[j]->MomentumAtPoint( 0 ) ; // need to clarify which momentum is it.
	    SaveRecoTrack( track_f[j], reco_path );
	    
	    // Get track based variables
	    std::vector< art::Ptr<recob::Hit> > hit_f        = findHits.at(track_f[j]->ID()); 
	    std::vector< art::Ptr<anab::Calorimetry> > cal_f = findCalorimetry.at(track_f[j]->ID());
	    std::vector< art::Ptr<anab::ParticleID> > pid_f  = findPID.at(track_f[j]->ID());

	    //Loop over PID associations 
	    for ( unsigned int k = 0 ; k < pid_f.size() ; ++k ){
	      if( !pid_f[k] ) continue ;
	      if( !pid_f[k]->PlaneID().isValid) continue ;
	      // only look at collection plane for dEdx information
	      if( pid_f[k]->PlaneID().Plane != 2 ) continue ;

	      //Loop over calo information also in collection plane
	      for ( unsigned int n = 0 ; n < cal_f.size() ; ++n ) {
		if( !cal_f[n] ) continue ;
		if( !cal_f[n]->PlaneID().isValid) continue ;
		if( cal_f[n]->PlaneID().Plane != 2 ) continue ;
		// save information 
		rnu_hits  = hit_f.size() ;
		r_chi2_mu = pid_f[k]->Chi2Muon() ;
		r_chi2_pi = pid_f[k]->Chi2Pion() ;
		r_chi2_p  = pid_f[k]->Chi2Proton() ;
		r_PIDA    = pid_f[k]->PIDA();
		r_missenergy = pid_f[k]->MissingE();
		r_KineticEnergy = cal_f[n]->KineticEnergy();
		dEdx = cal_f[n]->dEdx();
		r_Range = cal_f[n]->Range();
	
	      }

	      if( event_id == 35 ) {
		for ( unsigned int t_hit = 0 ; t_hit <  track_f[j]->LastValidPoint()+1; ++t_hit ){
		
		  tr_x =  track_f[j]->TrajectoryPoint( t_hit ).position.X();
		  tr_y =  track_f[j]->TrajectoryPoint( t_hit ).position.Y();
		  tr_z =  track_f[j]->TrajectoryPoint( t_hit ).position.Z();
		  // also need time
		  // dEdx
		  recoTrackInfo_tree -> Fill();
		}
	      }
 
	    }// close pip loop
	  }// close loop reco track
	}  // end if
    } // for
  } // if

  // FILL TREES
  event_tree      -> Fill();
  mcparticle_tree -> Fill();
  recotrack_tree  -> Fill();
  
}


// Functions

void TrackID::MyAnalysis::SaveMCTrack( simb::MCTrajectory const & mc_track, std::string const & path ) const {

  // This should be a function called Print : need to understand what it is printing exactly. Obviously a track, but is it the primary? 
  
  TH3D *h_track = new TH3D("h_track", " Particle Track ", mc_track.size(), mc_track.X(0), mc_track.X(mc_track.size()), mc_track.size(), mc_track.Y(0), mc_track.Y(mc_track.size()), mc_track.size(), mc_track.Z(0), mc_track.Z(mc_track.size()));
  
  for( unsigned int i = 0; i < mc_track.size(); ++i ){
    h_track-> Fill(mc_track.X(i), mc_track.Y(i), mc_track.Z(i), mc_track.E(i));
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
  c->SaveAs((path+".root").c_str());
  c->Clear();
}

void TrackID::MyAnalysis::SaveRecoTrack( art::Ptr<recob::Track> const & reco_track, std::string const & path ) const {
  
  TH3D *h_recotrack = new TH3D("h_recotrack", "Reco Particle Track ", int(reco_track->CountValidPoints()/10), reco_track->Vertex().X(), reco_track->End().X(), int(reco_track->CountValidPoints()/10), reco_track->Vertex().Y(), reco_track->End().Y(), int(reco_track->CountValidPoints()/10), reco_track->Vertex().Z(), reco_track->End().Z() );
  
  for ( unsigned int t_hit = 0 ; t_hit < reco_track->LastValidPoint()+1; ++t_hit ){
    h_recotrack -> Fill( reco_track->TrajectoryPoint( t_hit ).position.X(), reco_track->TrajectoryPoint( t_hit ).position.Y(), reco_track->TrajectoryPoint( t_hit ).position.Z(), reco_track->MomentumAtPoint( t_hit ));
  }
  
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
	    
  TCanvas *c2 = new TCanvas();
  h_recotrack->SetLineColor(2);
  h_recotrack->GetXaxis()->SetTitle("X");
  h_recotrack->GetYaxis()->SetTitle("Y");
  h_recotrack->GetXaxis()->SetTitle("Z");
  h_recotrack->Draw("hist");
  h_recotrack->Draw("BOX2Z");
  c2->SaveAs( (path+".root").c_str() );
  c2->Clear();
}

void TrackID::MyAnalysis::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of required member function here.
  // Here you add an external fcl file to change configuration
}

void TrackID::MyAnalysis::beginJob( )
{
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
  fDaughter_mu = 0 ;
  fDaughter_pi = 0 ;
  fDaughter_e = 0 ;
  fDaughter_p = 0 ;
  fDaughter_n = 0 ;
  fDaughter_photon = 0 ;
  fDaughter_other = 0 ;
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
  r_n_daughters = 0 ;
  r_photon_daughters = 0 ;
  r_other_daughters = 0 ;
  rVertex_x = -999. ;
  rVertex_y = -999. ;
  rVertex_z = -999. ;
  rEnd_x = -999. ;
  rEnd_y = -999. ;
  rEnd_z = -999. ;
  rLength = -999. ;
  rMomentum = -999. ;
  rnu_hits  = -999 ;
  r_chi2_mu = -999. ;
  r_chi2_pi = -999. ;
  r_chi2_p  = -999. ;
  r_PIDA    = -999. ;
  r_missenergy = -999. ;
  r_KineticEnergy = -999. ;
  r_Range = -999. ;

  //Track example information to test future Track class. This is the information needed for the Track class. 
  // Saving it into a root file to work offline...
  tr_x = -999. ;
  tr_y = -999. ;
  tr_z = -999. ;
  tr_dEdx = -999. ;
  tr_dQdx = -999. ;


  // Declare trees and branches
  event_tree      = new TTree( "event_tree",           "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree( "mcparticle_tree",      "MC tree:    True Particle track information");
  recotrack_tree  = new TTree( "recoparticle_tree",    "Reco tree: reconstructed information of the tracks, hit level included");

  recoTrackInfo_tree  = new TTree( "recoTrack_tree",    "Reco Track Info tree: Contains information needed for the fitter offline ");

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
  mcparticle_tree -> Branch( "Daughter_mu",            &fDaughter_mu,        "Daughter_mu/I");
  mcparticle_tree -> Branch( "Daughter_pi",            &fDaughter_pi,        "Daughter_pi/I");
  mcparticle_tree -> Branch( "Daughter_e",             &fDaughter_e,         "Daughter_e/I");
  mcparticle_tree -> Branch( "Daughter_p",             &fDaughter_p,         "Daughter_p/I");
  mcparticle_tree -> Branch( "Daughter_n",             &fDaughter_n,         "Daughter_n/I");
  mcparticle_tree -> Branch( "Daughter_photon",        &fDaughter_photon,    "Daughter_photon/I");
  mcparticle_tree -> Branch( "Daughter_other",         &fDaughter_other,     "Daughter_other/I");
  mcparticle_tree -> Branch( "MC_Length",              &fMCLength,           "Length/D");

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
  recotrack_tree  -> Branch( "n_daughters",       &r_n_daughters,     "n_daughters/I");
  recotrack_tree  -> Branch( "photon_daughters",  &r_photon_daughters,"photon_daughters/I");
  recotrack_tree  -> Branch( "other_daughters",   &r_other_daughters, "others_daughters/I");
  recotrack_tree  -> Branch( "Vertex_x",          &rVertex_x,         "rVertex_x/D");
  recotrack_tree  -> Branch( "Vertex_y",          &rVertex_y,         "rVertex_y/D");
  recotrack_tree  -> Branch( "Vertex_z",          &rVertex_z,         "rVertex_z/D");
  recotrack_tree  -> Branch( "End_x",             &rEnd_x,            "rEnd_x/D");
  recotrack_tree  -> Branch( "End_y",             &rEnd_y,            "rEnd_y/D");
  recotrack_tree  -> Branch( "End_z",             &rEnd_z,            "rEnd_z/D");
  recotrack_tree  -> Branch( "Length",            &rLength,           "rLength/F");
  recotrack_tree  -> Branch( "Momentum",          &rMomentum,         "rMomentum/F");
  recotrack_tree  -> Branch( "nu_hits",           &rnu_hits,          "rnu_hits/I");
  recotrack_tree  -> Branch( "r_chi2_mu",         &r_chi2_mu,         "r_chi2_mu/D");
  recotrack_tree  -> Branch( "r_chi2_pi",         &r_chi2_pi,         "r_chi2_pi/D");
  recotrack_tree  -> Branch( "r_chi2_p",          &r_chi2_p,          "r_chi2_p/D");
  recotrack_tree  -> Branch( "r_PIDA",            &r_PIDA,            "r_PIDA/D");
  recotrack_tree  -> Branch( "r_missing_energy",  &r_missenergy,      "r_missenergy/D");
  recotrack_tree  -> Branch( "r_KineticEnergy",   &r_KineticEnergy,   "r_KineticEnergy/D");
  recotrack_tree  -> Branch( "r_Range",           &r_Range,           "r_Range/D");
  

  recoTrackInfo_tree -> Branch( "tr_x",           &tr_x,              "tr_x/D");
  recoTrackInfo_tree -> Branch( "tr_y",           &tr_y,              "tr_y/D");
  recoTrackInfo_tree -> Branch( "tr_z",           &tr_z,              "tr_z/D");
  recoTrackInfo_tree -> Branch( "tr_dEdx",        &tr_dEdx,           "tr_dEdx/D");
  recoTrackInfo_tree -> Branch( "tr_dQdx",        &tr_dQdx,           "tr_dQdx/D");


  // Set directories
  event_tree->SetDirectory(0);
  mcparticle_tree->SetDirectory(0);
  recotrack_tree->SetDirectory(0);
  recoTrackInfo_tree -> SetDirectory( 0 );
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

  TFile fileTrack("output_track_ev35_mu.root", "RECREATE");
  recoTrackInfo_tree -> Write() ;
  fileTrack.Write();
  fileTrack.Close();

  delete event_tree ; 
  delete mcparticle_tree ; 
  delete recotrack_tree ; 
  delete recoTrackInfo_tree;
}

DEFINE_ART_MODULE(TrackID::MyAnalysis)
