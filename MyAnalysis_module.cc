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
#include "lardataobj/RecoBase/SpacePoint.h" 
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

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
  void StoreInformation( art::Event const & e, art::Handle< std::vector< recob::Track > > & trackHandle, art::Handle< std::vector< recob::Shower > > & showerHandle, art::FindManyP< recob::Track > & findTracks , int & part_id_f ) ;

private:

  // Declare member data here.

  // Labels
  std::string fTruthLabel, m_particleLabel, m_hitfinderLabel, m_recotrackLabel, m_recoshowerLabel, m_recoPIDLabel, m_recoCaloLabel;

  // Tree members
  TTree * event_tree, * mcparticle_tree, * recotrack_tree ; 
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
  recob::TrackTrajectory primary_trajectory ;
  int rnu_hits, rdQdx_size ;

  double r_chi2_mu[10], r_chi2_pi[10], r_chi2_p[10], r_PIDA[10] ;
  double r_missenergy[10], r_KineticEnergy[10], r_Range[10] ;
  float rLength ;
  float r_dQdx[100000]; // r_track_Q[100000];
  double r_track_x[100000], r_track_y[100000], r_track_z[100000], r_track_dQdx[100000];
  std::vector< float > r_dQdx_ID, r_dQdx_total ; // to save it ordered!
  std::vector< std::vector< float >  > r_dQdx_ID_all ; 
  lar_pandora::PFParticleMap particleMap ; 
  
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
  m_recoshowerLabel = p.get<std::string>("RecoShowerLabel","emshower");
  m_recoCaloLabel = p.get<std::string>("RecoCaloLabel","pandoraCalo");
  m_recoPIDLabel = p.get<std::string>("RecoPIDLabel","pandoraPid");

  this->reconfigure(p);
}

void TrackID::MyAnalysis::analyze(art::Event const & e)
{
  event_id = e.id().event();
  std::stringstream t_path, r_path ;
  t_path << "Histograms/eid_"<<event_id<<"_truth_track" ;
  r_path << "Histograms/eid_"<<event_id<<"_reco_track" ;
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
	  else if ( trueParticle.PdgCode() == 11   ) { ++fDaughter_e ;  } 
	  else if ( trueParticle.PdgCode() == 2212 ) { ++fDaughter_p ;  }
	  else if ( trueParticle.PdgCode() == 2112 ) { ++fDaughter_n ;  }
	  else if ( trueParticle.PdgCode() == 22  )  { ++fDaughter_photon ; }
	  else                                       { ++fDaughter_other ;  }
	}
      }
    }
  }

  /**************************************************************************************************
   *  RECO INFORMATION
   *************************************************************************************************/
  // Get PFParticle Handle
  art::Handle< std::vector< recob::PFParticle > > pfParticleHandle ;
  e.getByLabel(m_particleLabel, pfParticleHandle ) ;

  // Save map for hiearchy info
  lar_pandora::PFParticleVector pfplist; 
  
  //Find the hits associated to the reconstructed PFParticle
  art::Handle< std::vector< recob::Hit > > hitListHandle ;
  e.getByLabel(m_hitfinderLabel, hitListHandle);

  //Find the reco tracks
  art::Handle< std::vector< recob::Track > > trackHandle ;
  e.getByLabel(m_recotrackLabel, trackHandle ) ;

  //Find the reco showers
  art::Handle< std::vector< recob::Shower > > showerHandle ;
  e.getByLabel(m_recoshowerLabel, showerHandle ) ;
  
  // Get track associations with PFParticles from Pandora. Find all possible tracks associated to an event
  art::FindManyP< recob::Track > findTracks( pfParticleHandle, e, m_recotrackLabel );

  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid() ){
    art::fill_ptr_vector(pfplist, pfParticleHandle );
    lar_pandora::LArPandoraHelper::BuildPFParticleMap( pfplist, particleMap );
  }

  rLength = 0 ;
  rnu_hits = 0 ;
  rdQdx_size = 0 ;
  r_dQdx_ID.clear() ; // clean individual one

  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid()){
    
    for( unsigned int i = 0 ; i < pfParticleHandle->size(); ++i ){ // loop over pfParticleHandle to find Primary
      art::Ptr< recob::PFParticle > pfparticle( pfParticleHandle, i ) ; // Point to particle i 

      if( pfparticle->IsPrimary() == 1 ) { //found primary particle and starting point 
	if( pfparticle->NumDaughters() > 0 ) {
	  for( int j = 0 ; j < pfparticle->NumDaughters() ; ++j ){ // looping over daughters to read them in order 
	    int part_id_f = particleMap[ pfparticle->Daughters()[j] ] -> Self() ;
	    StoreInformation( e, trackHandle, showerHandle, findTracks, part_id_f ) ;

	    if( pfparticle->NumDaughters() > 0 ) {
	      for( int j2 = 0 ; j2 < pfparticle->NumDaughters() ; ++j2 ){ // looping over daughters to read them in order 
		
		if( particleMap[ pfparticle->Daughters()[j2] ] -> NumDaughters() > 0 ) { 
		  int id_2daughter = particleMap[ pfparticle->Daughters()[j2] ]->Daughters()[0] ;
		  StoreInformation( e, trackHandle, showerHandle, findTracks, id_2daughter ) ; 
		}
	      }
	    }
	  }
	}
      } // primary 
    } //pfparticleHandle
  } //valid handle
  
  // FILL TREES
  event_tree      -> Fill();
  mcparticle_tree -> Fill();
  recotrack_tree  -> Fill();
  
} // event 

void TrackID::MyAnalysis::StoreInformation( art::Event const & e, art::Handle< std::vector< recob::Track > > & trackHandle, art::Handle< std::vector< recob::Shower > > & showerHandle, art::FindManyP< recob::Track > & findTracks , int & part_id_f ) {
        
      // Save track info
      if ( findTracks.at( part_id_f ).size()!=0 ){
	std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(part_id_f);
	art::FindManyP< recob::Hit > findHits (  trackHandle, e, m_recotrackLabel ) ;
	art::FindManyP< anab::Calorimetry > findCalorimetry ( trackHandle, e, m_recoCaloLabel );
	art::FindManyP< anab::ParticleID > findPID ( trackHandle, e, m_recoPIDLabel );
	std::cout<< " hits - " << rnu_hits << std::endl;

	// Loop over tracks per event
	for( unsigned int n = 0 ; n < track_f.size() ; ++n ){
	  
	  rLength   += track_f[n]->Length() ;
	  // Get track based variables
	  std::vector< art::Ptr<recob::Hit> > hit_f        = findHits.at(track_f[n]->ID()); 
	  std::vector< art::Ptr<anab::Calorimetry> > cal_f = findCalorimetry.at(track_f[n]->ID());
	  std::vector< art::Ptr<anab::ParticleID> > pid_f  = findPID.at(track_f[n]->ID());
	  
	  //Loop over PID associations 
	  for ( unsigned int k = 0 ; k < pid_f.size() ; ++k ){
	    
	    if( !pid_f[k] ) continue ;
	    if( !pid_f[k]->PlaneID().isValid) continue ;
	    if( pid_f[k]->PlaneID().Plane != 2 ) continue ; // only look at collection plane for dEdx information
	    
	    //Loop over calo information also in collection plane
	    for ( unsigned int m = 0 ; m < cal_f.size() ; ++m ) {
	      if( !cal_f[m] ) continue ;
	      if( !cal_f[m]->PlaneID().isValid) continue ;
	      if( cal_f[m]->PlaneID().Plane == 2 ) {	    
		// save information -> Before instead of part_id_f it was j : check if right . 
		r_chi2_mu[part_id_f] = pid_f[k]->Chi2Muon() ;
		r_chi2_pi[part_id_f] = pid_f[k]->Chi2Pion() ;
		r_chi2_p[part_id_f]  = pid_f[k]->Chi2Proton() ;
		r_PIDA[part_id_f]    = pid_f[k]->PIDA();
		r_missenergy[part_id_f] = pid_f[k]->MissingE();
		r_KineticEnergy[part_id_f] = cal_f[m]->KineticEnergy();
		
		for( unsigned int l = 0 ; l < track_f[n]->LastValidPoint()+1 ; ++l ) {
		  r_track_x[l+rnu_hits] = track_f[n]->TrajectoryPoint( l ).position.X();
		  r_track_y[l+rnu_hits] = track_f[n]->TrajectoryPoint( l ).position.Y();
		  r_track_z[l+rnu_hits] = track_f[n]->TrajectoryPoint( l ).position.Z();
		}
		rnu_hits   += track_f[n]->LastValidPoint() + 1 ; // ?: +1 // valid hits
	      }// just collection plane 
	      
	      // calo information is stored in all planes:
	      
	      for( unsigned int l = 0 ; l < (cal_f[m]->dQdx()).size() ; ++l ) r_dQdx[l+rdQdx_size] = cal_f[m]->dQdx()[l];
	      for( unsigned int l = 0 ; l < (cal_f[m]->XYZ()).size() ; ++l ) {
		for( int t = 0 ; t < rnu_hits ; ++t ){
		  if( cal_f[m]->XYZ()[l].X() == r_track_x[t] && cal_f[m]->XYZ()[l].Y() == r_track_y[t] && cal_f[m]->XYZ()[l].Z() == r_track_z[t]) {
		    r_track_dQdx[t] = cal_f[m]->dQdx()[l] ; 
		  }
		}
	      }
	      
	      r_Range[part_id_f] = cal_f[m]->Range(); // check 
	      rdQdx_size += (cal_f[m]->dQdx()).size();
	    } //close calo
	  } //close pid
	} //close track  	  
      } else if( showerHandle.isValid() && showerHandle->size() ) { // if no track look in showers 
	art::FindManyP< recob::Hit > findHitShower( showerHandle, e, m_recoshowerLabel ) ;
	art::FindManyP< recob::SpacePoint > findSpacePoint( showerHandle, e, m_recoshowerLabel ) ;

	for( unsigned int y = 0 ; y < showerHandle->size() ; ++y ) {
	  art::Ptr< recob::Shower > shower_f( showerHandle, y ) ;
	  std::vector< art::Ptr<recob::Hit> > hit_sh_f = findHitShower.at(y) ; 
	  std::vector< art::Ptr<recob::SpacePoint> > spacepoint_f = findSpacePoint.at(y) ;
	  for( unsigned int l = 0 ; l < spacepoint_f.size(); ++l ) {
	      r_track_x[l+rnu_hits] = spacepoint_f[l]->XYZ()[0] ;
	      r_track_y[l+rnu_hits] = spacepoint_f[l]->XYZ()[1] ;
	      r_track_z[l+rnu_hits] = spacepoint_f[l]->XYZ()[2] ;
	    }
	    rnu_hits   += spacepoint_f.size() + 1 ; // ?: +1 // valid hits
	}	
      } // track vs shower
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
  rLength = -999. ;
  rnu_hits  = 0 ;
  
  // Declare trees and branches
  event_tree      = new TTree( "event_tree",           "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree( "mcparticle_tree",      "MC tree:    True Particle track information");
  recotrack_tree  = new TTree( "recoparticle_tree",    "Reco tree: reconstructed information of the tracks, hit level included");

  event_tree      -> Branch( "event_id",          &event_id, "event_id/I");

  /**
     MC PARTICLE TREE BRANCHES :
   */
  mcparticle_tree -> Branch( "event_id",                &event_id,            "event_id/I");
  mcparticle_tree -> Branch( "fTrack_ID",               &fTrack_ID,           "fTrack_ID/I");
  mcparticle_tree -> Branch( "ftrueEnergy",             &fTrueParticleEnergy, "fTrueParticleEnergy/F");
  mcparticle_tree -> Branch( "fPDG_Code",               &fPDG_Code,           "fPDG_Code/I");
  mcparticle_tree -> Branch( "fMass",                   &fMass,               "fMass/F");
  mcparticle_tree -> Branch( "fpx",                     &fpx,                 "fpx/F");
  mcparticle_tree -> Branch( "fpy",                     &fpy,                 "fpy/F");
  mcparticle_tree -> Branch( "fpz",                     &fpz,                 "fpz/F");
  mcparticle_tree -> Branch( "fpt",                     &fpt,                 "fpt/F");
  mcparticle_tree -> Branch( "fp",                      &fp,                  "fp/F");
  mcparticle_tree -> Branch( "fNumDaughters",           &fNumDaughters,       "fNumDaughters/I");
  mcparticle_tree -> Branch( "fDaughter_mu",            &fDaughter_mu,        "fDaughter_mu/I");
  mcparticle_tree -> Branch( "fDaughter_pi",            &fDaughter_pi,        "fDaughter_pi/I");
  mcparticle_tree -> Branch( "fDaughter_e",             &fDaughter_e,         "fDaughter_e/I");
  mcparticle_tree -> Branch( "fDaughter_p",             &fDaughter_p,         "fDaughter_p/I");
  mcparticle_tree -> Branch( "fDaughter_n",             &fDaughter_n,         "fDaughter_n/I");
  mcparticle_tree -> Branch( "fDaughter_photon",        &fDaughter_photon,    "fDaughter_photon/I");
  mcparticle_tree -> Branch( "fDaughter_other",         &fDaughter_other,     "fDaughter_other/I");
  mcparticle_tree -> Branch( "fMC_Length",              &fMCLength,           "fMCLength/D");

 /**
     RECONSTRUCTED PARTICLE TREE BRANCHES :
   */
  recotrack_tree  -> Branch( "event_id",            &event_id,          "event_id/I");
  recotrack_tree  -> Branch( "r_pdg_primary",       &r_pdg_primary,     "r_pdg_primary/I");
  recotrack_tree  -> Branch( "r_nu_daughters",      &r_nu_daughters,    "r_nu_daughters/I");
  recotrack_tree  -> Branch( "rLength",             &rLength,           "rLength/F");
  recotrack_tree  -> Branch( "rnu_hits",            &rnu_hits,          "rnu_hits/I");
  recotrack_tree  -> Branch( "r_chi2_mu",           &r_chi2_mu,         ("r_chi2_mu[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_chi2_pi",           &r_chi2_pi,         ("r_chi2_pi[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_chi2_p",            &r_chi2_p,          ("r_chi2_p[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_PIDA",              &r_PIDA,            ("r_PIDA[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_missing_energy",    &r_missenergy,      ("r_missenergy[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_KineticEnergy",     &r_KineticEnergy,   ("r_KineticEnergy[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_Range",             &r_Range,           ("r_Range[" + std::to_string(10)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_dQdx",              &r_dQdx,            ("r_dQdx[" + std::to_string(100000)+"]/F").c_str());
  recotrack_tree  -> Branch( "r_track_x",           &r_track_x,         ("r_track_x[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_track_y",           &r_track_y,         ("r_track_y[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_track_z",           &r_track_z,         ("r_track_z[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_track_dQdx",        &r_track_dQdx,      ("r_track_dQdx[" + std::to_string(100000)+"]/D").c_str());

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
