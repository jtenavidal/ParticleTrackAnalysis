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

private:

  // Declare member data here.

  // Labels
  std::string fTruthLabel, m_particleLabel, m_hitfinderLabel, m_recotrackLabel, m_recoPIDLabel, m_recoCaloLabel;

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
  int r_mu_daughters, r_pi_daughters, r_e_daughters, r_p_daughters, r_n_daughters, r_photon_daughters, r_other_daughters;
  recob::TrackTrajectory primary_trajectory ;
  double rVertex_x, rVertex_y, rVertex_z, rEnd_x, rEnd_y, rEnd_z;
  int rnu_hits, rnu_hits_size, rdQdx_size ;

  double r_chi2_mu, r_chi2_pi, r_chi2_p, r_PIDA, r_missenergy, r_KineticEnergy, r_Range ;
  float rLength ;
  float r_dQdx[100000]; // r_track_Q[100000];
  double r_track_x[100000], r_track_y[100000], r_track_z[100000];
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
	  else if ( trueParticle.PdgCode() == 11   ) { ++fDaughter_e ; } 
	  else if ( trueParticle.PdgCode() == 2212 ) { ++fDaughter_p ; }
	  else if ( trueParticle.PdgCode() == 2112 ) { ++fDaughter_n ; }
	  else if ( trueParticle.PdgCode() == 22  )  { ++fDaughter_photon ; }
	  else                                       { ++fDaughter_other ; }
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
  
  // Get track associations with PFParticles from Pandora. Find all possible tracks associated to an event
  art::FindManyP< recob::Track > findTracks( pfParticleHandle, e, m_recotrackLabel );


  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid() ){
    art::fill_ptr_vector(pfplist, pfParticleHandle );
    lar_pandora::LArPandoraHelper::BuildPFParticleMap( pfplist, particleMap );
  }

  r_mu_daughters = 0;
  r_pi_daughters = 0;
  r_e_daughters = 0;
  r_p_daughters = 0;
  r_n_daughters = 0; // expecting non reco 
  r_photon_daughters = 0; // should not find showers
  r_other_daughters = 0;
  rnu_hits = 0 ;
  rnu_hits_size = 0 ;
  rdQdx_size = 0 ;
  r_dQdx_ID.clear() ; // clean individual one

  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid()){

    for( unsigned int i = 0 ; i < pfParticleHandle->size(); ++i ){ // loop over pfParticleHandle to find Primary
      art::Ptr< recob::PFParticle > pfparticle( pfParticleHandle, i ) ; // Point to particle i 

      if( pfparticle->IsPrimary() == 1 ){
	if( pfparticle->NumDaughters() > 0 ) {
	  // Primary particle is now the first track :
	  r_pdg_primary = particleMap[ pfparticle->Daughters()[0] ] -> PdgCode() ;
	  r_nu_daughters = pfparticle->NumDaughters() - 1 ; // substracting the primary from daughters list

	  for( int j = 1 ; j < pfparticle->NumDaughters() ; ++j ){ // looping over daughters 

	    if        ( particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() == 13   ) { ++ r_mu_daughters ;
	    } else if ( particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() == 211  ) { ++ r_pi_daughters ;
	    } else if ( particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() == 11   ) { ++ r_e_daughters ;
	    } else if ( particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() == 2212 ) { ++ r_p_daughters ;
	    } else if ( particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() == 2112 ) { ++ r_n_daughters ;
	    } else if ( particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() == 22   ) { ++ r_photon_daughters ;
	    } else                                                                                { ++ r_other_daughters ; }

	  }

	  for( int j = 0 ; j < pfparticle->NumDaughters() ; ++j ){ // looping over daughters 
	    int part_id_f = particleMap[ pfparticle->Daughters()[j] ] -> Self() ;
	    std::cout<< "Even ID = " << event_id<< "    particle_id_f= " << part_id_f << std::endl;
	    if ( findTracks.at( part_id_f ).size()!=0 ){
	      std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(part_id_f);
	      art::FindManyP< recob::Hit > findHits (  trackHandle, e, m_recotrackLabel ) ;
	      art::FindManyP< anab::Calorimetry > findCalorimetry ( trackHandle, e, m_recoCaloLabel );
	      art::FindManyP< anab::ParticleID > findPID ( trackHandle, e, m_recoPIDLabel );

	      // Loop over tracks per event
	      for( unsigned int n = 0 ; n < track_f.size() ; ++n ){
	  
		rVertex_x = track_f[n]->Vertex( ).X() ;
		rVertex_y = track_f[n]->Vertex( ).Y() ;
		rVertex_z = track_f[n]->Vertex( ).Z() ;
		rEnd_x    = track_f[n]->End( ).X() ;
		rEnd_y    = track_f[n]->End( ).Y() ;
		rEnd_z    = track_f[n]->End( ).Z() ;
		rLength   = track_f[n]->Length() ;
		
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
		    if( cal_f[m]->PlaneID().Plane != 2 ) continue ;
		    
		    // save information 
		    r_chi2_mu = pid_f[k]->Chi2Muon() ;
		    r_chi2_pi = pid_f[k]->Chi2Pion() ;
		    r_chi2_p  = pid_f[k]->Chi2Proton() ;
		    r_PIDA    = pid_f[k]->PIDA();
		    r_missenergy = pid_f[k]->MissingE();
		    r_KineticEnergy = cal_f[m]->KineticEnergy();
		    
		    for( unsigned int l = 0 ; l < (cal_f[m]->dQdx()).size() ; ++l ) r_dQdx[l+rdQdx_size] = cal_f[n]->dQdx()[l];
		    r_Range = cal_f[m]->Range();
		    
		    for( unsigned int l = 0 ; l < track_f[n]->LastValidPoint()+1 ; ++l ) {
		      r_track_x[l+rnu_hits] = track_f[n]->TrajectoryPoint( l ).position.X();
		      r_track_y[l+rnu_hits] = track_f[n]->TrajectoryPoint( l ).position.Y();
		      r_track_z[l+rnu_hits] = track_f[n]->TrajectoryPoint( l ).position.Z();
		    }
		    //		for( unsigned int l = 0 ; l < hit_f.size() ; ++l ) r_track_Q[l+rnu_hits_size] = hit_f[l] -> Integral() ;
		    rnu_hits   += track_f[n]->LastValidPoint() + 1 ; // ?: +1
		    rnu_hits_size += hit_f.size() ;
		    rdQdx_size += (cal_f[m]->dQdx()).size();
		  } //close calo
		} //close pid
	      } //close track  	  
	    } // if find tracks 
	  } //loop j particle daughter
	} //if numdaugh>0
      }//if primary 
    } //pfparticleHandle
  } //valid handle

  // FILL TREES
  event_tree      -> Fill();
  mcparticle_tree -> Fill();
  recotrack_tree  -> Fill();
  
} // event 



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
  rnu_hits  = 0 ;
  rnu_hits_size = 0 ;
  r_chi2_mu = -999. ;
  r_chi2_pi = -999. ;
  r_chi2_p  = -999. ;
  r_PIDA    = -999. ;
  r_missenergy = -999. ;
  r_KineticEnergy = -999. ;
  r_Range = -999. ;
  r_dQdx_ID_all.clear();
  r_dQdx_total.clear();
  
  // Declare trees and branches
  event_tree      = new TTree( "event_tree",           "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree( "mcparticle_tree",      "MC tree:    True Particle track information");
  recotrack_tree  = new TTree( "recoparticle_tree",    "Reco tree: reconstructed information of the tracks, hit level included");

  event_tree      -> Branch( "event_id",          &event_id, "event_id/I");

  /**
     MC PARTICLE TREE BRANCHES :
   */
  mcparticle_tree -> Branch( "event_id",                &event_id,            "event_id/I");
  mcparticle_tree -> Branch( "fTrack_ID",               &fTrack_ID,           "Track_id/I");
  mcparticle_tree -> Branch( "ftrueEnergy",             &fTrueParticleEnergy, "TrueParticleEnergy/F");
  mcparticle_tree -> Branch( "fPDG_Code",               &fPDG_Code,           "PDG_Code/I");
  mcparticle_tree -> Branch( "fMass",                   &fMass,               "Mass/F");
  mcparticle_tree -> Branch( "fpx",                     &fpx,                 "px/F");
  mcparticle_tree -> Branch( "fpy",                     &fpy,                 "py/F");
  mcparticle_tree -> Branch( "fpz",                     &fpz,                 "pz/F");
  mcparticle_tree -> Branch( "fpt",                     &fpt,                 "pt/F");
  mcparticle_tree -> Branch( "fp",                      &fp,                  "p/F");
  mcparticle_tree -> Branch( "fNum_Daughters",          &fNumDaughters,       "num_d/I");
  mcparticle_tree -> Branch( "fDaughter_mu",            &fDaughter_mu,        "Daughter_mu/I");
  mcparticle_tree -> Branch( "fDaughter_pi",            &fDaughter_pi,        "Daughter_pi/I");
  mcparticle_tree -> Branch( "fDaughter_e",             &fDaughter_e,         "Daughter_e/I");
  mcparticle_tree -> Branch( "fDaughter_p",             &fDaughter_p,         "Daughter_p/I");
  mcparticle_tree -> Branch( "fDaughter_n",             &fDaughter_n,         "Daughter_n/I");
  mcparticle_tree -> Branch( "fDaughter_photon",        &fDaughter_photon,    "Daughter_photon/I");
  mcparticle_tree -> Branch( "fDaughter_other",         &fDaughter_other,     "Daughter_other/I");
  mcparticle_tree -> Branch( "fMC_Length",              &fMCLength,           "Length/D");

 /**
     RECONSTRUCTED PARTICLE TREE BRANCHES :
   */
  recotrack_tree  -> Branch( "event_id",            &event_id,          "event_id/I");
  recotrack_tree  -> Branch( "r_pdg_primary",       &r_pdg_primary,     "pdg_primary/I");
  recotrack_tree  -> Branch( "r_nu_daughters",      &r_nu_daughters,    "nu_daughters/I");
  recotrack_tree  -> Branch( "r_mu_daughters",      &r_mu_daughters,    "mu_daughters/I");
  recotrack_tree  -> Branch( "r_pi_daughters",      &r_pi_daughters,    "pi_daughters/I");
  recotrack_tree  -> Branch( "r_e_daughters",       &r_e_daughters,     "e_daughters/I");
  recotrack_tree  -> Branch( "r_p_daughters",       &r_p_daughters,     "p_daughters/I");
  recotrack_tree  -> Branch( "r_n_daughters",       &r_n_daughters,     "n_daughters/I");
  recotrack_tree  -> Branch( "r_photon_daughters",  &r_photon_daughters,"photon_daughters/I");
  recotrack_tree  -> Branch( "r_other_daughters",   &r_other_daughters, "others_daughters/I");
  recotrack_tree  -> Branch( "rVertex_x",           &rVertex_x,         "rVertex_x/D");
  recotrack_tree  -> Branch( "rVertex_y",           &rVertex_y,         "rVertex_y/D");
  recotrack_tree  -> Branch( "rVertex_z",           &rVertex_z,         "rVertex_z/D");
  recotrack_tree  -> Branch( "rEnd_x",              &rEnd_x,            "rEnd_x/D");
  recotrack_tree  -> Branch( "rEnd_y",              &rEnd_y,            "rEnd_y/D");
  recotrack_tree  -> Branch( "rEnd_z",              &rEnd_z,            "rEnd_z/D");
  recotrack_tree  -> Branch( "rLength",             &rLength,           "rLength/F");
  recotrack_tree  -> Branch( "rnu_hits",            &rnu_hits,          "rnu_hits/I");
  recotrack_tree  -> Branch( "rnu_hits_size",       &rnu_hits_size,     "rnu_hits_size/I");
  recotrack_tree  -> Branch( "r_chi2_mu",           &r_chi2_mu,         "r_chi2_mu/D");
  recotrack_tree  -> Branch( "r_chi2_pi",           &r_chi2_pi,         "r_chi2_pi/D");
  recotrack_tree  -> Branch( "r_chi2_p",            &r_chi2_p,          "r_chi2_p/D");
  recotrack_tree  -> Branch( "r_PIDA",              &r_PIDA,            "r_PIDA/D");
  recotrack_tree  -> Branch( "r_missing_energy",    &r_missenergy,      "r_missenergy/D");
  recotrack_tree  -> Branch( "r_KineticEnergy",     &r_KineticEnergy,   "r_KineticEnergy/D");
  recotrack_tree  -> Branch( "r_Range",             &r_Range,           "r_Range/D");
  recotrack_tree  -> Branch( "rdQdx_size",          &rdQdx_size,        "rdQdx_size/I");
  recotrack_tree  -> Branch( "r_dQdx",              &r_dQdx,            ("r_dQdx[" + std::to_string(100000)+"]/F").c_str());
  recotrack_tree  -> Branch( "r_track_x",           &r_track_x,         ("r_track_x[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_track_y",           &r_track_y,         ("r_track_y[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree  -> Branch( "r_track_z",           &r_track_z,         ("r_track_z[" + std::to_string(100000)+"]/D").c_str());
  //  recotrack_tree  -> Branch( "r_track_Q",           &r_track_Q,         ("r_track_Q[" + std::to_string(100000)+"]/F").c_str());


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
