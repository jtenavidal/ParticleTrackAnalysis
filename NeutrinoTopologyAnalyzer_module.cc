
////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoTopologyAnalyzer
// Plugin Type: analyzer (art v3_00_00)
// File:        NeutrinoTopologyAnalyzer_module.cc
//
// Generated at Mon Mar 25 06:08:16 2019 by Julia Tena Vidal using cetskelgen
// from cetlib version v3_04_00.
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
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
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
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/RecoUtils/ShowerUtils.h"

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
#include <map>

namespace test {
  class NeutrinoTopologyAnalyzer;
}


class test::NeutrinoTopologyAnalyzer : public art::EDAnalyzer {
public:
  explicit NeutrinoTopologyAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoTopologyAnalyzer(NeutrinoTopologyAnalyzer const&) = delete;
  NeutrinoTopologyAnalyzer(NeutrinoTopologyAnalyzer&&) = delete;
  NeutrinoTopologyAnalyzer& operator=(NeutrinoTopologyAnalyzer const&) = delete;
  NeutrinoTopologyAnalyzer& operator=(NeutrinoTopologyAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // -> Added functions ( * )

  void reconfigure(fhicl::ParameterSet const & p);
  //  void clearVariables() ;
  void beginJob() override;
  void endJob() override;
  void clearVariables() ; 
private:
  // Declare member data here ;
  // Labels
  std::string TruthLabel, G4Label, ParticleLabel, HitFinderLabel, RecoTrackLabel, RecoShowerLabel, RecoPIDLabel, RecoCaloLabel ; 

  // Detector information
  float DetectorHalfLengthX, DetectorHalfLengthY, DetectorHalfLengthZ, CoordinateOffSetX, CoordinateOffSetY, CoordinateOffSetZ, SelectedBorderX, SelectedBorderY, SelectedBorderZ ;

  // Tree members
  TTree * mcparticle_tree, * recoevent_tree ; 
  
  // Event Information
  int event_id ;

  // MC Event Information 
  // neutrino information
  int Tnu_PDG, T_interaction ;
  double t_vertex[3], t_momentum[3], t_vertex_energy ;
  bool is_cc ; 

  // Truth information ( = > T ) 
  int TPDG_Code, TNumDaughters, TDaughter_mu, TDaughter_pi, TDaughter_e, TDaughter_p, TDaughter_n, TDaughter_photon, TDaughter_other ;
  float TrueParticleEnergy, TMass;
  float Tpx, Tpy, Tpz, Tpt, Tp; 
  double TMCLength, TTrack_vertex_x, TTrack_vertex_y, TTrack_vertex_z, TTrack_vertex_t, TTrack_end_x, TTrack_end_y, TTrack_end_z, TTrack_end_t ;
  // vectorization
  int TPDG_Primary[10000], TNumDaughPrimary[10000];


  // Reco information
  lar_pandora::PFParticleMap particleMap ;   

  bool is_reconstructed, primary_vcontained, primary_econtained ;
  int r_pdg_primary[10000] ;
  int rnu_hits[10000] ;

  //    => neutrino vertex information
  bool nu_reconstructed , nu_vertex_contained ;
  int numb_nu_reco ;
  double nu_reco_vertex[3] ;
  double error_vertex_reco ;

  //  double r_chi2_mu[10], r_chi2_pi[10], r_chi2_p[10], r_PIDA[10] ;
  //double r_missenergy[10], r_KineticEnergy[10], r_Range[10] ;
  //float rLength ;
  //double r_track_x[100000], r_track_y[100000], r_track_z[100000], r_track_dEdx[100000];
  //std::vector< float > r_dEdx_ID, r_dEdx_total ; 
  //std::vector< std::vector< float >  > r_dEdx_ID_all ; 
  /*
  // -> Breakdown particles in event from Pandora
  int tr_id_energy, tr_id_charge, tr_id_hits;
  int pfps_truePDG[1000] ;
  int event_vcontained[1000], event_econtained[1000] ; 
  int pfps_hits[1000] , pfps_type[1000] ;
  float pfps_length[1000] ; 
  double pfps_dir_start_x[1000], pfps_dir_start_y[1000], pfps_dir_start_z[1000], pfps_dir_end_x[1000] , pfps_dir_end_y[1000] , pfps_dir_end_z[1000],pfps_start_x[1000] , pfps_start_y[1000] , pfps_start_z[1000] , pfps_end_x[1000], pfps_end_y[1000] , pfps_end_z[1000] ;
  */
};


test::NeutrinoTopologyAnalyzer::NeutrinoTopologyAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  TruthLabel = p.get<std::string>("TruthLabel");
  G4Label = p.get<std::string>("G4Label");
  ParticleLabel = p.get<std::string>("PFParticleModule","pandora");
  HitFinderLabel = p.get<std::string>("HitFinderModule","linecluster");
  RecoTrackLabel = p.get<std::string>("RecoTrackLabel","pandoraTrack");
  RecoShowerLabel = p.get<std::string>("RecoShowerLabel","emshower");
  RecoCaloLabel = p.get<std::string>("RecoCaloLabel","pandoraCalo");
  RecoPIDLabel = p.get<std::string>("RecoPIDLabel","pandoraPid");
  
  this->reconfigure(p);
}

void test::NeutrinoTopologyAnalyzer::analyze(art::Event const& e)
{
  clearVariables();
  // Implementation of required member function here.
  event_id = e.id().event();
  std::cout<< " Event ID = " << event_id <<std::endl;
  if( !e.isRealData()){
    /**************************************************************************************************
     *  MC INFORMATION
     *************************************************************************************************/
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(G4Label);
    art::Handle< std::vector< simb::MCTruth > > mct_handle;
    e.getByLabel(TruthLabel, mct_handle);
    // Get the GEANT information of the particles 
    art::FindManyP< simb::MCParticle > fmcp( mct_handle , e, G4Label );
    art::FindManyP< simb::GTruth > fmgt( mct_handle, e, TruthLabel );
    art::Ptr< simb::MCTruth > mct(mct_handle, 0);
    std::vector< art::Ptr<simb::MCParticle> > mcp_assn = fmcp.at(0);
    std::vector< art::Ptr<simb::GTruth> > mcgt_assn = fmgt.at(0);
    simb::MCNeutrino nu = mct->GetNeutrino();

    Tnu_PDG = nu.Nu().PdgCode() ;
    T_interaction = mcgt_assn[0]->fGint ;
    t_vertex[0] = nu.Nu().Vx();
    t_vertex[1] = nu.Nu().Vy();
    t_vertex[2] = nu.Nu().Vz();
    t_vertex_energy = nu.Nu().E();
    is_cc = nu.CCNC() == simb::curr_type_::kCC;
    
    if(mcParticles.isValid()){
      // Loop over true info
      for( unsigned int t = 0; t < mcParticles->size(); ++t ){
	const simb::MCParticle trueParticle = mcParticles->at(t) ;
	if( trueParticle.PdgCode() >= 1000018038 ) continue ; // Cut on PDG codes which refer to elements (Argon30 and above)
	if(trueParticle.Process() == "primary" ){
	  // updating to neutrino events : need to vectorize it all
	  TPDG_Primary[t] = trueParticle.PdgCode() ;
	  TNumDaughPrimary[t] = trueParticle.NumberDaughters() ;
	  TNumDaughters += 1 ;
	  
	  
	  TPDG_Code = trueParticle.PdgCode() ;
	  TrueParticleEnergy = trueParticle.E() ;
	  TMass = trueParticle.Mass() ;
	  Tpx = trueParticle.Px() ;
	  Tpy = trueParticle.Py() ;
	  Tpz = trueParticle.Pz() ;
	  Tpt = trueParticle.Pt() ;
	  Tp  = trueParticle.P() ;
	  TTrack_vertex_x = trueParticle.Trajectory().X( 0 ) ;
	  TTrack_vertex_y = trueParticle.Trajectory().Y( 0 ) ; 
	  TTrack_vertex_z = trueParticle.Trajectory().Z( 0 ) ; 
	  TTrack_vertex_t = trueParticle.Trajectory().T( 0 ) ;
	  TTrack_end_x = trueParticle.EndX() ;
	  TTrack_end_y = trueParticle.EndY() ;
	  TTrack_end_z = trueParticle.EndZ() ; 
	  TTrack_end_t = trueParticle.EndT() ;
	  TMCLength = trueParticle.Trajectory().TotalLength() ;
	  
	  //primary_vcontained = MCIsContained( trueParticle )[0] ; 
	  //primary_econtained = MCIsContained( trueParticle )[1] ; 
	} else { // secondary particle information : Just storing total number of each type
	  if      ( trueParticle.PdgCode() == 13   ) { ++TDaughter_mu ; }
	  else if ( trueParticle.PdgCode() == 211  ) { ++TDaughter_pi ; } 
	  else if ( trueParticle.PdgCode() == 11   ) { ++TDaughter_e ;  } 
	  else if ( trueParticle.PdgCode() == 2212 ) { ++TDaughter_p ;  }
	  else if ( trueParticle.PdgCode() == 2112 ) { ++TDaughter_n ;  }
	  else if ( trueParticle.PdgCode() == 22   ) { ++TDaughter_photon ; }
	  else                                       { ++TDaughter_other ;  }
	}
      }
    }
  }


  /**************************************************************************************************
   *  RECO INFORMATION
   *************************************************************************************************/
  // Get PFParticle Handle
  art::Handle< std::vector< recob::PFParticle > > pfParticleHandle ;
  e.getByLabel(ParticleLabel, pfParticleHandle ) ;

  // Save map for hiearchy info
  lar_pandora::PFParticleVector pfplist; 
  
  //Find the hits associated to the reconstructed PFParticle
  art::Handle< std::vector< recob::Hit > > hitListHandle ;
  e.getByLabel(HitFinderLabel, hitListHandle);

  //Find the reco tracks
  art::Handle< std::vector< recob::Track > > trackHandle ;
  e.getByLabel(RecoTrackLabel, trackHandle ) ;
  
  // Get track associations with PFParticles from Pandora. Find all possible tracks associated to an event
  art::FindManyP< recob::Track > findTracks( pfParticleHandle, e, RecoTrackLabel );
  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid() ){
    art::fill_ptr_vector(pfplist, pfParticleHandle );
    lar_pandora::LArPandoraHelper::BuildPFParticleMap( pfplist, particleMap );
  }

  //Find the reco showers
  art::Handle< std::vector< recob::Shower > > showerHandle ;
  e.getByLabel(RecoShowerLabel, showerHandle ) ;
  art::FindManyP< recob::Shower > findShowers( pfParticleHandle, e, ParticleLabel );

  if( pfParticleHandle.isValid() && pfParticleHandle->size() && hitListHandle.isValid() && trackHandle.isValid() ){
    is_reconstructed = true ; 
    for( unsigned int i = 0 ; i < pfParticleHandle->size(); ++i ){ // loop over pfParticleHandle to find Primary
      art::Ptr< recob::PFParticle > pfparticle( pfParticleHandle, i ) ; // Point to particle i 
      
      if( pfparticle->IsPrimary() == 1 ) { //found primary particle => NEUTRINO
	// Get vertex association
	numb_nu_reco += 1 ;
	if( !nu_reconstructed ) { // only looking for one neutrino event -> Check why there can be more than one  
	  nu_reconstructed = true ; 
	  
	  art::FindManyP< recob::Vertex  > fvtx( pfParticleHandle, e, ParticleLabel );
	  std::vector< art::Ptr<recob::Vertex> > vtx_assn = fvtx.at(pfparticle->Self()); 
	  
	  nu_reco_vertex[0] = vtx_assn[0]->position().X() ;
	  nu_reco_vertex[1] = vtx_assn[0]->position().Y() ;
	  nu_reco_vertex[2] = vtx_assn[0]->position().Z() ;

	  // Check if neutrino is contained : 
	  if( ( vtx_assn[0]->position().X() > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
	      || ( vtx_assn[0]->position().X() < (-CoordinateOffSetX + SelectedBorderX)) 
	      || ( vtx_assn[0]->position().Y() > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
	      || ( vtx_assn[0]->position().Y() < (-CoordinateOffSetY + SelectedBorderY)) 
	      || ( vtx_assn[0]->position().Z() > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
	      || ( vtx_assn[0]->position().Z() < (-CoordinateOffSetZ + SelectedBorderZ))) nu_vertex_contained = false ;
	  

	}
      }// end if primary
    }// end loop over particles
  }

  for ( int q = 0 ; q < 3 ; ++q ){
    error_vertex_reco += (t_vertex[q] - nu_reco_vertex[q]) * (t_vertex[q] - nu_reco_vertex[q]) ;
  }
  error_vertex_reco = sqrt( error_vertex_reco ) ;

  mcparticle_tree -> Fill();
  recoevent_tree -> Fill();


}


void test::NeutrinoTopologyAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of required member function here.
  // Here you add an external fcl file to change configuration
}


void test::NeutrinoTopologyAnalyzer::beginJob( )
{
  clearVariables();
  // Declare trees and branches
  mcparticle_tree   = new TTree( "mcparticle_tree",    "MC event tree: True event track information" ) ;
  recoevent_tree = new TTree( "recoevent_tree",  "Reco event tree: reconstructed information of event") ;

  //set branches
  mcparticle_tree -> Branch( "event_id",                &event_id,            "event_id/I");
  mcparticle_tree -> Branch( "Tnu_PDG",                 &Tnu_PDG,             "Tnu_PDG/I");
  mcparticle_tree -> Branch( "T_interaction",           &T_interaction,       "T_interaction/I");
  mcparticle_tree -> Branch( "t_vertex",                &t_vertex,            "tvertex[3]/D");
  mcparticle_tree -> Branch( "t_vertex_energy",         &t_vertex_energy,     "t_vertex_energy/D");
  mcparticle_tree -> Branch( "is_cc",                   &is_cc,               "is_cc/B");
  mcparticle_tree -> Branch( "TPDG_Primary",            &TPDG_Primary,        "TPDG_Primary[10000]/I");
  mcparticle_tree -> Branch( "TNumDaughPrimary",        &TNumDaughPrimary,    "TNumDaughPrimary[10000]/I");
  mcparticle_tree -> Branch( "TNumDaughters",           &TNumDaughters,       "TNumDaughters/I");


  recoevent_tree -> Branch( "event_id",                 &event_id,            "event_id/I");
  recoevent_tree -> Branch( "nu_reconstructed",         &nu_reconstructed,    "nu_reconstructed/B");
  recoevent_tree -> Branch( "numb_nu_reco",             &numb_nu_reco,        "numb_nu_reco/I");
  recoevent_tree -> Branch( "nu_reco_vertex",           &nu_reco_vertex,      "nu_reco_vertex[3]/D");
  recoevent_tree -> Branch( "error_vertex_reco",        &error_vertex_reco,   "error_vertex_reco/D");


  //set directory
  mcparticle_tree -> SetDirectory(0);
  recoevent_tree -> SetDirectory(0);
}


void test::NeutrinoTopologyAnalyzer::endJob( )
{
  // Implementation of required member function here.
  // Write and close files. Other comments also fit here 
  TFile file("output_eventtree.root" , "RECREATE" );
  mcparticle_tree ->Write();
  recoevent_tree ->Write();
  
  file.Write();
  file.Close();
  
  delete mcparticle_tree ;
  delete recoevent_tree ; 
}

void test::NeutrinoTopologyAnalyzer::clearVariables() {

  // Define default for parameters and create variables and trees
  // Detector Geometry
  DetectorHalfLengthX = 400 ;
  DetectorHalfLengthY = 400 ;
  DetectorHalfLengthZ = 500 ;
  CoordinateOffSetX = 200 ; 
  CoordinateOffSetY = 200 ; 
  CoordinateOffSetZ = 0 ; 
  SelectedBorderX = 16.5 ;
  SelectedBorderY = 15 ;
  SelectedBorderZ = 15 ; //47.5 ;

  // Clear variables
  event_id = 0 ;
  Tnu_PDG = 0 ;
  T_interaction = 0 ;
  t_vertex[0] = 0 ;
  t_vertex[1] = 0 ;
  t_vertex[2] = 0 ;
  t_momentum[0] = 0 ;
  t_momentum[1] = 0 ;
  t_momentum[2] = 0 ;
  t_vertex_energy = 0 ;
  is_cc = false ; 
  TPDG_Code = -999 ; 
  TNumDaughters = 0 ;
  TDaughter_mu = 0 ;
  TDaughter_pi = 0 ;
  TDaughter_e = 0 ;
  TDaughter_p = 0 ;
  TDaughter_n = 0 ;
  TDaughter_photon = 0 ;
  TDaughter_other =0 ;
  TrueParticleEnergy = 0 ;
  TMass = 0 ;
  Tpx = 0 ;
  Tpy = 0 ;
  Tpz = 0 ;
  Tpt = 0 ;
  Tp  = 0 ; 
  TMCLength = 0 ;
  TTrack_vertex_x = 0 ;
  TTrack_vertex_y = 0 ;
  TTrack_vertex_z = 0 ;
  TTrack_vertex_t = 0 ;
  TTrack_end_x = 0 ;
  TTrack_end_y = 0 ;
  TTrack_end_z = 0 ;
  TTrack_end_t = 0 ;

  for ( int i = 0 ; i < 10000 ; ++i ){
    TPDG_Primary[i] = 0 ;
    TNumDaughPrimary[i] = 0 ;
  }

  // RECO INFO

  is_reconstructed = false ;
  primary_vcontained = false ;
  primary_econtained = false ;
  for ( int i = 0 ; i < 10000 ; ++i ){
    r_pdg_primary[i] = 0 ;
    rnu_hits[i] = 0;
  }

  nu_reconstructed = false ;
  nu_vertex_contained = true ; 
  numb_nu_reco = 0 ;
  
  for( int i = 0 ; i < 3 ; ++i ){
    nu_reco_vertex[i] = 0 ;
  }

  error_vertex_reco = 0 ;
}

DEFINE_ART_MODULE(test::NeutrinoTopologyAnalyzer)
