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
  void StoreInformation( art::Event const & e, art::Handle< std::vector< recob::Track > > const & trackHandle, art::Handle< std::vector< recob::Shower > > const & showerHandle, art::FindManyP< recob::Track > const & findTracks,  std::map< int , std::vector< int > > & ShowerMothers, int const & part_id_f , int const & primary_daughter) ;
  bool IsMuonPionCandidateChi2( art::Ptr<anab::ParticleID> const & pid_f);
  bool IsMuonPionCandidatePIDA( art::Ptr<anab::ParticleID> const & pid_f);
  double EfficiencyCalo( art::Ptr<anab::ParticleID> const & pid_f , int const & true_pdg , std::string const & particle ) ;

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
  TTree * event_tree, * mcparticle_tree, * recoevent_tree ; 
  
  // Event Information
  int event_id ;
  bool is_reconstructed, has_reco_tracks, has_reco_showers ; 
  int has_reco_daughters ; 

  // MC Event Information 
  // Particle inventory service
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  std::map< int , std::vector< int > > ShowerMothers ;
  std::map<int, const simb::MCParticle*> trueParticles ; 

  // neutrino information
  int Tnu_PDG, T_interaction ;
  double t_vertex[3], t_momentum[3], t_vertex_energy ;
  bool is_cc ; 

  // Truth information of primary particles 
  std::map < int, int > mapTDaughters, mapTRescatter, mapTPrimary; 
  std::map < int, double > mapTLength;
  //  std::map < int, bool > mapTPrimary ;
  bool Thas_primary_mu, Thas_primary_pi, Tmuon_decay, Tdecay_e, Tdecay_nue, Tdecay_numu ; 

  // Reco information
  lar_pandora::PFParticleMap particleMap ;   

  std::map< int , int > mapMC_reco_pdg ;

  bool primary_vcontained, primary_econtained ;
  int r_pdg_primary[10000] ;
  int rnu_hits[10000] ;
  int rtrack_hits;

  //    => neutrino vertex information
  bool nu_reconstructed , nu_rvertex_contained ;
  int numb_nu_reco ;
  bool is_candidate ;
  double nu_reco_vertex[3] ;
  double error_vertex_reco ;

  // Efficiency calculation: just calorimetry information
  int reco_mu , reco_pi ;
  int true_mu, true_pi ;
  int signal_mu, signal_pi ;
  
  std::vector< int > daughters ;
  std::vector< std::vector< int > > daughter_hiearchy ; 

  // need to reconfigure 
  double r_chi2_mu[10], r_chi2_pi[10], r_chi2_p[10], r_PIDA[10] ;
  double r_missenergy[10], r_KineticEnergy[10], r_Range[10] ;
  float rLength ;
  double r_track_x[100000], r_track_y[100000], r_track_z[100000], r_track_dEdx[100000];
  std::vector< float > r_dEdx_ID, r_dEdx_total ; 
  std::vector< std::vector< float >  > r_dEdx_ID_all ; 
  
  // -> Breakdown particles in event from Pandora
  int tr_id_energy, tr_id_charge, tr_id_hits;
  int pfps_truePDG[1000] ;
  int event_vcontained[1000], event_econtained[1000] ; 
  int pfps_hits[1000] , pfps_type[1000] ;
  float pfps_length[1000] ; 
  double pfps_dir_start_x[1000], pfps_dir_start_y[1000], pfps_dir_start_z[1000], pfps_dir_end_x[1000] , pfps_dir_end_y[1000] , pfps_dir_end_z[1000],pfps_start_x[1000] , pfps_start_y[1000] , pfps_start_z[1000] , pfps_end_x[1000], pfps_end_y[1000] , pfps_end_z[1000] ;
  
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
  event_id = e.id().event();
  std::cout<< " Event ID = " << event_id <<std::endl;

  if( !e.isRealData()){
    /**************************************************************************************************
     *  MC INFORMATION
     *
     *  - > information stored : pdg, daughters, primarys id , length, rescatter
     *  - > Want to have a estimate of muon decay after neutrino interactions
     *************************************************************************************************/
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(G4Label);
    
    //List the particles in the event                                                                                                              
    const sim::ParticleList& particles = particleInventory->ParticleList(); // should change. Doing the same twice 
    
    //Make a map of Track id and pdgcode                                                                                                             
    for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
      const simb::MCParticle *particle = particleIt->second;
      trueParticles[particle->TrackId()] = particle;
    }

    for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
      const simb::MCParticle *particle_temp = particleIt->second;

      // Creating shower mother map:      
      if( particle_temp->Mother() != 0 ){
	if( trueParticles.find( particle_temp->TrackId() ) == trueParticles.end() || 
            trueParticles.find( particle_temp->Mother()) == trueParticles.end() ){ continue ; }

	while( particle_temp->Mother() != 0 && (TMath::Abs(trueParticles[particle_temp->Mother()]->PdgCode()) == 11 || 
               trueParticles[particle_temp->Mother()]->PdgCode() == 22 || trueParticles[particle_temp->Mother()]->PdgCode() == 111 )){
	  
	  particle_temp =  trueParticles[particle_temp->Mother()];
	  if( trueParticles.find(particle_temp->Mother()) == trueParticles.end()){ break; } // found mother
	}
      }
      if( ShowerMothers.find( particle_temp->TrackId() ) == ShowerMothers.end() 
          && (TMath::Abs(trueParticles[particle_temp->TrackId()]->PdgCode()) == 11 || trueParticles[particle_temp->TrackId()] -> PdgCode() == 22 
          || trueParticles[particle_temp->TrackId()] -> PdgCode() == 111 ) ) { 
	      ShowerMothers[particle_temp->TrackId()].push_back(particle_temp->TrackId()) ; 
      }
    }

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
	// Storing all particles in maps -> [ TrackId vs variable ]
	const simb::MCParticle trueParticle = mcParticles->at(t) ;
	if( trueParticle.PdgCode() >= 1000018038 ) continue ; // Cut on PDG codes which refer to elements (Argon30 and above)
	mapMC_reco_pdg[trueParticle.TrackId()] = trueParticle.PdgCode() ;
	mapTLength[trueParticle.TrackId()] = trueParticle.Trajectory().TotalLength() ;
	mapTDaughters[trueParticle.TrackId()] = trueParticle.NumberDaughters();
	mapTRescatter[trueParticle.TrackId()] = trueParticle.Rescatter();

	if(trueParticle.Process() == "primary" ) mapTPrimary[trueParticle.TrackId()] = 1 ;
	else mapTPrimary[trueParticle.TrackId()] = 2 ;

	if(trueParticle.Process() == "primary" && trueParticle.PdgCode() == 13  ) Thas_primary_mu = true ;
	if(trueParticle.Process() == "primary" && trueParticle.PdgCode() == 211 ) Thas_primary_pi = true ;
	
	if( trueParticle.Process() != "primary" && trueParticle.Mother() != 0 && TMath::Abs(mapMC_reco_pdg[trueParticle.Mother()]) == 13 ){
	  if( trueParticle.PdgCode() == 11 ) Tdecay_e = true ;
	  if( trueParticle.PdgCode() == -12 ) Tdecay_nue = true ;
	  if( trueParticle.PdgCode() == 14 ) Tdecay_numu = true ;
	  if( Tdecay_e && Tdecay_nue && Tdecay_numu ) Tmuon_decay = true ;
	}

      }
    }
    
  }  

  /**************************************************************************************************
   *  RECO INFORMATION
   *
   *  - > information stored : 
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
	      || ( vtx_assn[0]->position().Z() < (-CoordinateOffSetZ + SelectedBorderZ))) nu_rvertex_contained = false ;

	  for( int j = 0 ; j < pfparticle->NumDaughters() ; ++j ) { // loop over neutrino daughters
		int part_id_f = particleMap[ pfparticle->Daughters()[j] ] -> Self() ;
		is_candidate = false ; // reset for each pfparticle
		pfps_type[j] = particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() ; // this is the pandora pdg code
		StoreInformation( e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_f , j ) ;

		if( is_candidate == true ) { // just check possible muons and pions !
		  //		  daughters.push_back( particleMap[pfparticle->Daughters()[j] ] -> Self() ) ;
		  daughters.push_back( particleMap[pfparticle->Daughters()[j] ] -> NumDaughters() ) ;
		  if( particleMap[pfparticle->Daughters()[j] ] -> NumDaughters() == 1 ){ // at the moment just look more if there is one track. Splitted -> pion track 
		    for( int j2 = 0 ; j2 < particleMap[pfparticle->Daughters()[j] ] -> NumDaughters() ; ++j2 ) {
		      daughters.push_back(particleMap[ particleMap[ pfparticle->Daughters()[j] ]->Daughters()[j2] ] -> NumDaughters() ) ;
		    }
		  }
		  daughter_hiearchy.push_back(daughters);
		  daughters.clear();
		}
	    }
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

void test::NeutrinoTopologyAnalyzer::StoreInformation( 
      art::Event const & e, art::Handle< std::vector< recob::Track > > const & trackHandle, 
      art::Handle< std::vector< recob::Shower > > const & showerHandle, art::FindManyP< recob::Track > const & findTracks, 
      std::map< int , std::vector< int > > & ShowerMothers, int const & part_id_f , int const & primary_daughter) {

  // Save track info
  if ( findTracks.at( part_id_f ).size() != 0 ){
    std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(part_id_f);
    art::FindManyP< recob::Hit > findHits (  trackHandle, e, RecoTrackLabel ) ;
    art::FindManyP< anab::Calorimetry > findCalorimetry ( trackHandle, e, RecoCaloLabel );
    art::FindManyP< anab::ParticleID > findPID ( trackHandle, e, RecoPIDLabel );

    // Loop over tracks found for track_f
    for( unsigned int n = 0 ; n < track_f.size() ; ++n ){
      has_reco_tracks = true ; 
      rLength   += track_f[n]->Length() ; // add end of track
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
       	    // Get associated MCParticle ID using 3 different methods:
	    //    Which particle contributes the most energy to all the hits
	    //    Which particle contributes the reco charge to all the hits
	    //    Which particle is the biggest contributor to all the hits
	    tr_id_energy      = RecoUtils::TrueParticleIDFromTotalTrueEnergy(hit_f);
	    tr_id_charge      = RecoUtils::TrueParticleIDFromTotalRecoCharge(hit_f);
	    tr_id_hits        = RecoUtils::TrueParticleIDFromTotalRecoHits(hit_f);
	    // save the most common answer : 
	    if( tr_id_energy == tr_id_charge && tr_id_energy == tr_id_hits ) pfps_truePDG[primary_daughter] = mapMC_reco_pdg[tr_id_energy] ;
	    if( tr_id_energy == tr_id_charge && tr_id_energy != tr_id_hits ) pfps_truePDG[primary_daughter] = mapMC_reco_pdg[tr_id_energy] ;
	    if( tr_id_energy != tr_id_charge && tr_id_energy == tr_id_hits ) pfps_truePDG[primary_daughter] = mapMC_reco_pdg[tr_id_energy] ;
	    if( tr_id_energy != tr_id_charge && tr_id_charge == tr_id_hits ) pfps_truePDG[primary_daughter] = mapMC_reco_pdg[tr_id_charge] ;
	    if( tr_id_energy != tr_id_charge && tr_id_energy != tr_id_hits && tr_id_charge != tr_id_hits) pfps_truePDG[primary_daughter] = mapMC_reco_pdg[tr_id_hits] ;
	    //	    std::cout << " pdg -> " << pfps_truePDG[primary_daughter] << " Is candidate =  "<< IsMuonPionCandidateChi2( pid_f[k] ) << std::endl; 
	    //	    std::cout<< "efficiency muon " << EfficiencyCalo( pid_f[k] , pfps_truePDG[primary_daughter], "muon" ) << std::endl;
	    
	    if( IsMuonPionCandidateChi2( pid_f[k] ) == 0 ) {
	      is_candidate = false ;
	      continue ; // just read potential muon/pion tracks
	    } else is_candidate = true ;
	    // save information -- add pid methods here !
	    
	    /*
	      1) Find muon/pion candidates
	      2) Efficiency simple method
	      3) Apply topology consideretions -> New efficiency
	      4) Number of kinks -> Efficiency
	      5) Michael electrons -> Efficiency
	     */

	    r_chi2_mu[primary_daughter] = pid_f[k]->Chi2Muon() ;
	    r_chi2_pi[primary_daughter] = pid_f[k]->Chi2Pion() ;
	    r_chi2_p[primary_daughter]  = pid_f[k]->Chi2Proton() ;
	    r_PIDA[primary_daughter]    = pid_f[k]->PIDA();
	    r_missenergy[primary_daughter] = pid_f[k]->MissingE();
	    r_KineticEnergy[primary_daughter] = cal_f[m]->KineticEnergy();
	    
	    for( unsigned int l = 0 ; l < track_f[n]->LastValidPoint()+1 ; ++l ) {
	      r_track_x[l+rtrack_hits] = track_f[n]->TrajectoryPoint( l ).position.X();
	      r_track_y[l+rtrack_hits] = track_f[n]->TrajectoryPoint( l ).position.Y();
	      r_track_z[l+rtrack_hits] = track_f[n]->TrajectoryPoint( l ).position.Z();
	    }
	    
	    rtrack_hits   += track_f[n]->LastValidPoint() + 1 ;
	    pfps_hits[primary_daughter] = track_f[n]->LastValidPoint() + 1 ;
	    pfps_length[primary_daughter] = track_f[n]->Length() ;
	    pfps_dir_start_x[primary_daughter] = track_f[n]->StartDirection().X() ;
	    pfps_dir_start_y[primary_daughter] = track_f[n]->StartDirection().Y() ;
	    pfps_dir_start_z[primary_daughter] = track_f[n]->StartDirection().Z() ;
	    pfps_dir_end_x[primary_daughter] = track_f[n]->EndDirection().X() ;
	    pfps_dir_end_y[primary_daughter] = track_f[n]->EndDirection().Y() ;
	    pfps_dir_end_z[primary_daughter] = track_f[n]->EndDirection().Z() ;
	    pfps_start_x[primary_daughter] = track_f[n]->Start().X() ;
	    pfps_start_y[primary_daughter] = track_f[n]->Start().Y() ;
	    pfps_start_z[primary_daughter] = track_f[n]->Start().Z() ;
	    pfps_end_x[primary_daughter] = track_f[n]->End().X() ;
	    pfps_end_y[primary_daughter] = track_f[n]->End().Y() ;
	    pfps_end_z[primary_daughter] = track_f[n]->End().Z() ;
	    
	  }// just collection plane 
          
	  // calo information is stored in all planes:
          
	  for( unsigned int l = 0 ; l < (cal_f[m]->XYZ()).size() ; ++l ) {
	    for( int t = 0 ; t < rtrack_hits ; ++t ){
	      if( cal_f[m]->XYZ()[l].X() == r_track_x[t] && cal_f[m]->XYZ()[l].Y() == r_track_y[t] && cal_f[m]->XYZ()[l].Z() == r_track_z[t]) {
		r_track_dEdx[t] = cal_f[m]->dEdx()[l] ; 
	      }
	    }
	  }
	  r_Range[primary_daughter] = cal_f[m]->Range(); // check 
	} //close calo
      } //close pid
    } //close track    
  } else if( showerHandle.isValid() && showerHandle->size() != 0 ) { // if no track look into showers 
    //std::cout<< " is shower " << std::endl;
    has_reco_showers = true ; 
    art::FindManyP< recob::Hit > findHitShower( showerHandle, e, RecoShowerLabel ) ;
    art::FindManyP< recob::SpacePoint > findSpacePoint( showerHandle, e, RecoShowerLabel ) ;
    for( unsigned int y = 0 ; y < showerHandle->size() ; ++y ) {
      art::Ptr< recob::Shower > shower_f( showerHandle, y ) ;
      std::vector< art::Ptr<recob::Hit> > hit_sh_f = findHitShower.at(y) ; 
      std::vector< art::Ptr<recob::SpacePoint> > spacepoint_f = findSpacePoint.at(y) ;
      if( spacepoint_f.size() == 0 ) continue ; 

      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain( ShowerMothers, hit_sh_f , shower_f->best_plane() ) ;

      for( unsigned int l = 0 ; l < spacepoint_f.size() ; ++l ) {
	r_track_x[l+rtrack_hits] = spacepoint_f[l]->XYZ()[0] ;
	r_track_y[l+rtrack_hits] = spacepoint_f[l]->XYZ()[1] ;
	r_track_z[l+rtrack_hits] = spacepoint_f[l]->XYZ()[2] ;
      }

      rtrack_hits   += spacepoint_f.size() ;
      pfps_hits[primary_daughter] = spacepoint_f.size() ;

      if( shower_f->has_length() ) { 
	pfps_length[primary_daughter] = shower_f->Length() ;
      } else  {
	pfps_length[primary_daughter]  = pow(spacepoint_f[spacepoint_f.size()-1]->XYZ()[0] - spacepoint_f[0]->XYZ()[0], 2 ) ;
	pfps_length[primary_daughter] += pow(spacepoint_f[spacepoint_f.size()-1]->XYZ()[1] - spacepoint_f[0]->XYZ()[1], 2 ) ;
	pfps_length[primary_daughter] += pow(spacepoint_f[spacepoint_f.size()-1]->XYZ()[2] - spacepoint_f[0]->XYZ()[2], 2 ) ;
	pfps_length[primary_daughter]  = sqrt( pfps_length[primary_daughter] ) ;
      }
          
      pfps_dir_start_x[primary_daughter] = shower_f->Direction().X() ;
      pfps_dir_start_y[primary_daughter] = shower_f->Direction().Y() ;
      pfps_dir_start_z[primary_daughter] = shower_f->Direction().Z() ;
      // No end Direction
      pfps_start_x[primary_daughter] = shower_f->ShowerStart().X() ;
      pfps_start_y[primary_daughter] = shower_f->ShowerStart().Y() ;
      pfps_start_z[primary_daughter] = shower_f->ShowerStart().Z() ;
      // no end position
      if( ShowerTrackInfo.first ) { pfps_truePDG[primary_daughter] = mapMC_reco_pdg[ShowerTrackInfo.first] ; }
      else pfps_truePDG[primary_daughter] = 0 ;
      //std::cout<<" --> Shower pdg = " << pfps_truePDG[primary_daughter]  << std::endl; 

    }
  } // track vs shower
}


bool test::NeutrinoTopologyAnalyzer::IsMuonPionCandidateChi2( art::Ptr<anab::ParticleID> const & pid_f )
{
  if( ( pid_f ->Chi2Muon() < pid_f ->Chi2Proton() && pid_f ->Chi2Muon() < pid_f ->Chi2Kaon() ) || ( pid_f ->Chi2Pion() < pid_f ->Chi2Proton() && pid_f ->Chi2Pion() < pid_f ->Chi2Kaon() ) ) return true ;
	 
  return false ; 
}


bool test::NeutrinoTopologyAnalyzer::IsMuonPionCandidatePIDA( art::Ptr<anab::ParticleID> const & pid_f )
{
  /*
  // Muon : PIDA >= 5 && PIDA() < 9 
  //Pion : PIDA() >= 9 && PIDA() < 13 
  */
  if( pid_f ->PIDA() >= 5 && pid_f ->PIDA() < 13 ) return true ;
  return false ; 
}

double test::NeutrinoTopologyAnalyzer::EfficiencyCalo( art::Ptr<anab::ParticleID> const & pid_f , int const & true_pdg , std::string const & particle ) {
  double eff_mu, eff_pi, purity_mu, purity_pi ;
  // True information 
  if( true_pdg == 13 ) ++true_mu ;
  if( true_pdg == 211 ) ++true_pi ; 
  // reco information
  if( pid_f ->Chi2Muon() < pid_f ->Chi2Proton() && pid_f ->Chi2Muon() < pid_f ->Chi2Kaon() && pid_f ->Chi2Muon() < pid_f -> Chi2Pion() ) {
    ++reco_mu ;
    if( true_pdg == 13 ) ++signal_mu ;
  }

  if( pid_f ->Chi2Muon() < pid_f ->Chi2Proton() && pid_f ->Chi2Muon() < pid_f ->Chi2Kaon() && pid_f ->Chi2Pion() < pid_f -> Chi2Muon() ) {
    ++reco_pi ;
    if( true_pdg == 211 ) ++signal_pi ;
  }

  eff_mu = signal_mu / (double)true_mu ; 
  purity_mu = signal_mu / (double)reco_mu ; 
  eff_pi = signal_pi / (double)true_pi ; 
  purity_pi = signal_pi / (double)reco_pi ; 
  std::cout<<"reco mu " << reco_mu << " true mu " << true_mu << " signal mu " << signal_mu << std::endl;
  std::cout<<"reco pi " << reco_pi << " true pi " << true_pi << " signal pi " << signal_pi << std::endl;
  if( particle == "pion" ) return eff_pi ; 
  if( particle == "purity pion" ) return purity_pi ; 
  if( particle == "purity muon" ) return purity_mu ; 
  return eff_mu;
}

void test::NeutrinoTopologyAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of required member function here.
  // Here you add an external fcl file to change configuration
}


void test::NeutrinoTopologyAnalyzer::beginJob( )
{
  // don't initialize per event. 
  reco_mu = 0 ;
  reco_pi = 0 ;
  true_mu = 0 ;
  true_pi = 0 ;
  signal_mu = 0 ;
  signal_pi = 0 ; 

  clearVariables();
  // Declare trees and branches
  event_tree   = new TTree( "event_tree",    "Event tree: General information about the event" ) ;
  mcparticle_tree   = new TTree( "mcparticle_tree",    "MC event tree: True event track information" ) ;
  recoevent_tree = new TTree( "recoevent_tree",  "Reco event tree: reconstructed information of event") ;

  // Event tree
  event_tree      -> Branch( "event_id",                  &event_id,           "event_id/I");
  event_tree      -> Branch( "is_reconstructed",          &is_reconstructed,   "is_reconstructed/B");
  event_tree      -> Branch( "has_reco_daughters",        &has_reco_daughters, "has_reco_daughters/I");
  event_tree      -> Branch( "has_reco_tracks",           &has_reco_tracks,    "has_reco_tracks/B");
  event_tree      -> Branch( "has_reco_showers",          &has_reco_showers,   "has_reco_showers/B");

  // MC tree
  mcparticle_tree -> Branch( "event_id",                &event_id,            "event_id/I");
  mcparticle_tree -> Branch( "Tnu_PDG",                 &Tnu_PDG,             "Tnu_PDG/I");
  mcparticle_tree -> Branch( "T_interaction",           &T_interaction,       "T_interaction/I");
  mcparticle_tree -> Branch( "t_vertex",                &t_vertex,            "tvertex[3]/D");
  mcparticle_tree -> Branch( "t_vertex_energy",         &t_vertex_energy,     "t_vertex_energy/D");
  mcparticle_tree -> Branch( "is_cc",                   &is_cc,               "is_cc/B");
  mcparticle_tree -> Branch( "mapMC_reco_pdg",          "std::map<int,int>",   &mapMC_reco_pdg);
  mcparticle_tree -> Branch( "mapTDaughters",           "std::map<int,int>",   &mapTDaughters);
  mcparticle_tree -> Branch( "mapTRescatter",           "std::map<int,int>",   &mapTRescatter); 
  mcparticle_tree -> Branch( "mapTPrimary",             "std::map<int,int>",   &mapTPrimary);
  mcparticle_tree -> Branch( "mapTLength",              "std::map<int,double>",&mapTLength);
  mcparticle_tree -> Branch( "Tmuon_decay",             &Tmuon_decay,          "Tmuon_decay/B");
  mcparticle_tree -> Branch( "Thas_primary_mu",         &Thas_primary_mu,      "Thas_primary_mu/B");
  mcparticle_tree -> Branch( "Thas_primary_pi",         &Thas_primary_pi,      "Thas_primary_pi/B");


  // Reco tree
  recoevent_tree -> Branch( "event_id",                 &event_id,            "event_id/I");
  recoevent_tree -> Branch( "nu_reconstructed",         &nu_reconstructed,    "nu_reconstructed/B");
  recoevent_tree -> Branch( "nu_rvertex_contained",     &nu_rvertex_contained,"nu_rvertex_contained/B");
  recoevent_tree -> Branch( "numb_nu_reco",             &numb_nu_reco,        "numb_nu_reco/I");
  recoevent_tree -> Branch( "nu_reco_vertex",           &nu_reco_vertex,      "nu_reco_vertex[3]/D");
  recoevent_tree -> Branch( "error_vertex_reco",        &error_vertex_reco,   "error_vertex_reco/D");
  recoevent_tree -> Branch( "daughter_hiearchy" , "std::vector< std::vector< int > >", &daughter_hiearchy);


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

  // event tree variables 
  event_id = 0 ;
  is_reconstructed   =  false ; 
  has_reco_daughters =  0 ;
  has_reco_tracks    =  false ;
  has_reco_showers = false ; 

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
  mapTDaughters.clear() ;
  mapMC_reco_pdg.clear();
  mapTLength.clear();
  mapTPrimary.clear();
  mapTRescatter.clear();
  Tmuon_decay = false ; 
  Tdecay_e = false ; 
  Tdecay_nue = false ; 
  Tdecay_numu = false ; 
  Thas_primary_mu = false ;
  Thas_primary_pi = false ;  

  // RECO INFO

  primary_vcontained = false ;
  primary_econtained = false ;
  is_candidate = false ; 
  for ( int i = 0 ; i < 10000 ; ++i ){
    r_pdg_primary[i] = 0 ;
    rnu_hits[i] = 0;
  }

  rtrack_hits = 0 ;
  nu_reconstructed = false ;
  nu_rvertex_contained = true ; 
  numb_nu_reco = 0 ;
  
  for( int i = 0 ; i < 3 ; ++i ){
    nu_reco_vertex[i] = 0 ;
  }

  error_vertex_reco = 0 ;
  daughters.clear() ;
  daughter_hiearchy.clear() ; 


}

DEFINE_ART_MODULE(test::NeutrinoTopologyAnalyzer)
