/////////////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoTopologyAnalyzer                                       //
// Plugin Type: analyzer (art v3_00_00)                                        //
// File:        NeutrinoTopologyAnalyzer_module.cc                             //
//                                                                             //
// Generated at Mon Mar 25 06:08:16 2019 by Julia Tena Vidal using cetskelgen  //
// from cetlib version v3_04_00.                                               //
/////////////////////////////////////////////////////////////////////////////////


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
#include "THStack.h"
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
#include <iostream>
#include <fstream>

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
  void StoreInformation( art::Event const & e, art::Handle< std::vector< recob::Track > > const & trackHandle, art::Handle< std::vector< recob::Shower > > const & showerHandle, art::FindManyP< recob::Track > const & findTracks,  std::map< int , std::vector< int > > & ShowerMothers, int const & part_id_f ) ;
  bool IsMuonPionCandidateChi2( art::Ptr<anab::ParticleID> const & pid_f );
  bool IsMuonChi2( art::Ptr<anab::ParticleID> const & pid_f );
  bool IsPionChi2( art::Ptr<anab::ParticleID> const & pid_f );
  bool IsMuonPionCandidateChi2Proton( art::Ptr<anab::ParticleID> const & pid_f );
  bool IsMuonPionCandidatePIDA( art::Ptr<anab::ParticleID> const & pid_f);
  double EfficiencyCalo( art::Ptr<anab::ParticleID> const & pid_f , int const & true_pdg , std::string const & particle ) ;
  void IsReconstructed( int  const & best_id ) ;
  int FindBestMCID(art::Event const & e, art::Handle< std::vector< recob::Track > > const & trackHandle, 
		  art::Handle< std::vector< recob::Shower > > const & showerHandle, art::FindManyP< recob::Track > const & findTracks, 
		  std::map< int , std::vector< int > > & ShowerMothers, int const & part_id_f ) ;  
  std::vector< bool > MCIsContained( simb::MCParticle const & trueParticle ) ;
  std::vector< bool > IsContained( art::Event const & e, art::Handle< std::vector< recob::Track > > & trackHandle, art::Handle< std::vector< recob::Shower > > & showerHandle, art::FindManyP< recob::Track > & findTracks , int & part_id_f ) ;


  void reconfigure(fhicl::ParameterSet const & p);
  //  void clearVariables() ;
  void beginJob() override;
  void endJob() override;
  void clearVariables() ; 


private:
  // Declare member data here ;
  
  /******************************************
   *  LABELS                                *
   ******************************************/
  std::string TruthLabel, G4Label, ParticleLabel, HitFinderLabel, RecoTrackLabel, RecoShowerLabel, RecoPIDLabel, RecoCaloLabel ; 


  /******************************************
   *  DETECTOR INFORMATION                  *
   ******************************************/
  float DetectorHalfLengthX, DetectorHalfLengthY, DetectorHalfLengthZ, CoordinateOffSetX, CoordinateOffSetY, CoordinateOffSetZ, SelectedBorderX, SelectedBorderY, SelectedBorderZ ;

  /******************************************
   *  TREE DEFINITIONS                      *
   ******************************************/
  TTree * event_tree, * mcparticle_tree, * recoevent_tree ; 
  
  /******************************************
   *  GENERAL EVENT INFORMATION             *
   ******************************************/
  int event_id ;
  bool is_reconstructed, has_reco_tracks, has_reco_showers ; 
  int has_reco_daughters ; 

  /******************************************
   *  MC INFORMATION                        *
   ******************************************/
  
  // Particle inventory service
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  std::map< int , std::vector< int > > ShowerMothers ;
  std::map<int, const simb::MCParticle*> trueParticles ; 

  // neutrino information
  int Tnu_PDG, T_interaction ;
  double t_vertex[3], t_momentum[3], t_vertex_energy ;
  bool is_cc ; 

  // Truth information of primary particles 
  std::map < int, int >  map_IsReconstructed, mapTDaughters, mapTRescatter, mapTPrimary, mapTPiDaughterPdg; 
  std::map < int, double > mapTLength;
  //  std::map < int, bool > mapTPrimary ;
  bool Thas_primary_mu, Thas_primary_pi, Tmuon_decay, Tdecay_e, Tdecay_nue, Tdecay_numu ; 
  
  /******************************************
   *  RECONSTRUCTED INFORMATION             *   
   ******************************************/
  
  lar_pandora::PFParticleMap particleMap ;   

  // neutrino related quantities
  bool primary_vcontained, primary_econtained ;
  bool nu_reconstructed , nu_rvertex_contained ;
  int numb_nu_reco ;
  double nu_reco_vertex[3] ;
  double error_vertex_reco ;
  
  // reco - true id matching variables. tr_id_best -> best of the three methods
  int tr_id_best, tr_id_energy, tr_id_charge, tr_id_hits;

  // Particle True - Reco map 
  std::map< int , int > mapMC_reco_pdg ;
  
  // muon - pion candidates stored information :
  bool is_candidate ;
  std::map< int , int > map_MCID_RecoID ;
  std::map< int , std::vector<bool> > map_RecoContained ; 
  std::map< int , int > map_RecoHits, map_RecoPrimary, map_PandoraPDG ; 
  std::map< int , double > map_RecoLength, map_RecoKEnergy ; 
  std::map< int , std::vector< double > > map_RecoXPosition, map_RecoYPosition, map_RecoZPosition, map_RecodEdx ;
  std::map< int , std::vector< int > > map_RecoDaughters ; 
  std::map< int , int > map_RecoHiearchy ; 

  // Efficiency calculation: just calorimetry information
  int reco_primary_mu , reco_primary_pi ;
  int true_mu, true_pi, true_p ; // total pdg-particle in event. Includes secondaries 
  int true_primary_mu, true_primary_pi ;
  int signal_mu, signal_pi ;
  int bg_mu_pi, bg_mu_p, bg_mu_others ; 
  int bg_pi_mu, bg_pi_p, bg_pi_others ; 
  double eff_mu, eff_pi, purity_mu, purity_pi ;


  // Is reconstructed information ( look if particles are reconstructed by pandora )
  int reco_track_mu, reco_track_pi, reco_track_p ; 
  

  /******************************************
   * HISTOGRAMS FOR STUDIES                 *   
   ******************************************/
  TCanvas *c = new TCanvas();
  TLegend *l = new TLegend( 0.58, 0.68, 0.88, 0.88 );

  TH1D * TFinalStatePi = new TH1D("TFinalStatePi", " True Final State Pion " , 50 , 0 , 3000 ); // could be proton

  TH1D * Chi2p_Tmu = new TH1D("Chi2p_Tmu", " CHI2 under PROTON HYPOTESIS FOR MUONS"   , 50 , 0 , 200 );
  TH1D * Chi2p_Tpi = new TH1D("Chi2p_Tpi", " CHI2 under PROTON HYPOTESIS FOR PIONS"   , 50 , 0 , 200 );
  TH1D * Chi2p_Tp  = new TH1D("Chi2p_Tp",  " CHI2 under PROTON HYPOTESIS FOR PROTONS" , 50 , 0 , 200 );

  int total_reco_p , reco_p ; 

  // Information about TPC reconstructed particles
  TH1D * h_MCLength_mu_TPC_signal = new TH1D("MCLenght_mu_TPC_signal", "mu MC Length, TPC Signal" , 20, 0 , 250 ) ;
  TH1D * h_MCLength_pi_TPC_signal = new TH1D("MCLenght_pi_TPC_signal", "pi MC Length, TPC Signal" , 20, 0 , 200 ) ;
  TH1D * h_MCLength_p_TPC_signal  = new TH1D("MCLenght_p_TPC_signal" , "p MC Length, TPC Signal"  , 20, 0 , 200 ) ;
  TH1D * h_MCLength_mu_TPC_miss   = new TH1D("MCLenght_mu_TPC_miss"  , "mu MC Length, TPC Miss"   , 20, 0 , 250 ) ;
  TH1D * h_MCLength_pi_TPC_miss   = new TH1D("MCLenght_pi_TPC_miss"  , "pi MC Length, TPC Miss"   , 20, 0 , 200 ) ;
  TH1D * h_MCLength_p_TPC_miss    = new TH1D("MCLenght_p_TPC_miss"   , "p MC Length, TPC Miss"    , 20, 0 , 200 ) ;

  THStack * h_MCLength_mu_TPC = new THStack("h_MCLength_mu_TPC", "Mu MC Length TPC");
  THStack * h_MCLength_pi_TPC = new THStack("h_MCLength_pi_TPC", "Pi MC Length TPC");
  THStack * h_MCLength_p_TPC  = new THStack("h_MCLength_p_TPC" , "P MC Length TPC");

  // Hiearchy information for a candidate signature
  TH1D * h_recoDaughters_mu = new TH1D("recoDaughters_mu", " Total number of daughters for muons" , 4, -0.5, 3.5 ) ;
  TH1D * h_recoDaughters_pi = new TH1D("recoDaughters_pi", " Total number of daughters for pions" , 4, -0.5, 3.5 ) ;
  TH1D * h_recoDaughters_p  = new TH1D("recoDaughters_p", " Total number of daughters for protons" , 4, -0.5, 3.5 ) ;

  TH1D * h_reco3Daughters_mu = new TH1D("reco3Daughters_mu", " Total number of 3rd-generation daughters" , 4, -0.5, 3.5 ) ;
  TH1D * h_reco3Daughters_pi = new TH1D("reco3Daughters_pi", " Total number of 3rd-generation daughters" , 4, -0.5, 3.5 ) ;
  TH1D * h_reco3Daughters_p  = new TH1D("reco3Daughters_p", " Total number of 3rd-generation daughters" , 4, -0.5, 3.5 ) ;

  // reconstructed information
  // For the time being I will just fill the muon ones. No signal events for pions
  TH1D * h_recoLength_mu = new TH1D("recoLength_mu", " Candidate true Muon reconstructed length ", 50, 0, 400 ) ;
  TH1D * h_recoLength_pi = new TH1D("recoLength_pi", " Candidate true Pion reconstructed length ", 50, 0, 400 ) ;
  TH1D * h_recoLength_p  = new TH1D("recoLength_p",  " Candidate true Proton reconstructed length ", 50, 0, 400 ) ;

  TH1D * h_recoKE_mu = new TH1D("recoKE_mu", " Candidate true Muon reconstructed KE ", 50, 0, 10 ) ;
  TH1D * h_recoKE_pi = new TH1D("recoKE_pi", " Candidate true Pion reconstructed KE ", 50, 0, 10 ) ;
  TH1D * h_recoKE_p  = new TH1D("recoKE_p",  " Candidate true Proton reconstructed KE ", 50, 0, 10 ) ;

  TH1D * h_recoLength_mu_noD = new TH1D("recoLength_mu_noD", " Candidate true Muon reconstructed length with no daughters", 50, 0, 400 ) ;
  TH1D * h_recoLength_pi_noD = new TH1D("recoLength_pi_noD", " Candidate true Pion reconstructed length with no daughters", 50, 0, 400 ) ;
  TH1D * h_recoLength_p_noD  = new TH1D("recoLength_p_noD",  " Candidate true Proton reconstructed length with no daughters", 50, 0, 400 ) ;

  TH1D * h_recoKE_mu_noD = new TH1D("recoKE_mu_noD", " Candidate true Muon reconstructed KE with no daughters", 50, 0, 10 ) ;
  TH1D * h_recoKE_pi_noD = new TH1D("recoKE_pi_noD", " Candidate true Pion reconstructed KE with no daughters", 50, 0, 10 ) ;
  TH1D * h_recoKE_p_noD  = new TH1D("recoKE_p_noD",  " Candidate true Proton reconstructed KE with no daughters", 50, 0, 10 ) ;

  TH1D * h_recoLength_mu_wD = new TH1D("recoLength_mu_wD", " Candidate true Muon reconstructed length with daughters", 50, 0, 400 ) ;
  TH1D * h_recoLength_pi_wD = new TH1D("recoLength_pi_wD", " Candidate true Pion reconstructed length with daughters", 50, 0, 400 ) ;
  TH1D * h_recoLength_p_wD  = new TH1D("recoLength_p_wD",  " Candidate true Proton reconstructed length with daughters", 50, 0, 400 ) ;

  TH1D * h_recoKE_mu_wD = new TH1D("recoKE_mu_wD", " Candidate true Muon reconstructed KE with daughters", 50, 0, 10 ) ;
  TH1D * h_recoKE_pi_wD = new TH1D("recoKE_pi_wD", " Candidate true Pion reconstructed KE with daughters", 50, 0, 10 ) ;
  TH1D * h_recoKE_p_wD  = new TH1D("recoKE_p_wD",  " Candidate true Proton reconstructed KE with daughters", 50, 0, 10 ) ;
  
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
  std::cout<< " ***********************************************************" << event_id <<std::endl;
  std::cout<< " Event ID = " << event_id <<std::endl;
  std::cout<< " ***********************************************************" << event_id <<std::endl;

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
      //std::cout<< " MC INFORMATION : " << std::endl;
      //std::cout<< " ------------------" << std::endl;

      for( unsigned int t = 0; t < mcParticles->size(); ++t ){
	// Storing all particles in maps -> [ TrackId vs variable ]
	const simb::MCParticle trueParticle = mcParticles->at(t) ;
	if( trueParticle.PdgCode() >= 1000018038 ) continue ; // Cut on PDG codes which refer to elements (Argon30 and above)
	mapMC_reco_pdg[trueParticle.TrackId()] = trueParticle.PdgCode() ;
	mapTLength[trueParticle.TrackId()] = trueParticle.Trajectory().TotalLength() ;
	mapTDaughters[trueParticle.TrackId()] = trueParticle.NumberDaughters();
	mapTRescatter[trueParticle.TrackId()] = trueParticle.Rescatter();
	map_IsReconstructed[trueParticle.TrackId()] = 0 ; // creates the map empty  

	// Counting total amount of muons, pions and protons 
	if( TMath::Abs(trueParticle.PdgCode()) == 13 )    ++true_mu ;
	if( TMath::Abs(trueParticle.PdgCode()) == 211  )  ++true_pi ;
	if( trueParticle.PdgCode() == 2212 )              ++true_p  ;

	if(trueParticle.Process() == "primary" ) mapTPrimary[trueParticle.TrackId()] = 1 ;
	else mapTPrimary[trueParticle.TrackId()] = 2 ;

	// muons and pions
	if(trueParticle.Process() == "primary" && TMath::Abs(trueParticle.PdgCode()) == 13  ) { Thas_primary_mu = true ; ++true_primary_mu ; }
	if(trueParticle.Process() == "primary" && TMath::Abs(trueParticle.PdgCode()) == 211 ) { Thas_primary_pi = true ; ++true_primary_pi ;}

	// Study of mu^- decay. We expect a mu^- for numu CC interactions. If studying numubar CC interactions, change decay products.
	if( trueParticle.Process() != "primary" && trueParticle.Mother() != 0 && mapMC_reco_pdg[trueParticle.Mother()] == 13 ){
	  if( trueParticle.PdgCode() == 11 ) Tdecay_e = true ;
	  if( trueParticle.PdgCode() == -12 ) Tdecay_nue = true ;
	  if( trueParticle.PdgCode() == 14 ) Tdecay_numu = true ;
	  if( Tdecay_e && Tdecay_nue && Tdecay_numu ) Tmuon_decay = true ;
	}

	// Study of pi+ elastic and inelastic scattering : final products pdg
	if( trueParticle.Process() != "primary" && trueParticle.Mother() != 0 && TMath::Abs(mapMC_reco_pdg[trueParticle.Mother()]) == 211 ){
	  mapTPiDaughterPdg[trueParticle.PdgCode()] += 1 ;
	  TFinalStatePi -> Fill ( trueParticle.PdgCode() ) ;
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
      if( pfparticle->IsPrimary() == 1 && pfparticle->NumDaughters() != 0 ) { //found primary particle => NEUTRINO real candidate. Has daughers!
	// Get vertex association 
	numb_nu_reco += 1 ; // now checking number of neutrinos with daughters
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
	    int part_MCID_f = FindBestMCID(e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_f ) ; 
	    is_candidate = false ; // reset for each pfparticle daughter of the neutrino
	    StoreInformation( e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_f ) ;
	      
	    if( is_candidate == true ) { // store daugher information for muons and pion candidates
	      map_MCID_RecoID[part_id_f] = part_MCID_f ; 
	      map_RecoHiearchy[part_id_f] = 1 ; // reconstructed as primary 
	      map_PandoraPDG[part_id_f] = particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() ; // this is the pandora pdg code
	      map_RecoContained[part_id_f] = IsContained( e, trackHandle, showerHandle, findTracks, part_id_f ) ;
	    
	      for( int j2 = 0 ; j2 < particleMap[pfparticle->Daughters()[j] ] -> NumDaughters() ; ++j2 ) {
		// secondary particles 
		int part_id_2f = particleMap[ pfparticle->Daughters()[j] ]->Daughters()[j2];
		int part_MCID_2f = FindBestMCID(e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_2f ) ;
		map_MCID_RecoID[part_id_2f] = part_MCID_2f ; 
		map_RecoDaughters[ part_id_f ].push_back( part_id_2f ) ; 
		map_RecoHiearchy[part_id_2f] = 2 ; 
		map_PandoraPDG[part_id_2f] = particleMap[ pfparticle->Daughters()[j] ] -> PdgCode() ;
		map_RecoContained[part_id_2f] = IsContained( e, trackHandle, showerHandle, findTracks, part_id_2f ) ;
		StoreInformation( e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_2f ) ;
	
		for( int j3 = 0 ; j3 < particleMap[particleMap[ pfparticle->Daughters()[j] ]->Daughters()[j2]] -> NumDaughters() ; ++j3 ) {
		  //daugheter secondary particles
		  int part_id_3f = particleMap[ particleMap[pfparticle->Daughters()[j] ]->Daughters()[j2]]->Daughters()[j3];
		  int part_MCID_3f = FindBestMCID(e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_3f ) ;
		  map_MCID_RecoID[part_id_3f] = part_MCID_3f ; 
		  map_RecoDaughters[part_id_2f].push_back( part_id_3f ) ;
		  map_RecoHiearchy[part_id_3f] = 3 ; 
		  map_PandoraPDG[part_id_3f] = particleMap[ particleMap[pfparticle->Daughters()[j] ]->Daughters()[j2]] ->PdgCode() ; 
		  map_RecoContained[part_id_3f] = IsContained( e, trackHandle, showerHandle, findTracks, part_id_3f ) ;
		  StoreInformation( e, trackHandle, showerHandle, findTracks, ShowerMothers, part_id_3f ) ;
		}
	      }	
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

  // Loop the map to get the pdg of the reconstructed particles. This is per event  
  // it->second returns the number of tracks/showers reconstructed. Default 0 
  std::map<int,int>::iterator it ;
  for( it = map_IsReconstructed.begin(); it != map_IsReconstructed.end() ; ++it ){ 
    if ( TMath::Abs( mapMC_reco_pdg[ it -> first] ) == 13   && ( it -> second ) != 0 ) ++reco_track_mu ;
    if ( TMath::Abs( mapMC_reco_pdg[ it -> first] ) == 211  && ( it -> second ) != 0 ) ++reco_track_pi ;
    if ( mapMC_reco_pdg[ it -> first ] == 2212 && ( it -> second ) != 0 ) ++reco_track_p  ;
    //    if ( mapMC_reco_pdg[ it -> first ] == 13   ) reco_track_mu += ( it -> second ) ;
    //    if ( mapMC_reco_pdg[ it -> first ] == 211  ) reco_track_pi += ( it -> second ) ;
    //    if ( mapMC_reco_pdg[ it -> first ] == 2212 ) reco_track_p  += ( it -> second ) ;    
  }

  for( it = map_IsReconstructed.begin(); it != map_IsReconstructed.end() ; ++it ){ 
    // reconstructed in tpc 
    if ( TMath::Abs(mapMC_reco_pdg[ it -> first ]) == 13 && ( it -> second ) != 0 ) h_MCLength_mu_TPC_signal -> Fill( mapTLength[ it -> first ] );
    if ( TMath::Abs(mapMC_reco_pdg[ it -> first ]) == 211  && ( it -> second ) != 0 ) h_MCLength_pi_TPC_signal -> Fill( mapTLength[ it -> first ] );
    if ( mapMC_reco_pdg[ it -> first ] == 2212 && ( it -> second ) != 0 ) h_MCLength_p_TPC_signal -> Fill( mapTLength[ it -> first ] );
    // missed signal in tpc
    if ( TMath::Abs(mapMC_reco_pdg[ map_MCID_RecoID[it -> first] ]) == 13   && ( it -> second ) == 0 ) h_MCLength_mu_TPC_miss -> Fill( mapTLength[ map_MCID_RecoID[it -> first] ] );
    if ( TMath::Abs(mapMC_reco_pdg[ map_MCID_RecoID[it -> first] ]) == 211  && ( it -> second ) == 0 ) h_MCLength_pi_TPC_miss -> Fill( mapTLength[ map_MCID_RecoID[it -> first] ] );
    if ( mapMC_reco_pdg[ map_MCID_RecoID[it -> first] ] == 2212 && ( it -> second ) == 0 ) h_MCLength_p_TPC_miss  -> Fill( mapTLength[ map_MCID_RecoID[it -> first] ] );
  }
  // READING MAPS TO STUDY HIEARCHY OF FINAL STATE 
  for( it = map_RecoHiearchy.begin(); it != map_RecoHiearchy.end() ; ++it ) {
    if( it -> second == 1 && TMath::Abs(mapMC_reco_pdg[ map_MCID_RecoID[it -> first] ]) == 13 // check if primary and truth pdg code
	&& map_RecoContained[ it -> first ][0] == 1 && map_RecoContained[ it -> first ][1] == 1 ){ // check if contained 
      if( map_RecoDaughters.find( it -> first ) != map_RecoDaughters.end() ) { 
	if( map_RecoDaughters[it->first].size() > 2) h_recoDaughters_mu -> Fill( 3 ) ; // 3 means more than 2. 
	else h_recoDaughters_mu -> Fill( map_RecoDaughters[it->first].size() ) ;
      }
      else h_recoDaughters_mu -> Fill( 0 ) ;
     
      // Loop over ID secondaries
      for( unsigned int i2 = 0 ; i2 < map_RecoDaughters[it->first].size() ; ++i2 ){
	if( map_RecoDaughters.find(map_RecoDaughters[it->first][i2]) != map_RecoDaughters.end() ){
	  h_reco3Daughters_mu -> Fill( map_RecoDaughters[map_RecoDaughters[it->first][i2] ].size() ) ;}
        else { h_reco3Daughters_mu -> Fill(0);}
      }
    } 
    if( it -> second == 1 && TMath::Abs(mapMC_reco_pdg[ map_MCID_RecoID[it -> first] ]) == 211 
	&& map_RecoContained[ it -> first ][0] == 1 && map_RecoContained[ it -> first ][1] == 1 ){
      if( map_RecoDaughters.find( it -> first ) != map_RecoDaughters.end() ) { 
	if( map_RecoDaughters[it->first].size() > 2) h_recoDaughters_pi -> Fill( 3 ) ; // 3 means more than 2. 
	else h_recoDaughters_pi -> Fill( map_RecoDaughters[it->first].size() ) ;
      }
      else h_recoDaughters_pi -> Fill( 0 ) ;
      // Loop over ID secondaries
      for( unsigned int i2 = 0 ; i2 < map_RecoDaughters[it->first].size() ; ++i2 ){
	if( map_RecoDaughters.find(map_RecoDaughters[it->first][i2]) != map_RecoDaughters.end() ){
	  h_reco3Daughters_pi -> Fill( map_RecoDaughters[map_RecoDaughters[it->first][i2] ].size() ) ;}
        else { h_reco3Daughters_pi -> Fill(0);}
      }
    } 

    if( it -> second == 1 && mapMC_reco_pdg[ map_MCID_RecoID[it -> first] ] == 2212 
	&& map_RecoContained[ it -> first ][0] == 1 && map_RecoContained[ it -> first ][1] == 1 ){
      if( map_RecoDaughters.find( it -> first ) != map_RecoDaughters.end() ) { 
	if( map_RecoDaughters[it->first].size() > 2) h_recoDaughters_p -> Fill( 3 ) ; // 3 means more than 2. 
	else h_recoDaughters_p -> Fill( map_RecoDaughters[it->first].size() ) ;
      }
      else h_recoDaughters_p -> Fill( 0 ) ;
      // Loop over ID secondaries
      for( unsigned int i2 = 0 ; i2 < map_RecoDaughters[it->first].size() ; ++i2 ){
	if( map_RecoDaughters.find(map_RecoDaughters[it->first][i2]) != map_RecoDaughters.end() ){
	  h_reco3Daughters_p -> Fill( map_RecoDaughters[map_RecoDaughters[it->first][i2] ].size() ) ;}
        else { h_reco3Daughters_p -> Fill(0);}
      }
    } 
  }
  
  // Filling histograms lenght and KE vs number of daughters
  std::map<int,double>::iterator itD ;
  for( itD = map_RecoLength.begin() ; itD != map_RecoLength.end() ; ++itD ) {
    if( map_RecoDaughters[itD->first].size() == 0 ) {
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 13   ) h_recoLength_mu_noD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 211  ) h_recoLength_pi_noD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 2212 ) h_recoLength_p_noD  -> Fill( itD->second ) ; 
    } else {
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 13   ) h_recoLength_mu_wD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 211  ) h_recoLength_pi_wD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 2212 ) h_recoLength_p_wD  -> Fill( itD->second ) ; 
    }
  }

  for( itD = map_RecoKEnergy.begin() ; itD != map_RecoKEnergy.end() ; ++itD ) {
    if( map_RecoDaughters[itD->first].size() == 0 ) {
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 13   ) h_recoKE_mu_noD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 211  ) h_recoKE_pi_noD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 2212 ) h_recoKE_p_noD  -> Fill( itD->second ) ; 
    } else {
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 13   ) h_recoKE_mu_wD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 211  ) h_recoKE_pi_wD -> Fill( itD->second ) ; 
      if( mapMC_reco_pdg[ map_MCID_RecoID[ itD->first ] ] == 2212 ) h_recoKE_p_wD  -> Fill( itD->second ) ; 
    }
  }
  
 
  mcparticle_tree -> Fill();
  recoevent_tree -> Fill();

}

void test::NeutrinoTopologyAnalyzer::StoreInformation( 
      art::Event const & e, art::Handle< std::vector< recob::Track > > const & trackHandle, 
      art::Handle< std::vector< recob::Shower > > const & showerHandle, art::FindManyP< recob::Track > const & findTracks, 
      std::map< int , std::vector< int > > & ShowerMothers, int const & part_id_f ){//, int const & primary_daughter) {
  /*
    1) Find muon/pion candidates
    2) Efficiency simple method
    3) Apply topology consideretions -> New efficiency
    4) Number of kinks -> Efficiency
    5) Michael electrons -> Efficiency
  */
  
  // Save track info
  if ( findTracks.at( part_id_f ).size() != 0 ){
    std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(part_id_f);
    art::FindManyP< recob::Hit > findHits (  trackHandle, e, RecoTrackLabel ) ;
    art::FindManyP< anab::Calorimetry > findCalorimetry ( trackHandle, e, RecoCaloLabel );
    art::FindManyP< anab::ParticleID > findPID ( trackHandle, e, RecoPIDLabel );

    // Loop over tracks found for track_f
    for( unsigned int n = 0 ; n < track_f.size() ; ++n ){
      has_reco_tracks = true ; 
      // Get track based variables
      std::vector< art::Ptr<recob::Hit> > hit_f        = findHits.at(track_f[n]->ID()); 
      std::vector< art::Ptr<anab::Calorimetry> > cal_f = findCalorimetry.at(track_f[n]->ID());
      std::vector< art::Ptr<anab::ParticleID> > pid_f  = findPID.at(track_f[n]->ID());
      
      //Loop over PID associations 
      for ( unsigned int k = 0 ; k < pid_f.size() ; ++k ){
	
	if( !pid_f[k] ) continue ;
	if( !pid_f[k]->PlaneID().isValid) continue ;
	if( pid_f[k]->PlaneID().Plane != 2 ) continue ; // only look at collection plane for dEdx information

	// looks for muons /pions or daughters of muons /pions candidates
	// only checks for primary particle. If it's not a candidate we don't store the information for the primary neither the daughters
	if( is_candidate == false ) is_candidate = IsMuonPionCandidateChi2( pid_f[k] ) ; 
	
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

	    // save the most common answer in tr_id_best: 
	    if( tr_id_energy == tr_id_charge && tr_id_energy == tr_id_hits ) tr_id_best = tr_id_energy ;
	    if( tr_id_energy == tr_id_charge && tr_id_energy != tr_id_hits ) tr_id_best = tr_id_energy ;
	    if( tr_id_energy != tr_id_charge && tr_id_energy == tr_id_hits ) tr_id_best = tr_id_energy ;
	    if( tr_id_energy != tr_id_charge && tr_id_charge == tr_id_hits ) tr_id_best = tr_id_charge ;
	    if( tr_id_energy != tr_id_charge && tr_id_energy != tr_id_hits && tr_id_charge != tr_id_hits) {
	      if( tr_id_hits > 0 ) tr_id_best = tr_id_hits ;
	      else if( tr_id_energy > 0 ) tr_id_best = tr_id_energy ;
	      else if( tr_id_charge > 0 ) tr_id_best = tr_id_hits ;
	      else tr_id_best = -9999 ; // couldn't find id 
	    }

	    IsReconstructed( tr_id_best ) ; // check if MC particle is reconstructed and stores it in a map 

	    if( TMath::Abs(mapMC_reco_pdg[ tr_id_best ]) == 13   ) Chi2p_Tmu -> Fill( pid_f[k] ->Chi2Proton() ) ;
	    if( TMath::Abs(mapMC_reco_pdg[ tr_id_best ]) == 221  ) Chi2p_Tpi -> Fill( pid_f[k] ->Chi2Proton() ) ;
	    if( mapMC_reco_pdg[ tr_id_best ] == 2212 ) {
	      Chi2p_Tp -> Fill( pid_f[k] -> Chi2Proton() ) ;
	      ++total_reco_p ;
	      if( IsMuonPionCandidateChi2Proton( pid_f[k] ) == false )         ++reco_p ; 
	      if( pid_f[k] -> Chi2Proton() < 100 && pid_f[k] -> Chi2Proton() > 70 
		  && pid_f[k] -> Chi2Proton() < 2 * pid_f[k] -> Chi2Muon() 
		  && pid_f[k] -> Chi2Proton() < 2 * pid_f[k] -> Chi2Pion() )   ++reco_p ;
	    }

	    if( is_candidate == true ) {
	      // If the same track is broken into different ones, we add the information. 
	      // We do that only if Pandora knows about it. Not adding truth information here
	      map_RecoHits[ part_id_f ] += track_f[n]->LastValidPoint() + 1 ;
	      map_RecoLength[ part_id_f ] += track_f[n]->Length() ;
	      map_RecoKEnergy[ part_id_f ] += cal_f[m]->KineticEnergy(); 
	      for( unsigned int l = 0 ; l < track_f[n]->LastValidPoint() + 1 ; ++l ) {
		map_RecoXPosition[ part_id_f ].push_back( track_f[n]->TrajectoryPoint( l ).position.X() ) ; 
		map_RecoYPosition[ part_id_f ].push_back( track_f[n]->TrajectoryPoint( l ).position.Y() ) ; 
		map_RecoZPosition[ part_id_f ].push_back( track_f[n]->TrajectoryPoint( l ).position.Z() ) ; 
	      }
	      EfficiencyCalo( pid_f[k] , mapMC_reco_pdg[ tr_id_best ], "muon" ); 
	    }
	  }// just collection plane 
	  // calo information is stored in all planes. Need to read dEdx in an ordered way           
	  if( is_candidate == false ) continue ; 
	  
	  for( unsigned int l = 0 ; l < (cal_f[m]->XYZ()).size() ; ++l ) {
	    for( unsigned int t = 0 ; t < track_f[n] -> LastValidPoint() ; ++t ){
	      if( cal_f[m]->XYZ()[l].X() == map_RecoXPosition[ part_id_f ][t] 
		  && cal_f[m]->XYZ()[l].Y() == map_RecoYPosition[ part_id_f ][t] 
		  && cal_f[m]->XYZ()[l].Z() == map_RecoZPosition[ part_id_f ][t]  ){
		map_RecodEdx[ part_id_f ].push_back( cal_f[m]->dEdx()[l] ) ; 
	      }
	    }
	  }
	} //close calo

	if( is_candidate == true ) {
	  if( TMath::Abs(mapMC_reco_pdg[ tr_id_best ]) == 13   ) {
	    h_recoLength_mu -> Fill( map_RecoLength[ part_id_f ] ) ;
	    h_recoKE_mu -> Fill( map_RecoKEnergy[ part_id_f] ) ;
	  }
	  if( TMath::Abs(mapMC_reco_pdg[ tr_id_best ]) == 211  ) {
	    h_recoLength_pi -> Fill( map_RecoLength[ part_id_f ] ) ;
	    h_recoKE_pi -> Fill( map_RecoKEnergy[ part_id_f] ) ;
	  }
	  if( mapMC_reco_pdg[ tr_id_best ] == 2212 ){
	    h_recoLength_p  -> Fill( map_RecoLength[ part_id_f ] ) ;
	    h_recoKE_p -> Fill( map_RecoKEnergy[ part_id_f] ) ;
	  }
	}
      } //close pid
    } //close track
    
  } else if( showerHandle.isValid() && showerHandle->size() != 0 ) { // if no track look into showers 
    // Need to call with id! 
    has_reco_showers = true ; 
    art::FindManyP< recob::Hit > findHitShower( showerHandle, e, RecoShowerLabel ) ;
    art::FindManyP< recob::SpacePoint > findSpacePoint( showerHandle, e, RecoShowerLabel ) ;
    for( unsigned int y = 0 ; y < showerHandle->size() ; ++y ) {
      art::Ptr< recob::Shower > shower_f( showerHandle, y ) ;
      std::vector< art::Ptr<recob::Hit> > hit_sh_f = findHitShower.at(y) ; 
      std::vector< art::Ptr<recob::SpacePoint> > spacepoint_f = findSpacePoint.at(y) ;
      if( spacepoint_f.size() == 0 ) continue ; 

      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain( ShowerMothers, hit_sh_f , shower_f->best_plane() ) ;
      tr_id_best = ShowerTrackInfo.first ;
      if( !ShowerTrackInfo.first ) tr_id_best = -9999 ; // default no match 
     
      map_RecoHits[ part_id_f ] += spacepoint_f.size() ; 
      
      for( unsigned int l = 0 ; l < spacepoint_f.size() ; ++l ) {
	  map_RecoXPosition[ part_id_f ].push_back( spacepoint_f[l]->XYZ()[0] ) ; 
	  map_RecoYPosition[ part_id_f ].push_back( spacepoint_f[l]->XYZ()[1] ) ; 
	  map_RecoZPosition[ part_id_f ].push_back( spacepoint_f[l]->XYZ()[2] ) ;
	}
	          
      if( shower_f->has_length() ) { 
	map_RecoLength[ part_id_f ] += shower_f->Length(); //primary_daughter 
      } else  {
	map_RecoLength[ part_id_f ]  = pow(spacepoint_f[spacepoint_f.size()-1]->XYZ()[0] - spacepoint_f[0]->XYZ()[0], 2 ) ;
	map_RecoLength[ part_id_f ] += pow(spacepoint_f[spacepoint_f.size()-1]->XYZ()[1] - spacepoint_f[0]->XYZ()[1], 2 ) ;
	map_RecoLength[ part_id_f ] += pow(spacepoint_f[spacepoint_f.size()-1]->XYZ()[2] - spacepoint_f[0]->XYZ()[2], 2 ) ;
	map_RecoLength[ part_id_f ]  = sqrt( map_RecoLength[ part_id_f ] ) ;
      }

      IsReconstructed( tr_id_best ) ; // check if MC particle is reconstructed and stores it in a map 

    }
  } // track vs shower
}

int test::NeutrinoTopologyAnalyzer::FindBestMCID(art::Event const & e, art::Handle< std::vector< recob::Track > > const & trackHandle, 
      art::Handle< std::vector< recob::Shower > > const & showerHandle, art::FindManyP< recob::Track > const & findTracks, 
      std::map< int , std::vector< int > > & ShowerMothers, int const & part_id_f ) {
  if ( findTracks.at( part_id_f ).size() != 0 ){
    std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(part_id_f);
    art::FindManyP< recob::Hit > findHits (  trackHandle, e, RecoTrackLabel ) ;
    // Loop over tracks found for track_f
    for( unsigned int n = 0 ; n < track_f.size() ; ++n ){
        // Get track based variables
      std::vector< art::Ptr<recob::Hit> > hit_f        = findHits.at(track_f[n]->ID()); 
      // Get associated MCParticle ID using 3 different methods:
      //    Which particle contributes the most energy to all the hits
      //    Which particle contributes the reco charge to all the hits
      //    Which particle is the biggest contributor to all the hits
      tr_id_energy      = RecoUtils::TrueParticleIDFromTotalTrueEnergy(hit_f);
      tr_id_charge      = RecoUtils::TrueParticleIDFromTotalRecoCharge(hit_f);
      tr_id_hits        = RecoUtils::TrueParticleIDFromTotalRecoHits(hit_f);
      // save the most common answer in tr_id_best: 
      if( tr_id_energy == tr_id_charge && tr_id_energy == tr_id_hits ) tr_id_best = tr_id_energy ;
      if( tr_id_energy == tr_id_charge && tr_id_energy != tr_id_hits ) tr_id_best = tr_id_energy ;
      if( tr_id_energy != tr_id_charge && tr_id_energy == tr_id_hits ) tr_id_best = tr_id_energy ;
      if( tr_id_energy != tr_id_charge && tr_id_charge == tr_id_hits ) tr_id_best = tr_id_charge ;
      if( tr_id_energy != tr_id_charge && tr_id_energy != tr_id_hits && tr_id_charge != tr_id_hits) {
	if( tr_id_hits > 0 ) tr_id_best = tr_id_hits ;
	else if( tr_id_energy > 0 ) tr_id_best = tr_id_energy ;
	else if( tr_id_charge > 0 ) tr_id_best = tr_id_hits ;
	else tr_id_best = -9999 ; // couldn't find id 
      }
    }
  } else if( showerHandle.isValid() && showerHandle->size() != 0 ) { // if no track look into showers 
    art::FindManyP< recob::Hit > findHitShower( showerHandle, e, RecoShowerLabel ) ;
    art::FindManyP< recob::SpacePoint > findSpacePoint( showerHandle, e, RecoShowerLabel ) ;
    for( unsigned int y = 0 ; y < showerHandle->size() ; ++y ) {
      art::Ptr< recob::Shower > shower_f( showerHandle, y ) ;
      std::vector< art::Ptr<recob::Hit> > hit_sh_f = findHitShower.at(y) ; 
      std::vector< art::Ptr<recob::SpacePoint> > spacepoint_f = findSpacePoint.at(y) ;
      if( spacepoint_f.size() == 0 ) continue ; 
      
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain( ShowerMothers, hit_sh_f , shower_f->best_plane() ) ;
      tr_id_best = ShowerTrackInfo.first ;
      if( !ShowerTrackInfo.first ) tr_id_best = -9999 ; // default no match 
    }
  }// end shower
  return tr_id_best ;
}


bool test::NeutrinoTopologyAnalyzer::IsMuonChi2( art::Ptr<anab::ParticleID> const & pid_f ) {
  if( pid_f -> Chi2Proton() < 100 && pid_f -> Chi2Proton() > 70  && pid_f -> Chi2Proton() > 2 * pid_f -> Chi2Muon() ) return true ; 
  else if( IsMuonPionCandidateChi2Proton( pid_f ) == false ) return false ;
  else if( pid_f ->Chi2Muon() < pid_f ->Chi2Proton() && pid_f ->Chi2Muon() < pid_f ->Chi2Kaon() && pid_f ->Chi2Muon() < pid_f -> Chi2Pion() ) return true ; 
  return false ;
}

bool test::NeutrinoTopologyAnalyzer::IsPionChi2( art::Ptr<anab::ParticleID> const & pid_f ) {
  if( pid_f -> Chi2Proton() < 100 && pid_f -> Chi2Proton() > 70  && pid_f -> Chi2Proton() > 2 * pid_f -> Chi2Pion() ) return true ; 
  else if( IsMuonPionCandidateChi2Proton( pid_f ) == false ) return false ;
  else if( pid_f ->Chi2Pion() < pid_f ->Chi2Proton() && pid_f ->Chi2Pion() < pid_f ->Chi2Kaon() && pid_f ->Chi2Pion() < pid_f -> Chi2Muon() ) return true ; 
  return false ;
}

bool test::NeutrinoTopologyAnalyzer::IsMuonPionCandidateChi2( art::Ptr<anab::ParticleID> const & pid_f )
{

  if( IsMuonChi2( pid_f ) == true || IsPionChi2( pid_f ) == true ) return true ;
  return false ; 
}

bool test::NeutrinoTopologyAnalyzer::IsMuonPionCandidateChi2Proton( art::Ptr<anab::ParticleID> const & pid_f )
{
  if( pid_f->Chi2Proton() < 100 ) return false ; // R.J : 80  
  // If this cut is increased, many muons will be missreconstructed 
  // Keep the information for now 
	 
  return true ; 
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

void test::NeutrinoTopologyAnalyzer::IsReconstructed( int  const & best_id )  {
  std::map<int, int>::iterator it;
  for (it = map_IsReconstructed.begin(); it != map_IsReconstructed.end(); it++) {
    if ( it -> first == best_id ) ++map_IsReconstructed[best_id] ; // will return the number of tracks reconstructed! 
  }
}


double test::NeutrinoTopologyAnalyzer::EfficiencyCalo( art::Ptr<anab::ParticleID> const & pid_f , int const & true_pdg , std::string const & particle ) {
  // True information -> Calculated in the MC loop 
  if( IsMuonChi2( pid_f ) == true ) {
    ++reco_primary_mu ;
    if( TMath::Abs(true_pdg) == 13 ) { ++signal_mu ;
    } else if( TMath::Abs(true_pdg) == 211  ) { ++bg_mu_pi ;
    } else if( true_pdg == 2212 ) { ++bg_mu_p ;
    } else { ++bg_mu_others ; }
        
  } else if( IsPionChi2( pid_f ) == true ) {
    ++reco_primary_pi ;
    if( TMath::Abs(true_pdg) == 211 ) { ++signal_pi ;
    } else if( TMath::Abs(true_pdg) == 13  ) { ++bg_pi_mu ; // check for miss reco particles . was reco as pion is a muon
    } else if( true_pdg == 2212 ) { ++bg_pi_p ;
    } else { ++bg_mu_others ; }

  }

  eff_mu = signal_mu / (double)true_primary_mu ; 
  purity_mu = signal_mu / (double)reco_primary_mu ; 
  eff_pi = signal_pi / (double)true_primary_pi ; 
  purity_pi = signal_pi / (double)reco_primary_pi ; 

  if( particle == "pion" ) return eff_pi ; 
  if( particle == "purity pion" ) return purity_pi ; 
  if( particle == "purity muon" ) return purity_mu ; 
  return eff_mu;
}

std::vector< bool > test::NeutrinoTopologyAnalyzer::MCIsContained( simb::MCParticle const & trueParticle ) {
  // Checks if true track is contained in fiducial. Not looking for showers as we are interested in muon/pion candidates. Truth info does not contain showers
  bool vertex_contained = true , end_contained = true ;
  std::vector< bool > contained_info ; 

  if( ( trueParticle.Trajectory().X( 0 ) > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
      || ( trueParticle.Trajectory().X( 0 ) < (-CoordinateOffSetX + SelectedBorderX)) 
      || ( trueParticle.Trajectory().Y( 0 ) > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
      || ( trueParticle.Trajectory().Y( 0 ) < (-CoordinateOffSetY + SelectedBorderY)) 
      || ( trueParticle.Trajectory().Z( 0 ) > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
      || ( trueParticle.Trajectory().Z( 0 ) < (-CoordinateOffSetZ + SelectedBorderZ))) vertex_contained = false ;
  if( ( trueParticle.EndX() > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
      || ( trueParticle.EndX() < (-CoordinateOffSetX + SelectedBorderX)) 
      || ( trueParticle.EndY() > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
      || ( trueParticle.EndY() < (-CoordinateOffSetY + SelectedBorderY)) 
      || ( trueParticle.EndZ() > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
      || ( trueParticle.EndZ() < (-CoordinateOffSetZ + SelectedBorderZ))) end_contained = false ;
  contained_info.push_back( vertex_contained ) ; 
  contained_info.push_back( end_contained ) ; 
  return contained_info ; 
}


std::vector< bool > test::NeutrinoTopologyAnalyzer::IsContained( art::Event const & e, art::Handle< std::vector< recob::Track > > & trackHandle, art::Handle< std::vector< recob::Shower > > & showerHandle, art::FindManyP< recob::Track > & findTracks , int & part_id_f ) {
  bool vertex_contained = true , end_contained = true ;
  std::vector< bool > contained_info ; 

  if ( findTracks.at( part_id_f ).size()!=0 ){
    std::vector< art::Ptr<recob::Track> > track_f = findTracks.at(part_id_f);
    art::FindManyP< anab::ParticleID > findPID ( trackHandle, e, RecoPIDLabel );
    
    for( unsigned int n = 0 ; n < track_f.size() ; ++n ){      
      if( ( track_f[n]->Start().X() > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
	  || ( track_f[n]->Start().X() < (-CoordinateOffSetX + SelectedBorderX)) 
	  || ( track_f[n]->Start().Y() > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
	  || ( track_f[n]->Start().Y() < (-CoordinateOffSetY + SelectedBorderY)) 
	  || ( track_f[n]->Start().Z() > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
	  || ( track_f[n]->Start().Z() < (-CoordinateOffSetZ + SelectedBorderZ))) vertex_contained = false ;
      if( ( track_f[n]->End().X() > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
	  || ( track_f[n]->End().X() < (-CoordinateOffSetX + SelectedBorderX)) 
	  || ( track_f[n]->End().Y() > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
	  || ( track_f[n]->End().Y() < (-CoordinateOffSetY + SelectedBorderY)) 
	  || ( track_f[n]->End().Z() > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
	  || ( track_f[n]->End().Z() < (-CoordinateOffSetZ + SelectedBorderZ))) end_contained = false ;
      
    } //close track  	  
  } else if( showerHandle.isValid() && showerHandle->size() ) { // if no track look in showers 
    art::FindManyP< recob::SpacePoint > findSpacePoint( showerHandle, e, RecoShowerLabel ) ;

    for( unsigned int y = 0 ; y < showerHandle->size() ; ++y ) {
      art::Ptr< recob::Shower > shower_f( showerHandle, y ) ;
      std::vector< art::Ptr<recob::SpacePoint> > spacepoint_f = findSpacePoint.at(y) ;
      if( spacepoint_f.size() == 0 ) continue ;
      if( ( spacepoint_f[0]->XYZ()[0] > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
	  || ( spacepoint_f[0]->XYZ()[0] < (-CoordinateOffSetX + SelectedBorderX)) 
	  || ( spacepoint_f[0]->XYZ()[1] > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
	  || ( spacepoint_f[0]->XYZ()[1] < (-CoordinateOffSetY + SelectedBorderY)) 
	  || ( spacepoint_f[0]->XYZ()[2] > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
	  || ( spacepoint_f[0]->XYZ()[2] < (-CoordinateOffSetZ + SelectedBorderZ))) vertex_contained = false ;
      if( ( spacepoint_f[0]->XYZ()[ spacepoint_f.size()-1 ] > (DetectorHalfLengthX - CoordinateOffSetX - SelectedBorderX)) 
	  || ( spacepoint_f[spacepoint_f.size()-1]->XYZ()[0] < (-CoordinateOffSetX + SelectedBorderX)) 
	  || ( spacepoint_f[spacepoint_f.size()-1]->XYZ()[1] > (DetectorHalfLengthY - CoordinateOffSetY - SelectedBorderY)) 
	  || ( spacepoint_f[spacepoint_f.size()-1]->XYZ()[1] < (-CoordinateOffSetY + SelectedBorderY)) 
	  || ( spacepoint_f[spacepoint_f.size()-1]->XYZ()[2] > (DetectorHalfLengthZ - CoordinateOffSetZ - SelectedBorderZ)) 
	  || ( spacepoint_f[spacepoint_f.size()-1]->XYZ()[2] < (-CoordinateOffSetZ + SelectedBorderZ))) end_contained = false ;
    }	
  } 
  contained_info.push_back( vertex_contained ) ; 
  contained_info.push_back( end_contained ) ; 
  return contained_info ; 
}




void test::NeutrinoTopologyAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of required member function here.
  // Here you add an external fcl file to change configuration
}


void test::NeutrinoTopologyAnalyzer::beginJob( )
{
  // don't initialize per event. 

  true_mu = 0 ;
  true_pi = 0 ;
  true_p = 0 ; 
  reco_primary_mu = 0 ;
  reco_primary_pi = 0 ;
  true_primary_mu = 0 ;
  true_primary_pi = 0 ;
  signal_mu = 0 ;
  signal_pi = 0 ; 
  bg_mu_pi = 0 ;
  bg_mu_p = 0 ;
  bg_mu_others = 0 ; 
  bg_pi_mu = 0 ;
  bg_pi_p = 0 ;
  bg_pi_others = 0 ; 
  eff_mu = 0 ;
  eff_pi = 0 ;
  purity_mu = 0 ;
  purity_pi = 0 ;
  total_reco_p = 0 ; // number of reconstructed particles (p)
  reco_p =0 ; 
  reco_track_mu = 0 ;
  reco_track_pi = 0 ;
  reco_track_p = 0 ;

  
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
  mcparticle_tree -> Branch( "mapTPiDaughterPdg",       "std::map<int,int>",   &mapTPiDaughterPdg);

  // Reco tree
  recoevent_tree -> Branch( "event_id",                 &event_id,            "event_id/I");
  recoevent_tree -> Branch( "nu_reconstructed",         &nu_reconstructed,    "nu_reconstructed/B");
  recoevent_tree -> Branch( "nu_rvertex_contained",     &nu_rvertex_contained,"nu_rvertex_contained/B");
  recoevent_tree -> Branch( "numb_nu_reco",             &numb_nu_reco,        "numb_nu_reco/I");
  recoevent_tree -> Branch( "nu_reco_vertex",           &nu_reco_vertex,      "nu_reco_vertex[3]/D");
  recoevent_tree -> Branch( "error_vertex_reco",        &error_vertex_reco,   "error_vertex_reco/D");
  recoevent_tree -> Branch( "map_MCID_RecoID",          "std::map<int,int>",  &map_MCID_RecoID);
  recoevent_tree -> Branch( "mapMC_reco_pdg",           "std::map<int,int>",  &mapMC_reco_pdg);
  recoevent_tree -> Branch( "map_RecoContained",        "std::map<int,std::vector<bool>>",  &map_RecoContained);
  recoevent_tree -> Branch( "map_RecoHits",             "std::map<int,int>",  &map_RecoHits);
  recoevent_tree -> Branch( "map_RecoPrimary",          "std::map<int,int>",  &map_RecoPrimary);
  recoevent_tree -> Branch( "map_RecoDaughters",        "std::map<int,std::vector<int>>",  &map_RecoDaughters);
  recoevent_tree -> Branch( "map_PandoraPDG",           "std::map<int,int>",  &map_PandoraPDG);
  recoevent_tree -> Branch( "map_RecoLength",           "std::map<int,double>",  &map_RecoLength);
  recoevent_tree -> Branch( "map_RecoKEnergy",          "std::map<int,double>",  &map_RecoKEnergy);
  recoevent_tree -> Branch( "map_RecoXPosition",        "std::map<int,std::vector<double>>",  &map_RecoXPosition);
  recoevent_tree -> Branch( "map_RecoYPosition",        "std::map<int,std::vector<double>>",  &map_RecoYPosition);
  recoevent_tree -> Branch( "map_RecoZPosition",        "std::map<int,std::vector<double>>",  &map_RecoZPosition);
  recoevent_tree -> Branch( "map_RecodEdx",             "std::map<int,std::vector<double>>",  &map_RecodEdx);
  recoevent_tree -> Branch( "map_RecoHiearchy",        "std::map<int,int>",  &map_RecoHiearchy);
   
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
  recoevent_tree  ->Write();
  
  file.Write();
  file.Close();

  std::ofstream eff_chi2_file ; 
  
  eff_chi2_file.open("eff_chi2_information.txt");

  if( eff_chi2_file.is_open() ) {
    eff_chi2_file << " Primary particles study " << std::endl;
    eff_chi2_file << " ******************** MUONS CALORIMETRY STUDY *********************** \n" ;
    eff_chi2_file << " Number of reconstructed muons = " << reco_primary_mu << "\n" ;
    eff_chi2_file << " Number of true muons from MC = " << true_primary_mu << "\n" ;
    eff_chi2_file << " Signal events ( true & reco ) = " << signal_mu << "\n" ;
    eff_chi2_file << " Background events : \n " ;
    eff_chi2_file << "     -- > True MC ID - pions " << bg_mu_pi << "\n" ; 
    eff_chi2_file << "     -- > True MC ID - proton " << bg_mu_p << "\n" ; 
    eff_chi2_file << "     -- > True MC ID - other " << bg_mu_others << "\n" ; 
    eff_chi2_file << " EFFICIENCY =  " << eff_mu << "\n" ;
    eff_chi2_file << " PURITY =  " << purity_mu << "\n" ; 
    eff_chi2_file << " \n\n******************** PIONS CALORIMETRY STUDY *********************** \n" ;
    eff_chi2_file << " Number of reconstructed pions = " << reco_primary_pi << "\n" ;
    eff_chi2_file << " Number of true pions from MC = " << true_primary_pi << "\n" ;
    eff_chi2_file << " Signal events ( true & reco ) = " << signal_pi << "\n" ;
    eff_chi2_file << " Background events : \n " ;
    eff_chi2_file << "     -- > True MC ID - muons " << bg_pi_mu << "\n" ; 
    eff_chi2_file << "     -- > True MC ID - proton " << bg_pi_p << "\n" ; 
    eff_chi2_file << "     -- > True MC ID - other " << bg_pi_others << "\n" ; 
    eff_chi2_file << " EFFICIENCY =  " << eff_pi << "\n" ;
    eff_chi2_file << " PURITY =  " << purity_pi << "\n" ; 
    eff_chi2_file << " \n\n******************** RECONSTRUCT PROTONS STUDY *********************** \n" ;
    eff_chi2_file << " Number of reconstructed protons = " << total_reco_p << "\n" ;
    eff_chi2_file << " Number of reco p reconstructed as protons = " << reco_p << "\n" ;
    eff_chi2_file << " % reconstructed  = " << (double) reco_p / (double) total_reco_p * 100 << "\n" ;
    
  }


  if( eff_chi2_file.is_open() ) {
    eff_chi2_file << " \n\n******************** STUDY PANDORA RECONSTRUCTION ************************************* \n" ;
    eff_chi2_file << " * Looking at if particle is reconstructed, I don't care if we reconstruct them right  * \n" ;
    eff_chi2_file << " *************************************************************************************** \n" ;
    eff_chi2_file << "  # True muon leaves signal in TPC   : " << reco_track_mu << "   vs   #MC muons   = " << true_mu << " \n" ;  
    eff_chi2_file << "  # True pion leaves signal in TPC   : " << reco_track_pi << "   vs   #MC pions   = " << true_pi << " \n" ;  
    eff_chi2_file << "  # True proton leaves signal in TPC : " << reco_track_p  << "   vs   #MC protons = " << true_p  << " \n" ;  
  }

  eff_chi2_file.close();


  // PRINT HISTOGRAMS
  //**********************************************************************//  
  Chi2p_Tmu -> SetStats(0);
  Chi2p_Tmu ->GetXaxis()->SetTitle(" chi2 under proton hypotesis " ) ;
  Chi2p_Tmu ->GetYaxis()->SetTitle(" Events " ) ;
  
  Chi2p_Tp -> SetStats(0);
  Chi2p_Tp ->GetXaxis()->SetTitle(" chi2 under proton hypotesis " ) ;
  Chi2p_Tp ->GetYaxis()->SetTitle(" Events " ) ;
  
  l->SetBorderSize(0);
  l->AddEntry( Chi2p_Tmu,      " True Muon",    "l" );
  l->AddEntry( Chi2p_Tpi,      " True Pion",    "l" );
  l->AddEntry( Chi2p_Tp,       " True Proton ",           "l" );
    
  Chi2p_Tmu->SetLineColor(2);
  Chi2p_Tpi->SetLineColor(4);
  
  Chi2p_Tmu->Draw();
  Chi2p_Tpi->Draw("same");
  Chi2p_Tp->Draw("same");
  l->Draw();

  c->SaveAs("chi2_test.root");
  c->Clear();

  TFinalStatePi -> SetStats(0); 
  TFinalStatePi -> GetXaxis() -> SetTitle(" pdg rescatter" ) ; 
  TFinalStatePi -> Draw();

  c->SaveAs("pion_finalstate.root") ;
  c->Clear();
  //**********************************************************************//  
  h_recoLength_mu -> SetStats( 0 ) ;
  h_recoLength_pi -> SetStats( 0 ) ;
  h_recoLength_p -> SetStats( 0 ) ; 
  
  h_recoLength_mu -> SetLineColor( 1 ) ;
  h_recoLength_pi -> SetLineColor( 2 ) ;
  h_recoLength_p -> SetLineColor( 4 ) ;
  
  h_recoLength_mu -> GetXaxis() ->SetTitle( "Length[cm]" ) ;
  
  h_recoLength_mu -> Draw( ) ;
  h_recoLength_pi -> Draw( "same" ) ;
  h_recoLength_p -> Draw( "same" ) ;

  c->SaveAs("Length_recoCandidates.root");
  c->Clear();
  //**********************************************************************//  
  h_recoKE_mu -> SetStats( 0 ) ;
  h_recoKE_pi -> SetStats( 0 ) ;
  h_recoKE_p -> SetStats( 0 ) ; 
  
  h_recoKE_mu -> SetLineColor( 1 ) ;
  h_recoKE_pi -> SetLineColor( 2 ) ;
  h_recoKE_p -> SetLineColor( 4 ) ;
  
  h_recoKE_mu -> GetXaxis() ->SetTitle( "KE[GeV]" ) ;
  
  h_recoKE_mu -> Draw( ) ;
  h_recoKE_pi -> Draw( "same" ) ;
  h_recoKE_p -> Draw( "same" ) ;

  c->SaveAs("KE_recoCandidates.root");
  c->Clear();
  //**********************************************************************//  
  h_recoKE_mu_noD -> SetStats( 0 ) ;
  h_recoKE_pi_noD -> SetStats( 0 ) ;
  h_recoKE_p_noD -> SetStats( 0 ) ; 
  
  h_recoKE_mu_noD -> SetLineColor( 1 ) ;
  h_recoKE_pi_noD -> SetLineColor( 2 ) ;
  h_recoKE_p_noD -> SetLineColor( 4 ) ;
  
  h_recoKE_mu_noD -> GetXaxis() ->SetTitle( "KE[GeV]" ) ;
  
  h_recoKE_mu_noD -> Draw( ) ;
  h_recoKE_pi_noD -> Draw( "same" ) ;
  h_recoKE_p_noD -> Draw( "same" ) ;

  c->SaveAs("KE_recoCandidates_nodaughters.root");
  c->Clear();

  //**********************************************************************//  
  h_recoLength_mu_noD -> SetStats( 0 ) ;
  h_recoLength_pi_noD -> SetStats( 0 ) ;
  h_recoLength_p_noD -> SetStats( 0 ) ; 
  
  h_recoLength_mu_noD -> SetLineColor( 1 ) ;
  h_recoLength_pi_noD -> SetLineColor( 2 ) ;
  h_recoLength_p_noD -> SetLineColor( 4 ) ;
  
  h_recoLength_mu_noD -> GetXaxis() ->SetTitle( "Length[cm]" ) ;
  
  h_recoLength_mu_noD -> Draw( ) ;
  h_recoLength_pi_noD -> Draw( "same" ) ;
  h_recoLength_p_noD -> Draw( "same" ) ;

  c->SaveAs("Length_recoCandidates_noD.root");
  c->Clear();
  //**********************************************************************//  
  h_recoKE_mu_wD -> SetStats( 0 ) ;
  h_recoKE_pi_wD -> SetStats( 0 ) ;
  h_recoKE_p_wD -> SetStats( 0 ) ; 
  
  h_recoKE_mu_wD -> SetLineColor( 1 ) ;
  h_recoKE_pi_wD -> SetLineColor( 2 ) ;
  h_recoKE_p_wD -> SetLineColor( 4 ) ;
  
  h_recoKE_mu_wD -> GetXaxis() ->SetTitle( "KE[GeV]" ) ;
  
  h_recoKE_mu_wD -> Draw( ) ;
  h_recoKE_pi_wD -> Draw( "same" ) ;
  h_recoKE_p_wD -> Draw( "same" ) ;

  c->SaveAs("KE_recoCandidates_wD.root");
  c->Clear();

  //**********************************************************************//  
  h_recoLength_mu_wD -> SetStats( 0 ) ;
  h_recoLength_pi_wD -> SetStats( 0 ) ;
  h_recoLength_p_wD -> SetStats( 0 ) ; 
  
  h_recoLength_mu_wD -> SetLineColor( 1 ) ;
  h_recoLength_pi_wD -> SetLineColor( 2 ) ;
  h_recoLength_p_wD -> SetLineColor( 4 ) ;
  
  h_recoLength_mu_wD -> GetXaxis() ->SetTitle( "Length[cm]" ) ;
  
  h_recoLength_mu_wD -> Draw( ) ;
  h_recoLength_pi_wD -> Draw( "same" ) ;
  h_recoLength_p_wD -> Draw( "same" ) ;

  c->SaveAs("Length_recoCandidates_wD.root");
  c->Clear();
  //**********************************************************************//  
  h_MCLength_mu_TPC_signal -> SetStats(0);
  h_MCLength_mu_TPC_miss   -> SetStats(0);

  h_MCLength_mu_TPC_signal -> SetLineColor(1);
  h_MCLength_mu_TPC_miss   -> SetLineColor(2);
  
  h_MCLength_mu_TPC_signal -> GetXaxis() -> SetTitle( "Length[cm]");
  h_MCLength_mu_TPC_signal -> GetYaxis() -> SetTitle( "Entries[#]");
 
  h_MCLength_mu_TPC_miss     -> Draw();
  h_MCLength_mu_TPC_signal   -> Draw("same");
  
  c->SaveAs("MCLength_mu_signalTPC.root") ;
  c->Clear();
  //**********************************************************************//  
  h_MCLength_pi_TPC_signal -> SetStats(0);
  h_MCLength_pi_TPC_miss   -> SetStats(0);

  h_MCLength_pi_TPC_signal -> SetLineColor(1);
  h_MCLength_pi_TPC_miss   -> SetLineColor(2);
  
  h_MCLength_pi_TPC_signal -> GetXaxis() -> SetTitle( "Length[cm]");
  h_MCLength_pi_TPC_signal -> GetYaxis() -> SetTitle( "Entries[#]");
  
  h_MCLength_pi_TPC_miss -> Draw();
  h_MCLength_pi_TPC_signal -> Draw("same");
  
  c->SaveAs("MCLength_pi_signalTPC.root") ;
  c->Clear();
  //**********************************************************************//  
  h_MCLength_p_TPC_signal -> SetStats(0);
  h_MCLength_p_TPC_miss   -> SetStats(0);

  h_MCLength_p_TPC_signal -> SetLineColor(1);
  h_MCLength_p_TPC_miss   -> SetLineColor(2);

  h_MCLength_p_TPC_signal -> GetXaxis() -> SetTitle( "Length[cm]");
  h_MCLength_p_TPC_signal -> GetYaxis() -> SetTitle( "Entries[#]");

  h_MCLength_p_TPC_miss   -> Draw();
  h_MCLength_p_TPC_signal -> Draw("same");

  c->SaveAs("MCLength_p_signalTPC.root") ;
  c->Clear();
  //**********************************************************************//  
  h_recoDaughters_pi -> SetStats(0) ; 
  h_recoDaughters_p  -> SetStats(0) ;

  h_recoDaughters_mu -> SetLineColor(1) ; 
  h_recoDaughters_pi -> SetLineColor(2) ; 
  h_recoDaughters_p  -> SetLineColor(4) ;
    
  h_recoDaughters_mu -> GetXaxis() -> SetTitle( "#Daughters primary" ) ; 
  h_recoDaughters_mu -> GetYaxis() -> SetTitle( "#events" ) ; 

  h_recoDaughters_mu -> Draw() ; 
  h_recoDaughters_pi -> Draw("same") ; 
  h_recoDaughters_p  -> Draw("same") ; 

  c->SaveAs("Daughters_primary.root") ; 
  c->Clear();
  //**********************************************************************//  
  h_reco3Daughters_pi -> SetStats(0) ; 
  h_reco3Daughters_p  -> SetStats(0) ;

  h_reco3Daughters_mu -> SetLineColor(1) ; 
  h_reco3Daughters_pi -> SetLineColor(2) ; 
  h_reco3Daughters_p  -> SetLineColor(4) ;
    
  h_reco3Daughters_mu -> GetXaxis() -> SetTitle( "#Daughters primary" ) ; 
  h_reco3Daughters_mu -> GetYaxis() -> SetTitle( "#events" ) ; 

  h_reco3Daughters_mu -> Draw() ; 
  h_reco3Daughters_pi -> Draw("same") ; 
  h_reco3Daughters_p  -> Draw("same") ; 

  c->SaveAs("Daughters_3.root") ; 
  c->Clear();
  //**********************************************************************//  
  
  
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
  mapTPiDaughterPdg.clear();

  // RECO INFO
  tr_id_best = 0;
  is_candidate = false ; 
  
  nu_reconstructed = false ;
  nu_rvertex_contained = true ; 
  numb_nu_reco = 0 ;
  
  for( int i = 0 ; i < 3 ; ++i ){
    nu_reco_vertex[i] = 0 ;
  }

  error_vertex_reco = 0 ; 

  map_MCID_RecoID.clear();
  map_RecoContained.clear();
  map_PandoraPDG.clear();
  map_RecoDaughters.clear();
  map_RecoHiearchy.clear();
  map_RecoKEnergy.clear();
  map_RecoLength.clear();
  map_RecoHits.clear();
}

DEFINE_ART_MODULE(test::NeutrinoTopologyAnalyzer)
