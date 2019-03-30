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

private:
  // Declare member data here ;
  // Labels
  std::string TruthLabel, ParticleLabel, HitFinderLabel, RecoTrackLabel, RecoShowerLabel, RecoPIDLabel, RecoCaloLabel ; 

  // Tree members
  TTree * mcevent_tree, * recoevent_tree ; 
  
  // Event Information
  int event_id ;

  // MC Event Information 
  int MCPdgCode, MCNumDaughters ;
  
};


test::NeutrinoTopologyAnalyzer::NeutrinoTopologyAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  TruthLabel = p.get<std::string>("TruthLabel","largeant");
  ParticleLabel = p.get<std::string>("PFParticleModule","pandora");
  HitFinderLabel = p.get<std::string>("HitFinderModule","linecluster");
  RecoTrackLabel = p.get<std::string>("RecoTrackLabel","pandoraTrack");
  RecoShowerLabel = p.get<std::string>("RecoShowerLabel","emshower");
  RecoCaloLabel = p.get<std::string>("RecoCaloLabel","pandoraCalo");
  RecoPIDLabel = p.get<std::string>("RecoPIDLabel","pandoraPid");
  
  this->reconfiure(p);
}

void test::NeutrinoTopologyAnalyzer::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  event_id = e.id().event();
  std::cout<< " Event ID = " << event_id <<std::endl;

}


void TrackID::MyAnalysis::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of required member function here.
  // Here you add an external fcl file to change configuration
}


void TrackID::MyAnalysis::beginJob( )
{
  // Declare trees and branches
  mcevent_tree   = new TTree( "mcevent_tree",    "MC event tree: True event track information" ) ;
  recoevent_tree = new TTree( "recoevent_tree",  "Reco event tree: reconstructed information of event") ;

  //set branches

  //set directory
  mvevent_tree -> SetDirectory(0);
  recoevent_tree -> SetDirectory(0);
}


void TrackID::MyAnalysis::endJob( )
{
  // Implementation of required member function here.
  // Write and close files. Other comments also fit here 
  TFile file("output_eventtree.root" , "RECREATE" );
  mcevent_tree ->Write();
  recoevent_tree ->Write();
  
  file.Write();
  file.Close();
  
  delete mcevent_tree ;
  delete recoevent_tree ; 
}

DEFINE_ART_MODULE(test::NeutrinoTopologyAnalyzer)
