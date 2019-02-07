#include <iostream>
#include "../include/TrackFitter.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

void MainTest(){
unsigned int event =  12; // Breaks for pions: 10, 15
TrackFitter reco_track_all( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_new_mu.root") ;
reco_track_all.SaveStatisticsTrueEvent("stats");
// TrackFitter reco_track_all( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_pi.root") ;

// std::cout<< " i = " << event << " reco = " << reco_track_all.GetHits( event ) <<std::endl ;
/*for( unsigned int i = 1 ; i < 50 ; ++i ) {
  std::cout<< " i = " << i << " reco = " << reco_track_all.GetHits( i ) <<std::endl ;
}*/

// reco_track_all.SaveTrack( ("Track_testing_new_constructor_") , event ) ;
// reco_track_all.PlotLinearityTrack( " testing ", event ); // reco_track_all.GetHits(event)*6/10
//for( int i = 1 ; i < 8 ; ++i ) {
//  reco_track_all.SaveTrack( ("Track_testing_new_constructor_"+std::to_string(i)) , i ) ;
//}
//reco_track_all.PrintdQdx( "test_dqdx", event);
// reco_track_all.TruthParticles( event ) ;
// reco_track_all.StatisticsKinks( ) ;




//MUONS-> "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_mu.root"

//PIONS-> "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_pi.root"


} // MainTest()1
