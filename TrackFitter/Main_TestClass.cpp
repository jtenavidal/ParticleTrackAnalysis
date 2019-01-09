#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TrackFitter.h"
#include "TGraph.h"

void MainTest(){
/*  TrackFitter reco_track_ev1_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_example_pi.root" ) ;
  reco_track_ev1_pi.SaveTrack( "test_reco_track_ev1_pi" ) ;
  reco_track_ev1_pi.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_example" ) ;

  TrackFitter reco_track_ev2_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev2_pi.root" ) ;
  reco_track_ev2_pi.SaveTrack( "test_reco_track_ev2_pi" ) ;
  reco_track_ev2_pi.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev2_pi" ) ;

  TrackFitter reco_track_ev4_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev4_pi.root" ) ;
  reco_track_ev4_pi.SaveTrack( "test_reco_track_ev4_pi" ) ;
  reco_track_ev4_pi.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev4_pi" ) ;

  TrackFitter reco_track_ev7_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev7_pi.root" ) ;
  reco_track_ev7_pi.SaveTrack( "test_reco_track_ev7_pi" ) ;
  reco_track_ev7_pi.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev7_pi" ) ;

  TrackFitter reco_track_ev19_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev19_pi.root" ) ;
  reco_track_ev19_pi.SaveTrack( "test_reco_track_ev19_pi" ) ;
  reco_track_ev19_pi.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev19_pi" ) ;

  TrackFitter reco_track_ev10_mu( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev10_mu.root" ) ;
  reco_track_ev10_mu.SaveTrack( "test_reco_track_ev10_mu" ) ;
  reco_track_ev10_mu.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev10_mu" ) ;

  TrackFitter reco_track_ev11_mu( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev11_mu.root" ) ;
  reco_track_ev11_mu.SaveTrack( "test_reco_track_ev11_mu" ) ;
  reco_track_ev11_mu.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev11_mu" ) ;

  TrackFitter reco_track_ev15_mu( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev15_mu.root" ) ;
  reco_track_ev15_mu.SaveTrack( "test_reco_track_ev15_mu" ) ;
  reco_track_ev15_mu.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev15_mu" ) ;

  TrackFitter reco_track_ev17_mu( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev17_mu.root" ) ;
  reco_track_ev17_mu.SaveTrack( "test_reco_track_ev75_mu" ) ;
  reco_track_ev17_mu.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev17_mu" ) ;

  TrackFitter reco_track_ev35_mu( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev35_mu.root" ) ;
  reco_track_ev35_mu.SaveTrack( "test_reco_track_ev35_mu" ) ;
  reco_track_ev35_mu.PlotLinearity( 15, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_ev35_mu" ) ;
*/

TrackFitter reco_track_ev1_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_pi.root" , 34) ;//_pi
reco_track_ev1_pi.SaveTrack( "Track_testing_new_constructor" ) ;
reco_track_ev1_pi.PlotLinearity( 5, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_testing_new_constructor" ) ;


} // MainTest()
