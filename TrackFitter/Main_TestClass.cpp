#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TrackFitter.h"
#include "TGraph.h"

void MainTest(){
  TrackFitter reco_track( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_example_mu.root" ) ;
  reco_track.SaveTrack( "test_reco_track" ) ;
  reco_track.PlotLinearity( 10, "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/Histograms/output_track_example" ) ;
} // MainTest()
