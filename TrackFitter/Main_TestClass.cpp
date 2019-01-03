#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TrackFitter.h"
#include "TGraph.h"

void MainTest(){
  std::cout<<"WORKS ? " << std::endl;
  TrackFitter reco_track( "output_track_example.root" ) ;
  reco_track.SaveTrack( "test_reco_track" ) ;


} // MainTest()
