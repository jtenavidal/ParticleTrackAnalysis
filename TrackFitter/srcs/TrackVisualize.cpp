#include "../include/TrackFitter.h"
#include <iostream>
#include <string>
#include "TH1.h"
#include "TH3D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraph.h"
#include "TVector3.h"
#include "math.h"
#include "TArrow.h"


/**
* FUNCTIONS
* 1 - Get Properties
* 2 - To check or save information
*/

void TrackFitter::SaveTrack( std::string const & path , const unsigned int & event_id_track ) {
  std::vector< std::vector< double > > min_Linearity_position = FindMinimumLinearityPosition( event_id_track ) ;
  int bins = int(_event_hits[event_id_track-1]/10) ;
  double min_x, min_y, min_z, max_x, max_y, max_z ;
  std::string title ;

  min_x = FindMinCoordinate( event_id_track, 0 ) ;
  min_y = FindMinCoordinate( event_id_track, 1 ) ;
  min_z = FindMinCoordinate( event_id_track, 2 ) ;
  max_x = FindMaxCoordinate( event_id_track, 0 ) ;
  max_y = FindMaxCoordinate( event_id_track, 1 ) ;
  max_z = FindMaxCoordinate( event_id_track, 2 ) ;

  TH3D *h_track = new TH3D("h_track", " Particle Track ", bins, min_x, max_x, bins, min_y, max_y, bins, min_z, max_z );
  TH3D *h_track_MCvertex = (TH3D * ) h_track ->Clone();
  TH3D *h_track_kink = (TH3D * ) h_track ->Clone();
  h_track_kink->SetName("h_track_kink") ;

  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);

  for( int i = 0; i < _event_hits[event_id_track-1]; ++i ){
    h_track-> Fill(_event_tracks[event_id_track-1][0][i], _event_tracks[event_id_track-1][1][i], _event_tracks[event_id_track-1][2][i] ) ; }//, _event_tracks[event_id_track-1][3][i]);  }

  if( _TPDG_Code_Primary[event_id_track-1] == 13 ) title = "  #mu^{-} gun event : " ;
  if( _TPDG_Code_Primary[event_id_track-1] == 211 ) title = " #pi^{-} gun event : " ;
  if( _Tnu_mu[event_id_track-1] != 0 ) title += std::to_string(_Tnu_mu[event_id_track-1])+" #mu^{-} " ;
  if( _Tnu_pi[event_id_track-1] != 0 ) title += std::to_string(_Tnu_pi[event_id_track-1])+" #pi^{-} " ;
  if( _Tnu_p[event_id_track-1] != 0 ) title += std::to_string(_Tnu_p[event_id_track-1])+" p " ;
  if( _Tnu_e[event_id_track-1] != 0 ) title += std::to_string(_Tnu_e[event_id_track-1])+" #e^{-} " ;
  if( _Tnu_photon[event_id_track-1] != 0 ) title += std::to_string(_Tnu_photon[event_id_track-1])+" #gamma " ;
  if( _Tnu_others[event_id_track-1] != 0 ) title += std::to_string(_Tnu_others[event_id_track-1])+" others " ;
  std::cout<< _event_MC_vertex[event_id_track-1][0] << std::endl;
  h_track->SetTitle( title.c_str() ) ;
  TCanvas *c = new TCanvas();
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
  h_track->SetLineColor(2);
  h_track->GetXaxis()->SetTitle("X [cm]");
  h_track->GetYaxis()->SetTitle("Y [cm]");
  h_track->GetZaxis()->SetTitle("Z [cm]");
  leg->AddEntry( h_track, " Track 3D trajectory ");
  h_track->Draw("BOX2Z");

  h_track_MCvertex->Fill( _event_MC_vertex[event_id_track-1][0], _event_MC_vertex[event_id_track-1][1], _event_MC_vertex[event_id_track-1][2]) ;
  h_track_MCvertex->SetLineColor(2) ;
  h_track_MCvertex->SetLineWidth(3) ;
  h_track_MCvertex->SetFillColor(kRed);
  h_track_MCvertex->SetFillStyle(3001);
  h_track_MCvertex->Draw("BOX same") ;
  leg->AddEntry(h_track_MCvertex, "Truth vertex position" );

  if( min_Linearity_position.size() > 0 ) { // can remove
    for( unsigned int i = 0 ; i < min_Linearity_position.size() ; ++i ) {
      h_track_kink -> Fill( min_Linearity_position[i][0], min_Linearity_position[i][1], min_Linearity_position[i][2] ) ;
    }
    h_track_kink->SetLineColor(3) ;
    h_track_kink->SetLineWidth(3) ;
    h_track_kink->SetFillColor(kRed);
    h_track_kink->SetFillStyle(3004);
    h_track_kink->Draw("BOX same") ;
    leg->AddEntry(h_track_kink, "Identified kink/s position" );
  }

  leg->Draw();
  c->SaveAs((path+".root").c_str());
  //  c->Clear();
}

void TrackFitter::PrintdEdx( const std::string & path , const unsigned int & event_id_track ) const {

  TH1F * h_dEdx = new TH1F( "h_dEdx", "dEdx", int(_event_hits[event_id_track-1]), 0, _event_hits[event_id_track-1] );
  for ( int i = 0 ; i < _event_hits[event_id_track-1] ; ++i ){
    h_dEdx -> Fill ( i ,  _event_tracks[event_id_track-1][3][i] );//_reco_dEdx[i] ) ;//
  }

  TCanvas *c = new TCanvas() ;
  h_dEdx -> Draw("hist") ;
  c->SaveAs( (path+"_dEdx.root").c_str() ) ;

}


void TrackFitter::PlotLinearityData( const std::string & path , const unsigned int & event_id_track ) {//}, const std::vector< double > & Data_1, const std::vector< double > & Data_2  ) {
  std::vector< double > corrP = LinearityData( _event_tracks[event_id_track-1][0], _event_tracks[event_id_track-1][1] );//Data_1, Data_2) ;
  TH1F * h_Linearity = new TH1F( "h_Linearity", "Linearity", corrP.size(), 0,  corrP.size() );
  TLegend * legend = new TLegend(0.15,0.15,0.35,0.35) ;

  int window =  int( _event_hits[event_id_track-1] * 0.07 ) ;
  if( window < 5 ) window = 5 ;

  gStyle->SetOptStat(0);
  h_Linearity->SetLineColor(46);
  h_Linearity->SetLineWidth(2);
  h_Linearity->SetLineStyle(1);
  h_Linearity->GetYaxis()->SetRangeUser(0,1);
  legend->AddEntry( h_Linearity , "r", "l") ;
  gStyle->SetOptStat(0);

  for (unsigned int i = 0 ; i < corrP.size() ; ++i ){
    h_Linearity -> Fill ( i, corrP[i] ) ;
  }

  TCanvas *c = new TCanvas() ;
  h_Linearity -> Draw("HIST L") ;
  legend->Draw();
  c->SaveAs( (path+"_LinearityX.root").c_str() ) ;

}

void TrackFitter::PlotLinearityTrack( const std::string & path , const unsigned int & event_id_track ) {
  std::vector< double > corrP_XY = LinearityData( _event_tracks[event_id_track-1][0], _event_tracks[event_id_track-1][1] ) ;
  std::vector< double > corrP_XZ = LinearityData( _event_tracks[event_id_track-1][0], _event_tracks[event_id_track-1][2] ) ;
  std::vector< double > corrP_YZ = LinearityData( _event_tracks[event_id_track-1][2], _event_tracks[event_id_track-1][1] ) ;
  int window =  int( _event_hits[event_id_track-1] * 0.07 ) ;
  if( window < 5 ) window = 5 ;

  TH1F * h_Linearity = new TH1F( "h_Linearity", "Linearity", corrP_XY.size(), 0,  corrP_XY.size() );
  TH1F * h_LinearityXY = new TH1F( "h_LinearityXY", "Linearity", corrP_XY.size(), 0,  corrP_XY.size() );
  TH1F * h_LinearityXZ = new TH1F( "h_LinearityXZ", "Linearity", corrP_XZ.size(), 0,  corrP_XZ.size() );
  TH1F * h_LinearityYZ = new TH1F( "h_LinearityYZ", "Linearity", corrP_YZ.size(), 0,  corrP_YZ.size() );
  TLegend * legend = new TLegend(0.15,0.15,0.35,0.35) ;

  gStyle->SetOptStat(0);
  h_Linearity->SetLineColor(46);
  h_Linearity->SetLineWidth(2);
  h_Linearity->SetLineStyle(1);
  h_Linearity->GetYaxis()->SetRangeUser(0,1);
  legend->AddEntry( h_Linearity , "r", "l") ;
  gStyle->SetOptStat(0);
  h_LinearityXY->SetLineColor(2);
  h_LinearityXY->SetLineStyle(6);
//  h_LinearityXY->GetYaxis()->SetRangeUser(0,1);
  legend->AddEntry( h_LinearityXY , "r_{XY}", "l") ;
  h_LinearityXZ->SetLineColor(3);
  h_LinearityXZ->SetLineStyle(2);
  legend->AddEntry( h_LinearityXZ , "r_{XZ}", "l") ;
  h_LinearityYZ->SetLineColor(4);
  h_LinearityYZ->SetLineStyle(3);
  legend->AddEntry( h_LinearityYZ , "r_{YZ}", "l") ;

  for (unsigned int i = 0 ; i < corrP_XY.size() ; ++i ){
    h_LinearityXY -> Fill ( i, corrP_XY[i] ) ;
  }
  for (unsigned int i = 0 ; i < corrP_XY.size() ; ++i ){
    h_Linearity -> Fill ( i, corrP_XZ[i]*corrP_XY[i]*corrP_YZ[i] ) ;
  }
  for (unsigned int i = 0 ; i < corrP_XZ.size() ; ++i ){
    h_LinearityXZ -> Fill ( i, corrP_XZ[i] ) ;
  }
  for (unsigned int i = 0 ; i < corrP_YZ.size() ; ++i ){
    h_LinearityYZ -> Fill ( i, corrP_YZ[i] ) ;
  }

  TCanvas *c = new TCanvas() ;
  h_Linearity -> Draw("HIST L") ;
  h_LinearityXY -> Draw("HIST L SAME") ;
  h_LinearityXZ -> Draw("HIST L SAME") ;
  h_LinearityYZ -> Draw("HIST L SAME") ;
  h_Linearity->GetXaxis()->SetTitle("hits");
  h_Linearity->GetYaxis()->SetTitle("r");
  legend->Draw();
  c->SaveAs( (path+"_LinearityX.root").c_str() ) ;

}
