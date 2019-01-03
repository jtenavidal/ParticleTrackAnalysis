#include "TrackFitter.h"
#include <iostream>
#include <string>
#include "TH1.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraph.h"
#include "TVector3.h"
#include "math.h"

/**
* CONSTRUCTORS
* 1- Read from Track [X]
* 2- Read from root File [X] <- can work offline
* 3- Read from LArSoft [ ]
*/
TrackFitter::TrackFitter( const Track & p_track ){
  _hits = p_track.size();
  _vertex_position = p_track[0] ;// AccessVertex( p_track );
  _end_position = p_track[0] ;//AccessEnd( p_track );
  _particle_track = p_track;
  // Should add more info ...
}


TrackFitter::TrackFitter( const std::string & track_file_path ){
  TFile track_file( track_file_path.c_str() );
  TTree *recoTrack_tree   = (TTree*) track_file.Get("recoTrack_tree") ;
  // can add other trees later if needed

  TBranch *tr_x    = recoTrack_tree->GetBranch("tr_x");
  TBranch *tr_y    = recoTrack_tree->GetBranch("tr_y");
  TBranch *tr_z    = recoTrack_tree->GetBranch("tr_z");
  TBranch *tr_t    = recoTrack_tree->GetBranch("tr_t");
  TBranch *tr_dEdx = recoTrack_tree->GetBranch("tr_dEdx");
  TBranch *tr_dQdx = recoTrack_tree->GetBranch("tr_dQdx");
  // add here other info if needed

  _hits = recoTrack_tree->GetEntries();
  Hit_level track_hit ;
  for( unsigned int i = 0; i < _hits; ++i ){
    recoTrack_tree->GetEntry(i);
    track_hit.clear();
    track_hit.push_back( tr_x->GetLeaf("tr_x")->GetValue());
    track_hit.push_back( tr_y->GetLeaf("tr_y")->GetValue());
    track_hit.push_back( tr_z->GetLeaf("tr_z")->GetValue());
    track_hit.push_back( tr_t->GetLeaf("tr_t")->GetValue());
    _reco_dEdx.push_back( tr_dEdx->GetLeaf("tr_dEdx")->GetValue());
    _reco_dQdx.push_back( tr_dQdx->GetLeaf("tr_dQdx")->GetValue());
    _particle_track.push_back(track_hit);
  }
  _vertex_position = _particle_track[0] ;// AccessVertex( p_track );
  _end_position = _particle_track[_hits-1] ;//AccessEnd( p_track );
  // Should add more info ...
}

/**
* FUNCTIONS
* 1 - Get Properties
*/

unsigned int TrackFitter::GetHits( ){
  return _hits ;
}

Hit_level TrackFitter::AccessVertex( ) {
  return _vertex_position ;
}

Hit_level TrackFitter::AccessEnd( ) {
  return _end_position ;
}

Hit_level TrackFitter::GetdEdx( ){
  return _reco_dEdx ;
}

Hit_level TrackFitter::GetdQdx( ){
  return _reco_dQdx ;
}

Track TrackFitter::GetTrack( ){
  return _particle_track ;
}

/**
* FUNCTIONS
* 1 - Get Properties
* 2 - To check or save information
*/

void TrackFitter::SaveTrack( std::string const & path ) const {

  TH3D *h_track = new TH3D("h_track", " Particle Track ", int(_hits/10),
 _particle_track[0][0], _particle_track[_hits-1][0], int(_hits/10),
_particle_track[0][1], _particle_track[_hits-1][1], int(_hits/10),
_particle_track[0][2], _particle_track[_hits-1][2] );

  for( unsigned int i = 0; i < _hits; ++i ){
    h_track-> Fill(_particle_track[i][0], _particle_track[i][1], _particle_track[i][2]); // should add dEdX
  //  std::cout<<_particle_track[i][0]<< " "<<_particle_track[i][1]<< " " << _particle_track[i][2] <<std::endl;
  }

  TCanvas *c = new TCanvas();
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
  h_track->SetLineColor(2);
  h_track->GetXaxis()->SetTitle("X");
  h_track->GetYaxis()->SetTitle("Y");
  h_track->GetXaxis()->SetTitle("Z");
  h_track->Draw("hist");
  h_track->Draw("BOX2Z");
  c->SaveAs((path+".root").c_str());
  c->Clear();
}


Track TrackFitter::Straight( ) {
  Track track_hipotesis ;
  Hit_level v ;
  double distx, disty, distz ;

  track_hipotesis.reserve(_hits-1);
  distx = ( _end_position[0] - _vertex_position[0] )/ (_hits-1);
  disty = ( _end_position[1] - _vertex_position[1] )/ (_hits-1);
  distz = ( _end_position[2] - _vertex_position[2] )/ (_hits-1);

  for ( unsigned int i = 0 ; i < _hits ; ++i ){
    v.clear();
    v.push_back( _vertex_position[0]+i*distx );
    v.push_back( _vertex_position[1]+i*disty );
    v.push_back( _vertex_position[2]+i*distz );
    track_hipotesis.push_back(v);
  }
  return track_hipotesis ;
}


void TrackFitter::PrintHipotesis( const std::string & path ) {
  Track track_hipotesis = Straight ( ) ;

  TH3D *h_track = new TH3D("h_track", " Particle Track ", int(_hits/10),
   _particle_track[0][0], _particle_track[_hits-1][0], int(_hits/10),
  _particle_track[0][1], _particle_track[_hits-1][1], int(_hits/10),
  _particle_track[0][2], _particle_track[_hits-1][2] );

  for( unsigned int i = 0; i < _hits; ++i ){
    h_track-> Fill(_particle_track[i][0], _particle_track[i][1], _particle_track[i][2]); // should add dEdX
    //  std::cout<<_particle_track[i][0]<< " "<<_particle_track[i][1]<< " " << _particle_track[i][2] <<std::endl;
  }

  TCanvas *c = new TCanvas();
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
  h_track->SetLineColor(2);
  h_track->GetXaxis()->SetTitle("X");
  h_track->GetYaxis()->SetTitle("Y");
  h_track->GetXaxis()->SetTitle("Z");
  h_track->Draw("hist");
  h_track->Draw("BOX2Z");

  TH3D *h_hipotesis = new TH3D("h_hipotesis", "Reco hipotesis", int(_hits/10), _vertex_position[0], _end_position[0], int(_hits/10), _vertex_position[1], _end_position[1], int(_hits/10), _vertex_position[2], _end_position[2]);

  for ( unsigned int t_hit = 0 ; t_hit < track_hipotesis.size(); ++t_hit ){
      h_hipotesis -> Fill( track_hipotesis[t_hit][0], track_hipotesis[t_hit][1], track_hipotesis[t_hit][2]);
  }

  h_hipotesis->Draw("hist SAME ");
  c->SaveAs( (path+"_hipotesis.root").c_str() ) ;
  c->Clear();
}

double TrackFitter::FitToLine( ){
  /*
   * This function calculates the distance btw the hipotesis track and the reconstructed points, and returns the
   * cummulative distance.
   */

  TVector3 vertex( _vertex_position[0], _vertex_position[1], _vertex_position[2]);
  TVector3 direction, vertex_P ;
  std::vector< double > residual ;
  double residual_total = 0 ;
  direction.SetX(( _end_position[0] - _vertex_position[0])/ (_hits-1));
  direction.SetY(( _end_position[1] - _vertex_position[1])/ (_hits-1));
  direction.SetZ(( _end_position[2] - _vertex_position[2])/ (_hits-1));

  for( unsigned int i = 0; i < _hits; ++i ){
    vertex_P.SetX( _particle_track[i][0] - _vertex_position[0] );
    vertex_P.SetY( _particle_track[i][1] - _vertex_position[1] );
    vertex_P.SetZ( _particle_track[i][2] - _vertex_position[2] );
    residual.push_back((vertex_P.Cross( direction )).Mag()/direction.Mag()); // use this to find possible kinked trakcs: angle btw them?
    residual_total += residual[i] ;
  }

  return residual_total/_hits; // returns mean value
}


void TrackFitter::PrintdEdx( const std::string & path ) const {

  TH1F * h_dEdx = new TH1F( "h_dEdx", "dEdx", _reco_dEdx.size(), 0, _reco_dEdx.size() );
  for (unsigned int i = 0 ; i < _reco_dEdx.size() ; ++i ){
    h_dEdx -> Fill ( i , _reco_dEdx[i] ) ;
  }

  TCanvas *c = new TCanvas() ;
  h_dEdx -> Draw() ;
  c->SaveAs( (path+"_dEdx.root").c_str() ) ;

}


double TrackFitter::Linearity( ){
  return 0.;
}

std::vector< std::vector< double > > TrackFitter::MeanPosition( const int & window ){
  unsigned int starting_hit , end_hit;
  double muI_x, muI_y, muI_z;
  std::vector< double > mu_x, mu_y, mu_z ;
  std::vector< std::vector< double > > mu ;

  for( int i = 0; i < int(_hits); ++i ){
    muI_x = 0. ;
    muI_y = 0. ;
    muI_z = 0. ;
    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = (i+1) + window ;
    } else if ( i + window > int(_hits)) {
      starting_hit = i - window ;
      end_hit = _hits ;
    } else {
      starting_hit = i - window ;
      end_hit = (i+1) + window ;
    }

    for( unsigned int j = starting_hit ; j < end_hit -1 ; ++j ){
      muI_x += _particle_track[i][0];
      muI_y += _particle_track[i][1];
      muI_z += _particle_track[i][2];
      //std::cout<< j << std::endl;
    }
    mu_x.push_back( muI_x/( end_hit - starting_hit -1 ) );
    mu_y.push_back( muI_y/( end_hit - starting_hit -1 ) );
    mu_z.push_back( muI_z/( end_hit - starting_hit -1 ) );
    std::cout<< "hit number = "<< i << ":   "<<"muI_x= " << muI_x/( end_hit - starting_hit -1 )  <<" muI_y= " << muI_y/( end_hit - starting_hit -1 )  << " muI_z= "<< muI_z/( end_hit - starting_hit -1 )  << std::endl;
  }

  mu.push_back( mu_x );
  mu.push_back( mu_y );
  mu.push_back( mu_z );
  return mu;
}