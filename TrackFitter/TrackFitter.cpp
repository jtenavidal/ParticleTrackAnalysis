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

/* THIS CONSTRUCTOR WAS A TEST FOR A TEMPORARY TTREE : not removing it yet
TrackFitter::TrackFitter( const std::string & track_file_path ){
  TFile track_file( track_file_path.c_str() );
  TTree *recoTrack_tree   = (TTree*) track_file.Get("recoTrack_tree") ;
  // can add other trees later if needed

  TBranch *tr_x    = recoTrack_tree->GetBranch("tr_x");
  TBranch *tr_y    = recoTrack_tree->GetBranch("tr_y");
  TBranch *tr_z    = recoTrack_tree->GetBranch("tr_z");
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
    _reco_dEdx.push_back( tr_dEdx->GetLeaf("tr_dEdx")->GetValue());
    _reco_dQdx.push_back( tr_dQdx->GetLeaf("tr_dQdx")->GetValue());
    _particle_track.push_back(track_hit);
  }
  _vertex_position = _particle_track[0] ;// AccessVertex( p_track );
  _end_position = _particle_track[_hits-1] ;//AccessEnd( p_track );
  // Should add more info ...
}
*/
 // THIS CONSTRUCTOR WORKS WITH THE OUTPUT TREE ( event, mc, reco )

 TrackFitter::TrackFitter( const std::string & track_file_path ){
   TFile track_file( track_file_path.c_str() );
   // Event Tree
   TTree * event_tree   = (TTree*) track_file.Get("event_tree") ;
   TBranch * event_id = event_tree->GetBranch("event_id");
   // MC Tree
   TTree * mcparticle_tree   = (TTree*) track_file.Get("mcparticle_tree") ;
   TBranch * mc_event_id = mcparticle_tree->GetBranch("event_id");
   TBranch * Track_ID = mcparticle_tree->GetBranch("fTrack_ID");
   TBranch * trueEnergy = mcparticle_tree->GetBranch("ftrueEnergy");
   TBranch * PDG_Code = mcparticle_tree->GetBranch("fPDG_Code");
   TBranch * Mass = mcparticle_tree->GetBranch("fMass");
   TBranch * Px = mcparticle_tree->GetBranch("fpx");
   TBranch * Py = mcparticle_tree->GetBranch("fpy");
   TBranch * Pz = mcparticle_tree->GetBranch("fpz");
   TBranch * Pt = mcparticle_tree->GetBranch("fpt");
   TBranch * P = mcparticle_tree->GetBranch("fp");
   TBranch * Num_Daughters = mcparticle_tree->GetBranch("fNum_Daughters");
   TBranch * Daughter_mu = mcparticle_tree->GetBranch("fDaughter_mu");
   TBranch * Daughter_pi = mcparticle_tree->GetBranch("fDaughter_pi");
   TBranch * Daughter_e = mcparticle_tree->GetBranch("fDaughter_e");
   TBranch * Daughter_p = mcparticle_tree->GetBranch("fDaughter_p");
   TBranch * Daughter_n = mcparticle_tree->GetBranch("fDaughter_n");
   TBranch * Daughter_photon = mcparticle_tree->GetBranch("fDaughter_photon");
   TBranch * Daughter_other = mcparticle_tree->GetBranch("fDaughter_other");
   TBranch * MC_Lenght= mcparticle_tree->GetBranch("fMC_Lenght");
   // RECO Tree
   TTree * recoparticle_tree   = (TTree*) track_file.Get("recoparticle_tree") ;
   TBranch * reco_event_id = recoparticle_tree->GetBranch("event_id");
   TBranch * pdg_primary = recoparticle_tree->GetBranch("r_pdg_primary");
   TBranch * nu_daughters = recoparticle_tree->GetBranch("r_nu_daughters");
   TBranch * mu_daughters = recoparticle_tree->GetBranch("r_mu_daughters");
   TBranch * pi_daughters = recoparticle_tree->GetBranch("r_pi_daughters");
   TBranch * e_daughters = recoparticle_tree->GetBranch("r_e_daughters");
   TBranch * p_daughters = recoparticle_tree->GetBranch("r_p_daughters");
   TBranch * other_daughters = recoparticle_tree->GetBranch("r_other_daughters");
   TBranch * Length = recoparticle_tree->GetBranch("rLength");
   TBranch * Momentum = recoparticle_tree->GetBranch("rMomentum");
   TBranch * nu_hits = recoparticle_tree->GetBranch("rnu_hits");
   TBranch * r_chi2_mu = recoparticle_tree->GetBranch("r_chi2_mu");
   TBranch * r_chi2_pi = recoparticle_tree->GetBranch("r_chi2_pi");
   TBranch * r_chi2_p = recoparticle_tree->GetBranch("r_chi2_p");
   TBranch * r_PIDA = recoparticle_tree->GetBranch("r_PIDA");
   TBranch * r_missing_energy = recoparticle_tree->GetBranch("r_missing_energy");
   TBranch * r_KineticEnergy = recoparticle_tree->GetBranch("r_KineticEnergy");
   TBranch * rdEdx_size = recoparticle_tree->GetBranch("rdEdx_size");
   TBranch * rdQdx_size = recoparticle_tree->GetBranch("rdQdx_size");
   TBranch * r_track_x = recoparticle_tree->GetBranch("r_track_x");
   TBranch * r_track_y = recoparticle_tree->GetBranch("r_track_y");
   TBranch * r_track_z = recoparticle_tree->GetBranch("r_track_z");
   TBranch * r_dEdx = recoparticle_tree->GetBranch("r_dEdx");
   TBranch * r_dQdx = recoparticle_tree->GetBranch("r_dQdx");

   Hit_level track_hit;

   for( unsigned int i = 0; i < recoparticle_tree->GetEntries(); ++i ){
     _particle_track.clear();
     _reco_dEdx.clear();
     _reco_dQdx.clear();

     recoparticle_tree->GetEntry(i);
     _hits = nu_hits->GetLeaf("rnu_hits")->GetValue() ;
     _event_hits.push_back( _hits ) ;

     for( int j = 0; j < _hits; ++j ) {
       track_hit.clear();
       track_hit.push_back( r_track_x->GetLeaf("r_track_x")->GetValue(j));
       track_hit.push_back( r_track_y->GetLeaf("r_track_y")->GetValue(j));
       track_hit.push_back( r_track_z->GetLeaf("r_track_z")->GetValue(j));
       _particle_track.push_back(track_hit);
     }
     _event_tracks.push_back( _particle_track ) ;

     _dEdx_size = recoparticle_tree->GetLeaf("rdEdx_size")->GetValue() ;
     _dQdx_size = recoparticle_tree->GetLeaf("rdQdx_size")->GetValue() ;
     _event_dEdx_size.push_back( _dEdx_size ) ;
     _event_dQdx_size.push_back( _dQdx_size );

     for( int j = 0; j < _dEdx_size; ++j ) _reco_dEdx.push_back( recoparticle_tree->GetLeaf("r_dEdx")->GetValue(j));
     for( int j = 0; j < _dQdx_size; ++j ) _reco_dQdx.push_back( recoparticle_tree->GetLeaf("r_dQdx")->GetValue(j));
     _event_reco_dEdx.push_back( _reco_dEdx ) ;
     _event_reco_dQdx.push_back( _reco_dQdx ) ;

     _event_vertex.push_back( _event_tracks[i][0] ) ;// AccessVertex( p_track );
     _event_end.push_back( _event_tracks[i][_event_hits[i]-1] ) ;//AccessEnd( p_track );
     // Should add more info ...

   }

 }


 TrackFitter::TrackFitter( const std::string & track_file_path, const unsigned int & event_id_track ){
   TFile track_file( track_file_path.c_str() );
   // Event Tree
   TTree * event_tree   = (TTree*) track_file.Get("event_tree") ;
   TBranch * event_id = event_tree->GetBranch("event_id");
   // MC Tree
   TTree * mcparticle_tree   = (TTree*) track_file.Get("mcparticle_tree") ;
   TBranch * mc_event_id = mcparticle_tree->GetBranch("event_id");
   TBranch * Track_ID = mcparticle_tree->GetBranch("fTrack_ID");
   TBranch * trueEnergy = mcparticle_tree->GetBranch("ftrueEnergy");
   TBranch * PDG_Code = mcparticle_tree->GetBranch("fPDG_Code");
   TBranch * Mass = mcparticle_tree->GetBranch("fMass");
   TBranch * Px = mcparticle_tree->GetBranch("fpx");
   TBranch * Py = mcparticle_tree->GetBranch("fpy");
   TBranch * Pz = mcparticle_tree->GetBranch("fpz");
   TBranch * Pt = mcparticle_tree->GetBranch("fpt");
   TBranch * P = mcparticle_tree->GetBranch("fp");
   TBranch * Num_Daughters = mcparticle_tree->GetBranch("fNum_Daughters");
   TBranch * Daughter_mu = mcparticle_tree->GetBranch("fDaughter_mu");
   TBranch * Daughter_pi = mcparticle_tree->GetBranch("fDaughter_pi");
   TBranch * Daughter_e = mcparticle_tree->GetBranch("fDaughter_e");
   TBranch * Daughter_p = mcparticle_tree->GetBranch("fDaughter_p");
   TBranch * Daughter_n = mcparticle_tree->GetBranch("fDaughter_n");
   TBranch * Daughter_photon = mcparticle_tree->GetBranch("fDaughter_photon");
   TBranch * Daughter_other = mcparticle_tree->GetBranch("fDaughter_other");
   TBranch * MC_Lenght= mcparticle_tree->GetBranch("fMC_Lenght");
   // RECO Tree
   TTree * recoparticle_tree   = (TTree*) track_file.Get("recoparticle_tree") ;
   TBranch * reco_event_id = recoparticle_tree->GetBranch("event_id");
   TBranch * pdg_primary = recoparticle_tree->GetBranch("r_pdg_primary");
   TBranch * nu_daughters = recoparticle_tree->GetBranch("r_nu_daughters");
   TBranch * mu_daughters = recoparticle_tree->GetBranch("r_mu_daughters");
   TBranch * pi_daughters = recoparticle_tree->GetBranch("r_pi_daughters");
   TBranch * e_daughters = recoparticle_tree->GetBranch("r_e_daughters");
   TBranch * p_daughters = recoparticle_tree->GetBranch("r_p_daughters");
   TBranch * other_daughters = recoparticle_tree->GetBranch("r_other_daughters");
   TBranch * Length = recoparticle_tree->GetBranch("rLength");
   TBranch * Momentum = recoparticle_tree->GetBranch("rMomentum");
   TBranch * rnu_hits = recoparticle_tree->GetBranch("rnu_hits");
   TBranch * r_chi2_mu = recoparticle_tree->GetBranch("r_chi2_mu");
   TBranch * r_chi2_pi = recoparticle_tree->GetBranch("r_chi2_pi");
   TBranch * r_chi2_p = recoparticle_tree->GetBranch("r_chi2_p");
   TBranch * r_PIDA = recoparticle_tree->GetBranch("r_PIDA");
   TBranch * r_missing_energy = recoparticle_tree->GetBranch("r_missing_energy");
   TBranch * r_KineticEnergy = recoparticle_tree->GetBranch("r_KineticEnergy");
   TBranch * rdEdx_size = recoparticle_tree->GetBranch("rdEdx_size");
   TBranch * rdQdx_size = recoparticle_tree->GetBranch("rdQdx_size");
   TBranch * r_track_x = recoparticle_tree->GetBranch("r_track_x");
   TBranch * r_track_y = recoparticle_tree->GetBranch("r_track_y");
   TBranch * r_track_z = recoparticle_tree->GetBranch("r_track_z");
   TBranch * r_dEdx = recoparticle_tree->GetBranch("r_dEdx");
   TBranch * r_dQdx = recoparticle_tree->GetBranch("r_dQdx");

   Hit_level track_hit;

   recoparticle_tree->GetEntry( event_id_track -1 );
   _hits = rnu_hits->GetLeaf("rnu_hits")->GetValue();
   _particle_track.clear();

   for( int j = 0; j < _hits; ++j ) {
     track_hit.clear();
     track_hit.push_back( r_track_x->GetLeaf("r_track_x")->GetValue(j));
     track_hit.push_back( r_track_y->GetLeaf("r_track_y")->GetValue(j));
     track_hit.push_back( r_track_z->GetLeaf("r_track_z")->GetValue(j));
     _particle_track.push_back(track_hit);
   }
    _dEdx_size = recoparticle_tree->GetLeaf("rdEdx_size")->GetValue( );
    _dQdx_size = recoparticle_tree->GetLeaf("rdQdx_size")->GetValue( );
    for( int j = 0; j < _dEdx_size; ++j ) _reco_dEdx.push_back( recoparticle_tree->GetLeaf("r_dEdx")->GetValue(j));
    for( int j = 0; j < _dQdx_size; ++j ) _reco_dQdx.push_back( recoparticle_tree->GetLeaf("r_dQdx")->GetValue(j));
   _vertex_position = _particle_track[0] ;// AccessVertex( p_track );
   _end_position = _particle_track[_hits-1] ;//AccessEnd( p_track );
   // Should add more info ...

 }


//*/
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
/*
std::vector< std::vector< float > > TrackFitter::GetdEdx( ){
  return _reco_dEdx ;
}

std::vector< std::vector< float > > TrackFitter::GetdQdx( ){
  return _reco_dQdx ;
}*/

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

  for( int i = 0; i < _hits; ++i ){
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

  for ( int i = 0 ; i < _hits ; ++i ){
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

  for( int i = 0; i < _hits; ++i ){
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

  for( int i = 0; i < _hits; ++i ){
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

std::vector< std::vector< double > > TrackFitter::Linearity( const int & window ){
  // It calculates the linearity from Pearson correlation coefficient (corrP)
  unsigned int starting_hit , end_hit;
  std::vector< double > corrP_xy, corrP_xz, corrP_yz;
  std::vector<  std::vector< double > > corrP ;//[xy, xz, yz]
  std::vector< std::vector< double > > cov, dev ;
  cov = CovPosition( window ) ;
  dev = DevPosition( window ) ;

  for( int i = 0; i < int(_hits); ++i ){
    corrP_xy.push_back( sqrt(std::pow(cov[0][i]/(dev[0][i]*dev[1][i]),2)) ) ;
    corrP_xz.push_back( sqrt(std::pow(cov[1][i]/(dev[0][i]*dev[2][i]),2)) ) ;
    corrP_yz.push_back( sqrt(std::pow(cov[2][i]/(dev[2][i]*dev[1][i]),2)) ) ;
  }

  corrP.push_back( corrP_xy );
  corrP.push_back( corrP_xz );
  corrP.push_back( corrP_yz );

  return corrP; // corrP[x,y,z]
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
      end_hit = i + window ;
    } else if ( i + window >= int(_hits)) {
      starting_hit = i - window ;
      end_hit = _hits ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }

    for( unsigned int j = starting_hit ; j < end_hit ; ++j ){
      muI_x += _particle_track[j][0];
      muI_y += _particle_track[j][1];
      muI_z += _particle_track[j][2];
    }
    mu_x.push_back( muI_x/( end_hit - starting_hit ) );
    mu_y.push_back( muI_y/( end_hit - starting_hit ) );
    mu_z.push_back( muI_z/( end_hit - starting_hit ) );
  }
  mu.push_back( mu_x );
  mu.push_back( mu_y );
  mu.push_back( mu_z );
  return mu;
}


std::vector< std::vector< double > > TrackFitter::DevPosition( const int & window ){
  unsigned int starting_hit , end_hit;
  double devI_x, devI_y, devI_z;
  std::vector< double > dev_x, dev_y, dev_z ;
  std::vector< std::vector< double > > dev, mean ;
  mean = MeanPosition( window ) ;

  for( int i = 0; i < int(_hits); ++i ){
    devI_x = 0. ;
    devI_y = 0. ;
    devI_z = 0. ;
    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = i + window ;
    } else if ( i + window >= int(_hits)) {
      starting_hit = i - window ;
      end_hit = _hits ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }
    for( unsigned int j = starting_hit ; j < end_hit ; ++j ){
      devI_x += std::pow( _particle_track[j][0]-mean[0][i] ,2);
      devI_y += std::pow( _particle_track[j][1]-mean[1][i] ,2);
      devI_z += std::pow( _particle_track[j][2]-mean[2][i] ,2);
    }
    dev_x.push_back( std::sqrt( devI_x/( end_hit - starting_hit -1 ) ) );
    dev_y.push_back( std::sqrt( devI_y/( end_hit - starting_hit -1 ) ) );
    dev_z.push_back( std::sqrt( devI_z/( end_hit - starting_hit -1 ) ) );
  }

  dev.push_back( dev_x );
  dev.push_back( dev_y );
  dev.push_back( dev_z );
  return dev; // dev[x,y,z]
}


std::vector< std::vector< double > > TrackFitter::CovPosition( const int & window ){
  unsigned int starting_hit , end_hit;
  double covI_xy, covI_xz, covI_yz;
  std::vector< double > cov_xy, cov_xz, cov_yz ;
  std::vector< std::vector< double > > cov, mean ;
  mean = MeanPosition( window ) ;
  for( int i = 0; i < int(_hits); ++i ){

    covI_xy = 0. ;
    covI_xz = 0. ;
    covI_yz = 0. ;

    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = i + window ;
    } else if ( i + window >= int(_hits)) {
      starting_hit = i - window ;
      end_hit = _hits ;
    } else if ( i - window < 0  && i + window >= int(_hits)) {
      starting_hit = 0 ;
      end_hit = _hits ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }

    for( unsigned int j = starting_hit ; j < end_hit ; ++j ){
      covI_xy += (_particle_track[j][0]-mean[0][i])*(_particle_track[j][1]-mean[1][i]) ;
      covI_xz += (_particle_track[j][0]-mean[0][i])*(_particle_track[j][2]-mean[2][i]) ;
      covI_yz += (_particle_track[j][2]-mean[2][i])*(_particle_track[j][1]-mean[1][i]) ;
    }

    cov_xy.push_back( covI_xy/( end_hit - starting_hit -1 ) );
    cov_xz.push_back( covI_xz/( end_hit - starting_hit -1 ) );
    cov_yz.push_back( covI_yz/( end_hit - starting_hit -1 ) );
  }

  cov.push_back( cov_xy );
  cov.push_back( cov_xz );
  cov.push_back( cov_yz );
  return cov; // cov[xy,xz,yz]
}


void TrackFitter::PlotLinearity( const int & window, const std::string & path ) {
  std::vector< std::vector< double > > corrP = Linearity( window ) ;
  TH1F * h_Linearity = new TH1F( "h_Linearity", "Linearity", corrP[0].size(), 0,  corrP[0].size() );
  TH1F * h_LinearityXY = new TH1F( "h_LinearityXY", "Linearity", corrP[0].size(), 0,  corrP[0].size() );
  TH1F * h_LinearityXZ = new TH1F( "h_LinearityXZ", "Linearity", corrP[1].size(), 0,  corrP[1].size() );
  TH1F * h_LinearityYZ = new TH1F( "h_LinearityYZ", "Linearity", corrP[2].size(), 0,  corrP[2].size() );
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

  for (unsigned int i = 0 ; i < corrP[0].size() ; ++i ){
    h_LinearityXY -> Fill ( i, corrP[0][i] ) ;
  }
  for (unsigned int i = 0 ; i < corrP[0].size() ; ++i ){
    h_Linearity -> Fill ( i, corrP[0][i]*corrP[1][i]*corrP[2][i] ) ;
  }
  for (unsigned int i = 0 ; i < corrP[1].size() ; ++i ){
    h_LinearityXZ -> Fill ( i, corrP[1][i] ) ;
  }
  for (unsigned int i = 0 ; i < corrP[2].size() ; ++i ){
    h_LinearityYZ -> Fill ( i, corrP[2][i] ) ;
  }

  TCanvas *c = new TCanvas() ;
  h_Linearity -> Draw("HIST L") ;
  h_LinearityXY -> Draw("HIST L SAME") ;
  h_LinearityXZ -> Draw("HIST L SAME") ;
  h_LinearityYZ -> Draw("HIST L SAME") ;
  legend->Draw();
  c->SaveAs( (path+"_LinearityX.root").c_str() ) ;

}
