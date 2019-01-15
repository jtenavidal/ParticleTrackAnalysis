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
   TBranch * nu_hits = recoparticle_tree->GetBranch("rnu_hits");
   TBranch * r_chi2_mu = recoparticle_tree->GetBranch("r_chi2_mu");
   TBranch * r_chi2_pi = recoparticle_tree->GetBranch("r_chi2_pi");
   TBranch * r_chi2_p = recoparticle_tree->GetBranch("r_chi2_p");
   TBranch * r_PIDA = recoparticle_tree->GetBranch("r_PIDA");
   TBranch * r_missing_energy = recoparticle_tree->GetBranch("r_missing_energy");
   TBranch * r_KineticEnergy = recoparticle_tree->GetBranch("r_KineticEnergy");
   TBranch * rdQdx_size = recoparticle_tree->GetBranch("rdQdx_size");
   TBranch * r_track_x = recoparticle_tree->GetBranch("r_track_x");
   TBranch * r_track_y = recoparticle_tree->GetBranch("r_track_y");
   TBranch * r_track_z = recoparticle_tree->GetBranch("r_track_z");
   TBranch * r_dQdx = recoparticle_tree->GetBranch("r_dQdx");

   Hit_level track_hit_x, track_hit_y, track_hit_z;

   for( unsigned int i = 0; i < recoparticle_tree->GetEntries(); ++i ){
     _particle_track.clear();
     track_hit_x.clear();
     track_hit_y.clear();
     track_hit_z.clear();
     _reco_dQdx.clear();

     recoparticle_tree->GetEntry(i);
     _hits = nu_hits->GetLeaf("rnu_hits")->GetValue() ;
     _event_hits.push_back( _hits ) ;

     for( int j = 0; j < _hits; ++j ) {
       track_hit_x.push_back( r_track_x->GetLeaf("r_track_x")->GetValue(j));
       track_hit_y.push_back( r_track_y->GetLeaf("r_track_y")->GetValue(j));
       track_hit_z.push_back( r_track_z->GetLeaf("r_track_z")->GetValue(j));
     }
     _particle_track.push_back(track_hit_x);
     _particle_track.push_back(track_hit_y);
     _particle_track.push_back(track_hit_z);

     _event_tracks.push_back( _particle_track ) ;

     _dQdx_size = recoparticle_tree->GetLeaf("rdQdx_size")->GetValue() ;
     _event_dQdx_size.push_back( _dQdx_size );

     for( int j = 0; j < _dQdx_size; ++j ) _reco_dQdx.push_back( recoparticle_tree->GetLeaf("r_dQdx")->GetValue(j));
     _event_reco_dQdx.push_back( _reco_dQdx ) ;

     // HAVE TO FIX THIS:
    // _event_vertex.push_back( _event_tracks[i][0] ) ;// AccessVertex( p_track );
    // _event_end.push_back( _event_tracks[i][_event_hits[i]-1] ) ;//AccessEnd( p_track );
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
   TBranch * rnu_hits = recoparticle_tree->GetBranch("rnu_hits");
   TBranch * rnu_hits_size = recoparticle_tree->GetBranch("rnu_hits_size");
   TBranch * r_chi2_mu = recoparticle_tree->GetBranch("r_chi2_mu");
   TBranch * r_chi2_pi = recoparticle_tree->GetBranch("r_chi2_pi");
   TBranch * r_chi2_p = recoparticle_tree->GetBranch("r_chi2_p");
   TBranch * r_PIDA = recoparticle_tree->GetBranch("r_PIDA");
   TBranch * r_missing_energy = recoparticle_tree->GetBranch("r_missing_energy");
   TBranch * r_KineticEnergy = recoparticle_tree->GetBranch("r_KineticEnergy");
   TBranch * rdQdx_size = recoparticle_tree->GetBranch("rdQdx_size");
   TBranch * r_track_x = recoparticle_tree->GetBranch("r_track_x");
   TBranch * r_track_y = recoparticle_tree->GetBranch("r_track_y");
   TBranch * r_track_z = recoparticle_tree->GetBranch("r_track_z");
   TBranch * r_dQdx = recoparticle_tree->GetBranch("r_dQdx");

   Hit_level track_hit_x, track_hit_y, track_hit_z;
   int hits_size ; // number of hits
   recoparticle_tree->GetEntry( event_id_track -1 );
   _hits = rnu_hits->GetLeaf("rnu_hits")->GetValue(); // number of valid hits!
   _hits_size = recoparticle_tree->GetLeaf("rnu_hits_size")-> GetValue();
   _particle_track.clear();
   track_hit_x.clear();
   track_hit_y.clear();
   track_hit_z.clear();

   for( int j = 0; j < _hits; ++j ) {
     track_hit_x.push_back( r_track_x->GetLeaf("r_track_x")->GetValue(j));
     track_hit_y.push_back( r_track_y->GetLeaf("r_track_y")->GetValue(j));
     track_hit_z.push_back( r_track_z->GetLeaf("r_track_z")->GetValue(j));
   }

   _particle_track.push_back(track_hit_x);
   _particle_track.push_back(track_hit_y);
   _particle_track.push_back(track_hit_z);

    _dQdx_size = recoparticle_tree->GetLeaf("rdQdx_size")->GetValue( );
    for( int j = 0; j < _dQdx_size; ++j ) _reco_dQdx.push_back( recoparticle_tree->GetLeaf("r_dQdx")->GetValue(j));
    _vertex_position.push_back( _particle_track[0][0] ) ;// AccessVertex( p_track );
    _vertex_position.push_back( _particle_track[1][0] ) ;
    _vertex_position.push_back( _particle_track[2][0] ) ;
    _end_position.push_back( _particle_track[0][_hits-1] );
    _end_position.push_back( _particle_track[1][_hits-1] );
    _end_position.push_back( _particle_track[2][_hits-1] );

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

std::vector< float > TrackFitter::GetdQdx( ){
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
 _particle_track[0][0], _particle_track[0][_hits-1], int(_hits/10),
_particle_track[1][0], _particle_track[1][_hits-1], int(_hits/10), // need to define number of bins as a function of _hits to avoid bad memory allocation
_particle_track[2][0], _particle_track[2][_hits-1] );

  for( int i = 0; i < _hits; ++i ){
    h_track-> Fill(_particle_track[0][i], _particle_track[1][i], _particle_track[2][i]);  }

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
   _particle_track[0][0], _particle_track[0][_hits-1], int(_hits/10),
  _particle_track[1][0], _particle_track[1][_hits-1], int(_hits/10),
  _particle_track[2][0], _particle_track[2][_hits-1] );

  for( int i = 0; i < _hits; ++i ){
    h_track-> Fill(_particle_track[0][i], _particle_track[1][i], _particle_track[2][i]);
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
    vertex_P.SetX( _particle_track[0][i] - _vertex_position[0] );
    vertex_P.SetY( _particle_track[1][i] - _vertex_position[1] );
    vertex_P.SetZ( _particle_track[2][i] - _vertex_position[2] );
    residual.push_back((vertex_P.Cross( direction )).Mag()/direction.Mag()); // use this to find possible kinked trakcs: angle btw them?
    residual_total += residual[i] ;
  }

  return residual_total/_hits; // returns mean value
}

void TrackFitter::PrintdQdx( const std::string & path ) const {

  TH1F * h_dQdx = new TH1F( "h_dQdx", "dQdx", _dQdx_size, 0, _dQdx_size );
  for ( int i = 0 ; i < _dQdx_size ; ++i ){
    h_dQdx -> Fill ( i , _reco_dQdx[i] ) ;
  }

  TCanvas *c = new TCanvas() ;
  h_dQdx -> Draw() ;
  c->SaveAs( (path+"_dQdx.root").c_str() ) ;

}


// STATISTICS GENERAL FUNCTIONS
std::vector< double > TrackFitter::MeanData( const int & window, const std::vector<double> & data ){
  unsigned int starting_hit , end_hit;
  double muI_data;
  std::vector< double > mu_data ;

  for( int i = 0; i < int(data.size()); ++i ){
    muI_data = 0. ;
    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = i + window ;
    } else if ( i + window >= int(data.size())) {
      starting_hit = i - window ;
      end_hit = data.size() ;
    } else if ( i - window < 0  && i + window >= int(data.size())) {
      starting_hit = 0 ;
      end_hit = data.size() ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }

    for( unsigned int j = starting_hit ; j < end_hit ; ++j ){
      muI_data += data[j];
    }
    mu_data.push_back( muI_data/( end_hit - starting_hit ) );
  }
  return mu_data;
}

std::vector< double > TrackFitter::DevData( const int & window, const std::vector<double> & data ){
  unsigned int starting_hit , end_hit;
  double devI_data;
  std::vector< double > dev_data , mean ;
  mean = MeanData( window, data ) ;

  for( int i = 0; i < int(data.size()); ++i ){
    devI_data = 0. ;
    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = i + window ;
    } else if ( i + window >= int(data.size())) {
      starting_hit = i - window ;
      end_hit = data.size() ;
    } else if ( i - window < 0  && i + window >= int(data.size())) {
      starting_hit = 0 ;
      end_hit = data.size() ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }

    for( unsigned int j = starting_hit ; j < end_hit ; ++j ){
      devI_data += std::pow( data[j]-mean[i] ,2);
    }
    dev_data.push_back( std::sqrt( devI_data/( end_hit - starting_hit -1 ) ) );
  }

  return dev_data;
}

std::vector< double > TrackFitter::CovData( const int & window, const std::vector< double > & Data_1, const std::vector< double > & Data_2 ){
  unsigned int starting_hit , end_hit;
  double covI_12 ; // 1 - variable 1, 2 - second variable
  std::vector< double > cov_12 , mean1, mean2;
  mean1 = MeanData( window , Data_1 ) ; //variable 1
  mean2 = MeanData( window , Data_2 ) ; //variable 2

  if( Data_1.size() == Data_2.size() ) {
    for( int i = 0; i < int(Data_1.size()); ++i ) {

      covI_12 = 0. ;
      if ( i - window < 0 ) {
        starting_hit = 0 ;
        end_hit = i + window ;
      } else if ( i + window >= int(Data_1.size())) {
        starting_hit = i - window ;
        end_hit = Data_1.size() ;
      } else if ( i - window < 0  && i + window >= int(Data_1.size())) {
        starting_hit = 0 ;
        end_hit = Data_1.size() ;
      } else {
        starting_hit = i - window ;
        end_hit = i + window ;
      }

      for( unsigned int j = starting_hit ; j < end_hit ; ++j ){
        covI_12 += ( Data_1[j]-mean1[i])*( Data_2[j]-mean2[i]) ;
      }

      cov_12.push_back( covI_12/( end_hit - starting_hit -1 ) );
    }
  }

  return cov_12;
}

std::vector< double > TrackFitter::LinearityData( const int & window, const std::vector< double > & Data_1, const std::vector< double > & Data_2 ){
  // It calculates the linearity from Pearson correlation coefficient (corrP)
  unsigned int starting_hit , end_hit;
  std::vector< double > corrP_12, cov_12, dev1, dev2 ;
  cov_12 = CovData( window, Data_1, Data_2 ) ;
  dev1 = DevData( window, Data_1 ) ;
  dev2 = DevData( window, Data_2 ) ;

  if ( Data_1.size() == Data_2.size()) {
    for( int i = 0; i < int( Data_1.size() ); ++i ){
      corrP_12.push_back( sqrt(std::pow(cov_12[i]/(dev1[i]*dev2[i]),2)) ) ;
    }
  }
  return corrP_12;
}

void TrackFitter::PlotLinearityData( const int & window, const std::string & path ) {//}, const std::vector< double > & Data_1, const std::vector< double > & Data_2  ) {
  std::vector< double > corrP = LinearityData( window, _particle_track[0], _particle_track[1] );//Data_1, Data_2) ;
  TH1F * h_Linearity = new TH1F( "h_Linearity", "Linearity", corrP.size(), 0,  corrP.size() );
  TLegend * legend = new TLegend(0.15,0.15,0.35,0.35) ;

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

void TrackFitter::PlotLinearityTrack( const int & window, const std::string & path ) {
  std::vector< double > corrP_XY = LinearityData( window, _particle_track[0], _particle_track[1] ) ;
  std::vector< double > corrP_XZ = LinearityData( window, _particle_track[0], _particle_track[2] ) ;
  std::vector< double > corrP_YZ = LinearityData( window, _particle_track[2], _particle_track[1] ) ;

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
  legend->Draw();
  c->SaveAs( (path+"_LinearityX.root").c_str() ) ;

}


std::vector< TVector3 > TrackFitter::MeanDirectionData( const int & window ){
  // this will only be applied to the track information (x,y,z)
  unsigned int starting_hit , end_hit;
  TVector3 directionI ;
  std::vector< TVector3 > mean_direction ;
  for( int i = 0; i < int(_hits); ++i ){
    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = i + window ;
    } else if ( i + window >= int(_hits)) {
      starting_hit = i - window ;
      end_hit = _hits - 1 ;
    } else if ( i - window < 0  && i + window >= int(_hits)) {
      starting_hit = 0 ;
      end_hit = _hits - 1 ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }
    directionI.SetX( (_particle_track[0][end_hit]-_particle_track[0][starting_hit])/(end_hit-starting_hit) );
    directionI.SetY( (_particle_track[1][end_hit]-_particle_track[1][starting_hit])/(end_hit-starting_hit) );
    directionI.SetZ( (_particle_track[2][end_hit]-_particle_track[2][starting_hit])/(end_hit-starting_hit) );
    mean_direction.push_back( directionI );
  }
  return mean_direction;
}


std::vector< double > TrackFitter::AngleTrackDistribution( const int & window ) {
  // this will only be applied to the track information (x,y,z)
  unsigned int starting_hit , end_hit;
  std::vector< double > angle_distribution ;
  std::vector< TVector3 > mean_direction = MeanDirectionData( window ) ;
  TVector3 test1(1, 1, 0);
  TVector3 test2(1, -1, 0);
  for( int i = 0; i < int(_hits) - 1 ; ++i ){

  if( mean_direction[i].Angle(mean_direction[i+1]) > 1 ) {
        std::cout<< "i = "<< i << " angle : " << (180/TMath::Pi())*mean_direction[i].Angle(mean_direction[i+1]) << std::endl;
  }
    angle_distribution.push_back( mean_direction[i].Angle(mean_direction[i+1]) ) ;
  }
  //std::cout<< " angle test = " << (180/TMath::Pi())*test1.Angle(test2)<<std::endl;
  return angle_distribution;
}
