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

void TrackFitter::TruthParticles( unsigned int & event_id_track ){
  if ( event_id_track == 0 ) event_id_track = 1 ;
  std::cout<< "Primary particle PDG = " << _TPDG_Code_Primary[ event_id_track-1 ] << std::endl;
  std::cout<< "#mu = " << _Tnu_mu[event_id_track-1] << " #pi = " << _Tnu_pi[event_id_track-1] << " #p = " << _Tnu_p[event_id_track-1] << " #e = " << _Tnu_e[event_id_track-1]<< " #photon = " << _Tnu_photon[event_id_track-1] << " #others = " << _Tnu_others[event_id_track-1]<< std::endl;
}


/**
* FUNCTIONS
* 1 - Get Properties
* 2 - To check or save information
*/

void TrackFitter::SaveTrack( std::string const & path , const unsigned int & event_id_track ) {
  std::vector< std::vector< double > > min_Linearity_position = FindMinimumLinearityPosition( event_id_track ) ;
  int bins = int(_event_hits[event_id_track-1]/10) ;
  int window =  int( _event_hits[event_id_track-1] * 0.07 ) ;
  if( window < 5 ) window = 5 ;

  TH3D *h_track = new TH3D("h_track", " Particle Track ", bins,
 _event_tracks[event_id_track-1][0][0], _event_tracks[event_id_track-1][0][_event_hits[event_id_track-1]-1], bins,
_event_tracks[event_id_track-1][1][0], _event_tracks[event_id_track-1][1][_event_hits[event_id_track-1]-1], bins, // need to define number of bins as a function of _hits to avoid bad memory allocation
_event_tracks[event_id_track-1][2][0], _event_tracks[event_id_track-1][2][_event_hits[event_id_track-1]-1] );
  std::cout<< " x = " << _event_tracks[event_id_track-1][0][0] << " Xf= "<<  _event_tracks[event_id_track-1][0][_event_hits[event_id_track-1]-2] << " hits " << _event_hits[event_id_track-1] << std::endl;
  TH3D *h_track_kink = (TH3D * ) h_track ->Clone();
  h_track_kink->SetName("h_track_kink") ;

  TLegend *leg = new TLegend(0.9,0.7,0.48,0.9);

  for( int i = 0; i < _event_hits[event_id_track-1]; ++i ){
    h_track-> Fill(_event_tracks[event_id_track-1][0][i], _event_tracks[event_id_track-1][1][i], _event_tracks[event_id_track-1][2][i] ) ; }//, _event_tracks[event_id_track-1][3][i]);  }

  TCanvas *c = new TCanvas();
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
  h_track->SetLineColor(2);
  h_track->GetXaxis()->SetTitle("X");
  h_track->GetYaxis()->SetTitle("Y");
  h_track->GetZaxis()->SetTitle("Z");
  leg->AddEntry( h_track, " Track 3D trajectory ");
  h_track->Draw("BOX2Z");

  if( min_Linearity_position.size() > 0 ) { // can remove
    for( unsigned int i = 0 ; i < min_Linearity_position.size() ; ++i ) {
      h_track_kink -> Fill( min_Linearity_position[i][0], min_Linearity_position[i][1], min_Linearity_position[i][2] ) ;
    }
    h_track_kink->SetLineColor(3) ;
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

// STATISTICS GENERAL FUNCTIONS
std::vector< double > TrackFitter::MeanData( const std::vector<double> & data ){
  unsigned int starting_hit , end_hit;
  double muI_data;
  std::vector< double > mu_data ;
  int window =  int( data.size() * 0.05 ) ;
  if( window < 5 ) window = 5 ;

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

std::vector< double > TrackFitter::DevData( const std::vector<double> & data ){
  unsigned int starting_hit , end_hit;
  double devI_data;
  std::vector< double > dev_data , mean ;
  int window =  int( data.size() * 0.05 ) ;
  if( window < 5 ) window = 5 ;
  mean = MeanData( data ) ;

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

std::vector< double > TrackFitter::CovData( const std::vector< double > & Data_1, const std::vector< double > & Data_2 ){
  unsigned int starting_hit , end_hit;
  double covI_12 ; // 1 - variable 1, 2 - second variable
  std::vector< double > cov_12 , mean1, mean2;
  mean1 = MeanData( Data_1 ) ; //variable 1
  mean2 = MeanData( Data_2 ) ; //variable 2
  int window =  int( Data_1.size() * 0.05 ) ;
  if( window < 5 ) window = 5 ;


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

std::vector< double > TrackFitter::LinearityData( const std::vector< double > & Data_1, const std::vector< double > & Data_2 ){
  // It calculates the linearity from Pearson correlation coefficient (corrP)
  unsigned int starting_hit , end_hit;
  std::vector< double > corrP_12, cov_12, dev1, dev2 ;
  cov_12 = CovData( Data_1, Data_2 ) ;
  dev1 = DevData( Data_1 ) ;
  dev2 = DevData( Data_2 ) ;
  int window =  int( Data_1.size() * 0.05 ) ;
  if( window < 5 ) window = 5 ;


  if ( Data_1.size() == Data_2.size()) {
    for( int i = 0; i < int( Data_1.size() ); ++i ){
      corrP_12.push_back( sqrt(std::pow(cov_12[i]/(dev1[i]*dev2[i]),2)) ) ;
    }
  }
  return corrP_12;
}

std::vector< TVector3 > TrackFitter::MeanDirectionData( const unsigned int & event_id_track ){
  // this will only be applied to the track information (x,y,z)
  unsigned int starting_hit , end_hit;
  TVector3 directionI ;
  std::vector< TVector3 > mean_direction ;
  int window =  int( _event_hits[event_id_track-1] * 0.07 ) ;
  if( window < 5 ) window = 5 ;

  for( int i = 0; i < int(_event_hits[event_id_track-1]); ++i ){
    if ( i - window < 0 ) {
      starting_hit = 0 ;
      end_hit = i + window ;
    } else if ( i + window >= int(_event_hits[event_id_track-1])) {
      starting_hit = i - window ;
      end_hit = _hits - 1 ;
    } else if ( i - window < 0  && i + window >= int(_event_hits[event_id_track-1])) {
      starting_hit = 0 ;
      end_hit = _event_hits[event_id_track-1] - 1 ;
    } else {
      starting_hit = i - window ;
      end_hit = i + window ;
    }
    directionI.SetX( (_event_tracks[event_id_track-1][0][end_hit]-_event_tracks[event_id_track-1][0][starting_hit])/(end_hit-starting_hit) );
    directionI.SetY( (_event_tracks[event_id_track-1][1][end_hit]-_event_tracks[event_id_track-1][1][starting_hit])/(end_hit-starting_hit) );
    directionI.SetZ( (_event_tracks[event_id_track-1][2][end_hit]-_event_tracks[event_id_track-1][2][starting_hit])/(end_hit-starting_hit) );
    mean_direction.push_back( directionI );
  }
  return mean_direction;
}


std::vector< double > TrackFitter::AngleTrackDistribution( const unsigned int & event_id_track ) {
  // this will only be applied to the track information (x,y,z)
  unsigned int starting_hit , end_hit;
  std::vector< double > angle_distribution ;
  std::vector< TVector3 > mean_direction = MeanDirectionData( event_id_track ) ;
  TVector3 test1(1, 1, 0);
  TVector3 test2(-1, -1, 0);
  for( int i = 0; i < int(_event_hits[event_id_track-1]) - 1 ; ++i ){

  if( mean_direction[i].Angle(mean_direction[i+1]) > 1 ) {
        std::cout<< "i = "<< i << " angle : " << (180/TMath::Pi())*mean_direction[i].Angle(mean_direction[i+1]) << std::endl;
  }
    angle_distribution.push_back( mean_direction[i].Angle(mean_direction[i+1]) ) ;
  }
  std::cout<< " angle test = " << (180/TMath::Pi())*test1.Angle(test2)<<std::endl;
  return angle_distribution;
}

std::vector< std::vector< double > > TrackFitter::FindMinimumLinearityPosition( const unsigned int & event_id_track ){
    std::vector< std::vector< double > > min_Linearity_position ;
    std::vector< double > position ;
    std::vector< double > corrP_XY = LinearityData( _event_tracks[event_id_track-1][0], _event_tracks[event_id_track-1][1] ) ;
    std::vector< double > corrP_XZ = LinearityData( _event_tracks[event_id_track-1][0], _event_tracks[event_id_track-1][2] ) ;
    std::vector< double > corrP_YZ = LinearityData( _event_tracks[event_id_track-1][2], _event_tracks[event_id_track-1][1] ) ;
    std::vector< double > corrP ;
    int window =  int( _event_hits[event_id_track-1] * 0.07 ) ;
    if( window < 5 ) window = 5 ;
    double linearity_min = 2 ;
    unsigned int min_hit = 0 ;

    for( unsigned int i = 0 ; i < corrP_XY.size() ; ++i ) {
        corrP.push_back(corrP_XY[i]*corrP_XZ[i]*corrP_YZ[i] ) ;
    }

    for( unsigned int i = 0 ; i < corrP_XY.size() - 1 ; ++i ) {
        if( corrP[i+1] < corrP[i] && corrP[i] < linearity_min ) {
          linearity_min = corrP[i+1] ;
          min_hit = i+1 ;
        }

        if( linearity_min < 0.9 && i == min_hit + int( window / 2 ) && corrP[i] > corrP[ min_hit ] ) {
          // reseting : looking for other minima
          position.push_back( _event_tracks[event_id_track-1][0][min_hit] ) ;
          position.push_back( _event_tracks[event_id_track-1][1][min_hit] ) ;
          position.push_back( _event_tracks[event_id_track-1][2][min_hit] ) ;
          min_Linearity_position.push_back( position ) ;
          // std::cout<< " x hit  = " << _event_tracks[event_id_track-1][0][min_hit] << " y hit = " << _event_tracks[event_id_track-1][1][min_hit] << " z hit = " << _event_tracks[event_id_track-1][2][min_hit] << std::endl;

          position.clear();
          linearity_min = 2 ;
        }
    }

    return min_Linearity_position;
}

void TrackFitter::StatisticsKinks( ){
  int has_kink = 0 , is_straight = 0 , has_one_kink = 0 , has_more_kinks = 0 ;
  for ( unsigned int i = 1 ; i < _event_tracks.size() +1 ; ++i ){
    if ( FindMinimumLinearityPosition( i ).size() > 0 ) {
      ++has_kink ;
      if ( FindMinimumLinearityPosition( i ).size() == 1 ) {
          ++has_one_kink ;
        }
      if ( FindMinimumLinearityPosition( i ).size() > 1 ) {
          ++has_more_kinks ;
          std::cout<< " Event ID = " << i << std::endl;
          SaveTrack( ("Track_testing_new_constructor_"+std::to_string(i)) , i ) ;
          unsigned int event = i ;
          PlotLinearityTrack( " testing "+std::to_string(i) , event );
        }
    }
    else { ++is_straight ; }
  }

  std::cout<< " % kinked tracks = " << has_kink*100/(has_kink+is_straight) << std::endl;
  std::cout<< " --->   % 1  kinked tracks = " << has_one_kink*100/has_kink << std::endl;
  std::cout<< " kinq 1 " << has_one_kink << std::endl;

  std::cout<< " kinq 1 " << has_more_kinks << std::endl;
  std::cout<< " --->   % >1 kinked tracks = " << has_more_kinks*100/has_kink << std::endl;
  std::cout<< " % straight tracks = " << is_straight*100/(has_kink+is_straight) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
// Truth information functions :                                             //
///////////////////////////////////////////////////////////////////////////////
void TrackFitter::SaveStatisticsTrueEvent( const std::string & path ) {
/**
  * -> Badruns muons and pions
  * -> Muons with Michel e
  * -> Muons no Michel e
  * -> Breakdown muons and other particles
  * -> particle and number of kinks
  * -> Lenght primary particle
  * -> Ratio lenght primary / 1-daughter
  */
  TCanvas *c = new TCanvas();
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
  TLegend *leg = new TLegend(0.9,0.7,0.48,0.9);
  THStack * h_mu_stats = new THStack( "h_mu_stats", "Muon Statistics" ) ;
  TH1D *h_mu_rdaughter = new TH1D("h_mu_rdaughter", " Muon Statistics: reconstructed daughters ", 3, 0, 2 );
  TH1D *h_mu_charged = new TH1D("h_mu_charged", " Muon Statistics: Number of total charged daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_e = new TH1D("h_mu_charged_e", " Muon Statistics: Number of total e daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_mu = new TH1D("h_mu_charged_mu", " Muon Statistics: Number of total mu daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_pi = new TH1D("h_mu_charged_pi", " Muon Statistics: Number of total pi daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_p = new TH1D("h_mu_charged_p", " Muon Statistics: Number of total p daughters ", 3, 0, 2 );

  for( unsigned int j = 0 ; j < _rnu_daughters.size() ; ++j ){
    if( _TPDG_Code_Primary[j] != 13 ) continue ;
    if( _rnu_daughters[j] == 0 ) h_mu_rdaughter->Fill(0) ;
    if( _rnu_daughters[j] == 1 ) h_mu_rdaughter->Fill(1) ;
    if( _rnu_daughters[j] > 1 ) h_mu_rdaughter->Fill(2) ; // needs breakdown to understand
  }

  for( unsigned int j = 0 ; j < _TPDG_Code_Primary.size() ; ++j ){
    if ( _TPDG_Code_Primary[j] != 13 ) continue ;
    if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(0) ;
    } else if ( _Tnu_pi[j] == 1  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(1) ; h_mu_charged_pi->Fill(1) ;
    } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 1  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(1) ; h_mu_charged_mu->Fill(1) ;
    } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 1  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(1) ; h_mu_charged_e->Fill(1) ;
    } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 1  ) { h_mu_charged->Fill(1) ; h_mu_charged_p->Fill(1) ;
    } else { h_mu_charged->Fill(2) ; }
  }

  h_mu_charged_mu->SetFillColor(kRed);
  h_mu_charged_pi->SetFillColor(kBlue);
  h_mu_charged_e->SetFillColor(kGreen);
  h_mu_charged_p->SetFillColor(kYellow);
  h_mu_charged_mu->SetFillStyle(3004);
  h_mu_charged_pi->SetFillStyle(3004);
  h_mu_charged_e->SetFillStyle(3004);
  h_mu_charged_p->SetFillStyle(3004);
  leg->AddEntry(h_mu_charged_mu, "True muon daughters ") ;
  leg->AddEntry(h_mu_charged_pi, "True pion daughters ") ;
  leg->AddEntry(h_mu_charged_e, "True electron daughters ") ;
  leg->AddEntry(h_mu_charged_p, "True proton daughters ") ;
  h_mu_charged_mu->Scale(1/h_mu_charged->GetEntries());
  h_mu_charged_pi->Scale(1/h_mu_charged->GetEntries());
  h_mu_charged_e->Scale(1/h_mu_charged->GetEntries());
  h_mu_charged_p->Scale(1/h_mu_charged->GetEntries());
  h_mu_stats->Add(h_mu_charged_mu) ;
  h_mu_stats->Add(h_mu_charged_pi) ;
  h_mu_stats->Add(h_mu_charged_e) ;
  h_mu_stats->Add(h_mu_charged_p) ;
  h_mu_stats->Draw("hist") ;
  /////////////////////////////////////////////////////////////////////////////
  h_mu_rdaughter->SetLineColor(kPink);
  h_mu_rdaughter->SetLineWidth(2);
  h_mu_rdaughter->SetFillStyle(3003);
  h_mu_rdaughter->SetFillColor(kPink);
  h_mu_rdaughter->GetXaxis()->SetTitle("Number reconstructed particles");
  h_mu_rdaughter->GetYaxis()->SetTitle("%");
  h_mu_rdaughter->Scale(1/h_mu_rdaughter->GetEntries());
  h_mu_rdaughter->GetYaxis()->SetRangeUser(0.,1.);
  leg->AddEntry(h_mu_rdaughter, "Reconstructed secondary tracks ") ;
  h_mu_rdaughter->Draw("hist SAME");

  leg->Draw();
  c->SaveAs((path+"_muon_daughters_study.root").c_str());
  leg->Clear();
  c->Clear();

  /*c->Clear();

  TH1D *h_mu_Tdaughter = new TH1D("h_mu_Tdaughter", " Muon Statistics: true daughters ", 3, 0, 2 );
  for( unsigned int j = 0 ; j < _TPDG_Code_Primary.size() ; ++j ){
    if( _TPDG_Code_Primary[j] != 13 ) continue ;
    if( _Tnu_daughters[j] == 0 ) h_mu_Tdaughter->Fill(0) ;
    if( _Tnu_daughters[j] == 1 ) h_mu_Tdaughter->Fill(1) ;
    if( _Tnu_daughters[j] > 1 ) h_mu_Tdaughter->Fill(2) ; // needs breakdown to understand
  }

  h_mu_Tdaughter->SetLineColor(2);
  h_mu_Tdaughter->GetXaxis()->SetTitle("Number daughter electrons");
  h_mu_Tdaughter->GetYaxis()->SetTitle("%");
//  h_mu_Tdaughter->Scale(1/h_mu_Tdaughter->GetEntries());
  h_mu_Tdaughter->Draw("hist");
  c->SaveAs((path+"_muon_rdaughter.root").c_str());
  //  c->Clear();

*/
  // Histogram length second track

// HISTOGRAM WIHT NUMBER DAUGHTERS FOR PION> ADD ON TOP BREAKDOWN . FIX TTREE
THStack * h_pi_stats = new THStack( "h_pi_stats", "Muon Statistics" ) ;
TH1D *h_pi_rdaughter = new TH1D("h_pi_rdaughter", " Muon Statistics: reconstructed daughters ", 3, 0, 2 );
TH1D *h_pi_charged = new TH1D("h_pi_charged", " Muon Statistics: Number of total charged daughters ", 3, 0, 2 );
TH1D *h_pi_charged_e = new TH1D("h_pi_charged_e", " Muon Statistics: Number of total e daughters ", 3, 0, 2 );
TH1D *h_pi_charged_mu = new TH1D("h_pi_charged_mu", " Muon Statistics: Number of total mu daughters ", 3, 0, 2 );
TH1D *h_pi_charged_pi = new TH1D("h_pi_charged_pi", " Muon Statistics: Number of total pi daughters ", 3, 0, 2 );
TH1D *h_pi_charged_p = new TH1D("h_pi_charged_p", " Muon Statistics: Number of total p daughters ", 3, 0, 2 );

for( unsigned int j = 0 ; j < _rnu_daughters.size() ; ++j ){
  if( _TPDG_Code_Primary[j] != 211 ) continue ;
  if( _rnu_daughters[j] == 0 ) h_pi_rdaughter->Fill(0) ;
  if( _rnu_daughters[j] == 1 ) h_pi_rdaughter->Fill(1) ;
  if( _rnu_daughters[j] > 1 ) h_pi_rdaughter->Fill(2) ; // needs breakdown to understand
}

for( unsigned int j = 0 ; j < _TPDG_Code_Primary.size() ; ++j ){
  if ( _TPDG_Code_Primary[j] != 211 ) continue ;
  if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(0) ; std::cout<< "CUU" << std::endl;
  } else if ( _Tnu_pi[j] == 1  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(1) ; h_pi_charged_pi->Fill(1) ;
  } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 1  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(1) ; h_pi_charged_mu->Fill(1) ;
  } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 1  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(1) ; h_pi_charged_e->Fill(1) ;
  } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 1  ) { h_pi_charged->Fill(1) ; h_pi_charged_p->Fill(1) ;
  } else { h_pi_charged->Fill(2) ; }
}
///////////////////////////////////////////////////////////////////////////////
h_pi_charged->SetLineColor(1);
leg->AddEntry(h_pi_stats, "Total charged particles") ;
h_pi_charged->Draw("hist SAME") ;
/////////////////////////////////////////////////////////////////////////////
h_pi_charged_mu->SetFillColor(kRed);
h_pi_charged_pi->SetFillColor(kBlue);
h_pi_charged_e->SetFillColor(kGreen);
h_pi_charged_p->SetFillColor(kYellow);
h_pi_charged_mu->SetFillStyle(3004);
h_pi_charged_pi->SetFillStyle(3004);
h_pi_charged_e->SetFillStyle(3004);
h_pi_charged_p->SetFillStyle(3004);
leg->AddEntry(h_pi_charged_mu, "True muon daughters ") ;
leg->AddEntry(h_pi_charged_pi, "True pion daughters ") ;
leg->AddEntry(h_pi_charged_e, "True electron daughters ") ;
leg->AddEntry(h_pi_charged_p, "True proton daughters ") ;
//h_pi_charged_mu->Scale(1/h_pi_charged->GetEntries());
//  h_pi_charged_pi->Scale(1/h_pi_charged->GetEntries());
//h_pi_charged_e->Scale(1/h_pi_charged->GetEntries());
//h_pi_charged_p->Scale(1/h_pi_charged->GetEntries());
h_pi_stats->Add(h_pi_charged_mu) ;
h_pi_stats->Add(h_pi_charged_pi) ;
h_pi_stats->Add(h_pi_charged_e) ;
h_pi_stats->Add(h_pi_charged_p) ;
h_pi_stats->Draw("hist SAME") ;
/////////////////////////////////////////////////////////////////////////////
h_pi_rdaughter->SetLineColor(kPink);
h_pi_rdaughter->SetLineWidth(2);
h_pi_rdaughter->SetFillStyle(3003);
h_pi_rdaughter->SetFillColor(kPink);
h_pi_rdaughter->GetXaxis()->SetTitle("Number reconstructed particles");
h_pi_rdaughter->GetYaxis()->SetTitle("%");
//h_pi_rdaughter->Scale(1/h_pi_rdaughter->GetEntries());
//h_pi_rdaughter->GetYaxis()->SetRangeUser(0.,1.);
leg->AddEntry(h_pi_rdaughter, "Reconstructed secondary tracks ") ;
h_pi_rdaughter->Draw("hist SAME");
leg->Draw();
c->SaveAs((path+"_pion_daughters_study.root").c_str());
//c->Clear();
}
