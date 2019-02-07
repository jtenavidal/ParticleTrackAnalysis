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
