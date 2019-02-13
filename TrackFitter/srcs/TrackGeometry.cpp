#include "../include/TrackFitter.h"
#include <iostream>
#include "math.h"

//*/
/**
* FUNCTIONS
* 1 - Get Properties
*/

unsigned int TrackFitter::GetHits( const unsigned int & event_id_track ){
  return _event_hits[event_id_track-1] ;
}

Hit_level TrackFitter::AccessVertex( const unsigned int & event_id_track ) {
  return _event_tracks[event_id_track-1][0] ;
}

Hit_level TrackFitter::AccessEnd( const unsigned int & event_id_track ) {
  return _event_tracks[event_id_track-1][_event_hits[event_id_track]-1] ;
}

std::vector< double > TrackFitter::GetdEdx( const unsigned int & event_id_track ){
  return _event_tracks[event_id_track-1][3];
}

Track TrackFitter::GetTrack( const unsigned int & event_id_track ){
  return _event_tracks[event_id_track-1] ;
}

bool TrackFitter::IsContained( const unsigned int & event_id_track ){
  bool contained = true ;
  for ( int i = 0 ; i < _has_reco_daughters[event_id_track-1] + 1 ; ++i ){
    if( _event_pfps_vcontained[event_id_track-1][i] == 1 && _event_pfps_econtained[event_id_track-1][i] == 1 ) { contained = true ; }
    else contained = false ;
  }
  return contained ;
}
double TrackFitter::FindMaxCoordinate( const unsigned int & event_id_track , const int & coordinate_id ) {
  double max_coordinate = -600 ;
  for( int i = 0 ; i < _event_hits[event_id_track-1] ; ++i ){
    if ( _event_tracks[event_id_track-1][coordinate_id][i] > max_coordinate ) max_coordinate = _event_tracks[event_id_track-1][coordinate_id][i] ;
  }
  return max_coordinate ;
}

double TrackFitter::FindMinCoordinate( const unsigned int & event_id_track , const int & coordinate_id ) {
  double min_coordinate = 600 ;
  for( int i = 0 ; i < _event_hits[event_id_track-1] ; ++i ){
    if ( _event_tracks[event_id_track-1][coordinate_id][i] < min_coordinate ) min_coordinate = _event_tracks[event_id_track-1][coordinate_id][i] ;
  }
  return min_coordinate ;
}


std::map< int, std::vector<int> > TrackFitter::FindHiearchy( const unsigned int & event_track_id ) {
  int primary_id = -999 ;
  std::map< int, std::vector<int> > map_hiearchy ;
  std::vector<int> daughter;
  double dist_MC_vertex, dist_MC_end ;
  double min_dist_vertex = 99999 ;

  for( int i = 0 ; i < _has_reco_daughters[event_track_id-1] + 1 ; ++i ){
    // if just one track, take that as primary
    if( _has_reco_daughters[event_track_id-1] == 0 ) { primary_id = 0 ; }
    // if more than one track, find primary true vertex :
    // --> Joins are reconstructed as neutrino vertex. So the vertices may be the end of the real primary track, and the end may be the real vertex
    dist_MC_vertex =  pow(_event_MC_vertex[event_track_id-1][0] - _event_pfps_end_x[event_track_id-1][i] , 2 ) ;
    dist_MC_vertex += pow(_event_MC_vertex[event_track_id-1][1] - _event_pfps_end_y[event_track_id-1][i] , 2 ) ;
    dist_MC_vertex += pow(_event_MC_vertex[event_track_id-1][2] - _event_pfps_end_z[event_track_id-1][i] , 2 ) ;
    dist_MC_vertex =  sqrt( dist_MC_vertex ) ;

    if( dist_MC_vertex < min_dist_vertex ) {
      min_dist_vertex = dist_MC_vertex ;
      primary_id = i ;
    }
  }

  for( int i = 0 ; i < _has_reco_daughters[event_track_id-1] + 1 ; ++i ){
    if( i == primary_id ) continue ;
    // Look for daughters :
    dist_MC_end =  pow( _event_pfps_start_x[event_track_id-1][primary_id] - _event_pfps_start_x[event_track_id-1][i] , 2 ) ;
    dist_MC_end += pow( _event_pfps_start_y[event_track_id-1][primary_id] - _event_pfps_start_y[event_track_id-1][i] , 2 ) ;
    dist_MC_end += pow( _event_pfps_start_z[event_track_id-1][primary_id] - _event_pfps_start_z[event_track_id-1][i] , 2 ) ;
    dist_MC_end =  sqrt( dist_MC_end ) ;

    // some events seem to have the right order!! Ex 6
    // std::cout<< dist_MC_end << " dist x " << _event_pfps_start_x[event_track_id-1][i] << " dist y " << _event_pfps_start_y[event_track_id-1][i]<< " dist z " << _event_pfps_start_z[event_track_id-1][i] << std::endl;
    // std::cout<< dist_MC_end << " dist x " << _event_pfps_end_x[event_track_id-1][primary_id] << " dist y " << _event_pfps_end_y[event_track_id-1][primary_id]<< " dist z " << _event_pfps_end_z[event_track_id-1][primary_id] << std::endl;
    // std::cout<< _event_pfps_start_x[event_track_id-1][primary_id] - _event_pfps_start_x[event_track_id-1][i] << " , " << _event_pfps_start_y[event_track_id-1][primary_id] - _event_pfps_start_y[event_track_id-1][i] << " , " << _event_pfps_start_z[event_track_id-1][primary_id] - _event_pfps_start_z[event_track_id-1][i] << std::endl;
    if( dist_MC_end < 5 ) { // if min < cut lenght push back into vector. Cut out photons
      daughter.push_back( i ) ;
    }
    // std::cout<< " i = " << i << " track id = " << event_track_id -1 << "dist vertex = " << dist_MC_end  << " sixe = " << daughter.size() << std::endl ;
  ///     std::cout<< " min dist end " << min_dist_end  << " id - " << secondary_id << std::endl;

  }
  map_hiearchy[primary_id] = daughter ;
  return map_hiearchy ;
}
