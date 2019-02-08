#include "../include/TrackFitter.h"
#include <iostream>


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
