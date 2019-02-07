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
