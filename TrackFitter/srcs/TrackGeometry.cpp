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

unsigned int TrackFitter::GetDaughters( const unsigned int & event_id_track ){
  return _has_reco_daughters[event_id_track-1] ;
}

float TrackFitter::GetTotalRecoLength( const unsigned int & event_id_track ){
  return _event_RLenght[event_id_track-1] ;
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

bool TrackFitter::TruePrimaryVertexContained( const unsigned int & event_id_track ){
  bool contained = true ;
  if( _TPrimary_vcontained[event_id_track-1] == 1 ) { contained = true ; }
  else contained = false ;
  return contained ;
}

bool TrackFitter::TruePrimaryEndContained( const unsigned int & event_id_track ){
  bool contained = true ;
  if( _TPrimary_econtained[event_id_track-1] == 1 ) { contained = true ; }
  else contained = false ;
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
  std::vector<int> daughter, daughter2 ;
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

  // std::cout << " primary is - " << primary_id << std::endl;

  for( int i = 0 ; i < _has_reco_daughters[event_track_id-1] + 1 ; ++i ){
    if( i == primary_id ) continue ;
    // Look for daughters :
    // std::cout<< " **** i ***** " << i  << " vs " << primary_id << std::endl;
    // std::cout<< " x - " << _event_pfps_start_x[event_track_id-1][i] << " end " << _event_pfps_end_x[event_track_id-1][i]  << std::endl;
    // std::cout<< " x - primary" << _event_pfps_start_x[event_track_id-1][primary_id] << " end " << _event_pfps_end_x[event_track_id-1][primary_id]  << std::endl;
    // std::cout<< " y - " << _event_pfps_start_y[event_track_id-1][i] << " end " << _event_pfps_end_y[event_track_id-1][i]  << std::endl;
    // std::cout<< " y - primary" << _event_pfps_start_y[event_track_id-1][primary_id] << " end " << _event_pfps_end_y[event_track_id-1][primary_id]  << std::endl;
    // std::cout<< " z - " << _event_pfps_start_z[event_track_id-1][i] << " end " << _event_pfps_end_z[event_track_id-1][i]  << std::endl;
    // std::cout<< " z - primary" << _event_pfps_start_z[event_track_id-1][primary_id] << " end " << _event_pfps_end_z[event_track_id-1][primary_id]  << std::endl;
    dist_MC_end =  pow( _event_pfps_start_x[event_track_id-1][primary_id] - _event_pfps_start_x[event_track_id-1][i] , 2 ) ;
    dist_MC_end += pow( _event_pfps_start_y[event_track_id-1][primary_id] - _event_pfps_start_y[event_track_id-1][i] , 2 ) ;
    dist_MC_end += pow( _event_pfps_start_z[event_track_id-1][primary_id] - _event_pfps_start_z[event_track_id-1][i] , 2 ) ;
    dist_MC_end =  sqrt( dist_MC_end ) ;
    // std::cout<< " dist_MC_end "<< dist_MC_end << std::endl;
    // some events seem to have the right order! If cut not satisfied check if that is the case.
    if( dist_MC_end < 3 ) { // if min < cut lenght push back into vector. Cut out photons
      daughter.push_back( i ) ; //std::cout<< " secondary = " << i << std::endl;
      continue ;
    } else {
      dist_MC_end =  pow( _event_pfps_end_x[event_track_id-1][primary_id] - _event_pfps_start_x[event_track_id-1][i] , 2 ) ;
      dist_MC_end += pow( _event_pfps_end_y[event_track_id-1][primary_id] - _event_pfps_start_y[event_track_id-1][i] , 2 ) ;
      dist_MC_end += pow( _event_pfps_end_z[event_track_id-1][primary_id] - _event_pfps_start_z[event_track_id-1][i] , 2 ) ;
      dist_MC_end =  sqrt( dist_MC_end ) ;
      if( dist_MC_end < 3 ) { // if min < cut lenght push back into vector. Cut out photons. Check if this is the case ( look for example )
        daughter.push_back( i ) ; //std::cout<< " secondary = " << i << std::endl; std::cout<< " x - " << _event_pfps_start_x[event_track_id-1][i] << " end " << _event_pfps_end_x[event_track_id-1][i]  << std::endl;
      // For some cases it seems to be empty -> cannot find mother.
      } else {
        for ( int j = 0 ; j < _has_reco_daughters[event_track_id-1] + 1 ; ++j ) {
          if( _event_pfps_start_x[event_track_id-1][i] == 0 && _event_pfps_start_x[event_track_id-1][i] ==0 && _event_pfps_start_x[event_track_id-1][i] == 0 ) continue ;
          if( j == primary_id || i == j ) continue ;
          dist_MC_end =  pow( _event_pfps_start_x[event_track_id-1][i] - _event_pfps_end_x[event_track_id-1][j] , 2 ) ;
          dist_MC_end += pow( _event_pfps_start_y[event_track_id-1][i] - _event_pfps_end_y[event_track_id-1][j] , 2 ) ;
          dist_MC_end += pow( _event_pfps_start_z[event_track_id-1][i] - _event_pfps_end_z[event_track_id-1][j] , 2 ) ;
          dist_MC_end =  sqrt( dist_MC_end ) ;
          if( dist_MC_end < 3 ) { // if min < cut lenght push back into vector. Cut out photons
            daughter2.push_back( i );
            map_hiearchy[j] = daughter2;
            // std::cout<< " found daughter " << std::endl;
          } else {
            dist_MC_end =  pow( _event_pfps_start_x[event_track_id-1][i] - _event_pfps_start_x[event_track_id-1][j] , 2 ) ;
            dist_MC_end += pow( _event_pfps_start_y[event_track_id-1][i] - _event_pfps_start_y[event_track_id-1][j] , 2 ) ;
            dist_MC_end += pow( _event_pfps_start_z[event_track_id-1][i] - _event_pfps_start_z[event_track_id-1][j] , 2 ) ;
            dist_MC_end =  sqrt( dist_MC_end ) ;

            if( dist_MC_end < 3 ) { // if min < cut lenght push back into vector. Cut out photons
              daughter2.push_back( i );
              map_hiearchy[j] = daughter2;
              // std::cout<< " found daughter: " << i << " mother is " << j  << std::endl;
            }
          }
        }
      }
    }
  }
  map_hiearchy[primary_id] = daughter ;
  return map_hiearchy ;
}
