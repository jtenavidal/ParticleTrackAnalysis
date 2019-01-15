#ifndef TRACK_H
#define TRACK_H

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TBranch.h"

// Typedef for the map
typedef std::vector< std::vector<double> > Track;
typedef std::vector< double > Hit_level;
/**
 * @brief  Event class
 */
 class TrackFitter{
    public :
    /**
    * Constructor
    */
    TrackFitter( const Track & p_track ); // reads from easier and simpler file
    TrackFitter( const std::string & track_file_path ) ; // reads all events from track_file_path
    // -> Need to find a way to adapt the functions to this matrix!
    TrackFitter( const std::string & track_file_path, const unsigned int & event_id_track ) ; // reads only the event corresponding to event_id_track

    /**
    * Functions to get properties
    */
    unsigned int GetHits( ) ;
    Hit_level AccessVertex( ) ;
    Hit_level AccessEnd( ) ;
    std::vector< float > GetdQdx( ) ;
    Track GetTrack( ) ;

    /**
    * Functions to check and save information
    */
    void SaveTrack( const std::string & path ) const ;
    void PrintHipotesis( const std::string & path ) ;
    void PrintdQdx( const std::string & path ) const ;
    void PlotLinearityTrack( const int & window, const std::string & path ) ;
    void PlotLinearityData( const int & window, const std::string & path );//, const std::vector< double > & Data_1, const std::vector< double > & Data_2  ) ;

    /**
    * Functions to guess about the geometry of the track
    */
    double FitToLine( ) ;

  private :

  // Object information
  int _hits, _hits_size, _dQdx_size ;
  Hit_level _vertex_position, _end_position ;
  std::vector< float > _reco_dQdx ;
  Track _particle_track ;
  // If loading all file
  std::vector< int > _event_hits, _event_dQdx_size ;
  std::vector< std::vector<double> > _event_vertex, _event_end ;
  std::vector< std::vector< float > > _event_reco_dQdx ;
  std::vector< Track > _event_tracks;  // maps event and track
  /**
  * Functions to guess about the geometry of the track
  */
  Track Straight( ) ;


  /**
  * Functions to guess about the geometry of the track
  */

  std::vector< double > MeanData( const int & window, const std::vector< double > & data ) ;
  std::vector< double > DevData( const int & window, const std::vector< double > & dat ) ;
  std::vector< double > CovData( const int & window, const std::vector< double > & Data_1, const std::vector< double > & Data_2 ) ;
  std::vector< double > LinearityData( const int & window, const std::vector< double > & Data_1, const std::vector< double > & Data_2 ) ;

  }; // Event


#endif
