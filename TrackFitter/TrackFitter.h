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
    TrackFitter( const Track & p_track );
    TrackFitter( const std::string & track_file_path ) ;

    /**
    * Functions to get properties
    */
    unsigned int GetHits( ) ;
    Hit_level AccessVertex( ) ;
    Hit_level AccessEnd( ) ;
    Hit_level GetdEdx( ) ;
    Hit_level GetdQdx( ) ;
    Track GetTrack( ) ;

    /**
    * Functions to check and save information
    */
    void SaveTrack( const std::string & path ) const ;
    void PrintHipotesis( const std::string & path ) ;
    void PrintdEdx( const std::string & path ) const ;
    void PlotLinearity( const int & window, const std::string & path ) ;

    /**
    * Functions to guess about the geometry of the track
    */
    double FitToLine( ) ;
    std::vector< std::vector< double > > Linearity( const int & window ) ;
    std::vector< std::vector< double > > MeanPosition( const int & window ) ;
    std::vector< std::vector< double > > DevPosition( const int & window ) ;
    std::vector< std::vector< double > > CovPosition( const int & window ) ;

  private :

  // Object information
  unsigned int _hits ;
  Hit_level _vertex_position, _end_position, _reco_dEdx, _reco_dQdx ;
  Track _particle_track ;
  /**
  * Functions to guess about the geometry of the track
  */
  Track Straight( ) ;

  }; // Event


#endif
