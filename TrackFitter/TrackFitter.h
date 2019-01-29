#ifndef TRACK_H
#define TRACK_H

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TVector3.h"

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

    /**
    * Functions to get properties
    */
    // For constructor 1
    unsigned int GetHits( const unsigned int & event_id_track ) ;
    void TruthParticles( unsigned int & event_id_track ) ;
    Hit_level AccessVertex( const unsigned int & event_id_track ) ;
    Hit_level AccessEnd( const unsigned int & event_id_track ) ;
    std::vector< double > GetdQdx( const unsigned int & event_id_track ) ;
    Track GetTrack( const unsigned int & event_id_track ) ;

    /**
    * Functions to check and save information
    */
    void SaveTrack( const std::string & path , const unsigned int & event_id_track ) const ;
    void PrintdQdx( const std::string & path , const unsigned int & event_id_track ) const ;
    void SaveTrack( const std::string & path ) const ;
    //void PrintHipotesis( const std::string & path , const unsigned int & event_id_track ) const ) ;
    void PlotLinearityTrack( const int & window, const std::string & path , const unsigned int & event_id_track ) ;
    void PlotLinearityData( const int & window, const std::string & path , const unsigned int & event_id_track );//, const std::vector< double > & Data_1, const std::vector< double > & Data_2  ) ;
    void SaveStatisticsTrueEvent( const std::string & path );
    /**
    * Functions to guess about the geometry of the track
    */
    //double FitToLine( ) ;
    std::vector< TVector3 > MeanDirectionData( const int & window , const unsigned int & event_id_track) ;
    std::vector< double > AngleTrackDistribution( const int & window , const unsigned int & event_id_track) ;

  private :
  // Object Truth information
  std::vector< double > _event_TLenght ;
  std::vector< int > _nu_daughters, _TPDG_Code_Primary, _Tnu_mu, _Tnu_pi, _Tnu_p, _Tnu_e, _Tnu_n, _Tnu_photon, _Tnu_others ; // truth information pdg hiearchy

  // Object RECO information
  int _hits ;
  std::vector< double > _event_RLenght ;
  Hit_level _vertex_position, _end_position ;
  std::vector< float > _reco_dQdx ;
  Track _particle_track ;
  // If loading all file
  std::vector< int > _event_hits ;
  std::vector< std::vector<double> > _event_vertex, _event_end, _event_chi2_mu, _event_chi2_pi, _event_chi2_p, _event_PIDA  ;
  std::vector< std::vector< float > > _event_reco_dQdx ;
  std::vector< Track > _event_tracks;  // maps event and track

  /**
  * Functions to guess about the geometry of the track
  */
  Track Straight( const unsigned int & event_id_track ) ;

  /**
  * Functions to guess about the geometry of the track
  */
  std::vector< double > MeanData( const int & window, const std::vector< double > & data ) ;
  std::vector< double > DevData( const int & window, const std::vector< double > & dat ) ;
  std::vector< double > CovData( const int & window, const std::vector< double > & Data_1, const std::vector< double > & Data_2 ) ;
  std::vector< double > LinearityData( const int & window, const std::vector< double > & Data_1, const std::vector< double > & Data_2 ) ;

  }; // Event


#endif
