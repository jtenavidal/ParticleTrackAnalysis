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
    TrackFitter( const std::string & track_file_path ) ; // reads all events from track_file_path

    /**
    * Functions to get properties
    */
    // For constructor 1
    unsigned int GetHits( const unsigned int & event_id_track ) ;
    void TruthParticles( unsigned int & event_id_track ) ;
    Hit_level AccessVertex( const unsigned int & event_id_track ) ;
    Hit_level AccessEnd( const unsigned int & event_id_track ) ;
    std::vector< double > GetdEdx( const unsigned int & event_id_track ) ;
    Track GetTrack( const unsigned int & event_id_track ) ;

    /**
    * Functions to check and save information
    */
    void SaveTrack( const std::string & path , const unsigned int & event_id_track ) ;
    void PrintdEdx( const std::string & path , const unsigned int & event_id_track ) const ;
    void SaveTrack( const std::string & path ) const ;
    //void PrintHipotesis( const std::string & path , const unsigned int & event_id_track ) const ) ;
    void PlotLinearityTrack( const std::string & path , const unsigned int & event_id_track ) ;
    void PlotLinearityData( const std::string & path , const unsigned int & event_id_track );//, const std::vector< double > & Data_1, const std::vector< double > & Data_2  ) ;
    void SaveStatisticsTrueEvent( const std::string & path );
    /**
    * Functions to guess about the geometry of the track
    */
    //double FitToLine( ) ;
    std::vector< TVector3 > MeanDirectionData( const unsigned int & event_id_track) ;
    std::vector< double > AngleTrackDistribution( const unsigned int & event_id_track) ;
    void StatisticsKinks( ) ;

  private :
  // Event tree information
  std::vector< bool >  _is_reconstructed, _has_reco_tracks, _has_reco_showers ;
  std::vector< int > _has_reco_daughters ;
  // Object Truth information
  std::vector< bool > _TPrimary_vcontained, _TPrimary_econtained ;
  std::vector< double > _event_TLenght, _event_TPrimaryE, _event_TPrimaryMass ;
  std::vector< int > _Tnu_daughters, _TPDG_Code_Primary, _Tnu_mu, _Tnu_pi, _Tnu_p, _Tnu_e, _Tnu_n, _Tnu_photon, _Tnu_others ; // truth information pdg hiearchy
  std::vector< std::vector<double> > _event_MC_vertex, _event_MC_end ;

  // Object RECO information
  int _hits ;
  std::vector< double > _event_RLenght ;
  Hit_level _vertex_position, _end_position ;
  std::vector< float > _reco_dEdx ;
  Track _particle_track ;
  // If loading all file
  std::vector< int > _event_hits, _rnu_daughters ;
  std::vector< std::vector<double> > _event_vertex, _event_end, _event_chi2_mu, _event_chi2_pi, _event_chi2_p, _event_PIDA  ;
  std::vector< std::vector< float > > _event_reco_dEdx ;
  std::vector< Track > _event_tracks;  // maps event and track
  // Keeping hiearchy information per event_tree
  std::vector< std::vector< int > > _event_pfps_hits, _event_pfps_type ;
  std::vector< std::vector< int > > _event_pfps_vcontained, _event_pfps_econtained ;
  std::vector< std::vector< float > > _event_pfps_length ;
  std::vector< std::vector< double > > _event_pfps_dir_start_x, _event_pfps_dir_start_y, _event_pfps_dir_start_z, _event_pfps_dir_end_x,  _event_pfps_dir_end_y, _event_pfps_dir_end_z ;
  std::vector< std::vector< double > > _event_pfps_start_x, _event_pfps_start_y, _event_pfps_start_z, _event_pfps_end_x,  _event_pfps_end_y, _event_pfps_end_z ;


  double FindMaxCoordinate( const unsigned int & event_id_track , const int & coordinate_id ) ;
  double FindMinCoordinate( const unsigned int & event_id_track , const int & coordinate_id ) ;

  /**
  * Functions to guess about the geometry of the track
  */
  std::vector< double > MeanData( const std::vector< double > & data ) ;
  std::vector< double > DevData( const std::vector< double > & dat ) ;
  std::vector< double > CovData( const std::vector< double > & Data_1, const std::vector< double > & Data_2 ) ;
  std::vector< double > LinearityData( const std::vector< double > & Data_1, const std::vector< double > & Data_2 ) ;
  std::vector< std::vector< double > > FindMinimumLinearityPosition( const unsigned int & event_id_track );

  }; // Event


#endif
