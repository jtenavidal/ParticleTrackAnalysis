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

 TrackFitter::TrackFitter( const std::string & track_file_path ){
   TFile track_file( track_file_path.c_str() );
   // Event Tree
   TTree * event_tree   = (TTree*) track_file.Get("event_tree") ;
   TBranch * event_id = event_tree->GetBranch("event_id");
   TBranch * is_reconstructed = event_tree->GetBranch("is_reconstructed");
   TBranch * has_reco_daughters = event_tree->GetBranch("has_reco_daughters");
   TBranch * has_reco_tracks = event_tree->GetBranch("has_reco_tracks");
   TBranch * has_reco_showers = event_tree->GetBranch("has_reco_showers");
   // MC Tree
   TTree * mcparticle_tree   = (TTree*) track_file.Get("mcparticle_tree") ;
   TBranch * mc_event_id = mcparticle_tree->GetBranch("event_id");
   TBranch * Track_ID = mcparticle_tree->GetBranch("fTrack_id");
   TBranch * trueEnergy = mcparticle_tree->GetBranch("ftrueEnergy");
   TBranch * PDG_Code = mcparticle_tree->GetBranch("fPDG_Code");
   TBranch * Mass = mcparticle_tree->GetBranch("fMass");
   TBranch * Px = mcparticle_tree->GetBranch("fpx");
   TBranch * Py = mcparticle_tree->GetBranch("fpy");
   TBranch * Pz = mcparticle_tree->GetBranch("fpz");
   TBranch * Pt = mcparticle_tree->GetBranch("fpt");
   TBranch * P = mcparticle_tree->GetBranch("fp");
   TBranch * Num_Daughters = mcparticle_tree->GetBranch("fNumDaughters");
   TBranch * Daughter_mu = mcparticle_tree->GetBranch("fDaughter_mu");
   TBranch * Daughter_pi = mcparticle_tree->GetBranch("fDaughter_pi");
   TBranch * Daughter_e = mcparticle_tree->GetBranch("fDaughter_e");
   TBranch * Daughter_p = mcparticle_tree->GetBranch("fDaughter_p");
   TBranch * Daughter_n = mcparticle_tree->GetBranch("fDaughter_n");
   TBranch * Daughter_photon = mcparticle_tree->GetBranch("fDaughter_photon");
   TBranch * Daughter_other = mcparticle_tree->GetBranch("fDaughter_other");
   TBranch * MC_Lenght= mcparticle_tree->GetBranch("fMC_Lenght");
   // RECO Tree
   TTree * recoparticle_tree   = (TTree*) track_file.Get("recoparticle_tree") ;
   TBranch * reco_event_id = recoparticle_tree->GetBranch("event_id");
   TBranch * primary_vcontained = recoparticle_tree->GetBranch("primary_vcontained");
   TBranch * primary_econtained = recoparticle_tree->GetBranch("primary_econtained");
   TBranch * r_pdg_primary = recoparticle_tree->GetBranch("r_pdg_primary");
   TBranch * nu_daughters = recoparticle_tree->GetBranch("r_nu_daughters");
   TBranch * Length = recoparticle_tree->GetBranch("rLength");
   TBranch * rnu_hits = recoparticle_tree->GetBranch("rnu_hits");
   TBranch * r_chi2_mu = recoparticle_tree->GetBranch("r_chi2_mu");
   TBranch * r_chi2_pi = recoparticle_tree->GetBranch("r_chi2_pi");
   TBranch * r_chi2_p = recoparticle_tree->GetBranch("r_chi2_p");
   TBranch * r_PIDA = recoparticle_tree->GetBranch("r_PIDA");
   TBranch * r_missing_energy = recoparticle_tree->GetBranch("r_missing_energy");
   TBranch * r_KineticEnergy = recoparticle_tree->GetBranch("r_KineticEnergy");
   TBranch * r_track_x = recoparticle_tree->GetBranch("r_track_x");
   TBranch * r_track_y = recoparticle_tree->GetBranch("r_track_y");
   TBranch * r_track_z = recoparticle_tree->GetBranch("r_track_z");
   TBranch * r_track_dEdx = recoparticle_tree->GetBranch("r_track_dEdx");
   TBranch * r_dEdx = recoparticle_tree->GetBranch("r_dEdx");
   TBranch * pfps_hits = recoparticle_tree->GetBranch("pfps_hits");
   TBranch * pfps_type = recoparticle_tree->GetBranch("pfps_type");
   TBranch * event_vcontained = recoparticle_tree->GetBranch("event_vcontained");
   TBranch * event_econtained = recoparticle_tree->GetBranch("event_econtained");
   TBranch * pfps_length = recoparticle_tree->GetBranch("pfps_length");
   TBranch * pfps_dir_start_x = recoparticle_tree->GetBranch("pfps_dir_start_x");
   TBranch * pfps_dir_start_y = recoparticle_tree->GetBranch("pfps_dir_start_y");
   TBranch * pfps_dir_start_z = recoparticle_tree->GetBranch("pfps_dir_start_z");
   TBranch * pfps_dir_end_x = recoparticle_tree->GetBranch("pfps_dir_end_x");
   TBranch * pfps_dir_end_y = recoparticle_tree->GetBranch("pfps_dir_end_y");
   TBranch * pfps_dir_end_z = recoparticle_tree->GetBranch("pfps_dir_end_z");
   TBranch * pfps_start_x = recoparticle_tree->GetBranch("pfps_start_x");
   TBranch * pfps_start_y = recoparticle_tree->GetBranch("pfps_start_y");
   TBranch * pfps_start_z = recoparticle_tree->GetBranch("pfps_start_z");
   TBranch * pfps_end_x = recoparticle_tree->GetBranch("pfps_end_x");
   TBranch * pfps_end_y = recoparticle_tree->GetBranch("pfps_end_y");
   TBranch * pfps_end_z = recoparticle_tree->GetBranch("pfps_end_z");


   /**
     * General event information
     */
     for( unsigned int i = 0; i < event_tree->GetEntries(); ++i ){
       event_tree->GetEntry(i);
       _is_reconstructed.push_back(is_reconstructed->GetLeaf("is_reconstructed")->GetValue());
       _has_reco_daughters.push_back(has_reco_daughters->GetLeaf("has_reco_daughters")->GetValue());
       _has_reco_tracks.push_back(has_reco_tracks->GetLeaf("has_reco_tracks")->GetValue());
       _has_reco_showers.push_back(has_reco_showers->GetLeaf("has_reco_showers")->GetValue());
      }

   /**
     * MC information
     */
    for( unsigned int i = 0; i < mcparticle_tree->GetEntries(); ++i ){
      mcparticle_tree->GetEntry(i);
      //_event_TLenght.push_back( MC_Lenght->GetLeaf("fMCLength")->GetValue() ) ; // still breaks
      _Tnu_daughters.push_back(Num_Daughters->GetLeaf("fNumDaughters")->GetValue());
      _TPDG_Code_Primary.push_back( PDG_Code->GetLeaf("fPDG_Code")->GetValue() ) ;
      _Tnu_mu.push_back( Daughter_mu->GetLeaf("fDaughter_mu")->GetValue() ) ;
      _Tnu_pi.push_back( Daughter_pi->GetLeaf("fDaughter_pi")->GetValue() ) ;
      _Tnu_p.push_back( Daughter_p->GetLeaf("fDaughter_p")->GetValue() ) ;
      _Tnu_e.push_back( Daughter_e->GetLeaf("fDaughter_e")->GetValue() ) ;
      _Tnu_n.push_back( Daughter_n->GetLeaf("fDaughter_n")->GetValue() );
      _Tnu_photon.push_back( Daughter_photon->GetLeaf("fDaughter_photon")->GetValue() );
      _Tnu_others.push_back( Daughter_other->GetLeaf("fDaughter_other")->GetValue() );
     }

/**
  * Reco information
  */
   Hit_level track_hit_x, track_hit_y, track_hit_z, track_hit_dEdx;
   std::vector< double > chi2_mu, chi2_pi, chi2_p, pida ;
   std::vector< int > vcontained, econtained ;
   std::vector< int > pfps_hits_i, pfps_type_i ;
   std::vector< float > pfps_length_i ;
   std::vector< double > pfps_dir_start_i_x, pfps_dir_start_i_y, pfps_dir_start_i_z, pfps_dir_end_i_x, pfps_dir_end_i_y, pfps_dir_end_i_z ;
   std::vector< double > pfps_start_i_x, pfps_start_i_y, pfps_start_i_z, pfps_end_i_x, pfps_end_i_y, pfps_end_i_z ;


   for( unsigned int i = 0; i < recoparticle_tree->GetEntries(); ++i ){
     _particle_track.clear();
     track_hit_x.clear();
     track_hit_y.clear();
     track_hit_z.clear();
     track_hit_dEdx.clear();
     _reco_dEdx.clear();
     chi2_mu.clear();
     chi2_pi.clear();
     chi2_p.clear();
     pida.clear();
     vcontained.clear() ;
     econtained.clear() ;
     pfps_hits_i.clear() ;
     pfps_type_i.clear() ;
     pfps_length_i.clear() ;
     pfps_dir_start_i_x.clear() ;
     pfps_dir_start_i_y.clear() ;
     pfps_dir_start_i_z.clear() ;
     pfps_dir_end_i_x.clear() ;
     pfps_dir_end_i_y.clear() ;
     pfps_dir_end_i_z.clear() ;
     pfps_start_i_x.clear() ;
     pfps_start_i_y.clear() ;
     pfps_start_i_z.clear() ;
     pfps_end_i_x.clear() ;
     pfps_end_i_y.clear() ;
     pfps_end_i_z.clear() ;


     recoparticle_tree->GetEntry(i);
     _hits = rnu_hits->GetLeaf("rnu_hits")->GetValue() ;
     _event_hits.push_back( _hits ) ;
     _event_primamry_vcontained.push_back( primary_vcontained->GetLeaf("primary_vcontained")->GetValue() ) ;
     _event_primary_econtained.push_back( primary_econtained->GetLeaf("primary_econtained")->GetValue() ) ;

     if ( _hits != 0 ) {
       _event_RLenght.push_back( Length->GetLeaf("rLength")->GetValue() );
       for( int j = 0; j < _hits; ++j ) {
         track_hit_x.push_back( r_track_x->GetLeaf("r_track_x")->GetValue(j));
         track_hit_y.push_back( r_track_y->GetLeaf("r_track_y")->GetValue(j));
         track_hit_z.push_back( r_track_z->GetLeaf("r_track_z")->GetValue(j));
         track_hit_dEdx.push_back( r_track_dEdx->GetLeaf("r_track_dEdx")->GetValue(j));
       }
       _particle_track.push_back(track_hit_x);
       _particle_track.push_back(track_hit_y);
       _particle_track.push_back(track_hit_z);
       _particle_track.push_back(track_hit_dEdx);

       _event_tracks.push_back( _particle_track ) ;

       for( int j = 0; j < _hits; ++j ) _reco_dEdx.push_back( recoparticle_tree->GetLeaf("r_dEdx")->GetValue(j));
       _event_reco_dEdx.push_back( _reco_dEdx ) ;
       _rnu_daughters.push_back( nu_daughters->GetLeaf("r_nu_daughters")->GetValue() ) ;
       for( int j = 0; j < nu_daughters->GetLeaf("r_nu_daughters")->GetValue() + 1 ; ++j){ // loop over all particles from track
         chi2_mu.push_back( r_chi2_mu->GetLeaf("r_chi2_mu")->GetValue(j));
         chi2_pi.push_back( r_chi2_pi->GetLeaf("r_chi2_pi")->GetValue(j));
         chi2_p.push_back ( r_chi2_p ->GetLeaf("r_chi2_p")->GetValue(j));
         pida.push_back   ( r_PIDA   ->GetLeaf("r_PIDA")->GetValue(j));
         vcontained.push_back( event_vcontained->GetLeaf("event_vcontained")->GetValue(j) ) ;
         econtained.push_back( event_econtained->GetLeaf("event_econtained")->GetValue(j) ) ;
         pfps_hits_i.push_back( pfps_hits->GetLeaf("pfps_hits")->GetValue(j) ) ;
         pfps_type_i.push_back( pfps_type->GetLeaf("pfps_type")->GetValue(j) ) ;
         pfps_length_i.push_back( pfps_length->GetLeaf("pfps_length")->GetValue(j) ) ;
         pfps_dir_start_i_x.push_back( pfps_dir_start_x->GetLeaf("pfps_dir_start_x")->GetValue(j) ) ;
         pfps_dir_start_i_y.push_back( pfps_dir_start_y->GetLeaf("pfps_dir_start_y")->GetValue(j) ) ;
         pfps_dir_start_i_z.push_back( pfps_dir_start_z->GetLeaf("pfps_dir_start_z")->GetValue(j) ) ;
         pfps_dir_end_i_x.push_back( pfps_dir_end_x->GetLeaf("pfps_dir_end_x")->GetValue(j) ) ;
         pfps_dir_end_i_y.push_back( pfps_dir_end_y->GetLeaf("pfps_dir_end_y")->GetValue(j) ) ;
         pfps_dir_end_i_z.push_back( pfps_dir_end_z->GetLeaf("pfps_dir_end_z")->GetValue(j) ) ;
         pfps_start_i_x.push_back( pfps_start_x->GetLeaf("pfps_start_x")->GetValue(j) ) ;
         pfps_start_i_y.push_back( pfps_start_y->GetLeaf("pfps_start_y")->GetValue(j) ) ;
         pfps_start_i_z.push_back( pfps_start_z->GetLeaf("pfps_start_z")->GetValue(j) ) ;
         pfps_end_i_x.push_back( pfps_end_x->GetLeaf("pfps_end_x")->GetValue(j) ) ;
         pfps_end_i_y.push_back( pfps_end_y->GetLeaf("pfps_end_y")->GetValue(j) ) ;
         pfps_end_i_z.push_back( pfps_end_z->GetLeaf("pfps_end_z")->GetValue(j) ) ;
       }
       _event_chi2_mu.push_back(chi2_mu);
       _event_chi2_pi.push_back(chi2_pi);
       _event_chi2_p.push_back(chi2_p);
       _event_PIDA.push_back(pida);
       _event_pfps_hits.push_back(pfps_hits_i);
       _event_pfps_type.push_back(pfps_type_i);
       _event_pfps_vcontained.push_back(vcontained);
       _event_pfps_econtained.push_back(econtained);
       _event_pfps_length.push_back(pfps_length_i);
       _event_pfps_dir_start_x.push_back(pfps_dir_start_i_x);
       _event_pfps_dir_start_y.push_back(pfps_dir_start_i_y);
       _event_pfps_dir_start_z.push_back(pfps_dir_start_i_z);
       _event_pfps_dir_end_x.push_back(pfps_dir_end_i_x);
       _event_pfps_dir_end_y.push_back(pfps_dir_end_i_y);
       _event_pfps_dir_end_z.push_back(pfps_dir_end_i_z);
       _event_pfps_start_x.push_back(pfps_start_i_x);
       _event_pfps_start_y.push_back(pfps_start_i_y);
       _event_pfps_start_z.push_back(pfps_start_i_z);
       _event_pfps_end_x.push_back(pfps_end_i_x);
       _event_pfps_end_y.push_back(pfps_end_i_y);
       _event_pfps_end_z.push_back(pfps_end_i_z);

   } else {
     // seting to zero for empty events
     track_hit_x.push_back ( 0 ) ;
     track_hit_y.push_back ( 0 ) ;
     track_hit_z.push_back ( 0 ) ;
     track_hit_dEdx.push_back( 0 ) ;
     _particle_track.push_back(track_hit_x) ;
     _particle_track.push_back(track_hit_y) ;
     _particle_track.push_back(track_hit_z) ;
     _particle_track.push_back(track_hit_dEdx) ;
     _reco_dEdx.push_back( 0 ) ;
     _TPDG_Code_Primary.push_back( 0 ) ;
     _Tnu_mu.push_back( 0 ) ;
     _Tnu_pi.push_back( 0 ) ;
     _Tnu_p.push_back( 0 ) ;
     _Tnu_e.push_back( 0 ) ;
     _Tnu_n.push_back( 0 );
     _Tnu_photon.push_back( 0 );
     _Tnu_others.push_back( 0 );
     _event_tracks.push_back( _particle_track ) ;
     _event_reco_dEdx.push_back( _reco_dEdx ) ;

   }
 }

 _TPDG_Code_Primary.clear();
 _Tnu_mu.clear();
 _Tnu_pi.clear();
 _Tnu_p.clear();
 _Tnu_e.clear();
 _Tnu_n.clear();
 _Tnu_photon.clear();
 _Tnu_others.clear();
 _Tnu_daughters.clear();

}
