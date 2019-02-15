#include <iostream>
#include "../include/TrackFitter.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

void MainTest(){
/**
  * DEFINING OBJECTS
  */
TrackFitter reco_track_all_pi( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_new_pi.root") ;
TrackFitter reco_track_all_mu( "/home/jtenavidal/TrackFitter/ParticleTrackAnalysis/output_trees_trackID_new_mu.root") ;

/**
  * Containment study
  * 1 - Is the whole event contained ?
  * 2 - Is the primary contained ?
  * 3 - Contained / scaping vs event lenght
  */

TH1F * h_contained_pi = new TH1F( "h_contained", "Contained pion event ", 2, -0.5, 1.5 );
TH1F * h_contained_mu = new TH1F( "h_contained", "Contained muon event ", 2, -0.5, 1.5 );

for( unsigned int event = 1 ; event < 100 ; ++ event ){
  if ( reco_track_all_pi.GetHits( event ) < 35 ) continue ;
  if( reco_track_all_pi.IsContained(event) == true ) { h_contained_pi -> Fill ( 1 ) ;
  } else h_contained_pi ->Fill ( 0 ) ;
}
for( unsigned int event = 1 ; event < 100 ; ++ event ){
  if ( reco_track_all_mu.GetHits( event ) < 35 ) continue ;
  if( reco_track_all_mu.IsContained(event) == true ) { h_contained_mu -> Fill ( 1 ) ;
  } else h_contained_mu ->Fill ( 0 ) ;
}

TCanvas *c = new TCanvas() ;
TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);

h_contained_pi->SetLineColor(kRed);
h_contained_mu->SetLineColor(kBlue);
h_contained_pi->SetLineStyle(2);
h_contained_mu->SetLineStyle(1);
h_contained_pi->SetTitle("");
h_contained_mu->SetTitle("");
h_contained_pi->GetXaxis()->SetTitle( "containment" ) ;
h_contained_pi->GetYaxis()->SetTitle( "Events" ) ;
h_contained_mu->GetXaxis()->SetTitle( "containment" ) ;
h_contained_mu->GetYaxis()->SetTitle( "Events" ) ;
h_contained_pi -> Draw("hist") ;
h_contained_mu -> Draw("hist same") ;
h_contained_pi->SetStats(0) ;
leg->AddEntry(h_contained_mu, " Contained #mu^{-} event") ;
leg->AddEntry(h_contained_pi, " Contained #pi^{-} event") ;
leg->Draw();
c->SaveAs( "event_contaiment.root" ) ;


TH1F * h_contained_primary_pi = new TH1F( "h_contained_primary_pi", "Contained primary pions in event ", 2, -0.5, 1.5 );
TH1F * h_contained_primary_mu = new TH1F( "h_contained_primary_mu", "Contained primary muons in event ", 2, -0.5, 1.5 );
std::map<int, std::vector<int> >::iterator it ;

for( unsigned int event = 1 ; event < 100 ; ++ event ){
  if ( reco_track_all_pi.GetHits( event ) < 35 ) continue ;
  for( unsigned int i = 0 ; i < reco_track_all_pi.GetDaughters(event) + 1 ; ++i ){
    it = (reco_track_all_pi.FindHiearchy( event )).begin() ;
    if ( i != it->first ) continue ; // just looking for primary track
    if ( reco_track_all_pi.TruePrimaryVertexContained(event) == false ) continue ; // vertex must be always contained
    if( reco_track_all_pi.IsContained(event) == true ) { h_contained_primary_pi -> Fill ( 1 ) ;
    } else h_contained_primary_pi ->Fill ( 0 ) ; // scapes!
  }
}

for( unsigned int event = 1 ; event < 100 ; ++ event ){
  if ( reco_track_all_mu.GetHits( event ) < 35 ) continue ;
  for( unsigned int i = 0 ; i < reco_track_all_mu.GetDaughters(event) + 1 ; ++i ){
    it = (reco_track_all_mu.FindHiearchy( event )).begin() ;
    if ( i != it->first ) continue ; // just looking for primary track
    if ( reco_track_all_mu.TruePrimaryVertexContained(event) == false ) continue ;
    if( reco_track_all_mu.IsContained(event) == true ) { h_contained_primary_mu -> Fill ( 1 ) ;
    } else h_contained_primary_mu ->Fill ( 0 ) ;
  }
}

TCanvas *c2 = new TCanvas() ;
TLegend *leg2 = new TLegend(0.1,0.7,0.48,0.9);

h_contained_primary_pi->SetLineColor(kRed);
h_contained_primary_mu->SetLineColor(kBlue);
h_contained_primary_pi->SetLineStyle(2);
h_contained_primary_mu->SetLineStyle(1);
h_contained_primary_pi->SetTitle("");
h_contained_primary_mu->SetTitle("");
h_contained_primary_pi->GetXaxis()->SetTitle( "containment" ) ;
h_contained_primary_pi->GetYaxis()->SetTitle( "Events" ) ;
h_contained_primary_mu->GetXaxis()->SetTitle( "containment" ) ;
h_contained_primary_mu->GetYaxis()->SetTitle( "Events" ) ;
h_contained_primary_pi -> Draw("hist") ;
h_contained_primary_mu -> Draw("hist same") ;
h_contained_primary_pi->SetStats(0) ;
leg2->AddEntry(h_contained_primary_mu, " Contained #mu^{-} event") ;
leg2->AddEntry(h_contained_primary_pi, " Contained #pi^{-} event") ;
leg2->Draw();
c2->SaveAs( "primary_contaiment.root" ) ;



TH1F * h_Lcontained_pi = new TH1F( "h_Lcontained_pi", "Contained pion event ", 50, -0.5, 300 );
TH1F * h_Lcontained_mu = new TH1F( "h_Lcontained_pi", "Contained muon event ", 50, -0.5, 500 );
TH1F * h_Lescape_pi = new TH1F( "h_Lescape_pi", "Escaping pion event ", 30, -0.5, 300 );
TH1F * h_Lescape_mu = new TH1F( "h_Lescape_pi", "Escaping muon event ", 30, -0.5, 500 );

for( unsigned int event = 1 ; event < 100 ; ++ event ){
  if ( reco_track_all_pi.GetHits( event ) < 35 ) continue ;
  for( unsigned int i = 0 ; i < reco_track_all_pi.GetDaughters(event) + 1 ; ++i ){
    it = (reco_track_all_pi.FindHiearchy( event )).begin() ;
    if ( i != it->first ) continue ; // just looking for primary track
    if ( reco_track_all_pi.TruePrimaryVertexContained(event) == false ) continue ; // vertex must be always contained
    if( reco_track_all_pi.IsContained(event) == true ) { h_Lcontained_pi -> Fill ( reco_track_all_pi.GetTotalRecoLength( event ) ) ;
    } else h_Lescape_pi -> Fill ( reco_track_all_pi.GetTotalRecoLength( event ) ) ; // escapes!
  }
}


for( unsigned int event = 1 ; event < 100 ; ++ event ){
  if ( reco_track_all_mu.GetHits( event ) < 35 ) continue ;
  for( unsigned int i = 0 ; i < reco_track_all_mu.GetDaughters(event) + 1 ; ++i ){
    it = (reco_track_all_mu.FindHiearchy( event )).begin() ;
    if ( i != it->first ) continue ; // just looking for primary track
    if ( reco_track_all_mu.TruePrimaryVertexContained(event) == false ) continue ;
    if( reco_track_all_mu.IsContained(event) == true ) { h_Lcontained_mu -> Fill ( reco_track_all_mu.GetTotalRecoLength( event ) ) ;
    } else h_Lescape_mu -> Fill ( reco_track_all_mu.GetTotalRecoLength( event ) ) ;
  }
}


TCanvas *c3 = new TCanvas() ;
TLegend *leg3 = new TLegend(0.9,0.7,0.48,0.9);

h_Lcontained_pi->SetLineColor(kBlue);
h_Lcontained_mu->SetLineColor(kBlue);
h_Lescape_pi->SetLineColor(kRed);
h_Lescape_mu->SetLineColor(kRed);
h_Lcontained_pi->SetLineStyle(1);
h_Lcontained_mu->SetLineStyle(1);
h_Lescape_pi->SetLineStyle(2);
h_Lescape_mu->SetLineStyle(2);
h_Lcontained_pi->SetTitle("");
h_Lcontained_mu->SetTitle("");
h_Lescape_pi->SetTitle("");
h_Lescape_mu->SetTitle("");
h_Lcontained_pi->GetXaxis()->SetTitle( "Length [cm]" ) ;
h_Lcontained_pi->GetYaxis()->SetTitle( "Events" ) ;
h_Lcontained_mu->GetXaxis()->SetTitle( "Length [cm]" ) ;
h_Lcontained_mu->GetYaxis()->SetTitle( "Events" ) ;
h_Lcontained_mu -> Draw("hist") ;
h_Lescape_mu -> Draw("hist same") ;
h_Lcontained_pi->SetStats(0) ;
h_Lcontained_mu->SetStats(0) ;
leg3->AddEntry(h_Lcontained_mu, " Contained #mu^{-} event") ;
leg3->AddEntry(h_Lescape_mu, " Escaping #mu^{-} event") ;
leg3->Draw() ;

TCanvas *c4 = new TCanvas() ;
TLegend *leg4 = new TLegend(0.9,0.7,0.48,0.9);
h_Lcontained_pi -> Draw("hist") ;
h_Lescape_pi -> Draw("hist same") ;
leg4->AddEntry(h_Lcontained_pi, " Contained #pi^{-} event") ;
leg4->AddEntry(h_Lescape_pi, " Escaping #pi^{-} event") ;
leg4->Draw();
c4->SaveAs( "event_contained_length.root" ) ;



/**
  * Reconstructed particles
  * 1 -> Reconstructed primary
  * 2 -> Number of daughters
  * 3 -> Truth info about reco particles
  * 4 -> Lenght primary particle
  * 5 -> Lenght secondary particle
  * 6 -> Primary vs secondary particle lenght study
  * 7 -> Angle between primary vs secondary
  * 8 -> Miss reco e- 
  */






// TESTING FUNCTIONS HERE
// unsigned int event = 45 ; // Breaks for pions: 10, 15

// // it = (reco_track_all.FindHiearchy( event )).begin() ;
// std::cout<< "primary particle id -> " << it->first << " has #daughters = " << (it->second).size() << std::endl;
// reco_track_all.TruthParticles(event);


// reco_track_all.SaveTrack( ("Track_testing_new_constructor_") , event ) ;
/*for( int i = 1 ; i < 100 ; ++i ) {
   if( reco_track_all.GetHits( i ) > 35 && reco_track_all.IsContained( i ) == 1 ) {
     std::cout << i << std::endl;
     it = (reco_track_all.FindHiearchy( i )).begin() ;
     std::cout<< "primary particle id -> " << it->first << " has #daughters = " << (it->second).size() << std::endl;
     // reco_track_all.SaveTrack( ("Track_testing_new_constructor_"+std::to_string(i)) , i ) ;
   }
}*/
/*  it = (reco_track_all_pi.FindHiearchy( event )).begin() ;
  if ( i != it->first ) continue ; // just looking for primary track
*/
} // MainTest()1
