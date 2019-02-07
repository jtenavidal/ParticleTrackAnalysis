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




void TrackFitter::TruthParticles( unsigned int & event_id_track ){
  if ( event_id_track == 0 ) event_id_track = 1 ;
  std::cout<< "Primary particle PDG = " << _TPDG_Code_Primary[ event_id_track-1 ] << std::endl;
  std::cout<< "#mu = " << _Tnu_mu[event_id_track-1] << " #pi = " << _Tnu_pi[event_id_track-1] << " #p = " << _Tnu_p[event_id_track-1] << " #e = " << _Tnu_e[event_id_track-1]<< " #photon = " << _Tnu_photon[event_id_track-1] << " #others = " << _Tnu_others[event_id_track-1]<< std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// Truth information functions :                                             //
///////////////////////////////////////////////////////////////////////////////
void TrackFitter::SaveStatisticsTrueEvent( const std::string & path ) {
/**
  * -> Badruns muons and pions
  * -> Muons with Michel e
  * -> Muons no Michel e
  * -> Breakdown muons and other particles
  * -> particle and number of kinks
  * -> Lenght primary particle
  * -> Ratio lenght primary / 1-daughter
  */
  TCanvas *c = new TCanvas();
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);
  TLegend *leg = new TLegend(0.9,0.7,0.48,0.9);
  THStack * h_mu_stats = new THStack( "h_mu_stats", "Muon Statistics" ) ;
  TH1D *h_mu_rdaughter = new TH1D("h_mu_rdaughter", " Muon Statistics: reconstructed daughters ", 3, 0, 2 );
  TH1D *h_mu_charged = new TH1D("h_mu_charged", " Muon Statistics: Number of total charged daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_e = new TH1D("h_mu_charged_e", " Muon Statistics: Number of total e daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_mu = new TH1D("h_mu_charged_mu", " Muon Statistics: Number of total mu daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_pi = new TH1D("h_mu_charged_pi", " Muon Statistics: Number of total pi daughters ", 3, 0, 2 );
  TH1D *h_mu_charged_p = new TH1D("h_mu_charged_p", " Muon Statistics: Number of total p daughters ", 3, 0, 2 );

  for( unsigned int j = 0 ; j < _rnu_daughters.size() ; ++j ){
    if( _TPDG_Code_Primary[j] != 13 ) continue ;
    if( _rnu_daughters[j] == 0 ) h_mu_rdaughter->Fill(0) ;
    if( _rnu_daughters[j] == 1 ) h_mu_rdaughter->Fill(1) ;
    if( _rnu_daughters[j] > 1 ) h_mu_rdaughter->Fill(2) ; // needs breakdown to understand
  }

  for( unsigned int j = 0 ; j < _TPDG_Code_Primary.size() ; ++j ){
    if ( _TPDG_Code_Primary[j] != 13 ) continue ;
    if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(0) ;
    } else if ( _Tnu_pi[j] == 1  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(1) ; h_mu_charged_pi->Fill(1) ;
    } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 1  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(1) ; h_mu_charged_mu->Fill(1) ;
    } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 1  & _Tnu_p[j] == 0  ) { h_mu_charged->Fill(1) ; h_mu_charged_e->Fill(1) ;
    } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 1  ) { h_mu_charged->Fill(1) ; h_mu_charged_p->Fill(1) ;
    } else { h_mu_charged->Fill(2) ; }
  }

  h_mu_charged_mu->SetFillColor(kRed);
  h_mu_charged_pi->SetFillColor(kBlue);
  h_mu_charged_e->SetFillColor(kGreen);
  h_mu_charged_p->SetFillColor(kYellow);
  h_mu_charged_mu->SetFillStyle(3004);
  h_mu_charged_pi->SetFillStyle(3004);
  h_mu_charged_e->SetFillStyle(3004);
  h_mu_charged_p->SetFillStyle(3004);
  leg->AddEntry(h_mu_charged_mu, "True muon daughters ") ;
  leg->AddEntry(h_mu_charged_pi, "True pion daughters ") ;
  leg->AddEntry(h_mu_charged_e, "True electron daughters ") ;
  leg->AddEntry(h_mu_charged_p, "True proton daughters ") ;
  h_mu_charged_mu->Scale(1/h_mu_charged->GetEntries());
  h_mu_charged_pi->Scale(1/h_mu_charged->GetEntries());
  h_mu_charged_e->Scale(1/h_mu_charged->GetEntries());
  h_mu_charged_p->Scale(1/h_mu_charged->GetEntries());
  h_mu_stats->Add(h_mu_charged_mu) ;
  h_mu_stats->Add(h_mu_charged_pi) ;
  h_mu_stats->Add(h_mu_charged_e) ;
  h_mu_stats->Add(h_mu_charged_p) ;
  h_mu_stats->Draw("hist") ;
  /////////////////////////////////////////////////////////////////////////////
  h_mu_rdaughter->SetLineColor(kPink);
  h_mu_rdaughter->SetLineWidth(2);
  h_mu_rdaughter->SetFillStyle(3003);
  h_mu_rdaughter->SetFillColor(kPink);
  h_mu_rdaughter->GetXaxis()->SetTitle("Number reconstructed particles");
  h_mu_rdaughter->GetYaxis()->SetTitle("%");
  h_mu_rdaughter->Scale(1/h_mu_rdaughter->GetEntries());
  h_mu_rdaughter->GetYaxis()->SetRangeUser(0.,1.);
  leg->AddEntry(h_mu_rdaughter, "Reconstructed secondary tracks ") ;
  h_mu_rdaughter->Draw("hist SAME");

  leg->Draw();
  c->SaveAs((path+"_muon_daughters_study.root").c_str());
  leg->Clear();
  c->Clear();

  /*c->Clear();

  TH1D *h_mu_Tdaughter = new TH1D("h_mu_Tdaughter", " Muon Statistics: true daughters ", 3, 0, 2 );
  for( unsigned int j = 0 ; j < _TPDG_Code_Primary.size() ; ++j ){
    if( _TPDG_Code_Primary[j] != 13 ) continue ;
    if( _Tnu_daughters[j] == 0 ) h_mu_Tdaughter->Fill(0) ;
    if( _Tnu_daughters[j] == 1 ) h_mu_Tdaughter->Fill(1) ;
    if( _Tnu_daughters[j] > 1 ) h_mu_Tdaughter->Fill(2) ; // needs breakdown to understand
  }

  h_mu_Tdaughter->SetLineColor(2);
  h_mu_Tdaughter->GetXaxis()->SetTitle("Number daughter electrons");
  h_mu_Tdaughter->GetYaxis()->SetTitle("%");
//  h_mu_Tdaughter->Scale(1/h_mu_Tdaughter->GetEntries());
  h_mu_Tdaughter->Draw("hist");
  c->SaveAs((path+"_muon_rdaughter.root").c_str());
  //  c->Clear();

*/
  // Histogram length second track

// HISTOGRAM WIHT NUMBER DAUGHTERS FOR PION> ADD ON TOP BREAKDOWN . FIX TTREE
THStack * h_pi_stats = new THStack( "h_pi_stats", "Muon Statistics" ) ;
TH1D *h_pi_rdaughter = new TH1D("h_pi_rdaughter", " Muon Statistics: reconstructed daughters ", 3, 0, 2 );
TH1D *h_pi_charged = new TH1D("h_pi_charged", " Muon Statistics: Number of total charged daughters ", 3, 0, 2 );
TH1D *h_pi_charged_e = new TH1D("h_pi_charged_e", " Muon Statistics: Number of total e daughters ", 3, 0, 2 );
TH1D *h_pi_charged_mu = new TH1D("h_pi_charged_mu", " Muon Statistics: Number of total mu daughters ", 3, 0, 2 );
TH1D *h_pi_charged_pi = new TH1D("h_pi_charged_pi", " Muon Statistics: Number of total pi daughters ", 3, 0, 2 );
TH1D *h_pi_charged_p = new TH1D("h_pi_charged_p", " Muon Statistics: Number of total p daughters ", 3, 0, 2 );

for( unsigned int j = 0 ; j < _rnu_daughters.size() ; ++j ){
  if( _TPDG_Code_Primary[j] != 211 ) continue ;
  if( _rnu_daughters[j] == 0 ) h_pi_rdaughter->Fill(0) ;
  if( _rnu_daughters[j] == 1 ) h_pi_rdaughter->Fill(1) ;
  if( _rnu_daughters[j] > 1 ) h_pi_rdaughter->Fill(2) ; // needs breakdown to understand
}

for( unsigned int j = 0 ; j < _TPDG_Code_Primary.size() ; ++j ){
  if ( _TPDG_Code_Primary[j] != 211 ) continue ;
  if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(0) ; std::cout<< "CUU" << std::endl;
  } else if ( _Tnu_pi[j] == 1  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(1) ; h_pi_charged_pi->Fill(1) ;
  } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 1  & _Tnu_e[j] == 0  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(1) ; h_pi_charged_mu->Fill(1) ;
  } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 1  & _Tnu_p[j] == 0  ) { h_pi_charged->Fill(1) ; h_pi_charged_e->Fill(1) ;
  } else if ( _Tnu_pi[j] == 0  & _Tnu_mu[j] == 0  & _Tnu_e[j] == 0  & _Tnu_p[j] == 1  ) { h_pi_charged->Fill(1) ; h_pi_charged_p->Fill(1) ;
  } else { h_pi_charged->Fill(2) ; }
}
///////////////////////////////////////////////////////////////////////////////
h_pi_charged->SetLineColor(1);
leg->AddEntry(h_pi_stats, "Total charged particles") ;
h_pi_charged->Draw("hist SAME") ;
/////////////////////////////////////////////////////////////////////////////
h_pi_charged_mu->SetFillColor(kRed);
h_pi_charged_pi->SetFillColor(kBlue);
h_pi_charged_e->SetFillColor(kGreen);
h_pi_charged_p->SetFillColor(kYellow);
h_pi_charged_mu->SetFillStyle(3004);
h_pi_charged_pi->SetFillStyle(3004);
h_pi_charged_e->SetFillStyle(3004);
h_pi_charged_p->SetFillStyle(3004);
leg->AddEntry(h_pi_charged_mu, "True muon daughters ") ;
leg->AddEntry(h_pi_charged_pi, "True pion daughters ") ;
leg->AddEntry(h_pi_charged_e, "True electron daughters ") ;
leg->AddEntry(h_pi_charged_p, "True proton daughters ") ;
//h_pi_charged_mu->Scale(1/h_pi_charged->GetEntries());
//  h_pi_charged_pi->Scale(1/h_pi_charged->GetEntries());
//h_pi_charged_e->Scale(1/h_pi_charged->GetEntries());
//h_pi_charged_p->Scale(1/h_pi_charged->GetEntries());
h_pi_stats->Add(h_pi_charged_mu) ;
h_pi_stats->Add(h_pi_charged_pi) ;
h_pi_stats->Add(h_pi_charged_e) ;
h_pi_stats->Add(h_pi_charged_p) ;
h_pi_stats->Draw("hist SAME") ;
/////////////////////////////////////////////////////////////////////////////
h_pi_rdaughter->SetLineColor(kPink);
h_pi_rdaughter->SetLineWidth(2);
h_pi_rdaughter->SetFillStyle(3003);
h_pi_rdaughter->SetFillColor(kPink);
h_pi_rdaughter->GetXaxis()->SetTitle("Number reconstructed particles");
h_pi_rdaughter->GetYaxis()->SetTitle("%");
//h_pi_rdaughter->Scale(1/h_pi_rdaughter->GetEntries());
//h_pi_rdaughter->GetYaxis()->SetRangeUser(0.,1.);
leg->AddEntry(h_pi_rdaughter, "Reconstructed secondary tracks ") ;
h_pi_rdaughter->Draw("hist SAME");
leg->Draw();
c->SaveAs((path+"_pion_daughters_study.root").c_str());
//c->Clear();
}
