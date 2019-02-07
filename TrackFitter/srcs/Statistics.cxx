// #include "TrackFitter.h"
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

/**
  * STATISTICS : save information about event
  * 1 - Primary: vertex contained, end contained, scapes vs lenght
  * 2 - Length primary vs daughter
  * 3 - Reconstructed events vs empty ones
  * 4 - Showers/track events 
  */
