## Track Analyzer: Developed to look for specific features characteristic for muon and pion events -> develop a new PID method


OUTPUT: RooT file containging three trees,

Event Tree contains general information of the event regardless of its nature
   [X] event id
   [ ] time now
   [ ] mother particle pdg code (particle beam , sanity check )
   [ ] starting point
   [ ] End point 
   [ ] number of particles (including daughter particles) per type 
   [ ] number of tracks
   [ ] number of showers 

MC Particle tree contains information at the truth level
    Must have:
    [X] ID
    [X] PDG code
    [X] Mass
    [X] Momentum 
    [ ] direction -> Can have a TLorenzVector
    [X] Energy 
    [ ] Process
    [ ] Keep hiearchy information : possible at this level?
       [X] NumberDaughters() {int type}
       [ ] Daughter( ) {int type}
    [?] True hit level information -> detector effects, smearing (future?)
       [ ] simb::MCTrajectory Trajectory 


Reco track tree contains information of the reconstructed track, including hit level information
    [ ] Chi2 hipotesis for: muons, pions, kaons and protons (the last two should not be seen but just in case)
    [ ] dE/dx (must be able to split it if needed)
    [ ] Residual range
    [ ] Lenght 
    [ ] Kinetic energy and missing energy 
    [ ] Hit information: number of hits, hit_tpc, hit_plane
    [ ] Precise hit information: 3D track reconstructed ( position, time and calorimetry )
    [ ] Keep hiearchy information
  