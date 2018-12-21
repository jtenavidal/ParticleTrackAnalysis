# Track Analyzer: 

Developed to look for specific features characteristic for muon and pion events -> develop a new PID method

## OUTPUT file: 

RooT file containging three trees,

## Information in each tree: 

### Event Tree contains general information of the event regardless of its nature

- [X] event id
- [ ] time now
- [ ] mother particle pdg code (particle beam , sanity check )
- [ ] starting point
- [ ] End point 
- [ ] number of particles (including daughter particles) per type 
- [ ] number of tracks
- [ ] number of showers 

### MC Particle tree contains information at the truth level
- [X] ID
- [X] PDG code
- [X] Mass
- [X] Momentum 
- [ ] direction -> Can have a TLorenzVector
- [X] Energy 
- [ ] Process

Keep hiearchy information : possible at this level?
- [X] NumberDaughters() {int type}
- [ ] Daughter( ) {int type}

True hit level information
- [X] simb::MCTrajectory Trajectory -> Positioning and energy available
- [X] Enable plotting  ( 3D HIST with E as weight ) : SHOULD SET AS A FUNCTION
- [X] MCLenght information 

### Reco track tree contains information of the reconstructed track, including hit level information
- [X] Event_id
- [X] Chi2 hipotesis for: muons, pions, kaons and protons (the last two should not be seen but just in case) <- From hit based information
- [X] dE/dx (must be able to split it if needed) -> Ready to access, also dQdx
- [X] Residual range
- [X] Lenght 
- [X] Kinetic energy and missing energy 
- [X] Hit information: number of hits, hit_tpc [],  hit_plane
- [X] Precise hit information: 3D track reconstructed ( position, time and calorimetry ): TrajectoryPoint( i )
- [ ] Keep hiearchy information: need to understand why primary tracks are mostly empty?
- [ ] Need plots to further understand this: Chi2 plots?, dEdx, dQdx
  
 # To do List:
 - [ ] Access all information for MC particles
 - [ ] Access all information for RECO particles: PFParticles, Track, Hit information, also Calo information?
 - [ ] Add usefull features and functions for fitter 
 - [ ] Develop fitter itself
 - [ ] Check differences between muons and pions
 - [ ] Implement purity cuts

# To consider:
 -  Save track information in tree?
    -> Requires reading and programing all from scratch outside analyzer
 -  Analyze the track information in analyzer and save in tree number of kinks, DE kink, position? 
    ->  more compact information 
  ->  Requires running analyzer each time 
