# Track Analyzer: 

Developed to look for specific features characteristic for muon and pion events -> develop a new PID method

## OUTPUT file: 

RooT file containging three trees,

## Information in each tree: 

### Event Tree contains general information of the event regardless of its nature ( needed ?? )

- [X] event id
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
- [X] Keep hiearchy information: number and pdg of particles stored
- [ ] Need plots to further understand this: Chi2 plots?, dEdx, dQdx


# TrackFitter Class

## Method 1: using local linearity
- Calculates Pearson coefficient per hit regarding the surrounding hits within a specified window ( eg: taking 15 hits preceeding and after hit_i )
- r < 0.9 defined as a clear kink. The algorithm works for simple cases
- Only finds clear kinks: if the kinked track is close to the "mother" track, no deviation is found

## Method 2: should consider direction of track 

  
 # To do List:
 - [X] Access all information for MC particles
 - [X] Access all information for RECO particles: PFParticles, Track, Hit information, also Calo information?
 - [ ] Add usefull features and functions for fitter 
 - [ ] IN DEVELOPEMENT -> Develop fitter itself
 - [ ] Check differences between muons and pions
 - [ ] Implement purity cuts
