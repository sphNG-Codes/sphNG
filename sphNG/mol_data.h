c
c   Written by S. Glover, AMNH, 2004-2005, AIP, 2006
c
c H2 cooling function data -- based on data made available by 
c J. Le Bourlot, G. Pineau des Forets and  D. R. Flower 
c (see Le Bourlot et al, 1999, MNRAS, 305, 802; and also 
c  http://ccp7.dur.ac.uk/cooling_by_h2/).  
c
c The LTE rate is in units of erg s^-1 molecule^-1; the low density 
c rates, on the other hand, are written in terms of erg s^-1 cm^3.
c
c Data in the tables is for log(T) = 2.0 -> 4.0, at intervals of 0.05 dex
c
      integer nh2data
      parameter (nh2data = 41)
      REAL h2_lte(nh2data), h2_h_rate(nh2data), h2_h2_rate(nh2data)
      REAL h2_temp(nh2data)
c
c Low density limit -- almost pure H2 -- H/H2 = 1e-8
c
      DATA h2_h2_rate /-27.604, -27.340, -27.082, -26.856, -26.639, 
     $                 -26.440, -26.254, -26.078, -25.915, -25.756,
     $                 -25.610, -25.466, -25.330, -25.195, -25.065,
     $                 -24.937, -24.812, -24.689, -24.569, -24.453,
     $                 -24.338, -24.230, -24.124, -24.022, -23.923,
     $                 -23.828, -23.735, -23.645, -23.559, -23.475,
     $                 -23.395, -23.316, -23.243, -23.171, -23.103,
     $                 -23.038, -22.976, -22.918, -22.861, -22.809,
     $                 -22.758/
c
c Low density limit -- almost pure H -- H/H2 = 1e6
c
      DATA h2_h_rate /-28.151, -27.865, -27.593, -27.342, -27.107, 
     $                -26.882, -26.669, -26.467, -26.272, -26.085, 
     $                -25.905, -25.730, -25.561, -25.395, -25.235, 
     $                -25.078, -24.923, -24.772, -24.624, -24.478, 
     $                -24.334, -24.192, -24.050, -23.910, -23.772, 
     $                -23.634, -23.497, -23.359, -23.222, -23.086, 
     $                -22.952, -22.820, -22.691, -22.565, -22.446, 
     $                -22.332, -22.225, -22.123, -22.029, -21.943, 
     $                -21.943/ 
c
c LTE rate -- independent of H2/H ratio
c
      DATA h2_lte /-25.440, -25.107, -24.792, -24.503, -24.239, -23.995,
     $             -23.766, -23.550, -23.341, -23.136, -22.933, -22.731,
     $             -22.529, -22.330, -22.134, -21.940, -21.748, -21.555,
     $             -21.350, -21.122, -20.863, -20.575, -20.282, -20.008,
     $             -19.764, -19.552, -19.365, -19.201, -19.056, -18.927,
     $             -18.813, -18.712, -18.624, -18.549, -18.484, -18.429,
     $             -18.382, -18.341, -18.307, -18.278, -18.253/
