# TurbDNSPost
Post-processing code for Turbulence DNS:
Cf_old.f90  →  for Skin Friction Coefficient computation
FavAve.f90  →  for Favre average (density weighted average)
KED.f90     →  for Turbulent Kinetic Energy Transportation Equation
ave-channel-2D.f90  →  read averaged data for turbulent flow field
read-channel-3D.f90  →  read instantaneous data for turbulent flow field


DNS data for Channel or Flat-plate turbulence:
Channel-M0.8-R3,000  →  Ma=0.8; Re=3,000 for developed channel turbulence
Channel-M1.5-R3,000  →  Ma=1.5; Re=3,000 for developed channel turbulence
Plate-M0.8-R635,000  →  Ma=0.8; Re=635,000 for developing Flat-plate turbulence
Plate-M2.25-R635,000  →  Ma=2.25; Re=635,000 for developing Flat-plate turbulence
