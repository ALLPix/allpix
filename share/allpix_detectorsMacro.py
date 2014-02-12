from math import *

nDetectorsPerDisc = 18; # N detector per disc
nDiscs = 40;

############################  
# all offsets and vars

radius = 100.    # mm;
angleStep = 20.  # deg;
selfAngle = 90.  # deg; // startup angle
posAngle = 0.    # deg; // startup angle
zShift = 10.     # mm;  
posx = 0.   #*mm;
posy = 0.   #*mm;
posz = -200. #*mm;

for detItr in range(0, nDetectorsPerDisc*nDiscs):
    print("# Detector %d"%detItr)
    
    posAngle = detItr*angleStep;
    posx = radius*cos(posAngle*pi/180);
    posy = radius*sin(posAngle*pi/180);
    
#     // traslation
    print("/allpix/det/setPosition %f %f %f mm"%(posx, posy, posz))  # one detector positioned here
#     // rotation
    print("/allpix/det/setRotation %f %f %f deg"%(0, 0, selfAngle))
    print("/allpix/det/setLowTHL 5. keV")

    selfAngle -= angleStep;
    
    if (detItr+1)%nDetectorsPerDisc == 0 : # change disk
       posz += zShift
       selfAngle = 90.
       posAngle = 0.
       
