function [underSideDistance,beamDistance] = InsideDistance(point)
    %Using dimensions from subtrack 1

    % constants of the beam

    topHeight=2*.412;
    middleWidth=.313;
    distanceToCenter=23.4;
       
    %Distances to beam
    underSideDistance=(-1)*(distanceToCenter+topHeight)-point(3);
    beamDistance=abs(point(2))-(middleWidth/2);

          