%laserScanner.m
%Author: Patrick McKeen
%March 9th, 2016
%For OpenLoop Alliance in 2015-2016 SpaceX Hyperloop Competition

%Free to use and modify, so long as this line and the lines above it
%(lines 1 to 8) are not removed or changed, or the author grants
%permission.

%This works inside the Hyperloop simulation
%Assumes that the laser scanners central scan point:
%XY: [0,1,0]
%XZ: [0,0,-1]
%YZ: [0,0,-1]
function [angle,distance]=laserScanner(sensorPos,plane,podPos,rot,range,scans)
%sensorPos is sensor position in local frame
%plane is the plane that the sensor scans in, in the local frame. XY is
%[1,1,0], YZ is [0,1,1], and XZ is [1,0,1]. 
%podPos is pod position in global frame
%rot is rotation matrix
%range is angle of the of the arc of the laser's max and min sweep, in
%degrees.
%scans is the number of scans in this arc. The default is for each degree

%angle is the list of angles scanned, from the central beam, in degrees
%distance is distance to nearest object for each scan, at the angles in
%angle

%optional input
if nargin==5
    scans=range+1;
end

%angles to scan
angle=linspace(-range/2,range/2,scans);

if cross(plane,[0,0,1])==0 %XY plane
    f=@(x) [sin(x),cos(x) 0]';
    vec=arrayfunc(f,angle);
elseif cross(plane, [0,1,0])==0 %XZ plane
    f=@(x) [sin(x),0,-cos(x)]';
    vec=arrayfunc(f,angle);
elseif cross(Plane,[1,0,0])==0 %YZ plane
    f=@(x) [0,sin(x),-cos(x)]';
    vec=arrayfunc(f,angle);
else
    error('This is not a valid plane. Try [1,1,0] for XY, [0,1,1] for YZ, or [1,0,1] for XZ')
end
wrapper=@(x)arbitraryDistance(sensorPos,x,podPos,rot)';
distance=arrayfun(@(x) [0 1 0 0 0 0 0 0 0 0 0 0;]*wrapper(x),vec);
