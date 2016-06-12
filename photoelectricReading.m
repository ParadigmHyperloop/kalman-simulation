function output = photoelectricReading(sensorPosition,sensorOrientation, podPosition, podRotation, tube) 

%pull out variables of pod state
rx=podPosition(1);
ry=podPosition(2);
rz=podPosition(3);
q1=podRotation(1);  
q2=podRotation(2);
q3=podRotation(3);
q0=podRotation(4);

%generate rotation matrix
Rot=[1-2*q2^2-2*q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);...
 2*(q1*q2+q0*q3) 1-2*q1^2-2*q3^2 2*(q2*q3-q0*q1);...
 2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*q1^2-2*q2^2;];

%find sensor positions in global frame
sx=(Rot(1,:)*sensorPosition + rx);
sy=(Rot(2,:)*sensorPosition + ry);
sz=(Rot(3,:)*sensorPosition + rz);

%useful re-used variables
g=(((sz-tube.tubeCenterToTopOfRail)*Rot(3,:))+(sy*Rot(2,:)));
m=((dot(Rot(2,:),sensorOrientation))^2+(dot(Rot(3,:),sensorOrientation))^2);
tau = ((-dot(g,sensorOrientation))+sqrt((dot(g,sensorOrientation))^2-(m*((sz-tube.tubeCenterToTopOfRail)^2-tube.Radius^2+sy^2))))/m;

%x-coordinate of center of photoelectric sensor detection area on tube wall
alpha=sx+dot(Rot(1,:),sensorOrientation*tau);

dist=dot((Rot*sensorOrientation),tau);

%x-distance of sensor orientation in global frame
sensorOrientationx=dot(Rot(1,:),(sensorOrientation/norm(sensorOrientation)));

%radius of detection areaof PE sensor
detectionRad=dist*sin(tube.angleOfPESensitivity/2);

%detection ranges for alpha to be in, for PE to detect edge of strips
stripLeadingEdgeDistance=alpha-tube.stripDistances-tube.stripWidth/2;
stripTrailingEdgeDistance=alpha-tube.stripDistances+tube.stripWidth/2;

%PE reading
output=(tube.maxBrightness./tube.stripWidth).*sum((stripTrailingEdgeDistance>(-detectionRad)./sensorOrientationx).*(stripTrailingEdgeDistance+detectionRad./sensorOrientationx) - (stripLeadingEdgeDistance>(-detectionRad)./sensorOrientationx).*(stripLeadingEdgeDistance+detectionRad./sensorOrientationx) - (stripTrailingEdgeDistance>(detectionRad)./sensorOrientationx).*(stripTrailingEdgeDistance-detectionRad./sensorOrientationx) + (stripLeadingEdgeDistance>(detectionRad)./sensorOrientationx).*(stripLeadingEdgeDistance-detectionRad./sensorOrientationx));
        