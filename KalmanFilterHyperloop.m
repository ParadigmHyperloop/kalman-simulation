function [predictedState, predictedCovariance] = KalmanFilterHyperloop( prevState, prevCovariance, IMUData, sensorData, execution, timestep )

globals = globalData();
pod = podData();
tube = tubeData();

sensorPositions = zeros(7,3);
sensorPositions(1,:) = [0,0,-0.5];
sensorDirections = zeros(7,3);
sensorDirections(1,:) = [0,0,-1];
thicknessOfRail = .05;
tubeCenterToTopOfRail=.4;
maxBrightness=1;
angleOfPESensitivity=10*pi/180;

stripWidth=2*0.0254;

%PE sensor formula assumes the reflectivities above, and
%that the photodiodes are sensitive in a cone of angle angleOfPESensitivity
%around the normal. It assumes that the signal received by each photodiode
%is, when considering the intersection of the tube and the cone traceed out
%by the diode's sensitivity (which is treated as a plane), that the signal is equal to the reflectivty of
%the tape times the area of the intersection that is tape, plus the
%reflectivity of the tube, times the area of the intersection that is tube.
%
%Note that this treats the strips as going the full 2pi radians around the
%tube and ignores color differences.
%
%This also assumes that the cone traced by the photodiode is small enough
%that it is impossible for the photodiode to recieve reflections from two
%strips at once, in any realistic position of the pod.
%
%It also plays fast and loose with the geometry of cones.
%
%This is a believable, but very ideal, set of assumptions. A more correct
%formula will be instituted after sensor tests of the photodiodes,
%including the use of two photodiodes and their difference as the signal,
%hence the placeholder variable distanceBetweenPE;

p1 = sensorPositions(1,:)'; %Position of the sensor with respect to the center of mass of the pod
p2 = sensorPositions(2,:)'; %Position of the sensor with respect to the center of mass of the pod
p3 = sensorPositions(3,:)'; %Position of the sensor with respect to the center of mass of the pod
%tck = thicknessOfRail/2;
p4 = sensorPositions(4,:)'; %Position of the sensor with respect to the center of mass of the pod
p5 = sensorPositions(4,:)'; %Position of the sensor with respect to the center of mass of the pod
p6 = sensorPositions(4,:)'; %Position of the sensor with respect to the center of mass of the pod
p7 = sensorPositions(4,:)'; %Position of the sensor with respect to the center of mass of the pod


b4 = sensorDirections(4,:)'; %Ray indicating direction sensor points in local coordinates
b5 = sensorDirections(5,:)'; %Ray indicating direction sensor points in local coordinates
b6 = sensorDirections(6,:)'; %Ray indicating direction sensor points in local coordinates
b7 = sensorDirections(7,:)'; %Ray indicating direction sensor points in local coordinates


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing state, kalman gain etc%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Currently using IMU as control to make the state prediction and so
%initializing that

%%%%%%%%%%%%%%%% ----------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xkk=prevState; %This is the previous state and in the case of k=1 the initial state
uk=IMUData; %Control vector, currently IMU being used
Pkk=prevCovariance;

%
axL=uk(1);  %accelaration in x from IMU
ayL=uk(2);  %acceleration in y from IMU
azL=uk(3);  %acceleration in z from IMU
omegaX=uk(4);   %angular velocity in x from IMU
omegaY=uk(5);   %angular velocity in y from IMU
omegaZ=uk(6);   %angular velocity in z from IMU
%

%quaternions from the state
q1=xkk(7);  
q2=xkk(8);
q3=xkk(9);
q0=xkk(10);

%rotation matrix from the quaternions
Rot=[1-2*q2^2-2*q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);...
     2*(q1*q2+q0*q3) 1-2*q1^2-2*q3^2 2*(q2*q3-q0*q1);...
     2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*q1^2-2*q2^2;];


 %Fk matrix used for predicting the next state and the error in the
 %next state
%  Fk = [1 0 0 globals.kalmanTimestep 0 0 0 0 0 0;...
%        0 1 0 0 globals.kalmanTimestep 0 0 0 0 0;...
%        0 0 1 0 0 globals.kalmanTimestep 0 0 0 0;...
%        0 0 0 1 0 0 0 0 0 0;...
%        0 0 0 0 1 0 0 0 0 0;...
%        0 0 0 0 0 1 0 0 0 0;...
%        0 0 0 0 0 0 1 0 0 0;...
%        0 0 0 0 0 0 0 1 0 0;...
%        0 0 0 0 0 0 0 0 1 0;...
%        0 0 0 0 0 0 0 0 0 1;];

   Fk = [1 0 0 globals.kalmanTimestep 0 0 0 0 0 0;...
       0 1 0 0 globals.kalmanTimestep 0 0 0 0 0;...
       0 0 1 0 0 globals.kalmanTimestep 0 0 0 0;...
       0 0 0 1 0 0 globals.kalmanTimestep*(2*q2*ayL+2*q3*azL) globals.kalmanTimestep*(-4*q2*axL+2*q1*ayL+2*q0*azL) globals.kalmanTimestep*(-4*q3*axL-2*q0*ayL+2*q1*azL) globals.kalmanTimestep*(-2*q3*ayL+2*q2*azL);...
       0 0 0 0 1 0 globals.kalmanTimestep*(2*q2*axL-4*q1*ayL-2*q0*azL) globals.kalmanTimestep*(2*q1*axL+2*q3*azL) globals.kalmanTimestep*(2*q0*axL-4*q3*ayL+2*q2*azL) globals.kalmanTimestep*(2*q3*axL-2*q1*azL);...
       0 0 0 0 0 1 globals.kalmanTimestep*(2*q3*axL+2*q0*ayL-4*q1*azL) globals.kalmanTimestep*(-2*q0*axL+2*q3*ayL-4*q2*azL) globals.kalmanTimestep*(2*q1*axL+2*q2*ayL) globals.kalmanTimestep*(-2*q2*axL+2*q1*ayL);...
       0 0 0 0 0 0 1 globals.kalmanTimestep*omegaZ/2 -globals.kalmanTimestep*omegaY/2 globals.kalmanTimestep*omegaX/2;...
       0 0 0 0 0 0 -globals.kalmanTimestep*omegaZ/2 1 globals.kalmanTimestep*omegaX/2 globals.kalmanTimestep*omegaY/2;...
       0 0 0 0 0 0 globals.kalmanTimestep*omegaY/2 -globals.kalmanTimestep*omegaX/2 1 globals.kalmanTimestep*omegaZ/2;...
       0 0 0 0 0 0 -globals.kalmanTimestep*omegaX/2 -globals.kalmanTimestep*omegaY/2 -globals.kalmanTimestep*omegaZ/2 1;];
 



 %Bk matrix, multiplied to the control vector (currently IMU) in order
 %to also predict the next state of the pod
 Bk=[ 0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      Rot(1,1) Rot(1,2) Rot(1,3) 0 0 0;...
      Rot(2,1) Rot(2,2) Rot(2,3) 0 0 0;...
      Rot(3,1) Rot(3,2) Rot(3,3) 0 0 0;... %drop gravity term (constants=0 for this Jacobian)
      0 0 0 q0/2 -q3/2 q2/2;...
      0 0 0 q3/2 q0/2 -q1/2;...
      0 0 0 -q2/2 q1/2 q0/2;...
      0 0 0 -q1/2 -q2/2 -q3/2;]*globals.kalmanTimestep; 

%This is the portion for Wk and this is the same as patricks E80 code 
%for rn. We will be expirimentally determining this next semester%

%OLDuk=Control(:,k);
OmegaMatrix=[0 omegaZ -omegaY omegaX;...
            -omegaZ 0 omegaX omegaY;...
            omegaY -omegaX 0 omegaZ;...
            -omegaX -omegaY -omegaZ 0;];
% dqdt=(1/2).*(OmegaMatrix*xkk(7:10));
% dq1dt=dqdt(1);
% dq2dt=dqdt(2);
% dq3dt=dqdt(3);
% dq0dt=dqdt(4);
% domegadt=(NEXTuk(4:6)-OLDuk(4:6))./(2*globals.kalmanTimestep);
% domegaXdt=domegadt(1);
% domegaYdt=domegadt(2);
% domegaZdt=domegadt(3);


% wk=(0.5*globals.kalmanTimestep^2).*([(Rot*[axL;ayL;azL;])' ([1-4*q2*dq2dt-4*q3*dq3dt 2*(q1*dq2dt+q2*dq1dt-q0*dq3dt-q3*dq0dt) 2*(q1*dq3dt+q3*dq1dt+q0*dq2dt+q2*dq0dt);... %([1-4*q2*dq2dt-4*q3*dq3dt 2*(q1*dq2dt+q2*dq1dt+q0*dq3dt+q3*dq0dt) 2*(-q0*dq2dt-q2*dq0dt+q1*dq3dt+q3*dq1dt);%2*(-q0*dq3dt-q3*dq0dt+q1*dq2dt+q2*dq1dt) 1-4*q1*dq1dt-4*q3*dq3dt 2*(q2*dq3dt+q3*dq2dt+q0*dq1dt+q1*dq0dt); 2*(q1*dq3dt+q3*dq1dt+q0*dq2dt+q2*dq0dt) 2*(-q0*dq1dt-q1*dq0dt+q2*dq3dt+q3*dq2dt) 1-4*q2*dq2dt-4*q1*dq1dt;]*[axL;ayL;azL;]+Rot*([NEXTaxL-OLDaxL;NEXTayL-OLDayL;NEXTazL-OLDazL;]./(2*globals.kalmanTimestep)))'  
%                  2*(q1*dq2dt+q2*dq1dt+q0*dq3dt+q3*dq0dt) 1-4*q1*dq1dt-4*q3*dq3dt 2*(q2*dq3dt+q3*dq2dt-q0*dq1dt-q1*dq0dt);...
%                  2*(q1*dq3dt+q3*dq1dt-q0*dq2dt-q2*dq0dt) 2*(q2*dq3dt+q3*dq2dt+q0*dq1dt+q1*dq0dt) 1-4*q2*dq2dt-4*q1*dq1dt;]*[axL;ayL;azL;]+Rot*([NEXTaxL-OLDaxL;NEXTayL-OLDayL;NEXTazL-OLDazL;]./(2*globals.kalmanTimestep)))' ((1/2)*(OmegaMatrix*dqdt+...
%       [0 domegaZdt -domegaYdt domegaXdt;...
%       -domegaZdt 0 domegaXdt domegaYdt;...
%       domegaYdt -domegaXdt 0 domegaZdt;...
%       -domegaXdt -domegaYdt -domegaZdt 0;]*xkk(7:10)))'
%     ]').^2;
% Wk=(wk*wk')*(k~=1)+diag([wk])*(k==1);%diag(wk);
%Experimental determination bit ends%

Wk=zeros(10,10);


%The derivative of the rotation matrix in terms of quaternions

factor = 0.1;
Qk = factor*eye(6);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREDICTION STEP

xkp1k=xkk+globals.kalmanTimestep*[  xkk(4:6);...
                       (Rot*uk(1:3))-[0;0;globals.gravity;];...
                       (1/2).*(OmegaMatrix*xkk(7:10));];  %prediction step of the state 

Pkp1k=Fk*Pkk*Fk'+Bk*Qk*Bk'+Wk; %prediction step of the error, needs to be experimentally determined 


normQuat=sqrt(sum((xkp1k(7:10)).^2));
xkp1k(7:10)=xkp1k(7:10)./normQuat;


%%%%%%%%%%%%%%%%%%%----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECTIVE STEPS

% CORRECTIVE STEP 1, XZ LASER SCANNER
if execution(1)==1
    
    z1kp1 = sensorData(:,1);
    
    rz1=xkp1k(3);
    q11=xkp1k(7);  
    q12=xkp1k(8);
    q13=xkp1k(9);
    q10=xkp1k(10);
    Rot1=[1-2*q12^2-2*q13^2 2*(q11*q12-q10*q13) 2*(q11*q13+q10*q12);...
     2*(q11*q12+q10*q13) 1-2*q11^2-2*q13^2 2*(q12*q13-q10*q11);...
     2*(q11*q13-q10*q12) 2*(q12*q13+q10*q11) 1-2*q11^2-2*q12^2;];
    sz1=(Rot1(3,:)*p1+rz1+0.72136);
    sq1=sqrt(Rot1(2,2)^2+Rot1(1,2)^2);
    
    dd1dq1=(sq1*([2*q13 2*q10 -4*q11]*p1)-sz1*(Rot1(2,2)*(-4*q11) + Rot1(1,2)*2*q12)/sq1)/(sq1^2);
    dd1dq2=(sq1*([-2*q10 2*q13 -4*q12]*p1)-sz1*(Rot1(2,2)*0 + Rot1(1,2)*2*q11)/sq1)/(sq1^2);
    dd1dq3=(sq1*([2*q11 2*q12 0]*p1)-sz1*(Rot1(2,2)*(-4*q13) + Rot1(1,2)*-2*q10)/sq1)/(sq1^2);
    dd1dq0=(sq1*([-2*q12 2*q11 0]*p1)-sz1*(Rot1(2,2)*0 + Rot1(1,2)*-2*q13)/sq1)/(sq1^2);
    
    dtheta1dq1=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*((-4*q11)*Rot1(2,1)+2*q12*Rot1(1,2)-0*Rot1(2,2)-2*q12*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*(-4*q11) + Rot1(1,2)*2*q12)/sq1)/(sq1^2);
    dtheta1dq2=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*(0*Rot1(2,1)+2*q11*Rot1(1,2)-q12*-4*Rot1(2,2)-2*q11*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*0 + Rot1(1,2)*2*q11)/sq1)/(sq1^2);
    dtheta1dq3=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*(-4*q13*Rot1(2,1)+2*q10*Rot1(1,2)-q13*-4*Rot1(2,2)-q10*-2*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*(-4*q13) + Rot1(1,2)*-2*q10)/sq1)/(sq1^2);
    dtheta1dq0=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*(0*Rot1(2,1)+2*q13*Rot1(1,2)-0*Rot1(2,2)-q13*-2*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*0 + Rot1(1,2)*-2*q13)/sq1)/(sq1^2);
    
    
    H1kp1=[0 0 1/sq1 0 0 0 dd1dq1 dd1dq2 dd1dq3 dd1dq0;...
            0 0 0 0 0 0 dtheta1dq1 dtheta1dq2 dtheta1dq3 dtheta1dq0;];
%     S1kp1=diag([0.01,1000]);%XZScannerCovariance; %Experimentally determined
%     laserFactor = 100;
    S1kp1 = globals.laserCovariance*eye(2);
%     disp('--------------------------------------')
%     disp(Pkp1k);
%     disp(H1kp1);
%     disp(H1kp1*Pkp1k*H1kp1');
%     disp(H1kp1*Pkp1k*H1kp1'+S1kp1);
    K1kp1=Pkp1k*H1kp1'/(H1kp1*Pkp1k*H1kp1'+S1kp1);


    h1kp1=[sz1/sq1;...
           acos((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1);];
       
    x1kp1kp1=xkp1k+K1kp1*(z1kp1-h1kp1);
    P1kp1kp1=(eye(10,10)-K1kp1*H1kp1)*Pkp1k;
    
    normQuat1=sqrt(sum((x1kp1kp1(7:10)).^2));
    x1kp1kp1(7:10)=x1kp1kp1(7:10)./normQuat1;
else
    x1kp1kp1=xkp1k;
    P1kp1kp1=Pkp1k;
end



% CORRECTIVE STEP 2, YZ LASER SCANNER
if execution(2)==1
    
    z2kp1 = sensorData(2,:);
        
    rz2=x1kp1kp1(3);
    q21=x1kp1kp1(7);  
    q22=x1kp1kp1(8);
    q23=x1kp1kp1(9);
    q20=x1kp1kp1(10);
    Rot2=[1-2*q22^2-2*q23^2 2*(q21*q22-q20*q23) 2*(q21*q23+q20*q22);...
     2*(q21*q22+q20*q23) 1-2*q21^2-2*q23^2 2*(q22*q23-q20*q21);...
     2*(q21*q23-q20*q22) 2*(q22*q23+q20*q21) 1-2*q21^2-2*q22^2;];
    sz2=(Rot2(3,:)*p2+rz2);
    sq2=sqrt(Rot2(2,1)^2+Rot2(1,1)^2);
    
    
    dd2dq1=(sq2*([2*q23 2*q20 -4*q21]*p2)-sz2*(Rot2(2,1)*(2*q22) + Rot2(1,1)*0)/sq2)/(sq2^2);
    dd2dq2=(sq2*([-2*q20 2*q23 -4*q22]*p2)-sz2*(Rot2(2,1)*2*q21 + Rot2(1,1)*-4*q22)/sq2)/(sq2^2);
    dd2dq3=(sq2*([2*q21 2*q22 0]*p2)-sz2*(Rot2(2,1)*(2*q20) + Rot2(1,1)*-4*q23)/sq2)/(sq2^2);
    dd2dq0=(sq2*([-2*q22 2*q21 0]*p2)-sz2*(Rot2(2,1)*2*q23 + Rot2(1,1)*0)/sq2)/(sq2^2);
    
    dtheta2dq1=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*((-4*q21)*Rot2(2,1)+2*q22*Rot2(1,2)-0*Rot2(2,2)-2*q22*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,2)*(Rot2(2,1)*(2*q22) + Rot2(1,1)*0)/sq2)/(sq2^2));
    dtheta2dq2=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*(0*Rot2(2,1)+2*q21*Rot2(1,2)-q22*-4*Rot2(2,2)-2*q21*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,1)*2*q21 + Rot2(1,1)*-4*q22)/sq2)/(sq2^2);
    dtheta2dq3=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*(-4*q23*Rot2(2,1)+2*q20*Rot2(1,2)-q23*-4*Rot2(2,2)-q20*-2*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,1)*(2*q20) + Rot2(1,1)*-4*q23)/sq2)/(sq2^2);
    dtheta2dq0=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*(0*Rot2(2,1)+2*q23*Rot2(1,2)-0*Rot2(2,2)-q23*-2*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,1)*2*q23 + Rot2(1,1)*0)/sq2)/(sq2^2);
    
    
    H2kp1=[0 0 1/sq2 0 0 0 dd2dq1 dd2dq2 dd2dq3 dd2dq0;...
            0 0 0 0 0 0 dtheta2dq1 dtheta2dq2 dtheta2dq3 dtheta2dq0;];
    S2kp1=YZScannerCovariance; %Experimentally determined
    K2kp1=P1kp1kp1*H2kp1'/(H2kp1*P1kp1kp1*H2kp1'+S2kp1);

    
    h2kp1=[sz2/sq2;...
           acos((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2);];
        
    x2kp1kp1=x1kp1kp1+K2kp1*(z2kp1-h2kp1);
    P2kp1kp1=(eye(10,10)-K2kp1*H2kp1)*P1kp1kp1;
    
    normQuat2=sqrt(sum((x2kp1kp1(7:10)).^2));
    x2kp1kp1(7:10)=x2kp1kp1(7:10)./normQuat2;
else
    x2kp1kp1=x1kp1kp1;
    P2kp1kp1=P1kp1kp1;
end



% CORRECTIVE STEP 3, XY LASER SCANNER
if execution(3)==1
    
    z3kp1 = sensorData(3,:);
        
    ry3=x2kp1kp1(2);
    q31=x2kp1kp1(7);  
    q32=x2kp1kp1(8);
    q33=x2kp1kp1(9);
    q30=x2fkp1kp1(10);
    Rot3=[1-2*q32^2-2*q33^2 2*(q31*q32-q30*q33) 2*(q31*q33+q30*q32);...
     2*(q31*q32+q30*q33) 1-2*q31^2-2*q33^2 2*(q32*q33-q30*q31);...
     2*(q31*q33-q30*q32) 2*(q32*q33+q30*q31) 1-2*q31^2-2*q32^2;];
    sy3=(Rot3(2,:)*p3+ry3);
    sq3=sqrt(Rot3(3,3)^2+Rot3(1,3)^2);
    m3=Rot3(1,1)*Rot3(3,3)-Rot3(1,3)*Rot3(3,1);
    
    
    dd3dq1=-(thicknessOfRail/4)*((sq3)^(-3))*(2*Rot3(3,3)*-4*q31+2*Rot3(1,3)*2*q33)-(sq3*([2*q32 -4*q31 -2*q30]*p3)-0.5*(sy3*(2*Rot3(3,3)*-4*q31+2*Rot3(1,3)*2*q33)/sq3))/(sq3^2);
    dd3dq2=-(thicknessOfRail/4)*((sq3)^(-3))*(2*Rot3(3,3)*-4*q32+2*Rot3(1,3)*2*q30)-(sq3*([2*q31 0 2*q33]*p3)-0.5*(sy3*(2*Rot3(3,3)*-4*q32+2*Rot3(1,3)*2*q30)/sq3))/(sq3^2);
    dd3dq3=-(thicknessOfRail/4)*((sq3)^(-3))*(2*Rot3(3,3)*0+2*Rot3(1,3)*2*q31)-(sq3*([2*q30 -4*q33 2*q32]*p3)-0.5*(sy3*(2*Rot3(3,3)*0+2*Rot3(1,3)*2*q31)/sq3))/(sq3^2);
    dd3dq0=-(thicknessOfRail/4)*((sq3)^(-3))*(2*Rot3(3,3)*0+2*Rot3(1,3)*2*q32)-(sq3*([2*q33 0 -2*q31]*p3)-0.5*(sy3*(2*Rot3(3,3)*-2*q31+2*Rot3(1,3)*2*q32)/sq3))/(sq3^2);
    
    dtheta3dq1=-(((Rot3(1,1)*-4*q31+1*-4*q31*Rot3(3,3))-(Rot3(1,3)*2*q33+2*q33*Rot3(3,1)))*sq3-((0.5/sq3)*(2*Rot3(3,3)*-4*q31+2*Rot3(1,3)*2*q33))*m3)/((sq3^2)*sqrt(1-(m3/sq3)^2));
    dtheta3dq2=-(((Rot3(1,1)*-4*q32+1*-4*q32*Rot3(3,3))-(Rot3(1,3)*-2*q30+2*q30*Rot3(3,1)))*sq3-((0.5/sq3)*(2*Rot3(3,3)*-4*q32+2*Rot3(1,3)*2*q30))*m3)/((sq3^2)*sqrt(1-(m3/sq3)^2));
    dtheta3dq3=-(((Rot3(1,1)*0+1*-4*q33*Rot3(3,3))-(Rot3(1,3)*2*q31+2*q31*Rot3(3,1)))*sq3-((0.5/sq3)*(2*Rot3(3,3)*0+2*Rot3(1,3)*2*q31))*m3)/((sq3^2)*sqrt(1-(m3/sq3)^2));
    dtheta3dq0=-(((Rot3(1,1)*0+0*Rot3(3,3))-(Rot3(1,3)*-2*q32+2*q32*Rot3(3,1)))*sq3-((0.5/sq3)*(2*Rot3(3,3)*0+2*Rot3(1,3)*2*q32))*m3)/((sq3^2)*sqrt(1-(m3/sq3)^2));
    
    
    H3kp1=[0 1/sq3 0 0 0 0 dd3dq1 dd3dq2 dd3dq3 dd3dq0;...
            0 0 0 0 0 0 dtheta3dq1 dtheta3dq2 dtheta3dq3 dtheta3dq0;];
    S3kp1=XYScannerCovariance; %Experimentally determined
    K3kp1=P2kp1kp1*H3kp1'/(H3kp1*P2kp1kp1*H3kp1'+S3kp1);

    
    h3kp1=[((thicknessOfRail/2)-sy3)/sq3;
           acos((Rot3(1,1)*Rot3(3,3)-Rot3(1,3)*Rot3(3,1))/sq3);];
        
    x3kp1kp1=x2kp1kp1+K3kp1*(z3kp1-h3kp1);
    P3kp1kp1=(eye(10,10)-K3kp1*H3kp1)*P2kp1kp1;
    
    normQuat3=sqrt(sum((x3kp1kp1(7:10)).^2));
    x3kp1kp1(7:10)=x3kp1kp1(7:10)./normQuat3;
else
    x3kp1kp1=x2kp1kp1;
    P3kp1kp1=P2kp1kp1;
end


% CORRECTIVE STEP 4, PITOT TUBE
if execution(4)==1
    
    z4kp1 = sensorData(4,:);
    
    vPod4 = x3kp1kp1(4:6);
    q41=x3kp1kp1(7);  
    q42=x3kp1kp1(8);
    q43=x3kp1kp1(9);
    q40=x3fkp1kp1(10);
    Rot4=[1-2*q42^2-2*q43^2 2*(q41*q42-q40*q43) 2*(q41*q43+q40*q42);...
     2*(q41*q42+q40*q43) 1-2*q41^2-2*q43^2 2*(q42*q43-q40*q41);...
     2*(q41*q43-q40*q42) 2*(q42*q43+q40*q41) 1-2*q41^2-2*q42^2;];
 
    dRot4dq1=[0 2*q42 2*q43;...
     2*q42 -4*q41 -2*q40;...
     2*q43 2*q40 -4*q41;];
 
    dRot4dq2=[-4*q42 2*q41 2*q40;...
     2*q41 0 2*q43;...
     -2*q40 2*q43 -4*q42;];
 
    dRot4dq3=[-4*q43 -2*q40 2*q41;...
     2*q40 -4*q43 2*q42;...
     2*q41 2*q42 0;];
 
    dRot4dq0=[0 -2*q43 2*q42;...
     2*q43 0 -2*q41;...
     -2*q42 2*q41 0;];
    
    
    ddP4dq1=densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot(vPod4,(dRot4dq1*b4)));
    ddP4dq2=densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot(vPod4,(dRot4dq2*b4)));
    ddP4dq3=densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot(vPod4,(dRot4dq3*b4)));
    ddP4dq0=densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot(vPod4,(dRot4dq0*b4)));
    
    ddP4dvx = densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot((Rot4*b4),[1 0 0]'));
    ddP4dvy = densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot((Rot4*b4),[0 1 0]'));
    ddP4dvz = densityOfAir*(dot(vPod4,(Rot4*b4)))*(dot((Rot4*b4),[0 0 1]'));
    
    
    H4kp1=[0 0 0 ddP4dvx ddP4dvy ddP4dvz ddP4dq1 ddP4dq2 ddP4dq3 ddP4dq0];
    S4kp1=PitotCovariance; %Experimentally determined
    K4kp1=P3kp1kp1*H4kp1'/(H4kp1*P3kp1kp1*H4kp1'+S4kp1);

    
    h4kp1=0.5*densityOfAir*(dot(vPod4,b4))^2;
        
    x4kp1kp1=x3kp1kp1+K4kp1*(z4kp1-h4kp1);
    P4kp1kp1=(eye(10,10)-K4kp1*H4kp1)*P3kp1kp1;
    
    normQuat4=sqrt(sum((x4kp1kp1(7:10)).^2));
    x4kp1kp1(7:10)=x4kp1kp1(7:10)./normQuat4;
else
    x4kp1kp1=x3kp1kp1;
    P4kp1kp1=P3kp1kp1;
end



% CORRECTIVE STEP 5, 12 o'clock Photoelectric
if execution(5)==1
    % getting sensor data
    z5kp1 = sensorData(5,:);
    
    % Getting terms from the previous state prediction to use in calculation
    rx5=x4kp1kp1(1);
    ry5=x4kp1kp1(2);
    rz5=x4kp1kp1(3);
    q51=x4kp1kp1(7);  
    q52=x4kp1kp1(8);
    q53=x4kp1kp1(9);
    q50=x4kp1kp1(10);
    %Calculating the rotation matrix
    Rot5=[1-2*q52^2-2*q53^2 2*(q51*q52-q50*q53) 2*(q51*q53+q50*q52);...
     2*(q51*q52+q50*q53) 1-2*q51^2-2*q53^2 2*(q52*q53-q50*q51);...
     2*(q51*q53-q50*q52) 2*(q52*q53+q50*q51) 1-2*q51^2-2*q52^2;];
    sx5=(Rot5(1,:)*p5 + rx5);
    sy5=(Rot5(2,:)*p5 + ry5);
    sz5=(Rot5(3,:)*p5 + rz5);
    
    %derivatives of the rotation matrix w.r.t. each quaternion to be used to calculate the jacobian
    dRot5dq1=[0 2*q52 2*q53;...
     2*q52 -4*q51 -2*q50;...
     2*q53 2*q50 -4*q51;];
 
    dRot5dq2=[-4*q52 2*q51 2*q50;...
     2*q51 0 2*q53;...
     -2*q50 2*q53 -4*q52;];
 
    dRot5dq3=[-4*q53 -2*q50 2*q51;...
     2*q50 -4*q53 2*q52;...
     2*q51 2*q52 0;];
 
    dRot5dq0=[0 -2*q53 2*q52;...
     2*q53 0 -2*q51;...
     -2*q52 2*q51 0;];
     
    %mostly intermediate terms which are used over and over again so its easier to define them
    g5=(((sz5-tubeCenterToTopOfRail)*Rot5(3,:))+(sy5*Rot5(2,:)));
    dg5dq1 = (sz5-tubeCenterToTopOfRail)*dRot5dq1(3,:) + (dRot5dq1(3,:)*p5)*Rot5(3,:) + sy5*dRot5dq1(2,:) + (dRot5dq1(2,:)*p5)*Rot5(2,:);
    dg5dq2 = (sz5-tubeCenterToTopOfRail)*dRot5dq2(3,:) + (dRot5dq2(3,:)*p5)*Rot5(3,:) + sy5*dRot5dq2(2,:) + (dRot5dq2(2,:)*p5)*Rot5(2,:);
    dg5dq3 = (sz5-tubeCenterToTopOfRail)*dRot5dq3(3,:) + (dRot5dq3(3,:)*p5)*Rot5(3,:) + sy5*dRot5dq3(2,:) + (dRot5dq3(2,:)*p5)*Rot5(2,:);
    dg5dq0 = (sz5-tubeCenterToTopOfRail)*dRot5dq0(3,:) + (dRot5dq0(3,:)*p5)*Rot5(3,:) + sy5*dRot5dq0(2,:) + (dRot5dq0(2,:)*p5)*Rot5(2,:);
    %more intermediate terms and their partial derivative terms
    m5=((dot(Rot5(2,:),b5))^2+(dot(Rot5(3,:),b5))^2);
    dm5dq1 = 2*dot(b5,((dot(b5,Rot5(2,:)))*dRot5dq1(2,:) + (dot(b5,Rot5(3,:)))*dRot5dq1(3,:)));
    dm5dq2 = 2*dot(b5,((dot(b5,Rot5(2,:)))*dRot5dq2(2,:) + (dot(b5,Rot5(3,:)))*dRot5dq2(3,:)));
    dm5dq3 = 2*dot(b5,((dot(b5,Rot5(2,:)))*dRot5dq3(2,:) + (dot(b5,Rot5(3,:)))*dRot5dq3(3,:)));
    dm5dq0 = 2*dot(b5,((dot(b5,Rot5(2,:)))*dRot5dq0(2,:) + (dot(b5,Rot5(3,:)))*dRot5dq0(3,:)));
 
    tau5 = ((-dot(g5,b5))+sqrt((dot(g5,b5))^2-(m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))))/m5;
    %partial derivative of tau w.r.t. the position
    dtau5drx = 0;
    dtau5dry = (-dot(Rot5(2,:),b5) + ((dot(g5,b5))*(dot(b5,Rot5(2,:)))-m5*sy5)/sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2)))/m5;
    dtau5drz = (-dot(Rot5(3,:),b5) + ((dot(g5,b5))*(dot(b5,Rot5(3,:)))-m5*(sz5-tubeCenterToTopOfRail))/sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2)))/m5;
    %partial derivative of tau w.r.t. the quaternions
    dtau5dq1=m5*(-dot(b5, ((dot(p5,dRot5dq1(3,:))*Rot5(3,:) +(sz5-tubeCenterToTopOfRail)*dRot5dq1(3,:)+(dot(p5,dRot5dq1(2,:)))*Rot5(2,:)+sy5*dRot5dq1(2,:))) + (2*dot(g5,b5)*dot(dg5dq1,b5)-(2*m5*(sz5-tubeCenterToTopOfRail)*(dRot5dq1(3,:)*p5) + ((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))*dm5dq1))/(2*sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))))-((-(dot(g5,b5))+sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2)))*dm5dq1)/(m5^2);
    dtau5dq2=m5*(-dot(b5, ((dot(p5,dRot5dq2(3,:))*Rot5(3,:) +(sz5-tubeCenterToTopOfRail)*dRot5dq2(3,:)+(dot(p5,dRot5dq2(2,:)))*Rot5(2,:)+sy5*dRot5dq2(2,:))) + (2*dot(g5,b5)*dot(dg5dq2,b5)-(2*m5*(sz5-tubeCenterToTopOfRail)*(dRot5dq2(3,:)*p5) + ((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))*dm5dq2))/(2*sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))))-((-(dot(g5,b5))+sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2)))*dm5dq2)/(m5^2);
    dtau5dq3=m5*(-dot(b5, ((dot(p5,dRot5dq3(3,:))*Rot5(3,:) +(sz5-tubeCenterToTopOfRail)*dRot5dq3(3,:)+(dot(p5,dRot5dq3(2,:)))*Rot5(2,:)+sy5*dRot5dq3(2,:))) + (2*dot(g5,b5)*dot(dg5dq3,b5)-(2*m5*(sz5-tubeCenterToTopOfRail)*(dRot5dq3(3,:)*p5) + ((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))*dm5dq3))/(2*sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))))-((-(dot(g5,b5))+sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2)))*dm5dq3)/(m5^2);
    dtau5dq0=m5*(-dot(b5, ((dot(p5,dRot5dq0(3,:))*Rot5(3,:) +(sz5-tubeCenterToTopOfRail)*dRot5dq0(3,:)+(dot(p5,dRot5dq0(2,:)))*Rot5(2,:)+sy5*dRot5dq0(2,:))) + (2*dot(g5,b5)*dot(dg5dq0,b5)-(2*m5*(sz5-tubeCenterToTopOfRail)*(dRot5dq0(3,:)*p5) + ((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))*dm5dq0))/(2*sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2))))-((-(dot(g5,b5))+sqrt((dot(g5,b5))^2-m5*((sz5-tubeCenterToTopOfRail)^2-tubeRadius^2+sy5^2)))*dm5dq0)/(m5^2);
   
 
    alpha5=sx5+dot(Rot5(1,:),b5*tau5);
    %partial derivative of alpha w.r.t. the quaternions
    dalpha5drx = 1+dot(Rot5(1,:),b5*dtau5drx);
    dalpha5dry = dot(Rot5(1,:),b5*dtau5dry);
    dalpha5drz = dot(Rot5(1,:),b5*dtau5dry);
    %partial derivative of alpha w.r.t. the quaternions
    dalpha5dq1 = dRot5dq1(1,:)*p5 + dot(Rot5(1,:),b5*dtau5dq1) + dot(dRot5dq1(1,:),b5*tau5);
    dalpha5dq2 = dRot5dq2(1,:)*p5 + dot(Rot5(1,:),b5*dtau5dq2) + dot(dRot5dq2(1,:),b5*tau5);
    dalpha5dq3 = dRot5dq3(1,:)*p5 + dot(Rot5(1,:),b5*dtau5dq3) + dot(dRot5dq3(1,:),b5*tau5);
    dalpha5dq0 = dRot5dq0(1,:)*p5 + dot(Rot5(1,:),b5*dtau5dq0) + dot(dRot5dq0(1,:),b5*tau5);
    
    dist5=dot((Rot5*b5),tau5);
    
    bx5=dot(Rot5(1,:),(b5/norm(b5)));
    %partial derivative of b w.r.t. the position
    dbx5drx=0;
    dbx5dry=0;
    dbx5drz=0;
    %partial derivative of b w.r.t. the quaternions
    dbx5dq1=dot(dRot5dq1(1,:),(b5/norm(b5)));
    dbx5dq2=dot(dRot5dq2(1,:),(b5/norm(b5)));
    dbx5dq3=dot(dRot5dq3(1,:),(b5/norm(b5)));
    dbx5dq0=dot(dRot5dq0(1,:),(b5/norm(b5)));
    
    detectionRad5=dist5*sin(angleOfPESensitivity/2);
    %partial derivative of detection radius w.r.t. the position
    ddetRad5drx=sin(angleOfPESensitivity/2)*(dot((dRot5drx*b5),tau5)+dot(dtau5drx,(Rot5*b5)));
    ddetRad5dry=sin(angleOfPESensitivity/2)*(dot((dRot5dry*b5),tau5)+dot(dtau5dry,(Rot5*b5)));
    ddetRad5drz=sin(angleOfPESensitivity/2)*(dot((dRot5drz*b5),tau5)+dot(dtau5drz,(Rot5*b5)));
    %partial derivative of detection radius w.r.t. the quaternions
    ddetRad5dq1=sin(angleOfPESensitivity/2)*(dot((dRot5dq1*b5),tau5)+dot(dtau5dq1,(Rot5*b5)));
    ddetRad5dq2=sin(angleOfPESensitivity/2)*(dot((dRot5dq2*b5),tau5)+dot(dtau5dq2,(Rot5*b5)));
    ddetRad5dq3=sin(angleOfPESensitivity/2)*(dot((dRot5dq3*b5),tau5)+dot(dtau5dq3,(Rot5*b5)));
    ddetRad5dq0=sin(angleOfPESensitivity/2)*(dot((dRot5dq4*b5),tau5)+dot(dtau5dq4,(Rot5*b5)));
    
    stripLeadingEdgeDistance=alpha5-stripDistances-stripThickness/2;
    stripTrailingEdgeDistance=alpha5-stripDistances+stripThickness/2;
       
    %Finally calculating the actual terms for the Jacobian
    df5drx=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5drx+(bx5*ddetRad5drx-detectionRad5*dbx5drx)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5drx-(bx5*ddetRad5drx-detectionRad5*dbx5drx)/(bx5^2))));
    df5dry=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5dry+(bx5*ddetRad5dry-detectionRad5*dbx5dry)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5dry-(bx5*ddetRad5dry-detectionRad5*dbx5dry)/(bx5^2))));
    df5drz=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5drz+(bx5*ddetRad5drz-detectionRad5*dbx5drz)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5drz-(bx5*ddetRad5drz-detectionRad5*dbx5drz)/(bx5^2))));
    
    df5dq1=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5dq1+(bx5*ddetRad5dq1-detectionRad5*dbx5dq1)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5dq1-(bx5*ddetRad5dq1-detectionRad5*dbx5dq1)/(bx5^2))));
    df5dq2=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5dq2+(bx5*ddetRad5dq2-detectionRad5*dbx5dq2)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5dq2-(bx5*ddetRad5dq2-detectionRad5*dbx5dq2)/(bx5^2))));
    df5dq3=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5dq3+(bx5*ddetRad5dq3-detectionRad5*dbx5dq3)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5dq3-(bx5*ddetRad5dq3-detectionRad5*dbx5dq3)/(bx5^2))));
    df5dq0=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad5)./bx5)- (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(dalpha5dq0+(bx5*ddetRad5dq0-detectionRad5*dbx5dq0)/(bx5^2)) - ((stripTrailingEdgeDistance>(detectionRad5)./bx5)-(stripLeadingEdgeDistance>(detectionRad5)./bx5)).*(dalpha5dq0-(bx5*ddetRad5dq0-detectionRad5*dbx5dq0)/(bx5^2))));
    
    
    H5kp1=[df5drx df5dry df5drz 0 0 0 df5dq1 df5dq2 df5dq3 df5dq0];
    S5kp1=VerticalPECovariance; %Experimentally determined
    K5kp1=P4kp1kp1*H5kp1'/(H5kp1*P4kp1kp1*H5kp1'+S5kp1); %Kalman Gain

    
    h5kp1=(maxBrightness./stripThickness).*sum((stripTrailingEdgeDistance>(-detectionRad5)./bx5).*(stripTrailingEdgeDistance+detectionRad5./bx5) - (stripLeadingEdgeDistance>(-detectionRad5)./bx5).*(stripLeadingEdgeDistance+detectionRad5./bx5) - (stripTrailingEdgeDistance>(detectionRad5)./bx5).*(stripTrailingEdgeDistance-detectionRad5./bx5) + (stripLeadingEdgeDistance>(detectionRad5)./bx5).*(stripLeadingEdgeDistance-detectionRad5./bx5));
    % The next step compares the data from the sensors with the predicted state and alters
    %the state prediction by a factor determined by the Kalman Gain
    x5kp1kp1=x4kp1kp1+K5kp1*(z5kp1-h5kp1);
    P5kp1kp1=(eye(10,10)-K5kp1*H5kp1)*P4kp1kp1;
    
    normQuat5=sqrt(sum((x5kp1kp1(7:10)).^2));
    x5kp1kp1(7:10)=x5kp1kp1(7:10)./normQuat5;
else
    x5kp1kp1=x4kp1kp1;
    P5kp1kp1=P4kp1kp1;
end




% CORRECTIVE STEP 6, 10:30 Photoelectric
if execution(6)==1
    % getting sensor data
    z6kp1 = sensorData(6,:);
    
    % Getting terms from the previous state prediction to use in calculations
    rx6=x5kp1kp1(1);
    ry6=x5kp1kp1(2);
    rz6=x5kp1kp1(3);
    q61=x5kp1kp1(7);  
    q62=x5kp1kp1(8);
    q63=x5kp1kp1(9);
    q60=x5kp1kp1(10);
    %Calculating the rotation matrix
    Rot6=[1-2*q62^2-2*q63^2 2*(q61*q62-q60*q63) 2*(q61*q63+q60*q62);...
     2*(q61*q62+q60*q63) 1-2*q61^2-2*q63^2 2*(q62*q63-q60*q61);...
     2*(q61*q63-q60*q62) 2*(q62*q63+q60*q61) 1-2*q61^2-2*q62^2;];
    sx6=(Rot6(1,:)*p6 + rx6);
    sy6=(Rot6(2,:)*p6 + ry6);
    sz6=(Rot6(3,:)*p6 + rz6);
    
    %derivatives of the rotation matrix w.r.t. each quaternion to be used to calculate the jacobian
    dRot6dq1=[0 2*q62 2*q63;...
     2*q62 -4*q61 -2*q60;...
     2*q63 2*q60 -4*q61;];
 
    dRot6dq2=[-4*q62 2*q61 2*q60;...
     2*q61 0 2*q63;...
     -2*q60 2*q63 -4*q62;];
 
    dRot6dq3=[-4*q63 -2*q60 2*q61;...
     2*q60 -4*q63 2*q62;...
     2*q61 2*q62 0;];
 
    dRot6dq0=[0 -2*q63 2*q62;...
     2*q63 0 -2*q61;...
     -2*q62 2*q61 0;];
    %intermediate terms and their partial derivative terms
    g6=(((sz6-tubeCenterToTopOfRail)*Rot6(3,:))+(sy6*Rot6(2,:)));
    dg6dq1 = (sz6-tubeCenterToTopOfRail)*dRot6dq1(3,:) + (dRot6dq1(3,:)*p6)*Rot6(3,:) + sy6*dRot6dq1(2,:) + (dRot6dq1(2,:)*p6)*Rot6(2,:);
    dg6dq2 = (sz6-tubeCenterToTopOfRail)*dRot6dq2(3,:) + (dRot6dq2(3,:)*p6)*Rot6(3,:) + sy6*dRot6dq2(2,:) + (dRot6dq2(2,:)*p6)*Rot6(2,:);
    dg6dq3 = (sz6-tubeCenterToTopOfRail)*dRot6dq3(3,:) + (dRot6dq3(3,:)*p6)*Rot6(3,:) + sy6*dRot6dq3(2,:) + (dRot6dq3(2,:)*p6)*Rot6(2,:);
    dg6dq0 = (sz6-tubeCenterToTopOfRail)*dRot6dq0(3,:) + (dRot6dq0(3,:)*p6)*Rot6(3,:) + sy6*dRot6dq0(2,:) + (dRot6dq0(2,:)*p6)*Rot6(2,:);
    %more intermediate terms and their partial derivative terms
    m6=((dot(Rot6(2,:),b6))^2+(dot(Rot6(3,:),b6))^2);
    dm6dq1 = 2*dot(b6,((dot(b6,Rot6(2,:)))*dRot6dq1(2,:) + (dot(b6,Rot6(3,:)))*dRot6dq1(3,:)));
    dm6dq2 = 2*dot(b6,((dot(b6,Rot6(2,:)))*dRot6dq2(2,:) + (dot(b6,Rot6(3,:)))*dRot6dq2(3,:)));
    dm6dq3 = 2*dot(b6,((dot(b6,Rot6(2,:)))*dRot6dq3(2,:) + (dot(b6,Rot6(3,:)))*dRot6dq3(3,:)));
    dm6dq0 = 2*dot(b6,((dot(b6,Rot6(2,:)))*dRot6dq0(2,:) + (dot(b6,Rot6(3,:)))*dRot6dq0(3,:)));
 
    %more intermediate terms and their partial derivative terms
    tau6 = ((-dot(g6,b6))+sqrt((dot(g6,b6))^2-(m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))))/m6;
    
    dtau6drx = 0;
    dtau6dry = (-dot(Rot6(2,:),b6) + ((dot(g6,b6))*(dot(b6,Rot6(2,:)))-m6*sy6)/sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2)))/m6;
    dtau6drz = (-dot(Rot6(3,:),b6) + ((dot(g6,b6))*(dot(b6,Rot6(3,:)))-m6*(sz6-tubeCenterToTopOfRail))/sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2)))/m6;
    
    dtau6dq1=m6*(-dot(b6, ((dot(p6,dRot6dq1(3,:))*Rot6(3,:) +(sz6-tubeCenterToTopOfRail)*dRot6dq1(3,:)+(dot(p6,dRot6dq1(2,:)))*Rot6(2,:)+sy6*dRot6dq1(2,:))) + (2*dot(g6,b6)*dot(dg6dq1,b6)-(2*m6*(sz6-tubeCenterToTopOfRail)*(dRot6dq1(3,:)*p6) + ((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))*dm6dq1))/(2*sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))))-((-(dot(g6,b6))+sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2)))*dm6dq1)/(m6^2);
    dtau6dq2=m6*(-dot(b6, ((dot(p6,dRot6dq2(3,:))*Rot6(3,:) +(sz6-tubeCenterToTopOfRail)*dRot6dq2(3,:)+(dot(p6,dRot6dq2(2,:)))*Rot6(2,:)+sy6*dRot6dq2(2,:))) + (2*dot(g6,b6)*dot(dg6dq2,b6)-(2*m6*(sz6-tubeCenterToTopOfRail)*(dRot6dq2(3,:)*p6) + ((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))*dm6dq2))/(2*sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))))-((-(dot(g6,b6))+sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2)))*dm6dq2)/(m6^2);
    dtau6dq3=m6*(-dot(b6, ((dot(p6,dRot6dq3(3,:))*Rot6(3,:) +(sz6-tubeCenterToTopOfRail)*dRot6dq3(3,:)+(dot(p6,dRot6dq3(2,:)))*Rot6(2,:)+sy6*dRot6dq3(2,:))) + (2*dot(g6,b6)*dot(dg6dq3,b6)-(2*m6*(sz6-tubeCenterToTopOfRail)*(dRot6dq3(3,:)*p6) + ((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))*dm6dq3))/(2*sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))))-((-(dot(g6,b6))+sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2)))*dm6dq3)/(m6^2);
    dtau6dq0=m6*(-dot(b6, ((dot(p6,dRot6dq0(3,:))*Rot6(3,:) +(sz6-tubeCenterToTopOfRail)*dRot6dq0(3,:)+(dot(p6,dRot6dq0(2,:)))*Rot6(2,:)+sy6*dRot6dq0(2,:))) + (2*dot(g6,b6)*dot(dg6dq0,b6)-(2*m6*(sz6-tubeCenterToTopOfRail)*(dRot6dq0(3,:)*p6) + ((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))*dm6dq0))/(2*sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2))))-((-(dot(g6,b6))+sqrt((dot(g6,b6))^2-m6*((sz6-tubeCenterToTopOfRail)^2-tubeRadius^2+sy6^2)))*dm6dq0)/(m6^2);
   
    %more intermediate terms and their partial derivative terms
    alpha6=sx6+dot(Rot6(1,:),b6*tau6);
    
    dalpha6drx = 1+dot(Rot6(1,:),b6*dtau6drx);
    dalpha6dry = dot(Rot6(1,:),b6*dtau6dry);
    dalpha6drz = dot(Rot6(1,:),b6*dtau6dry);
   
    dalpha6dq1 = dRot6dq1(1,:)*p6 + dot(Rot6(1,:),b6*dtau6dq1) + dot(dRot6dq1(1,:),b6*tau6);
    dalpha6dq2 = dRot6dq2(1,:)*p6 + dot(Rot6(1,:),b6*dtau6dq2) + dot(dRot6dq2(1,:),b6*tau6);
    dalpha6dq3 = dRot6dq3(1,:)*p6 + dot(Rot6(1,:),b6*dtau6dq3) + dot(dRot6dq3(1,:),b6*tau6);
    dalpha6dq0 = dRot6dq0(1,:)*p6 + dot(Rot6(1,:),b6*dtau6dq0) + dot(dRot6dq0(1,:),b6*tau6);
    
    dist6=dot((Rot6*b6),tau6);
    %more intermediate terms and their partial derivative terms
    bx6=dot(Rot6(1,:),(b6/norm(b6)));
    
    dbx6drx=0;
    dbx6dry=0;
    dbx6drz=0;
    
    dbx6dq1=dot(dRot6dq1(1,:),(b6/norm(b6)));
    dbx6dq2=dot(dRot6dq2(1,:),(b6/norm(b6)));
    dbx6dq3=dot(dRot6dq3(1,:),(b6/norm(b6)));
    dbx6dq0=dot(dRot6dq0(1,:),(b6/norm(b6)));
    
    detectionRad6=dist6*sin(angleOfPESensitivity/2);
    %partial derivative of detection radius w.r.t. the position
    ddetRad6drx=sin(angleOfPESensitivity/2)*(dot((dRot6drx*b6),tau6)+dot(dtau6drx,(Rot6*b6)));
    ddetRad6dry=sin(angleOfPESensitivity/2)*(dot((dRot6dry*b6),tau6)+dot(dtau6dry,(Rot6*b6)));
    ddetRad6drz=sin(angleOfPESensitivity/2)*(dot((dRot6drz*b6),tau6)+dot(dtau6drz,(Rot6*b6)));
    %partial derivative of detection radius w.r.t. the quaternions
    ddetRad6dq1=sin(angleOfPESensitivity/2)*(dot((dRot6dq1*b6),tau6)+dot(dtau6dq1,(Rot6*b6)));
    ddetRad6dq2=sin(angleOfPESensitivity/2)*(dot((dRot6dq2*b6),tau6)+dot(dtau6dq2,(Rot6*b6)));
    ddetRad6dq3=sin(angleOfPESensitivity/2)*(dot((dRot6dq3*b6),tau6)+dot(dtau6dq3,(Rot6*b6)));
    ddetRad6dq0=sin(angleOfPESensitivity/2)*(dot((dRot6dq4*b6),tau6)+dot(dtau6dq4,(Rot6*b6)));
    
    stripLeadingEdgeDistance=alpha6-stripDistances-stripThickness/2;
    stripTrailingEdgeDistance=alpha6-stripDistances+stripThickness/2;
       
    %Finally calculating the actual terms for the Jacobian
    df6drx=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6drx+(bx6*ddetRad6drx-detectionRad6*dbx6drx)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6drx-(bx6*ddetRad6drx-detectionRad6*dbx6drx)/(bx6^2))));
    df6dry=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6dry+(bx6*ddetRad6dry-detectionRad6*dbx6dry)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6dry-(bx6*ddetRad6dry-detectionRad6*dbx6dry)/(bx6^2))));
    df6drz=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6drz+(bx6*ddetRad6drz-detectionRad6*dbx6drz)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6drz-(bx6*ddetRad6drz-detectionRad6*dbx6drz)/(bx6^2))));
    
    df6dq1=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6dq1+(bx6*ddetRad6dq1-detectionRad6*dbx6dq1)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6dq1-(bx6*ddetRad6dq1-detectionRad6*dbx6dq1)/(bx6^2))));
    df6dq2=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6dq2+(bx6*ddetRad6dq2-detectionRad6*dbx6dq2)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6dq2-(bx6*ddetRad6dq2-detectionRad6*dbx6dq2)/(bx6^2))));
    df6dq3=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6dq3+(bx6*ddetRad6dq3-detectionRad6*dbx6dq3)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6dq3-(bx6*ddetRad6dq3-detectionRad6*dbx6dq3)/(bx6^2))));
    df6dq0=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad6)./bx6)- (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(dalpha6dq0+(bx6*ddetRad6dq0-detectionRad6*dbx6dq0)/(bx6^2)) - ((stripTrailingEdgeDistance>(detectionRad6)./bx6)-(stripLeadingEdgeDistance>(detectionRad6)./bx6)).*(dalpha6dq0-(bx6*ddetRad6dq0-detectionRad6*dbx6dq0)/(bx6^2))));
    
    
    H6kp1=[df6drx df6dry df6drz 0 0 0 df6dq1 df6dq2 df6dq3 df6dq0];
    S6kp1=at1030PECovariance; %Experimentally determined
    K6kp1=P5kp1kp1*H6kp1'/(H6kp1*P5kp1kp1*H6kp1'+S6kp1); %Kalman Gain

    
    h6kp1=(maxBrightness./stripThickness).*sum((stripTrailingEdgeDistance>(-detectionRad6)./bx6).*(stripTrailingEdgeDistance+detectionRad6./bx6) - (stripLeadingEdgeDistance>(-detectionRad6)./bx6).*(stripLeadingEdgeDistance+detectionRad6./bx6) - (stripTrailingEdgeDistance>(detectionRad6)./bx6).*(stripTrailingEdgeDistance-detectionRad6./bx6) + (stripLeadingEdgeDistance>(detectionRad6)./bx6).*(stripLeadingEdgeDistance-detectionRad6./bx6));
    
    % The next step compares the data from the sensors with the predicted state and alters
    %the state prediction by a factor determined by the Kalman Gain        
    x6kp1kp1=x5kp1kp1+K6kp1*(z6kp1-h6kp1);
    P6kp1kp1=(eye(10,10)-K6kp1*H6kp1)*P5kp1kp1;
    
    normQuat6=sqrt(sum((x6kp1kp1(7:10)).^2));
    x6kp1kp1(7:10)=x6kp1kp1(7:10)./normQuat6;
else
    x6kp1kp1=x5kp1kp1;
    P6kp1kp1=P5kp1kp1;
end



% CORRECTIVE STEP 7, 1:30 Photoelectric
if execution(7)==1
    % getting sensor data
    z7kp1 = sensorData(7,:);
    
    % Getting terms from the previous state prediction to use in calculations
    rx7=x6kp1kp1(1);
    ry7=x6kp1kp1(2);
    rz7=x6kp1kp1(3);
    q71=x6kp1kp1(7);  
    q72=x6kp1kp1(8);
    q73=x6kp1kp1(9);
    q70=x6kp1kp1(10);
    %Calculating the rotation matrix
    Rot7=[1-2*q72^2-2*q73^2 2*(q71*q72-q70*q73) 2*(q71*q73+q70*q72);...
     2*(q71*q72+q70*q73) 1-2*q71^2-2*q73^2 2*(q72*q73-q70*q71);...
     2*(q71*q73-q70*q72) 2*(q72*q73+q70*q71) 1-2*q71^2-2*q72^2;];
    sx7=(Rot7(1,:)*p7 + rx7);
    sy7=(Rot7(2,:)*p7 + ry7);
    sz7=(Rot7(3,:)*p7 + rz7);
    
    %derivatives of the rotation matrix w.r.t. each quaternion to be used to calculate the jacobian
    dRot7dq1=[0 2*q72 2*q73;...
     2*q72 -4*q71 -2*q70;...
     2*q73 2*q70 -4*q71;];
 
    dRot7dq2=[-4*q72 2*q71 2*q70;...
     2*q71 0 2*q73;...
     -2*q70 2*q73 -4*q72;];
 
    dRot7dq3=[-4*q73 -2*q70 2*q71;...
     2*q70 -4*q73 2*q72;...
     2*q71 2*q72 0;];
 
    dRot7dq0=[0 -2*q73 2*q72;...
     2*q73 0 -2*q71;...
     -2*q72 2*q71 0;];
    %intermediate terms and their partial derivative terms
    g7=(((sz7-tubeCenterToTopOfRail)*Rot7(3,:))+(sy7*Rot7(2,:)));
    dg7dq1 = (sz7-tubeCenterToTopOfRail)*dRot7dq1(3,:) + (dRot7dq1(3,:)*p7)*Rot7(3,:) + sy7*dRot7dq1(2,:) + (dRot7dq1(2,:)*p7)*Rot7(2,:);
    dg7dq2 = (sz7-tubeCenterToTopOfRail)*dRot7dq2(3,:) + (dRot7dq2(3,:)*p7)*Rot7(3,:) + sy7*dRot7dq2(2,:) + (dRot7dq2(2,:)*p7)*Rot7(2,:);
    dg7dq3 = (sz7-tubeCenterToTopOfRail)*dRot7dq3(3,:) + (dRot7dq3(3,:)*p7)*Rot7(3,:) + sy7*dRot7dq3(2,:) + (dRot7dq3(2,:)*p7)*Rot7(2,:);
    dg7dq0 = (sz7-tubeCenterToTopOfRail)*dRot7dq0(3,:) + (dRot7dq0(3,:)*p7)*Rot7(3,:) + sy7*dRot7dq0(2,:) + (dRot7dq0(2,:)*p7)*Rot7(2,:);
    %more intermediate terms and their partial derivative terms
    m7=((dot(Rot7(2,:),b7))^2+(dot(Rot7(3,:),b7))^2);
    dm7dq1 = 2*dot(b7,((dot(b7,Rot7(2,:)))*dRot7dq1(2,:) + (dot(b7,Rot7(3,:)))*dRot7dq1(3,:)));
    dm7dq2 = 2*dot(b7,((dot(b7,Rot7(2,:)))*dRot7dq2(2,:) + (dot(b7,Rot7(3,:)))*dRot7dq2(3,:)));
    dm7dq3 = 2*dot(b7,((dot(b7,Rot7(2,:)))*dRot7dq3(2,:) + (dot(b7,Rot7(3,:)))*dRot7dq3(3,:)));
    dm7dq0 = 2*dot(b7,((dot(b7,Rot7(2,:)))*dRot7dq0(2,:) + (dot(b7,Rot7(3,:)))*dRot7dq0(3,:)));
    %more intermediate terms and their partial derivative terms
    tau7 = ((-dot(g7,b7))+sqrt((dot(g7,b7))^2-(m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))))/m7;
    
    dtau7drx = 0;
    dtau7dry = (-dot(Rot7(2,:),b7) + ((dot(g7,b7))*(dot(b7,Rot7(2,:)))-m7*sy7)/sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2)))/m7;
    dtau7drz = (-dot(Rot7(3,:),b7) + ((dot(g7,b7))*(dot(b7,Rot7(3,:)))-m7*(sz7-tubeCenterToTopOfRail))/sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2)))/m7;
    
    dtau7dq1=m7*(-dot(b7, ((dot(p7,dRot7dq1(3,:))*Rot7(3,:) +(sz7-tubeCenterToTopOfRail)*dRot7dq1(3,:)+(dot(p7,dRot7dq1(2,:)))*Rot7(2,:)+sy7*dRot7dq1(2,:))) + (2*dot(g7,b7)*dot(dg7dq1,b7)-(2*m7*(sz7-tubeCenterToTopOfRail)*(dRot7dq1(3,:)*p7) + ((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))*dm7dq1))/(2*sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))))-((-(dot(g7,b7))+sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2)))*dm7dq1)/(m7^2);
    dtau7dq2=m7*(-dot(b7, ((dot(p7,dRot7dq2(3,:))*Rot7(3,:) +(sz7-tubeCenterToTopOfRail)*dRot7dq2(3,:)+(dot(p7,dRot7dq2(2,:)))*Rot7(2,:)+sy7*dRot7dq2(2,:))) + (2*dot(g7,b7)*dot(dg7dq2,b7)-(2*m7*(sz7-tubeCenterToTopOfRail)*(dRot7dq2(3,:)*p7) + ((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))*dm7dq2))/(2*sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))))-((-(dot(g7,b7))+sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2)))*dm7dq2)/(m7^2);
    dtau7dq3=m7*(-dot(b7, ((dot(p7,dRot7dq3(3,:))*Rot7(3,:) +(sz7-tubeCenterToTopOfRail)*dRot7dq3(3,:)+(dot(p7,dRot7dq3(2,:)))*Rot7(2,:)+sy7*dRot7dq3(2,:))) + (2*dot(g7,b7)*dot(dg7dq3,b7)-(2*m7*(sz7-tubeCenterToTopOfRail)*(dRot7dq3(3,:)*p7) + ((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))*dm7dq3))/(2*sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))))-((-(dot(g7,b7))+sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2)))*dm7dq3)/(m7^2);
    dtau7dq0=m7*(-dot(b7, ((dot(p7,dRot7dq0(3,:))*Rot7(3,:) +(sz7-tubeCenterToTopOfRail)*dRot7dq0(3,:)+(dot(p7,dRot7dq0(2,:)))*Rot7(2,:)+sy7*dRot7dq0(2,:))) + (2*dot(g7,b7)*dot(dg7dq0,b7)-(2*m7*(sz7-tubeCenterToTopOfRail)*(dRot7dq0(3,:)*p7) + ((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))*dm7dq0))/(2*sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2))))-((-(dot(g7,b7))+sqrt((dot(g7,b7))^2-m7*((sz7-tubeCenterToTopOfRail)^2-tubeRadius^2+sy7^2)))*dm7dq0)/(m7^2);
   
    %more intermediate terms and their partial derivative terms
    alpha7=sx7+dot(Rot7(1,:),b7*tau7);
    
    dalpha7drx = 1+dot(Rot7(1,:),b7*dtau7drx);
    dalpha7dry = dot(Rot7(1,:),b7*dtau7dry);
    dalpha7drz = dot(Rot7(1,:),b7*dtau7dry);
   
    dalpha7dq1 = dRot7dq1(1,:)*p7 + dot(Rot7(1,:),b7*dtau7dq1) + dot(dRot7dq1(1,:),b7*tau7);
    dalpha7dq2 = dRot7dq2(1,:)*p7 + dot(Rot7(1,:),b7*dtau7dq2) + dot(dRot7dq2(1,:),b7*tau7);
    dalpha7dq3 = dRot7dq3(1,:)*p7 + dot(Rot7(1,:),b7*dtau7dq3) + dot(dRot7dq3(1,:),b7*tau7);
    dalpha7dq0 = dRot7dq0(1,:)*p7 + dot(Rot7(1,:),b7*dtau7dq0) + dot(dRot7dq0(1,:),b7*tau7);
    
    dist7=dot((Rot7*b7),tau7);
    %more intermediate terms and their partial derivative terms
    bx7=dot(Rot7(1,:),(b7/norm(b7)));
    
    dbx7drx=0;
    dbx7dry=0;
    dbx7drz=0;
    
    dbx7dq1=dot(dRot7dq1(1,:),(b7/norm(b7)));
    dbx7dq2=dot(dRot7dq2(1,:),(b7/norm(b7)));
    dbx7dq3=dot(dRot7dq3(1,:),(b7/norm(b7)));
    dbx7dq0=dot(dRot7dq0(1,:),(b7/norm(b7)));
    
    detectionRad7=dist7*sin(angleOfPESensitivity/2);
    %partial derivative of detection radius w.r.t. the position
    ddetRad7drx=sin(angleOfPESensitivity/2)*(dot((dRot7drx*b7),tau7)+dot(dtau7drx,(Rot7*b7)));
    ddetRad7dry=sin(angleOfPESensitivity/2)*(dot((dRot7dry*b7),tau7)+dot(dtau7dry,(Rot7*b7)));
    ddetRad7drz=sin(angleOfPESensitivity/2)*(dot((dRot7drz*b7),tau7)+dot(dtau7drz,(Rot7*b7)));
    %partial derivative of detection radius w.r.t. the quaternions
    ddetRad7dq1=sin(angleOfPESensitivity/2)*(dot((dRot7dq1*b7),tau7)+dot(dtau7dq1,(Rot7*b7)));
    ddetRad7dq2=sin(angleOfPESensitivity/2)*(dot((dRot7dq2*b7),tau7)+dot(dtau7dq2,(Rot7*b7)));
    ddetRad7dq3=sin(angleOfPESensitivity/2)*(dot((dRot7dq3*b7),tau7)+dot(dtau7dq3,(Rot7*b7)));
    ddetRad7dq0=sin(angleOfPESensitivity/2)*(dot((dRot7dq4*b7),tau7)+dot(dtau7dq4,(Rot7*b7)));
    
    stripLeadingEdgeDistance=alpha7-stripDistances-stripThickness/2;
    stripTrailingEdgeDistance=alpha7-stripDistances+stripThickness/2;
       
    %Finally calculating the actual terms for the Jacobian
    df7drx=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7drx+(bx7*ddetRad7drx-detectionRad7*dbx7drx)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7drx-(bx7*ddetRad7drx-detectionRad7*dbx7drx)/(bx7^2))));
    df7dry=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7dry+(bx7*ddetRad7dry-detectionRad7*dbx7dry)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7dry-(bx7*ddetRad7dry-detectionRad7*dbx7dry)/(bx7^2))));
    df7drz=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7drz+(bx7*ddetRad7drz-detectionRad7*dbx7drz)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7drz-(bx7*ddetRad7drz-detectionRad7*dbx7drz)/(bx7^2))));
    
    df7dq1=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7dq1+(bx7*ddetRad7dq1-detectionRad7*dbx7dq1)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7dq1-(bx7*ddetRad7dq1-detectionRad7*dbx7dq1)/(bx7^2))));
    df7dq2=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7dq2+(bx7*ddetRad7dq2-detectionRad7*dbx7dq2)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7dq2-(bx7*ddetRad7dq2-detectionRad7*dbx7dq2)/(bx7^2))));
    df7dq3=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7dq3+(bx7*ddetRad7dq3-detectionRad7*dbx7dq3)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7dq3-(bx7*ddetRad7dq3-detectionRad7*dbx7dq3)/(bx7^2))));
    df7dq0=(maxBrightness./stripThickness).*sum(((stripTrailingEdgeDistance>(-detectionRad7)./bx7)- (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(dalpha7dq0+(bx7*ddetRad7dq0-detectionRad7*dbx7dq0)/(bx7^2)) - ((stripTrailingEdgeDistance>(detectionRad7)./bx7)-(stripLeadingEdgeDistance>(detectionRad7)./bx7)).*(dalpha7dq0-(bx7*ddetRad7dq0-detectionRad7*dbx7dq0)/(bx7^2))));
    
    
    H7kp1=[df7drx df7dry df7drz 0 0 0 df7dq1 df7dq2 df7dq3 df7dq0];
    S7kp1=at130PECovariance; %Experimentally determined
    K7kp1=P6kp1kp1*H7kp1'/(H7kp1*P6kp1kp1*H7kp1'+S7kp1); %Kalman Gain

    
    h7kp1=(maxBrightness./stripThickness).*sum((stripTrailingEdgeDistance>(-detectionRad7)./bx7).*(stripTrailingEdgeDistance+detectionRad7./bx7) - (stripLeadingEdgeDistance>(-detectionRad7)./bx7).*(stripLeadingEdgeDistance+detectionRad7./bx7) - (stripTrailingEdgeDistance>(detectionRad7)./bx7).*(stripTrailingEdgeDistance-detectionRad7./bx7) + (stripLeadingEdgeDistance>(detectionRad7)./bx7).*(stripLeadingEdgeDistance-detectionRad7./bx7));
    
    % The next step compares the data from the sensors with the predicted state and alters
    %the state prediction by a factor determined by the Kalman Gain 
    x7kp1kp1=x6kp1kp1+K7kp1*(z7kp1-h7kp1);
    P7kp1kp1=(eye(10,10)-K7kp1*H7kp1)*P6kp1kp1;
    
    normQuat7=sqrt(sum((x7kp1kp1(7:10)).^2));
    x7kp1kp1(7:10)=x7kp1kp1(7:10)./normQuat7;
else
    x7kp1kp1=x6kp1kp1;
    P7kp1kp1=P6kp1kp1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

predictedState = x7kp1kp1;
predictedCovariance= P7kp1kp1;
    
    
    
    
    
    
    
    

