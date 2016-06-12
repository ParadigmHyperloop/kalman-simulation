classdef podData
    properties
        mass
        COM
        tensor
        length
        width
        height
        collisionPoints
        leftRailWheels
        rightRailWheels
        airskate
        skateSegmentLength
        eBrakeForce
        dragCoef
    end
    methods
        function pod = podData()
            pod.mass = 750;
            pod.COM = [-1.5904;0;0]; %%%%%% THIS IS WRONG, BUT WE DON'T HAVE A CONTROL ALGORITHM YET
            pod.length = 6;
            pod.width = 1;
            pod.height = 1;
            pod.tensor=[   1.0/12*pod.mass*(pod.height^2 + pod.width^2) 0 0; ...
                            0 1.0/12*pod.mass*(pod.height^2 + pod.length^2) 0; ...
                            0 0 1.0/12*pod.mass*(pod.length^2 + pod.width^2)];
                        
            pod.collisionPoints =[[1.04305	-0.08234	-0.06705]'...
                                  [1.04305	 0.08276	-0.06705]'...
                                  [1.00203	-0.09498	-0.09674]'...
                                  [1.00203	 0.09498	-0.09674]'...
                                  [-4.31659	-0.0951     -0.09669]'...
                                  [-4.31659	 0.0951     -0.09669]'...
                                  [0.25835	-0.4064 	-0.11263]'...
                                  [0.25835	 0.4064 	-0.11263]'...
                                  [-3.45005	-0.4064 	-0.11263]'...
                                  [-3.45005	 0.4064 	-0.11263]'...
                                  [0.25835	 0.0762 	-0.11263]'...
                                  [0.25835	-0.0762 	-0.11263]'...
                                  [-3.45005	 0.0762 	-0.11263]'...
                                  [-3.45005 -0.0762 	-0.11263]'];

            pod.skateSegmentLength = 1.14;
            pod.dragCoef = 2.7;
                    
            % GENERATE SKATE POINTS

            
            
            skateCoordinates = [
                [-3.3917382	-0.3556	-0.113919]'...
                [-2.16408	-0.3556	-0.113919]'...
                [-0.9364218	-0.3556	-0.113919]'...
                [-3.3917382	 0.127	-0.113919]'...
                [-2.16408	 0.127	-0.113919]'...
                [-0.9364218	 0.127	-0.113919]'];

            pod.airskate = zeros(6,4,3); 
    
            for i = 1:6
               % each segment has 10x46 holes 1" apart - simulate as fewer
               % holes for increases speed
               segment = zeros(4,3);
               for x = 0:45:45
                   for y = 0:9:9
                       segment(x/45*2 + y/9+1,:) = skateCoordinates(:,i) + .0254*[x;y;0];
%                        display(segment(x*10 + y+1,:))
                   end
                   
               end
               pod.airskate(i,:,:) = segment;
            end

            pod.eBrakeForce = 0.8*9.8*pod.mass;
            
        end
        
    end
end

    