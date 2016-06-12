classdef globalData

    properties
        
        %%% TIMING %%%
        timestep = .001 %s
        kalmanTimestep = .01 %s
        runtime = 5 %s
        numSteps

        %%% SIMULATION OPTIONS %%%
        randomNoise=false;
        noiseModifier=0.00000001;
            
        %%% TUBE SPECS %%%
        pusherForce = 17640 % newtons
        pusherDistance = 273 % m
        
        %%% NATURAL CONSTANTS %%%
        gravity = 9.81 % m/s^2
        TEMPERATURE = 350 % K
        TUBE_PRESSURE = 1000 % Pa
        GAS_CONSTANT = 287.058 % Pa * m^3 / K
        airDensity
        
        %%% Covariances 'n stuff %%%
        laserCovariance = .0000000001
        IMUCovariance = .0000000001
        
    end
    methods
        function globals = globalData()
            globals.numSteps = globals.runtime / globals.timestep;
            globals.airDensity = globals.TUBE_PRESSURE...
                /(globals.GAS_CONSTANT*globals.TEMPERATURE);
        end        
    end
end
