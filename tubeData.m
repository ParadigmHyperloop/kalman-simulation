classdef tubeData

    properties
        length = 1600 % m
        pressure = 1 % kPa
        radius = .889 % m
        stripDistances=[100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 ...
            4280 4280+(1/3) 4280+(2/3) 4281 4281+(1/3) 4281+(2/3) 4282 4281+(1/3) 4281+(2/3) 4282 4282+(1/3) 4282+(2/3) 4283 4283+(1/3) 4283+(2/3) 4284 4284+(1/3) 4284+(2/3) 4285 4285+(1/3) ...
            4300 4400 4500 4600 4700 ...
            4780 4780+(1/3) 4780+(2/3) 4781 4781+(1/3) 4781+(2/3) 4782 4782+(1/3) 4782+(2/3) 4783 ...
            4800 4900 5000 5100 5200]*12*0.0254; % m
        stripWidth=2*0.0254;% m
  
        tubeCenterHeight = .72 %m
        
        % rail height is total height minus height of the track
        railHeight = .127 - .0127 %m
        railTopThickness = .0104 %m
        topRailWidth = .127 %m
        midRailWidth = .00795 %m 
        maxBrightness = 1 %who the fuck knows, maybe candelas?
        tubeCenterToTopOfRail = 0.4 %m
        angleOfPESensitivity = 10*pi/180 %rad
        
    end
end
