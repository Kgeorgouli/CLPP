function t_point = out_of_sample(point, mapping)
%OUT_OF_SAMPLE Performs out-of-sample extension of new datapoints
%
%   t_points = out_of_sample(points, mapping)
%
% Performs out-of-sample extension of the new datapoints in points. The 
% information necessary for the out-of-sample extension is contained in 
% mapping (this struct can be obtained from COMPUTE_MAPPING).
% The function returns the coordinates of the transformed points in t_points.

%      Reference:
%                Georgouli, K., Del Rincon, J.M. and Koidis, A., 2017. Continuous statistical
%                modelling for rapid detection of adulteration of extra virgin olive oil using
%                mid infrared and Raman spectroscopic data. Food chemistry, 217, pp.735-742.
 
% (C) Konstantia Georgouli, Queen's University of Belfast

%

       [nSmp,nFea] = size(point);
       
        point = (point - repmat(mapping.mean,nSmp,1));          
        t_point= point * mapping.M;

end

