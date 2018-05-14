function [mappedX, mapping] = CLPP( data,corresp,no_dims,nc,ns,smoothSize,kratio,beta )
%CLPP Perform the CLPP algorithm
%
%   [mappedX, mapping] = CLPP( data,clpp_labels,no_dims,nc,ns,smoothSize,kratio,beta 
%         Input:
%         
%             data - Data matrix
%             corresp - vector indicating the identity of series 
%             no_dims - number of dimensions of the latent space
%             nc - number of continuous neighbours
%             ns - number of similar neighbours who define a fragment (repetition)
%             smoothSize - smoothing factor
%             beta -balance between the graphs
% 
%         Output:
%             mappedX - the coordinates of the low-dimensional data.            
%             mapping  - information on the mapping.
%
%         Reference:
%                Georgouli, K., Del Rincon, J.M. and Koidis, A., 2017. Continuous statistical
%                modelling for rapid detection of adulteration of extra virgin olive oil using
%                mid infrared and Raman spectroscopic data. Food chemistry, 217, pp.735-742.
 
% (C) Konstantia Georgouli, Queen's University of Belfast

 %% Execution of CLPP
               
       
        stdAcceptance=1; %1
        isCyclic =0; %0 if no cilci activity (end point conects with first point)
        

         
        
         
         % Y input data MxN WITH m NUMBER OF SAMPELS AND n NUMBER OF FEATURES OF
         % EACH CSAMPLE
        options = [];
        options.Metric = 'Euclidean';
        options.NeighborMode = 'KNN';
        options.PCARatio =kratio; %kratio 1
        
        
%         options.NeighborMode = 'Supervised';
%       options.gnd = classes_training;
%       options.bLDA = 1;
            
        
        
        
        resW1 = constructW1(data,ns,stdAcceptance,smoothSize);  %Affinity matrix      % stdAcceptance,smoothSize take values in the function
        W=resW1.G;
        sW = constructW2(data,resW1.DX,nc,corresp,isCyclic);           % isCyclic takes values in the function

        options.ReducedDim=no_dims;
%       %create clpp space
        [mappedX, mapping]=LPP(W,sW,options,data,beta); 
       

end

