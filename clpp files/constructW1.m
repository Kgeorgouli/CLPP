function [ res ] = constructW1( X, tw,stdAcceptancee,smoothSizee )
%CONSTRUCTW1 Summary of this function goes here
%  X data
%  tw number of neighbours who define a fragment (repetition)
%   Detailed explanation goes here
%% Construct similarity graph 
%     disp('Constructing repetition neighbourhood (temporal) graph...');
%     Konstantia
     stdAcceptance=stdAcceptancee; 
     smoothSize=smoothSizee;
     mintab=[];
    
    [P,d]=size(X); 
    DX = KG_clppL2(X',X',1);  %it is the same with TLE, I changed only the name
    
    if ~exist('DTW', 'var')
        [DTW1,DTW2] = ML_tleGetMatrixOfDTW(X,tw,P);   %calculates matrix of DTW distances
        %DTW1=sparse(DTW1);DTW2=sparse(DTW2);
        
    else
          DTW1=sparse(DTW.DTW1);
          DTW2=sparse(DTW.DTW2);
    end

        %%%%%%%%%%%%
        DTW.DTW1 = DTW1;
        DTW.DTW2 = DTW2;
    %%%%%%%%%%%
    
    
    [G]=ML_tleGetDTW_CTN(X,DX,tw,stdAcceptance,DTW1,DTW2,P,smoothSize,mintab); %constructs repetition neighbourhood graph
    
    res.G=G;
    res.DX=DX;

end

