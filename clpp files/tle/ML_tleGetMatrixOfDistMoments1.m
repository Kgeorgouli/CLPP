function [DF,DB] = ML_tleGetMatrixOfDistMoments1(DD1,DD2)
% ML_tleGetMatrixOfDistMoments - combines matrices of DTW distances
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    if(size(DD1,1)~=size(DD1,2))
        DF=ML_L2(DD1',DD1',1);
        %DF=sparse(DF);
    else
        DF=DD1;
    end
    if(size(DD2,1)~=size(DD2,2))
        DB=ML_L2(DD2',DD2',1);  
        %DB=sparse(DB);
    else
        DB=DD2;
    end
    
end