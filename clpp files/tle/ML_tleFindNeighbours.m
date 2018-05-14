function [mintab,DTW] = ML_tleFindNeighbours(DTW1,DTW2,tvwin,P,stdAcceptance,smoothSize)
% ML_tleFindNeighbours - finds repetition neighbours for each point
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk
if size(DTW2,1)==0
    DTW=DTW1;
else
    [DTW] = ML_tleGetMatrixOfDistMoments(DTW1,DTW2,tvwin,P,smoothSize);
end
    mintab=cell([P,1]); 
    for i=1:P
        [iS,iE] = ML_tleGetWindowIndex(i,P,tvwin);
        [mintab{i}]=ML_tleAdoptedToleranceExtremeExHum(DTW(i,:),iS,iE,P,stdAcceptance);
    end
end