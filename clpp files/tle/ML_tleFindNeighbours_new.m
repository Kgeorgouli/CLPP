function [mintab] = ML_tleFindNeighbours(DTW,tvwin,P,stdAcceptance,smoothSize)
% ML_tleFindNeighbours - finds repetition neighbours for each point
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    mintab=cell([P,1]); 
    for i=1:P
        [iS,iE] = ML_tleGetWindowIndex(i,P,tvwin);
        [mintab{i}]=ML_tleAdoptedToleranceExtremeExHum(DTW(i,:),iS,iE,P,stdAcceptance);
    end
end