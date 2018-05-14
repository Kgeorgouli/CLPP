function [mintab,maxtab]=ML_tleAdoptedToleranceExtremeExHum(v,iS,iE,P,stdAcceptance)
% ML_tleAdoptedToleranceExtremeExHum - finds minima in v
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

vBase = v(:);
x = (1:length(vBase))';

P=[0,P];
minTabAll=[];
maxTabAll=[];
tvwin=(iE-iS+1)/2;
tvwin=5;

minFind=0;
while(~minFind)
        
    maxTabT = [];
    minTabT = [];
    for m=2:length(P)

        maxtab = [];
        mintab = [];
        v=vBase(P(m-1)+1:P(m));  %print for nice function
%     clf(figure(1));
%     plot(v);
    
        sd = stdAcceptance*std(v(v~=inf));

        mn = Inf; mx = -Inf;
        mnpos = NaN; mxpos = NaN;
        currminpos = 1;
        lookformax = v(1)<v(2);

        for i=1:length(v)
          this = v(i);
          if this > mx, mx = this; mxpos = x(i); end
          if this < mn, mn = this; mnpos = x(i); end
          if lookformax
            if (isinf(this)==0 && (this < mx-sd || (i==length(v) && mxpos<length(v) && mxpos+tvwin<=length(v) && v(mxpos+1)<v(mxpos) && size(maxtab,1)>0 && ismember(mxpos,maxtab(:,1))==0))) %prawe zbocze maximum
                ismax=1;            
                if(ismax==1)
                    ismax=0;
                    if(size(mintab,1)==0)
                        lastmin=1;
                    else 
                        lastmin = mintab(end,1);
                    end    
                    for j=mxpos:-1:lastmin  
                        if(v(j)<mx-sd || (size(mintab,1)==0 && v(j)<mx && j>tvwin))
                            ismax=1;
                            break;
                        end
                    end
                    if(ismax==1)
                        maxtab = [maxtab ; mxpos mx];
                        mn = this; mnpos = x(i);
                    end
                end
                lookformax = 0;
            end  
          else
            if (isinf(this)==0 && (this > mn+sd || (i==length(v) && mnpos<length(v) && mnpos+tvwin<=length(v) && v(mnpos+1)>v(mnpos) && size(mintab,1)>0 && ismember(mnpos,mintab(:,1))==0)))  %prawe zbocze minimum
                ismin=0;
                if(size(maxtab,1)==0)
                    lastmax=1;
                else 
                    lastmax = maxtab(end,1);
                end
                for j=mnpos:-1:lastmax 
                    if(v(j)>mn+sd || (size(maxtab,1)==0 && v(j)>mn && j>tvwin))
                        ismin=1;
                        break;
                    end
                end
                if(ismin==1)
                    if(P(m-1)+mnpos < iS || P(m-1)+mnpos >iE)
                        mintab = [mintab ; mnpos mn];
                    else
                        currminpos = size(mintab,1)+1;
                    end
                    mx = this; mxpos = x(i);
                end
                lookformax = 1;
            else
                if(mnpos >= iS && mnpos <=iE && currminpos==1)
                    currminpos = size(mintab,1)+1;
                end
            end
          end

        end

        if(isempty(mintab)==0)
            xmax = sort(maxtab(:,2),'descend');
            xmin = sort(mintab(:,2),'ascend'); 

            mmm=mean([xmax(1),xmin(1)]);

            mintab(:,1) = mintab(:,1) + repmat(P(m-1),size(mintab,1),1);
            maxtab(:,1) = maxtab(:,1) + repmat(P(m-1),size(maxtab,1),1);    

            minTabAll=[minTabAll;mintab];
            maxTabAll=[maxTabAll;maxtab];
            minTabT=[minTabT;mintab];
            maxTabT=[maxTabT;maxtab];
        end
    end
    if(isempty(minTabT)==1)
        stdAcceptance=stdAcceptance-0.1;
        %warning('Standard deviation was reduced!!!');  %%%Konstantia
        if(stdAcceptance<0.1)
            tvwin=tvwin-1;
            warning('Standard deviation too small!!!');
            if(tvwin<=0)
                error('Standard deviation and tvwin too small!!!');
            end
        end
    else
        minFind=1;
    end
end
mintab=minTabAll;
maxtab=maxTabAll;

end