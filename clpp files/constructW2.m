function [ sG ] = constructW2( X,DX, tk, Pindex ,isCyclicc)
%CONSTRUCTW2 Summary of this function goes here
% X data
% DX Euclidean distance matrix
% tk number of temporal neighbours
%   Detailed explanation goes here
%% Construct Continuous graph 
     % disp('Constructing adjacent graph...');    Konstantia  Continuous
     % graph
     isCyclic =isCyclicc;  
     [P,d]=size(X); 
     n=P;
     if length(Pindex)==0;
        [c,Pindex] = histc(1:1:n,[1,P]);
     end
     sG=zeros(n,n);
     s=1;
     for i=1:n
        sG(i,i) = 1;
        if(i>P(s))
           s=s+1;
        end
        for j=1:tk
           if(i-j>0 && Pindex(i-j)==Pindex(i))
              sG(i,i-j) = 1;
           elseif(isCyclic==1)
              sG(i,P(s)+j-tk) = 1;
           end
           if(i+j<=n && Pindex(i+j)==Pindex(i))
               sG(i,i+j) = 1;
           elseif(isCyclic==1)
               LP = P(s);
               if(s>1)
                  LP = LP -  P(s-1);
               end
               sG(i,i-LP+j) = 1;
           end
        end
     end
 
     DX = DX - diag(diag(DX)); 
     for i=1:n
        sI = find(sG(i,:)>0);
        sG(i,sI) =DX(i,sI)';
        sG(i,i) = 0;
     end
     sG(1:size(sG,1)+1:end) = 1e-7;
     sG = max(sG, sG');
     sG = sG .^ 2;

end

