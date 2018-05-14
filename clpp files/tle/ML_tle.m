function [mappedX,mapping] = ML_tle(X,no_dims,tk,tw,stdAcceptance,isCyclic,smoothSize,mintab)
% ML_tle - reduced dimensionality of data using Temporal Laplacian Eigenmaps algorithm
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    if ~exist('stdAcceptance', 'var')
        stdAcceptance = 1.0; % std difference between reptitions
    end
    if ~exist('tk', 'var')
        tk=1;  % neighnours in time
    end
    if ~exist('tw', 'var')
        tw=20;  %half of the length of the fragment
    end
    if ~exist('isCyclic', 'var')
        isCyclic = 0;
    end 
    if(~exist('smoothSize','var'))
        smoothSize=30;
    end    
    if ~exist('mintab', 'var')
        mintab = [];
    end
    
    [P, d] = size(X);    
    DX = ML_tleL2(X',X',1);    
    [DTW1,DTW2] = ML_tleGetMatrixOfDTW(X,tw,P);    
    [G]=ML_tleGetDTW_CTN(X,DX,tw,stdAcceptance,DTW1,DTW2,P,smoothSize,mintab);
    
    n=P;
    [c,Pindex] = histc(1:1:n,[1,P]);
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
        sG(i,sI) = DX(i,sI)';
        sG(i,i) = 0;
    end
    sG(1:size(sG,1)+1:end) = 1e-7;
    sG = max(sG, sG');
    sG = sG .^ 2;
    
    mG = max(max(G));
    msG = max(max(sG));
    mG = max([mG msG]);
    
    G = G ./ mG;
    sG = sG ./ mG;
    
    G(G ~= 0) = exp(-G(G ~= 0));
    sG(sG ~= 0) = exp(-sG(sG ~= 0));    
    
    D = diag(sum(G, 2));
    L = D - G;
    L(isnan(L)) = 0; D(isnan(D)) = 0;
	L(isinf(L)) = 0; D(isinf(D)) = 0;
    
    sD = diag(sum(sG, 2));
    sL = sD - sG;    
    sL(isnan(sL)) = 0; sD(isnan(sD)) = 0;
	sL(isinf(sL)) = 0; sD(isinf(sD)) = 0;

    M = L+sL;
    D = D+sD;
    
	M(isnan(M)) = 0;
	M(isinf(M)) = 0;
    M = max(M, M');  
    
    D(isnan(D)) = 0;
	D(isinf(D)) = 0;
    D = max(D, D');  
    
    options.disp = 0;
    options.isreal = 1;
    options.issym = 1;
    options.v0=M(:,1);
    try
        tol=0;
        [eigenvecs, eigenvals] = eigs(M,D, max(no_dims) + 1, tol, options);
    catch
        warning('Exception in eigs, we try different tol');
        tol=1e-5;
        [eigenvecs, eigenvals] = eigs(M,D, max(no_dims) + 1, tol, options);
    end
    
    mappedX = cell(length(no_dims),1); 
    [eigenvals, ind] = sort(diag(eigenvals), 'ascend');
    eigenvals = eigenvals(2:no_dims + 1);
    
    for d=1:length(no_dims)
        mappedX{d} = real(eigenvecs(:,ind(2:no_dims(d) + 1)));
    end
    
    mapping.X = X';
    mapping.val = eigenvals;
    mapping.M = M;  
end