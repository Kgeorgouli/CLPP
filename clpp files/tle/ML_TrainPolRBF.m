function [Binv,Bfor,C,sigma] = ML_TrainPolRBF(X,Y,Nc,rbfType,clusterType,C)
%rbf - 0 thin plate, 1 - gaussian sigma=1, 2 - gaussian sigma unique, 3 - gaussian sigma differents, 4-real gaussian
%clusterType - 0 Kmeans, 1-EM (only for rbftype=4)

[N,d]=size(X);
if ~exist('C', 'var') || isempty(C)
    if ~exist('clusterType', 'var') || isempty(clusterType)
        opts = statset('MaxIter',200);
        [IDX,C] = kmeans(X,Nc,'EmptyAction','singleton','Replicates',200,'Options',opts);
    else
        if clusterType==0
            opts = statset('MaxIter',200);
            [IDX,C] = kmeans(X,Nc,'EmptyAction','singleton','Replicates',200,'Options',opts);
        else
            %EM
             estS = gmmb_fj(X, 'Cmax', Nc, 'thr', 1e-3,'maxloops',300);
            C = estS.mu';
            sigma = estS.sigma;
        end
    end
end
Nc = size(C,1);

if (rbfType == 1)
    sigma=1;
elseif (rbfType == 2)
    sigma = sum(sum(ML_tleL2(C',C',1)))/Nc;
elseif (rbfType == 3)
    for j=1:Nc
        values = X(IDX==j,:);
        values = (values-ones(size(values,1),1)*C(j,:)).^2;
        sigma(j,:) = sqrt( sum(values,1)/size(values,1) );
    end
    
else
    if ~exist('sigma', 'var')
        cont=1;
        for j=1:Nc
            values = X(IDX==j,:);
            if size(values,1)==1
                C(cont,:)=values;
                sigma(:,:,cont)=eye(d);
                cont=cont+1;
            elseif size(values,1)>0
                C(cont,:)=mean(values,1);
                sigma(:,:,cont)=cov(values);
                cont=cont+1;
            end

        %            [R,p]=chol(sigma(:,:,j));
        %             if p>0
        %                 sigma(:,:,j)=sqrtm(sigma(:,:,j));
        %             else
        %                 sigma(:,:,j)=R;
        %             end
        end
    end
end

H=ones(N,Nc);
tol=0.000000000001;
for i=1:N
    i
    for j=1:Nc
        if (rbfType == 3)
            H(i,j)= sum( ( ( (X(i,:) - C(j,:)) ./ sigma(j,:) ).^2 ) );% /2 );
        elseif (rbfType == 4)
            H(i,j)=  ((X(i,:) - C(j,:))*pinv(sigma(:,:,j))*(X(i,:) - C(j,:))'); %/2;
        else
            H(i,j)= norm(X(i,:) - C(j,:));
        end
    end
end
H(abs(H)<tol)=0;
if(rbfType==0)
    %r^2log(r)
    H(H>0) = H(H>0).*H(H>0).*log(H(H>0));
elseif (rbfType==1)||(rbfType==2)
    %exp(-r^2/sigma)
    H(H>0) = exp(-H(H>0).*H(H>0)/sigma);
else
    %exp(-r^2/sigma)
    H(H>0) = exp(-H(H>0));
end

Px=[ones(N,1) X];
Pc=[ones(Nc,1) C];

phi = [H Px; Pc' zeros(d+1,d+1)];
R = [Y; zeros(d+1,size(Y,2))];
Binv = pinv(phi)*R;
Bfor = pinv(Binv');

end