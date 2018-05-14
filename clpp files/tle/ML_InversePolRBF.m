function [P] = ML_InversePolRBF(X,C,Binv,sigma,rbfType)
%rbf - 0 thin plate, 1 - gaussian
    
    Nc = size(C,1);
    N = size(X,1);
    H=ones(N,Nc);
    tol=0.000000000001;    
    for i=1:N
        for j=1:Nc
            if (size(sigma,1)==1)&&(size(sigma,2)==1)
                H(i,j)= norm(X(i,:) - C(j,:));
            elseif (size(sigma,3)==1)
                H(i,j)= sum( ( ( (X(i,:) - C(j,:)) ./ sigma(j,:) ).^2 ) );% /2 );
            else
                H(i,j)=  ((X(i,:) - C(j,:))*pinv(sigma(:,:,j))*(X(i,:) - C(j,:))') ; %/2;
            end
        end
    end
    H(abs(H)<tol)=0;
    if(rbfType==0)
        %r^2log(r)            
        H(H>0) = H(H>0).*H(H>0).*log(H(H>0));
    else
        %exp(-r^2/sigma)
        if (size(sigma,1)==1)&&(size(sigma,2)==1)
            H(H>0) = exp(-H(H>0).*H(H>0)/sigma);
        else
            %exp(-r^2/sigma)
            H(H>0) = exp(-H(H>0));
        end
    end
    
    Px=[ones(N,1) X];
    H = [H Px];
    
    for i=1:N
        P(i,:)=H(i,:)*Binv;
    end
    
end