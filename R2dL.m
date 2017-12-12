function L_mx = R2dL( R_mx,n_modes )
%REAL_TO_LEGENDRE Summary of this function goes here
%   Detailed explanation goes here
global dP_mx;
global wts;
[n_rpts,~]=size(R_mx);
L_mx=zeros(n_rpts,n_modes);
for i=2:n_modes
    for j=1:n_rpts
        temp=0;
        for k=1:16
            temp=temp+R_mx(j,k)*dP_mx(i,k)*wts(k);
        end
        L_mx(j,i)=temp*(2*i-1)/2/(i-1)/i;            
    end
end         
end

