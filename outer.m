function [resF,ak,u,v,resf_out,d2cc] = outer( r,a,a1,a2,b0,cc_n,D1,D2 )
%OUTER Summary of this function goes here
%   Detailed explanation goes here
global P_mx;
global dP_mx;
global absc;
global wts;

setup_modes;

[N,mc]=size(cc_n);
M=mc+1;
cc=[zeros(N,1),cc_n];
u=zeros(N,M);
v=zeros(N,M);
d2cc=zeros(N,M);
u(:,2)=b0*cc(1,2)/3*(1./r.^3-1);
v(:,2)=b0*cc(1,2)/3*(-1-1/2./r.^3);
ak=zeros(M,1);
ak(2)=b0*cc(1,2)/6;
for k=3:M
    ak(k)=b0*cc(1,k)/(6*k-8);
    u(:,k)=-k*(k-1)*ak(k)*(r.^(1-k)-r.^(-k-1));
    v(:,k)=-ak(k)*((3-k)*r.^(1-k)+(k-1)*r.^(-k-1));
end
for k=1:M
    d2cc(:,k)=(D2+diag(2./r)*D1-k*(k-1)*diag(1./r.^2))*cc(:,k);
end
resfR=d2cc*P_mx(1:M,:)-a*((u*P_mx(1:M,:)).*(D1*cc*P_mx(1:M,:))+...
    (v*dP_mx(1:M,:)).*(diag(1./r)*cc*dP_mx(1:M,:)));
resf=R2L(resfR,M);
resf_out=resf;
resf(1,1)=D1(1,:)*cc(:,1);
resf(N,1)=cc(N,1);
resf(1,2)=D1(1,:)*cc(:,2)-2*a2*b0*cc(1,2)-a1;
resf(N,2)=cc(N,2);
for k=3:M
    resf(1,k)=D1(1,:)*cc(:,k)-2*k*(k-1)*a2*ak(k);
    resf(N,k)=cc(N,k);
end
resF=reshape(resf(:,2:M),[],1);
end

