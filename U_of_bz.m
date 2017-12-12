function [U,cc,Psic,u,v,r,con] = U_of_bz( beta,zeta,alpha )
%U_OF_BZ Summary of this function goes here
%   Detailed explanation goes here
eps=1e-6;

a1=6*beta*(1+2*alpha)*sinh(zeta/2)-6*alpha*beta*zeta*cosh(zeta/2);
a2=4*alpha*sinh(zeta/4)^2;
b0=-4*sinh(zeta/4)^2/cosh(zeta/2);

rmin=1;
rmax=64;
N=128;
n_modes=6;
[x,Dx]=chebdif(N,n_modes);
ra=(rmin-rmax)/2;
r=rmin-ra+ra*x;
D1=Dx(:,:,1)/ra;
D2=Dx(:,:,2)/ra^2;

cc=zeros(N,n_modes);

for i=1:50
    [Fu,Psic,u,v]=outer(r,alpha,a1,a2,b0,cc,D1,D2);
    DF=zeros(N*n_modes);
    for rd=1:N
        for md=1:n_modes
            cc0=cc;
            cc0(rd,md)=cc0(rd,md)+eps;
            tempv=(outer(r,alpha,a1,a2,b0,cc0,D1,D2)-Fu)/eps;
            DF(:,(md-1)*N+rd)=tempv;
        end
    end
    du=-DF\Fu;
    cc=cc+reshape(du,N,[]);
    if norm(Fu)<1e-10
        disp('stop at residule: ');
        disp(num2str(norm(Fu)));
        con=1;
        break;
    end
    if i==50
        disp('not converging.');
        con=0;
    end
end
U=b0*cc(1,1)/3;
end

