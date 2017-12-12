zeta=[0.959827883033704,1.757981276986011,2.386436292058683];
ind=-2:0.25:1.5;
beta=10.^(ind);
bn=length(ind);
U=zeros(3,bn);
Us=zeros(3,1);
for i=1:3
    for j=1:bn
        [U(i,j),~,~,~,~,con]=U_of_bz(beta(j),zeta(i),0.2);
        if con==0
            Us(i)=j-1;
            break;
        end
    end
end
%save('UvsB.mat','U','beta','-v7.3');