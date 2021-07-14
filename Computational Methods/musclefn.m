% Temperature approximation of all finite elements within the bone tissue is 
%calculated by this file.
clear i ii width p r k Cp w
% For definitions and units of i, ii, p, Cp, k, o, pblood, Cpblood, width,
% r and w, refer to corefn. file. 
i=1:n; 
p=1085; 
k=0.51; 
Cp=3800; 
w=p*3/(100*60*1000);
pblood=1059; 
Cpblood=3850;
width=(Rmuscle-Rbone); 
r=(width.*i+Rbone); 
o=5.67*10^-8; 
Vmuscle=pi*L*(Rmuscle^2-Rbone^2);
for ii=1:1:n
 if ii==1
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tmuscle(ii+1,tt)-Tmuscle(ii,tt))-Rbone*(Tmuscle(ii,tt)-Tbone(n,tt)))/(r(ii)^2-Rbone^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tmuscle(ii,tt));
 C=684.*deltat;
 Rout=-(o*Tmuscle(ii,tt)^4*2*pi*r(ii)*L+o*Tmuscle(ii,tt)^4*2*pi*Rbone*L)/Vmuscle;
 
 Rin=(o*Tmuscle(ii+1,tt)^4*2*pi*r(ii)*L+o*Tbone(n,tt)^4*2*pi*Rbone*L)/Vmuscle;
 
 Tmuscle(ii,tt+1)=Tmuscle(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
 if 1<ii && ii<n
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tmuscle(ii+1,tt)-Tmuscle(ii,tt))-r(ii-1)*(Tmuscle(ii,tt)-Tmuscle(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tmuscle(ii,tt));
 C=684.*deltat;
 Rout=-(o*Tmuscle(ii,tt)^4*2*pi*r(ii)*L+o*Tmuscle(ii,tt)^4*2*pi*r(ii-1)*L)/Vmuscle;
 Rin=(o*Tmuscle(ii+1,tt)^4*2*pi*r(ii)*L+o*Tmuscle(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vmuscle;
 
 Tmuscle(ii,tt+1)=Tmuscle(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 if ii==n
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tct(1,tt)-Tmuscle(ii,tt))-r(ii-1)*(Tmuscle(ii,tt)-Tmuscle(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tmuscle(ii,tt));
 C=684.*deltat;
 Rout=-(o*Tmuscle(ii,tt)^4*2*pi*r(ii)*L+o*Tmuscle(ii,tt)^4*2*pi*r(ii-1)*L)/Vmuscle;
 Rin=(o*Tct(1,tt)^4*2*pi*r(ii)*L+o*Tmuscle(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vmuscle;
 
 Tmuscle(ii,tt+1)=Tmuscle(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
end
