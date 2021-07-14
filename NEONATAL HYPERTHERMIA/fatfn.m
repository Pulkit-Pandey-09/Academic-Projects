% Temperature approximation of all finite elements within the bone tissue is calculated by this file.
clear i ii deltar p r k Cp w
% For definitions and units of i, ii, p, Cp, k, o, pblood, Cpblood, deltar,
% r and w, refer to corefn. file. 
i=1:n; 
p=920; 
k=0.21; 
Cp=2300; 
w=p*0.5/(100*60*1000);
pblood=1059; 
Cpblood=3850; 
deltar=(Rfat-Rct)./n; 
r=(deltar.*i+Rct); 
o=5.67*10^-8; 
Vfat=pi*L*(Rfat^2-Rct^2);
for ii=1:1:n
 if ii==1
 
 A=(1/deltar)*deltat*2*k*(r(ii)*(Tfat(ii+1,tt)-Tfat(ii,tt))-Rfat*(Tfat(ii,tt)-Tct(n,tt)));
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tfat(ii,tt));
 C=300.*deltat;
 Rout=-(o*Tfat(ii,tt)^4*2*pi*r(ii)*L+o*Tfat(ii,tt)^4*2*pi*Rct*L)/Vfat;
 Rin=(o*Tfat(ii+1,tt)^4*2*pi*r(ii)*L+o*Tct(n,tt)^4*2*pi*Rct*L)/Vfat;
 
 Tfat(ii,tt+1)=Tfat(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
 if 1<ii && ii<n
 
 A=(1/deltar)*deltat*2*k*(r(ii)*(Tfat(ii+1,tt)-Tfat(ii,tt))-r(ii-1)*(Tfat(ii,tt)-Tfat(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tfat(ii,tt));
 C=300.*deltat;
 Rout=-(o*Tfat(ii,tt)^4*2*pi*r(ii)*L+o*Tfat(ii,tt)^4*2*pi*r(ii-1)*L)/Vfat;
 Rin=(o*Tfat(ii+1,tt)^4*2*pi*r(ii)*L+o*Tfat(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vfat;
 Tfat(ii,tt+1)=Tfat(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 if ii==n
 
 A=(1/deltar)*deltat*2*k*(r(ii)*(Tskin(1,tt)-Tfat(ii,tt))-r(ii-1)*(Tfat(ii,tt)-Tfat(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tfat(ii,tt));
 C=300.*deltat;
 Rout=-(o*Tfat(ii,tt)^4*2*pi*r(ii)*L+o*Tfat(ii,tt)^4*2*pi*r(ii-1)*L)/Vfat;
 Rin=(o*Tskin(1,tt)^4*2*pi*r(ii)*L+o*Tfat(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vfat;
 
 Tfat(ii,tt+1)=Tfat(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
end
