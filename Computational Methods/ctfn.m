% Temperature approximation of all finite elements within the bone tissue is calculated by this file.
clear i ii width p r k Cp w
% For definitions and units of i, ii, p, Cp, k, o, pblood, Cpblood, width,
% r and w, refer to corefn. file. 
i=1:n; 
p=1085; 
k=0.47; 
Cp=3200; 
w=p*0.5/(100*60*1000);
pblood=1059; 
Cpblood=3850; 
width=(Rct-Rmuscle); 
r=(width.*i+Rmuscle); 
o=5.67*10^-8;
Vct=pi*L*(Rct^2-Rmuscle^2);
for ii=1:1:n
 if ii==1
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tct(ii+1,tt)-Tct(ii,tt))-Rmuscle*(Tct(ii,tt)-Tmuscle(n,tt)));
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tct(ii,tt));
 C=369.*deltat;
 Rout=-(o*Tct(ii,tt)^4*2*pi*r(ii)*L+o*Tct(ii,tt)^4*2*pi*Rmuscle*L)/Vct;
 
Rin=(o*Tct(ii+1,tt)^4*2*pi*r(ii)*L+o*Tmuscle(n,tt)^4*2*pi*Rmuscle*L)/Vct;
 
 Tct(ii,tt+1)=Tct(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
 if 1<ii && ii<n
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tct(ii+1,tt)-Tct(ii,tt))-r(ii-1)*(Tct(ii,tt)-Tct(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tct(ii,tt));
 C=369.*deltat;
 Rout=-(o*Tct(ii,tt)^4*2*pi*r(ii)*L+o*Tct(ii,tt)^4*2*pi*r(ii-1)*L)/Vct;
 Rin=(o*Tct(ii+1,tt)^4*2*pi*r(ii)*L+o*Tct(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vct;
 
 Tct(ii,tt+1)=Tct(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 if ii==n
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tfat(1,tt)-Tct(ii,tt))-r(ii-1)*(Tct(ii,tt)-Tct(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tct(ii,tt));
 C=369.*deltat;
 Rout=-(o*Tct(ii,tt)^4*2*pi*r(ii)*L+o*Tct(ii,tt)^4*2*pi*r(ii-1)*L)/Vct;
 Rin=(o*Tfat(1,tt)^4*2*pi*r(ii)*L+o*Tct(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vct;
 
 Tct(ii,tt+1)=Tct(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
end
