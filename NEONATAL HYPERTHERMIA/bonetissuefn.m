% Temperature approximation of all finite elements within the bone tissue is calculated by this file. 
clear i ii width p r k Cp w
% For definitions and units of i, ii, p, Cp, k, o, pblood, Cpblood, width,
% r and w, refer to corefn. file. 

% 'p' is density in [kg/m^3]. 
%'Cp' is specific heat capacity [J/(kg-K)]. 
%'k' is thermal conductivity in [W/(m-K)]. 
%'w' is basal blood perfusion rate in [m^3 of blood/ (m^3 of tissue-s]. 
% Density and specific heat of blood are represented by 'pblood' in [kg/m^3] and 'Cpblood' in [J/(kg-K)], respectively. 
i=1:n; 
p=1520; 
k=0.65; 
Cp=1700; 
w=0; 
pblood=1059; 
Cpblood=3850;

% 'width' is the width of each finite element within the tissue group. 
%'r' is the radius of the outer circle that bounds the finite shell. 
width=(Rbone-Rcore);  
r=(width).*i+Rcore;

% 'o' is the Stefan-Boltzmann constant in [W/(m^2-K^4)]
o=5.67*10^-8; 

Vbone=pi*L*(Rbone-Rcore)^2; 
for ii=1:1:n                   
 if ii==1
 
 A=(1/width)*deltat*2*k* (r(ii)*(Tbone(ii+1,tt)-Tbone(ii,tt))-Rcore*(Tbone(ii,tt)-Tcore(n,tt)))(r(ii)^2-Rcore^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tbone(ii,tt));
 C=368.*deltat;
 Rout=-(o*Tbone(ii,tt)^4*2*pi*r(ii)*L + o*Tbone(ii,tt)^4*2*pi*Rcore*L)/Vbone;
 
 Rin=(o*Tbone(ii+1,tt)^4*2*pi*r(ii)*L+o*Tcore(n,tt)^4*2*pi*Rcore*L)/Vbone;
 
 Tbone(ii,tt+1)=Tbone(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 
 if 1<ii && ii<n
 
 A =(1/width)*deltat*2*k*(r(ii)*(Tbone(ii+1,tt)-Tbone(ii,tt))-r(ii-1)*(Tbone(ii,tt)-Tbone(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B =deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/ncore-Tbone(ii,tt));
 C=368.*deltat;
 Rout=-(o*Tbone(ii,tt)^4*2*pi*r(ii)*L+o*Tbone(ii,tt)^4*2*pi*r(ii-1)*L)/Vbone;
 Rin=(o*Tbone(ii+1,tt)^4*2*pi*r(ii)*L+o*Tbone(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vbone;
 
 Tbone(ii,tt+1)=Tbone(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
 end
 if ii==n
 
 A=(1/width)*deltat*2*k*(r(ii)*(Tmuscle(1,tt)-Tbone(ii,tt))-r(ii-1)*(Tbone(ii,tt)-Tbone(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tbone(ii,tt));
 C=368.*deltat;
 Rout=-(o*Tbone(ii,tt)^4*2*pi*r(ii)*L+o*Tbone(ii,tt)^4*2*pi*r(ii-1)*L)/Vbone;
 Rin=(o*Tmuscle(1,tt)^4*2*pi*r(ii)*L+o*Tbone(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vbone;
 
 Tbone(ii,tt+1)=Tbone(ii,tt)+(1./p).*(1./Cp).*(A+B+C+Rout+Rin);
 
end
end
