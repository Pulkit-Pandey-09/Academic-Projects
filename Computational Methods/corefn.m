% Temperature approximation of all finite elements within the core tissue is calculated by this file.
clear i ii width p r k Cp w
% 'i' is a spatial index that defines the position of the finite 
% element in the tissue type being evaluated.
i=1:n;

% 'p' is density in [kg/m^3]. 
%'Cp' is specific heat capacity [J/(kg-K)]. 
%'k' is thermal conductivity in [W/(m-K)]. 
%'w' is basal blood perfusion rate in [m^3 of blood/ (m^3 of tissue-s]. 
p=970; 
Cp=3650;
k=0.47; 

% 'o' is the Stefan-Boltzmann constant in [W/(m^2-K^4)]
o=5.67*10^-8; 

% Density and specific heat of blood are represented by 'pblood' in [kg/m^3]
% and 'Cpblood' in [J/(kg-K)], respectively. 
pblood=1059; 
Cpblood=3850; 

% 'width' is the width of each finite element within the tissue group. 'r' 
% is the radius of the outer circle that bounds the finite shell. 
width=Rcore;                                                                %************ DOUBT ************
r=width.*i; 

% The mass of blood entering the core tissue is the sum of blood leaving the peripheral tissue groups. 
%'btotal' is the mass flow rate of blood entering the core tissue in units of [kg/s]. The basal perfusion rate of 
% core tissue 'wtotal' is calculated by dividing 'btotal' by the volume of the core tissue. 'wcore' is the basal blood perfusion rate of the core
% The volume of the core tissue element is represented by the variable 'Vcore' in units of [m^3].
btotal=(Vmuscle*543*10^-6+Vct*90*10^-6+Vfat*77*10^-6+Vskin*362*10^-6);
Vcore=Rcore^2*pi*L;
wtotal=btotal/Vcore;

% The mass fraction of blood sourced from the bone, muscle, connective, 
% fat, and skin tissues are represented by xbone, xmuscle, xct, xfat, and 
% xskin. These values are used to calculate the mass weighted average
% temperature of blood entering the core tissue. 
xbone=0;
xmuscle=0.87;
xct=0.02;
xfat=0.04;
xskin=0.07;
% The following for loop solves the approximated energy balance around 
% each tissue group in each discrete element at the time point defined in the 
% mainfn. file. if statements are included to account for differences in heat 
%transfer for different finite elements. The following coefficients are used:
% A= Conduction term
% B= Convection due to blood flow in the microvasculature
% C= Net heat generated and consumed by basal metabolism 
% Rout= Radiation emitted by the finite element
% Rin= Radiation entering the finite element 

% 'width' is the width of each finite element within the tissue group. 
%'r' is the radius of the outer circle that bounds the finite shell. 
%'deltat' is the time interval over which change in temperature is approximated. 
for ii=1:n
 
 if ii==1
 A=(1/width)*deltat*2*k*r(ii)*(Tcore(ii+1,tt)-Tcore(ii,tt));
 
 B=deltat.*wtotal.*pblood.*Cpblood.*(xbone*sum(Tbone(:,tt))/n+xmuscle*sum(Tmuscle(:,tt))/n+xct*sum(Tct(:,tt))/n+xfat*sum(Tfat(:,tt))/n+xskin*sum(Tskin(:,tt))/1-Tcore(ii,tt)); 
 C=828.*deltat;
 Rout=-(o*Tcore(ii,tt)^4*2*pi*r(ii)*L)/Vcore;
 Rin=(o*Tcore(ii+1,tt)^4*2*pi*r(ii)*L)/Vcore;
 
 Tcore(ii,tt+1)=Tcore(ii,tt)+(1/p)*(1/Cp)*(A+B+C+Rin+Rout);
 end
 
 if 1<ii && ii<n
 A=(1/width)*deltat*2*k*(r(ii)*(Tcore(ii+1,tt)-Tcore(ii,tt))-r(ii-1)*(Tcore(ii,tt)-Tcore(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 
B=deltat.*wtotal.*pblood.*Cpblood.*(xbone*sum(Tbone(:,tt))/n+xmuscle*sum(Tmuscle(:,tt))/n+xct*sum(Tct(:,tt))/n+xfat*sum(Tfat(:,tt))/n+xskin*sum(Tskin(:,tt))/1-Tcore(ii,tt)); 
C=828.*deltat;
 Rout=-(o*Tcore(ii,tt)^4*2*pi*r(ii)*L+o*Tcore(ii,tt)^4*2*pi*r(ii-1)*L)/Vcore;
 Rin=(o*Tcore(ii+1,tt)^4*2*pi*r(ii)*L+o*Tcore(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vcore;
 
 Tcore(ii,tt+1)=Tcore(ii,tt)+(1/p)*(1/Cp)*(A+B+C+Rin+Rout);
 
 end
 
 if ii==n
 A=(1/width)*deltat*2*k*(r(ii)*(Tbone(1,tt)-Tcore(ii,tt))-r(ii-1)*(Tcore(ii,tt)-Tcore(ii-1,tt)))./(r(ii)^2-r(ii-1)^2);
 
B=deltat.*wtotal.*pblood.*Cpblood.*(xbone*sum(Tbone(:,tt))/n+xmuscle*sum(Tmuscle(:,tt))/n+xct*sum(Tct(:,tt))/n+xfat*sum(Tfat(:,tt))/n+xskin*sum(Tskin(:,tt))/1-Tcore(ii,tt)); 
 C=828.*deltat;
 Rout=-(o*Tcore(ii,tt)^4*2*pi*r(ii)*L+o*Tcore(ii,tt)^4*2*pi*r(ii-1)*L)/Vcore;
 Rin=(o*Tbone(1,tt)^4*2*pi*r(ii)*L+o*Tcore(ii-1,tt)^4*2*pi*r(ii-1)*L)/Vcore;
 
 Tcore(ii,tt+1)=Tcore(ii,tt)+(1/p)*(1/Cp)*(A+B+C+Rin+Rout);
 
 end
 
end