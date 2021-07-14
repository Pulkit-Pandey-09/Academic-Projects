% Temperature approximation of all finite elements within the bone tissue is calculated by this file.
clear i ii deltar p r k Cp w
% For definitions and units of ii, p, Cp, k, o, pblood, Cpblood, deltar,
% r and w, refer to corefn. file. 
% I is the spectral irradiance of blue light [W/m^2-nm]
I=736*10^-6*100^2; 

p=1085; 
k=0.47; 
Cp=3680; 
w=p*0.5/(100*60*1000);
pblood=1059; 
Cpblood=3850; 
deltar=(Rskin-Rfat); 
r=Rskin; 
o=5.67*10^-8; 
Vskin=pi*L*(Rskin^2-Rfat^2);
ii=1;

Twall=33+273;
A=(1/deltar)*deltat*2*k*(Rskin*(Twall-Tskin(ii,tt))/2-1*Rfat*(Tskin(ii,tt)-Tfat(n,tt)))./(Rskin^2-Rfat^2);
B=deltat.*w.*pblood.*Cpblood.*(sum(Tcore(:,tt))/n-Tskin(ii,tt));
C=368*deltat;
%0.5*Rskin*(33+273-Tskin(ii,tt))
% The following are factors that were not pertinent to heat transfer of
% other tissue groups:
% The coefficient 'D' represents convective heat transfer with the 
% surrounding air. 'Tair' is the temperature of air in [K] and 'h' is the 
% convective heat transfer coefficient in [W/m^2-K]. 
Tair= 33+273;
h= 4.94;
D=-h*BSA*(Tskin(ii,tt)-Tair).*deltat/2;
% 'E' is the blue light incident on and absorbed by the skin.
E=BSA*I*50.28*deltat/2;
% 'F' is the radiative heat emmitted by the surrounding wall and the 
% environment added to radiation emitted by the fat layer under the skin.
% 'Twall' is the temperature of the blanket and surrounding walls,
% equal to 33 degrees Celsius. 
F=BSA*o*Twall.^4*deltat+(o*Tfat(n,tt)^4*BSA*deltat);
% 'G' is the heat transferred out of the skin as radiation
G=-(2*pi*Rskin*L*o*Tskin(ii,tt).^4+2*pi*Rfat*L*o*Tskin(ii,tt).^4)*deltat;
% Heat loss secondary to diffusive evaporation of water out the skin is 
% described by the coefficient 'H'. 
H=-BSA*0.06*2.2*(10.^(8.07-1730.63/(233.43+Tskin(ii,tt)-273))-37.7)*4.94*deltat/2;
Tskin(ii,tt+1)=Tskin(ii,tt)+(1./p).*(1./Cp).*(A+B+C+(1/Vskin)*(D+E+F+G+H));
