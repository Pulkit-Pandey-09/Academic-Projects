clc
% Temperature model of a neonatal human undergoing blue light
% phototherapy for management of hyperbilirubinemia 

% The scalar 't' is the magnitude of the time domain in [s] and 'deltat' 
% is the time interval over which change in temperature is approximated.
t=24*60*60;            % 24 hrs
deltat=1;

% 'Vtotal' is the volume of the neonatal body in [m^3].
Vtotal=0.0031;
Vcore=0.152*Vtotal;
Vbone=Vtotal*0.12;
Vmuscle=Vtotal*0.441;
Vct=0.069*Vtotal;
Vfat=0.162*Vtotal;
Vskin=0.053*Vtotal;

% 'BSA' is the body surface are of the neonate, calculated using the Meban
% formula, in [m^2]
BSA=0.21;

% 'L' is the length of the approximated neonatal body in [m]
L=1.11;

% The radius of the central cylinder and of the outer radius the
% concentric shells representing bone, muscle, connective, fat, and skin
% tissue are 'Rcore', 'Rbonetissue', 'Rmuscle', 'Rct', 'Rfat', and 'Rskin', 
% respectively'.
Rcore=(Vcore/(pi*L))^(1/2);              %pi*r^2*L = Vcore
Rbone=((Vcore+Vbone)/(pi*L))^(1/2); 
Rmuscle=((Vcore+Vbone+Vmuscle)/(pi*L))^(1/2); 
Rct=((Vct+Vcore+Vbone+Vmuscle)/(pi*L))^(1/2); 
Rfat=((Vfat+Vct+Vcore+Vbone+Vmuscle)/(pi*L))^(1/2); 
Rskin=((Vskin+Vfat+Vct+Vcore+Vbone+Vmuscle)/(pi*L))^(1/2); 

% Each tissue element is divided into 'n' number of finite elements. 
% The variable ncore is equal to n and describes the number of discretized elements in the 
% core. Note that n is not used in the calculation of skin temperature. Skin is 
% treated as a lumped system. 
n=2;
ncore=n;

% The temperature of each finite elements in the core, bone, muscle, 
% connective, fat and skin tissue are given by the matrices 'Tcore', 
% 'Tbone', 'Tmuscle', 'Tct', 'Tfat', and 'Tskin', respectively, in units 
% of [K]. 
% The rows of each matrix defines the position of each finite 
% element within each tissue group. The columns of the temperature matrix is 
% the time axis and the time lapse from one column to the next is equivalent 
% to 'deltat'. 
% Zeros are used as place holders in each temperature matrix. 
% Calculated values of temperature are inputted by the other m.files.
Tcore=zeros(n,t); 
Tbone=zeros(n,t); 
Tmuscle=zeros(n,t); 
Tct=zeros(n,t); 
Tfat=zeros(n,t); 
Tskin=zeros(1,t); 

% Each finite element in each tissue group is set to an initial temperature
% of 310 [K], consistent with the assumption that initial temperature is
% uniform. 
Tcore(:,:)=35+273; 
Tbone(:,:)=35+273; 
Tmuscle(:,:)=35+273;
Tct(:,:)=35+273;
Tfat(:,:)=35+273;
Tskin(:,:)=35+273;

% Approximated temperature values are inputted into the temperature matrices 
% by the following m.files. There is a separate m.file for each tissue group. 
% These files are repeated over the domain of time to calculate 
% temperature approximations over the entire time domain. The variable 'tt'
% is the array of time points being approximated. 
for tt=1:t
 
corefn;
bonetissuefn;
musclefn;
ctfn;
fatfn;
skinfn;
end
hold off

x=1:1:t;
y=x/3600;
% The matrices CoreTemp, BoneTemp, MuscleTemp, CTTemp, FatTemp, and SkinTemp contain the mean value of temperature averaged over all discrete elements in each tissue group at every time point. 
CoreTemp=sum(Tcore(:,x))/n; 
BoneTemp=sum(Tbone(:,x))/n; 
MuscleTemp=sum(Tmuscle(:,x))/n;
CTTemp=sum(Tct(:,x))/n;  
FatTemp=sum(Tfat(:,x))/n;  
SkinTemp=Tskin(:,x); 

timeA = linspace(1,24,86400);  %1:86400:25;
 
plot(timeA,CoreTemp , 'LineWidth',1.2, timeA,BoneTemp ,'LineWidth',1.2 , timeA,MuscleTemp ,'LineWidth',1.2 , timeA,CTTemp, 'LineWidth',1.2 , timeA,FatTemp , 'LineWidth',1.2 , timeA,SkinTemp , 'LineWidth',1.2 );
%plot(timeA,BoneTemp,'LineWidth',1);
xlabel('TIME');
ylabel('TEMPERATURE');
legend('CoreTemp' , 'BoneTemp' , 'MuscleTemp' , 'CTTemp' , 'FatTemp' , 'SkinTemp'  );


