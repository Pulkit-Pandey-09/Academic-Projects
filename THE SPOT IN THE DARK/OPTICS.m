close all
clc
clear all

x = 100*[-500:500]/1000;
y = 100*[-500:500]/1000;

load -ascii Er5
load -ascii Ei5

I5 = Er5.^2 + Ei5.^2 + 1e-30;

IRGB(:,:,1) = 0*I5;
IRGB(:,:,2) = I5;
IRGB(:,:,3) = 0*I5;

m = min(min(min(IRGB)));
M = max(max(max(IRGB)));
figure(1)
im = (IRGB-m)/(M-m);
imagesc(x,y,im);
axis square;
%colorbar;
xlabel('x in mm');
ylabel('y in mm');

M = max(max(max(log10(IRGB))));
m = -2; 
figure(2)
im = (log10(IRGB)-m)/(M-m);
im = (im >= 0).*im;
imagesc(x,y,im);
axis square;
%colorbar;
xlabel('x in mm');
ylabel('y in mm');

m = min(min(min(sqrt(IRGB))));
M = max(max(max(sqrt(IRGB))));
figure(3)
im = (sqrt(IRGB)-m)/(M-m);
im = (im >= 0).*im;
imagesc(x,y,im);
axis square;
%colorbar;
xlabel('x in mm');
ylabel('y in mm');

m = min(min(min(IRGB.^(0.4))));
M = max(max(max(IRGB.^(0.4))));
figure(4)
im = (IRGB.^(0.4)-m)/(M-m);
im = (im >= 0).*im;
imagesc(x,y,im);
axis square;
%colorbar;
xlabel('x in mm');
ylabel('y in mm');

load -ascii Iofr5;
figure(5)
plot(Iofr5(:,1),Iofr5(:,2),'b-',-Iofr5(:,1),Iofr5(:,2),'b-');
legend( '\lambda=0.532');
xlabel('Radius in mm')
ylabel('Log Intensity')
