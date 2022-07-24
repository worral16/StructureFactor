% CHEG825, HW9, P1
% 11/7/20
clear; close all;clc;
%% Import Data
csgr0 = importfile('csgr1.0.txt', 1, 1000);

%% Adjustment for erroneoues data
starter = 1;
endval = 1;
rdata = csgr0(:,1);
grdata = csgr0(:,2);
% Note: The rdata from the text file is normalized by the particle diameter 180 nm
rdata = rdata(starter:end-endval);
grdata = grdata(starter:end-endval);

%% Plot of output g(r) code from Frenkel and Smit code
figure(1)
plot(rdata,grdata)
ylabel(['g(r)/', char(949)])
xlabel('r/\sigma')
title('Pair Correlation Function')
rho = 1;
sLow = .0001;
sInc = sLow;
sUpper = .02;
scatteringAngleS = sLow:sInc:sUpper;
structureFactor = zeros(1,length(scatteringAngleS));
% convert rdata back to angstrom, user should modify this as needed
rdata = rdata*1800;

%% Numerical Fourier Transform Using Trapezoidal Integration
for s = sLow:sInc:sUpper
    
    y = zeros(1,length(grdata));
    
    for i = 1:length(grdata)
        y(i) = (grdata(i)-1) * sin(s*rdata(i))/(s*rdata(i))*rdata(i)^2;
    end
    % 10^-10 factor is for handling units in angstrom
    structureFactor(round(s*1/sLow)) = 1+2*pi*rho*trapz(rdata,y)*10^-10;
end

figure(2)
plot(scatteringAngleS,structureFactor)
xlabel("Q(A^{-1})")
ylabel("H(q)")
title('Structure Factor')
