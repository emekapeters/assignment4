%%Emeka Peters - 100953293
%%ELEC 4700 - Assignment 4 - Question 2

G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G0 = 1/1000;
alf = 100;

C = 0.25;
L = 0.2;

G = [1, 0, 0, 0, 0, 0, 0;...
    -G1, G1+G2, -1, 0, 0, 0, 0;...
    0, 1, 0, -1, 0, 0, 0;...
    0, 0, -1, G3, 0, 0, 0;...
    0, 0, 0, 0, -alf, 1, 0;...
    0, 0, 0, G3, -1, 0, 0;...
    0, 0, 0, 0, 0, -G4, G4+G0];

Cm = zeros(7, 7);
Cm(2, 1) = -C;
Cm(2, 2) = C;
Cm(3, 3) = -L;

%Vm = [V1, V2, IL, V3, I3, V4, V0];

F = zeros(1, 7);

V = G\F.';

%time step
dt = 1/100;

A = (Cm / dt) + G;

%i. Unit Step Input Simulation

Vm = zeros(7, 1);
t = 0:dt:1;
voutmat = zeros(1, length(t));
vinmat = zeros(1, length(t));

for i = 1:length(t)
    if(t(i) >= 0.03)
        F(1) = 1;
    else
        F(1) = 0;
    end
    
    Vmm = A\((Cm * Vm/dt) + F');
    voutmat(i) = Vmm(7);
    vinmat(i) = Vmm(1);
    Vm = Vmm;
end

figure(1);
subplot(2,1,1);
plot(t,vinmat);
title('Unit Step Input: Vin vs Time');
grid on

subplot(2,1,2);
plot(t,voutmat);
title('Unit Step Input: Vout over Time');
grid on


fs = 1000;
fvout = fft(voutmat);
n = length(voutmat);

Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(2);
plot(fshift,powershift);
title('Vout Spectrum: Unit Step');


fvin = fft(vinmat);
n = length(vinmat);

Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(3)
plot(fshift,powershift)
title('Vin Spectrum: Unit Step')



%ii. Sine Funtion Input Simulation
f = 1/0.03;
vin = @(t) sin(2 * pi * 33 *t);
Vm = zeros(7, 1);
t = 0:dt:1;
voutmat = zeros(1, length(t));
vinmat = zeros(1, length(t));

for i = 1:length(t)
    F(1) = vin(t(i));
    Vmm = A\((Cm * Vm/dt) + F');
    voutmat(i) = Vmm(7);
    vinmat(i) = Vmm(1);
    Vm = Vmm;
end

figure(4)
subplot(2,1,1)
plot(t, vinmat)
title('Sine Function Vin over Time')
grid on
     
subplot(2,1,2)
plot(t, voutmat)
title('Sine Function Vout over Time')
grid on


%Frequency Spectrums
  
fs = 1000;

fvin = fft(vinmat);
n = length(vinmat);
Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(5)
plot(fshift, powershift)
title('Vin: Sine Frequency Spectrum')

fvout = fft(voutmat);
n = length(voutmat);
Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(6)
plot(fshift,powershift)
title('Vout: Sine Frequency Spectrum')


% iii. Gaussian Function
dt = 1/100;
vin = @(t) exp(-0.5 *((t - 0.06)/(0.03)) ^ 2);
Vm = zeros(7, 1);
t = 0:dt:1;
voutmat = zeros(1, length(t));
vinmat = zeros(1, length(t));

for i = 1:length(t)
    F(1) = vin(t(i));
    Vmm = A\((Cm * Vm/dt) + F');
    voutmat(i) = Vmm(7);
    vinmat(i) = Vmm(1);
    Vm = Vmm;
end


figure(7)
subplot(2,1,1)
plot(t,vinmat)
title('Gaussian Function: Vin vs Time')
grid on

subplot(2,1,2)
plot(t,voutmat)
title('Guassian Function: Vout vs Time')
grid on

   
fs = 1000;

fvin = fft(vinmat);
n = length(vinmat);
Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(8)
plot(fshift, powershift)
title('Vin: Gaussian Frequency Spectrum')

fvout = fft(voutmat);
n = length(voutmat);
Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(9)
plot(fshift,powershift)
title('Vout: Gaussian Frequency Spectrum')

