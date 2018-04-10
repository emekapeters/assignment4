G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G0 = 1/1000;
alf = 100;

C = 0.25;
L = 0.2;

%Vm = [V1, V2, IL, V3, I3, V4, V0, In];% In parameter added, increased C
%and G matrix size

G = [1, 0, 0, 0, 0, 0, 0, 0;...
    -G1, G1+G2, -1, 0, 0, 0, 0, 0;...
    0, 1, 0, -1, 0, 0, 0, 0;...
    0, 0, -1, G3, 0, 0, 0, 1;...
    0, 0, 0, 0, -alf, 1, 0, -alf;...
    0, 0, 0, G3, -1, 0, 0, 1;...
    0, 0, 0, 0, 0, -G4, G4+G0, 0;...
    0, 0, 0, 0, 0, 0, 0, 1];

%Cap 1 = 0.000001F
Cn = 0.01;
Cm = zeros(8, 8);
Cm(2, 1) = -C;
Cm(2, 2) = C;
Cm(3, 3) = -L;
Cm(4, 4) = Cn;
Cm(6, 4) = Cn;
Cm(5, 4) = -alf * Cn;

%I3 was replaced with Itot which is the total of the currents from the
%resistor, capacitor, and current source

F = zeros(1, 8);

%time step
dt = 1/1000;

A = (Cm / dt) + G;

dt = 1/1000;
vin = @(t) exp(-0.5 *((t - 0.06)/(0.03)) ^ 2);
Vm = zeros(8, 1);
voutmat = zeros(1, 1001);
vinmat = zeros(1, 1001);
t = 0:dt:1; 

for i = 1:1001
    In = randn * 0.001; %Random noise current
    F(1) = vin(t(i));
    F(8) = In;
    Vmm = A\((Cm * Vm/dt) + F');
    voutmat(i) = Vmm(7);
    vinmat(i) = Vmm(1);
    Vm = Vmm;
end


figure(1)
subplot(2,1,1)
plot(t,vinmat)
title('Noisy Gaussian Function: Vin vs Time')
grid on

subplot(2,1,2)
plot(t,voutmat)
title('Noisy Guassian Function: Vout vs Time C = 0.01F')
grid on

   
fs = 1000;

fvin = fft(vinmat);
n = length(vinmat);
Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(2)
plot(fshift, powershift)
title('Vin: Noisy Gaussian Frequency Spectrum')

fvout = fft(voutmat);
n = length(voutmat);
Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(3)
plot(fshift,powershift)
title('Vout: Noisy Gaussian Frequency Spectrum C = 0.01F')

%Cap 2 = 0.0001F
Cn = 0.0001;
Cm = zeros(8, 8);
Cm(2, 1) = -C;
Cm(2, 2) = C;
Cm(3, 3) = -L;
Cm(4, 4) = Cn;
Cm(6, 4) = Cn;
Cm(5, 4) = -alf * Cn;

%I3 was replaced with Itot which is the total of the currents from the
%resistor, capacitor, and current source

F = zeros(1, 8);

%time step
dt = 1/1000;

A = (Cm / dt) + G;

dt = 1/1000;
vin = @(t) exp(-0.5 *((t - 0.06)/(0.03)) ^ 2);
Vm = zeros(8, 1);
voutmat = zeros(1, 1001);
vinmat = zeros(1, 1001);
t = 0:dt:1; 

for i = 1:1001
    In = randn * 0.001; %Random noise current
    F(1) = vin(t(i));
    F(8) = In;
    Vmm = A\((Cm * Vm/dt) + F');
    voutmat(i) = Vmm(7);
    vinmat(i) = Vmm(1);
    Vm = Vmm;
end


figure(4)
subplot(2,1,1)
plot(t,vinmat)
title('Noisy Gaussian Function: Vin vs Time')
grid on

subplot(2,1,2)
plot(t,voutmat)
title('Noisy Guassian Function: Vout vs Time C = 0.0001F')
grid on

   
fs = 1000;

fvin = fft(vinmat);
n = length(vinmat);
Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(5)
plot(fshift, powershift)
title('Vin: Noisy Gaussian Frequency Spectrum')

fvout = fft(voutmat);
n = length(voutmat);
Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(6)
plot(fshift,powershift)
title('Vout: Noisy Gaussian Frequency Spectrum C = 0.0001F')

%Cap 3 = 0.001F
Cn = 0.001;
Cm = zeros(8, 8);
Cm(2, 1) = -C;
Cm(2, 2) = C;
Cm(3, 3) = -L;
Cm(4, 4) = Cn;
Cm(6, 4) = Cn;
Cm(5, 4) = -alf * Cn;

%I3 was replaced with Itot which is the total of the currents from the
%resistor, capacitor, and current source

F = zeros(1, 8);

%time step
dt = 1/1000;

A = (Cm / dt) + G;

dt = 1/1000;
vin = @(t) exp(-0.5 *((t - 0.06)/(0.03)) ^ 2);
Vm = zeros(8, 1);
voutmat = zeros(1, 1001);
vinmat = zeros(1, 1001);
t = 0:dt:1; 

for i = 1:1001
    In = randn * 0.001; %Random noise current
    F(1) = vin(t(i));
    F(8) = In;
    Vmm = A\((Cm * Vm/dt) + F');
    voutmat(i) = Vmm(7);
    vinmat(i) = Vmm(1);
    Vm = Vmm;
end


figure(7)
subplot(2,1,1)
plot(t,vinmat)
title('Noisy Gaussian Function: Vin vs Time')
grid on

subplot(2,1,2)
plot(t,voutmat)
title('Noisy Guassian Function: Vout vs Time C = 0.001F')
grid on

   
fs = 1000;

fvin = fft(vinmat);
n = length(vinmat);
Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(8)
plot(fshift, powershift)
title('Vin: Noisy Gaussian Frequency Spectrum')

fvout = fft(voutmat);
n = length(voutmat);
Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(9)
plot(fshift,powershift)
title('Vout: Noisy Gaussian Frequency Spectrum C = 0.001F')