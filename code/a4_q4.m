%%Emeka Peters - 100953293
%%ELEC 4700 - Assignment 4 - Question 4

G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G0 = 1/1000;
alf = 100;
bet = 500;
gam = 2500;
C = 0.25;
L = 0.2;

G = [1, 0, 0, 0, 0, 0, 0;...
    -G1, G1+G2, -1, 0, 0, 0, 0;...
    0, 1, 0, -1, 0, 0, 0;...
    0, 0, -1, G3, 0, 0, 0;...
    0, 0, 0, 0, 0, 1, 0;...
    0, 0, 0, G3, -1, 0, 0;...
    0, 0, 0, 0, 0, -G4, G4+G0];

Cm = zeros(7, 7);
Cm(2, 1) = -C;
Cm(2, 2) = C;
Cm(3, 3) = -L;

%Vm = [V1, V2, IL, V3, I3, V4, V0];

F = zeros(1, 7);

V = G\F.';
I3 = F(5);
J = zeros(7, 7);
J(5, 5) = -alf - (2*bet*I3) - (3 * gam * (I3^2));
B = zeros(1, 7);
%time step
dt = 1/1000;

A = (Cm / dt) + G;

f = 1/0.03;
vin = @(t) sin(2 * pi * 33 *t);
Vm = zeros(7, 1);
Vm = G\F.';
Vmm = zeros(7, 1);
t = 0:dt:1;
voutmat = zeros(1, 1001);
vinmat = zeros(1, 1001);
conv = 0;

for i = 1:1001
    F(1) = vin(t(i));
    conv = 0;
    dvf = Cm/dt + G - J;
    while(conv < 200)
        
        I3 = Vm(5);
        J = zeros(7, 7);
        J(5, 5) = -alf - (2*bet*I3) - (3 * gam * (I3^2));
        B(6) = (alf * I3) + (bet * (I3^2)) + (gam * (I3^3));
        
        F2 = ((Cm/dt + G) * Vm) - ((Cm/dt)*Vmm) - F' - B';
        
        H = Cm/dt + G - J;
        dVm = pinv(H) * F2;
        
        Vm  = Vm + dVm;
        
        if(max(abs(dVm)) <= 0.001)
            break;
        end
        
        conv = conv + 1;
    end
    
    voutmat(i) = Vmm(7);
    vinmat(i) = F(1);
    Vmm = Vm;
end

figure(1);
subplot(2,1,1);
plot(t,vinmat);
title('Non Linear Simulation: Vin vs Time');
grid on

subplot(2,1,2);
plot(t,voutmat);
title('Non Linear Simulation: Vout over Time');
grid on


fs = 1000;
fvout = fft(voutmat);
n = length(voutmat);

Y = fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(2);
plot(fshift,powershift);
title('Vout Spectrum: Non Linear Simulation');


fvin = fft(vinmat);
n = length(vinmat);

Y = fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y) .^ 2/n;     % zero-centered power range

figure(3)
plot(fshift,powershift)
title('Vin Spectrum: Non Linear Simulation')


