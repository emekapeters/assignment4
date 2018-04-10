%%Emeka Peters - 100953293
%%ELEC 4700 - Assignment 4 - Question 1

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

vin = 10;
%Vm = [V1, V2, IL, V3, I3, V4, V0];

F = zeros(1, 7);
F(1) = vin;

%DC Sweep
V = G\F.';

vinmat = -10: 0.5: 10;
voutmat = zeros(1, length(vinmat));
v3mat = zeros(1, length(vinmat));

for i = 1:length(vinmat)
    vin = vinmat(i);
    F(1) = vin;
    V = G\F.';
    voutmat(i) = V(7);
    v3mat(i) = V(4);
end

figure(1);
plot(vinmat, voutmat);
grid on
title('Vout vs Vin');
xlabel('Vin (V)');
ylabel('Vout (V)');

figure(2);
plot(vinmat, v3mat);
grid on
title('V3 vs Vin');
xlabel('Vin (V)');
ylabel('V3 (V)');

%AC Simulation
npoints = -7:0.1:10;
wmat = zeros(1, length(npoints));

for i = 1:length(npoints) %Creating omega values
    wmat(i) = 10 ^ (npoints(i));
end

gains = zeros(1, length(npoints));

vin = 10;
F(1, :) = 0;
F(1) = vin;

for i = 1:length(npoints)
    w = wmat(i);
    V = (G + 1i * w * Cm)\F.';
    gains(i) = 20 * log10(V(7)/V(1));
end

figure(3)
semilogx(wmat, gains);
title('Bode Plot of Gain Vo/V1')
grid on
xlabel('Gain (dB)');
ylabel('w(rad/s)');


%Monte Carlo Simulation

vin = 10;
F(1) = vin;

Cs = normrnd(0.25, 0.05, 1, 1000);

Cm = zeros(7, 7);
Cm(2, 1) = -C;
Cm(2, 2) = C;
Cm(3, 3) = -L;

voutmat = zeros(1, 1000);

for i = 1:1000
    C = Cs(i);
    Cm(2, 1) = -C;
    Cm(2, 2) = C;
    V = (G + 1i * pi * Cm)\F.';
    voutmat(i) = real(V(7));
end
gains = 20 .* log10(voutmat./vin);
figure(4)
hist(gains,50);
title('Histogram of Gain');
xlabel('Capacitance (F)');
ylabel('Gain (dB)');

