clear
close all
clc
addpath(strcat(fileparts(mfilename('fullpath')), '\code'));
format compact

%Universal constants
c.eRestMass = 9.109E-31; %kg
c.boltzmann = 1.381E-23; %J/K

%Setup data structures
system.x = 200E-9; %m
system.y = 100E-9; %m
system.Temp = 300; %K
system.Tau = 0.2E-12; %s

voltageX = 0.8; % V
voltageY = 0;

system.EfieldX = voltageX./system.x; %V/m
system.EfieldY = voltageY./system.y; %V/m

%Calculation to get "number of electrons" in system
density = 10.^19; %1/m^2
system.EDensity = density.*system.x.*system.y;

numOfParticles = 1000;
numLoops = 1000;


electron.effM = 0.26.*c.eRestMass;
electron.num = numOfParticles;
electron.x = zeros(1,numOfParticles);
electron.y = zeros(1,numOfParticles);
electron.vx = zeros(1, numOfParticles);
electron.vy = zeros(1, numOfParticles);

%Calculation for thermal velocity and mean free path
system.thermalV = sqrt(2.*c.boltzmann.*system.Temp./(electron.effM));
system.meanFreePath = system.thermalV.*system.Tau;

[electron.x, electron.y] = assignPosition(system.x, system.y, electron.num);

%Generate the time step as time it takes to travel 150th
%of the minimum dimension at speed of thermal velocity
timeStep = min([system.x system.y])./(100.*system.thermalV);

% figure(1);
% xlim([0 system.x]);
% ylim([0 system.y]);

%Part 1 of assignment
%meanFreePath(system, electron, numLoops, timeStep)

%Part 2 of assignment
[Ex, Ey] = highResist(10E-15, 1, voltageX, voltageY, 200, 100, 'none');


%Part 3 of assignment
bottleNeck(system, electron, numLoops, timeStep, 'specular', Ex, Ey)