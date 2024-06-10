% Matlab script developed using version R2023a for synthesis of array
% pattern using the radiaton pattern of a square patch antenna element in 
% Phi and Theta polarized components, obtained from ANSYS HFSS version 
% 2022 and saved to realized_gain.csv.
% 
% The input impedance of the patch elements are matched to 50ohm.
% 
% The configuration of the array is rectangular (along x-y plane).
% Main input parameters are: 
% Nx = # array elements in the x dimension 
% Ny = # array elements in the y dimension. 
% dx = Center-to-center element spacing along the x dimensions
% dy = Center-to-center element spacing along the y dimensions
% PhiScan_deg = Scan angle in the Phi dimension (in degree)
% ThetaScan_deg = Scan angle in the Theta dimension (in degree)
%
% Output parameters are:
% Array pattern saved in the same format as the called element pattern
% with filename "arraypattern.mat"
% Plots of element/array patterns in 3D as well as Phi and Theta cuts. 
%
% Developed for the HEU 6GTandem Project by
% Yuyan Cao and Buon Kiong Lau, Lund University, Sweden
%
% Gilles Callebaut: correct small typo in version of 20 Feb 2024 (PhiScan_deg and phiScan_deg were used)
% Current version: 10 Jun 2024

close all;
clear;
clc;

M = readmatrix ('realized_gain.csv');
M_dB = M;

M_dB(:,3) = 10*log10(M(:,3)); % by definition, stored E-field data from Ansys HFSS is abs(Ephi)^2 abs(Etheta)^2 
M_dB(:,4) = 10*log10(M(:,4));

phi_deg = M(:, 1); 
theta_deg = M(:,2);
FieldPhi= sqrt(M(:, 3)); % converting the data back to abs(Ephi) and abs(theta) from squared quantities
FieldTheta = sqrt(M(:, 4));

c = 3e8; % speed of light in m/s 
f = 145e9; % frequency in Hz
lambda = c/f; % vacuum wavelength in m
beta0 = 2*pi/lambda; % 

%___input parameters____________________________________________
Nx = 4; % number of antennas in the x dimension (rectangular array)
Ny = 4; % number of antennas in the y dimension
dx = 0.55*lambda; % inter-element spacing in x direction in m
dy = dx; % inter-element spacing in y direction in m
PhiScan_deg   = 0; % scan angle \phi in degree
ThetaScan_deg = 30; % scan angle \theta in degree
%________________________________________________________________
 
% Array factor in the look directions (φ and θ in degrees)  
PhiXFed_deg = rad2deg(beta0*dx*sind(ThetaScan_deg)*cosd(PhiScan_deg));
PhiYFed_deg = rad2deg(beta0*dy*sind(ThetaScan_deg)*sind(PhiScan_deg));



% Array Pattern synthesis
FieldThetaSyn = 0;
FieldPhiSyn   = 0;
PhiXY = 0; % initial phase

for i = 1:Nx
    for j = 1:Ny       
        
        PhiXArray_Rad = beta0*dx.*sind(theta_deg).*cosd(phi_deg) .*(i-1); % increased phase in x aixs 
        PhiYArray_Rad = beta0*dy.*sind(theta_deg).*sind(phi_deg) .*(j-1); % increased phase in y aixs
        PhiArray_Rad  = PhiXArray_Rad+PhiYArray_Rad; % total increased phase at each position
 
        PhiXY(i,j) = PhiXFed_deg*(i-1)+PhiYFed_deg*(j-1); %feeding phase at each position compared to reference antenna

        % element pattern * array factor
        FieldThetaSyn = FieldThetaSyn + FieldTheta.*exp( 1j.*PhiArray_Rad )*exp(-1j.*deg2rad(PhiXY(i,j)))/sqrt(Nx*Ny); % last exponential term applies conjugate beamforming
        FieldPhiSyn   = FieldPhiSyn   + FieldPhi.*exp( 1j.*PhiArray_Rad )*exp(-1j.*deg2rad(PhiXY(i,j)))/sqrt(Nx*Ny); % array weight needs to be normalized by 1/sqrt(#array elements)
    end
end




FieldPhiSyn_dB = 20*log10(abs(FieldPhiSyn));
FieldThetaSyn_dB = 20*log10(abs(FieldThetaSyn));

%___input parameters____________________________________________
arraypat=zeros(length(M(:,1)),6);
arraypat(:,1:2)=M(:,1:2); % \phi and \theta angles
arraypat(:,3)=FieldPhi; % Amplitude of normalized electric field (square root of realized gain) in \phi polarization for single element 
arraypat(:,4)=FieldTheta; % Amplitude of normalized electric field in \theta polarization for single element 
arraypat(:,5)=FieldPhiSyn; % Complex Amplitude of normalized electric field in \phi polarization for array 
arraypat(:,6)=FieldThetaSyn; % Complex Amplitude of normalized electric field in \theta polarization for array 
save arraypattern arraypat 
%________________________________________________________________

% Element and array pattern plots

figure(1); patternCustom(M_dB(:,3),M(:,2),M(:,1)); title('Figure 1: Element Realized Gain (Phi)');

figure(2); patternCustom(M_dB(:,4),M(:,2),M(:,1)); title('Figure 2: Element Realized Gain (Theta)');

figure(3); patternCustom(FieldPhiSyn_dB,M(:,2),M(:,1)); title('Figure 3: Array Realized Gain (Phi)');

figure(4); patternCustom(FieldThetaSyn_dB,M(:,2),M(:,1)); title('Figure 4: Array Realized Gain (Theta)');

figure(5); patternCustom(M_dB(:,3),M(:,2),M(:,1),CoordinateSystem="rectangular",Slice="phi",SliceValue =[0 90]);
hold on; patternCustom(FieldPhiSyn_dB,M(:,2),M(:,1),CoordinateSystem="rectangular",Slice="phi",SliceValue=[PhiScan_deg PhiScan_deg+90]);
title('Figure 5: Realized GainPhi Comparison');
legend1 = {'element phi=0', 'element phi=90',['array phi=', num2str(PhiScan_deg)],['array phi=', num2str(PhiScan_deg+90)]};
legend(legend1);

figure(6);
patternCustom(M_dB(:,4),M(:,2),M(:,1),CoordinateSystem="rectangular",Slice="phi",SliceValue=[0 90]);
hold on; patternCustom(FieldThetaSyn_dB,M(:,2),M(:,1),CoordinateSystem="rectangular",Slice="phi",SliceValue=[PhiScan_deg PhiScan_deg+90]);
title('Figure 6: Realized GainTheta Comparison');
legend2 = {'element phi=0', 'element phi=90',['array phi=', num2str(PhiScan_deg)],['array phi=', num2str(PhiScan_deg+90)]};
legend(legend2);

figure(7);
patternCustom(10*log10(abs(FieldThetaSyn).^2+abs(FieldPhiSyn).^2),M(:,2),M(:,1));
title('Figure 7: Array Realized Gain');

% Using the saved data in arraypattern.mat
% patternCustom(10*log10(abs(arraypat(:,3)).^2+abs(arraypat(:,4)).^2),arraypat(:,2),arraypat(:,1));
