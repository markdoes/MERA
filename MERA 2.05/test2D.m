% two-dimensional demonstration code for MERA 2.0

clear all 
close all

% Create Example 2D Data
% Two Components, Bi-exp in two dimensions


SNR = 1000; % signal to noise ratio

% amplitude and time constants for signal component 1
M1 = .4; Ta1 = .020; Tb1 = 0.35; 
% amplitude and time constants for signal component 2
M2 = .6; Ta2 = .0750; Tb2 = .25; 
% amplitude and time constants for signal component 3
M3 = .04; Ta3 = .0750; Tb3 = 1; 

% sampling times in first dimension (rows)
Na = 128;
ta = .001*(1:Na)'; 
% sampling times in second dim (columns),
Nb = 64;
tb = logspace(log10(0.02),log10(5),Nb)'; 

D = zeros(Na,Nb);

for i=1:length(tb)
    for j=1:length(ta) 
        D(j,i) = M1 * exp(-(tb(i))./Tb1) *  exp(-(ta(j))./Ta1) ...
            + M2 * exp(-(tb(i))./Tb2) *  exp(-(ta(j))./Ta2) ... 
            + M3 * exp(-(tb(i))./Tb3) *  exp(-(ta(j))./Ta3);
    end
end

SNR = 1000;

v = randn(size(D))./SNR;% Gaussian random noise
Dv = (D+v);% Noisy signal 
data.D = Dv;
data.t = ta;
data.t2 = tb;

%

%% Calling MERA
% profile -memory on
clear fitting
close all
fitting.twoD = 'y';
 analysis.interactive = 'y'
fitting.regtyp = 'me';
fitting.rangeT = [10e-3,.25];
fitting.rangeT2 = [20e-3,1];
fitting.numbergauss = 2
fitting.regadj = 'manual';
analysis.graph = 'y';
fitting.regweight = 0.002;
fitting.numberT = 25;
fitting.numberT2 = 25;
analysis.extract = 'auto'

    [out2D,fitting_out]=MERA(data,fitting,analysis);
