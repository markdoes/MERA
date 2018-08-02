% demonstration/test code for MERA 2.0
% 1D multiple exponential examples
%
clear variables,
close all
clc


Ntrials = 10; % number of independent noise realizations
SNR = 250; % signal to noise ratio at time t = 0



% create 1D signals with time constants T2 = [T2_1, T2_2, ..., T1_K]'
% and corresponding signal amplitudes f = [f_1, f_2, ..., f_K]'
T2 = [20e-3, 100e-3]';
f = [2 8]';

% mimic CPMG measurement in NMR with refocussing pulse angle theta (in
% degrees)
theta = 155;

% samples in time are uniformly space (not a requirement for MERA) with
% period TE
TE = 7.5e-3;
t = (1:50)'*TE;

% maketestdata1D returns a structure containing the signal data.D and
% independent variable, data.t
data = maketestdata1D(T2,f,t,theta,Ntrials,SNR);

% now set fitting and analysis options

% analysis.interactive = 'y' provide an interactive GUI. This is useful for
% testing out different fitting options, but will only process one signal
analysis.interactive = 'n'; 

% If your signal came from an NMR CPMG measurement (see literature
% references in code) setting fitting.B1fit = 'y' will fit out the
% refocussing pulse flip angle. This is useful when measuring NMR
% transverse relaxation with imaging, where the RF transmit field likely
% varies over the image field of view. This feature only works for 1D
% signals.
fitting.B1fit = 'y';
% If this option is used, you may choose a range and number of flip angles
% (theta) to test. If you do not test a sufficient number or range, you may
% not get the best fit. 
fitting.rangetheta = [130 180];
fitting.numbertheta = 10;


% Choose the type of spectral regularization. 'none' is obvious. 
% 'me' (min amplitude energy) and 'mc' (min curvature energy') are very common
% choices. 
% 'upen' is the "uniform penalty" regularization (see literature
% references in code) which allows regularization weighting to vary across
% the spectrum, which may be useful if varying degrees of smoothness is
% expected across the spectrum.'upen' can be very slow.
% 'mg' is not a conventional regularization, but rather uses non-linear
% regression to fit a spectrum of relaxation time constants with a number of 
% guassian-shaped components in the log-time-constant domain. Because it is
% a non-linear fitting, initial estimates for time-constants can be supplied
% (fitting.T0guass) or will be generated automatically. Also, to avoid
% local minima solutions, the fitting can be automatically repeated with
% varied intial estimates (fitting.numberT0) and the lowest error solution is
% selected
fitting.regtyp = 'me';
% % optional inputs for 'mg' fitting
% fitting.widthgauss = 5;
% fitting.numbergauss = 2;
% fitting.T0gauss = [0.02 0.07]';
% fitting.numberT0 = 5;

% for conventional regularization, the regulization weighting must be
% adjusted. This can be done manually (fitting.regadj = 'manual') or using
% one of two automatic methods ('gcv' and 'erinc') -- see literature
% references in code
fitting.regadj = 'manual';
fitting.regweight = 0.001;

% graph the results or not. This input is irrelevant if
% analysis.interactive = 'y'

analysis.graph = 'n';

% define the range of time constants to fit. Note that for 'mg' fitting,
% this is the full range of the spectrum, but the lowest and highest mean
% component time constants (echoed to display) cover a narrower domain
fitting.rangeT = [2.5e-3 0.25];

% define the number of time constants in the spectrum
fitting.numberT = 200;

% set the non-negative least square (NNLS) code. In most cases, the supplied
% nnlsmex fuction is the fastest NNLS solver. You may need to compile this
% for your system. Or you can specify the MATLAB fuction,
% fitting.nnlscode = 'lsqnonneg'. 

% You can automatically or manually extract a finite number of component
% amplitudes and time constants. Test this out with the interactive GUI.
analysis.extract = 'auto';
analysis.numberextract = 2;

% If automatic component extraction is used, you can also automatically 
% disregard signals outside a specified T-domain, analysis.rangeA, or those 
% below a specified signal fraction, analysis.tolextract

% analysis.rangeA = [10e-3,1];
% analysis.tolextract = 2;


data = maketestdata1D(T2,f,t,theta,Ntrials,SNR);

tic
[out1D,fittingout]=MERA(data,fitting,analysis);
toc






%% fit with a variety of regularization options


Ntrials = 100; % number of independent noise realizations
SNR = 250; % signal to noise ratio at time t = 0
T2 = [20e-3, 100e-3]';
f = [2 8]';
theta = 160;

data = maketestdata1D(T2,f,t,theta,Ntrials,SNR);


% loop through all regulatization types and fit data 
regtyp_list = {'none','me','mc','upen','mg'};
analysis.interactive = 'n';
fitting.B1fit = 'y';
analysis.graph = 'n';
fitting.numbergauss = 2;
fitting.widthgauss = 5; 
runtime = zeros(1,length(regtyp_list));
close all
clear FvM FvS TvM TvS AvM

for kr = 1:length(regtyp_list)
   clear Av Fv Tv
   fitting.regtyp = regtyp_list{kr};
   tic
   [out1D] = MERA(data,fitting,analysis); 
   runtime(kr) = toc;
   % calculate mean and SD of signal component fractions and T2      
   AvM{kr} = mean(out1D.Av,2);
   FvM{kr} = mean(out1D.Fv,2);
   FvS{kr} = std(out1D.Fv,[],2);
   TvM{kr} = mean(out1D.Tv,2);
   TvS{kr} = std(out1D.Tv,[],2);
   
   outtest{kr} = out1D;
end

% Test results with SNR = 250, Nt = 100
% T2 = [20e-3, 100e-3]';
% f = [2 8]';
% means and coef of variation, each column is a different 'regtyp'
% 
% show benchmarks for each regtyp
%     'none'      'me'      'mc'    'mg'      'upen'
% Mac Pro, Early 2008, 2.8 GHz quad core Intel Xenon
%   %bi-exponential, theta = 180, B1fit = 'n'
%    0.5380    1.7323    1.6742    3.7294    9.9274
%   %bi-exponential, theta = 160, B1fit = 'y'
%    16.5054   30.4711   30.8034   54.6984  143.4539
% iMac 21.5 inch, mid 2011, 2.8 GHz Intel quad core i7
%   %bi-exponential, theta = 180, B1fit = 'n'
%    0.5068    1.6009    1.5481    4.4234    8.5202
%   %bi-exponential, theta = 160, B1fit = 'y'
%    12.4144   25.9889   24.8830   45.8386  122.6775

