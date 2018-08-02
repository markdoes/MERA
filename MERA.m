function [output,fitting,analysis] = MERA(data,fitting,analysis)
% function [output] = MERA(DATA,FITTING,ANALYSIS)
%
% MERA version 2.06
% fits data in structure DATA with a distribution of decaying exponential
% functions, subject to a non-negative constraint and other options as
% defined in the structure FITTING. Output includes analysis of
% results as defined by options in the structure ANALYSIS.
%
% Required Inputs:
%
%   DATA.D: column vector of signal magnitudes. If DATA.D is a matrix, then
%     each colum is fitted as a separate signal unless FITTING.twoD = 'y'
%     If FITTING.twoD = 'y', DATA.D must be a matrix of dimensions
%     [length(DATA.t),length(DATA.t2)]
%
%   DATA.t: column vector of sample times (units of seconds)
%
%   DATA.t2: column vector of sample times (unit of seconds) for the second
%     dminension of DATA.D when FITTING.twoD = 'y'
%
% Optional Inputs: (defaults)
%
%   FITTING.regtyp: type of regularization:
%         'none' = no regularization
%         'me'   = minimum energy
%         ('mc') = minimum curvature
%         'upen' = uniform penalty
%         'mg'   = multiple Gauss shaped components
%   If FITTING.regtyp = 'mg'
%     FITTING.numbergauss: number of gauss-shaped components (2)
%     FITTING.widthgauss: width of gauss-shaped components (10), acceptable
%     range is 0.4 - 20
%     FITTING.T0gauss: initial guesses of T2 values (uniformly spaced)
%       Note that unlike the other options, 'mg' fitting is a form of
%       non-linear regression, so results may be sensitive to initital
%       guesses of the T2 values.
%     FITTING.numberT0: defines the number of varied initial guesses to
%     test. The guesses are randomly varied about the values T0guess (1)
%
%   FITTING.regadj: method of regularizer adjustment:
%           'manual' = regularizer weighted by input FITTING.mux
%           ('gcv')  = generalized cross validation
%           'erinc'  = regularized weighted to increase resnom by
%                       FITTING.pinc percent
%   If FITTING.regadj = 'manual'
%     FITTING.regweight: regularization weighting (0.1)
%       note that this is different from "mu" in Version 1.0 by a square
%       root
%   If FITTING.regadj = 'erinc'
%     FITTING.percentinc: regularization weighting increased until
%       the fit residual norm is increased by FITTING.percentinc (1.0)
%   If FITTING.regtyp = 'mg', FITTING.regadj is set to 'none'
%
%   FITTING.T: vector of exponential time constants (units of seconds)
%     (FITTING.numberT logspaced values spanning FITTING.rangeT)
%   if FITTING.T is not presecribed
%     FITTING.rangeT: range of exponential time constants
%       ([3*min(DATA.t)/4,max(DATA.t)*2])
%     FITTING.numberT: number of exponential time constants (100)
%     if FITTING.twoD = 'y'
%       FITTING.rangeT2: range of 2nd dimension exponential time constants
%       FITTING.numberT2: number of 2nd dimension exponential time constants
%
%   FITTING.B1fit: flag to turn on B1 fitting ('y') or 'n'
%   if FITTING.B1fit = 'y'
%     FITTING.rangetheta: range of flip angles to test ([130 180])
%     FITTING.numbertheta: number of flip angles to test (10)
%     FITTING.fixedT1: the value to fix T1 for the EPG analysis (1)
%     ** B1 fitting is not configured for two-dimensional fitting
%     Note this option is only relevant to specific NMR measurements -- see
%     references to background materials below.
%
%   FITTING.twoD: flag to signify two dimensional analysis 'y' or ('n')
%
%   ANALYSIS.graph: flag to turn on graphing of results 'y' or ('n')
%
%   ANALYSIS.interactive: flag to open GUI for interactive fitting 'y' or
%   ('n')
%
%   ANALYSIS.extract: sets component extraction method as ('auto') or
%   'manual'
%
%   If ANALYSIS.extract = 'auto'
%     ANALYSIS.rangeA: open interval of exponential time constants
%        to include in analysis ([min(FITTING.rangeT) max(FITTING.rangeT)])
%     ANALYSIS.tolextract: minimum signal component amplitude in percent of
%           total used in analysis (2)
%     If FITTING.twoD = 'y'
%       ANAYLSIS.rangeA2: openn interval of 2nd dimension exponential time
%       constants to include in analysis (([min(FITTING.rangeT2)
%       max(FITTING.rangeT2)])
%
%   If ANALYSIS.extract = 'manual'
%     ANALYSIS.numberextract: set number of component (2) to manually
%     extract
%
%   OUTPUT
%   The primary output is written into the structure OUTPUT
%
%   OUTPUT.S: The fitted spectrum (or spectra) of signal amplitudes as a
%   function of time constants
%
%   OUTPUT.T: The time constants associated with colums of OUTPUT.S.
%   OUTPUT.T2: The time constants associated with the rows of OUTPUT.S, if
%     FITTING.twoD = 'y'
%
%   OUTPUT.R: Residuals to the fit(s)
%   OUTPUT.P: Fit-predicted signal amplitudes
%
%   OUTPUT.mu: Weighting of the regularizing term in the cost function, if
%     FITTING.regtyp = 'none'
%
%   OUTPUT.offset: amplitude of the fitted dc-offset term
%
%   OUTPUT.Er: Fit error(s) calculated as the sum of the squared residuals
%   OUTPUT.Er0: Fit error(s) without regularization
%
%   OUTPUT.Tv: Time constants of signal components extracted from the fit(s)
%     If FITTING.twoD = 'y', then OUTPUT.Tv(:,1) contains time constants
%     associated with OUTPUT.T, and OUTPUT.Tv(:,2) contains time constants
%     associated with OUTPUT.T2
%   OUTPUT.Av: Amplitude(s) of the signal components extracted from the fit(s)
%   OUTPUT.Fv: Same values as in OUTPUT.Av, but scaled to be signal
%   fractions
%
%   OUTPUT.SNR: The estimated signal-to-noise ratio of the input data.D,
%   calculated as sum(S)./std(R)
%
%   OUTPUT.theta: Fitted flip angle if FITTING.B1fit = 'y'
%
%   For convenience, two fitting parameters are passed to the output
%   structure
%   OUTPUT.regtyp: The type of regulization used
%   OUTPUT.regadj: The type of regulization adjustment used
%
%   All fitting parameters, updated with any changes made by/during the
%   processing, are output in the FITTING structure. Likewise for the
%   output ANALYSIS stucture.
%
%   This code as been formatted to fit in an 80-column width using 2-space
%   tabs
%% License

% <MERA Toolbox Version 2.02 (1D,2D Multi-exponential Relaxation Analysis)>
% Copyright (C) <2014>  <Mark D. Does, Vanderbilt University>
% MERA is written and maitained by Mark D. Does, with contributions from
% numerous people, listed below.
% Contact: mark.does@vanderbilt.edu

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% See <http://www.gnu.org/licenses/>.
%
% Note that this code has been tested using MATLAB R2013b. For full
% functionality, this code makes use of three MATLAB product toolboxes:
% statistics, optimization, and image processing.
%
% Please acknowledge use of this product by references the website:
% https://github.com/markdoes/MERA
% formerly
% http://www.vuiis.vanderbilt.edu/~doesmd/MERA/MERA_Toolbox.html
%
% MERA was developed with financial support from the NSF 0448915, and
% NIH EB001744 and EB019980
%
% If you have a burning desire to thank me for developing, maintaining, and
% sharing this code, consider that I am a big fan of single malt whisky!

%% Version History
%
% Version 2.0, 15 Feb 2014
%
% Version 2.01, 26 Feb 2014
%   corrects auto-extract error that may have cased a short time-constant
%   component to be missed or inaccurately calculated
%
%   forces analysis.graph = 'y' if analysis.extract = 'manual'
%
% Version 2.02, 11 Mar 2014
%   corrects auto-extract error in 2D fitting -- Thanks to Alan Seifert at
%   UPENN for finding this problem
%
% Version 2.03, 11 Sept 2014
%   corrected the help comments to clarify the structure name for input data to
%   be "Data.D" not "Data.d"
%
% Version 2.04, 13 Oct 2017
%   modified the B1_fitting to use uncontrainted fits for finding the
%   optimal theta, then the selected regularization only for the final fit.
%   This change is for efficiency only and shouldn't alter results.
%
% Version 2.05, 29 July 2018
%   finally because the website that hosts MERA is going through a required
%   address change, I figured it was time to move it to github. For this, a
%   new version number!
%
% Version 2.06, 1 Aug 2018
%   As an exercise to learn how to make code changes through github, I've
%   added the option to constrain T1 when B1fit = 'y'. In cases of contrast
%   agent loaded tissue samples, this might make a small difference in the
%   fitted T2 spectrum, but in the vast majority of cases, it's unimportant
%   see parameter fitting.fixedT1
%
%   Also, %done display changed to steps of 10%

%% Acknowledgements/Contributions

% This code was developed with support from the NSF (grant #0448915), and
% the NIH (grants EB001744 and EB019980)


% Richard D. Dortch provided the MEX NNLS code, and wrote the original
%   autopeaks.m function (some aspects of which remain in the current code)
%   The original version of the NNLS code was developed by
%   Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
%   1973 JUN 15, and published in the book
%   "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
%   Revised FEB 1995 to accompany reprinting of the book by SIAM.

% Meher Juttukonda, developed the initial beta versions of MERA v1.0 from a
%   collection of old code of mine, je developed the first verion of the MERA
%   GUI in MERA v1.0, he wrote the v1.0 Documentation.docx file, and wrote
%   the html and php files for the first MERA website. Although little of
%   the code of other materials Meher wrote still exist in MERA v2.0, I
%   would still be working from dozens of loosely related bits of code if
%   it weren't for his efforts in getting MERA off the ground.

% Kathryn L. West has helped test and debug numerous developmental
%   versions of MERA 2.0

% Thomas Prasloski, shared code for B1 correction at the Iceland ISMRM
%   white matter study group workshop. That code and the related paper,
%   Prasloski, T., Mädler, B., Xiang, Q.-S., MacKay, A., & Jones, C. (2012)
%   Applications of stimulated echo correction to multicomponent T2 analysis.
%   Magnetic Resonance In Medicine, 67(6), 1803-1814. doi:10.1002/mrm.23157
%   was used to guide the B1-correction code added to MERA 2.0

% Rather than use the MATLAB function 'kron' in 2D processing, I used code
% from the function 'kronecker' available from MATLAB Central File Exchange.
%   Author: Bruno Luong <brunoluong@yahoo.com>
% This is also how I originally discovered the bsxfun function in MATLAB!

% I have also received helpful feedback from a number of users. Please
% continue to send me bug reports and requests new functionality. I promise
% not to wait so long between updates next time!

%% Literature References
%
% description of standard non-negative least-squares fitting of
% relaxation-time spectra (or T2 spectra)
% without smoothing (regtyp = 'none')
% with minimum amplitude energy regularization (regtyp = 'me')
% or with minimum amplitude curvature regularization (regtyp = 'mc')
%
% Whittall, K. P., & MacKay, A. L. (1989).
%  Quantitative Interpretation of NMR Relaxation Data.
%  Journal Of Magnetic Resonance, 84, 134-152.
%
% Some earlier works worth noting are
%
% Provencher, S. W. (1982a). A constrained regularization method for
%  inverting data represented by linear algebraic or integral equations.
%  Computer Physics Communications, 27(3), 213-227.
% Provencher, S. W. (1982b). CONTIN: a general purpose constrained
%  regularization program for inverting noisy linear algebraic and integral
%  equations. Computer Physics Communications, 27(3), 229-242.
% Kroeker, R. M., & Mark Henkelman, R. (1986). Analysis of biological NMR
%  relaxation data with continuous distributions of relaxation times.
%  Journal of Magnetic Resonance (1969), 69(2), 218-235.
%
% Provencher's well-known fortran program, CONTIN, is available at
%  http://s-provencher.com/pages/contin.shtml
%
%
% The Uniform Penalty method of regularizing the solution, (regtyp =
% 'upen')comes from
%
% Borgia, G. C., Brown, R. J. S., & Fantazzini, P. (1998).
%  Uniform-Penalty Inversion of Multiexponential Decay Data.
%  Journal Of Magnetic Resonance, 132(1), 65-77.
%
% and
%
% Borgia, G. C., Brown, R. J. S., & Fantazzini, P. (2000).
%  Uniform-Penalty Inversion of Multiexponential Decay Data.
%  Journal Of Magnetic Resonance, 147(2), 273-285.
%
%
% Fitting multiple Gaussian-shaped spectral components (regtyp = 'mg')
% is described in
%
% Stanisz, G., & Henkelman, R. (1998). Diffusional anisotropy of T-2
%  components in bovine optic nerve. Magnetic Resonance In Medicine,
%  40(3), 405-410.
%
% and later used and described in
%
% Dortch, R. D., Harkins, K. D., Juttukonda, M. R., Gore, J. C.,
%  & Does, M. D. (2013). Characterizing inter-compartmental water exchange
%  in myelinated tissue using relaxation exchange spectroscopy.
%  Magnetic Resonance In Medicine, 70(5), 1450?1459 doi:10.1002/mrm.24571
%
%
% The Generalized Cross Validation approach to adjusting the regulizer
% weighting, regadj = 'gcv', is derived and presented in
%
% Golub, G. H., Heath, M., & Wahba, G. (1979). Generalized cross-validation
%  as a method for choosing a good ridge parameter. Technometrics, 215-223.
%
%
% Adjusting the regularer weight empircally until the mean-square
% error is increased by a given percentage from the case with no
% regularization, regadj = 'erinc'
% is presented in
%
% Graham, S. J., Stanchev, P. L., & Bronskill, M. J. (2011).
%  Criteria for analysis of multicomponent tissue T2 relaxation data.
%  Magnetic Resonance In Medicine, 35(3), 370-378.
%
%
% Fitting the refocusing pulse flip angle along with the T2 spectrum is
% presented in
%
% Prasloski, T., Mädler, B., Xiang, Q.-S., MacKay, A., & Jones, C. (2012)
%   Applications of stimulated echo correction to multicomponent T2 analysis.
%   Magnetic Resonance In Medicine, 67(6)
%
% Other useful papers related to this topic are
%
% Lebel, R. M., & Wilman, A. H. (2010). Transverse Relaxometry with
%  Stimulated Echo Compensation. Magnetic Resonance In Medicine, 64(4),
%  1005-1014.
%
% and
%
% Hennig, J. (1988). Multiecho imaging sequences with low refocusing
%  flip angles. Journal of Magnetic Resonance (1969), 78(3), 397-407.
% Hennig, J. (1991). Echoes-how to generate, recognize, use or avoid them
%  in MR-imaging sequences. Part I: Fundamental and not so fundamental
%  properties of spin echoes. Concepts in Magnetic Resonance, 3(3), 125-143.
%
%
% The singular-value-decomposition of the linear problems, used in MERA for
% the 2D problems is described in
%
% Venkataramanan, L., Yi-Qiao Song, & Hürlimann, M. D. (2002).
%  Solving Fredholm integrals of the first kind with tensor product
%  structure in 2 and 2.5 dimensions. IEEE Transactions on Signal
%  Processing, 50(5), 1017-1026.

%% EXAMPLE CODE for processing image data
%
% given object IMG which is Nx x Ny x NE in size
% and object TE which is NE x 1
% In order to fit each (or some) point in the image (Nx,Ny) to a spectrum
% of relaxation time constants, this 3D data must be reformatted into a 2D
% obect, Data.D
%
% first, identify which image points that are of interest to analyze. For
% example, choose all points about some threshold amplitude
%  THRESHOLD = 1000; % FILL IN APPROPRIATE VALUE HERE
%  bw = roicolor(IMG(:,:,1),THRESHOLD,1e12);
% Nvox = length(find(bw));
% NE = length(TE);
%
% IMGr = reshape(IMG,Nx*Ny,NE);
% bwr = reshape(bw,Nx*Ny,1);
% X = IMGr(bwr,:)'; % NOW X is an NE x Nvox matrix of decay data
%
% data.D = X;
% data.t = TE;
%
% tic
% [out1D]=MERAv2beta(data,fitting,analysis);
% toc
%
% Extract some of the out1D data back into image space
%
% The fitted signal amplitude at t = 0
% M0 = zeros(Nx*Ny,1);
% M0(bwr) = sum(out1D.S);
% M0 = reshape(M0,Ny,Nx);
%
% The fraction of the signal over some domain of relaxation time constants
% MWF = zeros(Nx*Ny,1);
% MWF(bwr) = sum(out1D.S(Tmin:Tmax,:))./sum(out1D.S);
% MWF = reshape(MWF,Ny,Nx);

%% main code
% warning off % avoids warnings from various non-linear regression steps

% Interactive or non-interactive analysis
if exist('analysis','var')
  if ~isfield(analysis,'interactive');
    analysis.interactive = 'n';
  end
else
  analysis.interactive = 'n';
end

if ~exist('fitting','var')
  fitting = [];
end

% call input argument parsing
[data,fitting,analysis] = parseinput(data,fitting,analysis);

% start the GUI or the loop
if strcmpi(analysis.interactive,'y')
  [output,fitting,analysis] = MERA_GUI(data,fitting,analysis);
else
  [output,fitting,analysis] = fittingloop(data,fitting,analysis);
end
end

function [output,fitting,analysis] = fittingloop(data,fitting,analysis)

%% create interactive window or waitbar

if strcmpi(analysis.interactive,'y')
  cla(analysis.htdomain);
  axes(analysis.hTdomain)
  cla;
  w = text(0.25,0.75,'Calculating ...');
  drawnow;
else
  fprintf('%3.0f %% done...',0);
end
%% Extract data and normalize data

if strcmpi(fitting.twoD,'y')
  mD = max(abs(data.D(:)));
  Dn = data.D/mD;
  t = data.t;
  t2 = data.t2;
else
  mD = max(abs(data.D));
  Dn = bsxfun(@rdivide, data.D, mD);
  t = data.t;
  t2 = [];
end

%% Preallocate memory and build large objects

% some shorthand variables
K1 = fitting.numberT+1;
K2 = fitting.numberT2+1;
K = K1*K2;

N1 = fitting.numberechoes;
N2 = fitting.numberechoes2;
N = N1*N2;

% build design matrix for d = As linear model, where d is a column vector
% of signal observations and s is a column vector of signal amplitudees for
% all relaxation time constants.build A matrix only once, reuse for all
% data vectors

[Amat_all] = buildA(fitting.theta_vector);

% build regularization matrix which will augment the design matrix, such
% that the linear model becomes [d;0] = [A;H]s

switch fitting.regtyp
  case 'none'
    H = [];
  case 'me' % Minimize amplitude energy
    H = speye(K,K); H(end,end) = 0;
  case {'mc','upen'} % Min curvature energy
    if strcmpi(fitting.twoD,'y')
      Ha = diag(ones(K1,1)); Ha([1,end])=0;
      Hb = diag(ones(K1-1,1),1)+diag(ones(K1-1,1),-1) ...
        -4*diag(ones(K1,1)); Hb([1,end],:) = 0;
      Hr = [Ha Hb Ha zeros(K1,K-3*K1)];
      H = spalloc(K,K,5*(K1-3)*(K2-3));
      for kH = 1:(K2-2)
        H(kH*K1+1:(kH+1)*K1,:) = circshift(Hr,[0,(kH-1)*K1]);
      end
    else
      H = -2*eye(K1,K1)+diag(ones(K1-1,1),1) ...
        +diag(ones(K1-1,1),-1);
      H(1,:) = 0; H(end,:) = 0;
    end
  case {'mg'}
    H = [];
end
fitting.H = H;

% pre-allocate space for various parameters

S = zeros(K,fitting.numbertrains);
R = zeros(N,fitting.numbertrains);
P = zeros(N,fitting.numbertrains);
muxm = zeros(1,fitting.numbertrains);
Er = zeros(1,fitting.numbertrains);
Er0 = zeros(1,fitting.numbertrains);
offset = zeros(1,fitting.numbertrains);
Tv = []; Fv = []; Av = [];

theta_min = 180*ones(1,fitting.numbertrains);

if strcmpi(fitting.regtyp,'mg')
  
  % short hand variables
  ng = fitting.numbergauss;
  u = log(fitting.T);
  sig = range(u)/100*fitting.widthgauss/2;
  lb_fitmG = exp(u(1)+5*sig)*ones(ng,1);
  ub_fitmG = exp(u(end)-5*sig)*ones(ng,1);
  
  if strcmpi(fitting.twoD,'y')
    u2 = log(fitting.T2);
    sig2 = range(u2)/100*fitting.widthgauss/2;
    Tv = zeros(fitting.numbergauss,2,fitting.numbertrains);
    TvT0 = zeros(fitting.numbergauss,2,fitting.numberT0);
    lb_fitmG(:,2) = exp(u2(1)+5*sig2)*ones(ng,1);
    ub_fitmG(:,2) = exp(u2(end)-5*sig2)*ones(ng,1);
    Ttxt = 'Ta/Tb';
  else
    Tv = zeros(fitting.numbergauss,1,fitting.numbertrains);
    TvT0 = zeros(fitting.numbergauss,1,fitting.numberT0);
    Ttxt = 'T';
  end
  display(['lowest Gauss component mean ', Ttxt,' = ', ...
    num2str(lb_fitmG(1,:)*1e3,3),' ms'])
  display(['highest Gauss component mean ', Ttxt,' = ', ...
    num2str(ub_fitmG(1,:)*1e3,3),' ms'])
  
  Fv = zeros(fitting.numbergauss,fitting.numbertrains);
  Av = zeros(fitting.numbergauss,fitting.numbertrains);
  ErT0 = zeros(1,fitting.numberT0);
  fopmg = optimset(@fmincon);
  fopmg = optimset(fopmg,'TolX',1e-8,'TolFun',1e-8,'Display','off',...
    'Algorithm','active-set');
end

% optimization for SpectrumCalc
fopsILT = optimset(@fminbnd);
fopsILT = optimset(fopsILT,'TolX',1e-8,'Display','off');

%% loop through every vector of decay data

for i = 1:fitting.numbertrains
  % extract a column vector of data
  if strcmpi(fitting.twoD,'y')
    Dx = Dn(:,:,i);
  else
    Dx = Dn(:,i);
  end
  
%   % update waitbar
%   if strcmpi(analysis.interactive,'n')
%     waitbar(i/fitting.numbertrains,w,'Calculating ...');
%   end
  % update waitbar
  if strcmpi(analysis.interactive,'n')&&~mod(i,fitting.numbertrains/10)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f %% done...',...
      i/fitting.numbertrains*100);
  end
  
  % initialize fit errors to zero for all possible theta
  Erb = zeros(fitting.numbertheta,1);
  Er0b = zeros(fitting.numbertheta,1);
  
  % spectral fitting at different refocussing flip angles
  for kb = 1:fitting.numbertheta
    % MDD 2.04 update, fit B1 without regularization for speed
    if strcmpi(fitting.B1fit,'y')&& kb==1
      Htemp = fitting.H;
      regAtemp = fitting.regtyp;
      fitting.H = [];
      fitting.regtyp = 'none';
    end
    
    Amat = Amat_all{kb};
    Dr = Amat.Ur'*Dx*Amat.Ur2;
    
    if strcmpi(fitting.regtyp,'mg')
      % for regtype=mg, Multi-Gauss fitting
      for k0 = 1:fitting.numberT0
        T0 = fitting.T0gauss(:,:,k0);
        [TvT0(:,:,k0),ErT0(k0)] = fmincon(@fitmG,T0,[],[],[],[], ...
          lb_fitmG,ub_fitmG,[],fopmg);
      end
      [~,i0] = min(ErT0); Tv(:,:,i) = TvT0(:,:,i0);
      T0min = fitting.T0gauss(:,:,i0);
      [Erb(kb),R(:,i),S(:,i),P(:,i),Av(:,i),Fv(:,i),offset(i)] =  ...
        fitmG(Tv(:,:,i));
    else
      % use SpectrumCalc for regtype =  me, mc, upen, none
      [S(:,i),R(:,i),P(:,i),muxm(i),Erb(kb),Er0b(kb),offset(i)] = ...
        SpectrumCalc(Dx,Dr,Amat,fitting,fopsILT);
    end
    
  end
  
  if strcmpi(fitting.B1fit,'y')
    % MDD 2.04 update, reset regularizer
    fitting.H = Htemp;
    fitting.regtyp = regAtemp;
    ErInterp=interp1(fitting.theta_vector,Erb, ...
      fitting.theta_vector_fine,'spline');
    [~,imn] = min(ErInterp);
    theta_min(i) = fitting.theta_vector_fine(imn);
    fitting.theta = theta_min(i);
    Amat = cell2mat(buildA(theta_min(i)));
    Dr = Amat.Ur'*Dx; %*Ur2; %Ur2 is Identity for 1D fitting
    
    % Now re-do the spectral fit at the fitted theta value
    if strcmpi(fitting.regtyp,'mg') % mg fitting
      % Solve for component central T2 by nonlinear regression
      Tv(:,:,i) = fmincon(@fitmG,T0min,[],[],[],[], ...
        lb_fitmG,ub_fitmG,[],fopmg);
      [Er(i),R(:,i),S(:,i),P(:,i),Av(:,i),Fv(:,i),offset(i)] = ...
        fitmG(Tv(:,:,i));
    else
      [S(:,i),R(:,i),P(:,i),muxm(i),Er(i),Er0(i),offset(i)] = ...
        SpectrumCalc(Dx,Dr,Amat,fitting,fopsILT);

    end
    
  else
    
    Er(i) = Erb(kb);
    Er0(i) = Er0b(kb);
  end
end

%renormalize signal amplitudes

S = bsxfun(@times,S,mD);
P = bsxfun(@times,P,mD);
R = bsxfun(@times,R,mD);
offset = bsxfun(@times,offset,mD);
Er = bsxfun(@times,Er,mD.^2);
Er0 = bsxfun(@times,Er0,mD.^2);

if strcmpi(fitting.regtyp,'mg')
  Av = bsxfun(@times,Av,mD);
  if ~strcmpi(fitting.twoD,'y')
    Tv = reshape(Tv,ng,fitting.numbertrains);
  end
end

SNR= sum(S)./std(R);

% if strcmpi(analysis.interactive,'n')
%   close(w)
% end

% reshape 2D fitting into matrices
if strcmpi(fitting.twoD,'y')
  S = reshape(S,K1,K2,fitting.numbertrains);
  S = S(1:end-1,1:end-1,:);
  P = reshape(P,N1,N2,fitting.numbertrains);
  R = reshape(R,N1,N2,fitting.numbertrains);
else
  S = S(1:end-1,:);
end

fitting.regweight = muxm;

[Av,Fv,Tv,analysis] = peakprocess(S,P,R,Fv,Tv,Av,Er,data,fitting,analysis);

% format output
output = struct('S',S(1:fitting.numberT,:),'T',fitting.T, ...
  'R',R,'P',P,'mu',muxm,'offset',offset,'Er',Er,'Er0',Er0,'Tv',Tv, ...
  'Av',Av,'Fv',Fv,'SNR',SNR,'theta',theta_min, ...
  'regtype', fitting.regtyp, 'regadj',fitting.regadj);
if strcmpi(fitting.twoD,'y')
  output.T2 = fitting.T2;
  output = orderfields(output,[1:2,16,3:15]);
end

%% sub functions

  function [Amat] = buildA(theta_vector)
   % For 1D problems, build Amatrix and associated objects for every theta
   % value in advance
   % 2D problems are not compatible with B1fitting, so theta_vector = 180
    
    numbertheta = length(theta_vector);
    Amat = cell(numbertheta,1);
    TE = t(2)-t(1);
    A = zeros(N1,K1);
    
    for k = 1:numbertheta
      th = theta_vector(k);
      % build A-matrix for a given theta
      if th == 180
        A(:,1:(K1-1)) = exp(-t*(1./fitting.T'));
      else
        for kt = 1:(K1-1)
          A(:,kt)=EPGdecaycurve(th,fitting.T(kt),TE,N1,fitting.fixedT1);
        end
      end
      A(:,K1) = ones(N1,1);
      % SVD dimension reduction
      [U,Sig,V] = svd(A,'econ');
      s = rank(Sig);
      Amat{k}.A = A;
      Amat{k}.Ur = U(:,1:s);
      Amat{k}.Ar = Sig(1:s,1:s)*V(:,1:s)';
    end
    
    if strcmpi(fitting.twoD,'y')
      A2 = zeros(N2,K2);
      A2(:,1:(K2-1)) = exp(-t2*(1./fitting.T2'));
      A2(:,K2) = ones(N2,1);
      [U,Sig,V] = svd(A2,'econ');
      s2 = rank(Sig);
      Amat{1}.A2 = A2;
      Amat{1}.Ur2 = U(:,1:s2);
      Ar2 = Sig(1:s2,1:s2)*V(:,1:s2)';
      
      % compute the kron tensor product of Ar and Ar2
      Ar2 = reshape(Ar2,[1 s2 1 K2]);
      Amat{1}.Ar = reshape(Amat{1}.Ar,[s 1 K1 1]);
      Amat{1}.Ar = reshape(bsxfun(@times,Ar2,Amat{1}.Ar),[s*s2 K]);
    else
      for k = 1:numbertheta
        Amat{k}.A2 = 1;
        Amat{k}.Ur2 = 1;
      end
    end
    
    function EchoAmp = EPGdecaycurve(theta,T2,TE,N,T1)
      % Assumes a CPMG condition -- see Hennig refs above
      % Arbitrarily set T1 = 1 s
      % T1 = 1;
      % version 2.06, added fitting.fixedT1, so T1 now set outside this func
      Nx = 2;
      Np = floor(N/Nx);
      
      % Compute rotation matrix in terms of coherence state
      ra = cosd(theta/2)^2; rb = sind(theta/2)^2;
      rc = -1i*sind(theta);rd = -0.5i*sind(theta);
      R1 = [ra;rb;rd]; R2 = [rb;ra;-rd]; R3 = [rc;-rc;cosd(theta)];
      Rot = [R1,R2,R3];
      
      % Initialize vector to track echo amplitude
      EchoAmp=zeros(N,1);
      
      % Initialize magnetization phase state vector (MPSV) and set all
      % magnetization in the F1 state.
      M=zeros(3*Np,1);
      M(1,1)=exp(-(TE/2)/T2);
      
      % Compute relaxation and transition matrix, E
      % same as T*E in Prasloski et al.
      E2 = exp(-TE/T2);
      E22 = exp(-TE/2/T2);
      E1 = exp(-TE/T1);
      ke = (1:Np-1)';
      r = [1;3*ke-1;3*ke+1;3*ke];
      c = [2;3*ke+2;3*ke-2;3*ke];
      m = 3*Np;
      e = [E2;ones(Np-1,1)*E2;ones(Np-1,1)*E2;ones(Np-1,1)*E1];
      RelTrans = sparse(r,c,e,m,m);
      
      r = [(1:m)';(1:m)';(1:m)'];
      c = (1:3:(m-2))'; c = [c,c,c]'; c = c(:);
      c = [c;c+1;c+2];
      p1 = R1(:,ones(Np,1)); p2 = R2(:,ones(Np,1)); p3 = R3(:,ones(Np,1));
      p = [p1(:);p2(:);p3(:)];
      BigRot = sparse(r,c,p,m,m);
      RRT = BigRot*RelTrans;
      
      % Apply first refocusing pulse and get first echo amplitude
      M(1:3)=Rot*M(1:3);
      EchoAmp(1)= M(2,1);
      
      % Apply relaxation matrix
      M=RRT*M;
      
      % Perform flip-relax sequence ETL-1 times
      for ke=2:N
        % Record the magnitude of the population of F1* as the echo amplitude
        EchoAmp(ke)=M(2,1);
        % Allow time evolution of magnetization between pulses
        M=RRT*M;
      end
      EchoAmp = abs(EchoAmp)*E22;
    end
    
  end

  function [Er,R,S,P,Av,Fv,dc] = fitmG(Tm)
    G = zeros(K1,K2,fitting.numbergauss+1);
    G(end,end,end) = 1;
    um = log(Tm(:,1));
    
    for kg = 1:fitting.numbergauss
      gc = 1/sqrt(2*pi*sig^2)*exp(-(u-um(kg)).^2/2/sig^2);
      if strcmpi(fitting.twoD,'y')
        um2 = log(Tm(:,2));
        gc2 = 1/sqrt(2*pi*sig2^2)*exp(-(u2-um2(kg)).^2/2/sig2^2);
        G(1:K1-1,1:K2-1,kg) = gc*gc2';
      else
        G(1:K1-1,:,kg) = gc;
      end
    end
    
    G = reshape(G,K1*K2,fitting.numbergauss+1);
    f = nnlsLH(Amat.Ar*G,Dr(:),fitting.nnlscode);
    P = Amat.A*reshape(G*f,K1,K2)*Amat.A2'; % Calc model
    P = P(:);
    R = P - Dx(:);% Calc residual
    dc = f(end);
    S = G(:,1:fitting.numbergauss)*f(1:fitting.numbergauss);
    Av = sum(G(:,1:fitting.numbergauss))'.*f(1:fitting.numbergauss);
    Fv = Av./sum(Av);
    Er = norm(R).^2;
    
  end

end

function [S,R,P,mux,Er,Er0,dc] = SpectrumCalc(D,Dr,Amat,fitting,fopsILT)
%% main fitting function

H = fitting.H;
mux = fitting.regweight;
dQ = diff(fitting.T(1:2));
Dr = Dr(:);

% unregularized fit: seed for upen, min error for erinc, and final fit for
% regtyp = 'none'

S = nnlsLH(Amat.Ar,Dr,fitting.nnlscode); % spectrum

% P: predicted t-domain signal
P = Amat.A*reshape(S,size(Amat.A,2),size(Amat.A2,2))*Amat.A2';
P = P(:);
R = P-D(:); % residuals
Er = norm(R).^2; % sum square error
Er0 = Er; % sum square error of unregularized fit
dc = S(end); % dc-offset term

if strcmpi(fitting.regtyp,'none')
  mux = 0;
else
  Dm = [Dr;zeros(size(H,1),1)]; % any regularization needs this
  
  % autocalculating of mux, if regadj = 'gcv', or 'erinc'
  if ~strcmpi(fitting.regadj,'manual')
    [mux,~,exflg] = fminbnd(@smoothILT,0,100,fopsILT);
  end
  
  % calculate spectrum for set value of mux
  if strcmpi(fitting.regtyp,'upen') % using upen
    [S] = upen_iter(mux);
  else % using uniform me or mc
    Am = [Amat.Ar;sqrt(mux)*H];
    S = nnlsLH(Am,Dm,fitting.nnlscode);
  end
  % recalculate P, R, and Er
  P = Amat.A*reshape(S,size(Amat.A,2),size(Amat.A2,2))*Amat.A2';
  P = P(:);
  R = P-D(:); % residuals
  Er = norm(R).^2; % sum square error
  
  dc = S(end);
  
end
%% Sub Functions

  function er = smoothILT(lam)
    if strcmpi(fitting.regtyp,'upen')
      [~,Rx,Hx] = upen_iter(lam);
    else
      Am = [Amat.Ar;sqrt(lam)*H];
      Sx = nnlsLH(Am,Dm,fitting.nnlscode);
      Px = Amat.A*reshape(Sx,size(Amat.A,2),size(Amat.A2,2))*Amat.A2';
      Px = Px(:);
      Rx = Px-D(:);
      Hx = H;
    end
    
    if strcmpi(fitting.regadj,'gcv')
      n = size(Amat.Ar,1);
      denm = eye(n)-Amat.Ar*((Amat.Ar'*Amat.Ar+lam.*(Hx'*Hx))\Amat.Ar');
      den = trace(denm)^2/n^2;
      er = norm(Rx,2)^2./n./den;
    else
      Erx = norm(Rx,2).^2;
      er = abs((Erx-Er0)./Er0*100-fitting.percentinc);
    end
  end

  function [Sq,Rq,Hm] = upen_iter(lam)
    betaL = ([1e-5,.6,0.3]); % for now, hard code "compliance" factors
    beta0 = betaL(1); betaP = betaL(2); betaC = betaL(3);
    Sq = S; Rq = R;
    for k = 1:5
      p = [0;(Sq(3:end)-Sq(1:end-2)).^2;0];
      c = [0;diff(Sq,2).^2;0];
      cm = max([c,[0;c(1:end-1)],[c(2:end);0]],[],2);
      ck = (betaP*dQ^2*p+betaC*cm);
      Ck = sqrt(sum(cm)./(beta0+ck));
      Hm = H.*(Ck*ones(1,length(Ck)));
      Am = [Amat.Ar;sqrt(lam)*Hm];
      Sq = nnlsLH(Am,Dm,fitting.nnlscode);
      Pq = Amat.A*reshape(Sq,size(Amat.A,2),size(Amat.A2,2))*Amat.A2';
      Rq = Pq(:)-D(:);
    end
  end

end

function [Av,Fv,Tv,analysis] = peakprocess(S,P,R,Fv,Tv,Av,Er,...
  data,fitting,analysis)
FS = 14;

% plot spectrum and time domain data
if strcmpi(analysis.graph, 'y')
  [hTfig,hTax] = plotdata;
end

if ~strcmpi(fitting.regtyp,'mg')
  if strcmpi(fitting.twoD,'y')
    [Av,Fv,Tv,Nc] = peaks2D;
  else
    [Av,Fv,Tv,Nc] = peaks1D;
  end
else
  Nc = fitting.numbergauss;
end

if strcmpi(analysis.extract,'auto')
  % Remove spurious peaks w/ weigh less than tol
  indKe = Fv > analysis.tolextract/1e2;
  Av = Av.*indKe;
  Tv = bsxfun(@times,Tv,indKe);
  % Remove spurious peaks outside acceptable range
  if strcmpi(fitting.twoD,'y')
    indKe = (Tv(:,1) >= analysis.rangeA(1)) & ...
      (Tv(:,1) <= analysis.rangeA(2));
    indKe(:,2) = (Tv(:,2) >= analysis.rangeA2(1)) & ...
      (Tv(:,2) <= analysis.rangeA2(2));
    indKe = indKe(:,1)&indKe(:,2);
  else
    indKe = (Tv >= analysis.rangeA(1)) & (Tv <= analysis.rangeA(2));
  end
  
  Av = bsxfun(@times,Av,indKe);
  Tv = bsxfun(@times,Tv,indKe);
  Fv = bsxfun(@rdivide,Av,sum(Av,1));
end

% sort Fv and Av by ascending Tv values
[~,ix] = sort(Tv,1);
for kj = 1:fitting.numbertrains
  if strcmpi(fitting.twoD,'y')
    Tv = Tv(ix(:,1),:);
  else
    Tv(:,kj) = Tv(ix(:,kj),kj);
  end
  Fv(:,kj) = Fv(ix(:,kj),kj);
  Av(:,kj) = Av(ix(:,kj),kj);
end

% remove empty component
z = sum(Tv,2)>0;
Tv = Tv(z,:);Fv = Fv(z,:);Av = Av(z,:);
Nc = size(Fv,1);

% label plots
if strcmpi(analysis.graph, 'y')
  plotlabel
end


  function [hTfig,hTax] = plotdata
    
    % plot t-domain signals
    if isempty(analysis.htdomain)
      htfig = figure;
      htPosition = get(htfig,'Position');
      htax = gca;
    else
      htax = analysis.htdomain;cla(htax)
      htfig = gcf;
    end
    figure(htfig), axes(htax)
    
    if strcmpi(fitting.twoD,'y')
      [t2m,tm] = meshgrid(data.t2,data.t);
      plot3(t2m,tm,data.D,'LineStyle','-'); grid on;
      view([90-37.5,30]),shading interp
      xlabel('t_b (s)','FontSize',FS,'Units','Normalized', ...
        'Position',[.2,0.025 0.0]);
      ylabel('t_a (s)','FontSize',FS,'Units','Normalized',...
        'Position',[.8,0.025 0.0])
      zlabel('Amplitude (a.u.)','FontSize',FS);
      title('t-domain Signal','FontSize',FS+4,'FontWeight','Bold')
      if fitting.numbertrains==1
        text(0.7,1,0,['RMS fit error = ',num2str(Er,3)],'FontSize',FS,...
          'Units','Normalized');
      end
    else
      hpD = plot(data.t,data.D,'.r','MarkerSize',8); hold on
      hpP = plot(data.t,P,'LineWidth',1.5); hold on
      hpR = plot(data.t,R,'--','LineWidth',1); hold off
      xlabel('time (s)','FontSize',FS)
      ylabel('Amplitude (a.u.)','FontSize',FS)
      legend([hpD(1),hpP(1),hpR(1)],' data ',' fit ',' residuals ', ...
        'location','East')
      title('t-domain Signal, Fit, and Residuals', ...
        'FontSize',FS+4,'FontWeight','Bold');
      if fitting.numbertrains==1
        text(0.7,.9,0,['RMS fit error = ',num2str(Er,3)],'FontSize',FS,...
          'Units','Normalized','Backgroundcolor','w');
      end
    end
    
    % plot T-domain signals
    if isempty(analysis.hTdomain)
      hTfig = figure;
      set(hTfig,'Position',htPosition-[0 htPosition(4)*1.2,0 0]);
      hTax = gca;
    else
      hTax = analysis.hTdomain; cla(hTax)
      hTfig = gcf;
    end
    set(hTfig,'Units','Normalized');
    figure(hTfig), axes(hTax)
    if strcmpi(fitting.twoD,'y')
      hcon = surf(fitting.T2,fitting.T,S,'LineStyle','none');
      xlabel('T_b (s)','FontSize',FS)
      ylabel('T_a (s)','FontSize',FS)
      grid on, view([0 90]), shading interp
      set(gca,'FontSize',FS,'YScale','log','XScale','log');
      axis([min(fitting.T2) max(fitting.T2) min(fitting.T) max(fitting.T)]);
    else
      semilogx(fitting.T,S,'-','LineWidth',2)
      xlabel('Time Constant (s)','FontSize',FS)
      ylabel('Amplitude (a.u.)','FontSize',FS)
    end
    title('T-domain Spectrum','FontSize',FS+4,'FontWeight','Bold')
    
    set([htax,hTax],'FontSize',FS,'XGrid','on','YGrid','on','ZGrid','on');
  end

  function plotlabel
    % add labels to components on T-domain plot
    figure(hTfig), axes(hTax)
    
    if strcmpi(fitting.twoD,'y')
      for k = 1:Nc
        Ta_ms = 1e3*Tv(k,1); Tb_ms = 1e3*Tv(k,2);
        txt = sprintf('T_a = %s ms\nT_b = %s ms\nF = %s %%', ...
          num2str(Ta_ms,floor(log10(Ta_ms))+1), ...
          num2str(Tb_ms,floor(log10(Tb_ms))+1),num2str(100*Fv(k),3));
        xnorm = find(Tv(k,2)<=fitting.T2, 1)/fitting.numberT2-0.15;
        ynorm = find(Tv(k,1)<=fitting.T, 1)/fitting.numberT+0.15;
        
        htxt{k} = text(xnorm,ynorm,txt,...
          'VerticalAlignment','Top', ...
          'HorizontalAlignment','Left','Fontsize',FS,'Color','w',...
          'ButtonDownFcn',{@activateMove,hTfig},'UserData',false,...
          'Units',get(hTfig,'Units'),'FontWeight','Bold');
      end
      
    else
      
      % only label 1D components if only 1 spectrum is plotted
      if size(S,2) == 1
        v = axis;
        for k = 1:Nc
          txt = sprintf('T = %s s\nF = %s %%', ...
            num2str(Tv(k),2),num2str(100*Fv(k),3));
          xnorm = 1+log10(Tv(k))./diff(log10(v(1:2)));
          ynorm = (k-1)/Nc*0.7+0.2;
          
          htxt{k} = text(xnorm,ynorm,txt,...
            'VerticalAlignment','Top','HorizontalAlignment',...
            'Left','Backgroundcolor','w','Fontsize',FS,'Color','k',...
            'FontWeight','Bold','ButtonDownFcn',{@activateMove,hTfig}, ...
            'UserData',false,'Units',get(hTfig,'Units'));
        end
        
      end
      
    end
    
    % instruct that labels can be moved
    txt = 'click and drag to move component labels on spectrum plot';
    infodisplay(txt,analysis);
    grid on
    drawnow
    
    % peak label movement
    
    function activateMove(hobj, ~, fig)
      
      set(fig,'WindowButtonMotionFcn',{@moveTxtPosition, hobj});
      set(fig,'WindowButtonDownFcn',{@openTxtPosition, hobj});
      set(fig,'WindowButtonUpFcn',{@closeTxtPosition});
    end
    
    function openTxtPosition(fig, ~, hobj)
      
      set(fig,'WindowButtonMotionFcn',{@moveTxtPosition, hobj})
    end
    
    function moveTxtPosition(fig, ~, hobj)
      
      prevPos = get(fig,'UserData');
      if isempty(prevPos)
        set(fig,'UserData',get(fig,'CurrentPoint'));
        return;
      end
      
      curPos = get(fig,'CurrentPoint');
      posChange = curPos - prevPos;
      
      txtPos = get(hobj,'pos');
      set(hobj,'pos', [posChange + txtPos(1:2), 0]);
      
      set(fig,'UserData',curPos);
      
    end
    
    function closeTxtPosition(fig, ~, ~)
      
      set(fig, 'WindowButtonMotionFcn','');
      set(fig, 'WindowButtonDownFcn','');
      set(fig, 'WindowButtonUpFcn','');
      set(fig,'UserData',[]);
    end
    
  end

  function [Av,Fv,Tv,Nc] = peaks1D
    
    if strcmpi(analysis.extract,'auto')
      nx = cell(fitting.numbertrains,1);
      for kt = 1:fitting.numbertrains
        S1 = S(:,kt);
        % automatic 1D peak identification
        dy = diff(-S1);
        n = find(([dy' 0]<0) & ([0 dy']>=0));
        nx{kt} = n;
      end
      
      Nc = max(cellfun('length',nx))+1; % max number of components found
      
    else
      Nc = analysis.numberextract; % fixed Nc for manual extraction
      
      % instruct manual extraction
      txt = ['Identify ', int2str(Nc-1), ...
        ' component boundary(ies) from left to right with the mouse'];
      infodisplay(txt,analysis);
      
      % check to make sure selections are within axis
      v = axis;
      x = zeros(Nc,1); y = x; hline = zeros(Nc-1,1);n = hline;
      for kc = 1:Nc-1
        [x(kc),y(kc)] = ginput(1);
        if x(kc)<min(fitting.T) || x(kc)>max(fitting.T) || y(kc)<v(3) || ...
            y(kc)>v(4)
          txt = ['Error. Select boundaries within the domain and range ', ...
            'of the plotted spectrum'];
          infodisplay(txt,analysis)
          Tv=[];Av=[];Fv=[];
          return
        end
        hline(kc)=line([x(kc),x(kc)],[0,v(4)], ...
          'LineWidth',2,'LineStyle','--');
        drawnow
        n(kc) = find(fitting.T<x(kc), 1, 'last' );
      end
      pause(1); delete(hline)
      n = [1,n',fitting.numberT];
      
    end
    
    Tv = zeros(Nc,fitting.numbertrains); Av = Tv;
    
    for kt = 1:fitting.numbertrains
      if strcmpi(analysis.extract,'auto')
        n = [1,nx{kt},fitting.numberT];
      end
      for kc = 1:length(n)-1
        xp = fitting.T(n(kc):n(kc+1));
        yp = S(n(kc):n(kc+1),kt);
        % Find area, mean, and weight
        Av(kc,kt) = sum(yp);
        if Av(kc,kt)>0
          Tv(kc,kt) = exp(sum(yp.*log(xp))/sum(yp)); %geometric mean
        else
          Tv(kc,kt) = exp(sum(log(xp))/length(yp)); % avoid NaN where Av = 0
        end
      end
    end
    
    Fv = bsxfun(@rdivide,Av,sum(Av,1));
    
  end

  function [Av,Fv,Tv,Nc] = peaks2D
    
    if strcmpi(analysis.extract,'auto')
      % automatic 2D peak identification
      [r,c] = size(S);
      zz = find(S); xz = min(S(zz));
      %          zz = S; xz = min(S(zz));
      [Si,Nc] = bwlabel(im2bw(S./xz));
      Tv = zeros(Nc,2);
      Av = zeros(Nc,1);
      
      for kc = 1:Nc
        [rc,cc] = find(Si == kc);
        fc = diag(S(rc,cc));
        Tv(kc,1) = exp(sum(fc.*log(fitting.T(rc)))./sum(fc));
        Tv(kc,2) = exp(sum(fc.*log(fitting.T2(cc)))./sum(fc));
        Av(kc) = sum(fc);
      end
      Fv=Av/sum(Av);
      
    else
      Nc = analysis.numberextract;
      txt = ['identify ', int2str(Nc), ' components by drawing ', ...
        'regions. Left click to define a vertex, ', ...
        'right or double click to terminate a region' ];
      infodisplay(txt,analysis);
      
      % Enter range for each valley
      [nx,ny] = size(S);
      mask = zeros(nx,ny,Nc);
      v = axis;
      hold on
      for ii = 1:Nc
        if ii>1
          txt = 'Draw region around next component';
          if strcmpi(analysis.interactive,'y')
            infodisplay(txt,analysis);
          else
            display(txt)
          end
        end
        
        % Taking input from the mouse
        button = 1;
        ct = 1;
        indx = zeros(Nc,1); indy = zeros(Nc,1);
        xs = zeros(Nc,1); ys = zeros(Nc,1);
        while (button == 1)
          [xi,yi,button] =  ginput(1);
          if xi<v(1) || xi>v(2) || yi<v(3) || yi>v(4)
            Tv=[];Av=[];Fv=[];
            return
          end
          [~,indx(ct)] = min(abs(xi-fitting.T2));
          [~,indy(ct)] = min(abs(yi-fitting.T));
          xs(ct) = fitting.T2(indx(ct));
          ys(ct) = fitting.T(indy(ct));
          hline_dashed(ct) = plot(xs,ys,'m--','LineWidth',2.0); drawnow
          ct = ct + 1;
        end
        
        
        % Plot
        hline(ii) = plot(xs,ys,'m-','LineWidth',2.0);
        delete(hline_dashed)
        
        
        % Mask mask
        mask(:,:,ii) = poly2mask(indx,indy,nx,ny);
      end
      pause(1); delete(hline),
      mT2 = zeros(Nc,1);
      mT1 = zeros(Nc,1);
      Av = zeros(Nc,1);
      
      for kc = 1:Nc
        
        % Find x and y values for peak
        [rc,cc] = find(mask(:,:,kc));
        fc = diag(S(rc,cc));
        mT1(kc) = exp(sum(fc.*log(fitting.T(rc)))./sum(fc));
        mT2(kc) = exp(sum(fc.*log(fitting.T2(cc)))./sum(fc));
        Av(kc) = sum(fc);
        
      end
      
      Tv(:,1) = mT1;
      Tv(:,2) = mT2;
      
      % Weight of each peak
      Fv = Av./sum(Av);
      
    end
    
  end

end

function infodisplay(txt,analysis)

if strcmpi(analysis.interactive,'y')
  if exist('analysis.hpInfotxt','var')
    delete(analysis.hpInfotxt);
  end
  analysis.hpInfotxt = uicontrol('Style','text','Parent',analysis.hpInfo, ...
    'Units','Normalized','String',txt,'Position',[0.05 0.05,0.8,0.8],...
    'HorizontalAlignment','Left','Visible','on', 'ForegroundColor','b', ...
    'FontSize',14,'FontWeight','Bold','BackgroundColor',[0.9 0.85 0.65]);
else
  disp(txt)
end

end

function [data,fitting,analysis] = parseinput(data,fitting,analysis)

% Check and correct for invalid or incomplete inputs

if (~exist('data','var') || isempty(data))
  error('Input structure data is required');
end

dataflags = isfield(data,{'D','t','t2'});
if ~all(dataflags(1:2))
  error('Input structure ''data'' requires fields ''D'' and ''t''');
end

D = double(real(data.D));
t = double(real(data.t));
if dataflags(3)
  t2 = double(real(data.t2));
else
  t2 = t;
end

clear data

if (isempty(D)||isscalar(D)||isempty(t)||~isvector(t)||isempty(t2)|| ...
    ~isvector(t2))
  error(['Inputs data.D, data.t, and data.t2 (for 2D fitting)' ...
    ' must be real and non-empty. D must be a vector or a matrix.' ...
    ' t and t2 must be vectors.']);
end

if isrow(t), t = t'; end
if isrow(t2), t2 = t2'; end
if isrow(D), D = D'; end
if size(t,1) ~= size(D,1)
  error('data.D and data.t should have equal number of rows');
end

% fitting options

acceptablefittingfields = {'regtyp','regadj','numbergauss','widthgauss', ...
  'regweight','percentinc','rangeT','numberT','T','rangetheta',...
  'numbertheta','B1fit','nnlscode','twoD','rangeT2','numberT2','T2',...
  'theta_vector','theta_vector_fine','numberechoes',...
  'numberechoes2','H','theta','T0gauss','numberT0','numbertrains','fixedT1'};
defaultfittingvalue = {'mc','gcv',2,10,1e-1,1,[],100,[],[130 180],10, ...
  'n','nnlsmex','n',[],35,[],[],[],[],[],[],[],[],1,[],1};
defaultfitting =cell2struct(defaultfittingvalue,acceptablefittingfields,2);

if (exist('fitting','var')==0 || isempty(fitting))
  fitting.regtyp = defaultfitting.regtyp;
  fitting.regadj = defaultfitting.regadj;
  fitting.twoD = defaultfitting.twoD;
end

fittingfields = fieldnames(fitting);

for kfields = 1:length(fittingfields)
  if ~any(strcmp(fittingfields{kfields},acceptablefittingfields))
    disp(['Warning, fitting field ''', ...
      fittingfields{kfields}, ''' is not recognized'])
  end
end

fittingflags = isfield(fitting,acceptablefittingfields);

if fittingflags(14)
  twoD = lower(fitting.twoD);
  if ~strcmpi(twoD,'y')
    twoD = defaultfitting.twoD;
  end
else
  twoD = defaultfitting.twoD;
end

if fittingflags(1)
  regtyp = lower(fitting.regtyp);
else
  regtyp = defaultfitting.regtyp;
end

if ~any(strcmpi(regtyp,{'none','me','mc','mg','upen'}))
  regtyp = defaultfitting.regtyp;
end

if fittingflags(2)
  regadj = lower(fitting.regadj);
else
  if strcmpi(twoD,'y')
    regadj = 'manual';
  else
    regadj = defaultfitting.regadj;
  end
end

if ~any(strcmpi(regadj,{'manual','gcv','erinc'}))
  if strcmpi(twoD,'y')
    regadj = 'manual';
  else
    regadj = defaultfitting.regadj;
  end
end

if fittingflags(3)
  numbergauss = real(fitting.numbergauss);
  if ~isscalar(numbergauss)||~any(numbergauss == 1:10)
    numbergauss = defaultfitting.numbergauss;
  end
else
  numbergauss = defaultfitting.numbergauss;
end

if fittingflags(4)
  widthgauss = real(fitting.widthgauss);
  if ~isscalar(widthgauss)||(widthgauss < 0.4||widthgauss >40)
    widthgauss = defaultfitting.widthgauss;
  end
else
  widthgauss = defaultfitting.widthgauss;
end

if fittingflags(5)
  regweight = real(fitting.regweight);
  if ~isscalar(regweight)||(regweight < 0)
    regweight = defaultfitting.regweight;
  end
else
  regweight = defaultfitting.regweight;
end

if fittingflags(6)
  percentinc = real(fitting.percentinc);
  if (isempty(percentinc)||percentinc<=0)
    percentinc = defaultfitting.percentinc;
  end
else
  percentinc = defaultfitting.percentinc;
end

if fittingflags(7)
  rangeT = real(fitting.rangeT);
  if (length(rangeT)~=2||any(rangeT<=0))
    rangeT = [t(1) 2*t(end)];
  end
else
  rangeT = [t(1) 2*t(end)];
end

if fittingflags(8)
  numberT = round(real(fitting.numberT));
  if ~isscalar(numberT)||numberT<1|| ...
      (numberT>200&&strcmpi(analysis.interactive,'y'))
    if strcmpi(twoD,'y')
      numberT = defaultfitting.numberT2;
    else
      numberT = defaultfitting.numberT;
    end
  end
else
  if strcmpi(twoD,'y')
    numberT = defaultfitting.numberT2;
  else
    numberT = defaultfitting.numberT;
  end
end

if strcmpi(regtyp,'mg')
  regadj = 'none';
else
  if strcmpi(regadj,'none')
    regadj = 'manual';
  end
end

if fittingflags(9)
  T = real(fitting.T);
  if ~isvector(T) || strcmpi(regtyp,'mg')||strcmpi(analysis.interactive,'y')
    T = logspace(log10(rangeT(1)),log10(rangeT(2)),numberT);
  end
  rangeT = [T(1) T(end)];
  numberT = length(T);
else
  T = logspace(log10(rangeT(1)),log10(rangeT(2)),numberT);
end
if isrow(T), T = T'; end


if fittingflags(10)
  rangetheta = real(fitting.rangetheta);
  if length(rangetheta)~=2||any(rangetheta<90|rangetheta>180)
    rangetheta = defaultfitting.rangetheta;
  end
else
  rangetheta = defaultfitting.rangetheta;
end

if fittingflags(11)
  numbertheta = round(real(fitting.numbertheta));
  if ~isscalar(numbertheta)||numbertheta<1||numbertheta>100
    numbertheta = defaultfitting.numbertheta;
  end
else
  numbertheta = defaultfitting.numbertheta;
end

if fittingflags(12)
  B1fit = lower(fitting.B1fit);
  if ~strcmpi(B1fit,'y')
    B1fit = 'n'; rangetheta = [180 180]; numbertheta = 1;
  end
else
  B1fit = 'n'; rangetheta = [180 180]; numbertheta = 1;
end
theta = [];
theta_vector = linspace(rangetheta(1),rangetheta(2),numbertheta);
theta_vector_fine = linspace(rangetheta(1),rangetheta(2),100);

if fittingflags(27)
  fixedT1 = real(fitting.fixedT1);
  if ~isscalar(fixedT1)||(fixedT1 < 0)
    fixedT1 = defaultfitting.fixedT1;
  end
else
  fixedT1 = defaultfitting.fixedT1;
end


if fittingflags(13)
  nnlscode = fitting.nnlscode;
  if ~any(strcmpi(nnlscode,{'lsqnonneg','nnlsmex'}))
    nnlscode = 'nnlsmex';
  end
else
  nnlscode = 'nnlsmex';
end
if strcmpi(nnlscode,'nnlsmex') && ~(exist(['nnlsMEX.',mexext],'file')==3)
  nnlscode = 'lsqnonneg';
end

if fittingflags(15)
  rangeT2 = real(fitting.rangeT2);
  if (length(rangeT2)~=2||any(rangeT2<=0))
    rangeT2 = [0.75*t2(1) 2*t2(end)];
  end
else
  rangeT2 = [t2(1) 2*t2(end)];
end

if fittingflags(16)
  numberT2 = round(real(fitting.numberT2));
  if ~isscalar(numberT2)||numberT2<1|| ...
      (numberT2>200&&strcmpi(analysis.interactive,'y'))
    numberT2 = defaultfitting.numberT2;
  end
else
  numberT2 = defaultfitting.numberT2;
end

if fittingflags(17)
  T2 = real(fitting.T2);
  if ~isvector(T2)|| strcmpi(regtyp,'mg')||strcmpi(analysis.interactive,'y')
    T2= logspace(log10(rangeT2(1)),log10(rangeT2(2)),numberT2);
  end
  rangeT2 = [T2(1) T2(end)];
  numberT2 = length(T2);
else
  T2 = logspace(log10(rangeT2(1)),log10(rangeT2(2)),numberT2);
end
if isrow(T2), T2 = T2'; end

if ~strcmpi(twoD,'y')
  numberT2 = 0; T2 = []; rangeT2 = [];
end

if strcmpi(regtyp,'mg')
  
  if fittingflags(24)
    T0gauss = real(fitting.T0gauss);
    if size(T0gauss,1) ~= numbergauss
      T0gauss = [];
    end
  else
    T0gauss = [];
  end
  if isempty(T0gauss)
    T0gauss(:,1) = logspace(log10(rangeT(1)), ...
      log10(rangeT(2)),numbergauss+2)';
    
    if strcmpi(twoD,'y')
      T0gauss(:,2) = logspace(log10(rangeT2(1)), ...
        log10(rangeT2(2)),numbergauss+2)';
    end
    T0gauss = T0gauss(2:end-1,:);
  end
  
  if fittingflags(25)
    numberT0 = round(real(fitting.numberT0));
  else
    numberT0 = defaultfitting.numberT0;
  end
  if isnumeric(numberT0) && isscalar(numberT0) && ~isempty(numberT0)
    if numberT0<1
      disp('fitting.numberT0 set 1')
      numberT0 = 1;
    end
    if numberT0>100
      disp(['warning: large number of varied initial guesses for T0 ', ...
        'may result in long run time']);
    end
  else
    numberT0 = defaultfitting.numberT0;
  end
  
  T0gauss = bsxfun(@plus,bsxfun(@times,randn([size(T0gauss),numberT0]),...
    T0gauss*0.15),T0gauss);
else
  T0gauss = [];
  numberT0 = [];
end

% analysis options

if (exist('analysis','var')==0 || isempty(fitting))
  analysis.graph = 'n';
  analysis.range = [min(rangeT) max(rangeT)];
  analysis.tolextract = 2;
end

analysisflags = isfield(analysis,{'graph','rangeA','tolextract', ...
  'interactive','rangeA2','htdomain','hTdomain','extract',...
  'numberextract','hpInfo','hpInfotxt'});

if analysisflags(1)
  graph = lower(analysis.graph);
  if ~strcmpi(graph,'y')
    graph = 'n';
  end
else
  graph = 'n';
end


if analysisflags(2)
  rangeA = real(analysis.rangeA);
  if length(rangeA)~=2||min(rangeA)<min(rangeT)||max(rangeA)>max(rangeT)
    rangeA = [min(rangeT) max(rangeT)];
  end
else
  rangeA = [min(rangeT) max(rangeT)];
end
rangeA = rangeA+[-eps,eps]; % avoid prob with round-off errors in Tv calc

if analysisflags(5)
  rangeA2 = real(analysis.rangeA2);
  if length(rangeA2)~=2||min(rangeA2)<min(rangeT2)||max(rangeA2)>max(rangeT2)
    rangeA2 = [min(rangeT2) max(rangeT2)];
  end
else
  rangeA2 = [min(rangeT2) max(rangeT2)];
end
if ~isempty(rangeA2)
  rangeA2 = rangeA2+[-eps,eps]; % avoid prob with round-off errors in Tv calc
end

if analysisflags(3)
  tolextract = real(analysis.tolextract);
  if (isempty(tolextract)||tolextract<0)
    tolextract = 2;
  end
else
  tolextract = 2;
end

if analysisflags(4)
  interactive = analysis.interactive;
  if ~strcmpi(interactive,'y')
    interactive = 'n';
  else
    graph = 'y';
  end
else
  interactive = 'n';
end

if analysisflags(6)
  htdomain = analysis.htdomain;
else
  htdomain = [];
end

if analysisflags(7)
  hTdomain = analysis.hTdomain;
else
  hTdomain = [];
end

if analysisflags(8)
  extract = analysis.extract;
  if ~strcmpi(extract,'manual')
    extract = 'auto';
  end
else
  extract = 'auto';
end
if strcmpi(regtyp,'mg')&&strcmpi(extract,'manual')
  txt = ['Manual peak extraction not permitted when fitting to '....
    'multiple gaussian components' ];
  display(txt);
  extract = 'auto';
end
if strcmpi(extract,'manual')
  graph = 'y';
end

if analysisflags(9)
  numberextract = real(analysis.numberextract);
  if ~isscalar(numberextract)||~any(numberextract == 1:10)
    numberextract = 2;
  end
else
  numberextract = 2;
end

if analysisflags(10)
  hpInfo = analysis.hpInfo;
else
  hpInfo = [];
end

if analysisflags(11)
  hpInfotxt = analysis.hpInfotxt;
else
  hpInfotxt = [];
end

if strcmpi(twoD,'y')
  [numberechoes,numberechoes2] = size(D);
  numbertrains = 1; % 2D fitting only for 1 data set at a time, for now
else
  [numberechoes,numbertrains] = size(D);
  numberechoes2 = 1;
  if strcmpi(interactive,'y') && numbertrains > 1
    txt = ['WARNING. Interactive analysis proceeding ',...
      'for the first signal only'];
    disp(txt),beep,pause(1)
    numbertrains = 1;
    D = D(:,1);
  end
end

% send objects back in structures

fitting = struct('regtyp',regtyp,'regadj',regadj,'numbergauss',numbergauss, ...
  'widthgauss',widthgauss,'regweight',regweight, 'twoD',twoD, ...
  'percentinc',percentinc,'rangeT',rangeT,'numberT',numberT,'T',T,...
  'rangetheta',rangetheta,'numbertheta', numbertheta,'B1fit',B1fit,...
  'theta_vector',theta_vector,'theta_vector_fine',theta_vector_fine,...
  'nnlscode',nnlscode,'rangeT2',rangeT2,'numberT2',numberT2,'T2',T2,...
  'theta',theta,'T0gauss',T0gauss,'numberT0',numberT0,...
  'numberechoes',numberechoes,'numberechoes2',numberechoes2,...
  'numbertrains',numbertrains,'fixedT1',fixedT1);

data = struct('D',D,'t',t,'t2',t2);

analysis = struct('graph',graph,'rangeA',rangeA,'tolextract',tolextract, ...
  'interactive',interactive,'rangeA2',rangeA2,'htdomain',htdomain, ...
  'hTdomain',hTdomain,'extract',extract,'numberextract',numberextract, ...
  'hpInfo',hpInfo,'hpInfotxt',hpInfotxt);

end

function x = nnlsLH(A,b,nnlscode)
%NNLS    Non-negative least-sqaures
%   X = NNLS(A,B) solves the least-squares problem
%   A*X = B subject to X >= 0
%
%   INPUT:
%       A = m x n matrix
%       b = m x 1 array
%
%   OUTPUT:
%       x = n x 1 array
%
%   The original version of this code was developed by
%   Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
%   1973 JUN 15, and published in the book
%   "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
%   Revised FEB 1995 to accompany reprinting of the book by SIAM.
%
%   The c-code distributed by the ASYNCHRONOUS PARALLEL PATTERN SEARCH (APPS)
%   was modified so as to create a MEX-file.
%

switch nnlscode
  
  case 'nnlsmex'
    % Get array dims
    [m,n] = size(A);
    
    % Check array dims
    if ((size(b,1) ~= m) || (size(b,2) ~= 1))
      error('Invalid variable dimension')
    end
    if ((size(A,1) ~= m) || (size(A,2) ~= n))
      error('Invalid variable dimension')
    end
    
    % Cast MATLAB variables for MEX funtion
    A = full(A);
    m = int32(m);
    n = int32(n);
    b = full(b);
    rnorm = 0;
    W = zeros(n,1);
    ZZ = zeros(m,1);
    IDX = int32(zeros(m*1e2,1));
    mode = int8(0);
    
    % Call MEX function
    At = A + 0;
    bt = b + 0;
    x = nnlsMEX(At,m,m,n,bt,rnorm,W,ZZ,IDX,mode);
    
  case 'lsqnonneg'
    x = lsqnonneg(A,b);
end
end

function [output,fitting,analysis]= MERA_GUI(data,fitting,analysis)
%% Create GUI
% Define standard dimensions
ScreenSize = get(0,'ScreenSize');
FigWidth = round(0.75*ScreenSize(3));
FigHeight = round(FigWidth/(16/9));
FontSize1 = floor(FigWidth/100);
FontSize2 = floor(0.9*FontSize1);
bkcolor = [0.9 0.85 0.65];
f = figure('Visible','off', ...
  'Position',[ScreenSize(3)-FigWidth,ScreenSize(4)-FigHeight,...
  FigWidth,FigHeight], ...
  'MenuBar','none','Name','MERA INTERACTIVE','NumberTitle','off', ...
  'PaperOrientation','landscape','Toolbar','none', ...
  'PaperType','uslegal','Color',bkcolor,'Resize','off');


%% Construct the components.

% Plotting Axes
analysis.htdomain = axes('Position',[.27 0.55 0.49 0.35],'Color',[1 1 1 ]);
httitle = uicontrol('Style','text','Units','pixels',...
  'String','Time Domain','HorizontalAlignment','Center',...
  'Position',[0.45 0.905 0.2 0.02],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);

analysis.hTdomain = axes('Position',[.27 0.1 0.49 0.35],'Color',[1 1 1]);
hTtitle = uicontrol('Style','text','Units','pixels',...
  'String','Spectral Domain','HorizontalAlignment','center',...
  'Position',[0.45 0.455 0.2 0.02],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);

set([httitle, hTtitle], ...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);


%% Control Panels
% B1 Fit Toggle
VscaleB1 = 0.125;
NelemB1 = 3+2;
EscaleB1 = 1./NelemB1;
hpB1 = uipanel('BackgroundColor',bkcolor, ...
  'Title','B1 Fitting', 'Position',[0.01 0.75 0.2 VscaleB1]);

% Regularization Panel
VscaleReg = 0.2;
NelemReg = 3+2;
EscaleReg = 1./NelemReg;
hpRegularize = uipanel('BackgroundColor',bkcolor,'Title','Regularization',...
  'Position',[0.01 0.5 0.2 VscaleReg]);

% T-domain Panel
VscaleTD = 0.3;
NelemTD = 8+3;
EscaleTD = 1./NelemTD;
hpTdomain = uipanel('BackgroundColor',bkcolor, ...
  'Title','Number/Range of Time Constants','Position',[0.01 0.15 0.2 VscaleTD]);

% Exit Pushbutton
hexit = uicontrol('Style','pushbutton','Units','Normalized',...
  'String',{'Exit'},...
  'Position',[0.84 0.82 0.1 0.05],...
  'TooltipString','Push to Exit GUI and return Output', ...
  'Callback',{@exit_callback});

% Instruction Panel
VscaleReg = 0.2;
NelemReg = 3;
EscaleInfo = 1./NelemReg;
analysis.hpInfo = uipanel('BackgroundColor',bkcolor,'Title','Information',...
  'Position',[0.79 0.58 0.20 VscaleReg]);
txt = ' ';
analysis.hpInfotxt = uicontrol('Style','text','Parent',analysis.hpInfo, ...
  'Units','Normalized','String',txt,'Position',[0.05 0.05,0.8,0.8],...
  'HorizontalAlignment','Left','Visible','on', 'ForegroundColor','b', ...
  'FontSize',14,'FontWeight','Bold','BackgroundColor',[0.9 0.85 0.65]);

% Analysis Panel
VscaleAn = 0.2;
NelemAn = 5;
EscaleAn = 1./NelemAn;
hpAnalysis = uipanel('BackgroundColor',bkcolor,'Title','Analysis',...
  'Position',[0.79 0.28 0.20 VscaleAn]);

set([hpB1,hpRegularize,hpTdomain,analysis.hpInfo,hexit,hpAnalysis],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);


%% B1 Panel Items
hB1fit = uicontrol('Style','checkbox','Parent',hpB1,'Units','Normalized',...
  'String','Y/N','Position',[0.025 3.5*EscaleB1,0.25,EscaleB1],...
  'TooltipString','Fit Imperfect Refocussing using EPG', ...
  'Callback',{@B1fit_callback});
if strcmpi(fitting.B1fit,'y')
  set(hB1fit,'Value',true)
else
  set(hB1fit,'Value',false)
end

hB1fittxt = uicontrol('Style','text','Parent',hpB1,...
  'String',['B1 fit = ',num2str(fitting.theta,3),'°'],'Units','Normalized',...
  'Position',[0.35 3.5*EscaleB1 0.6 EscaleB1],'HorizontalAlignment','Left');

hnumberthetatxt = uicontrol('Style','text','Parent',hpB1,...
  'String','Number of Flip Angles:','Units','Normalized',...
  'Position',[0.025 2*EscaleB1 0.6 EscaleB1],'HorizontalAlignment','Left');

hnumbertheta = uicontrol('Style','edit','Parent',hpB1,'Units','Normalized',...
  'String',int2str(fitting.numbertheta),'Min',1,'Max',1, ...
  'SelectionHighlight','off','Position', ...
  [0.6 1.75*EscaleB1 0.12 1.5*EscaleB1],'TooltipString', ...
  'Number of flip angles to test to fit B1','Callback',{@B1fit_callback});

hrangethetatxt = uicontrol('Style','text','Parent',hpB1,...
  'String','Range of Flip Angles:','Units','Normalized',...
  'Position',[0.025 0.5*EscaleB1 0.6 EscaleB1],'HorizontalAlignment','Left');

hrangethetaL = uicontrol('Style','edit','Parent',hpB1,'Units','Normalized',...
  'String',int2str(fitting.rangetheta(1)),'Min',1,'Max',1,...
  'SelectionHighlight','off','Position', ...
  [0.6 0.25*EscaleB1 0.15 1.5*EscaleB1],'TooltipString', ...
  'lowest flip angle to test','Callback',{@B1fit_callback});

hrangethetaU = uicontrol('Style','edit','Parent',hpB1,'Units','Normalized',...
  'String',int2str(fitting.rangetheta(2)),'Min',1,'Max',1, ...
  'SelectionHighlight','off','Position', ...
  [0.8 0.25*EscaleB1 0.15 1.5*EscaleB1],'TooltipString',...
  'highest flip angle to test','Callback',{@B1fit_callback});

set([hpB1,hB1fit,hB1fittxt,hnumberthetatxt,hnumbertheta,hrangethetatxt,...
  hrangethetaL,hrangethetaU],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);
set([hnumbertheta,hrangethetaL,hrangethetaU],'FontSize',FontSize2);


%% Regularization Panel Items

hregtyp = uicontrol('Style','popupmenu','Parent',hpRegularize,...
  'Units','Normalized','String',{'None','Min Energy','Min Curvature',...
  'Uniform Penalty','Multi Gauss'},'Value',...
  find(strcmpi(fitting.regtyp,{'none','me','mc','upen','mg'})), ...
  'Position',[0.025 3.5*EscaleReg,0.95,EscaleReg],...
  'TooltipString','Select Regularization Type','Callback',{@regtyp_callback});

hregadj = uicontrol('Style','popupmenu','Parent',hpRegularize,...
  'Units','Normalized','String', ...
  {'Manual','Cross Validation','% MSE Increase','None'}, ...
  'Value',find(strcmpi(fitting.regadj,{'manual','gcv','erinc','none'})), ...
  'Position',[0.025 2*EscaleReg,0.95,EscaleReg],...
  'TooltipString','Select Adjustment of Regularizer Weight',...
  'Callback',{@regadj_callback});

hregweight = uicontrol('Style','slider','Parent',hpRegularize, ...
  'Units','Normalized','Min',0.0001,'Max',1,'Value',fitting.regweight,...
  'SliderStep', [0.001,.1],'Position',[0.05 0.5*EscaleReg,0.4,EscaleReg],...
  'TooltipString','Adjust Regularizer Weight','Callback',{@slider_callback});

hregweighttxt = uicontrol('Style','text','Parent',hpRegularize,...
  'Units','Normalized','String', ...
  ['Reg Weight = ',num2str(fitting.regweight,2)],...
  'Position',[0.525 0.5*EscaleReg,0.45,EscaleReg],'HorizontalAlignment','Left');

hnumbergauss = uicontrol('Style','edit','Parent',hpRegularize,...
  'Units','Normalized','String',int2str(fitting.numbergauss),'Min',1,'Max',1,...
  'Position',[0.3 2*EscaleReg,0.15,EscaleReg],...
  'Callback',{@slider_callback});

hnumbergausstxt = uicontrol('Style','text','Parent',hpRegularize,...
  'Units','Normalized','String','# Comp''s',...
  'Position',[0.05 1.8*EscaleReg,0.25,EscaleReg],...
  'HorizontalAlignment','Left',...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);

hnumberT0 = uicontrol('Style','edit','Parent',hpRegularize,...
  'Units','Normalized','String',int2str(fitting.numberT0),'Min',1,'Max',1,...
  'Position',[0.3 0.7*EscaleReg,0.15,EscaleReg],...
  'Callback',{@slider_callback});

hnumberT0txt = uicontrol('Style','text','Parent',hpRegularize,...
  'Units','Normalized','String','# varied initial T0 conditions',...
  'Position',[0.05 0.5*EscaleReg,0.25,EscaleReg],...
  'HorizontalAlignment','Left',...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);

hwidthgauss = uicontrol('Style','slider','Parent',hpRegularize,...
  'Units','Normalized','Min',0.4,'Max',20,'Value',...
  floor(fitting.widthgauss), 'SliderStep', [0.01,0.1],...
  'Position',[0.55 1.5*EscaleReg,0.40,EscaleReg],...
  'TooltipString','Adjust Width of Gaussian Components',...
  'Callback',{@slider_callback});

hwidthgausstxt = uicontrol('Style','text','Parent',hpRegularize,...
  'Units','Normalized','String',['Width = ',int2str(fitting.widthgauss),'%'],...
  'Position',[0.55 2.5*EscaleReg,0.4,EscaleReg]);

herinc = uicontrol('Style','slider','Parent',hpRegularize,...
  'Units','Normalized','Min',0.1,'Max',5,'Value',...
  fitting.percentinc, 'SliderStep', [0.01,0.1],...
  'Position',[0.05 0.25*EscaleReg,0.4,0.5*EscaleReg],...
  'TooltipString','Adjust Regularization Weight by Percent Error Increase', ...
  'Callback',{@slider_callback});

herinctxt = uicontrol('Style','text','Parent',hpRegularize,...
  'Units','Normalized','String',...
  ['% MSE Incr = ',num2str(fitting.percentinc,2)],...
  'Position',[0.025 0.75*EscaleReg,0.4,EscaleReg]);

% set common properties
set([hpRegularize,hregtyp,hregadj,hregweighttxt,hnumbergauss,...
  hnumbergausstxt,hwidthgauss,hwidthgausstxt,herinc,herinctxt,...
  hnumberT0,hnumberT0txt],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);
set([hregtyp,hregadj,hnumbergauss],'FontSize',FontSize2);

%% T-Domain Panel Items

hnumberT = uicontrol('Style','slider','Parent',hpTdomain,...
  'Units','Normalized','Min',10,'Max',250,'Value',...
  round(fitting.numberT), 'SliderStep', [0.01,0.1],...
  'Position',[0.05 8.5*EscaleTD,0.45,EscaleTD],...
  'TooltipString','Adjust Number of Time Constants to Fit', ...
  'Callback',{@slider_callback});

hnumberTtxt = uicontrol('Style','text','Parent',hpTdomain,...
  'Units','Normalized','String',['# TConst = ',int2str(fitting.numberT)],...
  'Position',[0.55 8.5*EscaleTD,0.4,EscaleTD],'HorizontalAlignment','Left');

hrangeTL = uicontrol('Style','slider','Parent',hpTdomain,...
  'Units','Normalized','Min',data.t(1)/10,'Max',data.t(end)/10,'Value',...
  fitting.rangeT(1), 'SliderStep', [0.01,0.1],'Position',...
  [0.05 6*EscaleTD,0.40,EscaleTD],'TooltipString',...
  'Minimum Time Constant to Fit','Callback',{@slider_callback});

hrangeTLtxt = uicontrol('Style','text','Parent',hpTdomain,...
  'Units','Normalized','String',...
  ['Min T = ',num2str(fitting.rangeT(1)*1e3,3),' ms'],'Position', ...
  [0.05 7*EscaleTD,0.4,EscaleTD]);

hrangeTU = uicontrol('Style','slider','Parent',hpTdomain,...
  'Units','Normalized','Min',data.t(end)/9,'Max',data.t(end)*10,...
  'Value',fitting.rangeT(2), 'SliderStep', [0.01,0.1],'Position',...
  [0.5 6*EscaleTD,0.40,EscaleTD],'TooltipString', ...
  'Minimum Time Constant to Fit','Callback',{@slider_callback});

hrangeTUtxt = uicontrol('Style','text','Parent',hpTdomain,...
  'Units','Normalized','String',...
  ['Max T = ',num2str(fitting.rangeT(2),3),' s'],'Position',...
  [0.5 7*EscaleTD,0.4,EscaleTD]);

% set common properties
set([hpTdomain,hnumberT,hnumberTtxt,hrangeTL,hrangeTLtxt,...
  hrangeTU,hrangeTUtxt],...
  'FontSize',FontSize2,'FontWeight','Bold','BackgroundColor',bkcolor);
set(hnumberTtxt,'FontSize',FontSize2)

if strcmpi(fitting.twoD, 'y')
  
  hrangeT1label = uicontrol('Style','text','Parent',hpTdomain,...
    'Units','Normalized','String','1st Dimension',...
    'Position',[0.1 9.5*EscaleTD,0.8,EscaleTD]);
  
  hrangeT2label = uicontrol('Style','text','Parent',hpTdomain, ...
    'Units','Normalized','String',['2nd Dimension'],...
    'Position',[0.1 4*EscaleTD,0.8,EscaleTD]);
  
  hnumberT2 = uicontrol('Style','slider','Parent',hpTdomain,...
    'Units','Normalized','Min',10,'Max',250,'Value',round(fitting.numberT2),...
    'SliderStep', [0.01,0.1],'Position',[0.05 3*EscaleTD,0.45,EscaleTD],...
    'TooltipString','Adjust Number of Time Constants to Fit', ...
    'Callback',{@slider_callback});
  
  hnumberT2txt = uicontrol('Style','text','Parent',hpTdomain,...
    'Units','Normalized','String', ...
    ['#Time Const = ',int2str(fitting.numberT2)],...
    'Position',[0.55 3*EscaleTD,0.4,EscaleTD],...
    'HorizontalAlignment','Left');
  
  hrangeT2L = uicontrol('Style','slider','Parent',hpTdomain,...
    'Units','Normalized','Min',data.t2(1)/10,'Max',data.t2(1)*10,...
    'Value',fitting.rangeT2(1), 'SliderStep', [0.01,0.1],...
    'Position',[0.05 0.5*EscaleTD,0.40,EscaleTD],...
    'TooltipString','Minimum Time Constant to Fit',...
    'Callback',{@slider_callback});
  
  hrangeT2Ltxt = uicontrol('Style','text','Parent',hpTdomain,...
    'Units','Normalized','String',...
    ['Min T = ',num2str(fitting.rangeT2(1)*1e3,3),' ms'],...
    'Position',[0.05 1.5*EscaleTD,0.4,EscaleTD]);
  
  hrangeT2U = uicontrol('Style','slider','Parent',hpTdomain,...
    'Units','Normalized','Min',data.t2(end)/10,'Max',data.t2(end)*10,...
    'Value',fitting.rangeT2(2), 'SliderStep', [0.01,0.1],...
    'Position',[0.5 0.5*EscaleTD,0.40,EscaleTD],...
    'TooltipString','Minimum Time Constant to Fit',...
    'Callback',{@slider_callback});
  
  hrangeT2Utxt = uicontrol('Style','text','Parent',hpTdomain,...
    'Units','Normalized','String',...
    ['Max T = ',num2str(fitting.rangeT2(2),3),' s'],...
    'Position',[0.5 1.5*EscaleTD,0.4,EscaleTD]);
  
  % set common properties
  set([hnumberT2,hnumberT2txt,hrangeT2L,hrangeT2Ltxt,hrangeT2U,hrangeT2Utxt, ...
    hrangeT1label,hrangeT2label],...
    'FontSize',FontSize2,'FontWeight','Bold','BackgroundColor',bkcolor);
  set([hnumberT2txt,hrangeT2Ltxt],'FontSize',FontSize2)
  set(hpB1,'Visible','off')
else
  hnumberT2 = []; hnumberT2txt = [];
end

%% Analysis Panel Items
hgrid = uicontrol('Style','checkbox','Parent',hpAnalysis,...
  'Units','Normalized','String','grid','Position',...
  [0.025 3.5*EscaleAn,0.2,EscaleAn],'TooltipString','Toggle Grid Lines',...
  'Value',true,'Callback',{@grid_callback});
hlog = uicontrol('Style','checkbox','Parent',hpAnalysis,'Units','Normalized',...
  'String','log scale','Position',[0.275 3.5*EscaleAn,0.45,EscaleAn],...
  'TooltipString','Toggle Log Scale', 'Value',false, ...
  'Callback',{@log_callback});
hextract = uicontrol('Style','popupmenu','Parent',hpAnalysis,...
  'Units','Normalized','String',{'Auto','Manual'},...
  'Value',find(strcmpi(analysis.extract,{'auto','manual'})), ...
  'Position',[0.025 2*EscaleAn,0.35,EscaleAn],...
  'TooltipString','Select Regularization Type','Callback',{@extract_callback});
hnumberextract = uicontrol('Style','edit','Parent',hpAnalysis,...
  'Units','Normalized','String',...
  int2str(analysis.numberextract),'Min',1,'Max',1,'Position',...
  [0.25 .7*EscaleAn,0.12,EscaleAn],'Visible','off','TooltipString',...
  'Number of Peaks to Manually Extract','Callback',{@extract_callback});
hnumberextracttxt = uicontrol('Style','text','Parent',hpAnalysis,...
  'Units','Normalized','String','#Peaks','Position',...
  [0.05 0.5*EscaleAn,0.2,EscaleAn],'HorizontalAlignment',...
  'Left','Visible','off','FontSize',FontSize1,'FontWeight','Bold',...
  'BackgroundColor',bkcolor);
htolextract = uicontrol('Style','edit','Parent',hpAnalysis,...
  'Units','Normalized','String',num2str(analysis.tolextract),...
  'Min',1,'Max',1,'Position',[0.25 .7*EscaleAn,0.12,EscaleAn],...
  'TooltipString','Ignore Peaks below this Percent Amplitude', ...
  'Callback',{@extract_callback});
htolextracttxt = uicontrol('Style','text','Parent',hpAnalysis,...
  'Units','Normalized','String','Tol (%)','Position', ...
  [0.05 0.5*EscaleAn,0.2,EscaleAn],'HorizontalAlignment','Left',...
  'Visible','on','FontSize',FontSize1,'FontWeight','Bold',...
  'BackgroundColor',bkcolor);
hrangeAL = uicontrol('Style','edit','Parent',hpAnalysis,'Units','Normalized',...
  'String',num2str(analysis.rangeA(1)),'Min',1,'Max',1, ...
  'Position',[0.5 1*EscaleAn,0.2,EscaleAn], ...
  'TooltipString','Range of Time-Constants for Peak Extraction', ...
  'Callback',{@extract_callback});
hrangeAU = uicontrol('Style','edit','Parent',hpAnalysis,'Units','Normalized',...
  'String',num2str(analysis.rangeA(2)),'Min',1,'Max',1, ...
  'Position',[0.7 1*EscaleAn,0.2,EscaleAn],...
  'TooltipString','Range of Time-Constants for Peak Extraction', ...
  'Callback',{@extract_callback});
hrangeAtxt = uicontrol('Style','text','Parent',hpAnalysis,...
  'Units','Normalized','String','Range (s)',...
  'Position',[0.6 2.25*EscaleAn,0.25,EscaleAn],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);
if strcmpi(analysis.extract,'auto');
  set([htolextract,htolextracttxt,hrangeAL,hrangeAU,hrangeAtxt],'Visible','on');
  set([hnumberextract,hnumberextracttxt],'Visible','off');
else
  set([hnumberextract,hnumberextracttxt],'Visible','on');
  set([htolextract,htolextracttxt,hrangeAL,hrangeAU,hrangeAtxt],...
    'Visible','off');
end


if strcmpi(fitting.twoD, 'y')
  hrangeAtxt2 = uicontrol('Style','text','Parent',hpAnalysis,...
    'Units','Normalized','String','lower   upper',...
    'Position',[0.5 1.95*EscaleAn,0.4,0.6*EscaleAn],...
    'FontSize',FontSize2,'FontWeight','Bold','BackgroundColor',bkcolor);
  hrangeAtxt3 = uicontrol('Style','text','Parent',hpAnalysis,...
    'Units','Normalized','String','1D',...
    'Position',[0.9 1.2*EscaleAn,0.07,0.6*EscaleAn],...
    'FontSize',FontSize2,'FontWeight','Bold','BackgroundColor',bkcolor);
  hrangeAtxt4 = uicontrol('Style','text','Parent',hpAnalysis,...
    'Units','Normalized','String','2D',...
    'Position',[0.9 0.2*EscaleAn,0.07,0.6*EscaleAn],...
    'FontSize',FontSize2,'FontWeight','Bold','BackgroundColor',bkcolor);
  hrangeA2L = uicontrol('Style','edit','Parent',hpAnalysis,...
    'Units','Normalized','String',num2str(analysis.rangeA2(1)),...
    'Min',1,'Max',1,'Position',[0.5 0.1*EscaleAn,0.2,EscaleAn],...
    'TooltipString','Range of Time-Constants for Peak Extraction', ...
    'Callback',{@extract_callback});
  hrangeA2U = uicontrol('Style','edit','Parent',hpAnalysis,...
    'Units','Normalized','String',num2str(analysis.rangeA2(2)),...
    'Min',1,'Max',1,'Position',[0.7 0.1*EscaleAn,0.2,EscaleAn],...
    'TooltipString','Range of Time-Constants for Peak Extraction', ...
    'Callback',{@extract_callback});
  set([hrangeA2L,hrangeA2U,hrangeAtxt2,hrangeAtxt3,hrangeAtxt4],...
    'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);
  
  if strcmpi(analysis.extract,'auto');
    set([hrangeA2L,hrangeA2U,hrangeAtxt2,hrangeAtxt3,hrangeAtxt4],...
      'Visible','on');
  else
    set([hrangeA2L,hrangeA2U,hrangeAtxt2,hrangeAtxt3,hrangeAtxt4],...
      'Visible','off');
  end
else
  hrangeA2U = [];hrangeA2L = []; hrangeAtxt2 = []; hrangeAtxt3 =[];
  hrangeAtxt4 = [];
  
end

set([hgrid,hlog,hextract,hnumberextracttxt,hnumberextract,htolextract,...
  htolextracttxt,hrangeAtxt,hrangeAL,hrangeAU],...
  'FontSize',FontSize1,'FontWeight','Bold','BackgroundColor',bkcolor);

%%
updategraphics();
set(f,'Visible','on');
[data,fitting,analysis] = parseinput(data,fitting,analysis);
[output,fitting,analysis] = fittingloop(data,fitting,analysis);
updategraphics();

uiwait(f);
close(f);
return
%% Sub Functions
  function updategraphics
    
    set(hregtyp,'Value',find(strcmpi(fitting.regtyp, ...
      {'none','me','mc','upen','mg'})));
    set(hregadj,'Value',find(strcmpi(fitting.regadj, ...
      {'manual','gcv','erinc','none'})));
    
    switch(fitting.regtyp)
      case 'none'
        set([hregadj,hregweight,hregweighttxt, ...
          hwidthgauss,hwidthgausstxt, hnumberT0, hnumberT0txt, ...
          hnumbergauss, hnumbergausstxt],'Visible','off')
        set([hnumberT,hnumberT2,hnumberTtxt,hnumberT2txt],'Visible','on')
      case 'mg'
        set([hregadj,hregweight,hregweighttxt,hnumberTtxt,hnumberT2txt ...
          herinctxt,herinc,hnumberT,hnumberT2],'Visible','off')
        set([hwidthgauss,hwidthgausstxt,hnumbergauss, ...
          hnumberT0, hnumberT0txt, hnumbergausstxt],'Visible','on')
        set(hnumberT0,'String',int2str(fitting.numberT0));
      case {'me','mc','upen'}
        set([hregadj,hregweight,hregweighttxt,hnumberT,hnumberT2,...
          hnumberTtxt,hnumberT2txt],'Visible','on')
        set([hwidthgauss,hwidthgausstxt,hnumberT0, hnumberT0txt, ...
          hnumbergauss, hnumbergausstxt],'Visible','off')
        
        switch(fitting.regadj)
          case 'manual'
            set([hregweight,hregweighttxt],'Visible','on');
            set([herinc,herinctxt],'Visible','off');
          case 'gcv'
            set(hregweighttxt,'Visible','on');
            set([herinc,herinctxt,hregweight],'Visible','off');
          case 'erinc'
            set([herinc,herinctxt,hregweighttxt],'Visible','on');
            set(hregweight,'Visible','off');
          case 'none'
            set([herinc,herinctxt,hregweight,hregweighttxt],'Visible','off');
            
        end
    end
    
    if strcmpi(fitting.B1fit,'y')
      set(hB1fittxt,'String',['B1 fit = ', num2str(fitting.theta,4),'°']);
      set(hnumbertheta,'String',int2str(fitting.numbertheta));
      set(hrangethetaU,'String',int2str(fitting.rangetheta(2)));
      set(hrangethetaL,'String',int2str(fitting.rangetheta(1)));
    else
      set(hB1fittxt,'String','B1 fit','ForegroundColor',[0.5 0.5 0.5])
      set([hnumberthetatxt,hnumbertheta,hrangethetatxt,hrangethetaL,...
        hrangethetaU],'ForegroundColor',[0.5 0.5 0.5])
    end
    
    set(hregweighttxt,'String', ...
      ['Reg Weight = ',num2str(fitting.regweight,2)]);
    set(hnumberTtxt,'String', ...
      ['#Time Const = ',int2str(fitting.numberT)]);
    
  end

  function log_callback(~,~)
    if get(hlog,'Value')
      if strcmpi(fitting.twoD,'y')
        set(analysis.htdomain,'ZScale','log');
        zlim(analysis.htdomain,[output.Er(1)/2,max(output.P(1))*1.2]);
      else
        set(analysis.htdomain,'YScale','log');
        ylim(analysis.htdomain,[output.Er(1)/2,max(output.P(1))*1.2])
      end
    else
      set(analysis.htdomain,'YScale','linear','ZScale','linear');
      ylim(analysis.htdomain,'auto'), zlim(analysis.htdomain,'auto');
    end
  end

  function grid_callback(~,~)
    if get(hgrid,'Value')
      set([analysis.htdomain,analysis.hTdomain], ...
        'XGrid','on','YGrid','on','Zgrid','off');
    else
      set([analysis.htdomain,analysis.hTdomain], ...
        'XGrid','off','YGrid','off','ZGrid','off');
    end
  end

  function extract_callback(~,~)
    % Determine the selected data set.
    str = get(hextract,'String');
    val = get(hextract,'Value');
    % Set current data to the selected data set.
    
    % get new values for Peaks, tol, rangeA
    analysis.numberextract = str2double(get(hnumberextract,'String'));
    analysis.tolextract = str2double(get(htolextract,'String'));
    analysis.rangeA(1) = str2double(get(hrangeAL,'String'));
    analysis.rangeA(2) = str2double(get(hrangeAU,'String'));
    
    if strcmpi(fitting.twoD,'y')
      analysis.rangeA2(1) = str2double(get(hrangeA2L,'String'));
      analysis.rangeA2(2) = str2double(get(hrangeA2U,'String'));
    end
    
    if strcmpi(fitting.regtyp,'mg')
      val = 1;
      set(hextract,'Value',val);
      % instruct no manual extraction with regtyp = 'mg'
      txt = ['Manual peak extraction not permitted when fitting to '....
        'multiple gaussian components' ];
      infodisplay(txt,analysis);
      beep, pause(2)
    end
    
    switch str{val};
      case 'Auto'
        analysis.extract = 'auto';
        set([htolextract,htolextracttxt,hrangeAL,hrangeAU,hrangeAtxt,...
          hrangeA2U,hrangeA2L,hrangeAtxt2,hrangeAtxt3,hrangeAtxt4],...
          'Visible','on');
        set([hnumberextract,hnumberextracttxt],'Visible','off');
        
      case 'Manual'
        analysis.extract = 'manual';
        set([hnumberextract,hnumberextracttxt],'Visible','on');
        set([htolextract,htolextracttxt,hrangeAL,hrangeAU,hrangeAtxt,...
          hrangeA2U,hrangeA2L,hrangeAtxt2,hrangeAtxt3,hrangeAtxt4],...
          'Visible','off');
        
    end
    updategraphics();
    [output.Av, output.Fv, output.Tv, analysis] = ...
      peakprocess(output.S,output.P,output.R,output.Fv,output.Tv,output.Av, ...
      output.Er,data,fitting,analysis);
    updategraphics();
  end

  function B1fit_callback(~,~)
    if get(hB1fit,'Value')
      fitting.B1fit = 'y'; fitting.theta_vector = [];
      if fitting.numbertheta == 1
        fitting.numbertheta = [];
        fitting.rangetheta = [];
      else
        xL = sscanf(get(hrangethetaL,'String'),'%d');
        xU = sscanf(get(hrangethetaU,'String'),'%d');
        fitting.rangetheta = [xL,xU];
        fitting.numbertheta = sscanf(get(hnumbertheta,'String'),'%d');
      end
      set(hB1fittxt,'String',['B1 fit = ', num2str(fitting.theta,4),'°'], ...
        'ForegroundColor',[0 0 0])
      set([hnumberthetatxt,hnumbertheta,hrangethetatxt,hrangethetaL,...
        hrangethetaU],'ForegroundColor',[0 0 0])
    else
      fitting.B1fit = 'n'; fitting.theta_vector = [];
      set(hB1fittxt,'String','B1 fit','ForegroundColor',[0.5 0.5 0.5])
      set([hnumberthetatxt,hnumbertheta,hrangethetatxt,hrangethetaL,...
        hrangethetaU],'ForegroundColor',[0.5 0.5 0.5])
    end
    
    [~,fitting,analysis] = parseinput(data,fitting,analysis);
    updategraphics();
    [output,fitting,analysis] = fittingloop(data,fitting,analysis);
    updategraphics();
    
  end

  function regtyp_callback(source,~)
    % Determine the selected data set.
    str = get(source,'String');
    val = get(source,'Value');
    % Set current data to the selected data set.
    switch str{val};
      case 'None'
        fitting.regtyp = 'none';
      case 'Min Energy'
        fitting.regtyp = 'me';
      case 'Min Curvature'
        fitting.regtyp = 'mc';
      case 'Uniform Penalty'
        fitting.regtyp = 'upen';
      case 'Multi Gauss'
        fitting.regtyp = 'mg';
    end
    updategraphics();
    [~,fitting,analysis] = parseinput(data,fitting,analysis);
    [output,fitting,analysis] = fittingloop(data,fitting,analysis);
    updategraphics();
  end

  function regadj_callback(source,~)
    % Determine the selected data set.
    str = get(source,'String');
    val = get(source,'Value');
    % Set current data to the selected data set.
    switch str{val};
      case 'Manual'
        fitting.regadj = 'manual';
      case 'Cross Validation'
        fitting.regadj = 'gcv';
      case '% MSE Increase'
        fitting.regadj = 'erinc';
      case 'None'
        fitting.regadj = 'none';
    end
    
    updategraphics();
    [~,fitting,analysis] = parseinput(data,fitting,analysis);
    [output,fitting,analysis] = fittingloop(data,fitting,analysis);
    updategraphics();
    
  end

  function slider_callback(~,~)
    fitting.numberT = get(hnumberT,'Value');
    fitting.rangeT(1) = get(hrangeTL,'Value');
    fitting.rangeT(2) = get(hrangeTU,'Value');
    fitting.percentinc = get(herinc,'Value');
    fitting.widthgauss = get(hwidthgauss,'Value');
    fitting.numbergauss = str2double(get(hnumbergauss,'String'));
    fitting.regweight = get(hregweight,'Value');
    fitting.numberT0 = str2double(get(hnumberT0,'String'));
    
    set(hnumberTtxt,...
      'String',['#Time Const = ',int2str(fitting.numberT)]);
    set(hrangeTLtxt,...
      'String',['Min T = ',num2str(fitting.rangeT(1)*1e3,3),' ms']);
    set(hrangeTUtxt,...
      'String',['Max T = ',num2str(fitting.rangeT(2),3),' s']);
    set(herinctxt,'String', ...
      ['% MSE Incr = ',num2str(fitting.percentinc,2)]);
    set(hwidthgausstxt,'String', ...
      ['Width = ',num2str(fitting.widthgauss,2),'%']);
    set(hregweighttxt,'String', ...
      ['Reg Weight = ',num2str(fitting.regweight,2)]);
    
    if strcmpi(fitting.twoD,'y')
      fitting.numberT2 = get(hnumberT2,'Value');
      fitting.rangeT2(1) = get(hrangeT2L,'Value');
      fitting.rangeT2(2) = get(hrangeT2U,'Value');
      set(hnumberT2txt,...
        'String',['#Time Const = ',int2str(fitting.numberT2)]);
      set(hrangeT2Ltxt, ...
        'String',['Min T = ',num2str(fitting.rangeT2(1)*1e3,3),' ms']);
      set(hrangeT2Utxt, ...
        'String',['Max T = ',num2str(fitting.rangeT2(2),3),' s']);
    end
    
    fitting = rmfield(fitting,{'T','T2'});
    
    [~,fitting,analysis] = parseinput(data,fitting,analysis);
    [output,fitting,analysis] = fittingloop(data,fitting,analysis);
    updategraphics();
  end

  function exit_callback(~,~)
    disp('Goodbye')
    uiresume(f);
  end

end
