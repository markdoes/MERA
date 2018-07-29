function [data] = maketestdata1D(T2,s,t,theta,Ntrials,SNR)

rng('shuffle')

N = length(t);
TE = diff(t(1:2));

noise_sd = 1/SNR;
noise = randn(N,Ntrials)*noise_sd;


K = size(T2,1);
A = zeros(N,K);

% build A-matrix for a given theta
if theta == 180
  A(:,1:K) = exp(-t*(1./T2'));
else
  for kx = 1:K
    A(:,kx)=EPGdecaycurve(theta,T2(kx),TE,N);
  end
end

yx = A*s;
data.t = t;
data.D = bsxfun(@plus,yx,noise);

end



function EchoAmp = EPGdecaycurve(theta,T2,TE,N)
% Assumes a CPMG condition
% Arbitrarily set T1 = 1 ... should work well for most situations, but I might
% want to make this a variable later.
T1 = 1;

% Compute rotation matrix in terms of coherence state
Rot=[cosd(theta/2)^2,sind(theta/2)^2,-1i*sind(theta);...
  sind(theta/2)^2,cosd(theta/2)^2,1i*sind(theta);...
  -0.5i*sind(theta),0.5i*sind(theta),cosd(theta)];

% Initialize vector to track echo amplitude
EchoAmp=zeros(N,1);

% Initialize magnetization phase state vector (MPSV) and set all
% magnetization in the F1 state.
M=zeros(3*N,1);
M(1,1)=exp(-(TE/2)/T2);

% Compute relaxation and transition matrix, E
% same as T*E in Prasloski et al.

RelTrans=spalloc(3*N,3*N,3*N); % 3N x 3N, with 3N entries
E2 = exp(-TE/T2);
E1 = exp(-TE/T1);
RelTrans(1,2)=E2;
for ke=1:N-1
  RelTrans(3*ke-1,3*ke+2)=E2;
  RelTrans(3*ke+1,3*ke-2)=E2;
  RelTrans(3*ke,3*ke)=E1;
end

BigRot=spalloc(3*N,3*N,9*N); % 3N x 3N, with 9N entries
for kn=1:N
  BigRot(3*kn-2:3*kn,3*kn-2:3*kn)=Rot;
end

% Apply first refocusing pulse and get first echo amplitude
M(1:3)=Rot*M(1:3);
EchoAmp(1)=abs(M(2,1))*exp(-(TE/2)/T2);

% Apply relaxation matrix
M=RelTrans*M;
% Perform flip-relax sequence ETL-1 times
for ke=2:N
  % Perform the flip
  M=BigRot*M;
  % Record the magnitude of the population of F1* as the echo amplitude
  % and allow for relaxation
  EchoAmp(ke)=abs(M(2,1))*exp(-(TE/2)/T2);
  % Allow time evolution of magnetization between pulses
  M=RelTrans*M;
end

end

