function Test(N, vec, snr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N                %length time series
% vec              % Martin weight matrix
%                  % Identity weight martix
% snr              %SNR (in dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed = 123;  
rng(seed);
%%%%%%%%%%%%%%%
%Time series 1
K0a = 2;                                 %no. of sinusoids
beta0a = [2 2]';                              %amplitude
% f0a = round(N*0.05)/N;                   %frequency
f0a = round(N*[0.05 0.42]')/N;      %frequency
phi0a = zeros(K0a,1);                    %phase
%AR noise
a0a = [1 -0.85]';
b0a = 1;
%%%%%%%%%%%%%%%
%Time series 2

K0b = 1;                                 %no. of sinusoids
beta0b = 2;                              %amplitude
f0b = round(N*0.40)/N;                   %frequency
phi0b = zeros(K0b,1);                    %phase

a0b = 1;
% b0b = 1;
b0b = [1 -0.85]';
%%%%%%%%%%%%%%%

snr0 = snr;
fname = strcat('./N',num2str(N),'SNR',num2str(snr0), '_', vec, '/');
if isfolder(fname)==0
    mkdir(fname);
end
[Phi_signala, Phi_noisea, c_pa, s20a] = all(N,beta0a,f0a,b0a,a0a,snr0);
plot_spectra(Phi_signala+Phi_noisea, c_pa, 1);
[Phi_signalb, Phi_noiseb, c_pb, s20b] = all(N,beta0b,f0b,b0b,a0b,snr0);
plot_spectra(Phi_signalb+Phi_noiseb, c_pb, 2);

% Weight Matrix
if strcmp(vec, 'Martin')
    vecw = 0:N/2;
else
    vecw = [0 ones(1, N/2)];
end

DistTrue = comp_dist(c_pa,c_pb,vecw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W = diag((1:N/2));
% Sigma = (1/N)*(pi^2/6)*[eye(N/2-1) zeros(N/2-1,1); zeros(1,N/2-1) 2];
% mu = c_pa(2:N/2+1)-c_pb(2:N/2+1);
% 
% 2*trace((W*Sigma)^2) + 4*mu'*W*Sigma*W*mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nr = 5000;
Ymata = Inf(N,Nr);
Ymatb = Inf(N,Nr);

for run=1:Nr
    Ymata(:,run) = geny(N,beta0a,f0a,phi0a,b0a,a0a,s20a);
    Ymatb(:,run) = geny(N,beta0b,f0b,phi0b,b0b,a0b,s20b);
end

fname1 = strcat(fname,'simdata.mat');
save(fname1,'Ymata','Ymatb','DistTrue');

windowed_per(fname, fname1, vecw)

nulling_cep(fname, fname1, vecw);

FDRFER_cep(fname, fname1, vecw, 0.01);

FDRFER_cep(fname, fname1, vecw, 0.05);
end %function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = geny(N,beta,f,phi,b,a,s2)

t = (0:N-1)';
x = cos(2*pi*kron(t,f') + kron(ones(size(t)),phi'))*beta;

e = filter(b,a,sqrt(s2)*randn(2000+N,1));
e = e(end-N+1:end);

y = x+e;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s20, snr] = comp_s20(snr0,beta0,f0,b0,a0)
%Computes the noise variance s20 from the min. snr snr0
%The output snr contains the local snr for all the harmonics

om0 = 2*pi*f0;
snr1 = (beta0.^2/2)'./(comp_spec(b0',om0')./comp_spec(a0',om0'));
[m, ~] = min(snr1);

snr00 = 10^(snr0/10);
s20 = m/snr00;
snr = 10*log10(snr1/s20)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  sp=comp_spec(a0,om0)
%Auxiliary function for computing spectrum

i = sqrt(-1);
mm = (0:length(a0)-1)';
sp = abs(a0*exp(-i*kron(mm,om0))).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Phi_noise = comp_Phi_noise(A,B,s2,N)
%input:
%   A        = AR part of the model
%   B        = MA part of the model
%   s2       = variance driven noise
%   N        = no. of measurements
%output:
%   Phi_p    = true spectrum noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = A(:);
B = B(:);
freq = (0:N-1)/N;

i = sqrt(-1);
%AR part
kk = length(A)-1;
Exa = kron(freq',(0:kk));
Exa = exp(-i*2*pi*Exa);

%MA part
kk = length(B)-1;
Exb = kron(freq',(0:kk));
Exb = exp(-i*2*pi*Exb);

%true spectrum
Phi_noise = s2.*(abs(Exb*B).^2)./(abs(Exa*A).^2); 
Phi_noise = Phi_noise(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Phi_signal = comp_Phi_signal(f0,beta0,N)
%input:
%   beta0  = amplitude
%   f0     = frequency
%   N      = no. of measurements
%output:
%   Phi_p  = true spectrum signal

freq = (0:N-1)/N;

Phi_signal = zeros(size(freq));
for ind=1:length(f0)
    pos = find(freq==f0(ind));
    if length(pos)~= 1
        fprintf('Error: computations for line spectra');
        Phi_signal = [];
        return
    else
        Phi_signal(pos) = beta0(ind)^2/2;
        Phi_signal(N-pos+2) = Phi_signal(pos);
    end
end
Phi_signal = Phi_signal(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_p = comp_true_cep(Phi_p) 
%input:
%   Phi_p  = true spectrum
%output:
%   c_p    = true cepstral coeff

c_p = ifft(log(Phi_p));
if max(abs(imag(c_p)))<10^(-10)
    c_p = real(c_p);
else
    fprintf('ERROR: in the evaluation of c_p.\n');
    c_p = [];
    return
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi_signal, Phi_noise, c_p, s20] = all(N,beta0,f0,b0,a0,snr0)

    %true spectrum signal
    Phi_signal = comp_Phi_signal(f0,beta0,N);

    %variance driven noise 
    [s20, ~] = comp_s20(snr0,beta0,f0,b0,a0);

    %true spectrum noise
    Phi_noise = comp_Phi_noise(a0,b0,s20,N);

    %true cepstral coeff.
    c_p = comp_true_cep(Phi_signal+Phi_noise);

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_spectra(Phi, c_p, ts)
%input:
%ts:    1 for TS1, 2 for TS2

    N = length(Phi);
    figure(ts)
    subplot(2,1,1)
    semilogy(0:N-1,Phi,'b-')
    tit = sprintf('TS%i: Log Spectrum',ts);
    title(tit);
    subplot(2,1,2)
    semilogy(1:N-1,abs(c_p(2:end)),'r-') %!!!!!!!!
    tit = sprintf('TS%i: Magnitudes Cepstral Coefficients',ts);
    title(tit);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = comp_dist(cepa,cepb,vecw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:

N = length(cepa);
if length(cepb)~=N
    fprintf('Different no. of cepstral coefficients for the two time series.\n');
    dist = [];
    return;
end
cepa = cepa(:);
cepb = cepb(:);

diff = abs(cepa-cepb);


dist = sqrt( vecw*(diff(1:N/2+1).^2) );

%dist = vecw*(diff(1:N/2+1).^2);

end %function





