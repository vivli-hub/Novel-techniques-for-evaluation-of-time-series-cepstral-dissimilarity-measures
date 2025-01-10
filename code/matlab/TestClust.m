function TestClust(N, NoTS, snr0, dist, flag_plot, clust)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N                 %length time series
% NoTS              %no. of TS in each cluster
% snr0              %SNR (in dB)
% dist              %'Martin' or 'Euclidean'
% flag_plot         %1=plot clusters 0=no plot
% clust             %The method of cluster
%                   %'K-medoids_eu'= K-medoids with the Euclidean distance
%                   %'K-medoids_sq'= K-medoids with the squared Euclidean distance
%                   %'K-means'= K-means with the squared Euclidean distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set seed
seed = 123;  
rng(seed);
%Set up a folder to store useful results.
fname = strcat('./N',num2str(N),'SNR',num2str(snr0),'noTs',num2str(NoTS),'met', clust, '/');
if isfolder(fname)==0
    mkdir(fname);
end
%Generate the simulated data of the two clusters
%Param. cluster a
K0a = 2;                    %no. of sinusoids
beta0a = [2 2]';            %amplitude
phi0a = zeros(K0a,1); %phase
%AR noise
a0a = [1 -0.85]';
b0a = 1;


%Param. cluster b
K0b = 1;                     %no. of sinusoids
beta0b = 2;                  %amplitude
phi0b = zeros(K0b,1);        %phase

%MA noise
a0b = 1;
b0b = [1 -0.85]';

for i=1:NoTS
    f0a = [0.02 + (0.08-0.02)*rand(1);  0.39 + (0.45-0.39)*rand(1)];
    [s20a, ~] = comp_s20(snr0,beta0a,round(N*f0a)/N,b0a,a0a);
    Ymata(:,i) = geny(N,beta0a,f0a,phi0a,b0a,a0a,s20a);
    f0b = 0.37 + (0.44-0.37)*rand(1); 
    [s20b, ~] = comp_s20(snr0,beta0b,round(N*f0b)/N,b0b,a0b);
    Ymatb(:,i) = geny(N,beta0b,f0b,phi0b,b0b,a0b,s20b);
end
%label the simulated data
labels(1:NoTS) = {'A'};
labels(NoTS+1:2*NoTS)= {'B'};
true_label = {1:NoTS, NoTS+1:2*NoTS};
%Save the simulated data in 'simdata.mat'
fname1 = strcat(fname,'simdata.mat');
save(fname1,'Ymata','Ymatb');

% Two kinds of distance between the time series
if strcmp(dist, 'Euclidean') 
    vecw = [0 ones(1,N/2)];
elseif strcmp(dist, 'Martin')
    vecw = (0:N/2); 
else
    fprintf('Error: Distance\n');
    return
end
%Calculate the cepstral coefficients applying the window periodogram and
%cepstral nulling 
[CEP_WP, CEP_nulling] = comp_exp_CEP(fname1);
%shuffled_labels = labels(random_indices);
%true_label = cell(1, 2);
%true_label{1} = find(strcmp(shuffled_labels, 'A')); 
%true_label{2} = find(strcmp(shuffled_labels, 'B'));
CEP = [CEP_WP; CEP_nulling];
nulling_cep_clust(fname, CEP, vecw, labels, flag_plot, NoTS);
nonzero_nulling = squeeze(sum(CEP_nulling(:, 1:N/2+1 , :) ~=0, 2));
fname2 = strcat(fname,'CEP.mat');
save(fname2, 'CEP_WP', 'CEP_nulling', 'nonzero_nulling');

[silhouette_ID, sim_ID] = sim_silhouette(CEP_WP, CEP_nulling, true_label, 1, 2, clust);
[silhouette_Martin, sim_martin] = sim_silhouette(CEP_WP, CEP_nulling, true_label, 2, 2, clust);
fname3 = strcat(fname,'results.mat');
save(fname3,'silhouette_ID','sim_ID','silhouette_Martin', 'sim_martin');
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
%Phi_p: true spectrum
%c_p:   true cepstral coeff
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate cepstral coefficients with windowed periodogram and cepstral
% nulling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%   fname1 : the name of the address which save the simulated data
%   
%Output:
%   CEP_WP   = cepstral ceofficients from order of 1,...,N with windowed
%   periodogram
%   CEP_nulling = cepstral ceofficients from order of 1,...,N applied
%   cepstral nulling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CEP_WP, CEP_nulling] = comp_exp_CEP(fname1)
load(fname1,'Ymata', 'Ymatb');
% Put all of the simulated the data together and shuffle the data
y = [Ymata Ymatb];
%random_indices = randperm(size(y, 2));
%y = y(:, random_indices); 
% compute the cepstral coefficients
Nr = size(y,2);
N = size(y,1);
% calculate the cepstral coefficients applied the windowed periodogram
% hanning and hamming
ytr = y(1:N,:);
CEP_WP = comp_CEP_WP(ytr);
% compute the cepstral coefficients without any windowed periodogram
% apply the cepstral nulling
% Thresholds including periodogram, BIC, KSF, MRI
N_new = size(ytr, 1);
CEP_nulling = Inf(8, N_new, Nr);
[CEP, ~] = comp_CEP(ytr, N_new);
CEP_nulling(1:4,:,:) = CEP;
% Thresholds including FDR and FER, respectively for alpha = 0.01 and alpha
% = 0.05
% FDR_0.01 and FER_0.01
    CEP = comp_CEPf(ytr, 0.01);
    CEP_nulling(5:6,:,:) = CEP;
% FDR_0.05 and FER_0.05
    CEP = comp_CEPf(ytr, 0.05);
    CEP_nulling(7:8,:,:) = CEP;

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the cepstral coefficients with the window periodogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CEP = comp_CEP_WP(Y)

[N, Nr] = size(Y);
CEP     = Inf(2,N,Nr);

for run=1:Nr
    [CEP(1,:,run), ~] = wp( Y(:,run), 'rect' );
    [CEP(2,:,run), ~] = wp( Y(:,run), 'hann' );
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_p, Phi_WP_mod] = wp(y, v)



 N = length(y);

 if strcmp(v, 'rect')
     [Phi_WP, freq] = periodogram(y, ones(1, N), (0:2*pi/N:2*pi*(N-1)/N));
 elseif strcmp(v, 'hann')
     [Phi_WP, freq] = periodogram(y, hann(N), (0:2*pi/N:2*pi*(N-1)/N));
 else
        fprintf('error \n ')
 end
 
 L = N;
Phi_WP_mod = flipud([Phi_WP(L/2+1:L); Phi_WP(1:L/2)]);
 c_p = comp_cep_wp(Phi_WP);

end %funtion



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_p = comp_cep_wp(Phi_p) 
%input:
%   Phi_p  = spectrum
%output:
%   c_p    = cepstral coeff

%checking
c_p = ifft(log(Phi_p));
if max(abs(imag(c_p)))<10^(-4)
    c_p = real(c_p);
else
    fprintf('ERROR: in the evaluation of c_p.\n');
    c_p = [];
    return
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_phi,c_e, Phi_PER_mod] = comp_periodogram(y)
%input:
%   y = data vector
%output
%   Phi_PER_mod = periodogram
%   log_phi = log periodogram
%   c_e = cepstral coefficients

    N = length(y);
    L = N;
    Phi_PER = periodogram(y, ones(1, N), (0:2*pi/N:2*pi*(N-1)/N));
    Phi_PER_mod = flipud([Phi_PER(L/2+1:L); Phi_PER(1:L/2)]);
    Phi_PER_mod = Phi_PER_mod(:);
    log_phi = log(Phi_PER);
    c_e = ifft(log(Phi_PER));
    c_e(1) = c_e(1) + 0.57721; % Euler constant
    if max(abs( c_e(2:N/2)-c_e(N:-1:N/2+2) ))>10^(-4)
        fprintf('ERROR in the estimated cepstral coeff.!\n');
        return;
    else
       c_e(N:-1:N/2+2) = c_e(2:N/2); 
    end
    if max(abs(imag(c_e)))<10^(-4)
    c_e = real(c_e);
else
    fprintf('ERROR: in the evaluation of c_p.\n');
    c_e = [];
    return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CEP, vec_KSF] = comp_CEP(Ymat, N)

% Nr = Number of runs
Nr      = size(Ymat,2);
M       = N/2+1;
Lk      = comp_Lk(M);
%Grid for \mu
grid_mu  = (1:0.1:10);

vec_KSF = Inf(Nr,1);
CEP     = Inf(4,N,Nr);

vec_mu = [    1+sqrt(log(M))       % BIC
              Inf                  % NULL KSF
              1+sqrt(2*log(M))     % MRI
              ];

for rr=1:Nr
    [log_phi,c_e, ~]    = comp_periodogram(Ymat(1:N,rr));
    CEP(1,:,rr) = log_phi; %log periodogram
    vec_KSF(rr) = comp_mu_KSF(c_e',N,M,Lk,grid_mu,0,0);
    vec_mu(2)   = vec_KSF(rr);
    for met=1:3
        CEP(met+1,:,rr) = comp_c_til(c_e,N,vec_mu(met));
    end %met
end %rr

end %function

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_til = comp_c_til(c_e,N,mu)
% Turns to zero the cepstral coeff. 
% which are deemed to be small
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%   c_e     = estimated cepstral coeff.
%   N       = no. of measurements
%   mu      = threshold
%output:
%   c_til   = new estimates of the cepstral coeff.
%             \check c <-> c_til
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cc = c_e;

if abs(cc(1)) < mu*pi/sqrt(3*N)
  cc(1) = 0;
end
for k=2:N/2
  if abs(cc(k)) < mu*pi/sqrt(6*N)
    cc(k) = 0;
  end
end
if abs(cc(N/2+1)) < mu*pi/sqrt(3*N)
  cc(N/2+1) = 0;
end
for k=N/2+2:N
  if abs(cc(k)) < mu*pi/sqrt(6*N)
    cc(k) = 0;
  end
end

c_til = cc;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CEP = comp_CEPf(ytr, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%   ytr     = data vector
%   alpha   = pre-speciﬁed FDR or FER value
%output:
%   CEP   = new estimates of the cepstral coeff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the cepstral ceofficients
N = size(ytr,1);
Nr = size(ytr, 2);
ce = zeros(N, Nr);
for i = 1:Nr
    [~,c_e,~] = comp_periodogram(ytr(:, i));
    ce(:,i) = c_e;
end
%Cepstral nulling
% Nr = Number of runs
[N, Nr] = size(ce);
M       = N/2+1;
pv      = [1-(alpha.*(1:M)'./M ) alpha./(M+1-(1:M)')];
qv      = norminv(pv);

CEP     = Inf(2,N,Nr);

for rr=1:Nr
    %first M entries
    cor = ce(1:M, rr);
    c = cor;
    c(2:end-1) = abs(cor(2:end-1))./(pi/sqrt(6)/sqrt(N));
    c([1 end]) = abs(cor([1 end]))./(pi/sqrt(3)/sqrt(N));
    [cs, ind] = sort(c,'descend');

    for met=1:2
        comp = (cs>=qv(:,met));
        pk = find(comp==1, 1, 'last');
        if numel(pk)>0
            cor(ind(pk+1:M)) = 0;
        end
        CEP(met,1:M,rr) = cor;
    end
    
end %rr
end

% The function computes \mu_{KSF} by using the KSF criterion
function mu_KSF = comp_mu_KSF(c_e,N,M,Lk,vec_mu,ex,flag_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%   c_e       = estimated cepstral coeff.
%   N         = no. of measurements
%   M         = M = N/2+1
%   Lk        = code length for the structure
%   vec_mu    = grid for \mu
%   ex        = no. of the current experiment (used only for plots)
%   flag_plot = if flag_plot=2, then the following plots are generated:
%               L versus \mu_{KSF} and k versus \mu_{KSF}
%output:
%   mu_KSF    = KSF threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = pi/sqrt(N)*[1/sqrt(3); 1/sqrt(6)*ones(N/2-1,1); 1/sqrt(3)];
x = c_e(1:M)./s;
[ax_ord,ind] = sort(abs(x),'descend');
x_ord = x(ind);

%\ell <-> h_slen
h_slen = 0;
mu_KSF = 1 + h_slen;
%For h_slen = 0, KSF function is Inf !!!!
Lbest = Inf;

if flag_plot==2,
     L_vec = [];
     mu_KSF_vec = [];
     k_vec = [];
end

vec_mu1 = vec_mu-1;
vec_h = vec_mu1(2:end);
for h_slen=vec_h,
    
    ff = find( ax_ord >= (1+h_slen) );
    k = length(ff);
    if k>0,
        z = x_ord(1:k);
        az = ax_ord(1:k);
        tilde_z = (1+2*h_slen) + floor((az-1-h_slen)./(2*h_slen))*(2*h_slen);
        tilde_z = sign(z).*tilde_z;
        %centroid
        qx_ord = [tilde_z; zeros(M-k,1)];
        %code length for \tilde z 
        nt = sum(tilde_z.^2);
        logCk = (k/2)*log(nt/k) + k/2 + log(k)/2;
        logQk = -(k/2)*log(2*pi) + k*log(2*h_slen);
        %code length for the signal
        L_signal = Lk(k) + logCk - logQk;
        %code length for noise (distortion)
        L_noise = sum((x_ord-qx_ord).^2)/2;
        %KSF
        L = L_signal + L_noise;
 
        if L<Lbest,
            Lbest = L;
            mu_KSF = 1 + h_slen;
        end       
        if flag_plot==2,
            L_vec = [L_vec; L];
            mu_KSF_vec = [ mu_KSF_vec; (1 + h_slen)];
            k_vec = [k_vec; k];
        end       
    else
        % h_len is too large; all coeff. are quantized to zero
        break;
    end %k>0
end % for

if flag_plot==2,
    plot_L_k_versus_mu(L_vec,k_vec,mu_KSF_vec,ex,N);
end

end

% Function used by cep.m. Computes the code length for the structure
function Lk = comp_Lk(M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%   M       = N/2+1 (N=no. of measurements)
%output:
%   Lk      = vector which contains the code length for the structure
%             k\in\{1,\ldots,M-1\}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lk = zeros(M,1);
for k=1:M-1,
    temp = (M+1/2)*log(M)-(k+1/2)*log(k)-(M-k+1/2)*log(M-k) -1/2*log(2*pi);
    Lk(k) = min(M*log(2),temp+log(k)+log(1+log(M-1)));
end
Lk(M) = Inf; %at least one coeff. is turned to zero

end


function nulling_cep_clust(fname, CEP, vecw, labels, flag_plot, No)


Nr = size(CEP, 3);



MatDist = Inf(10,Nr,Nr);
for met=1:10
    for i=1:Nr
        for j=1:i-1
            MatDist(met,i,j) = comp_dist(CEP(met,:,i),CEP(met,:,j),vecw, 'eu');
        end
    end
end

methods = {
    'Rectangle Window'
    'Hanning'
    'log-Periodogram' 
    'BIC'
    'KSF'
    'MRI'
    'FDR_0.01'
    'FER_0.01'
    'FDR_0.05'
    'FER_0.05'};
              

for met=1:10
     MatDist(met,:,:) = transf_mat_dist( squeeze(MatDist(met,:,:)) );
    
    if flag_plot==1 %plot clusters
        figure(1) 
        YY=cmdscale(squeeze(MatDist(met,:,:)),2);
        subplot(5,2,met);
        plot(YY(1:No,1),YY(1:No,2),'ko', YY(No+1:2*No,1),YY(No+1:2*No,2),'r*');
        title(methods(met));
    end
end

fname2 = strcat(fname,'results.mat');
save(fname2, 'MatDist', "methods");

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = transf_mat_dist(MD)

D           = MD;
[I,J]       = find(MD==Inf);
if sum(I>J)>0
    fprintf('Inf-values in incorrect positions in the matrix of distances.\n');
end
D(D==Inf)   =0;
D           = D+D';
    
end %function
