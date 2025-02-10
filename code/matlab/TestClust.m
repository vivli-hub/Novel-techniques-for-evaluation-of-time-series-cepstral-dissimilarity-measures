function sim_all_met = TestClust(Nruns, N, NoTS, snr0, wmat, clust_dist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nruns             =no. runs experiment
% N                 =length time series
% NoTS              =no. of TS in each cluster
% snr0              =SNR (in dB)
% wmat              =matrix: 'Martin' or 'Identity'
% clust_dist        =dist clustering: 'euclidean' or 'sqEuclidean'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set seed
seed = 123;  
rng(seed);

%Set up a folder to store useful results
fname = strcat('./N',num2str(N),'SNR', num2str(snr0), ...
    'noTs',num2str(NoTS), ...
    '_', wmat, ...
    '_', clust_dist, '/');
fname1 = strcat(fname,'similarity.mat');
if isfolder(fname)==0
    mkdir(fname);
end
sim_all_met = []; 
save(fname1, 'sim_all_met');

%Run experiments
for ex=1:Nruns
    sim_all_met = exp_clust(fname1, N, NoTS, snr0, wmat, clust_dist);
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sim_all_met = exp_clust(fname1, N, NoTS, snr0, wmat, clust_dist)
 
%Read results previous experiments
load(fname1, 'sim_all_met');

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

Ymata = Inf(N,NoTS);
Ymatb = Inf(N,NoTS);
for i=1:NoTS
    f0a = [0.02 + (0.08-0.02)*rand(1);  0.39 + (0.45-0.39)*rand(1)];
    [s20a, ~] = comp_s20(snr0,beta0a,round(N*f0a)/N,b0a,a0a);
    Ymata(:,i) = geny(N,beta0a,f0a,phi0a,b0a,a0a,s20a);
    f0b = 0.37 + (0.44-0.37)*rand(1); 
    [s20b, ~] = comp_s20(snr0,beta0b,round(N*f0b)/N,b0b,a0b);
    Ymatb(:,i) = geny(N,beta0b,f0b,phi0b,b0b,a0b,s20b);
end
true_label = {(1:NoTS), (NoTS+1:2*NoTS)};

%Calculate the cepstral coefficients applying the windowed periodogram and
%cepstral nulling 
[CEP_WP, CEP_nulling] = comp_exp_CEP(Ymata, Ymatb);

% % Two kinds of distance between the time series
if strcmp(wmat, 'Identity') 
    vecw = ones(1,N/2);
elseif strcmp(wmat, 'Martin')
    vecw = sqrt(1:N/2); 
else
    fprintf('Error: Distance\n');
    return
end

sim_vector = Inf(1,9);
for ii=1:2
    sim_vector(ii) = clust_kmedoids(squeeze(CEP_WP(ii,:,:)), vecw, clust_dist, true_label);
end
for ii=1:7
    sim_vector(ii+2) = clust_kmedoids(squeeze(CEP_nulling(ii,:,:)), vecw, clust_dist, true_label);
end
sim_all_met = [sim_all_met; sim_vector];
save(fname1, 'sim_all_met');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate cepstral coefficients with windowed periodogram and cepstral
% nulling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%   fname1 : the name of the file which contains simulated data
%   
%Output:
%   CEP_WP      = cepstral ceofficients 1,...,N (windowed periodogram)
%   CEP_nulling = cepstral ceofficients 1,...,N (cepstral nulling) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CEP_WP, CEP_nulling] = comp_exp_CEP(Ymata, Ymatb)

% Put all of the simulated data together
ytr = [Ymata Ymatb];
[N, Nr] = size(ytr);

% Calculate the cepstral coefficients (windowed periodogram)
% rectangular and hanning
CEP_WP = comp_CEP_WP(ytr);

%Cepstral nulling
CEP_nulling = Inf(7, N, Nr);
[CEP, ~] = comp_CEP(ytr, N);
CEP_nulling(1:3,:,:) = CEP;
% FDR_0.01 and FER_0.01
CEP = comp_CEPf(ytr, 0.01);
CEP_nulling(4:5,:,:) = CEP;
% FDR_0.05 and FER_0.05
CEP = comp_CEPf(ytr, 0.05);
CEP_nulling(6:7,:,:) = CEP;

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the cepstral coefficients (windowed periodogram)
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
     [Phi_WP, ~] = periodogram(y, ones(1, N), (0:2*pi/N:2*pi*(N-1)/N));
 elseif strcmp(v, 'hann')
     [Phi_WP, ~] = periodogram(y, hann(N), (0:2*pi/N:2*pi*(N-1)/N));
 else
     fprintf('error \n ');
 end
 
 L = N;
 Phi_WP_mod = flipud([Phi_WP(L/2+1:L); Phi_WP(1:L/2)]);
 c_p = comp_cep_wp(Phi_WP);

end %function



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
CEP     = Inf(3,N,Nr);

vec_mu = [    1+sqrt(log(M))       % BIC
              Inf                  % NULL KSF
              1+sqrt(2*log(M))     % MRI
              ];

for rr=1:Nr
    [~,c_e, ~]      = comp_periodogram(Ymat(1:N,rr));
    vec_KSF(rr)     = comp_mu_KSF(c_e',N,M,Lk,grid_mu,0,0);
    vec_mu(2)       = vec_KSF(rr);
    for met=1:3
        CEP(met,:,rr) = comp_c_til(c_e,N,vec_mu(met));
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

if flag_plot==2
     L_vec = [];
     mu_KSF_vec = [];
     k_vec = [];
end

vec_mu1 = vec_mu-1;
vec_h = vec_mu1(2:end);
for h_slen=vec_h
    
    ff = find( ax_ord >= (1+h_slen) );
    k = length(ff);
    if k>0
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
 
        if L<Lbest
            Lbest = L;
            mu_KSF = 1 + h_slen;
        end       
        if flag_plot==2
            L_vec = [L_vec; L];
            mu_KSF_vec = [ mu_KSF_vec; (1 + h_slen)];
            k_vec = [k_vec; k];
        end       
    else
        % h_len is too large; all coeff. are quantized to zero
        break;
    end %k>0
end % for

if flag_plot==2
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
for k=1:M-1
    temp = (M+1/2)*log(M)-(k+1/2)*log(k)-(M-k+1/2)*log(M-k) -1/2*log(2*pi);
    Lk(k) = min(M*log(2),temp+log(k)+log(1+log(M-1)));
end
Lk(M) = Inf; %at least one coeff. is turned to zero

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function similarity = cluster_similarity(C, C_prime)
    k = length(C);
    similarity_sum = 0;
    
    for i = 1:k
        max_sim = 0;
        for j = 1:k
            sim = 2 * length(intersect(C{i}, C_prime{j})) / (length(C{i}) + length(C_prime{j}));
            if sim > max_sim
                max_sim = sim;
            end
        end
        similarity_sum = similarity_sum + max_sim;
    end
    similarity = similarity_sum / k;
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sim_index = clust_kmedoids(c_e, vecw, clust_dist, true_label)
%Input:
%c_e        = cep. coeff. for all time series
%vecw       = vector weights dist.
%clust_dist = dist clustering: 'euclidean' or 'seuclidean'
%true_label = true labels clusters
%Output:
%sim_index  = similarity index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, Nr] = size(c_e); 
D = Inf(N/2, Nr);
for i = 1:Nr
    D(:,i) = vecw'.*c_e(2:N/2+1,i);
end

noc = 2; %no. of clusters
[idx, ~] = kmedoids(D', noc, 'Distance',clust_dist,'replicates',3);
index_1 = find(idx(:) == 1); 
index_2 = find(idx(:) == 2);
C = {index_1, index_2};
sim_index = cluster_similarity(true_label, C);

end % function
