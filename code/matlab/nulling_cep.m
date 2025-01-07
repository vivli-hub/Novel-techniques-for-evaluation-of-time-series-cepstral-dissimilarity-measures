function nulling_cep(fname, vecw, dis)

%input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fname           % folder name
% vecw             % Weight martix
% dis              % sq-squared Eulcidean distance
%                  % eu-squared Eulcidean distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the simulated data and 'true' distance
fname1 = strcat(fname,'simdata.mat');
load(fname1,'Ymata', 'Ymatb', 'DistTrue');
%Nr-No. of runs
%N-length of time series
[N, Nr] = size(Ymata);

%Calculate the estimated cepstral ceofficients
[CEPa, vec_KSFa] = comp_CEP(Ymata, N);
[CEPb, vec_KSFb] = comp_CEP(Ymatb, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculated the estimated cepstral distance between the time series
MatDist = Inf*ones(4,Nr);
for met=1:4
    for rr=1:Nr
        MatDist(met,rr) = comp_dist(CEPa(met,:,rr),CEPb(met,:,rr),vecw, dis);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResCepNulling = cell(4,3);
ResCepNulling{1,1} = 'log_Periodogram';
ResCepNulling{2,1} = 'BIC';
ResCepNulling{3,1} = 'KSF';
ResCepNulling{4,1} = 'MRI';

MD = MatDist-DistTrue;
for ind=1:4
    ResCepNulling{ind,2} = mean(MD(ind,:));
    ResCepNulling{ind,3} = var(MD(ind,:));
end

%Find out the bias of the estimated distance (mean and std.)
fname2 = strcat(fname,'resultsCN.mat');
save(fname2, 'CEPa', 'vec_KSFa', 'CEPb', 'vec_KSFb', 'MatDist',"ResCepNulling");

end %function



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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_phi, c_e, Phi_PER_mod] = comp_periodogram(y)
%input:
%   y = data vector
%output
%   Phi_PER_mod = periodogram

    N = length(y);
    L = N;
    %v = ones(size(y));
    Phi_PER = periodogram(y, ones(1, N), (0:2*pi/N:2*pi*(N-1)/N));
    log_phi = log(Phi_PER);
    Phi_PER_mod = flipud([Phi_PER(L/2+1:L); Phi_PER(1:L/2)]);
    Phi_PER_mod = Phi_PER_mod(:);
    c_e = ifft(log(Phi_PER));
    c_e(1) = c_e(1) + 0.57721; % Euler constant
    %check the estimated cepstral coefficients, keep symmetry
    if max(abs( c_e(2:N/2)-c_e(N:-1:N/2+2) ))>10^(-4)
        fprintf('ERROR in the estimated cepstral coeff.!\n');
        return;
    else
       c_e(N:-1:N/2+2) = c_e(2:N/2); 
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
    [log_phi, c_e, ~]    = comp_periodogram(Ymat(1:N,rr));
    CEP(1,:,rr) = log_phi; %log-periodogram
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
