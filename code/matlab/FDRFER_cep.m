function FDRFER_cep(fname, fname1, vecw, alpha)

%input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fname           % folder name
% vecw            % Weight martix
% alpha           % pre-speciﬁed FDR or FER value 
% dis             % sq-squared Eulcidean distance
%                 % eu-squared Eulcidean distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fname1,'Ymata', 'Ymatb', 'DistTrue');
[N, Nr] = size(Ymata);

CEPa = comp_CEPf(Ymata, alpha);
CEPb = comp_CEPf(Ymatb, alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MatDist = Inf*ones(2,Nr);
for met=1:2
    for rr=1:Nr
        MatDist(met,rr) = comp_dist(CEPa(met,:,rr),CEPb(met,:,rr),vecw);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResCepNulling = cell(2,3);
ResCepNulling{1,1} = strcat('FDR',num2str(alpha));
ResCepNulling{2,1} = strcat('FER',num2str(alpha));

MD = MatDist-DistTrue;
for ind=1:2
    ResCepNulling{ind,2} = mean(MD(ind,:));
    ResCepNulling{ind,3} = var(MD(ind,:));
end

fname2 = strcat(fname, num2str(alpha),'resultsFDR.mat');
save(fname2, 'CEPa',  'CEPb',  'MatDist',"ResCepNulling");

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CEP = comp_CEPf(ytr, alpha)
% compute the cepstral ceofficients
N = size(ytr,1);
Nr = size(ytr, 2);
ce = zeros(N, Nr);
for i = 1:Nr
    c_e = comp_periodogram(ytr(:, i));
    ce(:,i) = c_e;
end
%Cepstral nulling
% Nr = Number of runs
[N, Nr] = size(ce);
M       = N/2+1;
pv      = [1 - (alpha.*(1:M)'./M ) alpha./(M+1-(1:M)')];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_e, Phi_PER_mod] = comp_periodogram(y)
%input:
%   y = data vector
%output
%   Phi_PER_mod = periodogram

    N = length(y);function FDRFER_cep(fname, vecw, alpha, dis)

%Load the simulated data and 'true' distance
fname1 = strcat(fname,'simdata.mat');
load(fname1,'Ymata', 'Ymatb', 'DistTrue');
%Nr-No. of runs
%N-length of time series
[N, Nr] = size(Ymata);

%Calculate the estimated cepstral ceofficients
CEPa = comp_CEPf(Ymata, alpha);
CEPb = comp_CEPf(Ymatb, alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculated the estimated cepstral distance between the time series
MatDist = Inf*ones(2,Nr);
for met=1:2
    for rr=1:Nr
        MatDist(met,rr) = comp_dist(CEPa(met,:,rr),CEPb(met,:,rr),vecw,dis);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResCepNulling = cell(2,3);
ResCepNulling{1,1} = strcat('FDR',num2str(alpha));
ResCepNulling{2,1} = strcat('FER',num2str(alpha));

MD = MatDist-DistTrue;
for ind=1:2
    ResCepNulling{ind,2} = mean(MD(ind,:));
    ResCepNulling{ind,3} = var(MD(ind,:));
end

%Find out the bias of the estimated distance (mean and std.)
fname2 = strcat(fname, num2str(alpha),'resultsFDR.mat');
save(fname2, 'CEPa',  'CEPb',  'MatDist',"ResCepNulling");

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CEP = comp_CEPf(ytr, alpha)
% compute the cepstral ceofficients
N = size(ytr,1);
Nr = size(ytr, 2);
ce = zeros(N, Nr);
for i = 1:Nr
    c_e = comp_periodogram(ytr(:, i));
    ce(:,i) = c_e;
end
%Cepstral nulling
% Nr = Number of runs
[N, Nr] = size(ce);
M       = N/2+1;
pv      = [1 - (alpha.*(1:M)'./M ) alpha./(M+1-(1:M)')];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_e, Phi_PER_mod] = comp_periodogram(y)
%input:
%   y = data vector
%output
%   Phi_PER_mod = periodogram

    N = length(y);
    L = N;
    v = ones(size(y));
    Phi_PER = periodogram(y, ones(1, N), (0:2*pi/N:2*pi*(N-1)/N));
    Phi_PER_mod = flipud([Phi_PER(L/2+1:L); Phi_PER(1:L/2)]);
    Phi_PER_mod = Phi_PER_mod(:);
    c_e = ifft(log(Phi_PER));
    c_e(1) = c_e(1) + 0.57721; % Euler constant
    if max(abs( c_e(2:N/2)-c_e(N:-1:N/2+2) ))>10^(-4)
        fprintf('ERROR in the estimated cepstral coeff.!\n');
        return;
    else
       c_e(N:-1:N/2+2) = c_e(2:N/2); 
    end
end

end %function



