function windowed_per(fname, vecw, dis)
%input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fname            % folder name
% vecw             % Matrix
% dis              % 'sq' - squared distance
%                  % 'eu' - distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the simulated data
fname1 = strcat(fname,'simdata.mat');
load(fname1,'Ymata', 'Ymatb', 'DistTrue');

%Calculate the estimated cepstral coefficients
CEPaWP = comp_CEP_WP(Ymata);
CEPbWP = comp_CEP_WP(Ymatb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%No. of runs
Nr = size(Ymata,2);
%Calculated the estimated cepstral distance between the time series
MatDistWP = Inf*ones(2,Nr);
for met=1:2
    for rr=1:Nr
        MatDistWP(met,rr) = comp_dist(CEPaWP(met,:,rr),CEPbWP(met,:,rr),vecw, dis);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResWP = cell(2,3);
ResWP{1,1} = 'Rectangle window';
ResWP{2,1} = 'Hann';


%Mean and var estimation errors
MD = MatDistWP-DistTrue;
for ind=1:2
    ResWP{ind,2} = mean(MD(ind,:));
    ResWP{ind,3} = var(MD(ind,:),1);
end

% Save all useful results in the 'resultsWP.mat'
fname2 = strcat(fname,'resultsWP.mat');
save(fname2, 'CEPaWP', 'CEPbWP', 'MatDistWP',"ResWP");


end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the estimated copstral coefficients with the windowed peridogram
function CEP = comp_CEP_WP(Y)

[N, Nr] = size(Y);
CEP     = Inf(2,N,Nr);

for run=1:Nr
    [CEP(1,:,run), ~] = wp( Y(:,run), 'rect' ); % rectangular window
    [CEP(2,:,run), ~] = wp( Y(:,run), 'hann' ); % Hann window
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
        fprintf('error \n ')
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

% Check the estimated cepstral coefficients
c_p = ifft(log(Phi_p));
 if max(abs(imag(c_p)))<10^(-4)
    c_p = real(c_p);
 else
    fprintf('ERROR: in the evaluation of c_p %14.13f\n', max(abs(imag(c_p))));
    c_p = [];
    return
 end

end %function




