function windowed_per(fname, fname1, vecw)

load(fname1,'Ymata', 'Ymatb', 'DistTrue');

CEPaWP = comp_CEP_WP(Ymata);
CEPbWP = comp_CEP_WP(Ymatb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nr = size(Ymata,2);
MatDistWP = Inf*ones(2,Nr);
for met=1:2
    for rr=1:Nr
        MatDistWP(met,rr) = comp_dist(CEPaWP(met,:,rr),CEPbWP(met,:,rr),vecw);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResCepNulling = cell(2,3);
ResCepNulling{1,1} = 'Rectangle window';
%  ResCepNulling{2,1} = 'Hamming';
ResCepNulling{2,1} = 'Hann';
% ResCepNulling{4,1} = 'hann with half overlapping';


MD = MatDistWP-DistTrue;
for ind=1:2
    ResCepNulling{ind,2} = mean(MD(ind,:));
    ResCepNulling{ind,3} = var(MD(ind,:),1);
end

fname2 = strcat(fname,'resultsWP.mat');
save(fname2, 'CEPaWP', 'CEPbWP', 'MatDistWP',"ResCepNulling");


end %function

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
     
      %Phi_WP = periodogramse(y,ones(1, N),N);
       [Phi_WP, freq] = periodogram(y, ones(1, N), (0:2*pi/N:2*pi*(N-1)/N));
 elseif strcmp(v, 'hann')
     %Phi_WP = periodogramse(y,hann(N),N);
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

c_p = ifft(log(Phi_p));
 if max(abs(imag(c_p)))<10^(-4)
    c_p = real(c_p);
 else
    fprintf('ERROR: in the evaluation of c_p %14.13f\n', max(abs(imag(c_p))));
    c_p = [];
    return
 end

end %function


