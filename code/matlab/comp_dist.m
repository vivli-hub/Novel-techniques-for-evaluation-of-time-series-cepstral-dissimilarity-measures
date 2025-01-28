%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = comp_dist(cepa,cepb,vecw,dis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cepa             % cepstral coefficients for time series a
% cepb             % cepstral coefficients for time series b
% vecw             % weight martix
% dis              % sq-squared Eulcidean distance
%                  % eu-squared Eulcidean distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(cepa);
if length(cepb)~=N
    fprintf('Different no. of cepstral coefficients for the two time series.\n');
    dist = [];
    return;
end
cepa = cepa(:);
cepb = cepb(:);

diff = abs(cepa-cepb);

if strcmp(dis, 'sq')
    dist = vecw*(diff(1:N/2+1).^2);
else
    dist = sqrt(vecw*(diff(1:N/2+1).^2));
end

end %function

