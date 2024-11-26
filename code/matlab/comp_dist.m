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


%dist = sqrt( vecw*(diff(1:N/2+1).^2) );

dist = vecw*(diff(1:N/2+1).^2);

end %function
