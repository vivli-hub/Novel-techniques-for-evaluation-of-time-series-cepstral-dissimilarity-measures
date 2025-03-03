function plot_sim(NoTS, snr0, wmat, clust_dist)
%NoTS = 100;
% snr0 = -1;
%snr0 = 10;
%wmat = 'Identity';
%clust_dist = 'euclidean';
% NoTS              =no. of TS in each cluster
% snr0              =SNR (in dB)
% wmat              =matrix: 'Martin' or 'Identity'
% clust_dist        =dist clustering: 'euclidean' or 'sqEuclidean'
mm = [];
for N=[128 256 512 1024 2048]

fname = strcat('/Volumes/Newsmy 1/phd/New_version/Simulation/Test_Cluster/N',num2str(N),'SNR', num2str(snr0), ...
    'noTs',num2str(NoTS), ...
    '_', wmat, ...
    '_', clust_dist, '/');
fname1 = strcat(fname,'similarity.mat');
load(fname1, 'sim_all_met');
mm = [mm; mean(sim_all_met)];

end

figure(1)
for j=1:size(mm,2)
    plot(mm(:,j),'-o');
    hold on
end
legend('Periodogram','Hann','BIC','KSF','MRI','FRD001','FER001','FRD005','FER005');
hold off


mm_col = mm(:);

B = repmat([128; 256; 512; 1024; 2048], 9, 1);

categories = ["Rectangle window"; "Hann"; "BIC"; "KSF"; "MRI"; ...
              "FDR0.01"; "FER0.01"; "FDR0.05"; "FER0.05"];
C = repelem(categories, 5);

final_table = table(mm_col, B, C, 'VariableNames', {'sim_metric', 'N', 'Method'});
name = strcat('SNR',num2str(snr0),'_sim_id.xlsx')
writetable(final_table, name);
end
