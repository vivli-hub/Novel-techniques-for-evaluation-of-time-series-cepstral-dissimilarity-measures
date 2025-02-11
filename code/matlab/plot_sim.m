

NoTS = 100;
% snr0 = -1;
snr0 = -1;
wmat = 'Martin';
clust_dist = 'euclidean';

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

writetable(final_table, 'SNR-1_sim_martin.xlsx');