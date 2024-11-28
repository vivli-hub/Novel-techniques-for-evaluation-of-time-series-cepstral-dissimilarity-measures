%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the similarity metric and silhouette width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%   CEP_WP   = cepstral ceofficients from order of 1,...,N with windowed
%   periodogram
%   CEP_nulling = cepstral ceofficients from order of 1,...,N applied
%   cepstral nulling
%   true label = the 'true' cluster for each of recordings
%   weight = weight matrix
%            1 : Identity matrix
%            2 : diagnal matrix with elements from 1, ..., N/2
%Output:
%   silhouette = 10*1 cell 
%   sim =  10*1 cell 
%   cluster = Nr*10 matrix
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [silhouette, sim, cluster] = sim_silhouette(CEP_WP, CEP_nulling, true_label, weight, No)
N = size(CEP_WP, 3);
%Windowed periodogram
silhouette_WP = zeros(2,1);
sim_WP = zeros(2,1);
cluster_WP = zeros(N,2);
for i = 1:2
    [cluster_no, silhouetteValues, simValues] = cluster_kmedoids(squeeze(CEP_WP(i,:,:)), 5, weight, No, true_label);
    similarity = mean(simValues);
    silhouette_WP(i) = mean(silhouetteValues);
    sim_WP(i) = similarity;
end
%Cepstral nulling

silhouette_nulling = zeros(7,1);
sim_nulling = zeros(7,1);
cluster_nulling = zeros(N,7);
for i = 1:7
    [cluster_no, silhouetteValues, simValues] = cluster_kmedoids(squeeze(CEP_nulling(i+1,:,:)), 5, weight, No, true_label);
    similarity = mean(simValues);
    silhouette_nulling(i) = mean(silhouetteValues);
    sim_nulling(i) = similarity;
end
 cluster = [cluster_WP, cluster_nulling];
 silhouette = [silhouette_WP; silhouette_nulling];
 sim = [sim_WP; sim_nulling];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_no, silhouetteValues, simValues] = cluster_kmedoids(c_e, numRuns, weight,k, true_label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
% c_e         = cepstral coefficients     
% numRuns     = number of Runs for kmedoids
% weight      = weight matrix
%            1 : Identity matrix
%            2 : diagnal matrix with elements from 1, ..., N/2
% k           = number of clusters
%Output:
% cluster_no   = matrix distances
%             size of cluster_no = method * Nr
%             method = number of methods for cepstral nulling
%             Nr = sample size for each run
% silhouette_values = the maximum silhouette index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(c_e, 1);
Nr = size(c_e, 2);
cluster_no = zeros(Nr,1);
diagonal_elements = arrayfun(@sqrt, 1:N/2);
diagonal_matrix = diag(diagonal_elements);
identity_matrix = ones(1, N/2);
if weight == 1
    weight_matrix = identity_matrix;
elseif weight == 2
    weight_matrix = diagonal_matrix;
end
D = zeros(N/2, Nr);
for i = 1:Nr
    D(:,i) = weight_matrix*c_e(2:N/2+1,i);
end

    % Define the number of clusters k
    allIdx = zeros(size(D', 1), numRuns);
    silhouetteValues = zeros(numRuns, 1);
    simValues = zeros(numRuns, 1);
    allMedoids = cell(numRuns, 1);
    for i = 1:numRuns
        rng(i); 
        % Perform clustering based on midpoint partitioning
        [idx, medoids] = kmedoids(D', k, 'Distance','euclidean' ,'replicates',10);
        allIdx(:, i) = idx;
        allMedoids{i} = medoids;
        index_1 = find(idx(:) == 1); 
        index_2 = find(idx(:) == 2);
        index_3 = find(idx(:) == 3);
        index_4 = find(idx(:) == 4);
        index_5 = find(idx(:) == 5);
        index_6 = find(idx(:) == 6);
        C = {index_1, index_2, index_3, index_4, index_5, index_6};
        % Silhouette width
        silhouetteValues(i) = mean(silhouette(D', idx));
        simValues(i) = cluster_similarity(true_label, C);
    end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end 
