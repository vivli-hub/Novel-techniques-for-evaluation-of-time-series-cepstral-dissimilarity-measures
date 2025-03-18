# Experiment
Evaluate the cepstral distance between simulated time series (sine-waves in additive noise)

# Authors
Miaotian Li, Ciprian Doru Giurcaneanu

## Data
1. Simulated Data: The results reported for simulated data are obtained by using the functions Test.m and Testcluster.m. More detailed information can be found in the Code section.
2. ECG signals: The ECG dataset utilized in this study was sourced from the EDG dataset available at PhySioNet [DataBase](https://physionet.org/content/ecg-fragment-high-risk-label/1.0.0/). This database comprises a collection of 2-second ECG signal segments with 360HZ sample rating exhibiting rhythm disturbances, categorized into distinct classes based on the severity of the threat to the patient’s life. For our experiments, we focused on two specific groups from this dataset: the first group consisted of 169 time series representing 2-second ECG recordings from individuals experiencing high rate ventricular tachycardia (VTHR) (The data is stored under the 3_Threatening_VT folder in the PhySioNet Database.), while the second group included 107 time series representing 2-second ECG recordings from healthy individuals with normal sinus rhythm (N) (The data is stored under the 6_Sinus_rhythm folder in the PhySioNet Database.). 

## Code
- [Test.m](#Test)
- [TestCluster.m](#TestCluster)
- [readdat.m](#readdat)
- [exp_ecg.m](#exp_ecg)

### Test
This is the main function used to generate the simulated data and calculate the mean and variance of the estimation errors for the cepstral distance. 
The function can be called as Test(N, vec, snr).

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%N`**    | The length of each time series                                                              |
| **`%vec`** | The matrix in the formula of the distance:<br> 'Martin' = Martin matrix<br> 'Identity' = Identity matrix |
| **`%snr`**   | Signal-to-Noise Ratio (in dB)   
| **`%dis`**   | Options Distance:<br> 'sq' = squared distance<br> 'eu' = distance |

**Output**
Test produces five different files. All of them are stored in one folder, whose name is given by `strcat('./N',num2str(N),'SNR',num2str(snr), '_', vec, '_', dis, '/')`. All of files have the extension .mat. 

`simdata.mat`

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`Ymata`**    | N×5000 matrix, consisting of 5000 simulated time series of length `%N`, generated based on the model of Time Series 1|                                                              |
| **`Ymatb`**    | N×5000 matrix, consisting of 5000 simulated time series of length `%N`, generated based on the model of Time Series 2|    
| **`DistTrue`** | The true distance between the Time Series 1 and Time Series 2   |

`resultsWP.mat`

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`CEPaWP`**    | (2*N)×5000 matrix, , which contains the cepstral coefficients computed from 'Ymata' in 'simdata.mat'. It includes results obtained by applying the rectangle window and Hann window.|                                                              |
| **`CEPbWP`**    | Similar to `CEPaWP`, the only difference is that it is based on the computation results of `Ymatb`|    
| **`MatDistWP`** | 2x5000 matrix, the estimated distance between the Time Series 1 and Time Series 2 based on cepstral coefficients (with the rectangle window and Hann window) |
| **`ResWP`** | 2x3 cell, the three columns represent the method name, the mean and the variance for the estimation errors |

 The contents of `resultsCN.mat`, `0.05resultsFDR.mat` and `0.01resultsFDR.mat` are essentially the same as in `resultsWP.mat`, with the only difference being the use of different methods to obtain the cepstral coefficients. The cepstral coefficients in `resultsCN.mat` were obtained using cepstral nulling with thresholds: BIC, KSF, and MRI. `0.05resultsFDR.mat` and `0.01resultsFDR.mat` are the cepstral coefficients obtained by applying FDR and FER thresholds at pre-specified FDR or FER values of 0.01 and 0.05, respectively.

### TestCluster

This is the main function used to generate the simulated data and cluster these time series based on the estimated cepstral distance. It calculates the similarity index to evaluate the performance of each method. For example, the function can be called as TestClust(Nruns, N, NoTS, snr0, wmat, clust_dist).

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%Nruns`**    | Number of runs in the experiment                                                              |
| **`%N`**    | The length of each time series                                                              |
| **`%NoTS`** | Number of Time Series in each cluster |
| **`%snr0`**   | Signal-to-Noise Ratio in dB                                                        |
| **`%dist`**   | The weighted matrix:<br> 'Martin' = Martin weighted matrix<br> 'Identity' = Identity weighted matrix |
| **`%clust_dist`**   | The distance used in the cluster:<br> 'euclidean' = distance <br> 'sqEuclidean' = squared distance |

**Output**
Test produces three different files. All of them are stored in one folder, whose name is given by `strcat('./N',num2str(N),'SNR', num2str(snr0), 'noTs', num2str(NoTS), '_', wmat, '_', clust_dist, '/')`. 

`similarity.mat`

The `%Nruns`×9 matrix records the clustering similarity. The 9 columns correspond to the following methods for computing cepstral coefficients: Rectangle window, Hann window, cepstral nulling with BIC, KSF, MRI, FDR with a pre-specified level of 0.01, FER with a pre-specified level of 0.01, FDR with a pre-specified level of 0.05, and FER with a pre-specified level of 0.05.

### readdat

Data can be downloaded from: [DataBase](https://physionet.org/content/ecg-fragment-high-risk-label/1.0.0/). In this experiment, we only need the data from the 3_Threatening_VT and 6_Sinus_rhythm folders. Function `readdat.m` reads the data from these two folders and creates a Mat file for each ECG signal in the folders N and VTHR.

### exp_ecg

Once you have these Mat files, you can run the main function `exp_ecg.m`. Additionally, the function requires `index.mat` to run. The index.mat file contains N and VTHR, each of which is a 100×100 matrix. These matrices store the indices of the 100 samples to be selected in 100 sampling iterations. In this function, process the ECG signals as follows: (i) Perform second-order differencing, (ii) Compute various types of estimated cepstral coefficients for the altered time series (including Hann window, Rectangle window, BIC threshold, KSF threshold, MRI threshold, and FDR and FER thresholds at pre-specified FDR or FER values of 0.01 and 0.05), (iii) Perform clustering using K-medoids with Distance, K-means with squared Distance, and K-medoids with squared Distance. The Distance is evaluated with the Identity matrix and the Martin matrix, (iv) Compute the similarity index to compare the clustering results with the ground truth.

**Input**

**`times`** represents the number of runs.

**Output**

**`results_ecg.mat`**: 

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`sim_vector_Identity_dis`**    | `times`x9 matrix (similarity index for Distance computed with Identity matrix, K-medoids clustering)|                                                             
| **`sim_vector_Identity_sqdis`**    | `times`x9 matrix (similarity index for squared Distance computed with Identity matrix, K-medoids clustering)|    
| **`sim_vector_Martin_dis`** | `times`x9 matrix (similarity index for Distance computed with Martin matrix, K-medoids clustering) |
| **`sim_vector_Martin_sqdis`** | `times`x9 matrix (similarity index for squared Distance computed with Martin matrix, K-medoids clustering)|
| **`sim_vector_Identity_kmeans`** | `times`x9 matrix (similarity index for squared Distance computed with Identity matrix, K-means clustering) |
| **`sim_vector_Martin_kmeans`** | `times`x9 matrix (similarity index for squared Distance computed with Martin matrix, K-means clustering)|

## Data management

- [make_table_sd.m](#make_table_sd)
- [plot_sim.m](#plot_sim)
- [make_table_ecg.m](#make_table_ecg)
- [dtw.Rmd](#dtw)
- [generate_graph.Rmd](#generate_graph)

### make_table_sd

Organize the output generated by the `Test.m` into an easy-to-read Excel file. For example, the function can be called as make_table_sd(specif). You can input specific characters to organize the files in all folders containing that character, such as "SNR". The output file has the extension .csv and is named 'mean_sd_eu'. The contents of the .mat files generated by Test.m are too large to upload. Therefore, only the CSV file processed by this function are included in the GitHub repository.

### plot_sim

The `plot_sim.m` function evaluates clustering performance across various signal lengths and methods for computing cepstral coefficients. It calculates similarity index between clustering results and ground truth, then plots these metrics and saves them to an Excel file. The contents of the .mat files generated by TestCluster.m are too large to upload. Therefore, only the xlsx file produced by this function are included in the GitHub repository.

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%NoTS`** | Number of Time Series in each cluster |
| **`%snr0`**   | Signal-to-Noise Ratio in dB                                                        |
| **`%dist`**   | The weighted matrix:<br> 'Martin' = Martin matrix<br> 'Identity' = Identity matrix |
| **`%clust_dist`**   | The distance used in the cluster:<br> 'euclidean' = distance <br> 'sqEuclidean' = squared distance |

**Output**
`SNR-10_sim_Identity.xlsx`
`SNR10_sim_Identity.xlsx`
`SNR-10_sim_Martin.xlsx`
`SNR10_sim_Martin.xlsx`

### make_table_ecg

Convert the data from `results_ecg.mat` to xlsx format, and add the type of cepstral coefficients, the type of distance, and the type of cluster corresponding to each result.

**Input**

**`results_ecg.mat`**

**Output**

**`ecg_table_100.xlsx`**: The results in 100 experiments

### dtw

We used the package `dtwclust` (https://cran.r-project.org/web/packages/dtwclust/index.html) to compute the Dynamic Time Warping distance between ECG signals. K-medoids was applied for clustering, and the similarity index was calculated.

**Input**

**`VTHR.csv`**: Convert the Mat files in the N folder to CSV format

**`N.csv`**: Convert the Mat files in the VTHR folder to CSV format

**`VTHR_index.csv`**: Convert VTHR from `index.mat` to CSV format

**`N_index.csv`**: Convert N from `index.mat` to CSV format

**Output**

**`true_label_1_100.xlsx`**: The true groupings of all ECG signals across 100 experiments

**`class_label_1_100.xlsx`**: The clustering results from 100 experiments


### generate_graph

1. Represent the results of the experiments done with `Test.m`. Use ggplot to aggregate and summarize the mean and standard deviation of the estimation errors

2. Represent the results of the experiments done with `TestCluster.m`. Use ggplot to aggregate and summarize the similarity index of the clustering based on the estimated cepstral distance

3. Represent the results of the experiments done with `exp_ecg.m`. Use ggplot to aggregate and summarize the mean and standard deviation of the similarity index

   

| **Type**   | **Files**                                                                                             |
|------------|-------------------------------------------------------------------------------------------------------|
| **Input**  | `simulation_mean_sd_SNR-1_.xlsx`<br>`simulation_mean_sd_SNR100_.xlsx`<br>`simulation_mean_sd_SNR-100_.xlsx`<br>`mean_sd_eu.csv`<br> `ecg_table_100.xlsx`|
| **Output** | `figrue2.pdf`<br>`mean_sd_SNR100.pdf`<br>`figure3.pdf`<br>`figure4.pdf` |

