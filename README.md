# Experiment
Evaluate the cepstral distance between simulated time series (sine-waves in additive noise)

# Authors
Miaotian Li, Ciprian Doru Giurcaneanu

## Data
The results reported for simulated data are obtained by using the functions Test.m and Testcluster.m. More detailed information can be found in the Code section.

## Code
- [Test.m](#Test)
- [TestCluster.m](#TestCluster)


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
| **`Ymata`**    | A N×5000 matrix, consisting of 5000 simulated time series of length `%N`, generated based on the model of Time Series 1|                                                              |
| **`Ymatb`**    | A N×5000 matrix, consisting of 5000 simulated time series of length `%N`, generated based on the model of Time Series 2|    
| **`DistTrue`** | The true distance between the Time Series 1 and Time Series 2   |

`resultsWP.mat`

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`CEPaWP`**    | A 2*N×5000 matrix, , which contains the cepstral coefficients computed from 'Ymata' in 'simdata.mat'. It includes results with applying the rectangle window and Hann window.|                                                              |
| **`CEPbWP`**    | Similar to `CEPaWP`, the only difference is that it is based on the computation results of `Ymatb`|    
| **`MatDistWP`** | A 2*5000 matrix, the estimated distance between the Time Series 1 and Time Series 2 based on cepstral coefficients with the rectangle window and Hann window |
| **`ResWP`** | A 2*3 cell, the three columns represent the method name, the mean of the bias for the estimated distance, and the variance of the estimated distance.|

 The contents among `resultsCN.mat`, `0.05resultsFDR.mat` and `0.01resultsFDR.mat` are essentially the same as in `resultsWP.mat`, with the only difference being the use of different methods to obtain the cepstral coefficients. The cepstral coefficients in `resultsCN.mat` were obtained using cepstral nulling with thresholds: BIC, KSF, and MRI. `0.05resultsFDR.mat` and `0.01resultsFDR.mat` are the cepstral coefficients obtained by applying FDR and FER thresholds at pre-specified FDR or FER values of 0.01 and 0.05, respectively.

### TestCluster

This is the main function used to generate the simulated data and cluster these time series based on the estimated cepstral distance. It represents the similarity index to evaluate the performance of each method. For example, the function can be called as TestClust(N, NoTS, snr0, dist, flag_plot).

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%N`**    | The length of each time series                                                              |
| **`%NoTS`** | Number of Time Series in each cluster |
| **`%snr0`**   | Signal-to-Noise Ration in dB                                                        |
| **`%dist`**   | The weighted matrix:<br> 'Martin' = Martin weighted matrix<br> 'Identity' = Identity weighted matrix |
| **`%flag_plot`**   | 0 = No plots<br> 1 = Plots for each method and each distance |

**Output**
Test produces three different files. All of them are stored in one folder, whose name is given by `strcat('./N',num2str(N),'SNR',num2str(snr0),'/')`. All of files have the extension .mat. 

`simdata.mat`

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`Ymata`**    | A N×NoTS matrix, consisting of 100 simulated time series of length `%N`, generated based on the model of time series in Cluster 1|                                                              |
| **`Ymatb`**    | A N×NoTS matrix, consisting of 100 simulated time series of length `%N`, generated based on the model of time series in Cluster 2|    

`results.mat`

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`silhouette_Martin`**    | A 9×1 matrix, clustering based on Martin Distance, resulting in the Silhouette coefficient. It includes a total of 9 different cepstral coefficients: Rectangle window, Hann window, BIC threshold, KSF threshold, MRI threshold, and FDR and FER thresholds at pre-specified FDR or FER values of 0.01 and 0.05.|                                                              |
| **`silhouette_ID`**    | A 9×1 matrix matrix, similar to `silhouette_Martin`, but it is based on clustering using Identity Distance.|    
| **`sim_Martin`**    | A 9×1 matrix matrix, clustering based on Martin Distance, resulting in the Similarity index.|    
| **`sim_ID`**    | A 9×1 matrix matrix, similar to `sim_Martin`, but it is based on clustering using Identity Distance.|  

`CEP.mat`
| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`CEP_WP`**    | A 2*N×2NoTS matrix, , which contains the cepstral coefficients computed from 'Ymata' and 'Ymatb' in 'simdata.mat'. It includes results with applying the rectangle window and Hann window.|                                                              |
| **`CEP_nulling`**    | A 8*N×2NoTS matrix, , which contains the cepstral coefficients computed from 'Ymata' and 'Ymatb' in 'simdata.mat'. It includes results with applying the rectangle window and cepstral nulling of thresholds: BIC, kSF, MRI, FDR and FER at pre-specified FDR or FER values of 0.01 and 0.05.|    

## Data management

- [make_table_sd.m](#make_table_sd)
- [make_table_cluster.m](#make_table_cluster)
- [plot_mean_sd.Rmd](#plot_mean_sd)
- [plot_cluster_similarity_index.Rmd](#plot_cluster_similarity_index)

### make_table_sd

Organize the output generated by the `Test.m` into an easy-to-read Excel file. For example, the function can be called as make_table_sd(specif). You can input specific characters to organize the files in all folders containing that character, such as "SNR". The output file has the extension .csv and is named 'mean_sd_eu'. The contents of the .mat files generated by Test.m are too large to upload. Therefore, only the CSV file processed by this function are included in the GitHub repository.

### make_table_cluster

Organize the output generated by the `TestCluster.m` into an easy-to-read Excel file. For example, the function can be called as make_table_cluster(specif). You can input specific characters to organize the files in all folders containing that character, such as "SNR-1". The output file has the extension .xlsx and is named 'similarity_index'.

### plot_mean_sd

Represent the results of the experiments. Use ggplot to aggregate and summarize the mean and standard deviation of the bias for the estimated cepstral distance.

| **Type**   | **Files**                                                                                             |
|------------|-------------------------------------------------------------------------------------------------------|
| **Input**  | `simulation_mean_sd_SNR-1_.xlsx`<br>`simulation_mean_sd_SNR100_.xlsx`<br>`simulation_mean_sd_SNR-100_.xlsx` |
| **Output** | `figure_2.pdf`<br>`figure_2'.pdf`<br>`mean_sd_snr100.pdf` |

### plot_cluster_similarity_index
Represent the results of the experiments. Use ggplot to aggregate and summarize the similarity index of the clustering based on the estimated cepstral distance.

| **Type**   | **Files**                                                                                             |
|------------|-------------------------------------------------------------------------------------------------------|
| **Input**  | `similarity_index.xlsx` |
| **Output** | `figure_3.pdf`|
