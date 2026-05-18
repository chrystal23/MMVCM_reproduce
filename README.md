# MMVCM_reproduce
Source code for reproducing the results in "Multi-dimensional multi-threshold varying coefficient model".

## Usage instructions

Three main folders correspond to the different components of the reproducibility materials:
- `MMVCM_functions`: Contains the core R functions and C++ implementations of our proposed MMVCM method.
- `Simulation`: Contains the R scripts and data needed to reproduce the results of multiple simulation examples.
- `RealData`: Contains the raw dataset, specific functions, and analysis R scripts to reproduce the real data application on infant sleep and maternal mental health.

A detailed list of the files within each folder is provided below.

### Folder: MMVCM_functions
|Code|Function|
|:-|:-|
|`simulation_function.R`|Core R functions used for generating the simulation data and implementing our MMVCM method.|
|`CPPlwls2d_s1.cpp`|C++ script implementing the local kernel smoothing for MMVCM Stage 1 (Estimation of the linear part).|
|`CPPlwls2d_s2.cpp`|C++ script implementing local kernel smoothing for estimating the smooth part in MMVCM Stage 2 (Estimation of the nonparametric part ).|
|`CPPlwls2d_s2_LR.cpp`|C++ script implementing one-sided local kernel smoothing for estimating the jumps in MMVCM Stage 2.|
|`compare.R`|R script used to implement the competing methods, including the traditional traditional semi-varying coefficient model (SVCoef), varying coefficient model (VCoef), and the threshold regression model with multiple change points (seqMS).|

### Folder: Simulation
|Code|Function|
|:-|:-|
|`exS_model.RData`|Saved R workspace containing the baseline model parameters and configuration for Example S.|
|`exM_model.RData`|Saved R workspace containing the baseline model parameters and configuration for Example M, M-Gauss, and M-Beta.|
|`Example_S.R`|Main R script to execute the simulation study for Example S.|
|`Example_M.R`|Main R script to execute the simulation study for Example M.|
|`Example_M-Gauss.R`|R script to execute the simulation study for Example M-Gauss, which follows the same model as Example M but generates index variables from multivariate Gaussian distribution.|
|`Example_M-Beta.R`|R script to execute the simulation study for Example M-Beta, which follows the same model as Example M but generates index variables from multivariate Beta distribution with Gaussian copula.|
|`Evaluation.R`|R script to evaluate, summarize, and visualize the final metrics from the simulation results.|

### Folder: RealData
|Code|Function|
|:-|:-|
|`Dataset_maternal_mental_health_infant_sleep.csv`|The raw dataset containing real-world maternal mental health and infant sleep records.|
|`Infant.R`|Main R script to run the application of the MMVCM method on the real dataset.|
|`infant_function.R`|R script containing helper functions specific to processing the infant sleep dataset and adapting our MMVCM method to the data characteristics.|
|`CPPlwls2d_s1_real.cpp`|C++ script implementing the local kernel smoothing for MMVCM Stage 1,tailored for the specific data structure of the infant sleep dataset.|
|`CPPlwls2d_s2_real.cpp`|C++ script implementing local kernel smoothing for estimating the smooth part in MMVCM Stage 2, tailored for the real dataset.|
|`CPPlwls2d_s2_LR_real.cpp`|C++ script implementing one-sided local kernel smoothing for estimating the jumps in MMVCM Stage 2, tailored for the real dataset.|


## Software MMVCM

We provide a R-based implementation of our MMVCM method for the users' convenience. For tutorials and other detials, please check [our repository](https://github.com/chrystal23) for the MMVCM R package.
