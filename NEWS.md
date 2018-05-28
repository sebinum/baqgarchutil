# baqgarchutil 0.3.0

The first release-version of`baqgarchutil`!

## Package structure

*   `baq_nifunction` - takes a `mGJR`-class object from `mgarchBEKK::mGJR` and 
    applies a news impact function function. It returns a `baq_nif`-class
    object containing the relevant parameters and the impacted covariances
    matrices
*   `baq_niplot` - takes a `baq_nif`-class object and creates a visual 
    interpretation of the results
*   `mv_ch_tests` - wrapper function for several tests to identify 
    conditional heteroscedasticity in multivariate time series.
*   `diag_mv_ch_model` is a wrapper function for several residual diagnostic 
    tests to determine model fit / adequacy
*   `diag_std_et` (and `diag_std_et_cnd`) - functions to transform a time 
    series (and it's conditional covariance matrices) to standardized series
    employed in `mv_ch_tests` and `diag_mv_ch_model`, also called in the
    following `diag_dufour_roy`-function
*   `diag_ljung_box` - the Ljung-Box Test Statistic for serial correlation in 
    uni- and multivariate time series (used in the above mentioned wrapper 
    functions)
*   `diag_dufour_roy` - the Rank-Based Test Statistic for serial correlation in 
    uni- and multivariate time series (used in the above mentioned wrapper 
    functions)
    
    



