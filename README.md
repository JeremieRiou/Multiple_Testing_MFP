# Multiple Testing MFP
The purpose of this space is to provide an R function that will allow you to apply the theory developed in the paper by Liquet, B. Riou, J. and Roux, M. The aim is to explain it, illustrate on a concrete example, and to make the results presented in this paper reproducible.


## Arguments of Resampling Function

    - **formula** : A formula for time to event outcome, e.g. Surv(time,event)~varcod+adjustment_variable
    
    - data: A dataframe 
    
    - varcod: name of the column corresponding to the variable of interest to be transformed
    
    - FP: a matrix with powers which are used for computing each fractional polynomial transformations. The FP argument needs to be a matrix. The number of rows correspond to the number of transformations tested, and the number of columns is the maximum number of degrees tested for a single transformation.
        
        Example:
        
        [1,]	-2    NA	NA	    NA
        [2,]	0.5	  1	    -0.5	2
        [3,]	-0.5  1	    NA	    NA

        - In this example, the user performs three transformations of the variable of interest. The first is a fractional polynomial transformation with one degree and a power of -2. The second transformation is a fractional polynomial transformation with four degrees and powers of 0.5,1,-0.5,and 2. The third transformation is a fractional polynomial transformation with two degrees and powers of -0.5, and 1.
        
    - N: the number of resampling that you want to do.
    
    - txcensure: censoring rate you expected
    
    - method: The resampling method you require:  "Permutation" or "Bootstrap"
    
    - alpha: FWER 

# Motivating Example

See Motivating_Example.md for all details about the use of the function in a real context
