This repository contains two NumPy files, Dbins.npy and PDFs.npy, i.e. the probability density functions (PDFs) of the log-scaled cosmic overdensity field across different smoothing scales and redshifts.

Dbins.npy
    Stores the overdensity bin edges.
    Shape: (5, 6, 100)
    Dimensions:
        First index (5 values): Smoothing scales R=[2,4,8,16,32] Mpc /h
        Second index (6 values): Redshifts Z=[127,3,2,1,0.5,0]
        Third index (100 values): Bins for the log-scaled overdensity field

PDFs.npy
    Stores the probability density function values corresponding to the bins in Dbins.npy.
    Shape: (5, 6, 99)
    Dimensions:
        First index (5 values): Smoothing scales R=[2,4,8,16,32] Mpc
        Second index (6 values): Redshifts Z=[127,3,2,1,0.5,0]
        Third index (99 values): PDF values for the corresponding bins (one less than Dbins.npy due to bin edges)
