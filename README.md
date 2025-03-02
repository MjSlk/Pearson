# Pearson

Modelising matter overdensity probablity density function (PDF) with Pearson differential equation

The differential equation is


   $$\frac{1}{f}\frac{df}{dx} = Q(x) = \frac{x-a_0}{a_1+a_2x+a_3x^2}$$

$a_i$ are free parameters (Pearson coefficients) and can be linked to the non-centered statistical moments via the recurrent equation:

$$a_0m_n+na_1m_{n-1}+(n+1)a_2m_n+(n+2)a_3m_{n+1}=-m_{n+1}$$
By rewriting the previous equation for n = 0, 1, 2, 3 and by taking into account $m_0 = 1$ (sum of probabilities = 1), we get a system of 4 equations linking each of the $a_i$ to $m_1, m_2, m_3, m_4$.

This modeling approach is applied to the PDF of the logarithmic overdensity field defined as
$$\tilde{\delta} = \ln (1+\delta)$$, and its corresponding PDF, $\tilde {P} (\tilde {\delta})$ is linked to the the PDF of $\delta$, $P(\delta)$, via the following relation 
$\tilde {P} (\tilde {\delta}) = P(\delta) (1+\delta)$

Once the fit is applied to $\tilde {P} (\tilde {\delta})$ we compute the corresponding fit to $P(\delta)$.
To quantify the goodness of fit we compute the Kullback-Leibler divergency


The data on which we test this approach are the density field of the fiducial cosmological model computed from Quijote N-Body simulations (https://quijote-simulations.readthedocs.io/)

Dependencies

The project was developed with the following library versions. Running with other versions may crash or produce incorrect results.

    Python 3.8.10
    numpy==1.19.5
    matplotlib==4.0

