# Pearson
Modelising matter density field PDF with Pearson differential equation

The differential equation is


   $$\frac{1}{f}\frac{df}{dx} = Q(x) = \frac{x-a_0}{a_1+a_2x+a_3x^2}$$

$a_i$ are free parameters and can be linked to the non-centered statistical moments via the recurrent equation:

$$a_0m_n+na_1m_{n-1}+(n+1)a_2m_n+(n+2)a_3m_{n+1}=-m_{n+1}$$
