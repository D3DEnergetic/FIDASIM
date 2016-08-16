title: Collisional Radiative Model

[TOC]

---

# Collisional Radiative Model
The collisions that a neutral particle experiences as it travels through a plasma changes the distribution of its energy level population.
The types of collisions that FIDASIM takes into account is as follows

* Spontaneous decay: \(A_{m \rightarrow n} / A_{n \rightarrow m}\)
* Electron-impact excitation/de-excitation: \(q^e_{m \rightarrow n} / q^e_{n \rightarrow m}\)
* Ion-impact excitation/de-excitation: \(q^i_{m \rightarrow n} / q^i_{n \rightarrow m}\)
* Impurity-impact excitation/de-excitation: \(q^Z_{m \rightarrow n} / q^Z_{n \rightarrow m}\)
* Electron-impact ionization: \(I^e_n\)
* Ion-impact ionization: \(I^i_n\)
* Impurity-impact ionization: \(I^Z_n\)
* Charge exchange with ions: \(X^i_n\)
* Charge exchange with impurities: \(X^Z_n\)

The above rate coefficients have units of \(cm^3/s\) and are calculated by averaging the respective collisional cross sections with a Maxwellian of the relevent species.
The population of the \(n^{th}\) energy level of a neutral atom, \(f_n\), can be described by the following time dependent differential equation
$$
\frac{df_n}{dt} = - \left ( \sum_{k=i,Z} f_n d_k X^k_n + \sum_{k=e,i,Z} f_n d_k I^k_n \right ) \\
+ \sum_{m>n} \left (f_m A_{m \rightarrow n} + \sum_{k=e,i,Z} (f_m d_k q^k_{m \rightarrow n} - f_n d_k q^k_{n \rightarrow m}) \right ) \\
+ \sum_{n>m} \left (-f_n A_{n \rightarrow m} + \sum_{k=e,i,Z} (f_m d_k q^k_{m \rightarrow n} - f_n d_k q^k_{n \rightarrow m}) \right )
$$
where the \(d_k\) are the respective target densities.

By rearranging the terms and letting \(q^k_{n \rightarrow m}\) represent excitation/de-excitation depending on the order of the indices we can get the following equation
$$
\frac{df_n}{dt} = C_{n,n} f_n + \sum_{m \ne n} C_{n,m} f_m
$$
where
$$
C_{n,n} = - \left [ \sum_{k=i,Z} d_k X^k_n + \sum_{k=e,i,Z} d_k I^k_n + \sum_{m \ne n} \left ( A_{n \rightarrow m} + \sum_{k=e,i,Z} d_k q^k_{n \rightarrow m} \right ) \right ]
$$
and 
$$
C_{n,m} = A_{m \rightarrow n} + \sum_{k=e,i,Z} d_k q^k_{m \rightarrow n}
$$

The system of differential equations can be compactly represented as a matrix multiplication.
$$
\frac{d \vec{f}}{dt} = C \cdot \vec{f}
$$

The system of differential equations can be solved analytically to give
$$\vec{f} (t) = S^{-1} \cdot \vec{f} (0) \cdot S \cdot \exp(\Lambda \,t)$$
where \(\vec{f}(t)\) is a vector of the neutral population flux [1/s] for each energy state at time \(t\), \(S\) is the matrix of the eigenvectors of \(C\) and \(\Lambda\) is a diagonal matrix containing the eigenvales of \(C\). 
The fractional flux of a neutral traveling through a uniform plasma is shown below.

![](|media|/neutral_attenuation.png)
{: style="text-align: center"}

As you can see the relative populations between states converges fairly quickly.

The number of neutrals in a given state after a time \(t\), \(\vec{n}(t)\), is given by
$$\vec{n}(t) = S^{-1} \cdot \vec{n}(0)\cdot S \cdot (\exp(\Lambda \,t) - 1)/\Lambda$$
If \(t\) represents the time spent inside a grid cell the neutral density can be calculated by dividing the above equation by \(V_{cell}\).
The total neutral density of a mc marker is shown below.

![](|media|/neutral_dens.png)
{: style="text-align: center"}

As you can see over time the total number of neutrals decreases exponentially.

#Fortran References
* [[colrad]]: Fortran implementation
* [[get_rate_matrix]]: Constructs rate matrix \(A\)
* [[AtomicRates]]: Derived type that stores populating and de-populating transitions
