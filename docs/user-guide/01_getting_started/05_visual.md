title: Visualization

[TOC]

---

These visualization tools can be used straight from the command line window.
They require an installation of Python 2.7 or higher.
The visualization scripts are located in `FIDASIM/lib/scripts/`.

#Visualizing Inputs

The following Python packages are required for plotting the inputs.
If you are computing on a cluster, this is probably already done for you.

* numpy
* h5py
* matplotlib
* os
* mpl_toolkits
* f90nml

Below is a brief step-by-step tutorial on how to visualize FIDASIM inputs inside an example directory `/u/garciaav/test` with run ID `330LT` and `LHD`.

The simplest command will plot of your inputs (plasma, fields, distribution function and geometry):
```bash
plot_inputs /u/garciaav/test/ 330LT
```

If you are interested in being more specific, then use optional arguments `-p`, `-f`, `-d` or `-g`.
Plotting the plasma from the equilibrium file:
```bash
plot_inputs /u/garciaav/test/ 330LT -p
```
![](|media|/visualize1.png){: width="200"}
{: style="text-align: center"}

and the fields:
```bash
plot_inputs /u/garciaav/test/ 330LT -f
```
![](|media|/visualize2.png){: width="200"}
{: style="text-align: center"}
Note: Similar to the FIDASIM code, plots of the plasma and fields use the mask variable to define the edge of the plasma.

Similarly, 2D projections of the distribution function can be plotted by executing:
```bash
plot_inputs /u/garciaav/test/ 330LT -d
```
![](|media|/visualize3.png){: width="200"}
{: style="text-align: center"}
Note: The R-Z plots of the distribution function will also use the mask variable from the equilibrium file if it is available.

The geometry file can be visualized interactively in a 3D rotatable plot:
```bash
plot_inputs /u/garciaav/test/ 330LT -g
```
![](|media|/visualize4.png){: width="200"}
{: style="text-align: center"}

If you are using many channels, then your legend is probably crowding the plot.
The following command hides the legend:
```bash
plot_inputs /u/garciaav/test/ 330LT -g -l
```
![](|media|/visualize5.png){: width="200"}
{: style="text-align: center"}
Several things to note:
* NBI centerline and the half-widths are plotted in black
* Beam grid boundaries are drawn in green
* Diagnostic line-of-sights are in autumn colors
* Plasma boundary is approximated with a purple torus

Lineouts for the 2D plots of the plasma, fields and distribution function are available.
For example, use the following command to plot lineouts at R = 153.4 cm:
```bash
plot_inputs /u/garciaav/test/ 330LT -p -f -d -rz 153.4
```
![](|media|/visualize6-1.png "Te rzlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize6-2.png "Br rzlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize6-3.png "Fzr rzlineout"){: width="400"}
{: style="text-align: center"}

and conversely for lineouts at constant Z = -10.8 cm:
```bash
plot_inputs /u/garciaav/test/ 330LT -p -f -d -zr -10.8
```
![](|media|/visualize7-1.png "Te zrlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize7-2.png "Br zrlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize7-3.png "Fzr zrlineout"){: width="400"}
{: style="text-align: center"}

Since the distribution function is defined with 2 more coordinates (Energy and Pitch), more lineouts are available.
See the example below:
```bash
plot_inputs /u/garciaav/test/ 330LT -d -pz .2 -pr -.25 -ez 80 -er 40 -ep 18
```
![](|media|/visualize8-1.png "Fzp pzlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize8-2.png "Fpr prlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize8-3.png "Fze ezlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize8-4.png "Fer erlineout"){: width="400"}
{: style="text-align: center"}
![](|media|/visualize8-5.png "Fpe eplineout"){: width="400"}
{: style="text-align: center"}
Note: The first letter corresponds to the variable for the constant line and the second letter corresponds to the variable of the axis that is parallel to the constant line.

If you provided 3D inputs, then specify the toroidal angle in radians to view a slice of your inputs.
```bash
plot_inputs /u/garciaav/test/ LHD -p -f -d -ph 0.157
```
![](|media|/visualize9-1.png "Te"){: width="200"} ![](|media|/visualize9-2.png "Br"){: width="200"} ![](|media|/visualize9-3.png "Fzr"){: width="200"}
{: style="text-align: center"}

You can change the number of contour levels used in the 2D plots with `-nl`

Always remember that you can execute `plot_inputs -h` to display the help message.
You are all set to explore FIDASIM inputs!


#Visualizing Outputs

The following Python packages are required for plotting the outputs.

* os
* re
* h5py
* argparse
* numpy
* matplotlib

Below is a brief step-by-step tutorial on how to visualize FIDASIM outputs inside example directories `/u/garciaav/test/` with run IDs `BAAE_1575`, `BAAE_1875`, `test_1` and `test_2` and `/u/garciaav/different/` with run ID `diff_1575`

The code executes in one of the following two cases:

* Plotting all available `*_spectra.h5` and `*_npa.h5` files located in a single directory
* Plotting selected files located anywhere

##Case 1: Plotting data from a single directory

The simplest command is to plot everything located in a single folder with the `-d` argument.
The code will search the provided directory and attempt to plot all available spectral and NPA data.

If the code finds `*_spectra.h5` files in the `-d` folder, you must indicate emission and channel information.

Use `-s` to plot all emission,
or select individually with `-fi` (FIDA), `-pf` (pFIDA), `-f` (Full), `-hf` (Half), `-t` (Third), `-b` (Bremsstrahlung), `-c` (Cold), `-hl` (Halo), and `-dc` (DCX).

Use `-as`, `-rs` or `-ls` to indicate all spectral channels, range of channels or list of channels.
Note: The first Python index is 0, but the code will add 1 to all subplots in order to avoid displaying “Channel 0”.

Similarly, if the code finds `*_npa.h5` files in the `-d` folder, you must indicate flux and channel information.
The syntax of the NPA arguments is similar to the spectral arguments.

Use `-n` to plot all emission,
or select individually with `-np` (NPA) and `-pn` (pNPA).

Use `-an`, `-rn` or `-ln` to indicate all NPA channels, range of channels or list of channels.

If you have many FIDASIM outputs in the `-d` folder, below are two ways you can filter what the visualization script plots.

Filter out the NPA files and plot only spectral data by using the `-os` argument.
On the other hand, filter out the spectral files and plot only NPA data by using the `-on` argument.

Also, you can filter files by defining the run IDs with the `-r` argument.

###Below is a summary of examples for Case 1:

`plot_outputs -d /u/garciaav/test/`

`-s -as -n -an`
Plots all spectral data, all spectral channels, all NPA data and all NPA channels.
![](|media|/visualize10-1.png "Spectra 1-9"){: width="200"}
{: style="text-align: center"}
![](|media|/visualize10-2.png "Spectra 10-18"){: width="200"}
{: style="text-align: center"}
![](|media|/visualize10-3.png "Spectra 19-24"){: width="200"}
{: style="text-align: center"}
![](|media|/visualize10-4.png "NPA 1-3"){: width="200"}
{: style="text-align: center"}

`-s -rs 11 19 -n -rn 1 2`
Plots all spectral data, range of spectral channels (11-19), all NPA data and range of NPA channels (1-2).
![](|media|/visualize11-1.png "Spectra 11-19"){: width="200"}
{: style="text-align: center"}
![](|media|/visualize11-2.png "NPA 1-2"){: width="200"}
{: style="text-align: center"}

`-fi -ls 1 14 22 -np -ln 2`
Plots active FIDA data, list of spectral channels (1, 14 & 19), active NPA data and list of NPA channels (2).
![](|media|/visualize12-1.png "Spectra 1, 14 & 22"){: width="200"}
{: style="text-align: center"}
![](|media|/visualize12-2.png "NPA 2"){: width="200"}
{: style="text-align: center"}

`-r BAAE_1875 -fi -ls 12 14 22 -os`
Plots run ID BAAE_1875, active FIDA data, list of spectral channels (12, 14 & 22), and ignores NPA files.
![](|media|/visualize13.png "FIDA runid 12, 14 & 22"){: width="200"}
{: style="text-align: center"}

##Case 2: Plotting data from user defined file paths

Alternatively, the code accepts specific file paths from the `-p` argument.

Similar to Case 1, specifying emission type and flux type is required when plotting data from `*_spectra.h5` and `*_npa.h5 files`, respectively.
Of course, channel information is also a requirement.

###Below is a summary of examples for Case 2:

`plot_outputs -p /u/garciaav/different/diff_1575_spectra.h5`

`-s -rs 10 18 -n -rn 2 3`
Plots all spectral data, range of spectral channels (10-18), all NPA data and range of NPA channels (2-3).
![](|media|/visualize14-1.png "Spectra 10-18"){: width="200"}
{: style="text-align: center"}
![](|media|/visualize14-2.png "NPA 2-3"){: width="200"}
{: style="text-align: center"}

`-fi -ls 1`
Plots active FIDA data and list of spectral channels (1).
![](|media|/visualize15.png "FIDA 1"){: width="200"}
{: style="text-align: center"}

##Changing subplot parameters (applies to Case 1 & 2)

You can change the x and y limits of the spectral plot with `-sx` and `-sy` arguments.
Also, you can use `-sl` to plot the y axis on a log-scale.

Similarly, change the NPA plots with the `-nx`, `-ny` and `-nl` arguments.

###Below are some examples:

`plot_outputs -p /u/garciaav/different/diff_1575_spectra.h5`

`-fi -ls 1 -sx 649 663`
Plots active FIDA data, list of spectral channels (1) and sets wavelength limits (649-663).
![](|media|/visualize16.png "FIDA xlim"){: width="200"}
{: style="text-align: center"}

`-fi -ls 1 -sx 649 663 -sy 3e11 3e16 -sl`
Plots active FIDA data, list of spectral channels (1), sets wavelength limits (649-663), sets radiance limits (3$\times$10$^11$-3$\times$10$^16$) and turns on the log scale.
![](|media|/visualize17.png "FIDA xlim ylim log"){: width="200"}
{: style="text-align: center"}

`plot_outputs -p /u/garciaav/test/test_1_npa.h5`

`-n -ln 1 -ny 4e8 2e11 -nl`
Plots NPA data, list of NPA channels (1), sets flux limits (4$\times$10$^8$-2$\times$10$^11$) and turns on the log scale.
![](|media|/visualize18.png "NPA ylim log"){: width="200"}
{: style="text-align: center"}

Lastly, execute `plot_inputs -h` to display the help message.
Congrats, you made it to the end of the tutorial!

If you have any questions or find a bug, please let us know on [GitHub](https://github.com/D3DEnergetic/FIDASIM/issues)
