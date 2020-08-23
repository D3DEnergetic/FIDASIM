title: Visualization

[TOC]

---

FIDASIM comes with visualization executables that make plotting inputs and outputs easy.
The code can be executed in the command line prompt or imported as a Python package.
This tutorial will focus on the former and provide many examples.
The scripts are configured to detect and use the default Python installation.
You can check the default path by running `which python` in the command line prompt.
If you do not have Python installed, please visit [here](https://www.python.org/downloads/) for more details.

You can skip the next part if the visualization scripts and data are located on your personal device.
Otherwise we have a few suggestions for users computing on a cluster (e.g. Iris or Portal).

*[NoMachine]( https://www.nomachine.com/download) is a useful application that will allow you to plot figures from within the cluster.
*If you insist on not using NoMachine, then you will need to set up X11 forwarding.
For example, Mac users can download [this application]( https://www.xquartz.org/).

The scripts depend upon the following packages

* numpy
* h5py
* matplotlib
* os
* mpl_toolkits
* f90nml
* re
* argparse

If you are plotting locally, use `conda` or `pip` to install dependencies.
If you are on a cluster, contact your Computer Systems Administrator to get the requirements downloaded.

#Visualizing Inputs

The inputs section assumes the following

* Example directory: `/u/garciaav/test`
* Run IDs: `330LT` and `LHD`.

Plot everything (plasma, fields, distribution function and geometry):
```bash
plot_inputs /u/garciaav/test/ 330LT
```

Plot the plasma
```bash
plot_inputs /u/garciaav/test/ 330LT -p
```
![](|media|/visualize1.png){: width="200"}
{: style="text-align: center"}

Plot the fields
```bash
plot_inputs /u/garciaav/test/ 330LT -f
```
![](|media|/visualize2.png){: width="200"}
{: style="text-align: center"}

Plot the distribution function
```bash
plot_inputs /u/garciaav/test/ 330LT -d
```
![](|media|/visualize3.png){: width="200"}
{: style="text-align: center"}

Plot the diagnostic and beam geometry
```bash
plot_inputs /u/garciaav/test/ 330LT -g
```
![](|media|/visualize4.png){: width="200"}
{: style="text-align: center"}

Hide the legend
```bash
plot_inputs /u/garciaav/test/ 330LT -g -l
```
![](|media|/visualize5.png){: width="200"}
{: style="text-align: center"}

Notes:

* Black: NBI centerline and the half-widths
* Green: Beam grid boundaries
* Autumn: Diagnostic line-of-sights
* Purple: Plasma boundary approximated as a torus

Plot lineouts at R = 153.4 cm
```bash
plot_inputs /u/garciaav/test/ 330LT -p -f -d -rz 153.4
```
![](|media|/visualize6-1.png "Te rzlineout"){: width="400"}
{: style="text-align: center"}

![](|media|/visualize6-2.png "Br rzlineout"){: width="400"}
{: style="text-align: center"}

![](|media|/visualize6-3.png "Fzr rzlineout"){: width="400"}
{: style="text-align: center"}

Plot lineouts at Z = -10.8 cm:
```bash
plot_inputs /u/garciaav/test/ 330LT -p -f -d -zr -10.8
```
![](|media|/visualize7-1.png "Te zrlineout"){: width="400"}
{: style="text-align: center"}

![](|media|/visualize7-2.png "Br zrlineout"){: width="400"}
{: style="text-align: center"}

![](|media|/visualize7-3.png "Fzr zrlineout"){: width="400"}
{: style="text-align: center"}

Plot Energy and Pitch lineouts for the distribution function
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

For 3D inputs, plot a slice by specifying the toroidal angle in radians
```bash
plot_inputs /u/garciaav/test/ LHD -p -f -d -ph 0.157
```
![](|media|/visualize9-1.png "Te"){: width="200"} ![](|media|/visualize9-2.png "Br"){: width="200"} ![](|media|/visualize9-3.png "Fzr"){: width="200"}
{: style="text-align: center"}

Lastly, below are a few more useful optional arguments for plotting inputs

* `-nl` number of contour levels
* `-si /u/garciaav/test/figs/` save important plots to a folder
* `-h` display help message

You are all set to explore FIDASIM inputs!

#Visualizing Outputs

The outputs section assumes the following

* Example directories: `/u/garciaav/test/` and `/u/garciaav/different/`
* Run IDs `BAAE_1575`, `BAAE_1875`, `test_1`, `test_2`, and `diff_1575`

Plotting spectral and NPA data can be done in two ways.
If all of your data is located in one folder, then scroll down to Case 1.
However, if your data is located in more than one folder or you want to plot individually selected items from a single folder, then follow the conventions of Case 2.

##Case 1: Plotting from a single directory

The simplest command instructs the code to look inside a single directory.

* `-d` path of single folder

If the code finds `*_spectra.h5` files in the folder, you must indicate emission and channel information.

* `-s` plot all emission,
* `-fi` FIDA
* `-pf` passive FIDA
* `-f` Full BES
* `-hf` Half BES
* `-t` Third BES
* `-b` Bremsstrahlung
* `-c` Cold D-alpha
* `-hl` Halo
* `-dc` DCX
* `-as` all spectral channels
* `-rs` range of channels
* `-ls` list of channels
Note: `-s` cannot be used with any other emission type argument.

Similarly, if the code finds any `*_npa.h5` files in the folder, you must indicate flux and channel information.
The syntax of the NPA arguments is similar to the spectral arguments.

* `-n` plot all NPA flux types
* `-np` NPA
* `-pn` passive NPA
* `-an` all NPA channels
* `-rn` range of channels
* `-ln` list of channels

If you have many output files in the folder, below are two filtering techniques

* `-os` Force only spectral plots
* `-on` Force only NPA plots
* `-r` Force a specific runid

Below are some examples for Case 1 capabilities

`plot_outputs -d /u/garciaav/test/ -s -as -n -an`
Plots all spectral data, all spectral channels, all NPA data and all NPA channels.

![](|media|/visualize10-1.png "Spectra 1-9"){: width="200"} ![](|media|/visualize10-2.png "Spectra 10-18"){: width="200"} ![](|media|/visualize10-3.png "Spectra 19-24"){: width="200"}
{: style="text-align: center"}

![](|media|/visualize10-4.png "NPA 1-3"){: width="200"}
{: style="text-align: center"}

Note: The first Python index is 0, but the code will add 1 to all subplots in order to avoid displaying “Channel 0”.

`plot_outputs -d /u/garciaav/test/ -s -rs 11 19 -n -rn 1 2`
Plots all spectral data, range of spectral channels (11-19), all NPA data and range of NPA channels (1-2).

![](|media|/visualize11-1.png "Spectra 11-19"){: width="200"}
{: style="text-align: center"}

![](|media|/visualize11-2.png "NPA 1-2"){: width="200"}
{: style="text-align: center"}

`plot_outputs -d /u/garciaav/test/ -fi -ls 1 14 22 -np -ln 2`
Plots active FIDA data, list of spectral channels (1, 14 & 19), active NPA data and list of NPA channels (2).

![](|media|/visualize12-1.png "Spectra 1, 14 & 22"){: width="200"}
{: style="text-align: center"}

![](|media|/visualize12-2.png "NPA 2"){: width="200"}
{: style="text-align: center"}

`plot_outputs -d /u/garciaav/test/ -r BAAE_1875 -fi -ls 12 14 22 -os`
Plots run ID BAAE_1875, active FIDA data, list of spectral channels (12, 14 & 22), and ignores NPA files.

![](|media|/visualize13.png "FIDA runid 12, 14 & 22"){: width="200"}
{: style="text-align: center"}

##Case 2: Plotting from anywhere

In this case, the single folder constraint is removed and the user must provide full file paths

* `-p` paths of selected files

Below are some examples for Case 2 capabilities

`plot_outputs -p /u/garciaav/different/diff_1575_spectra.h5 -s -rs 10 18 -n -rn 2 3`
Plots all spectral data, range of spectral channels (10-18), all NPA data and range of NPA channels (2-3).

![](|media|/visualize14-1.png "Spectra 10-18"){: width="200"}
{: style="text-align: center"}

![](|media|/visualize14-2.png "NPA 2-3"){: width="200"}
{: style="text-align: center"}

`plot_outputs -p /u/garciaav/different/diff_1575_spectra.h5 -fi -ls 1`
Plots active FIDA data and list of spectral channels (1).

![](|media|/visualize15.png "FIDA 1"){: width="200"}
{: style="text-align: center"}

##Changing subplot parameters

Below are the arguments to customize the subplots for spectral plots

* `-sx` change spectral x limits
* `-sy` change spectral y limits
* `-sl` spectral log plot

and similarly NPA plots
* `-nx` change NPA x limits
* `-ny` change NPA y limits
* `-nl` NPA log plot

Below are examples of subplot modifications

`plot_outputs -p /u/garciaav/different/diff_1575_spectra.h5 -fi -ls 1 -sx 649 663`
Plots active FIDA data, list of spectral channels (1) and sets wavelength limits (649-663).

![](|media|/visualize16.png "FIDA xlim"){: width="200"}
{: style="text-align: center"}

`plot_outputs -p /u/garciaav/different/diff_1575_spectra.h5 -fi -ls 1 -sx 649 663 -sy 3e11 3e16 -sl`
Plots active FIDA data, list of spectral channels (1), sets wavelength limits (649-663), sets radiance limits (\(3\times 10^11\),\(3\times 10^16\) ) and turns on the log scale.

![](|media|/visualize17.png "FIDA xlim ylim log"){: width="200"}
{: style="text-align: center"}

`plot_outputs -p /u/garciaav/test/test_1_npa.h5 -n -ln 1 -ny 4e8 2e11 -nl`
Plots NPA data, list of NPA channels (1), sets flux limits (\(4\times 10^8\),\(2\times 10^11\) ) and turns on the log scale.

![](|media|/visualize18.png "NPA ylim log"){: width="200"}
{: style="text-align: center"}

Lastly, below are a few more useful optional arguments for plotting outputs

* `-p` can be used with more than one argument (e.g. -p filepath1 filepath2 ...)
* `-ss -sn /u/garciaav/test/figs/` save spectral and NPA plots to a folder
* `-h` display help message

Congratulations, you made it to the end of the tutorial!

If you have any questions or find a bug, please let us know on [GitHub](https://github.com/D3DEnergetic/FIDASIM/issues)
