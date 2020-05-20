title: Visualization

[TOC]

---

These visualization tools can be used straight from the command line window.
They require an installation of Python 2.7 or higher.
The visualization scripts are located in `FIDASIM/lib/scripts/`.

#Visualizing Inputs

Please install the following modules with `pip install` or `conda install`.

* numpy
* h5py
* matplotlib
* os
* mpl_toolkits
* f90nml

Below is a brief step-by-step tutorial on how to visualize FIDASIM inputs inside an example directory `/u/garciaav/test` with runid `330LT` and `LHD`.

To plot all of your inputs (plasma, fields, distribution function and geometry) execute:
```bash
plot_inputs /u/garciaav/test/ 330LT
```

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

Visualizing your outputs can be done by executing `plot_outputs` found in `lib/scripts/`

Plotting data from your spectra or npa files, and printing the neutron rate can be done in various ways.
Let's start with the easiest method first.
Assuming all of your output files are located in the same folder, the following command will search that folder for all run IDs present.
Then, it will plot all of the spectra and npa data for each run ID on every FIDA and NPA channel.
Furthermore, it will also print out the neutron rate.
```bash
plot_outputs -d /u/garciaav/test/ -s -as -n -an
```
If you are interested in only one shot, e.g. 330LT, add `-r 330LT` to your command.

Let's now consider a more specific example.
What if your files are scattered around in different folders?
What if you are only interested in the active FIDA emission?
What if you only want to view channels 1 and 3?
The example below does this for two different files located in different folders.
```bash
plot_outputs -p /u/garciaav/test/330LT_spectra.h5 /different_path/test_2a_spectra.h5 -f -ls 1 3
```

Now that you have a grasp on how the script works, take a look at the help documentation `plot_outputs -h` to see what else the code is capable of.
