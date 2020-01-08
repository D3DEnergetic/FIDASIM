title: Visualization

#Visualizing Inputs and Outputs

[TOC]

---

##Visualization: Inputs

Visualizing your inputs can be done by executing `plot_inputs` found in `lib/scripts/`

Depending on what you wish to plot, your inputs, geometry, equilibrium and/or distribution files will need to be located in the same folder.
Below are brief descriptions and examples of what the script can currently handle:

To plot all of your inputs, simply indicate the directory and run ID.
```bash
plot_inputs /p/fida/lstagner/TEST/ test_1a
```

To plot only the beam and diagnostic geometry, append the optional argument -g
```bash
plot_inputs /p/fida/lstagner/TEST/ test_1a -g
```
In a similar fashion, append -p, -f and/or -d to plot the plasma, fields and/or distribution function inputs, respectively.

If you are plotting many FIDA or NPA line of sights, then it might be beneficial for you to append -l to remove the legend from the 3D plot.

If you wish to plot lineouts on your figures, simply indicate the value and dimension you want to cut through.
For example, if you are interested in seeing what the plasma lineout looks like at R = 170 cm along the z axis, execute the following command
```bash
plot_inputs /p/fida/lstagner/TEST/ test_1a -p -rz 170
```

There are many more possible lineouts that can be viewed, so run `plot_inputs -h` to look at the help documentation.

##Visualization: Outputs

Visualizing your outputs can be done by executing `plot_outputs` found in `lib/scripts/`

Plotting data from your spectra or npa files, and printing the neutron rate can be done in various ways.
Let's start with the easiest method first.
Assuming all of your output files are located in the same folder, the following command will search that folder for all run IDs present.
Then, it will plot all of the spectra and npa data for each run ID on every FIDA and NPA channel.
Furthermore, it will also print out the neutron rate.
```bash
plot_outputs -d /p/fida/lstagner/TEST/ -s -as -n -an
```
If you are interested in only one shot, e.g. test_1a, add `-r test_1a` to your command.

Let's now consider a more specific example.
What if your files are scattered around in different folders?
What if you are only interested in the active FIDA emission?
What if you only want to view channels 1 and 3?
The example below does this for two different files located in different folders.
```bash
plot_outputs -p /p/fida/lstagner/TEST/test_1a_spectra.h5 /different_path/test_2a_spectra.h5 -f -ls 1 3
```

Now that you have a grasp on how the script works, take a look at the help documentation `plot_outputs -h` to see what else the code is capable of.
