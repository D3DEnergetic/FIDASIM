title: Atomic Tables

[TOC]

---

# Atomic Cross Sections

As a neutral particle travels through a plasma it undergoes several different types of interactions

* Charge Exchange with Hydrogen and Impurities
* Excitation with Electrons, Hydrogen, and Impurities
* Ionization with Electrons, Hydrogen, and Impurities

These cross sections, as well as Maxwellian averaged reaction rates, are pre-computed over a range of logarithmically spaced collision energies and target temperatures.

# Approximate Hydrogen Charge Exchange Cross Sections 
Some of atomic transitions needed by FIDASIM are not available.
In particular, FIDASIM needs the n/m-resolved charge exchange cross sections.
While certain transitions are available through ADAS [4] others are not, as such, certain approximations are needed to fill out the table.

For instance, we use the equivalence principle (reversibility formula) to mirror the known ADAS cross sections.
$$\sigma(n_f \rightarrow n_i) = \frac{E_i}{E_f} \frac{n_i^2}{n_f^2} \sigma(n_i \rightarrow n_f)$$
This however is insufficient to completly fill out the table. 

Additionally, since the total cross sections for a transition from a given \(n\) to any \(m\) are given by Janev[2] we can then also assume that the probability of a transition from the \(n \rightarrow m\) state decreases exponentially with energy difference between the states we can "spread" the total cross section amoung the different m levels.

@note 
Total cross sections for \(n>4\) are not available so the \(n=4\) total cross sections are used.
Also we normalize the m levels to the Janev tables for consistancy.

A summary of the various approximations used in the charge exchange tables is given in the table below. (Spreading is done over m/rows)
<!-- Charge Exchange table made in http://www.tablesgenerator.com/html_tables -->
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;margin:0 auto;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;}
.tg .tg-amwm{font-weight:bold;text-align:center;vertical-align:top}
.tg .tg-fo0g{font-weight:bold;background-color:#009901;color:#ffffff;text-align:center;vertical-align:top}
.tg .tg-0mq6{font-weight:bold;background-color:#fe0000;color:#ffffff;text-align:center;vertical-align:top}
.tg .tg-bq31{font-weight:bold;background-color:#3166ff;color:#ffffff;text-align:center;vertical-align:top}
</style>
<table class="tg" >
  <caption>H-H Charge Exchange Data Source</caption>
  <tr>
    <th class="tg-amwm">n \ m</th>
    <th class="tg-amwm">1</th>
    <th class="tg-amwm">2</th>
    <th class="tg-amwm">3</th>
    <th class="tg-amwm">4</th>
    <th class="tg-amwm">5</th>
    <th class="tg-amwm">6</th>
    <th class="tg-amwm">Total</th>
  </tr>
  <tr>
    <td class="tg-amwm">1</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <th class="tg-amwm">Janev(n=1)</th>
  </tr>
  <tr>
    <td class="tg-amwm">2</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <th class="tg-amwm">Janev(n=2)</th>
  </tr>
  <tr>
    <td class="tg-amwm">3</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-0mq6">ADAS</td>
    <td class="tg-bq31">Spread</td>
    <th class="tg-amwm">ADAS/Janev(n=3)</th>
  </tr>
  <tr>
    <td class="tg-amwm">4</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <th class="tg-amwm">Janev(n=4)</th>
  </tr>
  <tr>
    <td class="tg-amwm">5</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <th class="tg-amwm">Janev(n=4)</th>
  </tr>
  <tr>
    <td class="tg-amwm">6</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-fo0g">Equivalence</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <td class="tg-bq31">Spread</td>
    <th class="tg-amwm">Janev(n=4)</th>
  </tr>
</table>

#Generating Tables
FIDASIM provides a routine, [[generate_tables(program)]] to generate the atomic tables. 
To generate the atomic tables with the default settings all you need to do is, from the install directory, run

```bash
make atomic_tables
```

@warning
This is computationally expensive so make sure you run this on a computer 
where you won't get angry emails for using up all the CPU's.
Optionally you can add `NTHREADS=#` to the command to set the number of threads.

The default settings should be appropriate for most devices but in some cases it may be necessary to generate custom tables.
For instance, the default tables are calculated assuming the main impurity is Carbon-6 so it would be inappropriate to the default tables if you have a different main impurity.
To generate custom tables from the tables directory run 
```bash
./generate_tables > table_settings.dat
```
to generate the namelist file that contains the default settings. Edit this file to change the settings.

After editing the namelist file run
```bash
./generate_tables table_settings.dat [NUM_THREADS] <-- NUM_THREADS is optional
```
and wait for a couple of hours depending on the number of threads used. 

#Atomic & Nuclear Data References
The atomic data is taken from a variety of sources [1-5]

1. [W.L. Wiese, M.W. Smith, and B.M. Glennon. *Atomic Transition Probabilities. Volume 1. Hydrogen through Neon*. National Bureau of Standards Washington DC Institute for Basic Standards, 1966.](http://www.dtic.mil/dtic/tr/fulltext/u2/634145.pdf)
2. [R.K. Janev, D. Reiter, and  U. Samm. *Collision processes in low-temperature hydrogen plasmas*. Forschungszentrum Jülich, Zentralbibliothek, 2003.](http://www.eirene.de/report_4105.pdf)
3. [M. O'Mullane. *Review of proton impact driven ionisation from the excited levels in neutral hydrogen beams*. ADAS note, 2009.](http://www.adas.ac.uk/notes/adas_c09-01.pdf)
4. [ADAS: Atomic Data and Analysis Structure](http://www.adas.ac.uk/)
5. [R.K. Janev and J.J. Smith. *Cross sections for collision processes of hydrogen atoms with electrons, protons and multiply charged ions.* Atomic and Plasma-Material Interaction Data for Fusion: Volume 4, 1993.](http://www-pub.iaea.org/books/IAEABooks/1839/Atomic-and-Plasma-Material-Interaction-Data-for-Fusion) 
6. [Reinhold, C. O., R. E. Olson, and W. Fritsch. *Excitation of atomic hydrogen by fully stripped ions.* Physical Review A 41.9 1990.](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.41.4837)
7. [Bosch, H-S., and G. M. Hale. *Improved formulas for fusion cross-sections and thermal reactivities.* !+ Nuclear fusion 32.4 1992.](http://iopscience.iop.org/article/10.1088/0029-5515/32/4/I07/meta)

A more precice references can be found in the Fortran function documentation. For example, [[p_cx_janev]]

#Relevent Namelist Settings
* `tables_file`: Location of atomic tables file

#Fortran References
* [[atomic_tables(module)]]: Module containing routines to calculate atomic tables and reaction rates
* [[generate_tables(program)]]: Program to create atomic tables file
* [[AtomicCrossSection]]: Defines a n/m-resolved atomic cross section table
* [[AtomicRates]]: Defines a n/m-resolved atomic cross section table
* [[AtomicTransitions]]: Defines an atomic table for populating and de-populating reaction rates
* [[AtomicTables]]: Atomic tables for various types of interactions need by FIDASIM
* [[NuclearRates]]: Defines nuclear fusion reaction rate table
* [[read_atomic_cross]]: Reads atomic cross section from file
* [[read_atomic_rate]]: Reads atomic rates from file
* [[read_atomic_transitions]]: Reads in a atomic transitions table from file
* [[read_nuclear_rates]]: Reads in nuclear reaction rates from file
* [[read_tables]]: Reads all cross sections and rates needed by FIDASIM
* [[m_spread]]: Spreads total n cross section amoung m states

#Hydrogen-Hydrogen Interactions

##\(H^+ + H(n) \rightarrow H(m) + H^+\)

![](|media|/H_H_cx_1_m.svg "H-H CX n=1 to m=1..6"){: width="400"} ![](|media|/H_H_cx_2_m.svg "H-H CX n=2 to m=1..6"){: width="400"}
{: style="text-align: center"}

![](|media|/H_H_cx_3_m.svg "H-H CX n=3 to m=1..6"){: width="400"} ![](|media|/H_H_cx_4_m.svg "H-H CX n=4 to m=1..6"){: width="400"}
{: style="text-align: center"}

##\(H^+ + H(n) \rightarrow H^+ + H(m), \, m \gt n\)

![](|media|/H_H_excit_1_m.svg "H-H Excitation n=1 to m=2..6"){: width="400"} ![](|media|/H_H_excit_2_m.svg "H-H Excitation n=2 to m=3..6"){: width="400"}
{: style="text-align: center"}

![](|media|/H_H_excit_3_m.svg "H-H Excitation n=3 to m=4..6"){: width="400"} ![](|media|/H_H_excit_4_m.svg "H-H Excitation n=4 to m=5..6"){: width="400"}
{: style="text-align: center"}

##\(H^+ + H(n) \rightarrow H^+ + H^+ + e\)

![](|media|/H_H_ioniz.svg "H-H Ionization"){: width="400"}
{: style="text-align: center"}

#Hydrogen-Electron Interactions

##\(e + H(n) \rightarrow e + H(m),\, m \gt n\)

![](|media|/H_e_excit_1_m.svg "H-e Excitation n=1 to m=2..6"){: width="400"} ![](|media|/H_e_excit_2_m.svg "H-e Excitation n=2 to m=3..6"){: width="400"}
{: style="text-align: center"}

![](|media|/H_e_excit_3_m.svg "H-e Excitation n=3 to m=4..6"){: width="400"} ![](|media|/H_e_excit_4_m.svg "H-e Excitation n=4 to m=5..6"){: width="400"}
{: style="text-align: center"}

##\(e + H(n) \rightarrow e + H^+ + e\)

![](|media|/H_e_ioniz.svg "H-H Ionization"){: width="400"}
{: style="text-align: center"}

#Hydrogen-Carbon₆ Interactions

##\(C^{6+} + H(n) \rightarrow C^{5+} + H^+\)

![](|media|/H_C6_cx.svg "H-H CX"){: width="400"}
{: style="text-align: center"}

##\(C^{6+} + H(n) \rightarrow C^{6+} + H(m), \, m \gt n\)

![](|media|/H_C6_excit_1_m.svg "H-C6 Excitation n=1 to m=2..6"){: width="400"} ![](|media|/H_C6_excit_2_m.svg "H-C6 Excitation n=2 to m=3..6"){: width="400"}
{: style="text-align: center"}

![](|media|/H_C6_excit_3_m.svg "H-C6 Excitation n=3 to m=4..6"){: width="400"} ![](|media|/H_C6_excit_4_m.svg "H-C6 Excitation n=4 to m=5..6"){: width="400"}
{: style="text-align: center"}

##\(C^{6+} + H(n) \rightarrow C^{6+} + H^+ + e\)

![](|media|/H_C6_ioniz.svg "H-C6 Ionization"){: width="400"}
{: style="text-align: center"}

