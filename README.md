# pimaimCageCF

pimaimCageCF is a tool to compute cage correlation functions from output files of pimaim. pimaim is a molecular dynamics simulation package with full ion polarization effects.

## Authors

The original pimaimCageCF was written by Mathieu Salanne, UPMC Univ Paris 06, UMR PECSA, Paris, France.
pimaimCageCF as proposed here has been partially rewritten by Maximilien Levesque, Postdoc under the supervision of Mathieu Salanne.

## Installation

pimaimCageCF has been written with linux in mind. It should work under MacOS, and certainly using Cygwin or other terminal emulator under Microsoft Windows.
1/ You need to have `scons` installed on your computer. `scons` is a `GNU make` made simple and intelligent. Search the web for more informations.

2/ Once `scons` is installed, just type in your terminal :
```
$ scons
```

## Licence

As we still don't have a clear policy for our software, you **HAVE** to ask all the authors their permission before using it.
The program comes with no warranty. It is intended for academic research, so you should verify it is correct (which also means correct for your application) before publishing anything with it.

## How to use pimaimCageCF

You'll need a text file named `cagecf.inpt`.

### Example file for cagecf.inpt

Order is important. You can keep comments.

    2           Number of species, 2 for B, O
    480         Number of ions of species 1. (B)
    320         Number of ions of species 2. (O)
    50          N(skip)
    43.659      Size of the cubic cell in Bohr (u.a.) i.e. in the natural unit of PIMAIM. May be found by a simple tail -n3 restart.dat
    1           How many files to be read from? pimaimCageCF should work with several position files, but it is not tested yet.
    /tmp/example/positions.out
    5000        Number of coordinate sets in file 1. If you have a MD of 1000000 steps of 41.341 a.u (ie you simulate 1ns) and print velocities every 10 steps, you have 100000 sets of coorinates in the file containing velocities.
    100         Time-span of the correlation function (in steps) You decide to compute your cage correlation function of how many MD steps ?
    3.5       cutoff for bond. First a first coordination shell cage correlation function, give the distance to the first minimum in the radial distribution fun.
    1           number of ions for change in cage
    41.341    MD timestep. The one introduced in runtime.inpt. Should be 41.341 if your MD step is 1femtosecond.
    1000        number of configs between MD dumps. Called "Number of steps inbetween periodic output (velocities etc)." in runtime.inpt


