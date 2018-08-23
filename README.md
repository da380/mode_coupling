# mode_coupling


############################################
#                  overview                #
############################################

Calculation of spectra is done in several stages,
and these will be described in turn. The basic
steps are:

1) Select a group of multiplets to use in the
coupling calculations;

2) Form the coupling matrices for a given 3D
earth model;

3) Calculate the spectra for a given source
and a number of seismic stations.


#############################################
#        Forming a list of multiplets       #
#############################################

The program "mdcpl" that calculates the coupling
matrices reads in a text file called modes.in.
This is a list of all the multiplets to be
used in the coupling calculations. The easiest
way to form such a list is to us the small
program "mode_list". This routine is given
a minimum and maximum frequency, and makes
a list of all multiplets with eigenfrequcies
within this range. Such a list can also be
made by hand, and the required format can be
seen by looking at an output from "mode_list".


#############################################
#      Building the coupling matrices       #
#############################################

This calculation is done by the program "mdcpl".
When this program is run it asks for following
inputs from the command line:

>> Input mins & maxs  eg:0 8 

This is the minimum and maximum structural degrees
to include from the given 3D model (by default S20RTS
at the momement). For example, if you want all the
variations in S20RTS, you would type

>> 0 20

but if you only wanted to include the degree two
heterogeneity you would input

>> 2 2

Next we are asked:

.. do you want include rotation? yes=1 no=0

which is clear

Then we get:

>> do you want include ellip.? yes=1 no=0

Here this means accounting for elliptical flattening
of the model's boundaries along with degree two perturbations
to the volumetric parametres due to internal level surfaces
also becoming elipsoidal.

The default should be to always include both rotation and
ellipticity.

Finally, we are asked:

>> afile = 

which is to be the name for the output binary file containing the
mode coupling matrices.

Note that a script "build_coupling" can be used to form a mode
list and build the coupling matrices more quickly. This is
run as:

>> build_coupling wmin wmax smin smax

and this includes all modes in the range wmin to wmax, and
all structures between smin and smax. As an example, we could
run:

>> build_coupling 0 3 0 20

which would be all modes below 3mHz, and include the
full S20RTS model



#############################################
#          calculating the spectra          #
#############################################

Once the coupling matrices have been made for a
given earth model, we are ready to compute
spectra. This is done for one earthquake at a
time, but many seismic stations can be added
within each calculation.

The source information is contained within a
simple text file. The following is an example
of its content (from the file "dspec.source"):

-13.8199997
 -67.25
  647.099976
 -7.58999992E+27
  7.75000019E+27
 -1.59999993E+26
 -2.50000004E+28
  4.20000007E+26
 -2.48000005E+27

The first line in the source latitude in degrees
The second line the the source longitude in degrees
The third linear is the source depth in km
The next six lines are the components of the
moment tensor in the same order as the CMT catalog.
I think the units here are in Nm. But this should not
be hard to check. For reference the above example is
for the 1994 Bolivia event.

Note that a step function in time is assumed for the
source at the moment. A more general source time function
can be added easily enough by convolving the results
appropriately. 

The station information is listed in a text file,
with the following ("dspec.rec") being an example:

30
49.1448   12.8803   -01.000   00.000   00.000  STAT1.Z
49.1448   12.8803    00.000   01.000   00.000  STAT1.N
49.1448   12.8803    00.000   00.000  -01.000  STAT1.E
49.1449   12.8782   -01.000   00.000   00.000  STAT2.Z


where I have only included the first five lines. The first
line contains the number of spectra to calculate, and there
is a line in the file for each one. Note that different
components are the same location are counted individually.
In each line there are six inputs:

latitude longitude cmp1 cpm2 cpm3 suffix

cmp1, cmp2, cmp3 are the components of a source polarisation
vector. The final entry is a string that will identify the
corresponding spectra and time series files. Note that
in this example the first three lines are the vertical (Z),
north (N), and east (E) components at the same station. And
so here you can see the appropriate polarisation vectors needed
to do this.

Once the source and receiver information is set, the spectra
are most easily calculated using the script "run_dspec"
the following is an annotated example:


================================================================================

#!/bin/bash

# set the number of processors
njob=4                              --- this is the number of processors used for the main calculations (can be equal to one)

dspec_pre << !
dspec.source                        --- name of the source file
dspec.rec                           ---  name of the receiver file
0.1                                 ---  minimum frequency for the output spectra (in mHz)
2                                   ---  maximum frequency for the output spectra (in mHz)
20                                  ---  time step for the time-series (in seconds)
128                                 ---  length of time series (in hours) 
0.01                                ---  fractional width of a filter in the frequency domain (this is a sensible value)
0.05                                ---  a parameter used in the IDSM (this is a sensible value)
0                                   ---  start time for a cosine bell window
128                                 ---  end time for a cosine bell window
$njob
matrix_parts.bin                    --- name of the coupling matrix file
!


for i in $(seq 1 $njob)
do
  dspec_cal $i > dspec_cal.out.$i &
done
wait

dspec_pro << !
$njob
dspec.out                           --- prefix for the spectra files (e.g. one might be dspec.out.STAT1.Z)
dspec.ts                            --- prefix for the time-series files (e.g. one might be dspec.ts.STAT1.Z)
!


================================================================================

Some comments will be useful:

1) the minimum frequency must always be non-zero. This is due to the form of the attenuation being
singular at zero frequency. In practice 0.1 mHz is below anything we actually observe seismically.

2) The spectra close to the ends of the frequency limits have been flitered slightly (this is what the
0.01 fractional filter does). Moreover, near the upper limit there might well be truncation effects to
be aware of -- in general to look at spectra in a given range, you should probably couple all modes
up to about 1mHz about the upper frequency of interest.

3) The calculations are very naively parallelised. This can be done properly on a real cluster as a batch job.
For large scale calculations this is useful to do.

4) A cosine bell window is applied in the time-domain prior to computing the final spectra. This is fairly standard
in normal mode work.

5) The output spectra files are text with four lines:

(frequency)  (real part of acceleration spectra) (complex part of acceleration spectra) (modulus of acceleration spectra)

The code can be made to output displacement or velocity spectra if needed -- its just an additional line of code as its all
done in the frequency domain anyway

