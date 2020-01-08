
This bundle contains the implementation of the Fermi-Hubbard model
solver via the determinant diagrammatic Monte Carlo method.
package is available from

<http://montecarlo.csi.cuny.edu/umass/fermion.html>.

The algorithm itself is detailed in the preprint cond-mat/060530,
available at

<http://arxiv.org/PS_cache/cond-mat/pdf/0605/0605350.pdf>,

thus only technical details are considered here.
Questions, comments, criticism, bug complaints and whatnot are
welcome at <burovski@physics.umass.edu>.

The code is written in the Fortran 90 language, and heavily uses
BLAS and LAPACK libraries. It also comes in two versions, MPI
and non-MPI.

*************** There are five source files: *********************

1. fermi-hubbard.f     --- the main source code file
2. fermi-hubbard-MPI.f --- MPI version
3. event_mod.f         --- supplementary tree search code
4. rndm_mumbers.f      --- random number generator code
5. det_n2.f            --- ``fast-update'' linear algebra routines


******************* Compiling the code: **************************

In order to compile the code one has to invoke BLAS and LAPACK
libraries, for which the general-purpose versions are available
from, e.g. www.netlib.org.

Nevertheless, it is highly advisable that you use an implementation
of these libraries optimized for your machine, for the performance
of the code crucially depends on the linear algebra library. It is
also highly advisable that you turn maximum code optimization
and do at least some profiling. Remember, the computational
complexity of a single update scales quadratically with number
of particles, thus if you, for example, have 10 fermions, the
simulation will be at least 100 times slower than what you would
have for bosons!

The general compile line is

<Fortran compiler> <optimization switches> <BLAS/LAPACK switches>
          det_n2.f event_mod.f rndm_mumbers.f fermi-hubbard.f

For example, on the Opteron machine with PGI compiler this
looks like:

pgf90 -fastsse -lacml det_n2.f event_mod.f rndm_mumbers.f
                     fermi-hubbard.f

On Xeon machines with PGI compiler use instead:

pgf90 -fast -llapack -lblas det_n2.f event_mod.f
              rndm_mumbers.f fermi-hubbard.f


and so on.

**************** Input parameter files: ***************************

The I/O protocol is devised for massively parallel runs
[few months of wallclock time on 200 processors, to make an
example :)]

First thing the code reads the file named infile. This should have
the following structure:

nfls
nc_1      par_XXX
nc_2      par_YYY
 ........
nc_nfls   par_ZZZ

Here nfls is the number of different parameter files, and nc_i is
the number of clones to run per i-th parameter set. The latter
number must match the corresponding entry in the parameter file
(see below).

For a parameter file par_XXX, _XXX is used as a suffix for all the
output: stat_XXX_#.dat, out_XXX_#.out etc, where # is the number of
a clone (1...nc_i).

For a non-MPI version simply keep nfls=nc_1=1 :).

Parameter files have the following structure (see the sample file
in the bundle):

3     ! Spatial dimension
4     ! Linear system size in x-, y-, and z-directions
4
4
0   ! 0 if new configuration, 1 if read one [from cnf_suffix_#.dat]
0   ! 0 if new statistics,    1 if read one [from stat_suffix_#.dat]
0   ! 0 if new rndm() seed [below], otherwise read
    ! the generator state [rndm_suffix_#.dat]
1   ! number of clones; must match the corresponding entry
    ! in the infile
-5.0d0              ! \mu, chemical potential; half filling
                    ! corresponds to U/2, the bottom of the
                    ! non-interacting band is at -6
7.915d0  0.015d0    ! - U, interaction and its initial value for
                    ! thermalization, which  starts with this
                    ! [presumably small] value of U, and gradually
                    ! adjusts it to avoid quenching metastable
                    ! configurations;
0.3d0               ! \beta, the inverse temperature
5.1d0               ! eta: G- vs Z-sector weight, adjust it manually
                    ! so that simulation spends roughly half of the
                    ! time in each sector
500                 ! Number of \tau points for GF tabulation mesh:
1.d2                ! tolerance level for determinant recalculation:
0  0.1d0            ! x- and \tau- radii for create/annihilate
                    ! ira&masha updates
1  0.1d0            ! x- and \tau- radii for leap_add/drop updates
50                  ! Number of \tau points for seeding the 'events'
                    ! in lead_add/drop updates
1d2                 ! number of 'sweeps' for thermalization,
                    ! a 'sweep' is defined as \beta*U*L**3
2.d5                ! # of MC steps for printing
6.d5                ! # of MC steps for writing to disk
1.d0                ! # of MC steps for measuring
5.7                 ! Wallclock time limit, hrs
------
0.0       ! address probablities for updates: add/drop a vertex
0.1       ! create/annihilate ira&masha
0.35      ! leap_add/drop
0.1       ! hop
------
4836 2738 ! rndm() seeds. The number of seed must match
4836 2748 ! the number of clones above.


********************* Output files: ******************************

Once again, the code has been written with massively parallel use
in mind, therefore no output goes directly to stdout. Instead, the
following files are written:

1. cnf_suffix_#.dat  --- MC configuration [not human-readable]
2. stat_suffix_#.dat --- MC statistics [not human-readable]
3. out_suffix_#      --- output [only human-readable]

These files are written by every process in the group [i.e.
if you request 10 clones in the parameter file, you will have
10 different statistics, configurations and outputs.]

In addition, you will have few 'service' files, one per group:

4. an_suffix  --- this will list the names of the statistics files
                  for the group. I feed this file to the service
                  program that merges statistics together.
5. par_suffix__ --- the 'rerun' parameter file. I.e. if you start
                    simulations from scratch, you have the
                    'initial' parameter file which says:
                    'don't read configuration/statistics and do
                    thermalization'. This 'rerun' file has these
                    switches off, i.e. if you want to continue
                    a simulation, simply replace your original
                    parameter file with this one.

Also, for the verification/debugging purposes, the code emits
three more files:

6. nmnm_suffix_#.dat --- the vertex number distribution
7. ct3d_suffix_#.dat --- the distribution of the temporal
                         distance between ira and masha
                         [integrated over the coordinates]
8. det_suffix_#.dat  --- log10( determinant ratio ) distribution
                         if you have a lot of large or small
                         det. ratios, it's a sign that something
                         might be wrong with the simulation.

***************** A test case ************************************

Also included in this bundle is the test case of a fairly small
system at fairly high temperature:

par_tst --- the parameter file to run
gim.dat --- ct3d_tst.dat file from your simulation
            must match this one
Vsmall.dat --- nmnm_tst.dat must match this one
res_tst.txt    --- this is what the simulation should converge to
