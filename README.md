Agent-Based Plasmid Simulations
===

This is the implementation of the evolutionary agent-based model of
hosts and plasmids used in the following research paper :

K. Kentzoglanakis, D. G. Lopez, S. P. Brown, and R. A. Goldstein. [The
Evolution of Collective Restraint: Policing and Obedience among
Non-conjugative
Plasmids.](http://dx.plos.org/10.1371/journal.pcbi.1003036) PLoS
Comput Biol, 9(4):e1003036+, 2013.


## OVERVIEW

The program `dps` that implements the actual model is written in C and
the source code is located in `src/`. Also provided are R scripts
located in `R` that deal with post-processing the results of `dps`,
plotting etc. The files `R/paper_*.R` plot the figures for specific
research articles (more about that below). 



## REQUIREMENTS

A simple `Makefile` is included for generating the `dps` program. Run
`make` in the dps/src directory in order to compile the code. `dps`
depends upon a few external libraries which should be present in your
system along with their respective header (development) files. On
Linux, these dependencies can be installed using the designated
package manager. On Mac systems, installation via
[homebrew](http://brew.sh) is advised.

The prerequisites for compiling `dps` are as follows:

1. [GNU GCC](http://gcc.gnu.org) >= 4.4.5

2. [GNU Scientific Library](http://www.gnu.org/software/gsl) >= 1.15

3. [HDF5](http://www.hdfgroup.org/HDF5) 

4. [GLib](https://developer.gnome.org/glib/) >= 2.24.2



## USAGE

The general usage pattern for `dps` is:

    $ dps [OPTIONS] out.h5

where `out.h5` is the name of the output file (see next section for
more details on that). Use `dps -h` for a list of options that can be
specified in the command line. The most important options include :

* `beta`, `kappa`, `alpha` : override the default initial value of one
  or more plasmid replication parameters. Defaults: β=0.05, κ=0, α=0.

* `mutate` : use any combination of `b`, `k`, `a` in order to activate
mutations on β, κ and α respectively. E.g. `--mutate bk --alpha 1`
sets α=1 and activates mutations on β and κ. Defaults: no mutations
and no copy number control (i.e. beta=0.05, kappa=0, alpha=0).

* `mu` : specifies the probability of mutation per plasmid replication
      event (default: 5e-3).

* `mut_rng` : specifies the width / 2 of the uniform distribution
        around a plasmid's current parameter values used for mutations
        (default: 0.05).

* `pconj` : specifies the probability of a successful horizontal
      transmission event per donor host at a given time step. Use
      `--pconj 0` (default) to switch off conjugation.

* `steps` : how many steps to run the simulation for. The `SIGINT`
      signal (C-c C-c) is taken to signify a user-requested premature
      end to the simulation and is handled gracefully by the program.

* `psize` : specifies the host population size (default: 1000)

* `load_from` : specify an input file (which could be the output of a
     separate simulation) from which the population will be loaded and
     usedas the initial population for the current simulation.

* `compete` : activates the competition mode with mutations turned
  off. This is used together with `--load_from` in order to set up a
  competition scenario by specifying an initial population that
  contains (exactly) two plasmid profiles. If one of the profiles goes
  extinct then the simulation stops, otherwise the simulation
  continues until the maximum number of steps is reached.

* `seg_type` : specifies the type of plasmid segregation during cell
  division, either "binomial" (default) or "perfect"

* `max_{beta,kappa,alpha}` : specifies the upper limits for the plasmid
  replication parameters

* `fparts` : how many distinct parts to split the histograms into, in
  order to capture frequencies at distinct stages of the simulation.


## OUTPUT


The simulation's output is stored in a file following the
[HDF format](http://www.hdfgroup.org/HDF5/); the file name (typically
ending in `.h5`) must be specified as the last argument to the `dps`
command. There exists a Java tool for viewing HDF files, called
[HDFView](http://www.hdfgroup.org/products/java/hdfview/index.html),
which can be useful for a quick inspection of a simulation's output. 

The HDF output file has four major sections :

1. `dynamics` : contains the evolutionary dynamics of the simulation
split in four groups

    * `counters` : records the numbers of various events (such as
      population size, total copy number, division events etc.) for
      each time step in the simulation and has the following columns:
	  * `n` : the number of hosts in the population
	  * `cn` : the total number of plasmids across all hosts
	  * `inf` : the number of plasmid-infected hosts
	  * `ptypes` : the number of distinct plasmid types in the population
	  * `loss` : the number of segregational losses across cell
        division events
      * `div.inf` : the number of division events of plasmid-infected
        hosts
      * `div` : the number of division events of all hosts
	  * `death` : the number of host deaths
	  * `rep` : the number of plasmid replication events
	  * `ht` : the number of horizontal transmission events
	  * `mut` : the number of mutation events

	
    * `intra`, `inter` and `global` : these contains descriptive
      statistics about evolutionary variables (i.e. means in subgroup
      `M`, variances in subgroup `V` and pairwise covariances in
      subgroup `C`) at three different levels: within hosts (intra),
      between hosts (inter) and across all plasmids regardless of
      hosts (global). The `M` and `V` subgroups have the same columns,
      while the `C` subgroup has column names of the form "$A.$B"
      where $A and $B enumerate the columns of the corresponding `M`
      group. For example, if the `M` subgroup has the columns "a", "b"
      and "c", then the `C` subgroup will have the columns "a.b",
      "a.c" and "b.c" with the corresponding covariances.

      The `M` subgroup of the `intra` and `global` members has the
      following columns:
	  * `cn` : average per-host copy number
	  * `beta`, `kappa`, `alpha` : average values of the plasmid
        replication parameters
      * `nr` : average number of plasmid replication events per host
	  * `ht` : average number of horizontal transmission events per
        host
      * `death` : average number of plasmid deaths per host
	  * `fitness` : average value of plasmid fitness
	  * `t{beta,kappa,alpha}` : transmission biases of the plasmid
        replication parameters

      The `M` subgroup of the `inter` member has the following
      columns:
	  * `domg` : the average value of ΔΩ across hosts
	  * `cn` : the average plasmid copy number per host
	  * `dev` : the average deviation of the host's copy number from
        the optimal copy number
      * `beta`, `kappa`, `alpha`, `nr`, `ht`, `death`, `fitness`
        `t{beta,kappa,alpha}` : see above (in this case per-host
        averages)

    * `relatedness` : contains the groups `wg` (for calculating
      whole-group relatedness) and `oo` (for calculating others-only
      relatedness). Each of these groups contains a `cov` (covariance)
      and a `var` (variance) data frame with members `beta`, `kappa`
      and `alpha`. In order to calculate the dynamics of, say, the
      others-only relatedness coefficient for β use:
      `r$dynamics$relatedness$oo$cov$beta /
      r$dynamics$relatedness$oo$var$beta`

2. `histograms` : contains the copy number and cell age (i.e. number
of simulation steps required for a host to divide) histograms. The
bins are stored in dataset `x`, whereas the counts in dataset `y`
which has dimensionality `(fparts x max_cn)`, where `fparts` and
`max_cn` are arguments to `dps`. `histograms` also contains all the
joint histograms between β, κ, and α in the respective groups `bk`,
`ba` and `ka`. In the case of the joint distributions, the bins are
stored in datasets `x` and `y`, whereas the counts in dataset `z`
which has dimensionality (`fparts` x `nbins`+1 x `nbins`+1), where
`nbins` is an argument to `dps`.

3. `settings` : provides access to the simulation's parameter values. 

4. `population` : contains the state of the plasmid (group `profiles`)
and host (group `hosts`) population at the end of the simulation as
nested-parentheses strings. The format of these strings are as
follows:

    * `profiles` : specifies the distinct plasmid profiles that exist
      in the population of which individual plasmids are
      incarnations. A profile is specified by the string `(ID X BETA
      KAPPA ALPHA X X X)`, where `ID` is a unique integer identifier,
      `BETA`, `KAPPA`, `ALPHA` are the values of the profiles
      replication parameters and `X` denotes fields we do not care
      about when loading the profile in a simulation using
      `--load_from` (these fields are determined automatically by the
      program). The string of all `profiles` in a population is of the
      form `(PROF1 PROF2 ... PROFN)`, i.e. a sequence of individual
      profile strings enclosed within parentheses. 

   * `hosts` : specifies all the hosts in the population along with
     their contained plasmids. A single host is specified by the
     string `(OMEGA OMEGA_0 AGE X X PLASMID_POOL)`, where the
     `PLASMID_POOL` enumerates all the distinct plasmid profiles in
     the host along with their corresponding copy number as follows:
     `((ID1 CN1) (ID2 CN2) ...)`, where `ID` refers to the profile's
     unique integer identifier (which should match the one in the
     `profiles` string) and `CN` is the number of profile copies
     within the host.

5. `competition` : this group is generated only when the competition
   mode between two plasmid profiles is activated (using `--compete`
   or `-c`). It contains the subgroups `contenders` (which gives
   information on the two competing profiles) and `frequencies` which
   tracks the total copy number of each plasmid profile over time.



## EXAMPLES

Run a NO-CNC simulation with no conjugation:

    dps --mutate b --kappa 0 --alpha 0 --pconj 0 --steps 50000 out0.h5

Run simulation where CNC evolves from NO-CNC with the migration
probability set to 0.01:

    dps --mutate bka --pconj 0.01 --steps 100000 out1.h5

Continue the previous simulation for another 100000 steps from where
it left off:

    dps --load_from out1.h5 --mutate bka --pconj 0.01 --steps 100000 out2.h5


## POST-PROCESSING

The R source file `R/dps.r` is provided as a convenient template for
viewing/processing/plotting the output of the simulation (HDF
file). Its requirements include `hdf5`, `multicore`, `tseries` and
`gplots`.

In R, source `dps.r` and use

    results <- dps.load("results.h5")

in order to load the simulation results located in file `results.h5`
into R object `results`. 

Having loaded the results, use

    dps.analyze(results)

to plot a basic view of the simulation. 

The files `R/paper_*.r` concern the plotting of results for specific
research publications 



