# Roadmap / to do

## General todo:
 - [ ] Fix so that `chromflock_gen` passes the gpseq file properly.

## Molecular dynamics, M-step

The M-step is where the model is updated. Initially the model is built from random positions, then at each following step:
 1. The contact indication matrix, W, for each structure is updated (A-step)
 2. New positions, X, are found in the structure that minimized the error functional (which includes W).

Our implementation of the M-step uses Simulated Annealing with a Brownian force, much like what Peter Cook and Davide Marenduzzo does. Alber does not provide any details on what they are doing. However they refer to Russel2012, but not much details can be found their either. Since they published the Nat. Prot. paper the details are at least available in their source code.

### TODO
 - [ ] a lua script should be the only way to set the dynamics parameters.
 - [ ] consider different stopping criteria, for example based on the radius of gyration.
 - [ ] Fix so that mflock crashes gracefully when fRad is specified but no `-r` file is supplied.
 - [ ] Append to log file if it already exists?
 - [ ] Warn if largest `max(R) > 1-r0`. Or convert `R' = R*(1-r0)`? Same goes of the output profiles, should `r` be output or `r*(1-r0)`? (`r0` is the bead radius).
 - [ ] Alber uses k=10 for consecutive beads (and k=1 for non consecutive).
 - [ ] Fix so that the supplied random seed is used.
 - [x] Put live view in main thread and optimization in child in order for SDL to be able to capture events on OSX.
 - [x] Implement ellipsoidal geometry.
 - [x] See if dInteraction = 4 is better than dInteraction = 3. -- seems slightly better. Changing default. See `fconf.dInteraction` in `md.c`.
 - [x] Chimera script to add extra annotations and save automatically -- see `util/chromflock_images.py`

 - [x] Removed any dependencies of GSL -- not necessary any more since only MD is used.
 - [x] Switched RNG and Normal Distribution generator to [A modified ziggurat algorithm for generating exponentially and normally distributed pseudorandom numbers](http://www.tandfonline.com/doi/abs/10.1080/00949655.2015.1060234).
( Previously used Zignor, copied from Jurgen A. Doornik (2005). "An Improved Ziggurat Method to Generate Normal Random Samples" (PDF). Nuffield College, Oxford, together with [PCG](http://www.pcg-random.org) which should be better than a Marsenne Twister in terms of statistical properties and has about the same speed.)
 - [x] Use 3D Gaussian for the Brownian force. (Used a force vector drawn uniformly from a 3D ball before).
 - [x] Long option names, i.e. not just single letter names for the parameters.
 - [x] Consider playing with the radius of the simulation (equivalent with the bead radius), one of the *tricks* used by Alber. -- Faster convergence when the volume quotient of the beads is smaller. However when doing that the beads does not use the full domain very well and the structures look very unrealistic.
 - [ ] Consider to change from double to single precision. The instructions might not be faster for SP, but the memory locality will be better. -- Probably not worthwhile.
 - [x] Any reason to use the varying `k_{con}` described on p. 64? -- Has possibly something todo with the *entropic forces* and convergence speed. -- Se below.
 - [x] Use a decreasing volume exclusion force during the simulations. By hypothesis this reduces the bias towards *old* contacts, i.e., contacts introduced at an early stage. To be specific, the old contacts are the ones closest to the diagonal. This means that the diagonal is well represented among the contacts in the structures but non the off-diagnoal elements. This can be seen as a lack of of sharpness at the edges of the blocks.
 - [x] Write error log to separate file, or use other mechanism to track crashes. -- Added `--joblog parallel.log` to `parallel`.
 - [x] Is there a bias in the brownian force? (that brings the structure to the side) -- there was. Issue resolved and unit test in place.
- [x] Implement Verlet integration. No significant changes in performance.
- [x] How many "brownian" steps and what falloff?
 - [x] Delta t, switched from 1 to 5 2019-04-17. What would be a good value?
 - [x] Implement live view of current state. -- available from SDL with the -a swithch (when compiled with -DSDL).
 - [x] Add a COM force (per chr).
 - [x] Switch to radial data in the "structure summary"
 - [x] Timing does not work in the -D mode.
 - [x] Use a label vector to generate cmms.
 - [x] Include L - label, both for cmm creation and for various per chromosome things.
 - [x] Integrate with GPSeq data.
 - [x] Add comparison of all `err` and all `grad` in the unit tests.
 - [x] Make bead size a parameter or set so that the volume occupancy is 20% (or some other number).
 - [x] Output coordinates for chimera.
 - [x] HASH table to find points in contact.
 - [x] Use a list of pairwise interactions rather than the full matrix.

And for the MD implemenation
 - [ ] Keep track of Ek, and potentially set it to a specific value as in IMP.
 - [ ] For the above point, use 1-diagonal only.
 - [x] Implement a force that brings chrs together to their centre of mass as Alber does.

### Parallelization

 - [x] Structure optimizations can be carried out in parallel, one thread per structure. The jobs are distributed with GNU parallel.

### Memory

## A-step
 - [ ] Consider re-running the m-step until the number of contacts falls below some threshold.
 - [ ] Consider writing coords as raw double data.
 - [ ] Add counting of failed contacts to the `-F` mode.
 - [ ] Add a flag to ignore inter contacts, then there won't be any need to generate a separate A-matrix when trying without inter contacts.
 - [x] Allow diploid experiments.
 - [x] Use zlib to compress the W matrices.
 - [x] Consider *even* assignment in the sense that each structure get the same number of contact initially.
 - [x] Implement re-assignment, and improve the assignment.
 - [x] Possible error when setting the bead radius from volume quotient?!
 - [x] Use zlib to write the cmms.
 - [x] Enable chromflock to use X from previous simulation with updated A.
 - [x] Implement the HiC-normalization from p. 53. The transform was actually already performed on the data!
 - [x] Implement the sphere contact probability matrix from p. 65.
 - [x] Determine W, see p. 66. Rather assign the contacts as in Tjong2016, i.e., the structures that already have the shortest distances will have to do the job.
And for flock:
 - [x] Read sphere contact probability matrix A.
 - [x] Initialize new experiment here or separate? -- in bash script
 - [x] Read and write the contact indicator tensor W.
 - [x] Write radial profile at -F

### Memory
 - Coordinates: Surprisingly not that much memory is required. Let's get an estimate:
 10,000 structure x 5,000 points x 3 dimensions = 1,200 Mb (64 bit floats). Currently aflock is using single precision to represent the coordinates. Using 10,000 structures x 3,300 points requires about .5 Gb.

 - Constraints: Each constraint requires 8 bytes `(uint_32, uint_32)`. By keeping separate lists for separate types of constraints there is no need to story the type.
 10,000 structures x 5,000 constraints = 400 MB.

### Parallelization
 - [x] The assignment process could be parallelized over each contact to be handed out, i.e., the qsorts. -- The generation of the contact distances is now parallelized using `pthread`. Initially all available cps are used.

## General and Random notes

 - [ ] Can we estimate the entropy based on the gradient norm?

 - [ ] Does the simulation tool yield end-to-end distances of L^{3/5} for unperturbated chains?

 - [ ] Update makefile and folder layout, have a look at [stackoverflow](https://stackoverflow.com/questions/7004702/how-can-i-create-a-makefile-for-c-projects-with-src-obj-and-bin-subdirectories)

 - [x] Make Hi-C look like TCC!

 - [ ] Estimate number of failed contacts after each `mflock`-run in order to decide if another round is needed of `aflock` before going to smaller thetas.

 - [ ] ~~libexpat for settings?~~

 - [x] lua for interactions/setting md/sa schedules?

 - [x] The chr4 trick?! -- No

 - [x] Rewrite `chromflock_run` so that it is easier to continue/restart failed jobs. New pipeline: `chromflock` generates `chromflock_gen` which should be edited and run. That script does in turn create `chromflock_run` which is a linear batch of jobs to be run. If anything in `chromflock_run` fails it should be easier to restart from that position.

 - [ ] Are beads of different type handled differently by Alber?

 - [ ] What are the results if only the structures with average number of contact restraints are studied?

 - [x] How many inter vs intra contacts are there when theta = 1, .2, ... , 0.01. (Is the algorithm agnostic to anything but the intra contacts? That would explain random orientation of the chromosomes).

 - [x] What if delta t is increased? -- I think that I got that parameter quite well.

 - PGP is only for diploid cells.

 - Maybe it makes sense to 'kill your darlings', i.e., open up for the possibility to not use all of the initial structures (discard or replace). We don't expect the nucleus to be completely random so random initializations might go wrong.

 - chromflock is snappier than chromoflock. chrof, crof.

 - Convergence critera shouldn't have to be that strict until final theta.

  - A very limited number of theta values, i.e., only few assignment steps (A) are used. This could potentially lead to the inclusion of conflicting constraints. However, it is extremely costly to increase this number...

 - The Hi-C maps does probably implicitly suggest radial positions. However what radial positions that are assigned to each loci depends on the reconstruction method. No reconstruction method is perfect, rather the opposite. This means that GPSeq adds extra constraints, i.e., reduces the space of probable states/configurations.

 - There is an upper limit on the number of neighbours each bead can have. Geometry tells us 12, see the wikipedia page on sphere packings.

## Further devolopment
 - [ ] Use local temperature like [Agrawal](https://www.biorxiv.org/node/96892.full)
 - [ ] Use proteins like [Cook and Marenduzzo](http://dx.doi.org/10.1093/nar/gkw135)
