Usage: arg-sample [OPTION]

  -s,--sites  <sites alignment>
    sequence alignment in sites format

  -f,--fasta  <fasta alignment>
    sequence alignment in FASTA format

  -o,--output  <output prefix>
    prefix for all output filenames (default='arg-sample')

  -a,--arg  <SMC file>
    initial ARG file (*.smc) for resampling (optional)

  --region  <start>-<end>
    sample ARG for only a region of the sites (optional)

  --maskmap  <sites mask>
    mask map file (optional)

Model parameters
  -N,--popsize  <population size>
    effective population size (default=1e4)

  -m,--mutrate  <mutation rate>
    mutations per site per generation (default=2.5e-8)

  -r,--recombrate  <recombination rate>
    recombination per site per generation (default=1.5e-8)

  -t,--ntimes  <ntimes>
    number of time points (default=20)

  --maxtime  <maxtime>
    maximum time point in generations (default=200e3)

  --time-step  <time>
    linear time step in generations (optional)

  --popsize-file  <popsize filename>
    file containing population sizes for each time span (optional)

  --times-file  <times filename>
    file containing time points (optional)

  -M,--mutmap  <mutation rate map file>
    mutation map file (optional)

  -R,--recombmap  <recombination rate map file>
    recombination map file (optional)

Sampling
  -n,--iters  <# of iterations>
    (default=1000)

  --resample-region  <start>-<end>
    region to resample of input ARG (optional)

  --resume
    resume a previous run

  --overwrite
    force an overwrite of a previous run

Miscellaneous
  -c,--compress-seq  <compression factor>
    alignment compression factor (default=1)

  --climb  <# of climb iterations>
    (default=0)

  --sample-step  <sample step size>
    number of iterations between steps (default=10)

  --no-compress-output
    do not use compressed output

  -x,--randseed  <random seed>
    seed for random number generator (default=current time)

  --gibbs
    use Gibbs sampling

Information
  -V,--verbose  <verbosity level>
    verbosity level 0=quiet, 1=low, 2=medium, 3=high

  -q,--quiet
    suppress logging to stderr

  -v,--version
    display version information

  -h,--help
    display help information

  --help-advanced
    display help information about advanced options
