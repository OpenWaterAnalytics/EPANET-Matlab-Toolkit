## Release Notes for EPANET Matlab Toolkit (EMT) v2.2.2

- Performance updates (getComputed..)
- Update [EX24_Parallel_computations.m](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/abe05e2327140c3b9599634295902e36ce43a117)
- Cleanup and some fixes
- When running `epanet('test.inp')` creates a temporary input file using the Matlab function `copyfile` instead of ENsaveinpfile the first time.
- Update [Movie Examples](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/tree/master/examples/movie-example)
- We have prepared a tutorial as part of the CCWI-WDSA 2022 conference. You can find the files in the `tutorial` folder at the following [GitHub repository](https://github.com/KIOS-Research/CCWI2022-EMT-Tutorial).
  - Basic functionality
  - Run first EPANET analysis
  - Run first MSX analysis
  - How to include uncertainties in simulations
  - How to analyze network graph
  - How to place pressure sensors
  - How to create leakage events
  - How to detect leakage events
  - How to place quality sensors
  - How to create contamination events (using MSX)
  - How to detect contamination events