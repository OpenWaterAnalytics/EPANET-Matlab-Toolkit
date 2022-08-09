## Release Notes for EPANET Matlab Toolkit (EMT) v2.2.2

- Performance updates (getComputed..)
- Update [EX24_Parallel_computations.m](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/abe05e2327140c3b9599634295902e36ce43a117)
- Cleanup and some fixes
- When running `epanet('test.inp')` creates a temporary input file using the Matlab function `copyfile` instead of ENsaveinpfile the first time.
- Update [Movie Examples](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/tree/master/examples/movie-example)
- Check out the repository https://github.com/KIOS-Research/CCWI2022-EMT-Tutorial in folder `tutorial`. You can find the short course <b>Indroduction To The EPANET-MATLAB Toolkit For Smart Water Networks Research</b> https://wdsa-ccwi2022.upv.es/introduction-to-the-epanet-matlab-toolkit-for-smart-water-networks-research/
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
- Check [ReleaseNotes2_2_1.md](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/ReleaseNotes2_2_1.md)