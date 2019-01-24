Chlorine Decay in EPANET Example Network Net2
=============================================

This example models bulk and wall chlorine demand in the standard EPANET
example network Net2. The EPANET input file (net2-cl2.inp) can be run
under EPANET to compare computed chlorine results against those produced
by EPANET-MSX. Because the MSX input file uses the turbulent mass transfer
coefficient equation for both turbulent and laminar flow regimes, the
results between the two models will be slightly different. The files
included are:

  net2-cl2.inp  --  the EPANET network input file with a specification
                    of the bulk and wall chlorine decay process
  net2-cl2.msx  --  the MSX reaction file with its own, equivalent
                    representation of bulk and wall chlorine decay
  net2-cl2.rpt  --  the report file containing the MSX results for
                    selected nodes
  net2-cl2.xls  --  an Excel file that compares the EPANET and MSX
                    solutions for selected nodes.

Chlorine Decay in EPANET Example Network Net3
=============================================

This example models bulk chlorine demand in the standard EPANET
example network Net3. The two sources have distinctly different
chlorine first order reaction rate constants, and these are
modeled in the blended water by an assumed linear relationship
between the overall reaction rate constant and the fraction of
water from both sources. The files included are:

  net3.inp      --  the EPANET network input file for hydraulic
                    calculations
  net3-cl2-compiler.msx  --  the MSX reaction file with its
                    representation of bulk chlorine decay
  net3-cl2-compiler.rpt  --  the report file containing the MSX results for
                    selected nodes
