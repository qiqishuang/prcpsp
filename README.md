Preemptive Resource Constrained Scheduling Problem (PRCPSP)
==================================================

A formulation for the PRCPSP which I came up with `prcpsp.py`. 

It uses the concept of events which reduces the number of variables significantly.
See `prcpsp.pdf` for more details.

Implemented here, `lcalg.py`, is a local constraint programming algorithm for the estimation of earliest and lastest start times of acitivities which allows the problem to be solved efficiently.
