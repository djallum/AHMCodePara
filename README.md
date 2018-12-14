# AHMCodePara
Modules for the parallel code

Diag.f90: Code used to determine the diagonalize the hamiltonian and calculate the DOS
DiagDIRECT.f90: Same as above but optimized for a direct ensemble
DiagIGNORE.f90: Not really useful. Allows for a quantity I don’t use anymore. Should be deleted.
DOSsetupPARA.f90: Code used to set up each system and sending the clusters to the Diag module
DOSsetupPARAIGNORE.f90: Like DiagIGNORE.f90, should be deleted.
Inputs.f90: Inputs/parameters module, this is not actually used. Updated version found in actual compile locations.
PreAnalysis.f90: Code used to analyze the systems before the cluster ensemble (FSC, DLB, etc)
Tools.f90: Subroutines that are just useful for coding and have nothing to do with physics (e.g. filenaming, data binning).

This version no longer “removes clusters” but instead treats them in order from left to right. There is no effective hopping so there is no point in setting them up in a “nice” way. Old version saved in AHMCopies.