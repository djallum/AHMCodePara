! MODULE INPUTS
!
! This module contains all of the input variables. In principle, only this file should be adjusted for new runs, also occasionally the calling program as well. The other pieces of code should only be opened to add new features or to fix bugs.
!
! ---------------------------
! Inputs :
! ---------------------------
! PHYSICAL PARAMETERS:
!      systemn : Number of systems in the ensemble
!      dim : Number of sites in each system
!      DELTA : Disorder strength
!      hop : hopping amplitude
!      uSite: Interaction Strength
!      ChemPot : Chemical potential

! CODING PARAMETERS:
!      bond_cutoff : Defines the bond strength seperating weak bonds and strong bonds
!      prune_cutoff : Defines the minimum difference in energy between ChemPot and the average energy of n.n. orbital pairs in order to ignore a bond
!      bins : number of elements in a number of arrays containing output data (aka number of bins in histograms)
!      ClusterMax : The maximum cluster size to be included in calculating the 'DOS' or 'Potential'
!      DOS_MaxCluster : Length of the first dimension of the DOS array.
!         If DOS_MaxCluster = ClusterMax. It is the length of the first dimension of DOS. ELSE: The first dimension of DOS has a length of 1.
!      DOS : Density of States Array. Histogram of the DOS.
!         DOS(i,:) DOS including one to i-site clusters. i ranges from 1 to DOS_MaxCluster
!      DOS_EMax : Maximum energy of a DOS contribution allowed
!      DOS_EMin : Minimum energy of a DOS contribution allowed
!      DroppedDos : Combined weight of DOS contributions outside of the range [DOS_EMin,DOS_EMax]
!      PrunedBonds : Number of n.n. pairs in which at least one orbital pair was ignored
!      StrongestBondsPruned : Number of pruned bonds
!      Potential: Distribution of site potentials
!      POT_EMax : Maximum site potential kept (Value is obvious but required for Bin_Data subroutine in Tools.f90)
!      POT_EMin : Minimum site potential kept
!      Bond_EMax : Maximum ln(Bond Strength) kept
!      Bond_EMin : Minimum ln(Bond Strength)
!      Diff_EMax : Maximum difference in energy 
!      Diff_EMin : Minimum difference in energy b/w n.n. sites

! Run This Part? Logical parameters to avoid running certain parts of the code. You can get DOS without disitribution of site potentials.
!      CalcFracSites : Calculate the fraction of sites per cluster size
!      CalcBondStrength : Calculate the bond strength disitribution
!      CalcEnergyDiff : Calculate the energy difference distibution
!      CalcDos : Calculate the Density of States
!      CalcPot : Calculate the distribution of site potentials
! See end of file for code overview


module inputs
  implicit none
  
  ! PHYSICAL PARAMETERS
  integer, parameter :: systemn = 2500000                      
  integer, parameter :: dim = 4                          
  real, parameter :: DELTA = 20                             
  real, parameter  :: hop = -1                              
  real, parameter :: uSite = 8                             
  real, parameter :: ChemPot = uSite/2                      
  
  ! CODING PARAMETERS
  real, parameter :: bond_cutoff = 0.5
  real, parameter :: prune_cutoff = 3
  integer, parameter :: bins = 1000
  integer, parameter :: ClusterMax = 7
  integer, parameter :: DOS_MaxCluster = ClusterMax
  real, allocatable :: DOS(:,:), DroppedDos(:) !DOS(DOS_MaxCluster,bins) ;; DroppedDos(DOS_MaxCluster)
  real, parameter :: DOS_EMax = 5, DOS_EMin = -DOS_EMax
  integer :: PrunedBonds 
  integer :: StrongestBondsPruned
  integer :: SitesMissed 
  real Potential(bins)
  real, parameter :: POT_EMax = Delta/2, POT_EMin = -Pot_EMax
  real, parameter :: Bond_EMax = 20, Bond_EMin = -20
  real, parameter :: Diff_EMax = 3*DELTA, Diff_EMin = -3*DELTA
  !real :: POT_EMax , POT_EMin, Bond_EMax, Bond_EMin, Diff_EMax, Diff_EMin
  ! RUN THIS PART?
  logical, private, parameter :: ff = .false., tt = .true.
  logical, parameter :: CalcFracSites = ff ! 
  logical, parameter :: CalcBondStrength = ff
  logical, parameter :: CalcEnergyDiff = ff
  logical, parameter :: CalcDos = tt
  logical, parameter :: CalcPot = ff
  logical, parameter :: GradientDOS = ff

  real :: times


  

end module inputs



! This file, along with the other .f90 files in the directory, are used to calculate important quantities related to a one-dimensional chain of atoms
! Quantities:
!    1. Fraction of sites per cluster size
!    2. Distribution of bond strengths
!    3. Energy difference between nearest neighbours distribution
!    4. Density of States (Most important)
!    5. Distribution of site potentials in clusters of X sites

! .f90 files:
!    1. Inputs.f90 (This file. Lists all system parameters and is accessed by every module to give these parameters artificial global scope)
!    2. Tools.f90 (Any routine that is strictly a tool for other modules to use. E.g. a subroutine that bins data or prints to file)
!    3. PreAnalysis.f90 (Subroutines related to creating each system and doing any pre-"RG" analysis. Quotes b/c not really RG)
!    4. DOSsetupPARA.f90 (All routines related to the cluster ensemble or "RG". Suffix "PARA" indicates that this version is parallelized)
!    5. Diag.f90 (Routines to set up diagonalization of clusters. Patrick's code that has been adapted for my purposes)

! Code was originally set up in one file but I have split it into the 5 modules. I keep the modules in this directory but the program file
! ,and output files, are found in a separate directory. If this directory is /AHMRG/CodePARA/ then one can change to the program directory
! with `cd ../Stage3/CodePARA/`. Stage3 is in reference to stages in my research but is currently (180821) the main folder for DOS calculations.

! In the program directory is a makefile which details how to build the executable. It is as simple as using the `make` command. Or the safer
! version (which guarantees the recompiling of every module) `make clean ; make`. Sometimes not everything is recompiled properly.

! Since this code is parallelized with OpenMPI we require the `mpirun` command to run the code and specify the nodes. Run the code with the command
! `mpirun -n 2 ./V3.e` `-n 2` specifies that I will run this on two processors. Otter/Salmon and Orca are all dual core and therefore 2 is the max.
! HOWEVER, oversubscribing cores with OpenMPI is a thing that exists. Feel free to google it if you really want more processes than processors.
! NOTE: `mpirun ./V3.e` runs the program on the number of processors that MPI detects. It is equivalent to the first run statement but is "safer" to
! use the former.

!**********************************

! This file:

! A variable declared in a module is available to any subroutine in the module. It is also accessible to any module or program that references this module.
! This allows my parameters to be accessed anywhere and avoids unnecessarily large subroutine call statements.
