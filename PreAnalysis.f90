Module PreAnalysis
  USE Inputs
  USE Tools
  Implicit none
  save
  Logical :: FirstRun = .true.

contains
  subroutine CorrectInputs( )
    implicit none
    character(len=1) :: Ans

    if ( .not. FirstRun ) return

    FirstRun = .false.
    
    print*, "Disorder Strength: ", DELTA, "Bond cutoff: ", bond_cutoff, "Interaction Strength", uSite
    print*, "System Size: ", dim, "Number of systems: ", systemn, "Hopping Potential: ", Hop
    print*, "Prune Cutoff: ", prune_cutoff; print*, "Drop cutoff: ", drop_cutoff
    

    if ( CalcDos ) then
       print*, "Density of States", CalcDos
       print*, "Maximum cluster size: ", ClusterMax, "Chemical Potential: ", ChemPot
    else
       print*, "Density of States", CalcDos
    end if

    if ( (DOS_MaxCluster .eq. 1) .and. ( .not. GradientDOS ) ) then
       print*, "Incompatible values of DOS_MaxCluster and GradientDOS"
       STOP
    end if
    
    if ( DOS_MaxCluster .eq. 1 ) then
       print*, "DOS per cluster size disabled"
    else if (DOS_MaxCluster .eq. ClusterMax) then
       print*, "DOS per cluster size enabled"
       if ( GradientDOS ) then
          print*, "DOS is cumulative"
       else
          print*, "DOS is not cumulative"
       end if
    else
       print*, "This value of DOS_MaxCluster is not handled: ", DOS_MaxCluster
       STOP
    end if


    
    

    if ( CalcPot ) then
       print*, "Distribution of Site Potentials", CalcPot
       print*, "Maximum Cluster size: ",  ClusterMax
    else
       print*, "Distribution of Site Potentials", CalcPot
    end if

    print*, "Fraction of Sites per Cluster Size", CalcFracSites

    print*, "Distribution of log(bonds): ", CalcBondStrength

    print*, "Distribution of \Delta\epsilon: ", CalcEnergyDiff

 
    return    
    do while ( .true. )
       print*, "Is this correct? (Y/N)"
       read(*,*) Ans

       SELECT CASE (Ans)
       CASE ('y', 'Y')
          return
       CASE ('n', 'N')
          print*, "Inputs incorrect"
          STOP
       CASE DEFAULT
          print*, "Invalid input, try again. Enter y/Y for yes, n/N for no"
       end SELECT
    end do
  end subroutine CorrectInputs

!
! This subroutine is from Eamonn's code and it defines a start point for the seed of the random number generator based on the time that the code is run
!
  subroutine init_random_seed(my_id)
    
    implicit none    
    integer, intent(in) :: my_id
    integer,dimension(:), allocatable :: seed !seed for random number generator
    integer  e,clock,i5  !clock is the system time in seconds, e is the size of seed of the random number generator
    call random_seed (size = e) !find out how big the seed is to be. This size of seed is e
    allocate (seed(e)) !allocate space for seed
    call system_clock (count= clock) !obtain the system time, in seconds
    seed = clock + 37*(/(i5 - 1,i5 = 1, e)/) !set the seed using the system time
    call random_seed (put=(seed+my_id)) !initialize the random number generator
    deallocate (seed) !deallocate the seed array
    return
  end subroutine init_random_seed
  
  !
  ! This subroutine creates the system, generates the site potentials, calculates bond strengths, stores hopping in array, might make it find all the clusters too
  ! Test002 shows that the site potentials are being selected correctly and the strongest bond has been chosen.
  !
  
  subroutine create_AHM( SitePotential, Hopping, Bonds, EDiff )
    implicit none
    
    !---------------------------Outputs------------------------------------------
    real, dimension(dim), intent(out) :: SitePotential, Hopping, Bonds
    real, dimension(dim), optional, intent(out) :: EDiff
    
    !---------------------------Coding Tools-------------------------------------
    real RandomNumber                                    ! random number for site potential calculation
    integer loop1                                        ! loop integer 
    integer IndexOfNeighbour                             ! this stores the site number of the right neighbour site in a bond. Used for
    ! boundary conditions
    
    !---------------------------Variable Allocations-----------------------------
    
    !---------------------------Initializing-------------------------------------
    RandomNumber = 0
    Hopping(1:dim) = 0 
    SitePotential(1:dim) = 0
    Bonds(1:dim) = 0

                !---------------------------Calculate Site Potentials/Hopping Amplitudes-----
    !This loop uses the RNG to determine the site potentials and initializes all the hopping potentials as the same.
    do loop1 = 1,dim
       call RANDOM_NUMBER(RandomNumber)                  ! Sets 'RandomNumber' to a random number between [0,1]
       !if ( mod(loop1,2) .eq. 0 ) then
          Hopping(loop1) = hop
       !end if
       SitePotential(loop1) = (0.5 - RandomNumber)*DELTA ! Shifts 'RandomNumber' from interval [0,1] to [-DELTA,DELTA]
    end do
    !print*, Hopping
                !---------------------------Calculate the 'Bonds' array----------------------
    do loop1 = 1,dim 
       IndexOfNeighbour = loop1+1                        ! Labels the neighbour site for the purposes of the boundary conditions

       if (loop1 .eq. dim) then                          ! If loop1 == number of sites then the right neighbour site is site 1
          IndexOfNeighbour = 1
       end if
       !---------------------------Calls the subroutine-----------------------------
       If ( Present(EDiff) ) then
          CALL Strongest_Bond( Bonds, SitePotential, Hopping, loop1, IndexOfNeighbour, EDiff )
       else
          CALL Strongest_Bond( Bonds, SitePotential, Hopping, loop1, IndexOfNeighbour )
       end If
       
       

    end do
   
    return
  end subroutine create_AHM

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  ! Strongest_Bond routine. Inputs: Site potentials, Hopping Amplitude, adjacent site labels. Output: Bonds array
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine Strongest_Bond( Bonds, SitePotential, Hopping, LeftNbr, RightNbr, EDiff )
    implicit none
                !---------------------------Output-------------------------------------------
    real, dimension(:), intent(out) :: Bonds
    real, dimension(:), optional, intent(out) :: EDiff

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), intent(in) :: SitePotential
    real, dimension(:), intent(in) :: Hopping
    integer, intent(in) :: LeftNbr                       ! Site label of the left neighbour site
    integer, intent(in) :: RightNbr                      ! Site label of the right neighbour site

                !---------------------------Programming Tools--------------------------------
    real ABond, HBond1, HBond2                           ! ABond: Anderson bond; HBond1: Hubbard bond with right neighbour UHO
                                                         ! HBond2: Hubbard bond with the left neighbour UHO
    real AAvg, HAvg, UAvg                                ! Average strength of LHOs (AAvg), of UHOs (UAvg) and of any LHO-UHO pair (HAvg)
    logical Pruned

    !---------------------------Calculate the three possible bond strengths------
    ABond = abs(Hopping(LeftNbr)/(SitePotential(LeftNbr) - SitePotential(RightNbr)))
    HBond1 = abs(Hopping(LeftNbr)/(SitePotential(LeftNbr) - (SitePotential(RightNbr) + uSite)))
    HBond2 = abs(Hopping(LeftNbr)/((SitePotential(LeftNbr) + uSite) - SitePotential(RightNbr)))

    !---------------------------Calculate the orbital energy average for all three bonds----------
    AAvg = (SitePotential(LeftNbr) + SitePotential(RightNbr))/2
    HAvg = (SitePotential(LeftNbr) + (SitePotential(RightNbr) + uSite))/2
    UAvg = ((SitePotential(LeftNbr) + uSite) + (SitePotential(RightNbr) + uSite))/2

    if ( Size(SitePotential) .eq. dim ) then
       CALL PruningStats( AAvg, HAvg, UAvg, ABond, HBond1, HBond2 )
    end if
    

                !---------------------------Determine which is strongest---------------------
    if ( (ABond .ge. HBond1) .and. &
         (ABond .ge. HBond2) ) then                      ! If the noninteracting bond is strongest, then store that bond
       
       Bonds(LeftNbr) = ABond
       
       If ( Present(EDiff) ) &
            EDiff(LeftNbr) = SitePotential(LeftNbr) - SitePotential(RightNbr)
       
    else if ( (HBond1 .ge. ABond) .and. &                ! If the bond between the LHO of the left site and the UHO of the right
         (HBond1 .ge. HBond2) ) then                     ! site is the strongest, then store that bond
       
       Bonds(LeftNbr) = HBond1
       
       If ( Present(EDiff) ) &
            EDiff(LeftNbr) = SitePotential(LeftNbr) - (SitePotential(RightNbr) + uSite)
       
    else if ( (HBond2 .ge. ABond) .and. &                ! If the bond between the UHO of the left site and the LHO of the right
         (HBond2 .ge. HBond1) ) then                     ! site is the strongest, then store that bond
       
       Bonds(LeftNbr) = HBond2
       
       If ( Present(EDiff) ) &
            EDiff(LeftNbr) = (SitePotential(LeftNbr) + uSite) - SitePotential(RightNbr)
       
    else if ( (IsNAN(ABond)) .or. (IsNAN(HBond1)) .or. (IsNAN(HBond2)) ) then
       Bonds(LeftNbr) = huge(0.0)
       
    else if ( (.not. ( Size(SitePotential) == 1 )) ) then
       print*, 'You belong in a museum, Error:001', SitePotential(LeftNbr), SitePotential(RightNbr)   ! Prints an error to screen in case there is an option not covered
                                                         ! Doesn't print when there is only a single site left in the system (bond
                                                         ! strength ends up being infinite)
       !print*, ABond, HBond1, HBond2
       Bonds(LeftNbr) = 0
    end if


    
       ! Note on above: only use 'greater than or equal to' when comparing bonds. This is because if they are equal it doesn't matter
       ! which is chosen
  end subroutine Strongest_Bond

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! PruningStats routine. Inputs: Bonds array, SitesRemoved. Output: WeakBonds array
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine PruningStats( AAvg, HAvg, UAvg, ABond, HBond1, HBond2 )
    implicit none
    real, intent(inout) :: ABond, HBond1, HBond2
    real, intent(in) :: AAvg, HAvg, UAvg
    real :: Bonds(3), Temp_Bonds(3)
    logical :: Pruned(3), IsOnePruned
    integer :: Strongest_index
    integer :: i !Looping
    real :: NextStrong
    
    Bonds(1) = ABond; Bonds(2) = HBond1; Bonds(3) = HBond2
    Pruned(1) = .false.; Pruned(2) = .false.; Pruned(3) = .false.
    IsOnePruned = .false.


    If ( (ABS(AAvg - ChemPot) .gt. prune_cutoff) .and. (ABS(UAvg - ChemPot) .gt. prune_cutoff) ) then
       IsOnePruned = .true.
       Pruned(1) = .true.
       ABond = 0
    end If

    if ( abs( HAvg - ChemPot ) .gt. prune_cutoff ) then
       IsOnePruned = .true.
       Pruned(2) = .true.; Pruned(3) = .true.
       HBond1 = 0;HBond2 = 0
    end if

    Strongest_index = MAXLOC(Bonds, dim = 1)

    !The following is an algorithm that determines whether or not the pruning of a bond is "useful"
    !A useful pruning is defined as a pruning in which the strongest bond is pruned and replaced by a bond below bond_cutoff
    !A useful pruning transforms a strong bond between a n.n. pair of sites into a weak bond
    !Only useful prunings increment "StrongestBondsPruned"
    if (Bonds(Strongest_index) .ge. bond_cutoff) then
       
       if ( Pruned(Strongest_index) ) then
          
          Temp_Bonds = Bonds
          
          Temp_Bonds(Strongest_index) = 0
          Strongest_index = MAXLOC(Temp_Bonds, dim = 1)

          if ( Pruned(Strongest_index) ) then
             
             Temp_Bonds(Strongest_index) = 0
             Strongest_index = MAXLOC(Temp_Bonds, dim = 1)

             if ( Pruned(Strongest_index) ) then
                StrongestBondsPruned = StrongestBondsPruned + 1
             else

                if ( Bonds(Strongest_index) .lt. bond_cutoff ) then
                   StrongestBondsPruned = StrongestBondsPruned + 1
                end if
                
             end if
             
          else
             
             if ( Bonds(Strongest_index) .lt. bond_cutoff ) then
                StrongestBondsPruned = StrongestBondsPruned + 1
             end if
             
          end if
          
       end if
       
    end if

    if ( IsOnePruned ) then
       PrunedBonds = PrunedBonds + 1
    end if
    

  end subroutine PruningStats
   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! CalcWeakBonds routine. Inputs: Bonds array, SitesRemoved. Output: WeakBonds array
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine CalcWeakBonds( Bonds, WeakBonds, SitesRemoved )
    implicit none
                !---------------------------Output-------------------------------------------
    integer, dimension(:), allocatable, intent(out) :: WeakBonds

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), intent(in) :: Bonds
    integer, optional, intent(in) :: SitesRemoved

                !---------------------------Programming Variables----------------------------
    integer NumWeakBonds                                 ! The number of weak bonds in the system at a given point in the RG
    integer loop1

                !---------------------------Initializing-------------------------------------
    NumWeakBonds = 0

    !---------------------------Count Weak Bonds---------------------------------
    if (Present(SitesRemoved)) then

       do loop1 = 1, ( dim - SitesRemoved )
          if ( Bonds(loop1) .lt. bond_cutoff ) then         ! Check if each bond is below the cutoff. If yes, increment NumWeakBonds
             NumWeakBonds = NumWeakBonds + 1
          end if
       end do       
       if ( NumWeakBonds .gt. 0 ) then                      ! If there are any weak bonds then calculate the array
          allocate(WeakBonds(NumWeakBonds))
          NumWeakBonds = 0
          do loop1 = 1,( dim - SitesRemoved ) 
             if ( Bonds(loop1) .lt. bond_cutoff ) then
                NumWeakBonds = NumWeakBonds + 1
                WeakBonds(NumWeakBonds) = loop1
                if ( WeakBonds(NumWeakBonds) .le. 0 ) then 
                   print*, 'You belong in a museum, Error:002', loop1
                   stop
                end if
                
             end if
          end do
       else if ( NumWeakBonds .eq. 0 ) then                 ! If there are no weak bonds then set WeakBonds as follows. This cannot happen
          ! under normal circumstances. Used to tell the code that there are no weakbond
          allocate(WeakBonds(1))
          WeakBonds(1) = 0
       end if


       
    else


       
       do loop1 = 1,dim
          if ( Bonds(loop1) .lt. bond_cutoff ) then         ! Check if each bond is below the cutoff. If yes, increment NumWeakBonds
             NumWeakBonds = NumWeakBonds + 1
          end if
       end do
       if ( NumWeakBonds .gt. 0 ) then                      ! If there are any weak bonds then calculate the array
          allocate(WeakBonds(NumWeakBonds))
          NumWeakBonds = 0

          do loop1 = 1,dim
             if ( Bonds(loop1) .lt. bond_cutoff ) then
                NumWeakBonds = NumWeakBonds + 1
                WeakBonds(NumWeakBonds) = loop1
                if ( WeakBonds(NumWeakBonds) .le. 0 ) then 
                   print*, 'You belong in a museum, Error:002', loop1
                   stop
                end if
                
             end if
          end do
       else if ( NumWeakBonds .eq. 0 ) then                 ! If there are no weak bonds then set WeakBonds as follows. This cannot happen
          ! under normal circumstances. Used to tell the code that there are no weakbond
          allocate(WeakBonds(1))
          WeakBonds(1) = 0
       end if


       
    end if
    
  end subroutine CalcWeakBonds

  subroutine AnalyzeClusters( )
    implicit none
    integer loop1
    real, dimension(dim) :: SitePotential, Hopping, Bonds
    real, dimension(dim) :: ClusterCount, EDiff
    integer, dimension(:), allocatable :: WeakBonds
    real, dimension(bins) :: HistoBonds, HistoEDiff
    

    ClusterCount = 0.0
    SitePotential = 0.0
    Hopping = 0.0
    Bonds = 0.0
    HistoBonds = 0.0
    HistoEDiff = 0.0
    
    do loop1=1,systemn

       CALL create_AHM( SitePotential, Hopping, Bonds, EDiff )
       CALL CalcWeakBonds( Bonds, WeakBonds )
       if ( CalcFracSites ) &
            CALL FracSites( WeakBonds, ClusterCount )
       
       if ( CalcBondStrength ) &
            CALL Bin_Data( HistoData = HistoBonds, Data = log(Bonds), Max = Bond_EMax, Min = Bond_EMin ) 

       if ( CalcEnergyDiff ) &
            CALL Bin_Data( HistoData = HistoEDiff, Data = EDiff, Max = Diff_EMax, Min = Diff_EMin )
       
    end do

    if ( CalcFracSites ) then
       CALL OpenFile( 100, "FracOfSites", "Fraction of Sites vs Cluster Size", &
            "Cluster Sizes", "Fraction of Sites" )
       CALL PrintData( 100, '(g12.5,g12.5)', 0.5, real(dim)+0.5, dim, ClusterCount )
       Close(100)
    end if
    if ( CalcBondStrength ) then
       CALL OpenFile( 200, "LogBonds", "Log(Distribtion of the log(bonds))", &
            "Right-edge of bins in log(bonds) space", "log(Distribution of log(bond strengths)" )
       CALL PrintData( 200, '(g12.5,g12.5)', Bond_EMin, Bond_EMax, bins, HistoBonds )
       Close(200)
    end if
    if ( CalcEnergyDiff ) then
       CALL OpenFile( 300, "EDIFF", "Distribution of energy differences", &
            "\Delta\epsilon", "Height of distribution at \Delta\epsilon" )
       CALL PrintData( 300, '(g12.5,g12.5)', Diff_EMin, Diff_EMax, bins, HistoEDiff )
       Close(300)
    end if


  end subroutine AnalyzeClusters
  
 
   
 subroutine FracSites(WeakBonds, ClusterCount)
    implicit none
    integer loop1, loop2!loop integer
    integer, dimension(:), intent(in) :: WeakBonds !An array where each element is 0 or 1. 1 means weak bond, 0 means strong bond
    real, intent(out) :: ClusterCount(dim) !Counts the number of clusters recorded for each cluster size
    integer ClusterSize !contains the size of the current cluster
    integer LbondLabel, RbondLabel !The number that labels the left and right weak bonds that neighbour each cluster

   

    if ( WeakBonds(1) .ne. 0 ) then
       do loop1 = 1,SIZE(WeakBonds) !cycle over the clusters in the system, number of weak bonds == number of clusters
          
          LbondLabel = WeakBonds(loop1) !Extract the label for the left bond for the current cluster
          if (loop1 .eq. SIZE(WeakBonds)) then !If the current cluster is the last cluster, then the second weak bond in this cluster is also the first weak bond in the system
             loop2 = 1
             LbondLabel = LbondLabel - dim !Also redefine the left label as an integer <= 0, the last bond on the chain (bond between last and first sites) is also the 0th bond in the cluster.
          else if (loop1 .lt. SIZE(WeakBonds)) then !If the current bond is not the last bond, then the label for the right neighbour bond is one element larger than the left neighbour bond in WeakBonds
             loop2 = loop1 + 1
          else !In case something funky happens. Can't see how but better safe than sorry
             print*, 'You belong in a museum: error 003' 
          end if
          RbondLabel = WeakBonds(loop2) !Extract the label for the right neighbour bond for the current cluster
          
          ClusterSize = abs(RbondLabel - LbondLabel) !this calculates the size of the cluster. If a cluster is between bond 3 and 4, then the cluster size is 4-3=1
          if (ClusterSize .eq. 0) then
             print*, 'Rbond', RbondLabel
             print*, 'Lbond', LbondLabel
             print*, 'Cluster Size', ClusterSize
             do loop2 = 1,5
                print*, 'WeakBonds', WeakBonds(loop2)
             end do
             print*, 'length of WeakBonds', size(WeakBonds)
             stop
          end if
          
          ClusterCount(ClusterSize) = ClusterCount(ClusterSize) + ClusterSize !Increase the number of clusters of size 'ClusterSize' by an amount equal to 'ClusterSize'
       end do
       
    else if ( WeakBonds(1) .eq. 0 ) then
       ClusterCount(dim) = ClusterCount(dim) + dim
    end if
       
    
  end subroutine FracSites

end Module PreAnalysis

  
  
  
