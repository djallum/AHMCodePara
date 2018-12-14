module DOSsetup
  USE Inputs
  USE Tools
  USE PreAnalysis
  USE Diag, only: DiagCluster
  USE mpi
  implicit none
  SAVE


   !This version is redesigned to allow for MPI
  !created: 180728; last edited: 180728
  !So far: added my_* variables to Full_DOS and system_DOS and getdos routines.
  ! Added an mpi_reduce statement, init and finalize.
  ! Current version: systemn = # of systems per process, rather than total number of systems.

contains

  subroutine GetPotential( SitePotential, weakL, weakR, Cl_Size, SitesRemoved )
    implicit none
    real(dp), dimension(:), intent(in) :: SitePotential
    integer, intent(in) :: Cl_Size, weakL, weakR, SitesRemoved
    real(dp), dimension(Cl_Size) :: Sites
    integer EndSite
    if ( weakL .gt. weakR) then
       EndSite = (dim - sitesremoved) - weakL
       Sites(1:EndSite) = SitePotential((weakL+1):(dim - sitesremoved))
       Sites((EndSite+1):Cl_Size) = SitePotential(1:weakR)
    else if ( weakL .eq. weakR) then
       Sites = SitePotential
    else
       Sites = SitePotential((weakL+1):weakR)
    end if


    CALL Bin_Data(HistoData = Potential, Data = Sites, Max = POT_EMax, Min = POT_EMin)

    

  end subroutine GetPotential

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Full_DoS routine. Inputs: Site Potentials, Hopping Amplitudes, Bond Strengths, Weak bonds array, Density of states parameters
  ! Outputs: Density of states for the system, number of contributions dropped from the density of states, Number of sites not
  ! renormalized
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine System_DoS( my_DOS, my_droppedDOS, my_SitesMissed, SitePotential, WeakBonds )
    implicit none

                !---------------------------Inputs-------------------------------------------
    real(dp), dimension(dim),intent(in) :: SitePotential             
    integer, dimension(:), intent(in) :: WeakBonds
    
                !---------------------------Outputs------------------------------------------                                                    
    real(dp), dimension(bins,DOS_MaxCluster), intent(inout) :: my_DOS
    real(dp), dimension(DOS_MaxCluster), intent(inout) :: my_droppedDOS
    integer, intent(inout) :: my_SitesMissed
    

                !---------------------------Programming Variables----------------------------
    integer numCluster, loop1
    integer nextWeak, weakL, weakR                                        ! The label for the weak bond on the "left" of the current cluster
    integer ClusterSize                                  ! Size of the current cluster

    real(dp), allocatable, dimension(:) :: Energy, Weight
    real(dp), dimension(:), allocatable :: Sites

    logical :: Periodic
                !---------------------------Initializing-------------------------------------
    
    weakL = 0
    weakR = 0

    Periodic = .false.

    
    do loop1 = 1,Size(WeakBonds)
       weakL = WeakBonds(Loop1)                                  ! Currently cluster has weakL as its left neighbour
       
       nextWeak = Loop1 + 1                                      ! Next bond label in WeakBonds creates a cluster with weakL
       if ( (Loop1 + 1) .gt. size(WeakBonds) ) then              ! If weakL is the last weakbond in the system then search the beginning of the WeakBonds array
          nextWeak = 1                                           
          weakR = WeakBonds(nextWeak)                            ! weakR is the first weak bond in the system
          ClusterSize = weakR - ( weakL - dim ) ! Have to calculate the cluster size in this way due to boundary conditions
       else                                                      ! For all other cases, the weakR label is just the next element in the weakbonds array
          weakR = WeakBonds(nextWeak)
          ClusterSize = weakR - weakL
       end if
       
       if (ClusterSize .gt. ClusterMax) then
          my_SitesMissed = my_SitesMissed + ClusterSize
          CYCLE
       end if
       
       

       allocate(Sites(ClusterSize))
       allocate(Energy(2*ClusterSize*(4**ClusterSize)))
       allocate(Weight(2*ClusterSize*(4**ClusterSize)))
	
       Energy = 0.0_dp
       Weight = 0.0_dp
       
       CALL DefineCluster( SitePotential, Sites, weakL, weakR, ClusterSize )

       CALL DiagCluster(ClusterSize,4**ClusterSize,Sites,Energy,Weight,ChemPot,uSite,hop,Periodic)
       
       CALL BinDOS( my_Dos, Energy, Weight, my_droppedDOS, ClusterSize )

       deallocate(Sites, Energy, Weight)
          
    end do

    if ( WeakBonds(1) .eq. 0 ) then
       ClusterSize = dim
       Periodic = .true.
       weakL = 0
       weakR = 1
       
       allocate(Sites(ClusterSize))
       allocate(Energy(2*ClusterSize*(4**ClusterSize)))
       allocate(Weight(2*ClusterSize*(4**ClusterSize)))
       
       Energy = 0.0_dp
       Weight = 0.0_dp
       
       CALL DefineCluster( SitePotential, Sites, weakL, weakR, ClusterSize )
       
       CALL DiagCluster(ClusterSize,4**ClusterSize,Sites,Energy,Weight,ChemPot,uSite,hop,Periodic)
       
       CALL BinDOS( my_Dos, Energy, Weight, my_droppedDOS, ClusterSize )
       
       deallocate(Sites, Energy, Weight)
    end if
       
        
    
  end subroutine System_DoS

  subroutine Full_Dos()
    implicit none
    real(dp), dimension(dim) :: SitePotential, Hopping, Bonds
    integer, dimension(:), allocatable :: WeakBonds
    
    integer SitesIgnored                                ! Keeps track of all the sites in clusters larger than this code can handle
    
    
    !----------------------------Coding Tools-------------------------------------
    integer i, j, k                                     ! Loop Integer
    real :: start, finish, TIME
    integer :: my_id, ierr, num_procs
    real(dp) :: my_DOS(bins,DOS_MaxCluster), my_droppedDOS(DOS_MaxCluster)
    integer :: my_SitesMissed
 
    !DOS = 0.0
    my_DOS = 0.d0
    my_droppedDOS = 0.d0
    my_SitesMissed = 0
    Potential = 0.d0
    PrunedBonds = 0
    




    CALL CPU_TIME(start)

    CALL mpi_init(ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)
   

    CALL init_random_seed(my_id)
    if (my_id .eq. 0) CALL CorrectInputs( )
    
    do i = 1,systemn
       if ( mod(real(i),0.1*systemn) == 0.0 ) then
          print*, int(i/real(systemn)*100), "%"
       end if
       SitePotential = 0.0_dp
       Hopping = 0.0_dp
       Bonds = 0.0_dp
       
       !---------------------------Create System------------------------------------
       call create_AHM( SitePotential, Bonds, Hopping)
       call CalcWeakBonds( Bonds, WeakBonds )
       if ( (WeakBonds(1) .eq. 0) .and. (dim .gt. ClusterMax) ) then
          !print*, "All Bonds Strong"
          CYCLE
       end if
       CALL system_Dos( my_DOS, my_droppedDOS, My_SitesMissed, SitePotential, WeakBonds )
    end do
    CALL mpi_Barrier(MPI_COMM_WORLD, ierr)
    
    allocate(DOS(bins,DOS_MaxCluster), DroppedDos(DOS_MaxCluster))
    DOS = 0.d0
    DroppedDos = 0.d0
    SitesMissed = 0
             

    do i=1,DOS_MaxCluster
       CALL mpi_reduce(my_DOS(:,i), DOS(:,i), bins, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr) 
    end do
    CALL mpi_reduce(my_droppedDOS(:), DroppedDos(:), DOS_MaxCluster, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(my_SitesMissed, SitesMissed, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)

    if ( my_id .ne. 0 ) deallocate(DOS, DroppedDos)
       
    if (my_id .eq. 0) then
      
       CALL CPU_TIME(finish)
       TIME = finish - start

       do i=1,ClusterMax
          !print*, "Cluster: ", i
          do j = 1,bins
             !print*, DOS(j,i)
             if (DOS(j,i) .ne. DOS(j,i)) then
                print*, "DOS NaN", j
             end if
             
          end do
       end do

       
       If ( CalcDos ) then
          
          
          CALL PrintDOS( TIME, num_procs )
          
       end If
       If ( CalcPot ) then
          CALL OpenFile(200, "POT", "Distribution of Site Potentials in counted clusters", "Site Potentials", &
               "Height of distribution", num_procs )
          write(200,*) "#Maximum cluster Size included: ", ClusterMax
          write(200,*) "#Time (s) = ", TIME
          CALL PrintData(200, '(g12.5,g12.5)', POT_EMin, POT_EMax, bins, Potential )
          Close(200)
       end If
    end if

    CALL mpi_finalize(ierr)
  end subroutine Full_Dos

  subroutine BinDOS( myDOS, Energy, Weight, myDropped, ClusterSize )
    implicit none
    integer :: i
    
    ! Inputs
    integer, intent(in) :: ClusterSize
    real(dp), dimension(2*ClusterSize*(4**ClusterSize)), intent(in) :: Energy, Weight

    ! Outputs
    real(dp), dimension(bins,DOS_MaxCluster), intent(inout) :: myDOS
    real(dp), dimension(DOS_MaxCluster), intent(inout) :: myDropped
    

    if ( .not. GradientDOS ) then
       
       CALL Bin_Data(myDOS(:,ClusterSize), Energy, Weight, myDropped(ClusterSize), DOS_EMax, DOS_EMin)
       
    else if (DOS_MaxCluster .eq. ClusterMax) then
       
       do i=ClusterSize,DOS_MaxCluster
          if (i .eq. ClusterMax) then
             CALL Bin_Data(myDOS(:,i), Energy, Weight, myDropped(i), DOS_EMax, DOS_EMin)
          else 
             CALL Bin_Data(HistoData = myDOS(:,i), Data = Energy, Weights = Weight &
                  , Max = DOS_EMax, Min = DOS_EMin )
          end if
       end do
       
    else if ( DOS_MaxCluster .eq. 1 ) then
       
       CALL Bin_Data(myDOS(:,1), Energy, Weight, myDropped(1), DOS_EMax, DOS_EMin)

    end if
    
    
  end subroutine BinDOS

  subroutine PrintDOS( TIME, num_procs )
    implicit none
    integer :: i,j
    real :: totalDOSweight
    
    ! Inputs
    real, intent(in) :: TIME
    integer, intent(in) :: num_procs
    if ( .not. GradientDOS ) then
       totalDOSweight = 0.0
       do i=1,ClusterMax
          totalDOSweight = totalDOSweight + Sum(DOS(:,i))+DroppedDos(i)
       end do
       do i=1,ClusterMax
          CALL OpenFile(100+i, "DOS_"//trim(str(i))//"Site_", "Density of States", "Energy", "Density of States", num_procs )
          write(100+i,400) "#Cluster Size included: ", i
          write(100+i,500) "#Fraction of sites missed: ", SitesMissed/real(dim*systemn*num_procs)
          write(100+i,500) "#Time (s) = ", TIME
          write(100+i,600) "#This is not a gradient DOS, part", i, "of", ClusterMax
          CALL PrintData(100+i, '(F15.8,F15.8)', DOS_EMin, DOS_EMax, bins, DOS(:,i), DroppedDos(i) )
          close(100+i)
       end do
    else if ( DOS_MaxCluster .eq. 1 ) then
       CALL OpenFile(100, "DOS", "Density of States", "Energy", "Density of States", num_procs )
       write(100,400) "#Maximum cluster Size included: ", ClusterMax
       write(100,500) "#Fraction of sites missed: ", SitesMissed/real(dim*systemn*num_procs)
       write(100,500) "#Time (s) = ", TIME
       write(100,'(A)') "#This is not a gradient DOS"
       CALL PrintData(100, '(F15.8,F15.8)', DOS_EMin, DOS_EMax, bins, DOS(:,1), DroppedDos(1))
       Close(100)
    else if ( DOS_MaxCluster .eq. ClusterMax ) then
              do i=1,ClusterMax
          CALL OpenFile(100+i, "DOSInc"//trim(str(i))//"_", "Density of States", "Energy", "Density of States", num_procs )
          write(100+i,400) "#Maximum cluster Size included: ", i
          print*, SitesMissed, real(dim*systemn*num_procs)
          write(100+i,500) "#Fraction of sites missed: ", SitesMissed/real(dim*systemn*num_procs)
          write(100+i,500) "#Time (s) = ", TIME
          write(100+i,600) "#This is a gradient DOS, part", i, "of", ClusterMax

          DOS(:,i) = DOS(:,i)/(Sum(DOS(:,DOS_MaxCluster)) + DroppedDos(DOS_MaxCluster))
          
          CALL PrintData(100+i, '(F15.8,F15.8)', DOS_EMin, DOS_EMax, bins, DOS(:,i))
          Close(100+i)
       end do
    end if

400 format(A,I2) !Text and 2-digit integer
500 format(A,F10.5) !Text and 9 digit real with a decimal point, with up to 5 decimal places
600 format(A,I2,A,I2) !Text, 2 digit integer, text, 2 digit integer
    
    
  end subroutine PrintDOS

  end module DOSsetup



  
