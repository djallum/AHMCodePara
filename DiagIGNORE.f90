! Log:
! 180320 -- Added 'Cell' type to enable calculating and storing the eigenvectors from each block simultaneously (they aren't a fixed size)
! 180401 -- Changed the neighbours routine to accurately represent the closed boundary conditions of a given cluster
! 180601 -- Added the 'CellI' type to enable calculating certain variables for all cluster sizes at the beginning of each ensemble



! In this module are the routines necessary to extract the DOS from a cluster of size 'nsite' that is sent to the subroutine 'DiagCluster' from the main program with suffix '-main.f90', the prefix is the same prefix of this file
! DiagCluster is called first and implements the other subroutines the determine the DOS of the cluster

module Diag
  use Inputs, only: drop_cutoff
  implicit none
  integer, dimension(:,:), allocatable :: fock_states       ! array that stores each FS represented in binary (see above)
  integer, dimension(:,:), allocatable :: PES_down, PES_up   ! lookup tables for PES of many-body groundstate (MBG) transformations (ex. c|Psi0>)
  integer, dimension(:,:), allocatable :: IPES_down, IPES_up ! lookup tables for PES of MBG transformations (ex. c^{dagger}|Psi0>)
  integer, dimension(:,:), allocatable :: phase_PES_down, phase_PES_up    ! to get anticommutation sign right
  integer, dimension(:,:), allocatable :: phase_IPES_down, phase_IPES_up  ! to get anticommutation sign right
  
  integer, dimension(:,:), allocatable :: msize       ! msize(i,j) is size of submatrix with n_up=i and n_dn=j
  integer, dimension(:,:), allocatable :: mblock      ! mblock(i,j) is index in fock_states array of first state with n_up=i,n_dn=j
  
  TYPE Cell
     real, allocatable :: comp(:)
  end type Cell  
  
contains
  !************************************************************************************
  !************************************************************************************
  subroutine PreSetUp(ClusterSize)
    implicit none
    integer, intent(in) :: ClusterSize

    if ( Allocated(fock_states) ) then
       deallocate( fock_states,PES_down,PES_up, &
             IPES_down,IPES_up,phase_PES_down, &
             phase_PES_up,phase_IPES_down,phase_IPES_up, &
             msize, mblock )
    end if
    allocate( fock_states(2,4**ClusterSize),PES_down(ClusterSize,4**ClusterSize),PES_up(ClusterSize,4**ClusterSize), &
         IPES_down(ClusterSize,4**ClusterSize),IPES_up(ClusterSize,4**ClusterSize),phase_PES_down(ClusterSize,4**ClusterSize), &
         phase_PES_up(ClusterSize,4**ClusterSize),phase_IPES_down(ClusterSize,4**ClusterSize), &
         phase_IPES_up(ClusterSize,4**ClusterSize),msize(0:ClusterSize,0:ClusterSize), mblock(0:ClusterSize,0:ClusterSize) )
    
    call num_sites(ClusterSize)
    call transformations(ClusterSize)
    call matrix_sizes(ClusterSize)

  end subroutine PreSetUp
  
  !************************************************************************************
  !************************************************************************************
  subroutine num_sites(ClusterSize)
	!    %----------------------------------------------------------------------------------------------%
	!    | Creates and orders the fock state basis. As an example the order for 2 site system is:       |
	!    |                                                                                              |
	!    |                                                 i=                                           |
	!    |                   | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |   |
	!    |   ----------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|   |
	!    |   fock_states(1,i)| 0 | 0 | 0 | 0 | 1 | 2 | 1 | 2 | 1 | 2  | 1  | 2  | 3	 | 3  | 3  | 3  |   |
	!    |   ----------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|   |
	!    |   fock_states(2,i)| 0 | 1 | 2 | 3 | 0 | 0 | 1 | 1 | 2 | 2  | 3  | 3  | 0  | 1  | 2  | 3  |   |
	!    |                                                                                              |
	!    %----------------------------------------------------------------------------------------------%
	!
	!    * details of the order are not particularily important and specifics can be found by working through source code
	!
	!    It first looks at only on integer (up electrons) and orders them by ne in increasing order. Then uses this ordered set
	!    of integers (variable is called order_states) to create the fock_states. The order in which each sub block of states with 
	!	 same number of electrons shows how it uses the list.
	!    
	!    Example for nsites=4
	!    	H00,H01,H02,H03,H04,H10,H11,H12,H13,H14,H20.......
	!
	!    It would move up by block of integers with n_up electrons (assigning them to fock_states(1,i)) then and cycle through the entire set
	!    (assigning fock_states(2,i) while keeping n_up in same block. (See example in the box above)
    implicit none
    integer :: i, j, istate, isite ! counters for loops
    integer, intent(in) :: ClusterSize
    integer :: max_electrons        ! maximum amount of spin-up electrons (same as nsites)
    integer :: tot_states_up        ! total number of fock states with n_dn = 0
    integer :: ne                   ! number of electrons
    integer, dimension(:), allocatable :: BlockIndex, temp_BlockIndex           ! BlockIndex(i) is index of the first state with i electrons (when only looking at states with n_dn=0)
    integer, dimension(:), allocatable :: nstates_up                  ! nstates_up(i) is the number of states with i up electrons and n_dn = 0   
    integer, dimension(:), allocatable :: states_order             ! orders the up part of FS and down part of FS before combining them

    allocate(BlockIndex(0:ClusterSize), temp_BlockIndex(0:ClusterSize), nstates_up(0:ClusterSize), &
         states_order(0:2**ClusterSize))
    tot_states_up = 2**ClusterSize     
    max_electrons = ClusterSize 
    ! enter in intial variables before starting loops
    BlockIndex = 0  
    nstates_up = 0
    BlockIndex(0) = 1
    nstates_up(0) = 1
    ! contruct the BlockIndex array
    do ne = 1, max_electrons
       nstates_up(ne) = choose(ClusterSize,ne)          ! number of fock states with n_up=ne (n_dn=0). Choose function defined at bottom of module
       BlockIndex(ne) = BlockIndex(ne-1) + nstates_up(ne-1)  ! start of ne block is the ne-1 block plus amount of ne-1 states there are
    end do
    temp_BlockIndex = BlockIndex                              ! temporary array that can be modified without destroying original block array
    
       ! %-----------------------------------------------------------------------%
       ! |	Order intergers by number of ones (electrons) they have in binary | 
       ! |                                                                       | 
       ! | Example for nsite=4:  0, 1,2,4,8, 3,5,6,9,10,12, 7,11,13,14, 15)      | 
       ! |                      |0||   1   ||      2      |     3     | 4 |      | 
       ! |                                                                       | 
       ! %-----------------------------------------------------------------------% 
    do istate=0,tot_states_up-1
       ne = 0
       do isite=1,ClusterSize
          ne = ne + ibits(istate,isite-1,1)     ! count the number of one of that binary number (up electrons of the state)
       end do
       states_order(temp_BlockIndex(ne)) = istate     ! place it at the lowest index of the block with the same n_up
       temp_BlockIndex(ne) = temp_BlockIndex(ne) + 1       ! raise the lowest index of that block by 1 (don't overwrite what you just put in)
    end do
    ! Uses the ordered list (states_order) created in last loop to make the fock_states
    
    istate = 1
    fock_states(:,:) = 0
    do ne=0,max_electrons
       do j=1,tot_states_up
          do i=BlockIndex(ne), BlockIndex(ne) + nstates_up(ne) - 1
             fock_states(1,istate) = states_order(i)
             fock_states(2,istate) = states_order(j)
             istate = istate + 1
          end do
       end do
    end do
    deallocate(BlockIndex, temp_BlockIndex, nstates_up, &
         states_order)

    
  end subroutine num_sites
  !************************************************************************************
  !************************************************************************************
  subroutine transformations(ClusterSize)
    !  %---------------------------------------------------------------------------------------------------%
    !  |  Creates the lookup tables needed to calculate c_{i,sigma}|Psi0> and c_{i,sigma}^{dagger}|Psi0>   |
    !  |       - Psi0 is the many-body ground state                                                        |
    !  |       - c_{i,sigma} is removable of electron of spin sigma from site i                            |
    !  |       - c_{i,sigma}^{dagger} is addition of electron of spin sigma on site i                      |
    !  |  Shows which fock state each of the basis vectors will transfer to after the specified removable  |
    !  |  or addition of up/dn electron from specified site.                                               |
    !  |                                                                                                   |
    !  |  PES_up(j,i) is the state (as in its index within the basis) that after removing an up electron   |
    !  |  from site j of it, will become state i.                                                          |
    !  |                                                                                                   |
    !  |  It does this by adding a '1' (binary) to the site then checking what fock state it now became    |
    !  |  by compairing the new integers to the FS basis. To get anticommutation sign right it counts the  |
    !  |  amount of up and dn electrons starting at the site the electron was added and going to           |
    !  |  site=nsites. This gives the number of anti-commutations that the creation operator would have    |
    !  |  to do.                                                                                           | 
    !  |                                                                                                   |
    !  |  The program then uses the PES tables to calculate the IPES tables since they are opposites       |
    !  %---------------------------------------------------------------------------------------------------%
    implicit none
    integer, intent(in) :: ClusterSize
    integer :: i, j                           ! counters for loops
    integer :: new_state(2)                   ! new fock state integers after the transition 
    integer :: new_index                      ! index of the new state in the array fock_states
    integer :: position, isite, istate                ! both counters for loops
    integer :: ne                             ! number of electrons
    
    !------------------Make the PES_up tables-------------------------

    istate = 4**ClusterSize
    
    
    PES_up(:,:) = 0; PES_down(:,:) = 0
    
    
    IPES_up(:,:) = 0; IPES_down(:,:) = 0
    
    
    phase_PES_up(:,:) = 0; phase_PES_down(:,:) = 0
    
    
    phase_IPES_up(:,:) = 0; phase_IPES_down(:,:) = 0
    
    
    do position = 1,ClusterSize                                         ! loop over all the sites (make PES_up for PE from each site)
       do i=1,istate                                             ! loop over each state
          ne = 0
          if (ibits(fock_states(1,i),position-1,1) == 1) then    ! if there is an electron on that site it can't be the result of PE
             PES_up(position,i) = 0                             ! zero everything then because it's not possible
             phase_PES_up(position,i) = 0                       ! zero everything then because it's not possible
          else
             do isite=position,ClusterSize                                      ! loop over all sites greater then site of PE (count number of anti-commutations)
                ne = ne + ibits(fock_states(1,i),isite-1,1)    ! count the number of up electrons it will have to commute with to be removed
             end do
             do isite=position,ClusterSize
                ne = ne + ibits(fock_states(2,i),isite-1,1)    ! count the number of dn electrons it will have to commute with to be removed
             end do
             if (MOD(ne,2) == 0) then 
                phase_PES_up(position,i) = 1                   ! if it had to do even number of anti-commutations it is positive
             else
                phase_PES_up(position,i) = -1                  ! if it had to do odd number of anti-commutations it is positive
             end if
             new_state(1) = ibset(fock_states(1,i),position-1)  ! add up electron to that site 
             new_state(2) = fock_states(2,i)                    ! the down portion remains the same
             do j=1,istate
                if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then 
                   new_index = j                              ! find the index of the new state by compairing it to the entire FS basis
                end if
             end do
             PES_up(position,i) = new_index                     ! record the state that will when PE will become state i
          end if
       end do
    end do
    !------------------Make the PES_dn tables-------------------------
    do position = 1,ClusterSize
       do i=1,istate
          ne = 0
          if (ibits(fock_states(2,i),position-1,1) == 1) then    ! if there is an electron on that site it can't be the result of PE
             PES_down(position,i) = 0                           ! zero everything then because it's not possible
             phase_PES_down(position,i) = 0                     ! zero everything then because it's not possible
          else
             do isite=position+1,ClusterSize                        ! +1 sicne the order is dn,up so wouldn't commute with the up electron on site=position
                ne = ne + ibits(fock_states(1,i),isite-1,1)    ! count the number of up electrons it will have to commute with to be removed
             end do
             do isite=position,ClusterSize
                ne = ne + ibits(fock_states(2,i),isite-1,1)    ! count the number of dn electrons it will have to commute with to be removed
             end do
             if (MOD(ne,2) == 0) then 
                phase_PES_down(position,i) = 1                 ! if it had to do even number of anti-commutations it is positive
             else
                phase_PES_down(position,i) = -1
             end if
             new_state(2) = ibset(fock_states(2,i),position-1)  ! add dn electron to that site 
             new_state(1) = fock_states(1,i)                    ! the up portion remains the same
             do j=1,istate
                if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
                   new_index = j                               ! find the index of the new state by compairing it to the entire FS basis
                end if
             end do
             PES_down(position,i) = new_index                    ! record the state that will when PE will become state i
          end if
       end do
    end do
    !-------Find the IPES tables----------
    do j=1,ClusterSize 
       do i=1,istate
          if (PES_down(j,i) /= 0) then
             phase_IPES_down(j,PES_down(j,i)) = phase_PES_down(j,i)
             IPES_down(j,PES_down(j,i)) = i
          end if
          if (PES_up(j,i) /= 0) then
             IPES_up(j,PES_up(j,i)) = i
             phase_IPES_up(j,PES_up(j,i)) = phase_PES_up(j,i)
          end if
       end do
    end do
    
    
    
  end subroutine transformations
  !************************************************************************************
  !************************************************************************************
  subroutine make_neighbours(nsites, neighbours)
    ! makes matrix that containes information about which sites are nearest neighbours
    ! neighbours(i,:) is a list of all the neighbours of site i. Each site has 2 nearest neighbours
    ! first and last site have neighbours with index 0 which will tell the program that it has no neighbour on that side
    implicit none
    integer, intent(in) :: nsites
    integer, intent(out), dimension(:,:) :: neighbours
    integer i
    
    neighbours = 0
    if ( nsites .ne. 1 ) then
       neighbours(1,2) = 2 
       neighbours(nsites,1) = nsites-1
    end if
    
    do i=2,nsites-1
       
       neighbours(i,1) = i-1
       
       
       neighbours(i,2) = i+1
       
    end do
    
  end subroutine make_neighbours
  !************************************************************************************
  !************************************************************************************
  subroutine Drop_Hopping(nsites, E, Hops, U)
    implicit none
    integer, intent(in) :: nsites
    real, intent(in) :: U, E(nsites)
    real, intent(inout) :: Hops(nsites-1)

    integer :: i     !Site
    integer :: j     !neighbour
    do i=1,nsites-1
       j=i+1

       if ( (E(i) - E(j) - U) .lt. drop_cutoff) then
          Hops(i) = 0.0
       else if ( (E(i) - E(j) + U) .lt. drop_cutoff) then
          Hops(i) = 0.0
       end if
    end do
    

  end subroutine Drop_Hopping
  
  !************************************************************************************
  !************************************************************************************
  subroutine matrix_sizes(ClusterSize)
    !  %------------------------------------------------------------------------------%
    !  |  This subroutine makes the array matrix_sizes which contains the dimensions  |
    !  |  of all the Hamiltonian submatrices                                          | 
    !  |                                                                              |
    !  |  matrix_sizes(i,j) is the size of the matrix for FS with i up electrons      |
    !  |  and j down electrons.                                                       |
    !  %------------------------------------------------------------------------------%
    implicit none
    integer, intent(in) :: ClusterSize
    integer :: n_up,n_dn
    
    
       
       msize(0:ClusterSize,0:ClusterSize) = 0; mblock(0:ClusterSize,0:ClusterSize) = 0
       
       do n_up=0,ClusterSize
          do n_dn=0,ClusterSize
             msize(n_up,n_dn) = choose(ClusterSize,n_up)*choose(ClusterSize,n_dn)
             if (n_dn == 0 .and. n_up == 0) then
                mblock(n_up,n_dn) = 1
             else if (n_dn == 0) then
                mblock(n_up,n_dn) = mblock(n_up-1,ClusterSize) + msize(n_up-1,ClusterSize)
             else 
                mblock(n_up,n_dn) = mblock(n_up,n_dn-1) + msize(n_up,n_dn-1)
             end if
          end do
       end do
    
    
  end subroutine matrix_sizes
  !************************************************************************************
  !************************************************************************************
  subroutine BuildHSUB(n_up, n_dn, HSUB, neighbours, U, E, Hops, nsites, nstates)
    implicit none
    integer, intent(in) :: n_up,n_dn
    real, dimension(:,:), intent(out) :: HSUB
    integer, intent(in) :: nsites, nstates
    integer, intent(in), dimension(:,:) :: neighbours
    real, intent(in) :: U, E(nsites)
    real, intent(in) :: Hops(nsites-1)       ! the hopping integral
    integer :: istate, isite, i, j, y        ! counters for loops
    integer :: inbr, hoplabel                ! site number of the neighbour
    integer :: new_state(2)                  ! new state after a hopping
    integer :: phase                         ! the number that keeps track of anti-commutations (+ or -)
    integer :: ne                            ! counts the number of electrons (used to calculate the phase)
    integer :: trans_site(2)                 ! holds the site number of the starting and ending site of a hopping (needed to find number of electrons)
    integer :: new_index                     ! the column index of the new FS after the hopping
    integer :: state_index                   ! the row index of the old FS before the hopping
    
    HSUB = 0.0
    
    do istate = mblock(n_up,n_dn),mblock(n_up,n_dn) + msize(n_up,n_dn)-1  ! loop over all the states in each submatrix
       do isite = 1,nsites                                               ! loop over each site of the state
          if (ibits(fock_states(1,istate),isite-1,1) == 1) then         ! check if there is a up electron on that site
             do y=1,size(neighbours,2)                                 ! if so loop over all the sites that nieghbour it
                new_state(1) = IBCLR(fock_states(1,istate),isite-1)   ! remove the up electron from current site
                inbr = neighbours(isite,y)                            ! this the site number of the neighbour
                if (inbr .eq. 0) CYCLE                                ! A neighbour index of 0 implies there is no neighbour (enforces closed boundary conditions)
                if (ibits(new_state(1),inbr-1,1) == 0) then           ! check if the neighbour site is empty (no up electron)
                   new_state(2) = fock_states(2,istate)              ! the state after the hopping has same down component
                   ne = 0                                            ! set the electron counter to zero (need to count them for anti-commutations)
                   trans_site(1) = inbr; trans_site(2) = isite       ! the site it the elelectron started at and finished at
                   do i=MINVAL(trans_site),MAXVAL(trans_site)-1
                      ne = ne + ibits(new_state(1),i-1,1)           ! count the number of up electrons between start and finish site
                      ne = ne + ibits(new_state(2),i-1,1)           ! count the number of down electrons between start and finish site
                   end do
                   if (MOD(ne,2) == 0) then                          ! if even amount of electrons between then anticommutation combine to 1
                      phase = 1
                   else
                      phase = -1                                    ! if odd amount of electrons between then anticommutation combine to -1
                   end if
                   new_state(1) = ibset(new_state(1),inbr-1)         ! find the up component of the new fock state
                   do j=1,nstates
                      if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
                         new_index = j                             ! search through the fock states to find the index of the new FS
                      end if
                   end do
                   state_index = istate + 1 - mblock(n_up,n_dn)      ! find the row and column to add the 't' to the Hamiltonian matrix 
                   new_index = new_index + 1 - mblock(n_up,n_dn)
                   hoplabel = isite + y - 2
                   HSUB(state_index,new_index) = Hops(hoplabel)*phase
                end if
             end do
          end if
          if (ibits(fock_states(2,istate),isite-1,1) == 1) then         ! Repeat identical process for hoppping of down electrons
             do y=1,size(neighbours,2)
                new_state(2) = IBCLR(fock_states(2,istate),isite-1)
                inbr = neighbours(isite,y)
                if ( inbr .eq. 0 ) CYCLE
                if (ibits(new_state(2),inbr-1,1) == 0) then
                   new_state(1) = fock_states(1,istate)
                   ne = 0
                   trans_site(1) = inbr; trans_site(2) = isite
                   do i=MINVAL(trans_site)+1,MAXVAL(trans_site)
                      ne = ne + ibits(new_state(1),i-1,1)     ! count the number of up electrons in that state
                      ne = ne + ibits(new_state(2),i-1,1)     ! count the number of down electrons in that state
                   end do
                   if (MOD(ne,2) == 0) then 
                      phase = 1
                   else
                      phase = -1
                   end if
                   new_state(2) = ibset(new_state(2),inbr-1)
                   do j=1,nstates
                      if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
                         new_index = j
                      end if
                   end do
                   state_index = istate + 1 - mblock(n_up,n_dn)
                   new_index = new_index + 1 - mblock(n_up,n_dn)
                   hoplabel = isite + y - 2
                   HSUB(state_index,new_index) = Hops(hoplabel)*phase
                end if
             end do
          end if
       end do
    end do
    do istate = 1,msize(n_up,n_dn)
       do isite=1,nsites                     ! loop over all the site of each state
          ne=0
          ne = ne + IBITS(fock_states(1,istate+mblock(n_up,n_dn)-1),isite-1,1)  ! check if up electron on that site of that state
          ne = ne + IBITS(fock_states(2,istate+mblock(n_up,n_dn)-1),isite-1,1)  ! check if down electron on that site of that state
          HSUB(istate,istate) = HSUB(istate,istate) + ne*E(isite)             ! add ne*(the site potential of that site)
          if (ne == 2) then 
             HSUB(istate,istate) = HSUB(istate,istate) + U                   ! if there is both an up and down electron add a U
          end if
       end do
    end do
  end subroutine BuildHSUB
  !************************************************************************************
  !************************************************************************************ 
  subroutine solve_hamiltonian1(E, U, mu, Hops, neighbours, e_ground, nsites, nstates)
    
    implicit none
    integer, intent(in) :: nsites, nstates
    real, intent(in) :: E(nsites), U, mu
    real, intent(in) :: Hops(nsites-1)                    ! the hopping integral
    integer, intent(in), dimension(:,:) :: neighbours
    real, intent(inout), dimension(0:nsites,0:nsites) :: e_ground       ! array of the lowest grand potential (Gpot) of each submatrix (e_ground(i,j) is lowest Gpot of Hij) 
    integer :: n_up,n_dn                     ! the number of up ad down electrons (used to loop over each submatrix)

    
    
    do n_up = 0,nsites             ! loop over all submatrices
       do n_dn = 0,nsites             ! loop over all submatrices
          Call GetEnergy(n_up,n_dn, neighbours, U, E, mu, Hops, nsites, nstates, e_ground)
       end do
    end do
    
  end subroutine solve_hamiltonian1
  !************************************************************************************
  !************************************************************************************
  subroutine GetEnergy(n_up,n_dn, neighbours, U, E, mu, Hops, nsites, nstates, e_ground)

    implicit none
    integer, intent(in) :: n_up, n_dn, nsites, nstates, neighbours(:,:)
    real, intent(inout), dimension(0:nsites,0:nsites) :: e_ground
    real, intent(in) :: E(nsites), U, mu, Hops(nsites-1)
    real :: HSUB(msize(n_up,n_dn),msize(n_up,n_dn)), VSUB(msize(n_up,n_dn),msize(n_up,n_dn))
    real :: WSUB(msize(n_up,n_dn))
    
    call BuildHSUB(n_up, n_dn, HSUB, neighbours, U, E, Hops, nsites, nstates)
    call ssyevr_lapack1(msize(n_up,n_dn),HSUB,WSUB,VSUB)
    e_ground(n_up,n_dn) = WSUB(1) - mu*(n_up+n_dn)
    
  end subroutine GetEnergy
  !************************************************************************************
  !************************************************************************************
  subroutine solve_hamiltonian2(E,U,mu,Hops,neighbours,nsites,nstates, &
       g_up,g_dn, eigenvectors, grand_potential)
    implicit none
    integer, intent(in) :: nsites, nstates
    real, intent(in) :: E(nsites), U, mu
    integer, intent(in) :: g_up,g_dn
    real, intent(in) :: Hops(nsites-1)                    ! the hopping integral
    integer, intent(in), dimension(:,:) :: neighbours
    real, intent(inout), dimension(:) :: grand_potential          ! grand potentials (eigenenergies - mu*number electrons)
    TYPE(cell), intent(inout), dimension(:) :: eigenvectors        ! the many body eigenvectors (MBE) only coefficients of basis states with same n_up,n_dn as it
    real :: HSUB(msize(g_up,g_dn),msize(g_up,g_dn)), VSUB(msize(g_up,g_dn),msize(g_up,g_dn)), WSUB(msize(g_up,g_dn))
    integer :: i
    do i=1,msize(g_up,g_dn)
       allocate(eigenvectors(i+mblock(g_up,g_dn)-1)%comp(1:msize(g_up,g_dn)))
    end do
    call BuildHSUB(g_up, g_dn, HSUB, neighbours, U, E, Hops, nsites, nstates)
    call ssyevr_lapack(msize(g_up,g_dn),HSUB,WSUB,VSUB)
    
    do i=1,msize(g_up,g_dn)
       grand_potential(i+mblock(g_up,g_dn)-1) = WSUB(i) - mu*(g_up+g_dn)  ! grand potentials
    end do
    
    do i=1,msize(g_up,g_dn)
       eigenvectors(i+mblock(g_up,g_dn)-1)%comp(1:msize(g_up,g_dn)) &
            = VSUB(1:msize(g_up,g_dn),i)  ! eigenvectors
    end do

    
    
  end subroutine solve_hamiltonian2
  !************************************************************************************
  !************************************************************************************
  subroutine DiagCluster(nsites,nstates,E,Energy,Weight,mu,U,Hops)
    implicit none
    integer, intent(in) :: nsites,nstates 
    real, intent(in) :: E(nsites), mu, U, Hops(nsites-1)
    real, dimension(nsites*2*nstates), intent(out) :: Energy, Weight
    
    
    real, dimension(nstates) :: grand_potential          ! grand potentials (eigenenergies - mu*number electrons)
    real :: grand_potential_ground=0.0                   ! the lowest grand ensemble energy
    real, dimension(0:nsites,0:nsites) :: e_ground       ! array of the lowest grand potential (Gpot) of each submatrix (e_ground(i,j) is lowest Gpot of Hij) 
    TYPE(cell), dimension(:), allocatable :: eigenvectors       ! the many body eigenvectors (MBE) only coefficients of basis states with same n_up,n_dn as it  
    integer, dimension(nsites,2) :: neighbours           ! neighbours(i,:) is all the site that are nearest neighbours to site i
    
    integer :: i, j, k                        ! counters for loops
    integer :: n_up, n_dn                        ! counters for loops (number of up or down electrons)
    integer :: g_up,g_dn                         ! number of electrons in the many-body ground state
    integer :: min_up, min_dn, max_up, max_dn    ! max/min number of electrons that a many-body eigenstate that can be transitioned to can have
    integer :: low, high                         ! range in array grand_potentials of states with a specified n_up, n_dn electrons
    integer :: groundloc(2)                      ! stores the location in array e_ground of minimum energy (this will find g_up, g_dn)
    integer :: location(1)                       ! stores the location in the grand_potential array of the lowest energy
    real, dimension(nstates) :: MBGvec           ! many-body ground state vector (MBG) written in fock basis (MBGvec(i) is coefficient of fock state "i")
    real, dimension(nsites,nstates) :: PESdn_MBG, PESup_MBG   ! MBG after a down or up photo emmision (PE) respectively (PESdn_MBG(i,:) is  c_{i,dn}|Psi0> )
    real, dimension(nsites,nstates) :: IPESdn_MBG, IPESup_MBG ! MBG after a down or up inverse photo emmision respectively (IPESdn_MBG(i,:) is  c^{dagger}_{i,dn}|Psi0> )
    real, dimension(nsites,2*nstates,2) :: LDOS      ! local density of states (LDOS(i,:) is LDOS of site i)
    real :: inner_prod_up, inner_prod_dn             ! inner products used when calculating weight of LDOS contributions (<Psi|PESdn_MBG> or <Psi|IPESdn_MBG>)

    allocate(eigenvectors(nstates))
    call make_neighbours(nsites, neighbours)
    MBGvec=0.0
    grand_potential_ground = 0.0
    LDOS = 0.0
    PESdn_MBG=0.0; PESup_MBG=0.0; IPESdn_MBG=0.0; IPESup_MBG=0.0
    grand_potential = 0
    e_ground = 0
    Energy = 0.0
    Weight = 0.0
    
    call solve_hamiltonian1(E, U, mu, Hops, neighbours, e_ground, nsites, nstates)  ! solve for the lowest grand potential of each hamiltonian sub-matrix
    
    
    grand_potential_ground = minval(e_ground)
    
    groundloc = minloc(e_ground)
    
    g_up = groundloc(1) - 1
    g_dn = groundloc(2) - 1
    
    min_up = MAX(0,g_up-1)             
    max_up = MIN(nsites,g_up+1)
    min_dn = MAX(0,g_dn-1)
    max_dn = MIN(nsites,g_dn+1)
    
    location = mblock(g_up,g_dn)       ! find the location of the lowest grand_potential in the array grand_potential_ground  
    
    call solve_hamiltonian2(E,U,mu,Hops,neighbours,nsites,nstates,&
         g_up,g_dn,eigenvectors,grand_potential)

    !---------Solve for eigenvectors and eigenvalues all sub hamiltonians with +/- 1 electron as MBG--------------------
    if (g_up /= 0) call solve_hamiltonian2(E,U,mu,Hops,neighbours,nsites,nstates,&
         g_up-1,g_dn,eigenvectors,grand_potential)
    
    if (g_dn /= 0) call solve_hamiltonian2(E,U,mu,Hops,neighbours,nsites,nstates,&
         g_up,g_dn-1,eigenvectors,grand_potential)
    
    if (g_up /= nsites) call solve_hamiltonian2(E,U,mu,Hops,neighbours,nsites,nstates,&
         g_up+1,g_dn,eigenvectors,grand_potential)
    
    if (g_dn /= nsites) call solve_hamiltonian2(E,U,mu,Hops,neighbours,nsites,nstates,&
         g_up,g_dn+1,eigenvectors,grand_potential)
    
    high = mblock(g_up,g_dn) + msize(g_up,g_dn) - 1                                 ! find range of indexs of fock states in the MBG's submatrix
    MBGvec(mblock(g_up,g_dn):high) = eigenvectors(location(1))%comp(1:msize(g_up,g_dn))   ! set MBGvec to the eigenvector corresponding to the lowest energy
    
    !------------------calculate PESdn_MBG, PESup_MBG, IPESup_MBG, IPESdn_MBG------------------------------------------
    
    do j=1,nsites
       do i=1,nstates
          if (PES_up(j,i)==0) then
             PESup_MBG(j,i) = 0.0
          else 
             PESup_MBG(j,i) = MBGvec(PES_up(j,i))*phase_PES_up(j,i)
          end if
          if (PES_down(j,i)==0) then
             PESdn_MBG(j,i) = 0.0
          else 
             PESdn_MBG(j,i) = MBGvec(PES_down(j,i))*phase_PES_down(j,i)
          end if
          if (IPES_up(j,i)==0) then
             IPESup_MBG(j,i) = 0.0
          else 
             IPESup_MBG(j,i) = MBGvec(IPES_up(j,i))*phase_IPES_up(j,i)
          end if
          if (IPES_down(j,i)==0) then
             IPESdn_MBG(j,i) = 0.0
          else 
             IPESdn_MBG(j,i) = MBGvec(IPES_down(j,i))*phase_IPES_down(j,i)
          end if
       end do
    end do
    
    
    !---------------------------calculate the LDOS for all the sites--------------------------------------------------
    
    do n_up=min_up,max_up                ! loop over all possible submatrices
       do n_dn=min_dn,max_dn
          if (n_up == min_up .and. n_dn == min_dn .and. g_up /= 0 .and. g_dn /= 0) CYCLE             ! not allowed (two electrons less then MBG)
          if (n_up == max_up .and. n_dn == max_dn .and. g_up /= nsites .and. g_dn /= nsites) CYCLE   ! not allowed (two electrons more then MBG)
          if (n_up == max_up .and. n_dn == min_dn .and. g_up /= nsites .and. g_dn /= 0) CYCLE        ! not allowed (two electrons different then MBG)
          if (n_up == min_up .and. n_dn == max_dn .and. g_up /= 0 .and. g_dn /= nsites) CYCLE        ! not allowed (two electrons different then MBG)
          if (n_up == g_up .and. n_dn == g_dn) CYCLE                                                 ! not allowed (same n_up and n_dn)
          low = mblock(n_up,n_dn)                               ! lowest value of range of fock states of the submatrix
          high = mblock(n_up,n_dn) + msize(n_up,n_dn) - 1       ! highest value of range of fock states of the submatrix
          do j=1,nsites                  ! loop over all the sites (PESdn_MBG for c_{j,sigma} with different j's)  
             do i=low,high              ! loop over only the states that will be non-zero (within range of submatrix)
                inner_prod_up = 0
                inner_prod_dn = 0
                if (n_up == min_up) then
                   inner_prod_up = (dot_product(PESup_MBG(j,low:high),eigenvectors(i)%comp(1:msize(n_up,n_dn))))**2
                end if
                if (n_dn == min_dn) then
                   inner_prod_dn =  (dot_product(PESdn_MBG(j,low:high),eigenvectors(i)%comp(1:msize(n_up,n_dn))))**2
                end if
                LDOS(j,i,1) = grand_potential_ground - grand_potential(i)              ! location of the peak
                LDOS(j,i,2) = (inner_prod_up + inner_prod_dn)*0.5                      ! weight of the peak (average up and down spin components)
               
                inner_prod_up = 0
                inner_prod_dn = 0
                if (n_up == max_up) then
                   inner_prod_up = (dot_product(IPESup_MBG(j,low:high),eigenvectors(i)%comp(1:msize(n_up,n_dn))))**2
                end if
                if (n_dn == max_dn) then
                   inner_prod_dn =  (dot_product(IPESdn_MBG(j,low:high),eigenvectors(i)%comp(1:msize(n_up,n_dn))))**2
                end if
                LDOS(j,i+nstates,1) = grand_potential(i) - grand_potential_ground       ! location of the peak
                LDOS(j,i+nstates,2) = (inner_prod_up + inner_prod_dn)*0.5               ! weight of the peak (average up and down spin components)
             end do
          end do
       end do
    end do
    
    k=1
    deallocate(eigenvectors)
    do i=1,nsites
       do j = 1,2*nstates
          Energy(k) = LDOS(i,j,1)
          Weight(k) = LDOS(i,j,2)
          k = k + 1
       end do
    end do

    
  end subroutine DiagCluster
  
  !************************************************************************************
  
  subroutine ssyevr_lapack1(dim,matrix,eigvalues,eigvectors)
    
    !  %-----------------------------------------------------------------%
    !  |  This subroutine calls the driver SSYEVR.                       |
    !  |                                                                 |
    !  |  Solves the only the lowest eigenvalues of the matrix.          |
    !  |      - Eigenvalues returned in eigvalues                        |
    !  %-----------------------------------------------------------------%
    
    
    implicit none
    
    integer, intent(in) :: dim                  ! dimension of the matrix
    real, intent(in) :: matrix(dim,dim)         ! the matrix to be solved
    real, intent(out) :: eigvalues(dim)         ! the outputed eigenvalues
    real, intent(out) :: eigvectors(dim,dim)    ! array contains no useful information this case
    
    integer, parameter :: NSELECT=1
    
    integer :: LDA, LDZ
    real :: VL,VU
    real :: ABSTOL = -1
    integer :: IL=1,IU=NSELECT
    
    integer :: M
    integer :: INFO
    integer :: LWORK, LIWORK
    integer, allocatable, dimension(:) :: ISUPPZ, IWORK
    real, allocatable,dimension(:) :: WORK
    eigvectors = 0.0
    eigvalues = 0.0   
    if (dim == 1) then
       eigvectors = 1
       eigvalues = matrix(1,1)
       return
    end if
    
    
    LDA = dim; LDZ = dim
    LWORK = -1
    LIWORK = -1
    
    allocate(ISUPPZ(2*dim))
    allocate(WORK(1),IWORK(1))
    WORK = 0.0; IWORK = 0
    call ssyevr('N','I','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    LWORK= int(WORK(1))
    LIWORK = IWORK(1)
    
    deallocate(WORK,IWORK)
    allocate(WORK(LWORK),IWORK(LIWORK))
    
    call ssyevr('N','I','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    deallocate(ISUPPZ,WORK,IWORK)
  end subroutine ssyevr_lapack1
  
  !************************************************************************************
  
  subroutine ssyevr_lapack(dim,matrix,eigvalues,eigvectors)
    
    !  %-----------------------------------------------------------------%
    !  |  This subroutine calls the driver SSYEVR.                       |
    !  |                                                                 |
    !  |  Solves the eigenvalues and eigenvectors of the matrix.         |
    !  |      - Eigenvalues returned in eigvalues                        |
    !  |      - Eigenvectors returned in eigvectors                      |
    !  %-----------------------------------------------------------------%
    
    implicit none
    
    integer, intent(in) :: dim
    real, intent(in) :: matrix(dim,dim)
    real, intent(out) :: eigvalues(dim)
    real, intent(out) :: eigvectors(dim,dim)
    
    integer, parameter :: NSELECT=1
    
    integer :: LDA, LDZ
    real :: VL,VU
    real :: ABSTOL = -1
    integer :: IL=1,IU=1
    
    integer :: M
    integer :: INFO
    integer :: LWORK, LIWORK
    integer, allocatable, dimension(:) :: ISUPPZ, IWORK
    real, allocatable,dimension(:) :: WORK
    
    eigvalues = 0.0
    eigvectors = 0.0
    
    if (dim == 1) then
       eigvectors = 1
       eigvalues = matrix(1,1)
       return
    end if
    
    
    LDA = dim; LDZ = dim
    LWORK = -1
    LIWORK = -1
    
    allocate(ISUPPZ(2*dim))
    allocate(WORK(1),IWORK(1))
    WORK = 0.0; IWORK = 0
    call ssyevr('V','A','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    LWORK= int(WORK(1)) * 2
    LIWORK = IWORK(1) * 2
    
    deallocate(WORK,IWORK)
    allocate(WORK(LWORK),IWORK(LIWORK))
    
    call ssyevr('V','A','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    deallocate(ISUPPZ,WORK,IWORK)
  end subroutine ssyevr_lapack
  !************************************************************************************
  integer function choose(j,k)
    ! %---------------------------------------------%
    ! |  simple choose function from statistics:    |
    ! |                                             |
    ! |                         i!                  |
    ! |       choose(i,j) =  --------               |
    ! |                      j!(i-j)!               |
    ! %---------------------------------------------%
    implicit none
    
    integer  :: j,k,i2
    integer  :: tj,tk,tchoose
    tchoose = 1
    tj = j
    tk = k
    do i2 = tj-tk+1,tj
       tchoose = tchoose * i2
    end do
    do i2 = 2, tk
       tchoose = tchoose / i2
    end do
    choose = tchoose
  end function choose
  !************************************************************************************
  
  
end module Diag


