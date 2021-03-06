module Tools
  USE Inputs
  implicit none
  SAVE

contains
  
  !     Bin_Data subroutine enables the binning of data into a histogram
  !
  !
  !     ------------------------------Arguments-------------------------------------
  !
  !     HistoData: A 1xN array, where N is the number of bins in the histogram. N is usually the parameter bins
  !
  !     Data: A 1xM array, where M is the number of contributions to be added to HistoData. Each element is the value of a contribution
  !
  !     Weights: (OPTIONAL) A 1xM array containing the weights of each contribution in Data. If not present, the weight of a contribution is 1
  !              
  !
  !     Dropped: (OPTIONAL) A variable of type integer that keeps track of the number of contributions outside of a specified range
  !              If not present, dropped contributions are ignored
  !
  !     Max/Min: Type real. The maximum (Max) and minimum (Min) values to be counted. Contributions outside are not counted
  
  subroutine Bin_Data( HistoData, Data, Weights, Dropped, Max, Min, bins )
    implicit none
    real(dp), dimension(:), intent(in) :: Data
    real(dp), dimension(:), Optional, intent(in) :: Weights
    integer, intent(in) :: bins
    real(dp), intent(inout) :: HistoData(bins)
    real, intent(in) :: Max
    real, intent(in) :: Min
    real(dp), optional, intent(inout) :: Dropped
    integer Loop1, Loop2 !Loop integers
 
    ! since both Dropped and Weights are optional there are 2^2 different versions of this subroutine.

    !First is if both Dropped and Weights are present
    If ( Present(Dropped) .and. Present(Weights) ) then
       
       do Loop1 = 1,size(Data)                                                !Loop through contributions
          if ( Weights(Loop1) .eq. 0 ) CYCLE                                  !If a contribution has a weight of 0, ignore.
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) then 
             Dropped = Dropped + Weights(Loop1)                               !If a contribution is \lt or \gt min/max increment Dropped
          else
             Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)          !Find bin number of contribution
             if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) then                !Had some precision errors near the edge of the band, this ignore those
                Dropped = Dropped + Weights(Loop1)
             end if
             
             HistoData(Loop2) = HistoData(Loop2) + Weights(Loop1)             !Increment Loop2 bin of HistoData by the Weight is the Loop2 element
             
          end if
          
       end do

    !Second is if Dropped is present and Weights is not
    else if ( Present(Dropped) .and. (.not. Present(Weights) ) ) then

       do Loop1 = 1,size(Data) 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) then
             Dropped = Dropped + 1
          else             
             Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
             if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE
             HistoData(Loop2) = HistoData(Loop2) + 1
             
          end if
          
       end do

    !Third is if Dropped is not present and Weights is present
    else if ( (.not. Present(Dropped)) .and. Present(Weights) ) then

       do Loop1 = 1,size(Data) 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) CYCLE
          Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
          if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE
          HistoData(Loop2) = HistoData(Loop2) + Weights(Loop1)
             
       end do
          
    !Last is if neither Dropped nor Weights are present.
    else if ( .not. ( Present(Dropped) .or. Present(Weights) ) ) then

       do Loop1 = 1,size(Data) 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) CYCLE
          Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
          if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE   
          HistoData(Loop2) = HistoData(Loop2) + 1
          
       end do

    else
       print*, "Bin_Data options exceeded"
       STOP
    end if
    
     

  end subroutine Bin_Data
  
  !     OpenFile creates a file with comments detailing system parameters and contents of the file
  !
  !
  !     ------------------------------Arguments-------------------------------------
  !
  !     Unum : Integer variable that contains the unit number for the file. Usually on the order of 100
  !
  !     Name : Character variable with any length that specificies the base file name.
  !
  !     Contents : Character variable with any length that names the data inside
  !
  !     Column1 : Character Variable with any length. Specifies the information in the first column
  !
  !     Column2 : Character Variable with any length. Specifies the information in the second column

  
  subroutine OpenFile( Unum, Name, Contents, Column1, Column2, num_procs )
    implicit none
    integer, intent(in) :: Unum
    character(len=*), intent(in) :: Name, Contents, Column1, Column2
    integer, intent(in), optional :: num_procs
    character(len=25) :: filename
    character(len=8) :: FormatP, FormatB
    filename = Name
    filename = FileNamer(filename)
    open (unit = Unum, file = "output/"//trim(filename)//"" ,status ='unknown' )

    FormatP = '(A,F6.1)'
    FormatB = '(A,F6.2)'
    
    write(Unum,'(A,A)') "#This file contains the following data type for the ensemble: ", Trim(Contents)
    write(Unum,'(A,A)') "#First column is: ", Trim(Column1)
    write(Unum,'(A,A)') "#Second column is: ", Trim(Column2)
    write(Unum,FormatP) "#Disorder Strength = ", DELTA
    write(Unum,FormatB) "#Bond Cutoff = ", bond_cutoff
    write(Unum,FormatB) "#Pruned bonds cutoff = ", prune_cutoff
    write(Unum,'(A)') "#Fraction of n.n. sites in which at least one orbital pair ignored = "
    write(Unum,'(A,F6.3)') "#", real(PrunedBonds)/real(dim*systemn)
    write(Unum,'(A,F6.3)') "#Fraction of prunings that are useful = ", real(StrongestBondsPruned)/real(PrunedBonds)
    write(Unum,FormatP) "#Hopping = ", hop
    write(Unum,FormatP) "#Interaction strength = ", uSite
    write(Unum,FormatP) "#Chemical Potential = ", ChemPot
    write(Unum,'(A,I4)') "#Dimensions = ", dim
    if ( Present(num_procs) ) then
       write(Unum,'(A,I4)') "#Number of processes = ", num_procs
       write(Unum,'(A,I10)') "#Number of systems = ", num_procs*systemn
    else if ( .not. Present(num_procs) ) then
       write(Unum,'(A,I10)') "#Number of systems = ", systemn
    end if

 
    

  end subroutine OpenFile
  
  !     PrintFile prints the data to the file created in OpenFile
  !
  !
  !     ------------------------------Arguments-------------------------------------
  !
  !     Unum : Integer variable that contains the unit number for the file. Usually on the order of 100
  !
  !     Form : Character variable of any length. Specifies the form the printed data will take in the file
  !
  !     Max/Min : real variable specifying the max and min values for the data being printed. Used for bin sizes
  !
  !     BinNum : Integer variable. Number of bins to be printed.
  !
  !     Data : Real rank-one array. Data to be printed.
  !
  !     Dropped : (OPTIONAL) real variable. Counts the combined weight that is less than min or greater than max for the purposes of normalization
  
  subroutine PrintData( Unum, Form, Min, Max, BinNum, Data_DOS, Data_dp, Data_int, Dropped )
    implicit none
    integer, intent(in) :: Unum, BinNum
    character(len=*), intent(in) :: Form
    real, intent(in) :: Min, Max
    real(dp), optional, intent(in) :: Dropped
    real(dp), optional, intent(in), dimension(:) :: Data_DOS, Data_dp
    integer(ip), optional, intent(in), dimension(:) :: Data_int
    real(dp) :: width
    integer :: Loop1
    integer(ip) :: sum1

    !This is Dropped is present
    if ( Present(Data_DOS) ) then
       if ( Present(Dropped) ) then
          width = (Max-Min)/BinNum
          do Loop1=1,BinNum
             write(Unum,(Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, Data_DOS(Loop1)/ &
                  ( (Sum(Data_DOS) + Dropped) * width )
          end do
       else
          !If not present this is used. Difference is minor but not dividing by Dropped or the sum of data.
          !Second part is to enable cumulative DOS (DOS with clustermax=1,2,3,4,5 normalized to the DOS with clustermax = 6)
          width = (Max-Min)/BinNum
          do Loop1=1,BinNum
             write(Unum,( Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, Data_DOS(Loop1)/( width )
          end do
       end if
    end if

    if ( Present(Data_dp) ) then
       width = (Max-Min)/BinNum
       do Loop1=1,BinNum
          write(Unum,(Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, Data_dp(Loop1)/( Sum(Data_dp)*width )
       end do
    end if

    if (Present(Data_int) ) then
       sum1=0
       do Loop1=1,BinNum
          sum1 = sum1 + Data_int(Loop1)
       end do
       width = (Max - Min)/real(BinNum)
       do Loop1=1,BinNum
          write(Unum,(Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, Data_int(Loop1)/(sum1*width)
       end do
    end if
    
    Close(Unum)

  end subroutine PrintData

  !     resize_array removes a number of elements from a rank-one array
  !
  !
  !     ------------------------------Arguments-------------------------------------
  !
  !     array : Real rank-one array that is to be resized
  !
  !     numRemove : Integer. The number of elements to remove starting from the element Start in array
  !
  !     Start : Integer. The first element to remove from array

  subroutine resize_array(array, numRemove, Start)
    real(dp), dimension(:), allocatable :: tmp_arr
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: numRemove                           ! Number of elements to remove from 'array'
    integer, intent(in) :: Start                               ! The element in which to start removing the array
    integer i, j                                               ! Looping integer
    
    allocate(tmp_arr(size(array) - numRemove))

    tmp_arr(:) = 0.d0

    j=1

    if ( (Start + numRemove) .le. size(array) ) then


       do i = 1, Start - 1

          tmp_arr(i) = array(i)
          j = j + 1

       end do


       do i = (Start + numRemove), Size(array)

          tmp_arr(j) = array(i)
          j = j + 1

       end do



    else if ( (Start + numRemove) .gt. size(array) ) then

       do i = (Start + numRemove - size(array)) , Start-1

          tmp_arr(j) = array(i)
          j = j + 1

       end do


    end if

    deallocate(array)
    allocate(array(size(tmp_arr)))
    
    array = tmp_arr
    
   
  end subroutine resize_array
  
  subroutine DefineCluster( SitePotential, Sites, weakL, weakR, ClusterSize )
    implicit none

    ! Inputs
    real(dp), dimension(dim), intent(in) :: SitePotential
    integer, intent(in) :: weakL, weakR
    integer, intent(in) :: ClusterSize

    ! Outputs
    real(dp), dimension(ClusterSize), intent(out) :: Sites

    ! Other
    integer :: EndSite
    
    if ( (weakL .ge. weakR) .and. (ClusterSize .gt. 1) ) then !2, 3 site clusters which wrap around the "edge" of the system. As well as 4-site clusters with a single weak bond
       EndSite = dim - weakL
       Sites(1:EndSite) = SitePotential((weakL+1):dim)
       Sites((EndSite+1):ClusterSize) = SitePotential(1:weakR)
    else if ( ClusterSize .eq. 1 ) then !1-site clusters
       Sites = SitePotential(weakR:weakR)
    else if ( (ClusterSize .eq. dim) .and. (weakL .eq. 0) ) then !4-site clusters with no weak bonds
       Sites = SitePotential
    else !All other clusters
       Sites = SitePotential((weakL+1):weakR)
    end if

  end subroutine DefineCluster    
       

  !     str(k) takes an integer k as input and returns a string with the same integer
  
  character(len=20) function str(k)                  
    integer, intent(in) :: k
    
    write (str, *) k
    str = adjustl(str)
  end function str
  
  !     FileNamer(Name) takes a base file name as input and adds a counter starting at 0 to the file name. If that exists, the counter is incremented by 1 (++)
  
  character(len=25) function FileNamer(Name)                                ! Increments filename labels Blah023.dat label:023
    character(len=25), intent(in) :: Name                                   ! Initial filename. Would be 'Blah' from above example
    integer num                                                             ! Keeps track of the label increment as this continues to check for name existence
    logical Exist                                                           ! Used to store if a filename exists
    
    Exist = .true.
    num = 0
    FileNamer = ""//trim(Name)//"000.text"                                  ! Takes Blah and puts it into the format 'Blah000.text'
    
    do while (Exist)                                                        ! Continue to change filename if the current one exists
       inquire( file = "output/"//trim(FileNamer)//"", exist = Exist )                           ! Checks the current folder (where this .f90 file is found) if the current filename exists
       if ( Exist ) then                                                    ! If the file name exists
          num = num + 1                                                     ! Increment the label
          
          if ( len(trim(str(num))) .eq. 1 ) then                            ! Labels can be from 000 to 999. To preserve filename length the leading zeros are required
             FileNamer = ""//trim(Name)//"00"//trim(str(num))//".text"      ! This if statement (and the one below) check the length of the label (1 has length 1, 23 has length 2, 450 has length 3)
          else if ( len(trim(str(num))) .eq. 2) then                        ! This determines the number of leading zeros to add to the filename label
             FileNamer = ""//trim(Name)//"0"//trim(str(num))//".text"
          else
             FileNamer = ""//trim(Name)//trim(str(num))//".text"
          end if
          
       else if ( .not. Exist ) then                                         ! If it does not exist. Note: Didn't increment 'num'
          if ( len(trim(str(num))) .eq. 1 ) then
             FileNamer = ""//trim(Name)//"00"//trim(str(num))//".text"
          else if ( len(trim(str(num))) .eq. 2) then
             FileNamer = ""//trim(Name)//"0"//trim(str(num))//".text"
          else
             FileNamer = ""//trim(Name)//trim(str(num))//".text"
          end if
          Exist = .false.
       else
          print*, "You belong in a museum. Error: FileNameErr"              ! Outputs error if something strange happens. Can't imagine what
       end if
    end do
  end function FileNamer
  
  
  
end module Tools
 
 
 
