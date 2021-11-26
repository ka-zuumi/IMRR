module ANALYSIS
use DOUBLE
implicit none

!This module is inteded to be changed according to what the user
!wants to analyze from a grid that was made.
!Analysis is done in roughly the order they appear here

!The name of this experiment
integer,parameter :: expfolder_length = 56
character(expfolder_length),parameter :: expfolder = "expcompareGradientsInGridtoNWChemGradients0_history3_20/"

!Set the number of grids to be analyzed; will start at 001 and increment
!If this number is larger than the number of grids in the folder,
!then it will default to the number of grids in the folder
integer,parameter :: Ngrid_cap = 1
!The number of grids we will end up using (never more than Ngrid_cap)
integer :: Ngrid_total
!Filechannels for each of those grids
integer,allocatable :: filechannels(:)

!Set the number of children cells to be checked
integer,parameter :: Norder_cap = 1

   !This takes much more time but you can force the checkState subroutine
   !to also check the rmsd of frames in adjacent cells
   logical :: force_Neighbors = .true.

   logical :: force_Permutations = .false.

   logical :: memory_flag = .true.
   integer :: single_index_max
   integer,allocatable :: populationbuffer2(:)

   !The point at which checkState stops checking neighbors
   !is determined here
   !The default (if none of these are set) is all zeros
   integer,parameter :: ssm_length = 2
   integer,dimension(ssm_length) :: ssm1 = (/ 24, 0 /)
   integer,dimension(ssm_length) :: ssm2 = (/ 00, 00 /)

      !Set the threshold RMSD to be used for any rejection method
      !real(dp),parameter :: !threshold_rmsd! = !.200100d0
      real(dp),parameter :: outer_threshold_SI = .15000d0
      real(dp),parameter :: inner_threshold_SI = .000000d0
      integer,parameter :: Ninterpolation_threshold = 0
      real(dp),parameter :: R1_threshold_SI = 3.0d3
      real(dp),parameter :: R2_threshold_SI = 3.0d3
      real(dp),parameter :: R1plusAlphaR2_threshold_SI = 3.0d3
      integer :: Nsort = 1

      !Set .true. to generate trajectories using md-calculated gradients
      !Otherwise, the program will use the above threshold as a rejection
      !method
      logical :: reject_flag = .false.
      logical :: readtrajectory_flag = .true.
      logical :: readtrajectory_onlyInitialQPandPDOT_flag = .true.
      character(len=*),parameter :: readtrajectory_path = ""
      integer,parameter :: filechannel5 = 78
      integer,parameter :: filechannel6 = 79

      !If reject_flag is false (and we are accepting frames) then
      !accept_first controls whether we use the first frame accepted
      !or use the frame that is closest in RMSD
      logical :: accept_first = .false.

      !If reject_flag is false (and we are accepting frames) then
      !accept_worst indicates to accept the frame with the
      !maximum RMSD (instead of the minimum RMSD)
      logical :: accept_worst = .false.

      !For diversity
      logical :: diversity_flag = .true.

      !For the selection strategy
      integer :: ssID = 2

      !For the replacement strategy
      integer :: rsID = 2

      !Set .true. if the real force calculations we do should be
      !added to the grid
      integer :: grid_addition = 0

         !Set .true. if interpolation should be used; that is to say
         !a weighted combination of acceptable frames are used to
         !calculate an approximate gradient
         logical :: interpolation_flag = .true.

         !Interpolation requires a scaling parameter for the weights
         !This is a positive, nonzero real number
         !(Now deprecated from introduction of 2nd order LS minimization)
         real(dp) :: interpolation_alpha1 = 2.0d0

         !For persistent data collection and analysis, a new file
         !is dedicated to interpolation results
         !Whether we record or not to this file is governed by the
         !gather_interpolation_flag
         logical :: gather_interpolation_flag = .true.
         character(17),parameter :: interpolationfile = "interpolation.dat"
         character(14),parameter :: interpolationfolder = "interpolation/"

            real(dp) :: alpha_ratio = 1.0d-1

            integer,parameter :: Ninterpolation_max = 30 !400
            integer,parameter :: Ninterpolation_cap = 30 !200

            integer,parameter :: Naccept_max = 2
            real(dp),allocatable :: trajRMSDbuffer(:,:)

character(len=:), allocatable :: gridpath4
character(len=:), allocatable :: gridpath5

integer :: Naccept
integer :: Nalpha_tries
logical :: temp_reject_flag
real(dp) :: Vprev
real(dp),allocatable :: qprev(:),pprev(:),pdotprev(:),ftprev(:)
real(dp),allocatable :: qhistory(:,:),phistory(:,:)
real(dp),allocatable :: pdothistory(:,:),fthistory(:,:)

integer :: MLmethodID = 1
integer :: RKHS_updateInverseKernel_counter = 0
integer :: RKHS_KernelSize
logical :: useLocalRKHS_flag = .false.

! The 'coordinates' to use for the IMRR
! 0 = the 3N xyz coordinates
! 1 = distance matrix jacobian product
integer :: coordinatesID = 0

! Related to IMRR use in a trajectory:
logical :: usePreviousFrameInTrajectory_flag = .true.
logical :: usePreviousFramesExactlyOnce_flag = .true.
integer :: Nhistory = 0
integer,parameter :: Nhistory_max = 3

logical :: checkNoState_OnlyHistory_flag = .false.

logical :: historyMayBeIMRR_flag = .false.
logical :: dontUseRecentHistoricalIMRR_flag = .false.
integer :: NremoveEvenMoreHistoricalIMRR = 0

contains


subroutine initializeTrajectory(current_Ntraj)
use PARAMETERS
implicit none

!COMMON/PRLIST/T,V,H

double precision :: T, V, H

character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(trajectory_text_length) :: Ntraj_text

character(200) :: readtrajectory_aline

integer,intent(in) :: current_Ntraj
integer :: cr
integer :: n,m

!Set the trajectory number
Ntraj = current_Ntraj

!Erase the memory buffer
!if (memory_flag) then
    populationbuffer2 = -1
!end if

!If moving on to a new trajectory,
!close the old trajectory data files
if (Ntraj > 1) then
    do m = 1, Ngrid_max
        close(filechannels(1+m))
    end do
end if

!And open the new trajectory data files
!if we are not done
write(variable_length_text,FMT=FMT5_variable)&
        trajectory_text_length
write(Ntraj_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        Ntraj

write(variable_length_text,FMT=FMT5_variable)&
        Ngrid_text_length
do m = 1, Ngrid_max
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            m
    open(filechannels(1+m),file=gridpath5//&
            Ngrid_text//"_"//Ntraj_text//".dat",&
            position="append")
!           Ngrid_text//"_"//Ntraj_text//".dat")
end do

Naccept = 0
Nalpha_tries = 1
temp_reject_flag = .true.

allocate(qprev(3*NATOMS),pprev(3*NATOMS),&
         pdotprev(3*NATOMS),ftprev(3*NATOMS))
allocate(qhistory(3*NATOMS,Nhistory_max),&
         phistory(3*NATOMS,Nhistory_max),&
         pdothistory(3*NATOMS,Nhistory_max),&
         fthistory(3*NATOMS,Nhistory_max))

if (readtrajectory_flag) then
    read(filechannel5,FMT="(A)") readtrajectory_aline
    open(filechannel6,file=trim(adjustl(readtrajectory_aline)))
end if

return

end subroutine initializeTrajectory

subroutine finalizeGridCheck()
use PARAMETERS
implicit none

integer :: m

do m = 1, Ngrid_max
    close(filechannels(1+m))
end do

if (gather_interpolation_flag) then
    close(filechannel3)
end if

if (.true.) then
    close(filechannel4)
end if

if (readtrajectory_flag) then
    close(filechannel5)
    close(filechannel6)
end if

return

end subroutine finalizeGridCheck

subroutine setSubcellSearchMax()
use PARAMETERS
use FUNCTIONS
implicit none

integer :: i,j
integer :: single_index
integer,dimension(Nvar) :: var_index

do i = 1, min(Norder_max+1,ssm_length)
    subcellsearch_max(i) = ssm1(i)
    subcellsearch_max1(i) = ssm1(i)
    subcellsearch_max2(i) = ssm2(i)
end do

!if (memory_flag) then

    !Initialize the memory buffer assuming
    !only subcellsearch_max1 and the first
    !order search are used
    j = subcellsearch_max(&
            Norder_order(1)+1)
    single_index_max = 0

    do i = 1, Nvar

        var_index = 0

        var_index(i) = j
        call getFlattened(Nvar,var_index,&
                single_index)
        single_index_max = max(single_index,&
                single_index_max)

        var_index(i) = -j
        call getFlattened(Nvar,var_index,&
                single_index)
        single_index_max = max(single_index,&
                single_index_max)

    end do

!end if

return

end subroutine setSubcellSearchMax

end module ANALYSIS



