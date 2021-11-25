
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               interactMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This module governs how a frame interacts with a grid
!
!               A frame can be added to the grid or it can be
!               checked alongside the grid (for the closest frame)
!
!               The grid has a particular file formatting system
!               which stores frames according to the variables they are
!               associated with (ex. r1, r2); this formatting is
!               initialized in PARAMETERS
!
!               The grid also has an internal file counting system
!               which keeps track of the size of files, whether they
!               have children or not, and which children have a file
!               associated with them; this system is initialized in
!               PARAMETERS and reset whenever a grid is made 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               FILECHANNELS(2:Ngrid_total+1)   WRITE
!               FILECHANNELS(1)                 OPEN, READ, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!		checkState                      vals                            intent(in),REAL(DP),DIM(Nvar)
!                                               coords                          intent(in),REAL(DP),DIM(3,Natoms)
!                                               gradient                        intent(inout),REAL(DP),DIM(3,Natoms)
!                                               min_rmsd                        intent(inout),REAL(DP)
!                                               filechannels                    intent(in),INTEGER,DIM(Ngrid_total+1)
!
!                                               number_of_frames                intent(out),INTEGER
!                                               order                           intent(out),INTEGER
!                                               neighbor_check                  intent(out),INTEGER
!
!		getNeighbors                    scaling                         intent(in),REAL
!                                               multiplier                      intent(in),REAL
!                                               FMTorder                        intent(in),CHARACTER(6)
!                                               index                           intent(in),INTEGER
!                                               filechannel                     intent(in),INTEGER
!
!                                               coords                          intent(in),REAL(DP),DIM(3,Natoms)
!                                               min_rmsd                        intent(inout),REAL(DP)
!                                               gradient                        intent(inout),REAL(DP),DIM(3,Natoms)
!                                               U                               intent(inout),REAL(DP),DIM(3,3)
!                                               number_of_frames                intent(inout),INTEGER
!
!		getRMSD_dp                      filechannel_thread              intent(in),INTEGER
!                                               subcell                         intent(in),CHARACTER(*)
!
!                                               coords                          intent(in),REAL(DP),DIM(3,Natoms)
!                                               population                      intent(out),INTEGER
!                                               min_rmsd                        intent(inout),REAL(DP)
!                                               gradient                        intent(inout),REAL(DP),DIM(3,Natoms)
!                                               U                               intent(inout),REAL(DP),DIM(3,3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               getVarMaxMin                    VARIABLES
!
!               getNeighbors                    interactMultipleGrids
!               getRMSD_dp                      interactMultipleGrids
!
!               rmsd                            ls_rmsd_original
!               matmul                          INTRINSIC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE      SUBROUTINE                      FMT
!
!               gridpath2//                     DAT             getRMSD_dp              FMT7,FMT3
!                 trim(subcell).dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module interactMultipleGrids
use ANALYSIS
!use VARIABLES
use PARAMETERS
implicit none

!For the RKHS
real(dp),allocatable :: RKHScoordsbuffer1(:,:,:)
real(dp),allocatable :: RKHSgradientbuffer1(:,:,:)

!Whereas the buffer1 (the array holding information on
!cells of order = 1) may be large

real(dp),allocatable :: valsbuffer1(:,:)
real(dp),allocatable :: coordsbuffer1(:,:,:)
real(dp),allocatable :: gradientbuffer1(:,:,:)
real(dp),allocatable :: potEbuffer1(:)

real(dp),allocatable :: Ubuffer1(:,:,:)
!real(dp),allocatable :: RMSDbuffer1(:)
!real(dp),allocatable :: CMdiffbuffer1(:)

integer,allocatable :: Ntrajbuffer1(:)
integer,allocatable :: varindexbuffer1(:)

!For the memory buffer

logical :: memory_available
integer,dimension(Nvar) :: previous_var_index
integer,allocatable :: vals_hash(:,:)

!integer,allocatable :: populationbuffer2(:)
real(dp),allocatable :: valsbuffer2(:,:,:)
integer,allocatable :: Ntrajbuffer2(:,:)
real(dp),allocatable :: coordsbuffer2(:,:,:,:)
real(dp),allocatable :: gradientbuffer2(:,:,:,:)
real(dp),allocatable :: potEbuffer2(:,:)

integer,allocatable :: temppopulationbuffer2(:)
real(dp),allocatable :: tempvalsbuffer2(:,:,:)
integer,allocatable :: tempNtrajbuffer2(:,:)
real(dp),allocatable :: tempcoordsbuffer2(:,:,:,:)
real(dp),allocatable :: tempgradientbuffer2(:,:,:,:)
real(dp),allocatable :: temppotEbuffer2(:,:)

!Arrays for interpolation

logical, allocatable :: acceptable_frame_mask(:)
real(dp),allocatable :: inputCLS(:,:)

!In the latter case, we need to keep track of how large
!the buffer gets, and increase it if it gets too large

integer :: buffer1_size
integer :: buffer2_size

!It is useful to estimate how many cells we may be
!looking at at one time

!integer,dimension(Norder_max+1) :: subcellsearch_max1 = (/ 0, 3 /)
!integer,dimension(Norder_max+1) :: subcellsearch_max2 = (/ 0, 3 /)

!integer,dimension(Norder_max+1) :: subcellsearch_max = (/ 0, 3 /)
integer,dimension(Norder_max+1) :: local_frame_count

!Assuming that each subcellsearch_max < 5
!This only approximates the number of cells in the outer
!shell of a Nvar-dimensional cube (not diamond, like
!used in the program)
!integer,parameter :: number_of_cells_max = &
!        2 * ((1 + 2 * (5))**(Nvar) - (1 + 2 * (5 - 1))**(Nvar))
!integer,parameter :: number_of_frames_max = &
!        number_of_cells_max * var_overcrowd(2)

!Store the index of the best candidate in approximation_index

integer :: number_of_cells
integer :: number_of_memory_cells
integer,allocatable :: approximation_index(:)

!Variables related to interpolation

integer :: Ninterpolation

!Other global variables to clean things up

real(dp),dimension(3,Natoms) :: candidate_gradient

integer :: Totalnumber_of_frames
integer :: Norder
integer :: var_filechannel
real(dp) :: candidate_rmsd

real(dp),dimension(25) :: qtest, pdottest
real(dp),dimension(25) :: qhold, pdothold


contains


subroutine checkState_PCM(vals,coords,gradient,&
                number_of_frames,order,neighbor_check)
use ls_rmsd_original
use ANALYSIS
!use VARIABLES
use PARAMETERS
use FUNCTIONS
use SIMILARITY
implicit none

integer :: i,j,k,l

!
! for GNU fortran
integer,intent(out),optional :: order,number_of_frames,neighbor_check
! for intel fortran
!integer*8,intent(out),optional :: order,number_of_frames,neighbor_check
!
!
integer :: population
integer :: chosen_index
integer :: OMP_GET_THREAD_NUM
logical :: stop_flag
real(dp), dimension(Nvar), intent(in) :: vals
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), dimension(3,Natoms) :: new_coords
real(dp), dimension(3,Natoms), intent(out) :: gradient
!real(dp), intent(inout) :: min_rmsd
real(dp) :: min_rmsd
!integer, dimension(1+Ngrid_total),intent(in) :: filechannels
character(Ngrid_text_length) :: Ngrid_text
character(5) :: variable_length_text
real(dp), dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error
real(dp),dimension(Natoms,Natoms) :: temp_CM

!For VENUS
min_rmsd = default_SIs(Nsort)

Ntrajbuffer1 = 0
var_filechannel = filechannels(1)
!var_coords = coords
!call getCoulombMatrix(Natoms,var_coords,charges,var_CM)

call setTarget(coords)

!For VENUS
interpolated_SIs = default_SIs

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
order = 0
number_of_cells = 0
neighbor_check = 0
Ninterpolation = 0
stop_flag = .false.

number_of_memory_cells = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrieve the index of each variable with
!respect to the grid and the real number
!(rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

call shiftMemory((var_index(:,Norder_order(1)+1) -&
                 previous_var_index))
previous_var_index = var_index(:,Norder_order(1)+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (usePreviousFrameInTrajectory_flag) then
if (.not.(usePreviousFramesExactlyOnce_flag)) then
do i = 1, Nhistory
  call addFrameToBufferNotInLibrary_PCM(&
          reshape(qhistory(1:3*Natoms,i),(/3,Natoms/)),&
          reshape(-pdothistory(1:3*Natoms,i)*0.529177249d0/(627.5095d0*0.04184d0) ,(/3,Natoms/)),&
          Vprev)
end do

if ((historyMayBeIMRR_flag).and.(dontUseRecentHistoricalIMRR_flag)) then
  if (Naccept > 0) call eraseBuffer(1,Naccept)
  Ninterpolation = Ninterpolation - Naccept

  ! This assumes the frames are in a specific order...
  !   namely that Naccept_max IMRR points are
  !   always accepted
  ! This may not work because the frames are sorted
  !   (because we did not use the 'force' method)
  !   and now they may be in a different order
  if (NremoveEvenMoreHistoricalIMRR > 0) then
    i = min(1+NremoveEvenMoreHistoricalIMRR,Ninterpolation)
    if (i > 1) call eraseBuffer(2,i)
    Ninterpolation = Ninterpolation + 1 - i
  end if
end if
end if
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Now, we start iterating over the grids
do k = 1, Ngrid_total
    
    !Streamline the process by storing the
    !path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    !The way we check, we check by order first
    do l = 1, Norder_max+1
    
        !The user can specify in what order to
        !check the orders through this array
        Norder = Norder_order(l)
        
        !Read the frames in the cells and
        !process their RMSDs

        !If 
        if (k == grid_addition) then
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_PCM(&
                        var_index(:,Norder+1),&
                        population)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if
            local_frame_count(Norder+1) = &
                    population
        else
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_PCM(&
                        var_index(:,Norder+1),&
                        population)
            else
            end if
        end if
        
        if (population >= var_overcrowd(Norder+1)) then
            if (Norder < Norder_max) cycle
        else
            Totalnumber_of_frames = &
                Totalnumber_of_frames + population
        end if
        
        !If the cell is unpopulated or a certain flag is true,
        !then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            !Integer i keeps track of how far away from the original
            !subcell we are; we look at cells on the 'diamond' surrounding
            !the original subcell
            do i = 1, subcellsearch_max(Norder+1)
        
                var_index_diff = 0
        
                call getRelativeIndex_PCM(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if
        

        !If we found a non-empty cell then our search has
        !been over this particular order and does not need
        !to go over other orders
        if (Totalnumber_of_frames > 0) then
                order = order + Norder
                exit
        end if
    
    end do

    !Record the RMSD encountered here for later analysis:
    !particularly, percent-RMSD graphs
    trajRMSDbuffer(k,Naccept+1) = SIbuffer1(Nsort,1)
    
end do

number_of_frames = Totalnumber_of_frames
neighbor_check = number_of_cells
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (usePreviousFrameInTrajectory_flag) then
if (usePreviousFramesExactlyOnce_flag) then
if ((historyMayBeIMRR_flag).and.(dontUseRecentHistoricalIMRR_flag)) then
  do i = Nhistory, Naccept+1, -1
    call addFrameToBufferNotInLibrary_force_PCM(&
            reshape(qhistory(1:3*Natoms,i),(/3,Natoms/)),&
            reshape(-pdothistory(1:3*Natoms,i)*0.529177249d0/(627.5095d0*0.04184d0) ,(/3,Natoms/)),&
            Vprev)
  end do

  ! This assumes the frames are in a specific order...
  !   namely that Naccept_max IMRR points are
  !   always accepted
  if (NremoveEvenMoreHistoricalIMRR > 0) then
    i = min(1+NremoveEvenMoreHistoricalIMRR,Ninterpolation)
    if (i > 1) call eraseBuffer(2,i)
    Ninterpolation = Ninterpolation + 1 - i
  end if
else
  do i = Nhistory, 1, -1
    call addFrameToBufferNotInLibrary_force_PCM(&
            reshape(qhistory(1:3*Natoms,i),(/3,Natoms/)),&
            reshape(-pdothistory(1:3*Natoms,i)*0.529177249d0/(627.5095d0*0.04184d0) ,(/3,Natoms/)),&
            Vprev)
  end do
end if
end if
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (Ninterpolation == 0) then
    candidate_rmsd = min_rmsd

    if (MLmethodID == 2) then
        call getRKHSGradient(new_coords,gradient)
        gradient = gradient * 0.529177249d0 / ( 627.5095d0 * 0.04184d0 )
        candidate_gradient = gradient
    end if

    return
end if

if ((interpolation_flag).and.(Ninterpolation > 1)) then
    select case(MLmethodID)
    case(1)
        call getInterpolatedGradient(new_coords,gradient)
    case(2)
        call getRKHSGradient(new_coords,gradient)
        gradient = gradient * 0.529177249d0 / ( 627.5095d0 * 0.04184d0 )
    end select

!   call getSIs(coords+new_coords,&
!               new_coords,new_U,new_SI)

    min_SIs = SIbuffer1(:,1)
    interpolated_SIs = sqrt(sum(new_coords(1:3,1:Natoms)**2)/Natoms)

else
    select case(MLmethodID)
    case(1)
        gradient = matmul(Ubuffer1(:,:,1),&
                          gradientbuffer1(:,:,1))
    case(2)
        call getRKHSGradient(new_coords,gradient)
        gradient = gradient * 0.529177249d0 / ( 627.5095d0 * 0.04184d0 )
    end select

!   largest_weighted_rmsd = min_rmsd
!   largest_weighted_rmsd2 = min_rmsd**2

!   interpolated_CMdiff = SIbuffer1(2,1)

    min_SIs = SIbuffer1(:,1)
    interpolated_SIs = SIbuffer1(:,1)
    largest_weighted_SIs = SIbuffer1(:,1)
    largest_weighted_SIs2 = SIbuffer1(:,1)**2

    print *, "RMSD:", SIbuffer1(Nsort,1), " vals:", &
                      valsbuffer1(:,1), " W:", 1.0d0
    
end if

min_rmsd = interpolated_SIs(Nsort)

candidate_rmsd = min_SIs(Nsort)
candidate_gradient = &
        matmul(Ubuffer1(:,:,1),&
        gradientbuffer1(:,:,1))

return

end subroutine checkState_PCM


recursive subroutine getRelativeIndex_PCM(currentVar,var_index,&
                                      N,var_index_diff,stop_flag)
use ANALYSIS
use PARAMETERS
implicit none

integer,intent(in) :: currentVar
integer,dimension(Nvar),intent(in) :: var_index
integer,intent(in) :: N

integer :: nextVar
integer :: j, k
integer :: population

integer,dimension(Nvar),intent(inout) :: var_index_diff
logical,intent(inout) :: stop_flag

if (currentVar == 1) then
!   stop_flag = .false.
    var_index_diff = 0
    k = 0
else
    k = sum(abs(var_index_diff(1:currentVar-1)))
end if

if (k == N) then
    do nextVar = currentVar, Nvar
        var_index_diff(nextVar) = 0
    end do

    call getRMSD_1_PCM(var_index + var_index_diff,population)

    if (population > 0) then
        Totalnumber_of_frames = &
        Totalnumber_of_frames + population

        stop_flag = .true.
    end if

else if (currentVar == Nvar) then

    var_index_diff(currentVar) = N - k

    call getRMSD_1_PCM(var_index + var_index_diff,population)

    if (population > 0) then
        Totalnumber_of_frames = &
        Totalnumber_of_frames + population

        stop_flag = .true.
    end if

    var_index_diff(currentVar) = -N + k

    call getRMSD_1_PCM(var_index + var_index_diff,population)

    if (population > 0) then
        Totalnumber_of_frames = &
        Totalnumber_of_frames + population

        stop_flag = .true.
    end if

else
    do j = -N + k, N - k
        var_index_diff(currentVar) = j
        nextVar = currentVar + 1
        call getRelativeIndex_PCM(nextVar,var_index,&
                              N,var_index_diff,stop_flag)
    end do
end if

end subroutine getRelativeIndex_PCM



subroutine checkNoState_OnlyHistory(coords,gradient)
use ls_rmsd_original
use ANALYSIS
!use VARIABLES
use PARAMETERS
use FUNCTIONS
use SIMILARITY
implicit none

integer :: i,j,k,l

!
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), dimension(3,Natoms) :: new_coords
real(dp), dimension(3,Natoms), intent(out) :: gradient
real(dp) :: min_rmsd
real(dp), dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

!For VENUS
min_rmsd = default_SIs(Nsort)

call setTarget(coords)

!For VENUS
interpolated_SIs = default_SIs

Ninterpolation = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((historyMayBeIMRR_flag).and.(dontUseRecentHistoricalIMRR_flag)) then
  do i = Nhistory, Naccept+1, -1
    call addFrameToBufferNotInLibrary_force_PCM(&
            reshape(qhistory(1:3*Natoms,i),(/3,Natoms/)),&
            reshape(-pdothistory(1:3*Natoms,i)*0.529177249d0/(627.5095d0*0.04184d0) ,(/3,Natoms/)),&
            Vprev)
  end do

  ! This assumes the frames are in a specific order...
  !   namely that Naccept_max IMRR points are
  !   always accepted
  if (NremoveEvenMoreHistoricalIMRR > 0) then
    i = min(1+NremoveEvenMoreHistoricalIMRR,Ninterpolation)
    if (i > 1) call eraseBuffer(2,i)
    Ninterpolation = Ninterpolation + 1 - i
  end if
else
  do i = Nhistory, 1, -1
    call addFrameToBufferNotInLibrary_force_PCM(&
            reshape(qhistory(1:3*Natoms,i),(/3,Natoms/)),&
            reshape(-pdothistory(1:3*Natoms,i)*0.529177249d0/(627.5095d0*0.04184d0) ,(/3,Natoms/)),&
            Vprev)
  end do
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (Ninterpolation == 0) then
    candidate_rmsd = min_rmsd
    return
end if

if ((interpolation_flag).and.(Ninterpolation > 1)) then
    select case(MLmethodID)
    case(1)
        call getInterpolatedGradient(new_coords,gradient)
    case(2)
        call getRKHSGradient(new_coords,gradient)
        gradient = gradient * 0.529177249d0 / ( 627.5095d0 * 0.04184d0 )
    end select

!   call getSIs(coords+new_coords,&
!               new_coords,new_U,new_SI)

    min_SIs = SIbuffer1(:,1)
    interpolated_SIs = sqrt(sum(new_coords(1:3,1:Natoms)**2)/Natoms)

else
    select case(MLmethodID)
    case(1)
        gradient = matmul(Ubuffer1(:,:,1),&
                          gradientbuffer1(:,:,1))
    case(2)
        call getRKHSGradient(new_coords,gradient)
        gradient = gradient * 0.529177249d0 / ( 627.5095d0 * 0.04184d0 )
    end select

!   largest_weighted_rmsd = min_rmsd
!   largest_weighted_rmsd2 = min_rmsd**2

!   interpolated_CMdiff = SIbuffer1(2,1)

    min_SIs = SIbuffer1(:,1)
    interpolated_SIs = SIbuffer1(:,1)
    largest_weighted_SIs = SIbuffer1(:,1)
    largest_weighted_SIs2 = SIbuffer1(:,1)**2

    print *, "RMSD:", SIbuffer1(Nsort,1), " vals:", &
                      valsbuffer1(:,1), " W:", 1.0d0
    
end if

min_rmsd = interpolated_SIs(Nsort)

candidate_rmsd = min_SIs(Nsort)
candidate_gradient = &
        matmul(Ubuffer1(:,:,1),&
        gradientbuffer1(:,:,1))

return

end subroutine checkNoState_OnlyHistory




subroutine getInterpolatedGradient(weighted_coords,weighted_gradient)
use ls_rmsd_original
use FUNCTIONS
use ANALYSIS
use PARAMETERS
use SIMILARITY
implicit none

!Information of the frame is held here
real(dp), dimension(3,Natoms), intent(out) :: weighted_coords
real(dp), dimension(3,Natoms), intent(out) :: weighted_gradient
real(dp) :: current_rmsd

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)

real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp),allocatable :: minimized_differences2(:,:)

real(dp),allocatable :: temp_valsbuffer1(:,:)
real(dp),allocatable :: temp_coordsbuffer1(:,:,:)
real(dp),allocatable :: temp_gradientbuffer1(:,:,:)
real(dp),allocatable :: temp_Ubuffer1(:,:,:)
integer,allocatable :: temp_approximation_index(:)
logical ,allocatable :: temp_acceptable_frame_mask(:)
real(dp),allocatable :: temp_inputCLS(:,:)

integer,allocatable :: RMSD_indexes(:)

integer :: Ninterpolation_true
real(dp) :: weight
real(dp) :: total_weight

!Incremental integers
integer :: i, j, k

integer :: int1, int2
real :: real1, real2

Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)

print *, "interpolation start"

allocate(frame_weights(Ninterpolation),&
         outputCLS(3*NATOMS+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         minimized_differences2(Ninterpolation,1))

do i = 1, Ninterpolation
    inputCLS(3*NATOMS+i,:) = 0.0d0
    inputCLS(3*NATOMS+i,i) = alpha_ratio * &
        SIbuffer1(Nsort,i)**2
end do

restraints = 1.0d0
restraint_values = 1.0d0
outputCLS(1:3*NATOMS) = 0.0d0
outputCLS(3*NATOMS+1:3*NATOMS+Ninterpolation) = 0.0d0

call CLS2(inputCLS(1:3*NATOMS+Ninterpolation,&
                   1:Ninterpolation),&
          3*NATOMS+Ninterpolation,&
          Ninterpolation,&
          restraints,1,restraint_values,&
          outputCLS,frame_weights)

if (all(frame_weights == 0.0d0)) then
    do i = 1, Ninterpolation
        outputCLS(i) = inputCLS(3*NATOMS+i,i)
    end do
    frame_weights(minloc(outputCLS(1:Ninterpolation))) = 1.0d0
end if

weighted_coords = 0.0d0
weighted_gradient = 0.0d0
total_weight = 0.0d0

largest_SIs = 0.0d0
largest_weighted_SIs = 0.0d0
largest_weighted_SIs2 = 0.0d0

do i = 1, Ninterpolation

    weight = frame_weights(i)
    total_weight = total_weight + weight

    do j = 1, NSIs
        largest_SIs(j) = max(largest_SIs(j),&
                             SIbuffer1(j,i))
        largest_weighted_SIs(j) = max(&
                largest_weighted_SIs(j),&
                abs(weight)*SIbuffer1(j,i))
        largest_weighted_SIs2(j) = &
                largest_weighted_SIs2(j) +&
                abs(weight)*((SIbuffer1(j,i)**2))
!               (weight*(SIbuffer1(j,i)**2))**2
    end do

    print *, "RMSD:", SIbuffer1(Nsort,i), " vals:", &
                      valsbuffer1(:,i), " W:", weight

    weighted_coords = weighted_coords + &
             weight * reshape(inputCLS(1:3*NATOMS,i),&
                                       (/3,Natoms/))
    
    weighted_gradient = weighted_gradient + &
             weight * matmul(Ubuffer1(:,:,i),&
                             gradientbuffer1(:,:,i))
end do

!largest_weighted_SIs2 = &
!        sqrt(largest_weighted_SIs2 / Natoms)

deallocate(frame_weights,outputCLS,&
           restraints,restraint_values,&
           minimized_differences2)

print *, "interpolation end"

return

end subroutine getInterpolatedGradient


subroutine getRKHSGradient(weighted_coords,weighted_gradient)
use ls_rmsd_original
use FUNCTIONS
use ANALYSIS
use PARAMETERS
use SIMILARITY
implicit none

!Information of the frame is held here
real(dp), dimension(3,Natoms), intent(out) :: weighted_coords
real(dp), dimension(3,Natoms), intent(out) :: weighted_gradient

real(dp), allocatable :: Omega(:,:),libUgradients(:,:)
real(dp), allocatable :: inputJacobians(:,:), Jacobian1(:,:)
real(dp), allocatable :: inputFeatures(:,:), Feature1(:)

real(dp),dimension(3,NAtoms) :: coords2
real(dp),dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

real(dp),dimension(3*NATOMS) :: outputF

integer :: Ncoords,Nfeatures
integer :: i, j, k

Ncoords = 3 * NATOMS
Nfeatures = (Natoms * (Natoms-1) ) / 2

allocate(libUgradients(Ncoords,Ninterpolation))

!allocate(inputJacobians(Nfeatures,Ncoords*Ninterpolation),&
!         Jacobian1(Nfeatures,Ncoords),&
!         inputFeatures(Nfeatures,Ninterpolation),&
!         Feature1(Nfeatures))
allocate(inputJacobians(Nfeatures,Ncoords*800),&
         Jacobian1(Nfeatures,Ncoords),&
         inputFeatures(Nfeatures,800),&
         Feature1(Nfeatures))

call getCoulombMatrixVector(reshape(&
         coords_target,(/Ncoords/)),Feature1)
call getCoulombMatrixGradient(reshape(&
         coords_target,(/Ncoords/)),Jacobian1)

!call sVVFwJ(inputFeatures(1:Nfeatures,1:Ninterpolation),&
!            inputJacobians(1:Nfeatures,1:Ncoords*Ninterpolation),&
!            Ncoords,Nfeatures,Ninterpolation,&
!            inputFeatures(1:Nfeatures,1),inputJacobians(1:Nfeatures,1:Ncoords),&
!            Omega(1:Ncoords,1:Ncoords*Ninterpolation))

if (RKHS_updateInverseKernel_counter < 1) then

! This part only looks at the inputs

if (useLocalRKHS_flag) then

RKHS_KernelSize = min(100,Ninterpolation)
do i = 1, RKHS_KernelSize
    call getSIs(coordsbuffer1(:,:,i),&
            coords2,new_U,new_SI)
!   libUgradients(:,i) = reshape(matmul(&
!           new_U,gradientbuffer1(:,:,i)),(/Ncoords/))

    RKHScoordsbuffer1(1:3,1:Natoms,i) = coords2
    RKHSgradientbuffer1(1:3,1:Natoms,i) = &
        matmul(new_U,gradientbuffer1(:,:,i))

    call getCoulombMatrixVector(&
        RKHScoordsbuffer1(1:3,1:Natoms,i),&
        inputFeatures(1:Nfeatures,i))
    call getCoulombMatrixGradient(&
        RKHScoordsbuffer1(1:3,1:Natoms,i),&
        inputJacobians(1:Nfeatures,&
            (i-1)*Ncoords+1:i*Ncoords))
end do
RKHS_updateInverseKernel_counter = 0 ! Always update (local)

call storeRKHSInverseKernel(inputFeatures(1:Nfeatures,1:RKHS_KernelSize),&
            inputJacobians(1:Nfeatures,1:Ncoords*RKHS_KernelSize),&
            Ncoords,Nfeatures,RKHS_KernelSize,&
            reshape(RKHSgradientbuffer1(1:3,1:Natoms,1:RKHS_KernelSize),(/Ncoords*RKHS_KernelSize/)))

else

!open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/training.xyz")

print *, "got here A"

! The old one:
!RKHS_KernelSize = 706
!open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/library700_uniformCV_entireRegion.xyz")

! The new one:
RKHS_KernelSize = 700
open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/newreactive-0.6mil_train.xyz")

print *, "got here B"
do i = 1, RKHS_KernelSize
    read(6666,FMT=*)
print *, "got here C", i
    do j = 1, Natoms
!       read(6666,FMT="(6(E12.5,1x))") &
        read(6666,FMT=*) &
                RKHScoordsbuffer1(1:3,j,i),&
                RKHSgradientbuffer1(1:3,j,i)
RKHScoordsbuffer1(1:3,j,i) = RKHScoordsbuffer1(1:3,j,i) * 1.8897259886d0
print *, &
                RKHScoordsbuffer1(1:3,j,i),&
                RKHSgradientbuffer1(1:3,j,i)
    end do

    call getCoulombMatrixVector(&
        RKHScoordsbuffer1(1:3,1:Natoms,i),&
        inputFeatures(1:Nfeatures,i))
    call getCoulombMatrixGradient(&
        RKHScoordsbuffer1(1:3,1:Natoms,i),&
        inputJacobians(1:Nfeatures,&
            (i-1)*Ncoords+1:i*Ncoords))
end do
close(6666)
print *, "got here D"
RKHS_updateInverseKernel_counter = 100000 ! Basically, never update (global)

call storeRKHSJacobian(inputFeatures(1:Nfeatures,1:RKHS_KernelSize),&
            inputJacobians(1:Nfeatures,1:Ncoords*RKHS_KernelSize),&
            Ncoords,Nfeatures,RKHS_KernelSize,&
            reshape(RKHSgradientbuffer1(1:3,1:Natoms,1:RKHS_KernelSize),(/Ncoords*RKHS_KernelSize/)))

!!!!!!!!!!!!!!!!!!!!

call getCoulombMatrixVector(reshape(&
         RKHScoordsbuffer1(1:3,1:Natoms,1),(/Ncoords/)),Feature1)
call getCoulombMatrixGradient(reshape(&
         RKHScoordsbuffer1(1:3,1:Natoms,1),(/Ncoords/)),Jacobian1)

!../alpha3667089.dat  ../alpha4710512.dat  ../alpha6499445.dat
!../alpha7288897.dat  ../alpha9421774.dat
!../alpha4385781.dat  ../alpha5685752.dat  ../alpha7154585.dat
!../alpha9165360.dat

open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha3667089.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha4710512.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha6499445.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha7288897.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha9421774.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha4385781.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha5685752.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha7154585.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do


open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha9165360.dat")
read(6666,FMT=*) RKHS_F(1:Ncoords*RKHS_KernelSize)
close(6666)
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

do i = 1, Natoms 
  print *, outputF(3*(i-1)+1) / RKHSgradientbuffer1(1,i,1)
  print *, outputF(3*(i-1)+2) / RKHSgradientbuffer1(2,i,1)
  print *, outputF(3*(i-1)+3) / RKHSgradientbuffer1(3,i,1)
end do

call getCoulombMatrixVector(reshape(&
         coords_target,(/Ncoords/)),Feature1)
call getCoulombMatrixGradient(reshape(&
         coords_target,(/Ncoords/)),Jacobian1)

!!!!!!!!!!!!!!!!!!!!

end if

print *, "got here E"

else
    if (Ninterpolation >= 100) then
        RKHS_updateInverseKernel_counter = &
                RKHS_updateInverseKernel_counter - 1
    end if
end if

if (useLocalRKHS_flag) then
!call getRKHSTarget(Ncoords,Nfeatures,Ninterpolation,&
call getRKHSTarget(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))
!           inputFeatures(1:Nfeatures,1),inputJacobians(1:Nfeatures,1:Ncoords),&

else
call getRKHSTarget_fromalpha(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))
end if

print *, "got here F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temporary (test)
!if (Ninterpolation > 10) then
if (.false.) then

RKHS_KernelSize = Ninterpolation
RKHS_updateInverseKernel_counter = 0 !100
call storeRKHSInverseKernel(inputFeatures(1:Nfeatures,1:RKHS_KernelSize),&
            inputJacobians(1:Nfeatures,1:Ncoords*RKHS_KernelSize),&
            Ncoords,Nfeatures,RKHS_KernelSize,&
            reshape(libUgradients(1:Ncoords,1:RKHS_KernelSize),(/Ncoords*RKHS_KernelSize/)))

call getRKHSTarget(Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            outputF(1:Ncoords))

print *, "getRKHSTarget input (training) printed!"
open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/RKHStraining.xyz")
do i = 1, Ninterpolation
print *, i
    call getSIs(coordsbuffer1(:,:,i),&
            coords2,new_U,new_SI)
    write(6666,FMT="(I3)") Natoms
    write(6666,FMT="(A,I3)") "training point", i
    k=1
    do j = 1, Natoms
        write(6666,FMT="(2(1x,3(F9.5,1x)))") &
            coords2(1:3,j), libUgradients(k:k+2,i)
        k=k+3
    end do
end do

weighted_gradient = reshape(outputF(1:Ncoords),(/3,NATOMS/))
print *, "getRKHSTarget (test) printed!"
k=1
write(6666,FMT="(I3)") Natoms
write(6666,FMT="(A,I3)") "test point RKHS output", 1
do i = 1, Natoms
    write(6666,FMT="(2(1x,3(F9.5,1x)))") &
        coords_target(1:3,i), outputF(k:k+2)
    k=k+3
end do
close(6666)

stop
end if

if (.false.) then
allocate(Omega(Ncoords,Ncoords*Ninterpolation))
call sVVFwJ(inputFeatures(1:Nfeatures,1:RKHS_KernelSize),&
            inputJacobians(1:Nfeatures,1:Ncoords*RKHS_KernelSize),&
            Ncoords,Nfeatures,RKHS_KernelSize,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),&
            Omega(1:Ncoords,1:Ncoords*RKHS_KernelSize))
weighted_gradient = reshape(matmul(Omega(1:Ncoords,1:Ncoords*RKHS_KernelSize),&
                            reshape(libUgradients(1:Ncoords,1:RKHS_KernelSize),&
                                    (/Ncoords*RKHS_KernelSize,1/))),&
                            (/3,NATOMS/))
deallocate(Omega)
print *, "sVVFwJ output:"
do i = 1, Natoms
    print *, weighted_gradient(:,i)
end do
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!weighted_coords = 0.0d0
!weighted_gradient = 0.0d0
!do i = 1, Ninterpolation
!    weighted_coords = weighted_coords + reshape(matmul(&
!            Omega(1:Ncoords,(i-1)*Ncoords+1:i*Ncoords),&
!            (inputCLS(1:Ncoords,i:i)+&
!             reshape(coords_target,(/Ncoords,1/)))      ),&
!                                                (/3,NATOMS/))
!    weighted_gradient = weighted_gradient + reshape(matmul(&
!            Omega(1:Ncoords,(i-1)*Ncoords+1:i*Ncoords),&
!            libUgradients(1:Ncoords,i:i)),(/3,NATOMS/))
!end do
!!R1 = sqrt(sum((  matmul(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),&
!!            inputMFN(1:Ncoords*Ninterpolation_current,&
!!                     1:1)) - inputMFN(1:Ncoords,1:1))**2)/Natoms)
!!R2 = sum(  matmul(abs(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)),&
!!           inputMFN(1:Ncoords*Ninterpolation_current,&
!!                    2:Ncoords*Ninterpolation_current+1))) / alpha_ratio

weighted_coords = 0.0d0
weighted_gradient = reshape(outputF(1:Ncoords),(/3,NATOMS/))

deallocate(libUgradients) !deallocate(Omega,libUgradients)
deallocate(inputJacobians,Jacobian1,&
           inputFeatures,Feature1)

return

end subroutine getRKHSGradient


subroutine addFrameToBufferNotInLibrary_PCM(coords,gradient,potE)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use SIMILARITY
implicit none

real(dp),intent(in) :: potE
real(dp),dimension(3,Natoms),intent(in) :: coords, gradient

integer :: index_switch

!Stores values temporarily
real(dp),dimension(Nvar) :: current_vals
real(dp),dimension(3,Natoms) :: temp_coords
real(dp),dimension(3,Natoms) :: current_coords,current_gradient
logical :: rejectNtraj_flag

!Outputs from the similarity identifier
real(dp),dimension(3,Natoms) :: new_coords
real(dp),dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI


!Incremental integers and iostate checking
integer :: n,i,j,k,iostate

do n = 1, Nindistinguishables

    BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)

    do j = 1, Natoms
         temp_coords(:,j) = coords(:,BOND_LABELLING_DATA(j))
    end do

    call getSIs(temp_coords,new_coords,new_U,new_SI)

    !If in the "accept best" method
    !and the RMSD is low enough:
    if ((new_SI(Nsort) < outer_threshold_SI).and.&
        (new_SI(Nsort) >= inner_threshold_SI)) then

!       !For diversity, accept best
!       rejectNtraj_flag = .false.
!       index_switch = 1
!       do
!           if (index_switch > Ninterpolation) exit
!           if (Ntrajbuffer2(population,single_index) == Ntrajbuffer1(index_switch)) then
!               if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
!                   rejectNtraj_flag = diversity_flag !.true.
!               else
!                   Ninterpolation = Ninterpolation - 1
!               end if
!   
!               exit
!           end if
!           index_switch = index_switch + 1
!       end do
!       if (rejectNtraj_flag) cycle
        index_switch = min(Ninterpolation + 1, Ninterpolation_max)

        if (accept_first) iostate = 1

        if (accept_worst) then
            do i = 1, index_switch
                if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
            end do
        else
            do i = 1, index_switch
                if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
            end do
        end if

        if (Ninterpolation < Ninterpolation_max) &
            Ninterpolation = Ninterpolation + 1

        call shiftBuffer(i+1,index_switch)

        do j = 1, Natoms
            coordsbuffer1(:,j,i) = coords(:,BOND_LABELLING_DATA(j))
            gradientbuffer1(:,j,i) = gradient(:,BOND_LABELLING_DATA(j))
        end do

        potEbuffer1(i) = potE
        valsbuffer1(:,i) = 0.00 !valsbuffer2(:,population,single_index)
        Ntrajbuffer1(i) = min(minval(Ntrajbuffer1(1:Ninterpolation)),-1)-1
        Ubuffer1(:,:,i) = new_U
        SIbuffer1(:,i) = new_SI

        ! Also put the vectorized coordinates into
        ! the constrained least squares input
        ! matrix called inputCLS, after substracting
        ! the target coordinates
        if (coordinatesID == 0) then
            inputCLS(1:3*NATOMS,i) =&
                   reshape(new_coords - &
                   coords_target,(/3*NATOMS/))
        else if (coordinatesID == 1) then
            call getCoulombMatrixJacobianProductForInputCLS(&
                coords_target,new_coords,inputCLS(1:3*NATOMS,i))
        end if
    end if

end do

return

end subroutine addFrameToBufferNotInLibrary_PCM

subroutine addFrameToBufferNotInLibrary_force_PCM(coords,gradient,potE)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use SIMILARITY
implicit none

real(dp),intent(in) :: potE
real(dp),dimension(3,Natoms),intent(in) :: coords, gradient

integer :: index_switch

!Stores values temporarily
real(dp),dimension(Nvar) :: current_vals
real(dp),dimension(3,Natoms) :: temp_coords
real(dp),dimension(3,Natoms) :: current_coords,current_gradient
logical :: rejectNtraj_flag

!Outputs from the similarity identifier
real(dp),dimension(3,Natoms) :: new_coords
real(dp),dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

real(dp),dimension(3,Natoms,Nindistinguishables) :: potentialCoords
real(dp),dimension(3,3,Nindistinguishables) :: potentialU
real(dp),dimension(NSIs,Nindistinguishables) :: potentialSIs


!Incremental integers and iostate checking
integer :: n,i,j,k,iostate

! This picks the best frame of all permutations

do n = 1, Nindistinguishables

    BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)

    do j = 1, Natoms
         temp_coords(:,j) = coords(:,BOND_LABELLING_DATA(j))
    end do

    call getSIs(temp_coords,new_coords,new_U,new_SI)

    potentialCoords(1:3,1:Natoms,n) = new_coords
    potentialU(1:3,1:3,n) = new_U
    potentialSIs(1:NSIs,n) = new_SI

end do

if (accept_worst) then
  n = maxloc(potentialSIs(Nsort,1:Nindistinguishables),1)
! n = maxloc(potentialSIs(Nsort,1:Nindistinguishables),1,&
!            mask=((potentialSIs(Nsort,1:Nindistinguishables)<outer_threshold_SI).and.&
!                  (potentialSIs(Nsort,1:Nindistinguishables)>=inner_threshold_SI)))
else
  n = minloc(potentialSIs(Nsort,1:Nindistinguishables),1)
! n = minloc(potentialSIs(Nsort,1:Nindistinguishables),1,&
!            mask=((potentialSIs(Nsort,1:Nindistinguishables)<outer_threshold_SI).and.&
!                  (potentialSIs(Nsort,1:Nindistinguishables)>=inner_threshold_SI)))
end if

new_coords(1:3,1:Natoms) = potentialCoords(1:3,1:Natoms,n)
new_U(1:3,1:3) = potentialU(1:3,1:3,n)
new_SI(1:NSIs) = potentialSIs(1:NSIs,n)

i = 1
index_switch = min(Ninterpolation + 1, Ninterpolation_max)

if (Ninterpolation < Ninterpolation_max) &
    Ninterpolation = Ninterpolation + 1

call shiftBuffer(i+1,index_switch)

BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
do j = 1, Natoms
    coordsbuffer1(:,j,i) = coords(:,BOND_LABELLING_DATA(j))
    gradientbuffer1(:,j,i) = gradient(:,BOND_LABELLING_DATA(j))
end do

potEbuffer1(i) = potE
valsbuffer1(:,i) = 0.00 !valsbuffer2(:,population,single_index)
Ntrajbuffer1(i) = min(minval(Ntrajbuffer1(1:Ninterpolation)),-1)-1
Ubuffer1(:,:,i) = new_U
SIbuffer1(:,i) = new_SI

! Also put the vectorized coordinates into
! the constrained least squares input
! matrix called inputCLS, after substracting
! the target coordinates
if (coordinatesID == 0) then
    inputCLS(1:3*NATOMS,i) =&
           reshape(new_coords - &
           coords_target,(/3*NATOMS/))
else if (coordinatesID == 1) then
    call getCoulombMatrixJacobianProductForInputCLS(&
        coords_target,new_coords,inputCLS(1:3*NATOMS,i))
end if

return

end subroutine addFrameToBufferNotInLibrary_force_PCM


subroutine getRMSD_1_PCM(var_index,population)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use SIMILARITY
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index
integer :: single_index, index_switch

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence
logical :: stop_flag

!Outputs from file reading
integer,intent(out) :: population

!Stores values temporarily
real(dp),dimension(Nvar) :: current_vals
integer :: current_Ntraj
real(dp) :: current_potE
real(dp),dimension(3,Natoms) :: temp_coords
real(dp),dimension(3,Natoms) :: current_coords,current_gradient
logical :: rejectNtraj_flag

!Outputs from the similarity identifier
real(dp),dimension(3,Natoms) :: new_coords
real(dp),dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

!Incremental integers and iostate checking
integer :: n,i,j,k,iostate
integer :: endpoint

if (Norder == Norder_order(1)) then
    call getFlattened(Nvar,&
            (var_index-previous_var_index),&
            single_index)
    population = populationbuffer2(single_index)
else
    single_index = single_index_max + 1
    population = -1
end if

if (population < 0) then

    write(var_filename,FMT=var_multipleFMT&
          (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
            var_index * multiplier(:,Norder+1)
    
    !Construct the subcell filename
    subcell = gridpath3//trim(var_filename)
    
    !See whether this subcell exists
    inquire(file=trim(subcell),exist=subcell_existence)
    
    if (.not. subcell_existence) then
        population = 0
        populationbuffer2(single_index) = 0
        return
    else
        number_of_cells = number_of_cells + 1
    end if

    !Open the file corresponding to the cell
    if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
    else
        open(var_filechannel,action="read",&
             file=trim(subcell))
    end if
    
    !Initialize a variable
    population = 0
    stop_flag = .true.
    do
    
        !Read the candidate frame
        if (unreadable_flag) then
    
            !In unformatted files, the first line are the variables
            !which do not need to be stored
            read(var_filechannel,iostat=iostate) &
                    (current_vals(i),i=1,Nvar)
 
            !If there are no more lines, stop; the population of the cell should be
            !one less the number of times that this portion of the loop was called
            if (iostate /= 0) exit

            population = population + 1

            valsbuffer2(:,population,single_index) = current_vals

            !The next line is trajectory number
            read(var_filechannel) &
                   Ntrajbuffer2(population,single_index)

            !The next line is the potential energy
            read(var_filechannel) &
                   potEbuffer2(population,single_index)

            !The next line is the coordinates
            read(var_filechannel) &
                   ((coordsbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
    
            !In unformatted files, the last line is the gradient
            read(var_filechannel) &
                    ((gradientbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
        else
    
            !In formatted files, everything (variables, coordinates, and gradient)
            !are stored in one line; FMT1 reads the variables
            read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                    (current_vals(i),i=1,Nvar)

            !If there are no more lines, stop; the population of the cell should be
            !one less the number of times that this portion of the loop was called
            if (iostate /= 0) exit

            population = population + 1

            valsbuffer2(:,population,single_index) = current_vals

            read(var_filechannel,FMT="(I5)",advance="no") &
                   Ntrajbuffer2(population,single_index)
            read(var_filechannel,FMT="(ES10.2)",advance="no") &
                   potEbuffer2(population,single_index)
    
            !In formatted files, FMT3 reads the coordinates
            read(var_filechannel,FMT=FMT3,advance="no") &
                   ((coordsbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
    
            !In formatted files, FMT3 reads the gradient as well
            read(var_filechannel,FMT=FMT3) &
                    ((gradientbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
        end if

        !For diversity, accept first
!       rejectNtraj_flag = .false.
!       do i = 1, Ninterpolation
!           if (Ntrajbuffer2(population,single_index) == Ntrajbuffer1(i)) then
!               rejectNtraj_flag = diversity_flag !.true.
!           end if
!       end do
!       if (rejectNtraj_flag) cycle

        do n = 1, Nindistinguishables
    
            BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
        
            do j = 1, Natoms
                 temp_coords(:,j) = &
                     coordsbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
            end do

            call getSIs(temp_coords,new_coords,new_U,new_SI)

            !If in the "accept best" method
            !and the RMSD is low enough:
            if ((new_SI(Nsort) < outer_threshold_SI).and.&
                (new_SI(Nsort) >= inner_threshold_SI)) then

            !For diversity, accept best
            rejectNtraj_flag = .false.
            index_switch = 1
            do
                if (index_switch > Ninterpolation) exit
                if (Ntrajbuffer2(population,single_index) == Ntrajbuffer1(index_switch)) then
                    if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
                        rejectNtraj_flag = diversity_flag !.true.
                    else
                        Ninterpolation = Ninterpolation - 1
                    end if
    
!                   index_switch = index_switch + 1
                    exit
                end if
                index_switch = index_switch + 1
            end do
            if (rejectNtraj_flag) cycle
!           index_switch = index_switch - 1
        
                if (accept_first) iostate = 1
        
                if (accept_worst) then
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
                    end do
                else
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
                    end do
                end if
        
                if (Ninterpolation < Ninterpolation_max) &
                    Ninterpolation = Ninterpolation + 1
        
                call shiftBuffer(i+1,index_switch)
        
                do j = 1, Natoms
                    coordsbuffer1(:,j,i) = &
                        coordsbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
                    gradientbuffer1(:,j,i) = &
                        gradientbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
                end do
        
                potEbuffer1(i) = potEbuffer2(population,single_index)
                valsbuffer1(:,i) = valsbuffer2(:,population,single_index)
                Ntrajbuffer1(i) = Ntrajbuffer2(population,single_index)
                Ubuffer1(:,:,i) = new_U
                SIbuffer1(:,i) = new_SI
        
                ! Also put the vectorized coordinates into
                ! the constrained least squares input
                ! matrix called inputCLS, after substracting
                ! the target coordinates
                if (coordinatesID == 0) then
                    inputCLS(1:3*NATOMS,i) =&
                           reshape(new_coords - &
                           coords_target,(/3*NATOMS/))
                else if (coordinatesID == 1) then
                    call getCoulombMatrixJacobianProductForInputCLS(&
                        coords_target,new_coords,inputCLS(1:3*NATOMS,i))
                end if
            end if
    
        end do
    
        if (population == buffer2_size) then
            stop_flag = .false.
            exit
        end if
    
    end do

    populationbuffer2(single_index) = population

    if (stop_flag) then
        close(var_filechannel)
        return
    end if

else if (population == 0) then
    return

else if (population < buffer2_size) then

    do k = 1, population

        !For diversity, accept first
!       rejectNtraj_flag = .false.
!       do i = 1, Ninterpolation
!           if (Ntrajbuffer2(k,single_index) == Ntrajbuffer1(i)) then
!               rejectNtraj_flag = diversity_flag !.true.
!           end if
!       end do
!       if (rejectNtraj_flag) cycle
    
        do n = 1, Nindistinguishables
    
            BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
        
            do j = 1, Natoms
                 temp_coords(:,j) = &
                     coordsbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
            end do

            call getSIs(temp_coords,new_coords,new_U,new_SI)

            !If in the "accept best" method
            !and the RMSD is low enough:
            if ((new_SI(Nsort) < outer_threshold_SI).and.&
                (new_SI(Nsort) >= inner_threshold_SI)) then

            !For diversity, accept best
            rejectNtraj_flag = .false.
            index_switch = 1
            do
                if (index_switch > Ninterpolation) exit
                if (Ntrajbuffer2(k,single_index) == Ntrajbuffer1(index_switch)) then
                    if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
                        rejectNtraj_flag = diversity_flag !.true.
                    else
                        Ninterpolation = Ninterpolation - 1
                    end if
    
!                   index_switch = index_switch + 1
                    exit
                end if
                index_switch = index_switch + 1
            end do
            if (rejectNtraj_flag) cycle
!           index_switch = index_switch - 1
        
                if (accept_first) iostate = 1
        
                if (accept_worst) then
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
                    end do
                else
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
                    end do
                end if
        
                if (Ninterpolation < Ninterpolation_max) &
                    Ninterpolation = Ninterpolation + 1
        
                call shiftBuffer(i+1,index_switch)
        
                do j = 1, Natoms
                    coordsbuffer1(:,j,i) = &
                        coordsbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
                    gradientbuffer1(:,j,i) = &
                        gradientbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
                end do
        
                potEbuffer1(i) = potEbuffer2(k,single_index)
                valsbuffer1(:,i) = valsbuffer2(:,k,single_index)
                Ntrajbuffer1(i) = Ntrajbuffer2(k,single_index)
                Ubuffer1(:,:,i) = new_U
                SIbuffer1(:,i) = new_SI
        
                ! Also put the vectorized coordinates into
                ! the constrained least squares input
                ! matrix called inputCLS, after substracting
                ! the target coordinates
                if (coordinatesID == 0) then
                    inputCLS(1:3*NATOMS,i) =&
                           reshape(new_coords - &
                           coords_target,(/3*NATOMS/))
                else if (coordinatesID == 1) then
                    call getCoulombMatrixJacobianProductForInputCLS(&
                        coords_target,new_coords,inputCLS(1:3*NATOMS,i))
                end if
            end if
    
        end do
    
    end do

    number_of_cells = number_of_cells + 1
else
    write(var_filename,FMT=var_multipleFMT&
          (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
            var_index * multiplier(:,Norder+1)
    
    subcell = gridpath3//trim(var_filename)

    number_of_cells = number_of_cells + 1
    population = 0

    !Open the file corresponding to the cell
    if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
    else
        open(var_filechannel,action="read",&
             file=trim(subcell))
    end if
end if

do

    !Read the candidate frame
    if (unreadable_flag) then

        !In unformatted files, the first line are the variables
        !which do not need to be stored
        read(var_filechannel,iostat=iostate) &
                (current_vals(i),i=1,Nvar)
 
        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) exit

        population = population + 1

        !The next line is trajectory number
        read(var_filechannel) &
               current_Ntraj

        !The next line is the potential energy
        read(var_filechannel) &
               current_potE

        !The next line is the coordinates
        read(var_filechannel) &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        !In unformatted files, the last line is the gradient
        read(var_filechannel) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)
    else

        !In formatted files, everything (variables, coordinates, and gradient)
        !are stored in one line; FMT1 reads the variables
        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                (current_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) exit

        population = population + 1

        read(var_filechannel,FMT="(I5)",advance="no") current_Ntraj
        read(var_filechannel,FMT="(ES10.2)",advance="no") current_potE

        !In formatted files, FMT3 reads the coordinates
        read(var_filechannel,FMT=FMT3,advance="no") &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        !In formatted files, FMT3 reads the gradient as well
        read(var_filechannel,FMT=FMT3) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)

    end if

    !For diversity, accept first
!   rejectNtraj_flag = .false.
!   do i = 1, Ninterpolation
!       if (current_Ntraj == Ntrajbuffer1(i)) then
!           rejectNtraj_flag = diversity_flag !.true.
!       end if
!   end do
!   if (rejectNtraj_flag) cycle

    do n = 1, Nindistinguishables

        BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
    
        do j = 1, Natoms
             temp_coords(:,j) = &
                 current_coords(:,BOND_LABELLING_DATA(j))
        end do
    
        call getSIs(temp_coords,new_coords,new_U,new_SI)

        !If in the "accept best" method
        !and the RMSD is low enough:
        if ((new_SI(Nsort) < outer_threshold_SI).and.&
            (new_SI(Nsort) >= inner_threshold_SI)) then

        !For diversity, accept best
        rejectNtraj_flag = .false.
        index_switch = 1
        do
            if (index_switch > Ninterpolation) exit
            if (current_Ntraj == Ntrajbuffer1(index_switch)) then
                if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
                    rejectNtraj_flag = diversity_flag !.true.
                else
                    Ninterpolation = Ninterpolation - 1
                end if

!               index_switch = index_switch + 1
                exit
            end if
            index_switch = index_switch + 1
        end do
        if (rejectNtraj_flag) cycle
!       index_switch = index_switch - 1
    
            if (accept_first) iostate = 1
    
            if (accept_worst) then
!               do i = 1, index_switch+1
                do i = 1, index_switch
                    if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
                end do
            else
!               do i = 1, index_switch+1
                do i = 1, index_switch
                    if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
                end do
            end if
    
            if (Ninterpolation < Ninterpolation_max) &
                Ninterpolation = Ninterpolation + 1
    
            call shiftBuffer(i+1,index_switch)
    
            do j = 1, Natoms
                coordsbuffer1(:,j,i) = &
                    current_coords(:,BOND_LABELLING_DATA(j))
                gradientbuffer1(:,j,i) = &
                    current_gradient(:,BOND_LABELLING_DATA(j))
            end do
    
            potEbuffer1(i) = current_potE
            valsbuffer1(:,i) = current_vals
            Ntrajbuffer1(i) = current_Ntraj
            Ubuffer1(:,:,i) = new_U
            SIbuffer1(:,i) = new_SI
    
            ! Also put the vectorized coordinates into
            ! the constrained least squares input
            ! matrix called inputCLS, after substracting
            ! the target coordinates
            if (coordinatesID == 0) then
                inputCLS(1:3*NATOMS,i) =&
                       reshape(new_coords - &
                       coords_target,(/3*NATOMS/))
            else if (coordinatesID == 1) then
                call getCoulombMatrixJacobianProductForInputCLS(&
                    coords_target,new_coords,inputCLS(1:3*NATOMS,i))
            end if
        end if

    end do

end do
close(var_filechannel)

end subroutine getRMSD_1_PCM




subroutine getRMSD_2(var_index,population)
use ANALYSIS
use PARAMETERS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence

!Outputs from file reading
integer,intent(out) :: population

real(dp),dimension(Nvar) :: dummy_vals
real(dp),dimension(3,Natoms) :: dummy_coords
integer :: dummy_Ntraj
real(dp) :: dummy_potE

!Incremental integers and iostate checking
integer :: i,j,k,iostate
integer :: endpoint

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Construct the subcell filename
subcell = gridpath3//trim(var_filename)

!See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

if (.not. subcell_existence) then
    population = 0
    return
else
    number_of_cells = number_of_cells + 1
end if

!Open the file corresponding to the cell
if (unreadable_flag) then
    open(var_filechannel,action="read",form="unformatted",&
         file=trim(subcell))
else
    open(var_filechannel,action="read",&
         file=trim(subcell))
end if

!Initialize a variable
population = 1

!Because the buffer may not be large enough to hold all the frames
!and, theoretically, we don't know how many times we need to
!increase the buffer, we need an overarching do loop that is capable
!of increasing the buffer without bounds
do
    !Read the candidate frame
    if (unreadable_flag) then

        !In unformatted files, the first line are the variables
        !which do not need to be stored
        read(var_filechannel,iostat=iostate) &
                (dummy_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) then
            population = population - 1
            exit
        end if

        !The next line is the trajectory number
        read(var_filechannel) &
               dummy_Ntraj

        !The next line is the potential energy
        read(var_filechannel) &
               dummy_potE

        !The next line is the coordinates
        read(var_filechannel) &
               ((dummy_coords(i,j),i=1,3),j=1,Natoms)

        !In unformatted files, the last line is the gradient
        read(var_filechannel) &
                ((dummy_coords(i,j),i=1,3),j=1,Natoms)
    else

        !In formatted files, everything (variables, coordinates, and gradient)
        !are stored in one line; FMT1 reads the variables
        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                (dummy_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) then
            population = population - 1
            exit
        end if

        read(var_filechannel,FMT="(I5)",advance="no") dummy_Ntraj
        read(var_filechannel,FMT="(ES10.2)",advance="no") dummy_potE

        !In formatted files, FMT3 reads the coordinates
        read(var_filechannel,FMT=FMT3,advance="no") &
               ((dummy_coords(i,j),i=1,3),j=1,Natoms)

        !In formatted files, FMT3 reads the gradient as well
        read(var_filechannel,FMT=FMT3) &
                ((dummy_coords(i,j),i=1,3),j=1,Natoms)
    end if

    !Increment the number of frames visited
    population = population + 1
end do
close(var_filechannel)

end subroutine getRMSD_2





subroutine addState_new(vals,coords,gradient,potE)
use PARAMETERS
implicit none

!Inputs for the file
real(dp),dimension(Nvar),intent(in) :: vals
real(dp),dimension(3,Natoms),intent(in) :: coords
real(dp),dimension(3,Natoms),intent(in) :: gradient
real(dp),intent(in) :: potE

!Keeps track of how many frames are in a cell
integer :: population

!An index used to uniquely identify the cell
integer,dimension(Nvar,Norder_max+1) :: var_index

!Character strings used to identify the file
character(50) :: var_filename

!Incremental integers
integer :: i,j, k,l

!print *, "additon flag 0"
!print *, coords
!print *, shape(coords)
!print *, gradient
!print *, shape(gradient)

!Retreive the index of each variable with respect to the grid
!and the real number (rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

do l = 1, Norder_max+1

    Norder = Norder_order(l)

    population = local_frame_count(Norder+1)

    if (population == var_overcrowd(Norder+1)) then

        call frameAddition(vals,coords,gradient,&
                           potE,&
                           var_index(:,Norder+1))

        valsbuffer1(:,population+1) = vals
        coordsbuffer1(:,:,population+1) = coords
        gradientbuffer1(:,:,population+1) = gradient

        headers(Norder+1) = headers(Norder+1) + 1

        if (Norder < Norder_max) then

            call divyUp(population+1)

        end if

        exit

    else if (population < var_overcrowd(Norder+1)) then

        if (population == 0) then

            if ((Norder == 1) .and. &
                (local_frame_count(1) <= &
                 var_overcrowd(1))) cycle

            Nfile = Nfile + 1

        end if

        call frameAddition(vals,coords,gradient,&
                           potE,&
                           var_index(:,Norder+1))

        exit

    else
        cycle
    end if

end do

end subroutine addState_new


subroutine frameAddition(vals,coords,gradient,&
                         potE,&
                         var_index,nolabel_flag)
use PARAMETERS
implicit none

!Inputs for file writing
real(dp),dimension(Nvar),intent(in) :: vals
real(dp),dimension(3,Natoms),intent(in) :: coords,gradient
real(dp),intent(in) :: potE
integer,dimension(Nvar),intent(in) :: var_index

!Character strings used to identify the file
character(50) :: var_filename

logical,optional :: nolabel_flag

!Incremental integers
integer :: i,j,k

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Pretty self-explanatory
if (unreadable_flag) then
    open(filechannel1,file=gridpath2//trim(var_filename),position="append",form="unformatted")
    if ((force_NoLabels).or.(present(nolabel_flag).and.(nolabel_flag))) then
        write(filechannel1) (vals(j),j=1,Nvar)
        write(filechannel1) Ntraj
        write(filechannel1) potE
        write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
    else
        write(filechannel1) (vals(j),j=1,Nvar)
        write(filechannel1) Ntraj
        write(filechannel1) potE
        write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
    end if
else
    open(filechannel1,file=gridpath2//trim(var_filename),position="append")
    if ((force_NoLabels).or.(present(nolabel_flag).and.(nolabel_flag))) then
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT="(I5)",advance="no") Ntraj
        write(filechannel1,FMT="(ES10.2)",advance="no") potE
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
    else
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT="(I5)",advance="no") Ntraj
        write(filechannel1,FMT="(ES10.2)",advance="no") potE
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
    end if
end if
close(filechannel1)

end subroutine frameAddition



subroutine divyUp(frames)
use PARAMETERS
use FUNCTIONS
implicit none

!The number of frames that are to be distributed
integer,intent(in) :: frames

!An index used to uniquely identify the cell
integer,dimension(Nvar) :: var_index

logical :: nolabel_flag = .true.

integer :: i, j

Norder = Norder + 1

do i = 1, frames

    do j = 1, Nvar
        var_index(j) = int(valsbuffer1(j,i) * divisor(j,Norder+1))
    end do

    call frameAddition(valsbuffer1(:,i),&
                       coordsbuffer1(:,:,i),&
                       gradientbuffer1(:,:,i),&
                       potEbuffer1(i),&
                       var_index,nolabel_flag)
end do

Norder = Norder - 1

end subroutine divyUp


subroutine setAllocations()
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    implicit none

    allocate(RKHScoordsbuffer1(3,Natoms,800),&
             RKHSgradientbuffer1(3,Natoms,800))

    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             potEbuffer1(buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             SIbuffer1(NSIs,buffer1_size))

    allocate(Ntrajbuffer1(buffer1_size),&
             varindexbuffer1(buffer1_size))

    allocate(acceptable_frame_mask(buffer1_size),&
             inputCLS(3*NATOMS+buffer1_size,buffer1_size))

!   if (memory_flag) then
        allocate(vals_hash(single_index_max,Nvar))

        allocate(populationbuffer2(&
                    single_index_max+1),&
                 valsbuffer2(Nvar,&
                    buffer2_size,&
                    single_index_max+1),&
                 Ntrajbuffer2(&
                    buffer2_size,&
                    single_index_max+1),&
                 coordsbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 gradientbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 potEbuffer2(&
                    buffer2_size,&
                    single_index_max+1))

        allocate(temppopulationbuffer2(&
                    single_index_max+1),&
                 tempvalsbuffer2(Nvar,&
                    buffer2_size,&
                    single_index_max+1),&
                 tempNtrajbuffer2(&
                    buffer2_size,&
                    single_index_max+1),&
                 tempcoordsbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 tempgradientbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 temppotEbuffer2(&
                    buffer2_size,&
                    single_index_max+1))
!   end if

    return

end subroutine setAllocations

subroutine unsetAllocations()
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    implicit none

    deallocate(RKHScoordsbuffer1,&
               RKHSgradientbuffer1)

    deallocate(valsbuffer1,&
               coordsbuffer1,&
               gradientbuffer1,&
               potEbuffer1,&
               Ubuffer1,&
               SIbuffer1)

    deallocate(Ntrajbuffer1,&
               varindexbuffer1)

    deallocate(acceptable_frame_mask,&
               inputCLS)

!   if (memory_flag) then
        deallocate(vals_hash)

        deallocate(populationbuffer2,&
                   valsbuffer2,&
                   Ntrajbuffer2,&
                   coordsbuffer2,&
                   gradientbuffer2,&
                   potEbuffer2)

        deallocate(temppopulationbuffer2,&
                   tempvalsbuffer2,&
                   tempNtrajbuffer2,&
                   tempcoordsbuffer2,&
                   tempgradientbuffer2,&
                   temppotEbuffer2)
!   end if

    return

end subroutine unsetAllocations

subroutine shiftBuffer(first_index,last_index)
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    use FUNCTIONS
    implicit none
    integer,intent(in) :: first_index,last_index
    integer :: j
 
!   do j = index_switch, i+1, -1
    do j = last_index, first_index, -1
        valsbuffer1(:,j) = valsbuffer1(:,j-1)
        Ntrajbuffer1(j) = Ntrajbuffer1(j-1)
        coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
        gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
        potEbuffer1(j) = potEbuffer1(j-1)
        Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
        SIbuffer1(:,j) = SIbuffer1(:,j-1)
        inputCLS(:,j) = inputCLS(:,j-1)
    end do

    return
        
end subroutine shiftBuffer

subroutine eraseBuffer(first_index,last_index)
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    use FUNCTIONS
    implicit none
    integer,intent(in) :: first_index,last_index
    integer :: delta_first_last_index
    integer :: j

    delta_first_last_index = 1 + last_index - first_index
 
    do j = first_index, Ninterpolation-delta_first_last_index
        valsbuffer1(:,j) = valsbuffer1(:,j+delta_first_last_index)
        Ntrajbuffer1(j) = Ntrajbuffer1(j+delta_first_last_index)
        coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j+delta_first_last_index)
        gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j+delta_first_last_index)
        potEbuffer1(j) = potEbuffer1(j+delta_first_last_index)
        Ubuffer1(:,:,j) = Ubuffer1(:,:,j+delta_first_last_index)
        SIbuffer1(:,j) = SIbuffer1(:,j+delta_first_last_index)
        inputCLS(:,j) = inputCLS(:,j+delta_first_last_index)
    end do

    return
        
end subroutine eraseBuffer

subroutine shiftMemory(delta_var_index)
    use PARAMETERS
    use ANALYSIS
    use FUNCTIONS
    implicit none
    integer,dimension(Nvar),intent(in) :: delta_var_index
    integer :: population
    integer :: single_index,xflattened
    integer,dimension(Nvar) :: xexpanded

    !Memory buffer may be turned off
    if (.not.(memory_flag)) &
        populationbuffer2 = -1

    if (all(delta_var_index == 0)) return

    temppopulationbuffer2 = populationbuffer2
    tempvalsbuffer2 = valsbuffer2
    tempNtrajbuffer2 = Ntrajbuffer2
    tempcoordsbuffer2 = coordsbuffer2
    tempgradientbuffer2 = gradientbuffer2
    temppotEbuffer2 = potEbuffer2

    populationbuffer2 = -1

    do single_index = 1, single_index_max

        population = temppopulationbuffer2(single_index)
        if (population == -1) cycle

        call getExpanded(Nvar,single_index,xexpanded)
        call getFlattened(Nvar,xexpanded-delta_var_index,xflattened)

        if (xflattened <= single_index_max) then
            populationbuffer2(xflattened) = population
            valsbuffer2(:,1:population,xflattened) = &
                tempvalsbuffer2(:,1:population,single_index)
            Ntrajbuffer2(1:population,xflattened) = &
                tempNtrajbuffer2(1:population,single_index)
            coordsbuffer2(:,:,1:population,xflattened) = &
                tempcoordsbuffer2(:,:,1:population,single_index)
            gradientbuffer2(:,:,1:population,xflattened) = &
                tempgradientbuffer2(:,:,1:population,single_index)
            potEbuffer2(1:population,xflattened) = &
                temppotEbuffer2(1:population,single_index)
        end if

    end do

    number_of_memory_cells = &
        sum(populationbuffer2(1:single_index_max),&
            mask=(populationbuffer2 < 0))

    !Memory buffer may be turned off
    if (.not.(memory_flag)) &
        populationbuffer2 = -1

    return

end subroutine shiftMemory






subroutine initializeGG()
use PARAMETERS
use ANALYSIS
use SIMILARITY

character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

integer :: cr
integer :: n,m

! Initialize PARAMETERS

print *, "Successful grid making initialization!"

!We need to initialize how we format the text
!This is big and unruly

write(Natom_text,FMT="(I0.6)") NATOMS
write(FMT2,FMT="(A22)") "("//Natom_text//"(6(1x,F14.10)))"
write(FMT3,FMT="(A22)") "("//Natom_text//"(3(1x,F14.10)))"

var_multipleFMT = ""
do m = 1, Norder_max+1
var_multipleFMT = trim(adjustl(var_multipleFMT))//"("

do n = 1, Nvar
        var_multipleFMT = trim(adjustl(var_multipleFMT))//&
        var_singleFMT(1+singleFMT_length*(m-1):&
                        singleFMT_length*m)

        if (n == Nvar) exit

        var_multipleFMT = trim(adjustl(var_multipleFMT))//&
                ',"_",'
end do

var_multipleFMT = trim(adjustl(var_multipleFMT))//',".dat")'
end do

!We also need to initialize a few scaling factors
!to speed up computation
multiplier(:,1) = var_spacing(:)
do m = 1, Norder_max
do n = 1, Nvar
        multiplier(n,m+1) = var_spacing(n) / product(var_scaling(n,1:m),DIM=1)
end do
end do

do m = 1, Norder_max+1
do n = 1, Nvar
        divisor(n,m) = 1.0 / multiplier(n,m)
end do
end do

!Let's make the actual directory
allocate(character(len(gridpath0)+&
        Ngrid_text_length+1) :: gridpath1)
allocate(character(len(gridpath0)+&
        Ngrid_text_length+1+5) :: gridpath2)
allocate(character(len(gridpath0)+&
        Ngrid_text_length+1+5) :: gridpath3)

!Let's make the directories we plan to populate
call system("mkdir "//gridpath0)
do Ngrid = 1, Ngrid_max
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") Ngrid

    call system("mkdir "//gridpath0//&
            Ngrid_text//"/")
    call system("mkdir "//gridpath0//&
            Ngrid_text//"/grid/")
end do

! Initialize ANALYSIS

print *, "Successful grid checking initialization!"

!Initialize the permutation variables
Nindistinguishables = 2
allocate(INDISTINGUISHABLES(&
        Nindistinguishables,NATOMS),&
        BOND_LABELLING_DATA(NATOMS),&
        charges(NATOMS),masses(NATOMS))
do n = 1, NATOMS
    BOND_LABELLING_DATA(n) = n
end do
INDISTINGUISHABLES(1,:) = BOND_LABELLING_DATA

INDISTINGUISHABLES(2,:) = BOND_LABELLING_DATA
INDISTINGUISHABLES(2,3) = 5
INDISTINGUISHABLES(2,5) = 3

!Initialize the subcell search variables
call setSubcellSearchMax()

!Initialize the clock
call system_clock(count_rate=cr)
!system_clock_rate = 1.0/real(cr)

!Initialize the experiment and data folders
allocate(character(len(gridpath0)+&
        expfolder_length) :: gridpath4)
allocate(character(len(gridpath0)+&
        expfolder_length+5) :: gridpath5)
gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//"data/"
call system("mkdir "//gridpath4)
call system("mkdir "//gridpath5)

!The grid number will uniquely identify one trajectory
!Open all these files under filechannels
allocate(filechannels(1+Ngrid_max))
allocate(trajRMSDbuffer(Ngrid_max,Naccept_max+1))
filechannels(1) = 1000
do m = 1, Ngrid_max
    filechannels(1+m) = 1000 + 69 * m
end do

!Which grids are we checking?
Ngrid_total = 1
Ngrid = 1

!Which grids are we adding to?
grid_addition = 0
write(variable_length_text,FMT=FMT5_variable)&
        Ngrid_text_length
write(Ngrid_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//&
        ")") grid_addition
gridpath1 = gridpath0//Ngrid_text//"/"
gridpath2 = gridpath0//Ngrid_text//"/grid/"

if (gather_interpolation_flag) &
    open(filechannel3,file=gridpath5//interpolationfile,&
                      position="append")
!   open(filechannel3,file=gridpath5//interpolationfile)

if (.true.) &
    open(filechannel4,file=gridpath5//"energydrift.dat",&
                      position="append")
!   open(filechannel4,file=gridpath5//"energydrift.dat")

if (readtrajectory_flag) then
    call system("cp "//&
                readtrajectory_path//&
                readtrajectoryfile//" "//&
                gridpath4//readtrajectoryfile)
    open(filechannel5,file=gridpath4//readtrajectoryfile)
end if

! Initialize SIMILARITY

allocate(coords_target(3,NATOMS),CM_target(NATOMS,NATOMS))
allocate(coords_transformed(3,NATOMS,NSIs))

!Buffers for later
buffer1_size = 2 + Ninterpolation_max
buffer2_size = 50
call setAllocations()

return
end subroutine initializeGG

end module interactMultipleGrids
