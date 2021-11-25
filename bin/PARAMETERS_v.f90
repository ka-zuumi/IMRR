module PARAMETERS

integer, parameter :: Natoms = 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      FILENAMES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!All intermediate files are kept in this folder
character(5) :: intermediatefolder = "data/"
!All readtrajectory trajectories are kept in this folder
character(5) :: readtrajectoryfolder = "traj/"

!File that holds all trajectories for readTrajectory
character(20),parameter :: readtrajectoryfile = "readtrajectories.txt"
!File that will keep the trajectory folder names
character(16),parameter :: trajectories = "trajectories.txt"
!File that writes the progress of the program
character(12),parameter :: progressfile = "progress.txt"
!File that writes the progress of a trajectory
character(14),parameter :: trajectoryfile = "trajectory.xyz"
!File that writes the rmsd retrieved from checkState every frame
character(14),parameter :: checkstatefile = "checkstate.dat"
!File that writes the progress of multiple trajectories
character(16),parameter :: trajectoriesfile = "trajectories.dat"
!File that has the parameters (self)
character(14),parameter :: parametersfile = "PARAMETERS.f90"
!File for any gnuplot scripting
character(11),parameter :: gnuplotfile = "gnuplotfile"
!File for storing trajectory data from across many grids
character(23),parameter :: cumulativefile = "cumulative_trajectories"
!File for storing the first and last frames of multiple trajectories
character(13),parameter :: timeslicefile = "timeslice.dat"
!File for storing information on the grid over multiple trajectories
character(15),parameter :: informaticsfile = "informatics.dat"
!File for storing scattering angle data over multiple trajectories
character(15),parameter :: SAfile = "SA_trajectories.dat"
!File for storing initial bonding data over multiple trajectories
character(26),parameter :: initialfile = "initbonds_trajectories.dat"
!File for storing TRV energy changes
character(16),parameter :: TRVfile = "TRV_trajectories.dat"
!File for both Eenrgy Decomposition and Scattering Angle over multiple trajectories
character(22),parameter :: SATRVfile = "SATRV_trajectories.dat"
character(28),parameter :: binnedSATRVfile = "binnedSATRV_trajectories.dat"
!File for checking error files
character(14),parameter :: errorcheckfile = "errorcheck.dat"

!File to write to for system calls
character(8),parameter :: temporaryfile1 = "tmp1.txt"
character(8),parameter :: temporaryfile2 = "tmp2.txt"
character(8),parameter :: temporaryfile3 = "tmp3.txt"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     FILE CHANNELS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Channel for the progressfile
integer,parameter :: progresschannel = 70
integer,parameter :: trajectorieschannel = 71
integer,parameter :: frameschannel = 72
integer,parameter :: filechannel1 = 73
integer,parameter :: filechannel2 = 74
integer,parameter :: filechannel3 = 75
integer,parameter :: filechannel4 = 76
integer,parameter :: gnuplotchannel = 77


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If true, the grid has unformatted files; thus, keep this constant
!The processed files and check trajectory files are still formatted though
logical,parameter :: unreadable_flag = .false.

character(14),parameter :: FMT1 = "(2(1x,F14.10))"          !For two variables
character(22) :: FMT2         !For six atoms, coords and gradient
character(22) :: FMT3         !For six atoms, coords or gradient
character(4),parameter :: FMT4 = "(I4)"                     !For subcell names
character(4),parameter :: FMT5 = "(I9)"                     !For number of substates
character(7),parameter :: FMT6 = "(F12.8)"                  !For RMSD
character(19),parameter :: FMT7 = "(30x,18(1x,F14.10))"     !For skipping the variables

character(4),parameter :: FMT5_variable = "(I5)"
character(4),parameter :: FMT8_counter = "(I8)"
character(6),parameter :: FMT6_pos_real0 = "(F6.5)"
character(6),parameter :: FMT6_pos_real1 = "(F6.4)"
character(6),parameter :: FMT6_neg_real1 = "(F6.3)"
character(6),parameter :: FMT6_pos_int = "(I0.6)"

!These formats depend on the number of atoms and bonds in the system
!So they cannot be parameters and must be changed later
!These are always changed in the MAIN program
character(6) :: Nbond_text
character(19) :: FMTinitial
character(6) :: Natom_text
character(19) :: FMTtimeslice

!Formatting for files that are processed AFTER grid making and checking
character(27),parameter :: FMTinformatics = "(2(F12.7),I6,2(I5),I8,F8.4)"
character(15),parameter :: FMTsa = "((F6.4),(F8.4))"
character(10),parameter :: FMTtrv = "(3(F11.6))"
character(32),parameter :: FMTdata = "((F6.4),(F8.4),2(F11.6),(F13.9))"
character(33),parameter :: FMTnow = "('Time: ',I2.2,':',I2.2,':',I2.2)"

!Because we want to not use trim and adjustl all the time, and also
!because we want nice-looking graphs, we want to store various numbers
!in strings (with formattings as specified above)
!But to adjust to different sizes of strings, we need to establish
!the length of certain strings
integer,parameter :: trajectory_text_length = 5
integer,parameter :: Ngrid_text_length = 3

!Gridpath is one variable that must be supplied by the user and is changed
!on the LOCAL version of a library; this is changed AUTOMATICALLY
!But within the library itself it should stay constant
character(len=*),parameter :: gridpath0="HBrCO2library/"

character(len=:), allocatable :: gridpath1
character(len=:), allocatable :: gridpath2
character(len=:), allocatable :: gridpath3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The number of atoms in the system
!Natoms is another variable that must be supplied by the user and is changed
!on the LOCAL version of a library; this must be changed MANUALLY
!integer,parameter :: Natoms = 5
!integer,parameter :: Ncoords = Natoms*3

!Permutation
integer :: Nindistinguishables
!integer,dimension(8,Natoms),parameter :: INDISTINGUISHABLES
integer,allocatable :: INDISTINGUISHABLES(:,:)
!integer,dimension(Natoms) :: BOND_LABELLING_DATA
integer,allocatable :: BOND_LABELLING_DATA(:)

!The number of variables in use
!Right now, we only support a two-variable grid (easier to script and visualize)
integer,parameter :: Nvar = 2
integer,parameter :: Nvar_next = Nvar + 1

!The charges and masses of each atom
!real,dimension(Natoms),parameter :: charges = &
!        (/ 1.0e0, 35.0e0, 8.0e0, 6.0e0, 8.0e0 /)
!real,dimension(Natoms),parameter :: masses = &
!        (/ 1.008e0, 79.904e0, 15.999e0, 12.011e0, 15.999e0 /)
real,allocatable :: charges(:)
real,allocatable :: masses(:)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     GRID PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Many of the variables below are supplied by the user and are changed
!on the LOCAL version of alibrary; this is done AUTOMATICALLY

!Norder_max controls how many children generation we will make
integer,parameter :: Norder_max = 1


!Some default values are good if we don't know at execution
!how many variables we plan to use

real,parameter,dimension(4) :: default_spacing = &
        (/ 0.005, 0.005, 0.005, 0.005 /)
!       (/ 0.0025, 0.0025, 0.0025, 0.0025 /)
real,parameter,dimension(4) :: default_maxvar = &
        (/ 12.0, 12.0, 10.0, 10.0 /)
!       (/ 10.0, 12.0, 10.0, 10.0 /)
integer,parameter,dimension(4,4) :: default_scaling = &
        reshape( &
                (/  4,  4,  4,  4, &
                    10, 10, 10, 10, &
                    10, 10, 10, 10, &
                    10, 10, 10, 10    /), &
                (/ 4, 4 /))
integer,parameter,dimension(4) :: default_overcrowd = &
        (/  10000,     10000,      50,       50 /)
!       (/     50,     10000,      50,       50 /)
integer,dimension(4),parameter :: default_Norder_order = &
        (/ 0, 1, 2, 3 /)
!       (/ 1, 0, 2, 3 /)

!After making these default values, we decide how many we
!can use, then simply truncate the extra variables off

real,parameter,dimension(Nvar) :: var_spacing = &
        default_spacing(1:Nvar)
real,parameter,dimension(Nvar) :: var_maxvar = &
        default_maxvar(1:Nvar)
integer,parameter,dimension(Nvar,Norder_max+1) :: var_scaling = &
        default_scaling(1:Nvar,1:Norder_max+1)
integer,parameter,dimension(Norder_max+1) :: var_overcrowd = &
        default_overcrowd(1:Norder_max+1)
integer,dimension(Norder_max+1) :: Norder_order = &
        default_Norder_order(1:Norder_max+1)

!We also need to define some default formats for cells
!of different orders;
!The singleFMT strings are formats for each variable, which
!will be combined later to formats for each cell
!
!One crucial thing here is that each string is of the length
!specified in singleFMT_length

integer,parameter :: singleFMT_length = 6
character(singleFMT_length*4),parameter :: default_singleFMT = &
        "(F0.4)"// &
        "(F0.6)"// &
        "(F0.7)"// &
        "(F0.8)"
!       "(F0.2)"// &
!       "(F0.4)"// &
!       "(F0.5)"// &
!       "(F0.6)"

character(singleFMT_length*(Norder_max+1)),parameter :: var_singleFMT = &
        default_singleFMT(1:singleFMT_length*(Norder_max+1))

integer,parameter :: multipleFMT_length = (Nvar) * singleFMT_length + &
                                          5 * (Nvar - 1) + 7 + 2
character(multipleFMT_length*(Norder_max+1)) :: var_multipleFMT


!The "bounds", or the upper limit to the indexing of
!the variables, is given a little padding (400 in this case)

integer,parameter,dimension(Nvar) :: var_bounds = &
        (/ (ceiling(var_maxvar(increment) / var_spacing(increment)) + 400, increment=1,Nvar) /)

!The "resolution" is how many new cells are created when a cell
!of a certain order has divyUp called on it

!integer,parameter,dimension(Norder_max+1) :: var_resolution = &
!        (/ (product(var_scaling(:,increment),DIM=1), increment=1,Norder_max+1) /)

!And then some scaling factors useful to save calculation time later

real,dimension(Nvar,Norder_max+1) :: multiplier
real,dimension(Nvar,Norder_max+1) :: divisor

!And we have some internal limit to how many cells we will search out

integer,parameter,dimension(4) :: default_subcellsearch_max = &
        (/ 0, 0, 0, 0 /)
integer,dimension(Norder_max+1) :: subcellsearch_max =&
        default_subcellsearch_max(1:Norder_max+1)
integer,dimension(Norder_max+1) :: subcellsearch_max1 =&
        default_subcellsearch_max(1:Norder_max+1)
integer,dimension(Norder_max+1) :: subcellsearch_max2 =&
        default_subcellsearch_max(1:Norder_max+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    MEMORY OVERHEAD COST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If we want to add duplicate copies of a state (but permuted labels)
!Then set this to .true.
!Otherwise, all frames are relabelled to a specific format
!Note: if true, this will fill up the grid much faster so a larger
!value of overcrowd is recommended
logical,parameter :: force_Duplicates = .false.

!If force_Duplicate is set to false there is an additonal feature here
!for the frame to be added to have no relabelling
!Note: this makes it so that two indistinguishable frames
!      possibly do not have the same RMSD or even cell index
logical,parameter :: force_NoLabels = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      COUNTERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,dimension(Norder_max+1) :: headers
integer,dimension(Norder_max+1) :: headers_old
integer,dimension(Norder_max+1) :: Norder_total

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   MULTI-GRID PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Ntraj_max, Ngrid_max, and Ngrid_check_min are three variables supplied by
!the user to change the grid-making procedure; this is updated AUTOMATICALLY

!The number of trajectories we add to a grid
integer,parameter :: Ntraj_max = 4000
!The maximum amount of time (seconds) we are willing to wait for a single trajectory to finish
real,parameter :: trajectory_CPU_time_max = 600.0
!The number of grids we will make
integer :: Ngrid_max = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  GRID FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Because Ngrid is a shared variable, but it is often manipulated within runTrajectory,
!we need each thread to have its own private copy of Ngrid
!(and consequently, gridpath1 and gridpath2)

integer :: Ngrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  TRAJECTORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: Ntraj,Ntraj_allowed,Nfile

end module PARAMETERS



