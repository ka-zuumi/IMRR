!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               analyzeScatteringAngleDistributions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!		This module processes outputs of trjaectories, gains information like the
!		scattering angle and energy change decomposition, and plots them
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
!               TRAJECTORIESCHANNEL             OPEN, WRITE, CLOSE
!               FILECHANNEL1                    OPEN, WRITE/READ, CLOSE
!               FILECHANNEL2                    OPEN, WRITE/READ, CLOSE
!               FILECHANNEL3                    OPEN, WRITE/READ, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!               getScatteringAngles1            prefix_filename         intent(in),character(*)
!                                               PNGfilename             intent(in),character(*)
!
!               getConvergenceImage             lowerlimit              intent(in),real(dp)
!                                               upperlimit              intent(in),real(dp)
!                                               SATRVcolumn             intent(in),integer
!                                               SATRVname               intent(in),character(*)
!
!               getScatteringAngles2            prefix_filename         intent(in),character(*)
!                                               PNGfilename             intent(in),character(*)
!
!               getInitialImages                prefix_filename         intent(in),character(*)
!                                               PNGfilename             intent(in),character(*)
!
!               getComparedScatteringAngles     lowerlimit              intent(in),real(dp)
!                                               upperlimit              intent(in),real(dp)
!                                               imagename               intent(in),character(*)
!                                               SATRVcolumn             intent(in),integer
!                                               SATRVname               intent(in),character(*)
!
!               postProcess                     prefix_filename         intent(in),character(*)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               SYSTEM                          INTRINSIC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the initial conditions for the set of
!                 initialfile                                                   trajectories associated with this prefix
!               gridpath0//prefix_filename//    DAT                             Stores the initial and final frames of a
!                 timeslicefile                                                 set of trajectories corresponding to the
!                                                                               specified prefix
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//prefix_filename//    DAT                             Same as the SATRV file but the values are
!                 binnedSATRVfile                                               stored as integers corresponding to some
!                                                                               bin in a specific binning
!               gridpath1//Initial//#traj//     DAT                             Stores the scattering angle and energy change
!                 binnedSATRVfile                                               decomposition for all trajectories in the library;
!                                                                               its is already binned
!               gridpath0//RMSD//SATRVname      DAT                             Stores the averages, devations, and other data
!                 cumulativefile                                                associated with a particular SATRV for all
!                                                                               trajectories in the library for some sampling
!               gridpath0//Adjusted//           DAT                             Similar to cumulative file but normalized and
!                 SATRVname                                                     dependent only on the number of sets and number
!                                                                               of trajectories per set
!               gridpath0//Convergence//        PNG                             The reference distribution of the library for
!                 SATRVname                                                     some sampling and some SATRV; also listed are
!                                                                               the RMSD, KRP, and running average
!               gridpath0//PNGfilename          PNG                             Generic format for image naming
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module analyzeScatteringAngleswithMultipleGrids
use DOUBLE
implicit none

!In the check-grid program we often want to analyze distributions, and to analyze those,
!we need to know the bounds of our variables
!Here we declare variables for the minimum and maximum of each variable
!This can be over-rided by user input

real :: max_absenergychange, min_absenergychange
real :: max_relenergychange, min_relenergychange
real :: max_rotenergychange, min_rotenergychange
real :: abs_energychange, rel_energychange, rot_energychange
real :: max_TranslationalEnergy

!In the check-grid program, we also need to figure out the number and sizes of
!the bins in our distribution

real(dp) :: sizeEnergyBin,sizeAngleBin
real(dp) :: sizeAbsEnergyBin,sizeRelEnergyBin,sizeRotEnergyBin
real(dp) :: sizeDeltaEnergyBin
integer,parameter :: energyBins = 100                  !This is for the SA heat map
integer,parameter :: angleBins = 200                   !Same
integer,parameter :: scatteringangleBins = 50          !This is for the regular histogram
                                                       !  (Note: make this a divisor of angleBins)
integer,parameter :: energychangeBins = 50             !Same

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getScatteringAngles1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine does two things:
!                 1) It reads in frames from the file corresponding to the prefix provided,
!                    proceeding then to bin them and put them in a differently-name filed
!                    with the same prefix
!                 2) It creates a scattering angle heatmap of the set of trajectories
!                    corresponding to the prefix provided
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
!               PNGfilename                     CHARACTER(*)                    The suffix we use for the output PNG
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntraj                           INTEGER                         The number of trajectories in this
!                                                                               set of trajectories, or the "sampling"
!               energychangeBins                INTEGER                         The number of bins used for distributions
!                                                                               involving energy change
!               energyBins                      INTEGER                         The number of bins used for distributions
!                                                                               involving the final translational energy
!               angleBins                       INTEGER                         The number of bins used for distributions
!                                                                               involving scattering angle
!               angleratio                      INTEGER                         The ratio between the size of the scattering
!                                                                               angle distribution and heatmap;
!                                                                               The heatmap bins should be SMALLER
!               angle_energy_bins               INTEGER,                        The array storing binned scattering angle
!                                               DIM(angleBins,energyBins)       data
!
!               sizeAbsEnergyBin                REAL(DP)                        The size of a bin for the distribution
!                                                                               of absolute translational energy change
!               sizeRelEnergyBin                REAL(DP)                        The size of a bin for the distribution
!                                                                               of relative translational energy change
!               sizeRotEnergyBin                REAL(DP)                        The size of a bin for the distribution
!                                                                               of rotational energy change
!               sizeDeltaEnergyBin              REAL(DP)                        The size of a bin for the distribution
!                                                                               of all kinetic energy changes
!
!               sizeEnergyBin                   REAL(DP)                        The size of a bin for the distribution
!                                                                               of final translational energy
!               sizeAngleBin                    REAL(DP)                        The size of a bin for the distribution
!                                                                               of scattering angles
!               angleslice                      REAL(DP)                        The dividing factor that slices the
!                                                                               scattering angle heatmap into a more
!                                                                               visually pleasing figure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//prefix_filename//    DAT                             Same as the SATRV file but the values are
!                 binnedSATRVfile                                               stored as integers corresponding to some
!                                                                               bin in a specific binning
!               gridpath0//prefix_filename//    PNG                             Scattering angle and energy change decompositon
!                 PNGfilename                                                   distributions for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//prefix_filename//    PNG                             Scattering angle heatmap for a set of
!                 Heatmap_//PNGfilename                                         trajectories corresponding to the specified
!                                                                               prefix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getScatteringAngles1(prefix_filename,PNGfilename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI

double precision :: c1,c2,c3,c4,c5,c6,c7,pi,halfpi,twopi

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(150) :: old_filename

!SCATTERING ANGLES
real(dp),dimension(3,3) :: coords_initial,velocities_initial,coords_final,velocities_final
real(dp) :: ScatteringAngle_real, TranslationalEnergy_real
integer :: ScatteringAngle, TranslationalEnergy
integer :: AbsEnergyChange, RelEnergyChange, RotEnergyChange

!HEAT MAP VARIABLES
integer,allocatable :: angle_energy_bins(:,:)
integer :: occurence_max
real :: bin_width
integer :: angle_ratio = angleBins / scatteringangleBins
integer :: angle_slice = 6

!I/O HANDLING
integer :: iostate

!Incremental Integers
integer :: i, j, k

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

!The plots are named starting with
!Ntraj (the number of trajectories)
write(variable_length_text,FMT=FMT5_variable)&
        trajectory_text_length
write(Ntraj_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        Ntraj

allocate(angle_energy_bins(angleBins,energyBins))
angle_energy_bins = 0

!The size of the bins, of course,
!depends on how many bins we have as
!well as the bounds of the distribution
sizeAbsEnergyBin = (max_absenergychange - &
                    min_absenergychange) /&
                    energychangeBins
sizeRelEnergyBin = (max_relenergychange - &
                    min_relenergychange) /&
                    energychangeBins
sizeRotEnergyBin = (max_rotenergychange - &
                    min_rotenergychange) /&
                    energychangeBins

sizeDeltaEnergyBin = (max(max_absenergychange,&
                          max_relenergychange,&
                          max_rotenergychange) - &
                      min(min_absenergychange,&
                          min_relenergychange,&
                          min_rotenergychange))/&
                          energychangeBins
sizeEnergyBin = max_TranslationalEnergy / &
        (energyBins)
!sizeAngleBin = pi / (angleBins)
!for VENUS:
sizeAngleBin = 3.14159265359d0 / (angleBins)

!If we are doing this for a comparison,
!we may have user-defined bounds
if (comparison_flag .and. &
        (comparison_upperlimit /= &
         comparison_lowerlimit)) then

if (trim(adjustl(comparison_SATRVname))&
        == "ScatteringAngle") then
    sizeAngleBin = (comparison_upperlimit - &
                    comparison_lowerlimit) /&
                    angleBins
else if (trim(adjustl(comparison_SATRVname))&
        == "RelativeEnergyChange") then
    sizeRelEnergyBin = (comparison_upperlimit - &
                        comparison_lowerlimit) /&
                        energychangeBins
else if (trim(adjustl(comparison_SATRVname))&
        == "AbsoluteEnergyChange") then
    sizeAbsEnergyBin = (comparison_upperlimit - &
                        comparison_lowerlimit) /&
                        energychangeBins
else if (trim(adjustl(comparison_SATRVname))&
        == "RotationalEnergyChange") then
    sizeRotEnergyBin = (comparison_upperlimit - &
                        comparison_lowerlimit) /&
                        energychangeBins
else
end if

end if

!Now we can actually bin them
!Here we open up the SATRVfile (which has
!the observables) and the binnedSATRVfile
!(which has the observables after binning)
open(filechannel1,file=gridpath5//&
        prefix_filename//SATRVfile)
open(filechannel2,file=gridpath5//&
        prefix_filename//binnedSATRVfile)
do i = 1, Ntraj

    !First, we read the line
    read(filechannel1,FMT=FMTdata,&
            iostat=iostate)&
            ScatteringAngle_real,&
            TranslationalEnergy_real,&
            abs_energychange,&
            rel_energychange,&
            rot_energychange

    !Then we figure out what bin it should be in
    ScatteringAngle = ceiling(ScatteringAngle_real/&
            sizeAngleBin)
    TranslationalEnergy = ceiling(TranslationalEnergy_real/&
            sizeEnergyBin)
    AbsEnergyChange = ceiling((abs_energychange)/&
            sizeAbsEnergyBin)
    RelEnergyChange = ceiling((rel_energychange)/&
            sizeRelEnergyBin)
    RotEnergyChange = ceiling((rot_energychange)/&
            sizeRotEnergyBin)

    !If the bin is out of bounds,
    !put it back in bounds
    if (ScatteringAngle > angleBins)&
            ScatteringAngle = angleBins
    if (ScatteringAngle == 0)&
            ScatteringAngle = 1

    if (TranslationalEnergy > energyBins)&
            TranslationalEnergy = energyBins
    if (TranslationalEnergy == 0)&
            TranslationalEnergy = 1

    if (AbsEnergyChange > energychangeBins)&
            AbsEnergyChange = energychangeBins
    if (AbsEnergyChange == 0)&
            AbsEnergyChange = 1

    if (RelEnergyChange > energychangeBins)&
            RelEnergyChange = energychangeBins
    if (RelEnergyChange == 0)&
            RelEnergyChange = 1

    if (RotEnergyChange > energychangeBins)&
            RotEnergyChange = energychangeBins
    if (RotEnergyChange == 0)&
            RotEnergyChange = 1

    !For the scattering angle heatmap, we
    !also store these values in another array
    angle_energy_bins(ScatteringAngle,&
                      TranslationalEnergy) = &
        angle_energy_bins(ScatteringAngle,&
                          TranslationalEnergy) + 1

    !We have two separate bin sizes for the
    !scattering angle heatmap and the scattering
    !angle distribution. The former have SMALLER
    !bins. Thus, we can get the binning of the
    !latter by dividing by some factor which I
    !I call "angle_ratio": the ratio between 
    !the sizes of these two binnings
    write(filechannel2,FMT=*)&
            ceiling((ScatteringAngle-0.5)/angle_ratio),&
            TranslationalEnergy,&
            AbsEnergyChange,&
            RelEnergyChange,&
            RotEnergyChange
end do
close(filechannel1)
close(filechannel2)


!This is the gnuplot code to make the scattering
!angle and energy change (SATRV) plots. This is
!only plotted if there is no comparison active

if (comparison_flag) then
    deallocate(angle_energy_bins)
    return
end if

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo size 1200,1200'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,FMT="(A)") 'set multiplot layout 4,1 margins 0.10,0.95,.1,.95 spacing 0,0.1'//&
                        'title "Scattering Angle and '//&
                        'Energy Change Distributions of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",18" offset 0,3'
write(gnuplotchannel,FMT="(A)") 'set style fill solid 1.0 noborder'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Occurence"'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'unset key'

write(gnuplotchannel,FMT="(A)") 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,FMT="(A)") 'pi = 3.14159265'
write(gnuplotchannel,FMT=*) 'box_width = ', sizeAngleBin*angle_ratio
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") 'set autoscale y'
write(gnuplotchannel,FMT="(A)") 'set xtics pi/2'
write(gnuplotchannel,FMT="(A)") "set format x '%.1P π'"
write(gnuplotchannel,FMT="(A)") 'set xrange [0:pi]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($1-0.5)):(1.0) smooth frequency with boxes'

write(gnuplotchannel,FMT="(A)") 'scaling = 1'

!write(gnuplotchannel,FMT="(A)") 'set xlabel "Absolute Translation Energy Change (meV)"'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Absolute Translation Energy Change (meV)"'
write(gnuplotchannel,FMT=*) 'min_E = scaling * ', min_absenergychange
write(gnuplotchannel,FMT=*) 'max_E = scaling * ', max_absenergychange
write(gnuplotchannel,FMT=*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,FMT="(A)") 'set xrange [min_E:max_E]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") "set format x '%.3f'"
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($3-0.5)+min_E):(1.0) smooth frequency w boxes'

!write(gnuplotchannel,FMT="(A)") 'set xlabel "Relative Translational Energy Change (meV)"'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Relative Translational Energy Change (kcal/mol)"'
write(gnuplotchannel,FMT=*) 'min_E = scaling * ', min_relenergychange
write(gnuplotchannel,FMT=*) 'max_E = scaling * ', max_relenergychange
write(gnuplotchannel,FMT=*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,FMT="(A)") 'set xrange [min_E:max_E]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($4-0.5)+min_E):(1.0) smooth frequency w boxes'

!write(gnuplotchannel,FMT="(A)") 'set xlabel "Rotational Energy Change (meV)"'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Rotational Energy Change (kcal/mol)"'
write(gnuplotchannel,FMT=*) 'min_E = scaling * ', min_rotenergychange
write(gnuplotchannel,FMT=*) 'max_E = scaling * ', max_rotenergychange
write(gnuplotchannel,FMT=*) 'Nbins = ', energychangeBins
write(gnuplotchannel,FMT="(A)") 'box_width = (max_E-min_E) / Nbins'
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") 'set xrange [min_E:max_E]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($5-0.5)+min_E):(1.0) smooth frequency w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

!For the scattering angle - final translational energy heatmap, all we need
!to do is write the values of angle_energy_bins onto a file and plot it
!For convenience, I also record the maximum frequency in the array
occurence_max = maxval(angle_energy_bins)
open(filechannel1,file=gridpath5//temporaryfile1)
do i = 1, angleBins/angle_slice+1
    do j = 1, energyBins
        write(filechannel1,FMT=*) (i-1)*sizeAngleBin,&
                (j-1)*sizeEnergyBin, angle_energy_bins(angleBins-i,j)
!               (j-1)*sizeEnergyBin, angle_energy_bins(i,j)
    end do
    write(filechannel1,FMT=*) ""
end do
close(filechannel1)
deallocate(angle_energy_bins)

!This is the gnuplot code to make the scattering angle - final translational energy heatmap
!This took A LOT of fine-tuning to make it look nice, so please don't change it
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo size 1200,1200'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//'HeatMap_'//&
                        PNGfilename//'.png"'
write(gnuplotchannel,FMT="(A)") 'set title "Scattering Angle Distribution" font ",18" offset 0,-20'
write(gnuplotchannel,FMT="(A)") 'set lmargin at screen 0.05'
write(gnuplotchannel,FMT="(A)") 'set rmargin at screen 0.85'
write(gnuplotchannel,FMT="(A)") 'set bmargin at screen 0.1'
write(gnuplotchannel,FMT="(A)") 'set pm3d map'
write(gnuplotchannel,FMT="(A)") 'set pm3d corners2color c1'
write(gnuplotchannel,FMT="(A)") 'unset key'
write(gnuplotchannel,FMT="(A)") 'set multiplot'
write(gnuplotchannel,FMT="(A)") 'unset border'
write(gnuplotchannel,FMT="(A)") 'unset xtics'
write(gnuplotchannel,FMT="(A)") 'unset ytics'
write(gnuplotchannel,FMT="(A)") 'set angles radian'
write(gnuplotchannel,FMT=*) 'r = ', max_TranslationalEnergy
write(gnuplotchannel,FMT="(A)") 'rmax = r*1.1'
write(gnuplotchannel,FMT="(A)") 'pi = 3.14159'
write(variable_length_text,FMT="(I5)") occurence_max
write(gnuplotchannel,FMT="(A)") 'set cbrange [0:'//trim(adjustl(variable_length_text))//']'
write(gnuplotchannel,FMT="(A)") 'set cblabel "Frequency"'
write(gnuplotchannel,FMT="(A)") 'set palette defined ( 0 0 1 0, 0.3333 0 0 1, 0.6667 1 0 0,\'
write(gnuplotchannel,FMT="(A)") '     1 1 0.6471 0 )'
write(gnuplotchannel,FMT="(A)") 'scaling_factor = 2.4'
write(gnuplotchannel,FMT="(A)") 'set size ratio -1'
write(gnuplotchannel,FMT="(A)") 'set xrange[0:scaling_factor*rmax]'
write(gnuplotchannel,FMT="(A)") 'set yrange[0:scaling_factor*rmax]'
write(gnuplotchannel,FMT="(A)") 'set colorbox user origin 0.9,0.1 size 0.03,0.5'
write(gnuplotchannel,FMT="(A)") 'set origin 0.0, 0.0'
write(gnuplotchannel,FMT="(A)") 'splot "'//gridpath5//temporaryfile1//'" u (scaling_factor*$2*cos($1)):'//&
                        '(scaling_factor*$2*sin($1)):3'
write(gnuplotchannel,FMT="(A)") 'set style line 11 lc rgb "black" lw 2'
write(gnuplotchannel,FMT="(A)") 'set samples 1000'
write(gnuplotchannel,FMT=*) 'phi_max = pi/', angle_slice
write(gnuplotchannel,FMT="(A)") 'dphi = phi_max/3'
write(gnuplotchannel,FMT="(A)") 'dR = rmax/4'
write(gnuplotchannel,FMT="(A)") 'set xr [0:rmax]'
write(gnuplotchannel,FMT="(A)") 'set yr [0:rmax]'
write(gnuplotchannel,FMT="(A)") 'set xtics out nomirror'
write(gnuplotchannel,FMT="(A)") 'unset ytics'
write(gnuplotchannel,FMT="(A)") 'set style line 42 lc rgb "black" dt 3'
write(gnuplotchannel,FMT="(A)") 'unset key'
write(gnuplotchannel,FMT="(A)") 'set origin 0.08, 0.09'
write(gnuplotchannel,FMT="(A)") 'set lmargin at screen 0.050'
write(gnuplotchannel,FMT="(A)") 'set rmargin at screen 0.77'
write(gnuplotchannel,FMT="(A)") 'set bmargin at screen 0.13'
write(gnuplotchannel,FMT="(A)") 'unset title'
write(gnuplotchannel,FMT="(A)") "set format x '%.3f'"
write(gnuplotchannel,FMT="(A)") 'set xtics 0, rmax/4, rmax'
!write(gnuplotchannel,FMT="(A)") 'set xlabel "Translational Energy (eV)"'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Translational Energy (kcal/mol)"'
write(gnuplotchannel,FMT="(A)") 'set for [i=0:phi_max/dphi+1] label at first (rmax*1.05)*cos(i*dphi),'//&
                        ' first (rmax*1.05)*sin(i*dphi) center sprintf(''%d^o'',i*dphi*1.02*180/pi)'
write(gnuplotchannel,FMT="(A)") 'plot for [i=1:ceil(rmax/dR)] "+" u (i*dR*cos($1*phi_max/rmax)):(i*dR*sin($1*phi_max/rmax))'//&
                        ' w l ls 42'! scale 0.8, 0.8'
write(gnuplotchannel,FMT="(A)") 'set origin 0.08, 0.09'
write(gnuplotchannel,FMT="(A)") 'plot for [i=0:phi_max/dphi+1] "+" u ($1*cos(i*dphi)):($1*sin(i*dphi))'//&
                        ' w l ls 42'! scale 0.8, 0.8'
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)


end subroutine getScatteringAngles1







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getConvergenceImage
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine takes some binning of the library's trajectories and finds
!               how it converges; it produces error bars along the way which are stored
!               in a separate file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               lowerlimit                      REAL(DP)                        The lower bound of the distribution
!               upperlimit                      REAL(DP)                        The upper bound of the distribution
!               SATRVcolumn                     INTEGER                         The column in the SATRV file this distribution
!                                                                               is recorded in
!               SATRVname                       CHARACTER(*)                    The name of the distribution (needs to be exact)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntesttraj                       INTEGER                         The number of trajectories per sample
!                                                                               set of trajectories, or the "sampling"
!               scatteringangleBins             INTEGER                         The number of bins used for distributions
!                                                                               (all of them, in this case)
!               bin_width                       REAL(DP)                        The size of the binning
!               Nsamples_max                    INTEGER                         The number of samples the library has
!                                                                               in total
!
!               binTotal                        INTEGER,DIM(                    The size of each bin of the distribution
!                                               scatteringangleBins,            for every sample
!                                               Nsamples_max)
!               binCumulative                   INTEGER,DIM(                    The size of each bin of the cumulative
!                                               scatteringangleBins,            distribution for every sample; this is
!                                               Nsamples_max)                   used in the Kolmogorov-Smirnov difference
!               sampleSize                      INTEGER,DIM(                    The number of trajectories per sample;
!                                               Nsamples_max)                   (should be constant)
!
!               binAverage                      REAL,DIM(                       The average per bin over all samples;
!                                               scatteringangleBins)            this is used in the final distribution
!               binSD                           REAL,DIM(                       The standard deviation per bin over all samples;
!                                               scatteringangleBins)            this is used in the final distribution
!
!               sampleKS                        REAL,DIM(                       The Kolmogorov-Smirnov difference per
!                                               Nsamples_max)                   sample
!               sampleRMSD                      REAL,DIM(                       The root mean square difference per
!                                               Nsamples_max)                   sample
!               sampleKRP                       REAL,DIM(                       The KRP per sample
!                                               Nsamples_max)
!
!               lambda_penalty                  REAL                            The value of lambda used to
!                                                                               calculate the KRP
!               minsd_penalty                   REAL                            The minimum standard deviation used to
!                                                                               calculate the KRP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath1//#traj//              DAT                             Stores the scattering angle and energy change
!                 binnedSATRVfile                                               decomposition for all trajectories in the library;
!                                                                               its is already binned
!               gridpath0//RMSD//SATRVname      DAT                             Stores the averages, devations, and other data
!                 cumulativefile                                                associated with a particular SATRV for all
!                                                                               trajectories in the library for some sampling
!               gridpath0//Adjusted//           DAT                             Similar to cumulative file but normalized and
!                 SATRVname                                                     dependent only on the number of sets and number
!                                                                               of trajectories per set
!               gridpath0//Convergence//        PNG                             The reference distribution of the library for
!                 SATRVname                                                     some sampling and some SATRV; also listed are
!                                                                               the RMSD, KRP, and running average
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getConvergenceImage(lowerlimit,upperlimit,SATRVcolumn,SATRVname)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

real,intent(in) :: lowerlimit, upperlimit
character(*),intent(in) :: SATRVname
integer,intent(in) :: SATRVcolumn
integer,dimension(5) :: SATRVdata

!SCATTERING ANGLE BINNING
integer :: Nsamples, Nsamples_max
real :: bin_width, binMean, binRMSD
real,allocatable :: binAverage(:), binSD(:),sampleKS(:), sampleRMSD(:), sampleKRP(:)
integer,allocatable :: binCumulative(:,:), binTotal(:,:),sampleSize(:)
integer :: binTally
real :: binThreshold = 1.0

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text

!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5

!INCREMENTAL INTEGER
integer :: i

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

!Some initialization
print *, "convergence start"
print *, "upperlimit:", upperlimit
print *, "lowerlimit:", lowerlimit
print *, "SAbins:", scatteringangleBins
bin_width = (upperlimit-lowerlimit) / scatteringangleBins
Nsamples_max = (Ngrid_max * Ntraj_max) / Ntesttraj
allocate(binCumulative(scatteringangleBins,Nsamples_max),&
         binTotal(scatteringangleBins,Nsamples_max),&
         binAverage(scatteringangleBins),binSD(scatteringangleBins),&
         sampleSize(Nsamples_max),sampleKS(Nsamples_max),&
         sampleKRP(Nsamples_max),sampleRMSD(Nsamples_max))

!The minimum standard deviation (used to calculate KRP) of a bin is
!set to be the theoretical standard deviation if a value for the
!observable occurs only once over the whole set. Thus:
!
! minBin = 1 / (Nsamples_max * Ntesttraj)
!
!  minSD = sqrt( SUM( (minBin - valueBin)^2 ) /
!                (Nsamples_max - 1)                   )
!
!        = sqrt( ((Nsamples_max - 1) * (minBin - 0)^2
!                        + (1) * (minBin - 1)^2 )) /
!                (Nsamples_max - 1)                   )
!        = sqrt( ((Nsamples_max - 1) * (1 / (Nsamples_max * Ntesttraj)^2)
!                        + ( 1 - 1/ (Nsamples_max * Ntesttraj))^2 ) /
!                (Nsamples_max - 1)                   )
!        = sqrt( ((Nsamples_max - 1) / (Nsamples_max * Ntesttraj)^2
!                        + ((Nsamples_max * Ntesttraj) - 1)^2 / (Nsamples_max * Ntesttraj)^2) /
!                (Nsamples_max - 1)                   )
!        = sqrt( ((Nsamples_max * Ntesttraj - 1)^2 + Nsamples_max - 1) /
!                (Nsamples_max - 1)                   ) / (Nsamples_max * Ntesttraj)
!
minsd_penalty = sqrt( ((Nsamples_max * Ntesttraj - 1)**2 + Nsamples_max - 1) * 1.0 / &
                      (Nsamples_max - 1) ) * 1.0 / (Nsamples_max * Ntesttraj)

!And the lambda used in the KRP is just set to be 2
!to mimic the RMSD
lambda_penalty = 2.0

!Now we want a reference distribution
!This will be the distribution of ALL trajectories
binTotal = 0
binCumulative = 0
sampleSize = 0
write(variable_length_text,FMT=FMT5_variable)&
        trajectory_text_length
write(Ntraj_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        Ngrid_max * Ntraj_max
write(variable_length_text,FMT=FMT5_variable)&
        Ngrid_text_length
write(Ngrid_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        Ngrid_max
open(filechannel1,file=gridpath5//&
        Ntraj_text//binnedSATRVfile,action="read")
do Nsamples = 1, Nsamples_max
    do i = 1, Ntesttraj
        read(filechannel1,FMT=*) SATRVdata
        binTotal(SATRVdata(SATRVcolumn),Nsamples) = &
                 binTotal(SATRVdata(SATRVcolumn),Nsamples) + 1
    end do
    sampleSize(Nsamples) = Nsamples
end do
close(filechannel1)

do i = 1,scatteringangleBins
    binAverage(i) = sum(binTotal(i,:)) * 1.0 / Nsamples_max
    binSD(i) = sqrt(sum((binTotal(i,:)*1.0 - binAverage(i))**2)/(Nsamples_max - 1))
end do

!Now we want error bars
!For this, we will get a running average distribution by adding Ntesttraj trajectories at a time
!When the running average converges then we stop and see the variance
binTally = 0
sampleKS = 0
sampleRMSD = 0
sampleKRP = 0
open(filechannel1,file=gridpath5//"RMSD"//SATRVname//cumulativefile//".dat",action="write")
do Nsamples = 1, Nsamples_max
    binRMSD = 0.0
    do i = 1, scatteringangleBins
        binRMSD = binRMSD + (sum(binTotal(i,1:Nsamples)) * 1.0 / Nsamples - binAverage(i))**2
        sampleRMSD(Nsamples) = sampleRMSD(Nsamples) + (binTotal(i,Nsamples) - binAverage(i))**2
        binCumulative(i,Nsamples) = sum(binTotal(1:i,Nsamples))
        sampleKS(Nsamples) = max(sampleKS(Nsamples),abs(binCumulative(i,Nsamples)- &
                                 sum(binAverage(1:i))))
        sampleKRP(Nsamples) = sampleKRP(Nsamples) + (abs(binTotal(i,Nsamples) - binAverage(i)) /&
                                                    max(minsd_penalty,binSD(i)))**lambda_penalty
    end do
    sampleRMSD(Nsamples) = sqrt(sampleRMSD(Nsamples) * 1.0 / (scatteringangleBins - 1.0))
    sampleKRP(Nsamples) = (sampleKRP(Nsamples) * 1.0 / (scatteringangleBins - 1.0))**(1.0/lambda_penalty)

    binRMSD = sqrt(binRMSD/scatteringangleBins)
    write(filechannel1,FMT=*) Nsamples*Ntesttraj, binRMSD, binThreshold,&
                              sampleKS(Nsamples), sampleRMSD(Nsamples), sampleKRP(Nsamples)

    if (binRMSD < binThreshold) then
        binTally = binTally + 1
    else
        binTally = 0
    end if

!   if (binTally == 5) exit
!   if (Nsamples == Nsamples_max) then
!       print *, "    No convergence for true scattering angle"
!       print *, ""
!       exit
!   end if
end do
close(filechannel1)

!Now we must do the laborious job of binning these
!Each bin has its own average and standard deviation based on how many samples we took
open(filechannel1,file=gridpath5//"Adjusted"//SATRVname//cumulativefile//".dat")
do i = 1, scatteringangleBins
    write(filechannel1,*) (i-0.5)*bin_width, binAverage(i), binSD(i)
end do
close(filechannel1)

!We have everything we need to draw the distribution with error bars
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") "set terminal pngcairo size 1200,1200"
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//'Convergence'//&
                        trim(adjustl(variable_length_text))//SATRVname//'.png"'
write(gnuplotchannel,FMT="(A)") 'set multiplot layout 2,2 columnsfirst'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,FMT="(A)") 'set title "Convergence of the '//trim(adjustl(variable_length_text))//&
                        ' '//SATRVname//' Distribution with '//trim(adjustl(Ngrid_text))//' Grids"'
write(gnuplotchannel,FMT="(A)") 'unset key'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Number of Trajectories"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "RMSD of the Distribution Over a Running Average"'
write(gnuplotchannel,FMT="(A)") 'set autoscale x'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,FMT="(A)") 'set xtics '//trim(adjustl(variable_length_text))
write(gnuplotchannel,FMT="(A)") 'set autoscale y'
write(variable_length_text1,FMT="(I5)") Ntesttraj
write(variable_length_text2,FMT="(F5.3)") binThreshold
write(gnuplotchannel,FMT="(A)") 'set label 1 "Convergence Threshold" at first '//&
                        variable_length_text1//','//variable_length_text2
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat" u 1:2 w lines, \'
write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat" u 1:3 w lines lc -1'
write(gnuplotchannel,FMT="(A)") 'unset label 1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,FMT="(A)") 'set title "Final '//SATRVname//' Distribution for '//&
                        trim(adjustl(variable_length_text))//' Trajectories"'
if (SATRVname == "ScatteringAngle") then
write(gnuplotchannel,FMT="(A)") 'scaling = 1'
write(gnuplotchannel,FMT="(A)") 'pi = 3.14159265'
write(gnuplotchannel,*) 'lowerlimit = ', lowerlimit
write(gnuplotchannel,*) 'upperlimit = ', upperlimit
write(gnuplotchannel,FMT="(A,I8,A)") 'set xtics lowerlimit, 10*(upperlimit-lowerlimit)/',&
                        scatteringangleBins,', upperlimit'
write(gnuplotchannel,FMT="(A)") "set format x '%.3P π'"
write(gnuplotchannel,FMT="(A)") 'set xrange [lowerlimit:upperlimit]'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Scattering Angle (rad)"'
else
write(gnuplotchannel,FMT="(A)") 'scaling = 1000'
write(gnuplotchannel,*) 'E_min = scaling * ', lowerlimit
write(gnuplotchannel,*) 'E_max = scaling * ', upperlimit
write(gnuplotchannel,FMT="(A)") 'set xrange [E_min:E_max]'
write(gnuplotchannel,FMT="(A,I8,A)") 'set xtics E_min, 10*(E_max-E_min)/',scatteringangleBins,', E_max'
write(gnuplotchannel,FMT="(A)") "set format x '%.3f'"
write(gnuplotchannel,FMT="(A)") 'set xlabel "Energy (meV)"'
end if
write(gnuplotchannel,FMT="(A)") 'set autoscale y'
write(gnuplotchannel,FMT="(A)") 'set ytics autofreq'
write(gnuplotchannel,FMT="(A)") 'set boxwidth'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Frequency"'
write(gnuplotchannel,FMT="(A)") 'set style fill solid 1.0 noborder'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//'.dat" u (scaling*($1)):2 w boxes, \'
write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
write(gnuplotchannel,FMT="(A)") 'set title "Distribution of the KRP Among\n'//&
                        trim(adjustl(variable_length_text))//' '//SATRVname//' Samplings Across '//&
                        trim(adjustl(Ngrid_text))//' Grids"'
write(gnuplotchannel,FMT="(A)") 'set xlabel "KRP"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'minKRP = ', max(0.0,minval(sampleKRP)-2*(maxval(sampleKRP)-minval(sampleKRP))/Nsamples_max)
write(gnuplotchannel,*) 'maxKRP = ', maxval(sampleKRP) + 2*(maxval(sampleKRP)-minval(sampleKRP))/Nsamples_max
write(gnuplotchannel,FMT="(A)") 'set xrange [minKRP:maxKRP]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,*) 'box_width = (maxKRP-minKRP) / ',(2*Nsamples_max)
write(gnuplotchannel,FMT="(A)") 'set xtics minKRP, (maxKRP-minKRP)/4, maxKRP'
write(gnuplotchannel,FMT="(A)") "set format x '%.2f'"
write(gnuplotchannel,FMT="(A)") 'set ytics 1'
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,FMT="(A)") 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat"'//&
                        'u (rounded($6)):(1.0) smooth frequency with boxes'
write(gnuplotchannel,FMT="(A)") 'set title "Distribution of the Root Mean Square Difference Among\n'//&
                        trim(adjustl(variable_length_text))//' '//SATRVname//' Samplings Across '//&
                        trim(adjustl(Ngrid_text))//' Grids"'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Root Mean Square Difference"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'minRMSD = ', max(0.0,minval(sampleRMSD)-2*(maxval(sampleRMSD)-minval(sampleRMSD))/Nsamples_max)
write(gnuplotchannel,*) 'maxRMSD = ', maxval(sampleRMSD) + 2*(maxval(sampleRMSD)-minval(sampleRMSD))/Nsamples_max
write(gnuplotchannel,FMT="(A)") 'set xrange [minRMSD:maxRMSD]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,*) 'box_width = (maxRMSD-minRMSD) / ',(2*Nsamples_max)
write(gnuplotchannel,FMT="(A)") 'set xtics minRMSD, (maxRMSD-minRMSD)/4, maxRMSD'
write(gnuplotchannel,FMT="(A)") "set format x '%.2f'"
write(gnuplotchannel,FMT="(A)") 'set ytics 1'
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,FMT="(A)") 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat"'//&
                        'u (rounded($5)):(1.0) smooth frequency with boxes'
close(gnuplotchannel)

deallocate(binAverage,binSD)
deallocate(binCumulative)
deallocate(binTotal)
deallocate(sampleKS,sampleRMSD,sampleSize,sampleKRP)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getConvergenceImage




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getScatteringAngles2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine reads in frames from the file corresponding to the prefix provided,
!               and then proceeds to visualize them in a distribution with a reference distribution
!               if provided with one
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
!               PNGfilename                     CHARACTER(*)                    The suffix we use for the output PNG
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntraj                           INTEGER                         The number of trajectories in this
!                                                                               set of trajectories, or the "sampling"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//Adjusted//           DAT                             Stores the averages and deviations for the 
!                 SATRVname                                                     distribution of some SATRV for all the
!                                                                               trajectories in the library for some sampling
!               gridpath0//PNGfilename          PNG                             Distributions of scattering angle and energy
!                                                                               change decomposition for some set of
!                                                                               trajectories
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getScatteringAngles2(prefix_filename,PNGfilename)
use PARAMETERS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

integer :: iostate
real :: speed_out, ScatteringAngle

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

inquire(file=gridpath5//'AdjustedScatteringAngle'//cumulativefile//'.dat',exist=grid_is_done)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo enhanced size 1200,1200'
write(gnuplotchannel,FMT="(A)") 'set encoding utf8'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,FMT="(A)") 'unset key'
write(gnuplotchannel,FMT="(A)") 'pi = 3.14159265'
write(gnuplotchannel,*) 'Nbins = ', energychangeBins
write(gnuplotchannel,FMT="(A)") 'scaling = 1'
write(gnuplotchannel,FMT="(A)") 'box_width = pi / Nbins'
write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
write(gnuplotchannel,FMT="(A)") 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,FMT="(A)") 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//variable_length_text//")") Ntraj
write(gnuplotchannel,FMT="(A)") 'set multiplot layout 4,1 margins 0.10,0.95,.1,.95 spacing 0,0.1'//&
                        'title "Scattering Angle and '//&
                        'Energy Change Distributions of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",18" offset 0,3'
write(gnuplotchannel,FMT="(A)") 'set style histogram clustered gap 1'
write(gnuplotchannel,FMT="(A)") 'set style fill solid 1.0 noborder'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Occurence"'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'

write(gnuplotchannel,FMT="(A)") 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,FMT="(A)") 'set xrange [0:pi]'
write(gnuplotchannel,FMT="(A)") 'set xtics pi/2'
write(gnuplotchannel,FMT="(A)") "set format x '%.1P π'"

if (grid_is_done) then
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$1>=pi?(pi-0.5*box_width):'//&
                                '(rounded($1))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedScatteringAngle'//&
                                cumulativefile//'.dat" u 1:2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedScatteringAngle'//&
                                cumulativefile//'.dat" u 1:2:3 w yerrorbars'
else
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$1>=pi?(pi-0.5*box_width):'//&
                                '(rounded($1))):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,FMT="(A)") 'scaling = 1000'

        write(gnuplotchannel,FMT="(A)") 'set xlabel "Absolute Translational Energy Change (meV)"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_absenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_absenergychange
        write(gnuplotchannel,*) 'Nbins = ', energychangeBins
        write(gnuplotchannel,FMT="(A)") 'box_width = (max_E-min_E) / Nbins'
        write(gnuplotchannel,FMT="(A)") 'set xrange [min_E:max_E]'
        write(gnuplotchannel,FMT="(A)") 'set xtics min_E, box_width*Nbins/5, max_E'
        write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
        write(gnuplotchannel,FMT="(A)") 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,FMT="(A)") 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,FMT="(A)") "set format x '%.3f'"
if (grid_is_done) then
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$3>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($3))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedAbsoluteEnergyChange'//&
                                cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedAbsoluteEnergyChange'//&
                                cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$3>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($3))):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,FMT="(A)") 'set xlabel "Relative Translational Energy Change (meV)"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_relenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_relenergychange
        write(gnuplotchannel,*) 'Nbins = ', energychangeBins
        write(gnuplotchannel,FMT="(A)") 'box_width = (max_E-min_E) / Nbins'
        write(gnuplotchannel,FMT="(A)") 'set xrange [min_E:max_E]'
        write(gnuplotchannel,FMT="(A)") 'set xtics min_E, box_width*Nbins/5, max_E'
        write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
        write(gnuplotchannel,FMT="(A)") 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,FMT="(A)") 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
if (grid_is_done) then
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$4>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($4))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedRelativeEnergyChange'//&
                                cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedRelativeEnergyChange'//&
                                cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$4>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($4))):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,FMT="(A)") 'set xlabel "Rotational Energy Change (meV)"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_rotenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_rotenergychange
        write(gnuplotchannel,*) 'Nbins = ', energychangeBins
        write(gnuplotchannel,FMT="(A)") 'box_width = (max_E-min_E) / Nbins'
        write(gnuplotchannel,FMT="(A)") 'set xrange [min_E:max_E]'
        write(gnuplotchannel,FMT="(A)") 'set xtics min_E, box_width*Nbins/5, max_E'
        write(gnuplotchannel,FMT="(A)") "set format x '%.3f'"
        write(gnuplotchannel,FMT="(A)") 'set boxwidth box_width'
        write(gnuplotchannel,FMT="(A)") 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,FMT="(A)") 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
if (grid_is_done) then
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$5>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($5))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedRotationalEnergyChange'//&
                                cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,FMT="(A)") '     "'//gridpath5//'AdjustedRotationalEnergyChange'//&
                                cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$5>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($5))):(1.0) smooth frequency w boxes'
end if
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getScatteringAngles2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getComparedScatteringAngles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine reads all selected set of trajectories and compares
!               them with respect to some number of trajectories and some SATRV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               lowerlimit                      REAL(DP)                        The lower bound of the distribution
!               upperlimit                      REAL(DP)                        The upper bound of the distribution
!               imagename                       CHARACTER(*)                    The name of the output image
!               SATRVcolumn                     INTEGER                         The column in the SATRV file this distribution
!                                                                               is recorded in
!               SATRVname                       CHARACTER(*)                    The name of the distribution (needs to be exact)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntesttraj                       INTEGER                         The number of trajectories per sample
!                                                                               set of trajectories, or the "sampling"
!               Nbins                           INTEGER                         The number of bins used for distributions
!               comparison_number               INTEGER                         The number of sets of trajectories being
!                                                                               compared
!
!               binTotal                        INTEGER,DIM(                    The size of each bin of the distribution
!                                               Nbins,comparison_number         for each set of trajectories
!
!               referenceBins                   REAL,DIM(Nbins)                 The x-value corresponding to a bin
!               referenceMeans                  REAL,DIM(Nbins)                 The y-value corresponding to a bin
!               referenceSDs                    REAL,DIM(Nbins)                 The y-value uncertainty corresponding to a bin
!
!               comparisonCDF                   REAL                            The cumulative distribution function of some bin
!               comparisonKS                    REAL                            The Kolmogorov-Smirnov difference of some bin
!               comparisonRMSD                  REAL                            The RMSD of some distribution
!               comparisonKRP                   REAL                            The KRP of some distribution
!
!               lambda_penalty                  REAL                            The value of lambda used to
!                                                                               calculate the KRP
!               minsd_penalty                   REAL                            The minimum standard deviation used to
!                                                                               calculate the KRP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//Adjusted//           DAT                             Stores the averages and deviations for the 
!                 SATRVname                                                     distribution of some SATRV for all the
!                                                                               trajectories in the library for some sampling
!               gridpath0//allprefixes(...)//   DAT                             Same as the SATRV file but the values are
!                 binnedSATRVfile                                               stored as integers corresponding to some
!                                                                               bin in a specific binning
!               gridpath0//imagename            PNG                             Compares distributions as selected in
!                                                                               allprefixes by some sampling and some SATRV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getComparedScatteringAngles(lowerlimit,upperlimit,imagename,SATRVcolumn,SATRVname)
use PARAMETERS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
integer,intent(in) :: SATRVcolumn
character(*), intent(in) :: SATRVname
integer,dimension(5) :: SATRVdata

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: imagename

!Upper and Lower Limits for the plot
real,intent(in) :: upperlimit, lowerlimit

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(8) :: difference_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

integer :: iostate,i,j
real :: speed_out, ScatteringAngle

real,allocatable :: referenceBins(:)
real,allocatable :: referenceMeans(:)
real,allocatable :: referenceSDs(:)
integer,allocatable :: binTotal(:,:)
real :: comparisonRMSD, comparisonKS, comparisonCDF, comparisonKRP
integer :: Nbins

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

!Initialization
if (SATRVname == "ScatteringAngle") then
    Nbins = scatteringangleBins
else
    Nbins = energychangeBins
end if

allocate(referenceBins(Nbins),&
         referenceMeans(Nbins),&
         referenceSDs(Nbins))
allocate(binTotal(Nbins,comparison_number))

open(filechannel1,file=gridpath5//"Adjusted"//&
        SATRVname//cumulativefile//".dat")
do j = 1, Nbins
    read(filechannel1,FMT=*) referenceBins(j),&
            referenceMeans(j), referenceSDs(j)
end do
close(filechannel1)

binTotal = 0
open(filechannel1,file=gridpath5//&
        allprefixes(1:alllengths(1)-1)//binnedSATRVfile)
do j = 1, Ntesttraj
    read(filechannel1,FMT=*) SATRVdata
    binTotal(SATRVdata(SATRVcolumn),1) = &
             binTotal(SATRVdata(SATRVcolumn),1) + 1
end do
close(filechannel1)

do i = 1, comparison_number-1
    open(filechannel1,file=gridpath5//&
            allprefixes(sum(alllengths(1:i))+1:&
            sum(alllengths(1:i+1))-1)//binnedSATRVfile)
    do j = 1, Ntesttraj
        read(filechannel1,FMT=*) SATRVdata
        binTotal(SATRVdata(SATRVcolumn),i+1) = &
                 binTotal(SATRVdata(SATRVcolumn),i+1) + 1
    end do
    close(filechannel1)
end do

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 1200,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath4//imagename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset xtics'
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I"//variable_length_text//")") Ntraj
write(gnuplotchannel,*) 'set multiplot layout ',comparison_number,&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"'//SATRVname//' Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",32" offset 0,-3'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set ytics font ",16"'

if (SATRVname == "ScatteringAngle") then
    write(gnuplotchannel,*) 'scaling = 1'
    write(gnuplotchannel,*) 'pi = 3.14159265'
else
    write(gnuplotchannel,*) 'scaling = 1000'
end if

write(gnuplotchannel,*) 'lowerlimit = scaling * ', lowerlimit
write(gnuplotchannel,*) 'upperlimit = scaling * ', upperlimit
write(gnuplotchannel,*) 'set xrange [lowerlimit:upperlimit]'
write(gnuplotchannel,*) 'box_width = (upperlimit - lowerlimit) / ', Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(scaling * x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (x - 0.5)'

comparisonRMSD = 0
comparisonKS = 0
comparisonCDF = 0
comparisonKRP = 0
do j = 1, Nbins
    comparisonCDF = comparisonCDF + 1.0*binTotal(j,1) - referenceMeans(j)
    comparisonKS = max(comparisonKS, abs(comparisonCDF))
    comparisonRMSD = comparisonRMSD + (1.0*bintotal(j,1) - referenceMeans(j))**2
    comparisonKRP = comparisonKRP + (abs(1.0*bintotal(j,1) - referenceMeans(j)) / &
                                    max(minsd_penalty,referenceSDs(j)))**lambda_penalty
end do

write(gnuplotchannel,*) 'set label 1 "'//allprefixes(1:alllengths(1)-1)//'" at graph 0.825, 0.9'
write(difference_text,FMT="(F8.4)") sqrt(comparisonRMSD * 1.0 / (Nbins - 1.0))
write(gnuplotchannel,*) 'set label 2" RMSD: '//difference_text//'" at graph 0.85,0.825'
write(difference_text,FMT="(F8.4)") (comparisonKRP * 1.0 / (Nbins - 1.0))**(1.0/lambda_penalty)
write(gnuplotchannel,*) 'set label 3" KRP: '//difference_text//'" at graph 0.85,0.750'
write(variable_length_text,FMT=FMT5_variable) SATRVcolumn
write(gnuplotchannel,*) 'plot "'//gridpath5//allprefixes(1:alllengths(1)-1)//&
                        binnedSATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
                        ')):(1.0) smooth frequency w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                        '.dat" u (scaling*($1)):2 w boxes'//&
                        ' fs transparent solid 0.5 noborder, \'
write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                        '.dat" u (scaling*($1)):2:3 w yerrorbars'
do i = 1, comparison_number-1
    comparisonRMSD = 0
    comparisonKS = 0
    comparisonCDF = 0
    comparisonKRP = 0
    do j = 1, Nbins
        comparisonCDF = comparisonCDF + 1.0*binTotal(j,i+1) - referenceMeans(j)
        comparisonKS = max(comparisonKS, abs(comparisonCDF))
        comparisonRMSD = comparisonRMSD + (1.0*binTotal(j,i+1) - referenceMeans(j))**2
        comparisonKRP = comparisonKRP + (abs(1.0*bintotal(j,i+1) - referenceMeans(j)) / &
                                        max(minsd_penalty,referenceSDs(j)))**lambda_penalty
    end do

    write(gnuplotchannel,*) 'set label ', 3*i+1, '"'//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1))-1)//&
                            '" at graph 0.825, 0.9'
    write(difference_text,FMT="(F8.4)") sqrt(comparisonRMSD * 1.0 / (Nbins - 1.0))
    write(gnuplotchannel,*) 'set label ',3*i+2,'" RMSD: '//difference_text//'" at graph 0.85,0.825'
    write(difference_text,FMT="(F8.2)") (comparisonKRP * 1.0 / (Nbins - 1.0))**(1.0/lambda_penalty)
    write(gnuplotchannel,*) 'set label ',3*i+3,'" KRP: '//difference_text//'" at graph 0.85,0.750'
    write(gnuplotchannel,*) 'unset label ', 3*i-2
    write(gnuplotchannel,*) 'unset label ', 3*i-1
    write(gnuplotchannel,*) 'unset label ', 3*i
    if (i == comparison_number-1) then
        if (SATRVname == "ScatteringAngle") then
            write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)" font ",24"'
            write(gnuplotchannel,*) 'set xtics lowerlimit, 10*(upperlimit-lowerlimit)/',&
                                    Nbins,', upperlimit font ",16"'
            write(gnuplotchannel,*) "set format x '%.3P π'"
        else
            write(gnuplotchannel,*) 'set xlabel "Energy Change (meV)" font ",24"'
            write(gnuplotchannel,*) 'set xtics lowerlimit, 10*box_width, upperlimit font ",16"'
            write(gnuplotchannel,*) "set format x '%.3f'"
        end if
    end if
    write(variable_length_text,FMT=FMT5_variable) SATRVcolumn
    write(gnuplotchannel,*) 'plot "'//gridpath5//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1))-1)//&
                            binnedSATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
                            ')):(1.0) smooth frequency w boxes, \'
    write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                            '.dat" u (scaling*($1)):2 w boxes'//&
                            ' fs transparent solid 0.5 noborder, \'
    write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                            '.dat" u (scaling*($1)):2:3 w yerrorbars'
end do
close(gnuplotchannel)

deallocate(referenceBins,referenceMeans,referenceSDs)
deallocate(binTotal)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getComparedScatteringAngles


subroutine postProcess(prefix_filename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

!FILENAMING
character(*),intent(in) :: prefix_filename

!Translational, Rotational, Vibrational Energies
real(dp) :: scatteringAngle, speed_in, speed_out
real(dp) :: bond_distance, rot_speed
real(dp) :: relTranslationalEnergy1, relTranslationalEnergy2
real(dp) :: absTranslationalEnergy1, absTranslationalEnergy2
real(dp) :: RotationalEnergy1, RotationalEnergy2
real(dp) :: TotalEnergy1, TotalEnergy2

!TRAJECTORY DATA
!real(dp),dimension(Nbonds,6) :: initial_bonding_data
real(dp),dimension(3,Natoms) :: coords_initial,velocities_initial
real(dp),dimension(3,Natoms) :: coords_final,velocities_final
real(dp),dimension(3) :: velocity_in, velocity_out
real(dp),dimension(3) :: velocity1, velocity2
real(dp),dimension(3) :: coords1, coords2
real(dp),dimension(3) :: velocity1_CM, velocity2_CM

!BONDING DATA INDEXES
integer :: atom1, atom2

!INCREMENTAL INTEGERS
integer :: i, j, k

    open(filechannel1,file=gridpath0//prefix_filename//SATRVfile)
    do k = 1, Ntraj
        read(filechannel1,FMTdata) scatteringAngle, speed_out, &
                                  abs_energychange, &
                                  rel_energychange, &
                                  rot_energychange

        max_TranslationalEnergy = max(speed_out,max_TranslationalEnergy)

        max_absenergychange = max(abs_energychange,max_absenergychange)
        min_absenergychange = min(abs_energychange,min_absenergychange)
        max_relenergychange = max(rel_energychange,max_relenergychange)
        min_relenergychange = min(rel_energychange,min_relenergychange)
        max_rotenergychange = max(rot_energychange,max_rotenergychange)
        min_rotenergychange = min(rot_energychange,min_rotenergychange)
    end do
    close(filechannel1)

end subroutine postProcess


end module analyzeScatteringAngleswithMultipleGrids
