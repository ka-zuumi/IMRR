!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               analyzeRMSDThresholdwithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This module plots the percentage of frames of a trajectory that are
!               below a certain threshold; these DAT files are prepared beforehand
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!		FILECHANNEL1			OPEN, WRITE, CLOSE
!		FILEHCANNEL2			OPEN, READ, CLOSE
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!               gridpath1//prefix_filename//    DAT                     Stores the RMSD of the frame accepted at
!                       _#traj                                          each timestep of some trajectory coresponding
!                                                                       to some prefix
!		gridpath0//percent_rmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!               gridpath0//JPGfilename          JPG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module analyzeRMSDThresholdwithMultipleGrids
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getRMSDThresholds1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks through the RMSD of the accepted frames over all trajectories
!               corresponding to prefix_filename, compares them to some threshold, and plots the
!               percentage of frames below the threshold as a distribution per number of grids used
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
!               JPGfilename                     CHARACTER(*)                    The suffix we use for the output JPG
!               selection                       INTEGER                         An integer describing which grids to plot
!                                                                               by representing "plot" as a 1 and "no
!                                                                               plot" as a 0 in binary
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntesttraj                       INTEGER                         The number of trajectories in this
!                                                                               set of trajectories
!               Ngrid_total                     INTEGER                         The number of grids used in this
!                                                                               set of trajectories
!               Ngrid_plotting                  INTEGER                         The number of grids to be plotted
!               grid_selection                  INTEGER,DIM(Ngrid_max)          An array with indexes of grids to plot
!               RMSD_Nbins                      INTEGER                         The number of bins for this distribution
!
!               threshold_rmsd                  REAL                            The RMSD threshold
!               frames                          INTEGER                         The number of frames read so far
!               total_threshold_rmsd            INTEGER                         The number of frames whose accepted frame
!                                                                               has an RMSD under the threshold so far
!               percent_threshold_rmsd          REAL                            The percentage of frames in a set of
!                                                                               trajectories whose accepted frame has an
!                                                                               RMSD under the threshold
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                     Stores the RMSD of the frame accepted at
!                       Ngrid_#traj                                     each timestep of some trajectory corresponding
!                                                                       to some prefix
!		gridpath5//percent_rmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!               gridpath4//JPGfilename          PNG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getRMSDThresholds1(prefix_filename,JPGfilename,selection)
use PARAMETERS
use ANALYSIS
use SIMILARITY
implicit none

!NUMBER OF TRAJECTORIES CHECKED
integer :: n_testtraj

!NUMBER OF GRIDS TO BE PLOTTED
integer :: Ngrid_plotting

!INTERNAL TALLY FOR RMSD BELOW A THRESHOLD
real :: percent_threshold_rmsd
integer :: total_threshold_rmsd, frames

!ARRAY HOLDING DATA FROM FILE
real(dp) :: min_rmsd, min_rmsd_prime

!ARRAY HOLDING SELECTION OF GRIDS TO PLOT
integer,optional,intent(in) :: selection
integer,dimension(Ngrid_max) :: grid_selection
integer :: selection_current

!FORMAT OF JPG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: JPGfilename
character(*), intent(in) :: prefix_filename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text

integer :: i1,i2,i3,i4
real(dp) :: r1,r2,r3,r4
integer :: count1, count2, total1, total2

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

!If we only want to plot a select number of grids,
!then that can be passed in by the user as the
!variable "selection"
if (present(selection)) then
    Ngrid_plotting = 0
    selection_current = selection
    do Ngrid = 1, Ngrid_max
        selection_current = selection_current&
                - 2 ** (Ngrid - 1)

        if (selection_current < 0) exit

        !If the binary representation of the
        !selection has a 1 in that digits place
        !then that is a signal by the user that
        !the grid should be plotted
        if (modulo(selection_current,2**Ngrid) == 0) then
            Ngrid_plotting = Ngrid_plotting + 1
            grid_selection(Ngrid_plotting) = Ngrid
        else
            selection_current = selection_current&
                    + 2 ** (Ngrid - 1)
        end if
    end do

!If not supplied, then we will just plot all
!grids up until Ngrid_total
else
    Ngrid_plotting = Ngrid_total
    do Ngrid = 1, Ngrid_plotting
        grid_selection(Ngrid) = Ngrid
    end do
end if


!Now we need to open up each grid's set of trajectories
do Ngrid = 1, Ngrid_plotting
write(variable_length_text,"(I5)")&
        Ngrid_text_length
write(Ngrid_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        grid_selection(Ngrid)

!We will bin data by GRID, not by trajectory
!So we uniquely name each output .dat and graph by the
!grid number
open(filechannel1,file=gridpath5//&
        "percent_rmsd"//Ngrid_text//".dat")
do n_testtraj = 1, Ntesttraj
    write(variable_length_text,"(I5)")&
            trajectory_text_length
    write(Ntraj_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            n_testtraj

    !Read the trajectory (which has the rmsd) across
    !all grids line-by-line
    frames = 0
    total_threshold_rmsd = 0
    open(filechannel2,file=gridpath0//&
            prefix_filename//Ngrid_text//&
            "_"//Ntraj_text//".dat")

    do
        !Tally how many frames we have in the trajectory
        read(filechannel2,FMT=*,iostat=iostate) min_rmsd
        if (iostate /= 0) exit
        frames = frames + 1

        !And if the RMSD is below the threshhold we
        !tally that separately
        if (.not.(accept_worst) .and. &
                 (min_rmsd < outer_threshold_SI))&
                 total_threshold_rmsd = &
                 total_threshold_rmsd + 1
        if ((accept_worst) .and. &
                (min_rmsd > 0.0d0))&
                 total_threshold_rmsd = &
                 total_threshold_rmsd + 1
    end do
    close(filechannel2)

    !We want the percentage of frames that has an
    !RMSD below the threshhold so we keep track of
    !the number of frames and divide by that
    percent_threshold_rmsd = &
            total_threshold_rmsd * 100.0 / frames
    write(filechannel1,FMT="(I6,1x,F7.3,1x,I8)")&
            n_testtraj, percent_threshold_rmsd, frames
end do
close(filechannel1)

end do

!Finally, plot the data
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term jpeg size 1200,1200'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//JPGfilename//'.jpg"'
write(gnuplotchannel,FMT="(A)") 'set tmargin 0'
write(gnuplotchannel,FMT="(A)") 'set bmargin 0'
write(gnuplotchannel,FMT="(A)") 'set lmargin 1'
write(gnuplotchannel,FMT="(A)") 'set rmargin 1'
write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I"//trim(adjustl(variable_length_text))//")") Ngrid_plotting
write(gnuplotchannel,FMT="(A)") 'set multiplot layout '//trim(adjustl(Ngrid_text))//&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0'&
                        //' title "Trajectory RMSD Distribution with '//expfolder//' method"'
write(gnuplotchannel,FMT="(A)") 'set title "Percentages of Trajectories with RMSD Below Threshold"'
write(gnuplotchannel,FMT="(A)") 'set style fill solid 1.0 noborder'
write(gnuplotchannel,FMT="(A)") 'scaling = 1'
write(gnuplotchannel,FMT="(A)") 'unset key'
write(gnuplotchannel,FMT="(A)") 'unset xtics'
write(gnuplotchannel,FMT="(A)") 'unset xlabel'
write(gnuplotchannel,FMT="(A)") 'ymax = 100.0'
write(variable_length_text,"(I5)") RMSD_Nbins
write(gnuplotchannel,FMT="(A)") 'Nbins = '//trim(adjustl(variable_length_text))
write(gnuplotchannel,FMT="(A)") 'bin_width = ymax / Nbins'
write(gnuplotchannel,FMT="(A)") 'set boxwidth bin_width'
write(gnuplotchannel,FMT="(A)") 'min(x,y) = (x < y) ? x : y'
write(gnuplotchannel,FMT="(A)") 'bin_number(x) = min(floor(x/bin_width),Nbins-1)'
write(gnuplotchannel,FMT="(A)") 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,FMT="(A)") 'set xrange [0:ymax]'
write(gnuplotchannel,FMT="(A)") 'set yrange [0:]'
write(gnuplotchannel,FMT="(A)") 'unset ylabel'

do Ngrid = 1, Ngrid_plotting
if (Ngrid == Ngrid_plotting) then
    write(gnuplotchannel,FMT="(A)") 'set label 1 "Occurence" at screen 0.01,0.45 rotate by 90'
    write(gnuplotchannel,FMT="(A)") 'set xtics'
    write(gnuplotchannel,FMT="(A)") 'set xlabel "Percentage of Frames with RMSD Below Threshold"'
end if

write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u (rounded($2)):(1.0/scaling) smooth frequency with boxes'
write(gnuplotchannel,FMT="(A)") 'unset title'
end do

write(gnuplotchannel,FMT="(A)") 'set title "Distribution of Trajectories, Percent RMSD vs. Trajectory Length"'
write(gnuplotchannel,FMT="(A)") 'scaling = 1000'
write(gnuplotchannel,FMT="(A)") 'set autoscale y'
write(gnuplotchannel,FMT="(A)") 'unset xtics'
write(gnuplotchannel,FMT="(A)") 'unset xlabel'
write(gnuplotchannel,FMT="(A)") 'set xrange [0:ymax]'
write(gnuplotchannel,FMT="(A)") 'unset ylabel'

write(variable_length_text,"(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_plotting
if (Ngrid == Ngrid_plotting) then
    write(gnuplotchannel,FMT="(A)") 'set label 2 "Length of Trajectory (Thousands of Frames)" at screen 0.51,0.40 rotate by 90'
    write(gnuplotchannel,FMT="(A)") 'set xtics'
    write(gnuplotchannel,FMT="(A)") 'set xlabel "Percentage of Frames with RMSD Below Threshold"'
end if

write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u 2:(($3)/scaling) with points'
write(gnuplotchannel,FMT="(A)") 'unset title'

end do

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

write(variable_length_text,"(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_plotting
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
call system("rm "//gridpath5//"percent_rmsd"//Ngrid_text//".dat")
end do


!frames = 0
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 1
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,I8)") i4, max(min_rmsd,.00001), i1
!frames = frames + 1
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 2
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,I8)") i4, max(min_rmsd_prime,.00001), i1
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!total1 = 0
!total2 = 0
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 3
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,F9.6,1x,I8)") i4, max(abs(min_rmsd - min_rmsd_prime),0.00001), &
!                                             min_rmsd - min_rmsd_prime, i1
!if ((min_rmsd-min_rmsd_prime) == 0.0d0) then
!else if ((min_rmsd-min_rmsd_prime) < 0.0d0) then
!        total1 = total1 + 1
!        total2 = total2 + 1
!else
!        total1 = total1 + 1
!end if
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!print *, "total1 (number of nonzero RMSD differences): ", total1
!print *, "total2 (number of negative RMSD differences): ", total2
!
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 4
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,I8)") i4, max(min_rmsd_prime/min(min_rmsd,1.0),.000001), i1
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!
!
!
!Ngrid_plotting = 4
!
!!Finally, plot the data
!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
!write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'.jpg"'
!write(gnuplotchannel,*) 'set tmargin 0'
!write(gnuplotchannel,*) 'set bmargin 0'
!write(gnuplotchannel,*) 'set lmargin 1'
!write(gnuplotchannel,*) 'set rmargin 1'
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I"//trim(adjustl(variable_length_text))//")") Ngrid_plotting
!write(gnuplotchannel,*) 'set multiplot layout '//trim(adjustl(Ngrid_text))//&
!                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0'&
!                        //' title "Trajectory RMSD Distribution with '//prefix_filename//' method"'
!write(gnuplotchannel,*) 'set title "RMSD Histogram Over the Trajectory"'
!write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
!write(gnuplotchannel,*) 'scaling = ', frames
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'unset xtics'
!write(gnuplotchannel,*) 'unset xlabel'
!!write(gnuplotchannel,*) 'ymax = 100.0'
!write(gnuplotchannel,*) 'ymax =', 0.00001
!write(variable_length_text,"(I5)") RMSD_Nbins
!write(gnuplotchannel,*) 'Nbins = '//trim(adjustl(variable_length_text))
!!write(gnuplotchannel,*) 'bin_width = ymax / Nbins'
!write(gnuplotchannel,*) 'bin_width = (2.0 - log10(ymax)) / Nbins'
!write(gnuplotchannel,*) 'set boxwidth bin_width'
!write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
!!write(gnuplotchannel,*) 'bin_number(x) = min(floor(x/bin_width),Nbins-1)'
!write(gnuplotchannel,*) 'bin_number(x) = min(floor(log10(x)/bin_width),Nbins-1)'
!write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
!!write(gnuplotchannel,*) 'set xrange [0:ymax]'
!write(gnuplotchannel,*) 'set xrange [log10(ymax):2.0]'
!write(gnuplotchannel,*) 'set yrange [0:0.2]'
!write(gnuplotchannel,*) 'unset ylabel'
!
!do Ngrid = 1, Ngrid_plotting
!if (Ngrid == Ngrid_plotting) then
!        write(gnuplotchannel,*) 'set label 1 "Occurence" at screen 0.01,0.45 rotate by 90'
!	write(gnuplotchannel,*) 'set xtics'
!	write(gnuplotchannel,*) 'set xlabel "log(RMSD)"'
!end if
!
!if (Ngrid == 3) then
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (rounded($2)):(1.0/scaling) lc rgb "blue" '//&
!                                'smooth frequency with boxes,\'
!        write(gnuplotchannel,*) '     "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u ($3>0?(rounded($3)):1/0):(1.0/scaling) lc rgb "red" '//&
!                                'smooth frequency with boxes'
!        write(gnuplotchannel,*) 'unset title'
!else
!        write(variable_length_text,"(I5)") Ngrid_text_length
!        !write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (rounded($2)):(1.0/scaling) smooth frequency with boxes'
!        write(gnuplotchannel,*) 'unset title'
!end if
!end do
!
!write(gnuplotchannel,*) 'set title "RMSD Distribution vs. The Number of Frames Checked"'
!write(gnuplotchannel,*) 'scaling = 100'
!write(gnuplotchannel,*) 'set autoscale y'
!write(gnuplotchannel,*) 'unset xtics'
!write(gnuplotchannel,*) 'unset xlabel'
!!write(gnuplotchannel,*) 'set xrange [0:ymax]'
!write(gnuplotchannel,*) 'set xrange [log10(ymax):2.0]'
!write(gnuplotchannel,*) 'unset ylabel'
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!do Ngrid = 1, Ngrid_plotting
!if (Ngrid == Ngrid_plotting) then
!        write(gnuplotchannel,*) 'set label 2 "Number of Frames Checked (Hundreds of Frames)" at screen 0.51,0.40 rotate by 90'
!	write(gnuplotchannel,*) 'set xtics'
!	write(gnuplotchannel,*) 'set xlabel "log(RMSD)"'
!end if
!!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
!!write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!!                        '" u 2:(($3)/scaling) with points'
!
!if (Ngrid == 3) then
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (log10($2)):(($4)/scaling) lc rgb "blue",\'
!        write(gnuplotchannel,*) '     "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u ($3>0?(log10($3)):1/0):(($4)/scaling) lc rgb "red"'
!        write(gnuplotchannel,*) 'unset title'
!else
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (log10($2)):(($3)/scaling) with points'
!        write(gnuplotchannel,*) 'unset title'
!end if
!end do
!
!close(gnuplotchannel)
!
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!do Ngrid = 1, Ngrid_plotting
!!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!!call system("rm "//gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!end do


end subroutine getRMSDThresholds1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getRMSDDifferences1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at the RMSD acceptance rates between two trials, particularly
!               for different subcell neighbor checks (subcellsearch_max)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
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
!               Nbins                           INTEGER                         The number of bins for this distribution
!               threshold_rmsd                  REAL                            The RMSD threshold
!               frames                          INTEGER                         The number of frames read so far
!               number_of_frames_accepted       INTEGER                         The number of frames whose accepted frame has
!                                                                               an RMSD under the threshold so far
!               min_rmsd                        REAL                            The RMSD of the accepted frame (trial 1)
!               min_rmsd_prime                  REAL                            The RMSD of the accepted frame (trial 2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath5//checkstatefile       DAT                     Stores the RMSD of the accepted frames of
!                                                                       each trial (and other data)
!		gridpath5//percent_rmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!               gridpath4//PNGfilename          PNG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRMSDDifferences1(PNGfilename)
use PARAMETERS
use ANALYSIS
use SIMILARITY
implicit none

!NUMBER OF FRAMES IN THE TRAJECTORY
integer :: frames

!RMSD
real(dp) :: min_rmsd, min_rmsd_prime
real(dp) :: min_min_rmsd
integer :: number_of_frames_accepted

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

!OTHER VARIABLES IN THE CHECKSTATE FILE
integer :: i1,i2,i3,i4
real(dp) :: r1,r2,r3,r4

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
integer :: Nbins
real(dp) :: bin_width

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

frames = 0
number_of_frames_accepted = 0
min_min_rmsd = default_SIs(1)

open(filechannel2,file=gridpath5//checkstatefile)
!read(filechannel1) number_of_frames,order,i&
!                   neighbor_check,steps,&
!                   min_rmsd,min_rmsd_prime,&
!                   vals(1),vals(2),U,KE
do 
    read(filechannel2,FMT=*,iostat=iostate)&
            i1, i2, i3, i4,&
            min_rmsd,min_rmsd_prime,&
            r1,r2,r3,r4
    if (iostate /= 0) exit

    min_min_rmsd = min(min_min_rmsd,&
            min_rmsd,min_rmsd_prime)
    if (min_rmsd_prime < outer_threshold_SI) &
        number_of_frames_accepted =&
        number_of_frames_accepted + 1
    frames = frames + 1
end do
close(filechannel2)


Nbins = 50
bin_width = 2 * log10(default_SIs(1)/&
        min_min_rmsd) / Nbins

write(variable_length_text,"(I5)")&
        Ngrid_text_length
write(Ngrid_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        4
open(filechannel1,file=gridpath5//&
        "percent_rmsd"//Ngrid_text//".dat")

open(filechannel2,file=gridpath5//checkstatefile)
!read(filechannel1) number_of_frames,order,&
!                   neighbor_check,steps,&
!                   min_rmsd,min_rmsd_prime,&
!                   vals(1),vals(2),U,KE
do 
    read(filechannel2,FMT=*,iostat=iostate)&
            i1, i2, i3, i4,&
            min_rmsd,min_rmsd_prime,&
            r1,r2,r3,r4
    if (iostate /= 0) exit

    write(filechannel1,FMT="(I6,1x,F9.6,1x,"//&
            "F9.6,1x,I8,1x,I8,1x,I8)") &
            i4, min_rmsd, min_rmsd_prime, i2, &
            nint(log10((min_rmsd_prime)/&
            (min_min_rmsd/default_SIs(1))) / (0.5*bin_width)),&
            nint(log10((min_rmsd/min_rmsd_prime)/&
            (min_min_rmsd/default_SIs(1))) / bin_width)
end do
close(filechannel2)
close(filechannel1)



open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "RMSD Encountered in Binary Check"'
write(gnuplotchannel,*) 'set key left top'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'min_min_rmsd = ', min_min_rmsd
write(gnuplotchannel,*) 'set xrange [min_min_rmsd:',default_SIs(1),']'
write(gnuplotchannel,*) 'set yrange [min_min_rmsd:',default_SIs(1),']'
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Encountered for the First Check"'
write(gnuplotchannel,*) 'set ylabel "RMSD (A) Encountered for the Second Check"'
write(gnuplotchannel,*) 'plot "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '($4==0?$2:1/0):3 w p lc rgb "red" title "Both Order 0",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '($4==1?$2:1/0):3 w p lc rgb "blue" title "One Order 1",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '($4==2?$2:1/0):3 w p lc rgb "green" title "Both Order 1"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'_1.png"'
write(gnuplotchannel,*) 'set title "Ratio of RMSD Encountered in Binary Check"'
write(gnuplotchannel,*) 'set key left top'
write(gnuplotchannel,*) 'scaling = ', frames
write(gnuplotchannel,*) 'xmin = ', min_min_rmsd / default_SIs(1)
write(gnuplotchannel,*) 'xmax = ', default_SIs(1) / min_min_rmsd
write(gnuplotchannel,*) 'set xlabel "RMSD (first check) / RMSD (second check)"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'

write(gnuplotchannel,*) 'Nbins = ', Nbins
write(gnuplotchannel,*) 'bin_width = ', bin_width
write(gnuplotchannel,*) 'set boxwidth 1'
write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
write(gnuplotchannel,*) 'set xrange [0:Nbins]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style fill transparent solid 0.5'

write(gnuplotchannel,*) 'set xtics ('//&
                                           '"10^{-5}" (-5-log10(xmin))/bin_width, '//&
                                           '"10^{-4}" (-4-log10(xmin))/bin_width, '//&
                                           '"10^{-3}" (-3-log10(xmin))/bin_width, '//&
                                           '"10^{-2}" (-2-log10(xmin))/bin_width, '//&
                                           '"10^{-1}" (-1-log10(xmin))/bin_width, '//&
                                               '"10^0" (0-log10(xmin))/bin_width, '//&
                                               '"10^1" (1-log10(xmin))/bin_width, '//&
                                               '"10^2" (2-log10(xmin))/bin_width, '//&
                                               '"10^3" (3-log10(xmin))/bin_width, '//&
                                               '"10^4" (4-log10(xmin))/bin_width, '//&
                                               '"10^5" (5-log10(xmin))/bin_width)'

write(gnuplotchannel,*) 'plot "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '6:($4==2?(1.0/scaling):0.0) '//&
                               'smooth frequency w boxes lc rgb "green" title "Both Order 1",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '6:($4==1?(1.0/scaling):0.0) '//&
                               'smooth frequency w boxes lc rgb "blue" title "One Order 1",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '6:($4==0?(1.0/scaling):0.0) '//&
                               'smooth frequency w boxes lc rgb "red" title "Both Order 0"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'_2.png"'
write(gnuplotchannel,*) 'set title "Ratio of RMSD Encountered for Binary Check"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'scaling = ', frames
write(gnuplotchannel,*) 'xmin = ', min_min_rmsd
write(gnuplotchannel,*) 'xmax = ', default_SIs(1)
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Encountered After the Second Check"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'

write(variable_length_text,"(F5.2)") number_of_frames_accepted * 100.0 / frames
write(gnuplotchannel,*) 'set label 1 "Acceptance Rate: '//variable_length_text//'%" at graph 0.1,0.8'

write(gnuplotchannel,*) 'Nbins = ', Nbins
write(gnuplotchannel,*) 'bin_width = ', 0.5*bin_width
write(gnuplotchannel,*) 'set boxwidth 1'
write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
write(gnuplotchannel,*) 'set xrange [-0.5:Nbins+.05]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style fill transparent solid 0.5'

write(gnuplotchannel,*) 'set xtics ('//&
                                           '"10^{-5}" (-5-log10(xmin))/bin_width, '//&
                                           '"10^{-4}" (-4-log10(xmin))/bin_width, '//&
                                           '"10^{-3}" (-3-log10(xmin))/bin_width, '//&
                                           '"10^{-2}" (-2-log10(xmin))/bin_width, '//&
                                           '"10^{-1}" (-1-log10(xmin))/bin_width, '//&
                                               '"10^0" (0-log10(xmin))/bin_width, '//&
                                               '"10^1" (1-log10(xmin))/bin_width, '//&
                                               '"10^2" (2-log10(xmin))/bin_width, '//&
                                               '"10^3" (3-log10(xmin))/bin_width, '//&
                                               '"10^4" (4-log10(xmin))/bin_width, '//&
                                               '"10^5" (5-log10(xmin))/bin_width)'

write(gnuplotchannel,*) 'plot "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '5:(1.0/scaling) '//&
                               'smooth frequency w boxes lc rgb "green"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)



end subroutine getRMSDDifferences1


subroutine getAlphaErrorDistribution(PNGfilename)
use PARAMETERS
use ANALYSIS
implicit none

integer :: frames
integer,dimension(Nalpha) :: alpha_flagging
integer,dimension(Nalpha) :: alpha_binning
real(dp),dimension(Nalpha) :: alpha_array

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

alpha_binning = 0
frames = 0

open(6666,file=gridpath5//alphafile)
do
    read(6666,FMT=*,iostat=iostate)&
            (alpha_flagging(n),n=1,Nalpha)

    if (iostate /= 0) exit

    frames = frames + 1
    alpha_binning = alpha_binning +&
            alpha_flagging

end do
close(6666)

do n = 1, Nalpha
    alpha_array(n) = alpha_start + (n-1-0.5) * &
        (alpha_end - alpha_start) / (Nalpha-1)
end do

if (logarithmic_alpha_flag) then
    do n = 1, Nalpha
        alpha_array(n) = &
             10.0d0 ** alpha_array(n)
    end do
end if

open(6666,file=gridpath5//"tmpA.dat")
do n = 1, Nalpha
    write(6666,FMT=*) alpha_array(n),&
            alpha_binning(n) * 1.0d2 / frames
end do
close(6666)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "Blue Point Occurence by Alpha Ratio"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set yrange [0:', &
        min(100.0d0,maxval(alpha_binning)*110.0d0/frames),']'

if (logarithmic_alpha_flag) then
    write(gnuplotchannel,*) 'set logscale x'
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set xrange [',&
            10.0d0 ** (alpha_start - (alpha_end-alpha_start)/Nalpha),&
            10.0d0 ** (alpha_end + (alpha_end-alpha_start)/Nalpha),']'
else
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set xrange [',&
            (alpha_start - (alpha_end-alpha_start)/Nalpha),&
            (alpha_end + (alpha_end-alpha_start)/Nalpha),']'
end if


write(gnuplotchannel,*) 'set xlabel "Alpha Ratio"'
write(gnuplotchannel,*) 'set ylabel "Percentage of Frames with IE > AE"'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'plot "'//gridpath5//'tmpA.dat" u 1:2 w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getAlphaErrorDistribution



subroutine getRMSDinterpolation(vals,delta_vals,PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real(dp),dimension(Nvar),intent(in) :: vals, delta_vals
character(5) :: vals_interpolation_text

!NUMBER OF FRAMES IN THE DATA
real(dp),dimension(7) :: min_rmsd_vals, max_rmsd_vals
!real(dp),dimension(7) :: min_rsv, max_rsv, rsv

!RMSD
real(dp) :: rmsd_x, rmsd_y, rmsd_z, rmsd_fx
real(dp) :: min_rmsd_x, max_rmsd_x
real(dp) :: min_rmsd_y, max_rmsd_y
real(dp) :: min_rmsd_z, max_rmsd_z
real(dp) :: min_rmsd_fx, max_rmsd_fx

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

!VARIABLES IN THE INTERPOLATION FILE
integer :: Ninterpolation
integer :: min_Ninterpolation, max_Ninterpolation
real :: vals1, vals2
!real(dp),dimension(2) :: temp_vals
real(dp),dimension(3) :: RMSDheatmap_coeff

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
integer :: Nbins

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

!frames = 0
!tally1 = 0
!tally2 = 0
!
!min_rmsd_x = default_rmsd
!min_rmsd_y = default_rmsd
!min_rmsd_z = default_rmsd
!min_rmsd_fx = 1.0d9
!min_Ninterpolation = 1000
!max_rmsd_x = 0.0d0
!max_rmsd_y = 0.0d0
!max_rmsd_z = 0.0d0
!max_rmsd_fx = 0.0d0
!max_Ninterpolation = 0
!
!write(vals_interpolation_text,FMT="(F7.3,'_',F7.3)") vals(1),vals(2)
!open(filechannel1,file=gridpath0//vals_interpolation_text//interpolationfile)
!open(filechannel2,file=gridpath0//interpolationfile)
!do 
!        read(filechannel2,FMT=*,iostat=iostate) vals1, vals2, Ninterpolation, &
!                                       rmsd_y, rmsd_z, &
!                                       rmsd_x_prime, rmsd_fx_prime, &
!                                       rmsd_x, rmsd_fx
!        if (iostate /= 0) exit
!
!        if ((abs(vals1-vals(1)) > delta_vals(1)).or.&
!            (abs(vals2-vals(2)) > delta_vals(2))) cycle
!
!        min_rmsd_x = min(min_rmsd_x, rmsd_x)
!        max_rmsd_x = max(max_rmsd_x, rmsd_x)
!        min_rmsd_y = min(min_rmsd_y, rmsd_y)
!        max_rmsd_y = max(max_rmsd_y, rmsd_y)
!        min_rmsd_z = min(min_rmsd_z, rmsd_z)
!        max_rmsd_z = max(max_rmsd_z, rmsd_z)
!        min_rmsd_fx = min(min_rmsd_fx, rmsd_fx)
!        max_rmsd_fx = max(max_rmsd_fx, rmsd_fx)
!        min_Ninterpolation = min(min_Ninterpolation, Ninterpolation)
!        max_Ninterpolation = max(max_Ninterpolation, Ninterpolation)
!
!        frames = frames + 1
!!       if (Ninterpolation == 1) tally1 = tally1 + 1
!!       if (Ninterpolation == 2) tally2 = tally2 + 1
!
!        write(filechannel1,FMT=*) vals1, vals2, Ninterpolation,&
!               rmsd_y, rmsd_z, rmsd_x, rmsd_fx
!
!end do
!close(filechannel1)
!close(filechannel2)
!
!if (frames == 0) return
!
!Nbins = 100
!min_Nbin = 10
!max_Nbin = 0
!bin_width = log10((max_rmsd_fx/min_rmsd_x) / (min_rmsd_fx/max_rmsd_x)) / Nbins
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 4
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//interpolationfile)
!do 
!        read(filechannel2,FMT=*,iostat=iostate) vals1, vals2, Ninterpolation, &
!                                       interpolation_alpha, threshold_rmsd_1, &
!                                       rmsd_x_prime, rmsd_fx_prime, &
!                                       rmsd_x, rmsd_fx
!        if (iostate /= 0) exit
!
!        if ((abs(vals1-vals(1)) > delta_vals(1)).or.&
!            (abs(vals2-vals(2)) > delta_vals(2))) cycle
!
!        Nbin = nint(log10((rmsd_fx/rmsd_x) / (min_rmsd_fx/max_rmsd_x)) / bin_width)
!
!        min_Nbin = min(min_Nbin,Nbin)
!        max_Nbin = max(max_Nbin,Nbin)
!
!        write(filechannel1,FMT=*) vals1, vals2, Nbin
!end do
!close(filechannel2)
!close(filechannel1)
!
!bin_width1 = log10(max_rmsd_z / min_rmsd_z) / Nbins
!bin_width2 = log10(max_rmsd_x / min_rmsd_x) / Nbins
!
!allocate(RMSDheatmap(Nbins,Nbins))
!RMSDheatmap = 1.0d-7
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 4
!open(filechannel2,file=gridpath0//interpolationfile)
!do 
!        read(filechannel2,FMT=*,iostat=iostate) vals1, vals2, Ninterpolation, &
!                                       rmsd_y, rmsd_z, &
!                                       rmsd_x_prime, rmsd_fx_prime, &
!                                       rmsd_x, rmsd_fx
!        if (iostate /= 0) exit
!
!        if ((abs(vals1-vals(1)) > delta_vals(1)).or.&
!            (abs(vals2-vals(2)) > delta_vals(2))) cycle
!
!        Nrmsd1 = nint(log10(rmsd_z/min_rmsd_z)/bin_width1)
!        if (Nrmsd1 < 1) Nrmsd1 = 1
!        Nrmsd2 = nint(log10(rmsd_x/min_rmsd_x)/bin_width2)
!        if (Nrmsd2 < 1) Nrmsd2 = 1
!
!        RMSDheatmap(Nrmsd1,Nrmsd2) = max(&
!                RMSDheatmap(Nrmsd1,Nrmsd2),rmsd_fx)
!end do
!close(filechannel2)
!
!Nheatmap = 0
!do Nrmsd1 = 1, Nbins
!        do Nrmsd2 = 1, Nbins
!                if (RMSDheatmap(Nrmsd1,Nrmsd2) > 1.0d-7) then
!                         Nheatmap = Nheatmap + 1
!                end if
!        end do
!end do
!
!allocate(A(Nheatmap,3),b(Nheatmap))
!Nheatmap = 0
!do Nrmsd1 = 1, Nbins
!        do Nrmsd2 = 1, Nbins
!                if (RMSDheatmap(Nrmsd1,Nrmsd2) > 1.0d-7) then
!                         Nheatmap = Nheatmap + 1
!                         A(Nheatmap,:) = (/ (Nrmsd1-0.5)*bin_width1,&
!                                 (Nrmsd2-0.5)*bin_width2,1.0d0/)
!                         b(Nheatmap) = RMSDheatmap(Nrmsd1,Nrmsd2)
!                end if
!        end do
!end do
!
!call LS(A,Nheatmap,3,b,RMSDheatmap_coeff)
!
!open(filechannel2,file=gridpath0//"heatmap_rmsd"//Ngrid_text//".dat")
!do Nrmsd1 = 1, Nbins
!        do Nrmsd2 = 1, Nbins
!                write(filechannel2,FMT=*) (Nrmsd1-0.5)*bin_width1,&
!                        (Nrmsd2-0.5)*bin_width2,RMSDheatmap(Nrmsd1,Nrmsd2),&
!                        max((Nrmsd1-0.5)*bin_width1*RMSDheatmap_coeff(1)+&
!                            (Nrmsd2-0.5)*bin_width2*RMSDheatmap_coeff(2)+&
!                                                    RMSDheatmap_coeff(3),1.0d-7)
!        end do
!        write(filechannel2,*) ""
!end do
!close(filechannel2)
!

Nbins = 100

call processInterpolationFile(vals,delta_vals,&
        RMSDheatmap_coeff,&
        min_rmsd_vals,max_rmsd_vals,&
        min_Ninterpolation,max_Ninterpolation)

write(vals_interpolation_text,FMT="(F7.3,'_',F7.3)") vals(1),vals(2)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,2400'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "RMSD Comparison of a Frame and Gradient with Interpolation"'
write(gnuplotchannel,*) 'set multiplot layout 3,1'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.1,0.900'
write(gnuplotchannel,FMT="(A,I6,A)") &
        'set label 2 "Ntraj = ', Ntraj, '" at screen 0.1,0.875'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_x
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_x
!write(gnuplotchannel,*) 'min_y = ', min_rmsd_fx
!write(gnuplotchannel,*) 'max_y = ', max_rmsd_fx
write(gnuplotchannel,*) 'min_x = ', min_rmsd_vals(2)
write(gnuplotchannel,*) 'max_x = ', max_rmsd_vals(2)
write(gnuplotchannel,*) 'min_y = ', min_rmsd_vals(6)
write(gnuplotchannel,*) 'max_y = ', max_rmsd_vals(6)
write(gnuplotchannel,*) 'min_cx = ', min_Ninterpolation
write(gnuplotchannel,*) 'max_cx = ', max_Ninterpolation
write(gnuplotchannel,*) 'set xrange [min_x:max_x]'
write(gnuplotchannel,*) 'set yrange [min_y:max_y]'
write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx+1]'
write(gnuplotchannel,*) 'set palette defined (min_cx "blue", max_cx "red")'
write(gnuplotchannel,*) 'set xlabel "Largest RMSD Between the Input Frame and the Interpolation Points"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between the Output Gradient and the Interpolation"'
write(gnuplotchannel,*) 'set cblabel "Number of Frames Used to Interpolate"'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'plot "'//gridpath5//vals_interpolation_text//interpolationfile//'" u '//&
                               '(($5)):7:3 w p lw 6 palette'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "SQRT of Weighted Largest RMSD Squared Between the Input Frame and the Interpolation Points"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between the Output Gradient and the Interpolation"'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_y
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_y
write(gnuplotchannel,*) 'min_x = ', min_rmsd_vals(1)
write(gnuplotchannel,*) 'max_x = ', max_rmsd_vals(1)
write(gnuplotchannel,*) 'set xrange [0.5*sqrt(min_x):2*sqrt(max_x)]'
write(gnuplotchannel,*) 'plot "'//gridpath5//vals_interpolation_text//interpolationfile//'" u '//&
                               '(sqrt($4)):7:3 w p lw 6 palette'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "Weighted RMSD Between the Input Frame and the Interpolation Points"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between the Output Gradient and the Interpolation"'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_z
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_z
write(gnuplotchannel,*) 'min_x = ', min_rmsd_vals(5)
write(gnuplotchannel,*) 'max_x = ', max_rmsd_vals(5)
write(gnuplotchannel,*) 'set xrange [0.5*min_x:2*max_x]'
write(gnuplotchannel,*) 'plot "'//gridpath5//vals_interpolation_text//interpolationfile//'" u '//&
                               '(($6)):7:3 w p lw 6 palette'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//"heatmap"//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "RMSD Comparison of a Frame and Gradient with Interpolation"'
write(gnuplotchannel,*) 'set multiplot layout 1,2'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.1,0.800 front'
write(gnuplotchannel,FMT="(A,I6,A)") &
        'set label 2 "Ntraj = ', Ntraj, '" at screen 0.1,0.775 front'
write(gnuplotchannel,*) 'xmin = ', min_rmsd_vals(1)
write(gnuplotchannel,*) 'xmax = ', max_rmsd_vals(1)
write(gnuplotchannel,*) 'ymin = ', min_rmsd_vals(5)
write(gnuplotchannel,*) 'ymax = ', max_rmsd_vals(5)
!write(gnuplotchannel,*) 'xmin = ', min_rmsd_z
!write(gnuplotchannel,*) 'xmax = ', max_rmsd_z
!write(gnuplotchannel,*) 'ymin = ', min_rmsd_x
!write(gnuplotchannel,*) 'ymax = ', max_rmsd_x
write(gnuplotchannel,*) 'min_cx = ', 1.0e-7
write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_vals(6),1.0e-7)
!write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_fx,1.0e-7)
write(gnuplotchannel,*) 'set cbrange [log10(.0000001):log10(.0001)]'
write(gnuplotchannel,*) 'set palette defined ('//&
                        'log10(.0000001) "white", '//&
                        'log10(.0000005) "yellow", '//&
                        'log10(.000001) "green", '//&
                        'log10(.000005) "cyan", '//&
                        'log10(.00001) "blue", '//&
                        'log10(.00005) "magenta", '//&
                        'log10(.0001) "red"'//&
                        ')'
write(gnuplotchannel,*) 'set xlabel "RMSD Variable 1 from Frame"'
write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2 from Frame"'
write(gnuplotchannel,*) 'set cblabel "Maximum RMSD Between the Output Gradient and the Interpolation"'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001/(xmin)), '//&
                                           '"1e-6" log10(0.000001/(xmin)), '//&
                                           '"1e-5" log10(0.00001/(xmin)), '//&
                                           '"1e-4" log10(0.0001/(xmin)), '//&
                                           '"1e-3" log10(0.001/(xmin)), '//&
                                           '"1e-2" log10(0.01/(xmin)), '//&
                                           '"1e-1" log10(0.1/(xmin)), '//&
                                               '"1e0" log10(1.0/(xmin)), '//&
                                               '"1e1" log10(10.0/(xmin)))'
write(gnuplotchannel,*) 'set ytics ('//&
                                           '"1e-7" log10(0.0000001/(ymin)), '//&
                                           '"1e-6" log10(0.000001/(ymin)), '//&
                                           '"1e-5" log10(0.00001/(ymin)), '//&
                                           '"1e-4" log10(0.0001/(ymin)), '//&
                                           '"1e-3" log10(0.001/(ymin)), '//&
                                           '"1e-2" log10(0.01/(ymin)), '//&
                                           '"1e-1" log10(0.1/(ymin)), '//&
                                               '"1e0" log10(1.0/(ymin)), '//&
                                               '"1e1" log10(10.0/(ymin)))'
write(gnuplotchannel,*) 'set cbtics ('//&
                                         '"1e-8" log10(.00000001), '//&
                                         '"5e-8" log10(.00000005), '//&
                                          '"1e-7" log10(.0000001), '//&
                                          '"5e-7" log10(.0000005), '//&
                                           '"1e-6" log10(.000001), '//&
                                           '"5e-6" log10(.000005), '//&
                                           '"1e-5"  log10(.00001), '//&
                                           '"5e-5"  log10(.00005), '//&
                                           '"1e-4"   log10(.0001), '//&
                                           '"5e-4"   log10(.0005), '//&
                                           '"1e-3"    log10(.001), '//&
                                           '"5e-3"    log10(.005), '//&
                                           '"1e-2"     log10(.01), '//&
                                           '"5e-2"     log10(.05), '//&
                                           '"1e-1"      log10(.1), '//&
                                           '"5e-1"      log10(.5), '//&
                                           ' "1.0"       log10(1), '//&
                                   ')'
write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_rmsd.dat" u '//&
                                '1:2:(log10($3)) w image palette'
write(gnuplotchannel,FMT="(A)") 'unset label 1'
write(gnuplotchannel,FMT="(A,F9.5,A,F9.5,A,F9.5,A)") &
        'set label 1 "f(x,y) = ', RMSDheatmap_coeff(1),'x + ',&
                                  RMSDheatmap_coeff(2),'y + ',&
                                  RMSDheatmap_coeff(3),'" at screen 0.7,0.800 front'
write(gnuplotchannel,FMT="(A)") 'unset label 2'
write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_rmsd.dat" u '//&
                                '1:2:(log10($4)) w image palette'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//"heatmap_freq"//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set multiplot layout 1,2'
write(gnuplotchannel,*) 'set title "RMSD Comparison of a Frame and Gradient with Interpolation"'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.1,0.800 front'
write(gnuplotchannel,FMT="(A,I6,A)") &
        'set label 2 "Ntraj = ', Ntraj, '" at screen 0.1,0.775 front'
write(gnuplotchannel,*) 'xmin = ', min(min_rmsd_vals(4),min_rmsd_vals(6))
write(gnuplotchannel,*) 'xmax = ', max(max_rmsd_vals(4),max_rmsd_vals(6))
write(gnuplotchannel,*) 'ymin = xmin'
write(gnuplotchannel,*) 'ymax = xmax'
write(gnuplotchannel,*) 'set xrange [0:log10(xmax/xmin)]'
write(gnuplotchannel,*) 'set yrange [0:log10(ymax/ymin)]'
write(gnuplotchannel,*) 'set autoscale cb'
!write(gnuplotchannel,*) 'xmin = ', min_rmsd_z
!write(gnuplotchannel,*) 'xmax = ', max_rmsd_z
!write(gnuplotchannel,*) 'ymin = ', min_rmsd_x
!write(gnuplotchannel,*) 'ymax = ', max_rmsd_x
!write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_fx,1.0e-7)
write(gnuplotchannel,*) 'set xlabel "RMSD Between Approximate and Real Gradient (N=1)"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between Approximate and Real Gradient (N>=1)"'
write(gnuplotchannel,*) 'set cblabel "Occurence"'
write(gnuplotchannel,*) 'set palette defined ('//&
                        '0 "white", '//&
                        '1 "yellow", '//&
                        '2 "green", '//&
                        '3 "cyan", '//&
                        '4 "blue", '//&
                        '5 "magenta", '//&
                        '6 "red"'//&
                        ')'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001/(xmin)), '//&
                                           '"1e-6" log10(0.000001/(xmin)), '//&
                                           '"1e-5" log10(0.00001/(xmin)), '//&
                                           '"1e-4" log10(0.0001/(xmin)), '//&
                                           '"1e-3" log10(0.001/(xmin)), '//&
                                           '"1e-2" log10(0.01/(xmin)), '//&
                                           '"1e-1" log10(0.1/(xmin)), '//&
                                               '"1e0" log10(1.0/(xmin)), '//&
                                               '"1e1" log10(10.0/(xmin)))'
write(gnuplotchannel,*) 'set ytics ('//&
                                           '"1e-7" log10(0.0000001/(ymin)), '//&
                                           '"1e-6" log10(0.000001/(ymin)), '//&
                                           '"1e-5" log10(0.00001/(ymin)), '//&
                                           '"1e-4" log10(0.0001/(ymin)), '//&
                                           '"1e-3" log10(0.001/(ymin)), '//&
                                           '"1e-2" log10(0.01/(ymin)), '//&
                                           '"1e-1" log10(0.1/(ymin)), '//&
                                               '"1e0" log10(1.0/(ymin)), '//&
                                               '"1e1" log10(10.0/(ymin)))'
write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_freq1_rmsd.dat" u '//&
                                '1:2:3 w image palette'
write(gnuplotchannel,*) 'unset pm3d'
write(gnuplotchannel,*) 'xmin = ', min_rmsd_vals(7)
write(gnuplotchannel,*) 'xmax = ', max_rmsd_vals(7)
write(gnuplotchannel,*) 'deltax = (log10(xmax/xmin)) / ', Nbins/5
write(gnuplotchannel,*) 'set boxwidth deltax'
write(gnuplotchannel,*) 'set style fill solid 1.0'
write(gnuplotchannel,*) 'set xrange [-deltax:log10(xmax/xmin)+deltax]'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set autoscale ymax'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-3" log10(0.001/(xmin)), '//&
                                           '"5e-3" log10(0.005/(xmin)), '//&
                                           '"1e-2" log10(0.01/(xmin)), '//&
                                           '"5e-2" log10(0.05/(xmin)), '//&
                                           '"1e-1" log10(0.1/(xmin)), '//&
                                           '"5e-1" log10(0.5/(xmin)), '//&
                                               '"1e0" log10(1.0/(xmin)), '//&
                                               '"5e0" log10(5.0/(xmin)), '//&
                                               '"1e1" log10(10.0/(xmin)), '//&
                                               '"5e1" log10(50.0/(xmin)), '//&
                                               '"1e2" log10(100.0/(xmin)), '//&
                                               '"5e2" log10(500.0/(xmin)))'
write(gnuplotchannel,*) 'set xlabel "Ratio of RMSD Between Approximate and Real Gradient (N=1) and (N>=1)"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'plot "'//gridpath5//'heatmap_freq2_rmsd.dat" u '//&
                               '1:2 w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)




!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
!write(gnuplotchannel,*) 'set output "'//gridpath0//PNGfilename//'_1.png"'
!write(gnuplotchannel,*) 'set title "Ratio of RMSD Between the Frame and the Gradient"'
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'scaling = ', frames
!write(gnuplotchannel,*) 'xmin = ', min_rmsd_fx / max_rmsd_x
!write(gnuplotchannel,*) 'xmax = ', max_rmsd_fx / min_rmsd_x
!write(gnuplotchannel,*) 'set xlabel "Ratio of RMSD"'
!write(gnuplotchannel,*) 'set ylabel "Occurence"'
!
!write(gnuplotchannel,*) 'Nbins = ', Nbins
!write(gnuplotchannel,*) 'bin_width = ', bin_width
!write(gnuplotchannel,*) 'set boxwidth 1'
!write(gnuplotchannel,*) 'set xrange [-0.5:Nbins+.05]'
!write(gnuplotchannel,*) 'set yrange [0:]'
!write(gnuplotchannel,*) 'set style fill transparent solid 0.5'
!
!write(gnuplotchannel,*) 'set xtics ('//&
!                                           '"1e-5" (-5-log10(xmin))/bin_width, '//&
!                                           '"1e-4" (-4-log10(xmin))/bin_width, '//&
!                                           '"1e-3" (-3-log10(xmin))/bin_width, '//&
!                                           '"1e-2" (-2-log10(xmin))/bin_width, '//&
!                                           '"1e-1" (-1-log10(xmin))/bin_width, '//&
!                                               '"1e0" (0-log10(xmin))/bin_width, '//&
!                                               '"1e1" (1-log10(xmin))/bin_width, '//&
!                                               '"1e2" (2-log10(xmin))/bin_width, '//&
!                                               '"1e3" (3-log10(xmin))/bin_width, '//&
!                                               '"1e4" (4-log10(xmin))/bin_width, '//&
!                                               '"1e5" (5-log10(xmin))/bin_width)'
!
!write(gnuplotchannel,*) 'plot "'//gridpath0//"percent_rmsd"//Ngrid_text//'.dat" u '//&
!                               '3:(1.0/scaling) '//&
!                               'smooth frequency w boxes lc rgb "green"'
!close(gnuplotchannel)
!
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)


!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
!write(gnuplotchannel,*) 'set output "'//gridpath0//PNGfilename//'_2.png"'
!write(gnuplotchannel,*) 'set title "RMSD Ratio Dependence on Variables"'
!write(gnuplotchannel,*) 'set pm3d map'
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_x
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_x
!write(gnuplotchannel,*) 'min_y = ', min_rmsd_fx
!write(gnuplotchannel,*) 'max_y = ', max_rmsd_fx
!write(gnuplotchannel,*) 'min_cx = ', min_Nbin
!write(gnuplotchannel,*) 'max_cx = ', max_Nbin
!write(gnuplotchannel,*) 'set xrange [0:',var_maxvar(1),']'
!write(gnuplotchannel,*) 'set yrange [0:',var_maxvar(2),']'
!write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx]'
!write(gnuplotchannel,*) 'set palette defined (min_cx "blue", max_cx "red")'
!write(gnuplotchannel,*) 'set xlabel "Var1 (A)"'
!write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
!write(gnuplotchannel,*) 'set cblabel "Ratio of RMSD"'
!
!write(gnuplotchannel,*) 'plot "'//gridpath0//"percent_rmsd"//Ngrid_text//'.dat" u '//&
!                               '1:2:3 w p lw 4 palette'
!close(gnuplotchannel)
!
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end subroutine getRMSDinterpolation

subroutine getRMSDErrorPlots(PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!NUMBER OF FRAMES IN THE DATA
real(dp),dimension(9) :: rsv

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename
integer, dimension(comparison_number,3) :: counters

integer,allocatable :: RMSDvsError_heatmap(:,:)
integer,allocatable :: CMDvsError_heatmap(:,:)
integer :: population_max
integer :: maxRMSDvsError_heatmap, maxCMDvsError_heatmap
integer :: RMSDbins, CMDbins, Errorbins
integer :: RMSDbin, CMDbin, Errorbin
real(dp) :: RMSDbinwidth, CMDbinwidth, Errorbinwidth
real(dp) :: RMSDlimit, CMDlimit, Errorlimit

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

call processMultipleInterpolationFiles(counters)

!    read(filechannel2,FMT=*,iostat=iostate) &
!            vals, Ninterpolation, &
!!           rsv1, rsv2, &
!!           rmsd_best, rmsd_interpolated, &
!!           min_CMdiff, interpolated_CMdiff, &
!!           error_best, error_interpolated
!            rsv(3), rsv(4), &
!            rsv(1), rsv(2), &
!            rsv(8), rsv(9), &
!            rsv(5), rsv(6)

!    tmp file: Ninterpolation, rsv

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_1.png"'
write(gnuplotchannel,*) 'set title "Error Convergence as Interpolations Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD (A) Between Target and Interpolated Frame" font ",24"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '3:7 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_2.png"'
write(gnuplotchannel,FMT="(A)") 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD (A) Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '2:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_2_linear.png"'
write(gnuplotchannel,FMT="(A)") 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD (A) Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '2:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_12.png"'
write(gnuplotchannel,FMT="(A)") 'set multiplot layout 1'//&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0,0.1 title '//&
                        '"Error Convergence as Candidates/Interpolations Get Closer" font ",36" offset 0,9'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD (A) Between Target and Candidate Frame"'//&
        ' font ",24" offset "0,3"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Candidate/Interpolated Gradient"'//&
        ' font ",24" offset -3,0'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set yrange [.001:3.0]'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'set ytics ('//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '2:6 w p'

write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                          '"1e-7" .0000001, '//&
                                           '"1e-6" .000001, '//&
                                           '"1e-5"  .00001, '//&
                                           '"1e-4"   .0001, '//&
                                           '"1e-3"    .001, '//&
                                           '"1e-2"     .01, '//&
                                           '"1e-1"      .1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD (A) Between Target and Interpolated Frame"'//&
        ' font ",24" offset "0,3"'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '3:7 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_2_linear.png"'
write(gnuplotchannel,FMT="(A)") 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD (A) Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xrange [0:]'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '2:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_CM2.png"'
write(gnuplotchannel,FMT="(A)") 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Coulomb Matrix Difference Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '9:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_CM2_linear.png"'
write(gnuplotchannel,FMT="(A)") 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A)") 'set xlabel "Coulomb Matrix Difference Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,FMT="(A)") 'set ylabel "Error (A/fs) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,FMT="(A)") 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '9:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)



!Now let's make a good old fashioned heatmap
RMSDbins = 30
CMDbins = 30
Errorbins = 30

RMSDlimit = 0.05d0
CMDlimit = 0.05d0
Errorlimit = 0.1d0

RMSDbinwidth = RMSDlimit / RMSDbins
CMDbinwidth = CMDlimit / CMDbins
Errorbinwidth = Errorlimit / Errorbins

allocate(RMSDvsError_heatmap(RMSDbins,Errorbins),&
          CMDvsError_heatmap(CMDbins,Errorbins))

RMSDvsError_heatmap = 0
CMDvsError_heatmap = 0
open(filechannel1,file=gridpath5//"tmp"//interpolationfile)
do
    read(filechannel1,iostat=iostate,FMT=*) n, rsv
    if (iostate /= 0) exit

    RMSDbin = floor(rsv(1) / RMSDbinwidth) + 1
    CMDbin = floor(rsv(8) / CMDbinwidth) + 1
    Errorbin = floor(rsv(5) / Errorbinwidth) + 1

    if (RMSDbin < 1) RMSDbin = 1
    if (CMDbin < 1) CMDbin = 1
    if (Errorbin < 1) Errorbin = 1

    if (Errorbin > Errorbins) cycle

    if (RMSDbin <= RMSDbins) &
        RMSDvsError_heatmap(RMSDbin,Errorbin) = &
        RMSDvsError_heatmap(RMSDbin,Errorbin) + 1
    if (CMDbin <= CMDbins) &
        CMDvsError_heatmap(CMDbin,Errorbin) = &
        CMDvsError_heatmap(CMDbin,Errorbin) + 1
end do
close(filechannel1)

maxRMSDvsError_heatmap = 0
open(filechannel1,file=gridpath5//temporaryfile1)
do RMSDbin = 1, RMSDbins
    do Errorbin = 1, Errorbins
        write(filechannel1,FMT=*) &
            RMSDbinwidth * (RMSDbin-0.5d0),&
            Errorbinwidth * (Errorbin-0.5d0),&
            RMSDvsError_heatmap(RMSDbin,Errorbin)
        maxRMSDvsError_heatmap = max(&
            maxRMSDvsError_heatmap,&
            RMSDvsError_heatmap(RMSDbin,Errorbin))
    end do
    write(filechannel1,FMT=*) ""
end do
close(filechannel1)

maxCMDvsError_heatmap = 0
open(filechannel1,file=gridpath5//temporaryfile2)
do CMDbin = 1, CMDbins
    do Errorbin = 1, Errorbins
        write(filechannel1,FMT=*) &
            CMDbinwidth * (CMDbin-0.5d0),&
            Errorbinwidth * (Errorbin-0.5d0),&
            CMDvsError_heatmap(CMDbin,Errorbin)
        maxCMDvsError_heatmap = max(&
            maxCMDvsError_heatmap,&
            CMDvsError_heatmap(CMDbin,Errorbin))
    end do
    write(filechannel1,FMT=*) ""
end do
close(filechannel1)

population_max = 30000
!population_max = max(maxRMSDvsError_heatmap,maxCMDvsError_heatmap)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGfilename//'_linear_heatmap.png"'
write(gnuplotchannel,FMT="(A)") 'set multiplot layout 1,2 '//&
        'margins screen 0.10,0.85,0.12,0.93 spacing 0,0.1 '//&
        'title "Error Convergence As Candidates Get Closer" '//&
        'font ",32" offset 0,-5'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set xtics font ",24"'
write(gnuplotchannel,*) 'set ytics font ",24"'
write(gnuplotchannel,*) 'set cbtics font ",24"'

write(gnuplotchannel,FMT="(A)") 'set colorbox user origin screen 0.87, screen 0.12 '//&
        'size screen 0.02, screen 0.81'
write(gnuplotchannel,*) 'population_max = ', population_max
write(gnuplotchannel,FMT="(A)") 'set palette defined ('//&
         '0 "white",'//&
         '0.5 "white",'//&
         '0.5 "yellow",'//&
         'population_max/5 "orange",'//&
         'population_max "red")'
write(gnuplotchannel,*) 'set cbrange [0:population_max]'
write(gnuplotchannel,FMT="(A)") 'set cblabel "Number of Candidate Frames"'//&
        ' font ",32" offset 6,0'

write(gnuplotchannel,FMT="(A)") 'set ylabel "Error Between Candidate '//&
        'and\nTarget Gradient (E_h/a_0)" font ",32" offset -5,0'
write(gnuplotchannel,*) 'ymax = ', Errorlimit
write(gnuplotchannel,*) 'set yrange [0:ymax]'

write(gnuplotchannel,FMT="(A)") 'set xlabel "RMSD Between Candidate '//&
        'and\nTarget Frame (A)" font ",32" offset 0,-1'
write(gnuplotchannel,*) 'xmax = ', RMSDlimit
write(gnuplotchannel,*) 'set xrange [0:xmax]'
write(gnuplotchannel,FMT="(A)") 'splot "'//gridpath5//temporaryfile1//&
                                '" u 1:2:3 w image palette'

write(gnuplotchannel,*) 'unset colorbox'
write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,FMT="(A)") 'set xlabel "CMD Between Candidate '//&
        'and\nTarget Frame (1/A)" font ",32" offset 0,-1'
write(gnuplotchannel,*) 'xmax = ', CMDlimit
write(gnuplotchannel,*) 'set xrange [0:xmax]'
write(gnuplotchannel,FMT="(A)") 'splot "'//gridpath5//temporaryfile2//&
                                '" u 1:2:3 w image palette'
close(gnuplotchannel)

deallocate(RMSDvsError_heatmap,CMDvsError_heatmap)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getRMSDErrorPlots



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               processInterpolationFile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at the processed interpolation file (truncated file)
!               and creates data files for plotting, excising anything outside of the
!               region of the phase space of interest
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               vals                            REAL,DIM(Nvar)                  The region of the phase space we are
!                                                                               interested in plotting
!               delta_vals                      REAL,DIM(Nvar)                  The size of the region of the phase
!                                                                               space of interest
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               RMSDheatmap_coeff               REAL,DIM(3)                     The coefficients of the linear fit
!                                                                               produced by this data
!               min_rsv                         REAL,DIM(7)                     Minimums of the RMSD special variables
!               max_rsv                         REAL,DIM(7)                     Maximums of the RMSD special variables
!               min_Ninterpolation              INTEGER                         Minimum of the Ninterpolation
!               max_Ninterpolation              INTEGER                         Maximum of the Ninterpolation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Nbins                           INTEGER                         The number of bins for the binned data
!               rsv_binwidth                    REAL,DIM(7)                     Width of the binning for each RSV
!               rsv_bin                         INTEGER,DIM(7)                  Temporary buffer for each RSV bin data
!               ERheatmap2D                     REAL(Nbins,Nins)                Heatmap that describes the occurence
!                                                                               of two RSV
!               ERheatmap1D                     REAL(Nbins/5)                   Heatmap that describes the occurence
!                                                                               of one RSV
!               IEheatmap2D                     REAL(Nbins,Nins)                Heatmap that describes the maximum RSV
!                                                                               (usually IE) observed given two RSV
!               Nheatmap                        INTEGER                         The number of unique nonzero entries
!                                                                               in IEheatmap2D
!               A                               REAL(Nheatmap,3)                The "independent variable matrix" for
!                                                                               IEheatmap2D used for linear fitting
!               b                               REAL(Nheatmap)                  The "dependent variable matrix" for
!                                                                               IEheatmap2D used for linear fitting
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath5//truncated//          DAT			Holds processed interpolation data
!		       interpolationfile
!               gridpath4//JPGfilename          PNG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine processInterpolationFile(&
                vals,delta_vals,&
                RMSDheatmap_coeff,&
                min_rsv,max_rsv,&
                min_Ninterpolation,max_Ninterpolation)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real(dp),dimension(Nvar),intent(in) :: vals, delta_vals

!NUMBER OF FRAMES IN THE DATA
real(dp),dimension(7),intent(out) :: min_rsv,max_rsv
real(dp),dimension(3),intent(out) :: RMSDheatmap_coeff
real(dp),dimension(7) :: rsv, rsv_binwidth
integer,dimension(7) :: rsv_bin

!FORMATTING OF PNG FILES
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5

!VARIABLES IN THE INTERPOLATION FILE
integer :: Ninterpolation
integer,intent(out) :: min_Ninterpolation
integer,intent(out) :: max_Ninterpolation
real(dp),dimension(Nvar) :: current_vals

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
integer :: Nbins
integer :: i,j,k,l,Nheatmap
real(dp),allocatable :: IEheatmap2D(:,:)
integer,allocatable :: ERheatmap2D(:,:),ERheatmap1D(:)
real(dp),allocatable :: A(:,:), b(:)

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

!The truncated interpolation file is more
!manageable
open(filechannel2,file=gridpath5//&
        "truncated"//interpolationfile)

!In particular, the first line describes the
!bounds of the data
read(filechannel2,FMT="(1x,3(I9,1x),2(I6,1x),"//&
        "2(7(G16.6,1x)))",iostat=iostate) &
        i,j,k,&
        min_Ninterpolation, max_Ninterpolation,&
        min_rsv, max_rsv

!The binning is done logarithmically and the
!number of bins is temporarily hardcoded
!The 2D heatmap looks better with higher
!resolution (compared to the 1D heatmap)
Nbins = 100
rsv_binwidth = log10(max_rsv/min_rsv)*1.0d0/Nbins
allocate(ERheatmap2D(Nbins,Nbins),&
         ERheatmap1D(Nbins/5),&
         IEheatmap2D(Nbins,Nbins))
ERheatmap2D = 0
ERheatmap1D = 0
IEheatmap2D = 1.0d-7

do 
    read(filechannel2,FMT=*,iostat=iostate) &
            current_vals, Ninterpolation, rsv
    if (iostate /= 0) exit

    !If the frame is not in the region of the phase
    !space of interest, then skip it
    if (any(abs(current_vals-vals) > &
            delta_vals)) cycle

    !Logarithms involve division by the minimum to
    !get a correct binning
    !The if statement is to make sure the edge cases
    !are counted in the distribution correctly
    rsv_bin = floor(log10(rsv/min_rsv)/rsv_binwidth) + 1
    do l = 1, 7
        if (rsv_bin(l) < 1) rsv_bin(l) = 1
        if (rsv_bin(l) > Nbins) rsv_bin(l) = Nbins
    end do

    !The variables chosen for the distributions are
    !not arbitrarily picked (but may be if wanted)
    ERheatmap2D(rsv_bin(5),rsv_bin(6)) =&
            ERheatmap2D(rsv_bin(5),rsv_bin(6)) + 1
    ERheatmap1D(rsv_bin(7)/5) =&
            ERheatmap1D(rsv_bin(7)/5) + 1

    !This last heatmap records the maximum of some
    !third value, not the total occurence (ideally
    !for maxmimum error checking)
    IEheatmap2D(rsv_bin(3),rsv_bin(4)) = max(&
            IEheatmap2D(rsv_bin(3),rsv_bin(4)),rsv(6))
end do
close(filechannel2)


!Then just bin the data normally

open(filechannel2,file=gridpath5//"ERheatmap2D.dat")
do i = 1, Nbins
    do j = 1, Nbins
        write(filechannel2,FMT=*)&
                (i-0.5)*rsv_binwidth(5),&
                (j-0.5)*rsv_binwidth(6),&
                ERheatmap2D(i,j)
    end do
    write(filechannel2,*) ""
end do
close(filechannel2)

open(filechannel2,file=gridpath5//"ERheatmap1D.dat")
do i = 1, Nbins/5
    write(filechannel2,FMT=*)&
            (i-0.5)*rsv_binwidth(7),&
            ERheatmap1D(i)
end do
close(filechannel2)

!For the third heatmap we are also interested in
!finding a fit for it, so we need to figure out
!how many nonzero entries are in it that we
!need to fit to

Nheatmap = 0
do i = 1, Nbins
    do j = 1, Nbins
        !We can see whether an entry is
        !nonzero by checking whether it is
        !greater than its initial value
        if (IEheatmap2D(i,j) > 1.0d-7) then
             Nheatmap = Nheatmap + 1
        end if
    end do
end do

!We store this data is matrices A and b
!A is the "independent data"
!b is the "dependent data"
!and we are trying to find a set of
!coefficients in matrix RMSDheatmap_coeff
!that, multiplied with A, produces a
!matrix close to b
allocate(A(Nheatmap,3),b(Nheatmap))
Nheatmap = 0
do i = 1, Nbins
    do j = 1, Nbins
        if (IEheatmap2D(i,j) > 1.0d-7) then
             Nheatmap = Nheatmap + 1

             !In this case, our fit is
             ! c2 * log10(rsv3) +
             !   c1 * log10(rsv4) +
             !                   c0 = log10(b)
             A(Nheatmap,:) = (/&
                     (i-0.5)*rsv_binwidth(3),&
                     (j-0.5)*rsv_binwidth(4),&
                     1.0d0/)
             b(Nheatmap) = log10(IEheatmap2D(i,j))
        end if
    end do
end do

!This fit is linear so can be obtained by doing
!a least squares regression on the above matrices
call LS(A,Nheatmap,3,b,RMSDheatmap_coeff)


!And then just plot the two sets of data (the
!interpolation and the fit) side by side

open(filechannel2,file=gridpath5//"IEheatmap2D.dat")
do i = 1, Nbins
    do j = 1, Nbins
        write(filechannel2,FMT=*)&
                (i-0.5)*rsv_binwidth(3),&
                (j-0.5)*rsv_binwidth(4),&
                IEheatmap2D(i,j),&
                10.0d0**(max((i-0.5)*rsv_binwidth(3)*RMSDheatmap_coeff(1)+&
                             (j-0.5)*rsv_binwidth(4)*RMSDheatmap_coeff(2)+&
                                            RMSDheatmap_coeff(3),-7.0d0))
    end do
    write(filechannel2,*) ""
end do
close(filechannel2)

end subroutine processInterpolationFile



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               processInterpolationFile2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at the interpolation file, gets minimums and
!               maximums, generates secondary variables, and writes the data down in
!               a separate truncated file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               vals                            REAL,DIM(Nvar)                  The region of the phase space
!                                                                               (described by collective variables)
!               rsv                             REAL,DIM(7)                     The variables in order are:
!                                                                               1.RMSD of accept best frame
!                                                                               2.RMSD of interpolated frame
!                                                                               3.Special Variable 1
!                                                                               4.Special Variable 2
!                                                                               5.Error of accept best gradient
!                                                                               6.Error of interpolated gradient
!                                                                               7.Error of accept best / interpolated
!               min_rsv                         REAL,DIM(7)                     Minimums of the RMSD special variables
!               max_rsv                         REAL,DIM(7)                     Maximums of the RMSD special variables
!               Ninterpolation                  INTEGER                         Number of frames used in the
!                                                                               interpolation
!               min_Ninterpolation              INTEGER                         Minimum of the Ninterpolation
!               max_Ninterpolation              INTEGER                         Maximum of the Ninterpolation
!               frames                          INTEGER                         Number of frames in the interpolation
!                                                                               file
!               neginfinity_counter             INTEGER                         Number of frames where RSV5 = 0
!               posinfinity_counter             INTEGER                         Number of frames where RSV6 = 0
!               firstliner                      CHARACTER()                     A string with minimums and maximums
!                                                                               of the RSV and other data appended as
!                                                                               the first line of the truncated file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath5//interpolationfile    DAT			Holds the raw interpolation data
!		gridpath5//truncated//          DAT			Holds processed interpolation data
!		       interpolationfile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine processInterpolationFile2()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
character(15) :: vals_interpolation_text
character(12) :: short_prefix_text

!VARIABLES IN THE INTERPOLATION FILE
real(dp),dimension(Nvar) :: vals
real(dp),dimension(9) :: rsv, min_rsv, max_rsv
integer :: min_Ninterpolation,max_Ninterpolation,Ninterpolation
integer :: Nswitch

integer :: frames,neginfinity_counter,posinfinity_counter
character(1+3*10+2*7+2*9*17) :: firstliner

!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

min_Ninterpolation = 1000
max_Ninterpolation = 0

min_rsv = 1.0d9
max_rsv = 0.0d0

frames = 0
neginfinity_counter = 0
posinfinity_counter = 0

!The data we read from the interpolation file, we do some
!minor processing and then record to the truncated interpolation
!file
open(filechannel1,file=gridpath5//"truncated"//interpolationfile)
open(filechannel2,file=gridpath5//interpolationfile)
do 
    !For each frame in the file, the data read is:
    !1.Its collective variables
    !2.Number of frames used for its interpolation
    !3.The RMSD special variables that results from it
    read(filechannel2,FMT="(2(F10.6,1x),I5,1x,8(ES10.2))",iostat=iostate) &
            vals, Ninterpolation, &
!           rsv1, rsv2, &
!           rmsd_best, rmsd_interpolated, &
!           min_CMdiff, interpolated_CMdiff, &
!           error_best, error_interpolated
            rsv(3), rsv(4), &
            rsv(1), rsv(2), &
            rsv(8), rsv(9), &
            rsv(5), rsv(6)
            
    if (iostate /= 0) exit

!   rsv(5) = rsv(5) * RU_energy / eV
!   rsv(6) = rsv(6) * RU_energy / eV

    !One of the processed variables that this
    !subroutine makes is the ratio between
    !RSV(5) and RSV(6) which means a few
    !conditionals need to be made to treat them

    !This may bias the data so these occurences
    !are recorded as negative and positive
    !infinity (since we will later take a
    !logarithm

    if ((rsv(5) == 0.0d0).and.(rsv(6)==0.0d0)) then
        rsv(7) = 1.0d0
        write(filechannel1,FMT="(2(F10.6,1x),I5,1x,9(ES10.2))") &
                vals, Ninterpolation, rsv

        min_Ninterpolation = min(min_Ninterpolation,&
                Ninterpolation)
        max_Ninterpolation = max(max_Ninterpolation,&
                Ninterpolation)

        do n = 1, 9
            max_rsv(n) = max(max_rsv(n),rsv(n))

            if (rsv(n) == 0.0d0) then
            else
                min_rsv(n) = min(min_rsv(n),rsv(n))
            end if
        end do
        
    else if (rsv(5) == 0.0d0) then
        neginfinity_counter = neginfinity_counter + 1

    else if (rsv(6) == 0.0d0) then
        posinfinity_counter = posinfinity_counter + 1

    else
        rsv(7) = rsv(5)/rsv(6)
        write(filechannel1,FMT="(2(F10.6,1x),I5,1x,9(ES10.2))") &
                vals, Ninterpolation, rsv

        min_Ninterpolation = min(min_Ninterpolation,&
                Ninterpolation)
        max_Ninterpolation = max(max_Ninterpolation,&
                Ninterpolation)

        do n = 1, 9
            max_rsv(n) = max(max_rsv(n),rsv(n))

            if (rsv(n) == 0.0d0) then
            else
                min_rsv(n) = min(min_rsv(n),rsv(n))
            end if
        end do
    
    end if

    !This subroutine also counts the total number of
    !frames in the file
    frames = frames + 1
end do
close(filechannel1)
close(filechannel2)

!All of this metadata on the bounds of the data,
!the positive/negative infinity bias, and the
!total number of frames is recorded in the first
!line of the file
write(firstliner,FMT="(A1,3(I9,1x),2(I6,1x),"//&
        "2(9(E16.6,1x)))") &
        "#",posinfinity_counter,&
        neginfinity_counter, frames, &
        min_Ninterpolation, max_Ninterpolation,&
        min_rsv, max_rsv

!Sed is used to deposit the first line onto the
!file
call system("sed -i '1i\"//firstliner//"' "//&
        gridpath5//"truncated"//&
        interpolationfile)

end subroutine processInterpolationFile2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               processMultipleInterpolationFiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at multiple processed interpolation files (truncated
!               files), calculates lower and upper bounds, bins the data, and writes down the
!               data in separate files for later plotting
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               counters                        INTEGER                         Number of frames and other counters
!                  (comparison_number,3)                                        in each interpolation file for
!                                                                               comparison
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               comparison_number               INTEGER                         The number of experiments to be
!                                                                               compared (not necessarily unique)
!               vals                            REAL,DIM(Nvar)                  The region of the phase space
!                                                                               (described by collective variables)
!               delta_vals                      REAL,DIM(Nvar)                  The size of the region of the phase
!                                                                               space of interest
!               rsv                             REAL,DIM(7)                     The variables in order are:
!                                                                               1.RMSD of accept best frame
!                                                                               2.RMSD of interpolated frame
!                                                                               3.Special Variable 1
!                                                                               4.Special Variable 2
!                                                                               5.Error of accept best gradient
!                                                                               6.Error of interpolated gradient
!                                                                               7.Error of accept best / interpolated
!               min_rsv                         REAL,DIM(7)                     Minimums of the RMSD special variables
!               max_rsv                         REAL,DIM(7)                     Maximums of the RMSD special variables
!               rsv_binwidth                    REAL,DIM(7)                     Width of the binning for each RSV
!               rsv_bin                         INTEGER,DIM(7)                  Temporary buffer for each RSV bin data
!               Ninterpolation                  INTEGER                         Number of frames used in the
!                                                                               interpolation
!               min_Ninterpolation              INTEGER                         Minimum of the Ninterpolation
!               max_Ninterpolation              INTEGER                         Maximum of the Ninterpolation
!               Ninterpolation_binwidth         REAL                            Width of the binning for
!                                                                               Ninterpolation
!               Ninterpolation_bin              INTEGER                         Temporary buffer for Ninterpolation
!                                                                               bin data
!               frames                          INTEGER                         Number of frames in the interpolation
!                                                                               file
!               frames_trials                   INTEGER                         Number of frames in each of the
!                  (comparison_number)                                          interpolation file for comparison
!               rsv_binning                     REAL,DIM(7,                     The RSV data binned (scaled so
!                                                 3*comparison_number,          to integrate to 1)
!                                                 Nbins)
!               Ninterpolation_binning          REAL,DIM(                       The Ninterpolation data binned
!                                                 3*comparison_number,          (scaled so to integrate to 1)
!                                                 Nbins)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath0//comparison_file      DAT			Holds the parameters that describe
!                                                                       what regions of the phase space to
!                                                                       excise data from
!		gridpath0//allprefixes()//      DAT			Holds processed interpolation data
!		       intermediatefolder//
!		       interpolationfile
!		gridpath5//tmp//                DAT			Holds slightly more processed
!		       interpolationfile                                interpolation data
!		gridpath5//SATRVname//          DAT			Holds binned interpolation data
!		       _binning                                         for later binning
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine processMultipleInterpolationFiles(counters)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real(dp),dimension(Nvar) :: vals, delta_vals, tmpvals
character(expfolder_length-1) :: otherexpfolder
character(70) :: longtext

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

integer,dimension(comparison_number,3),intent(out) :: counters
integer,dimension(comparison_number) :: frames_trials

real(dp),dimension(9) :: rsv, min_rsv, max_rsv
real(dp),dimension(comparison_number,9) :: lower_rsv, upper_rsv
real(dp),dimension(comparison_number,9) :: arithmetic_mean_rsv,geometric_mean_rsv
integer :: min_Ninterpolation,max_Ninterpolation,Ninterpolation
integer,dimension(comparison_number) :: lower_Ninterpolation,upper_Ninterpolation
real(dp),dimension(comparison_number) :: arithmetic_mean_Ninterpolation
real(dp),dimension(comparison_number) :: geometric_mean_Ninterpolation

!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
real(dp) :: Ninterpolation_binwidth
integer :: Ninterpolation_bin
real(dp), dimension(9) :: rsv_binwidth
integer, dimension(9) :: rsv_bin
integer :: frames, Nbins
real(dp),allocatable :: Ninterpolation_binning(:,:)
real(dp),allocatable :: rsv_binning(:,:,:)
integer :: starting_index
real(dp),allocatable :: zero_array(:)

real(dp),dimension(Nvar) :: vals_prev
integer,dimension(Nvar) :: deltavals
integer :: Nswitch
integer :: beyond_counter,ten_counter
integer,dimension(3) :: pos_counter,neg_counter

logical :: debug_flag = .false.

!INTEGER INCREMENTALS
integer :: n,m,l

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

counters = 0
frames_trials = 0
frames = 0

min_Ninterpolation = 1000000
max_Ninterpolation = 0
min_rsv = 1.0d19
max_rsv = 0.0d0

open(filechannel1,file=gridpath0//comparison_file)
open(filechannel3,file=gridpath5//"tmp"//interpolationfile)

!We open a file once per comparison number (even if
!it's a duplicate)
do n = 1, comparison_number

    !The comparison file has information regarding
    !what region of the phase space we want to plot
    read(filechannel1,FMT=*,iostat=iostate) vals, delta_vals
    if (iostate /= 0) then
        !If the amount of data on the comparison file does]
        !not match the comparison number than something
        !went wrong
        print *, ""
        print *, "Bad formatting in interpolation comparison arguments"
        print *, "   Must be 2*Nvar floats"
        print *, ""

        close(filechannel1)
        close(filechannel3)

        return
    end if

    !The string allprefixes has information regarding
    !what experiment we want to plot and the array
    !alllengths has information regarding how long
    !each string in it is
    if (n == 1) then
        starting_index = 0
    else
        starting_index = sum(alllengths(1:n-1))
    end if
    open(filechannel2,file=gridpath0//&
            allprefixes(starting_index+1:sum(alllengths(1:n)))//&
            intermediatefolder//"truncated"//interpolationfile)

    !The first line in the truncated interpolation file
    !has the bounds to the data
    read(filechannel2,FMT="(1x,3(I9,1x),2(I6,1x),"//&
            "2(9(G16.6,1x)))",iostat=iostate) &
            counters(n,1),counters(n,2),counters(n,3),&
            lower_Ninterpolation(n), upper_Ninterpolation(n),&
            lower_rsv(n,:),&
            upper_rsv(n,:)
!           (lower_rsv(n,i),i=1,7),&
!           (upper_rsv(n,i),i=1,7),&
!           (arithmetic_mean_rsv(n,i),i=1,7),&
!           (geometric_mean_rsv(n,i),i=1,7)

    pos_counter = 0
    neg_counter = 0
    beyond_counter = 0
    ten_counter = 0

    min_Ninterpolation = min(min_Ninterpolation,&
            lower_Ninterpolation(n))
    max_Ninterpolation = max(max_Ninterpolation,&
            upper_Ninterpolation(n))
    
    do m = 1, 9
        min_rsv(m) = min(min_rsv(m),lower_rsv(n,m))
        max_rsv(m) = max(max_rsv(m),upper_rsv(n,m))
    end do

    upper_Ninterpolation(n) = min_Ninterpolation
    upper_rsv(n,:) = min_rsv

    !We also excise any data that is not in the region
    !of the phase space of interest
    !All data is stored in the temporary interpolation
    !file; yes, there will be duplicates if the same
    !data is encountered across two trials
    do
        read(filechannel2,FMT="(2(F10.6,1x),I5,1x,9(ES10.2))",iostat=iostate)&
                tmpvals, Ninterpolation, rsv
        if (iostate /= 0) exit

        if (frames_trials(n) > 0) then
    
        deltavals = floor(abs(tmpvals-vals_prev)*divisor(:,1))
        deltavals = abs(floor(tmpvals*divisor(:,2)) - &
                        floor(vals_prev*divisor(:,2)))
        Nswitch = sum(deltavals)+1
    
        if ((Ninterpolation <= 9)) then
            ten_counter = ten_counter + 1
        else
!       if ((Nswitch < 4).and.(Ninterpolation > 9)) then
        if ((Nswitch < 4)) then
        if (rsv(7) > 1.0d0) then
            pos_counter(Nswitch) = pos_counter(Nswitch) + 1
        else
            neg_counter(Nswitch) = neg_counter(Nswitch) + 1
        end if
        else
            beyond_counter = beyond_counter + 1
        end if
        end if
    
        end if
    
        vals_prev = tmpvals

        tmpvals = abs(tmpvals-vals) - abs(delta_vals)

        if (all(delta_vals>0)) then
            if (any(tmpvals > 0)) cycle
        else
            if (all(tmpvals <= 0)) cycle
        end if

        upper_Ninterpolation(n) = max(upper_Ninterpolation(n),&
                Ninterpolation)
        do m = 1, 9
            upper_rsv(n,m) = max(upper_rsv(n,m),rsv(m))
        end do

        !Recording the number of frames in each
        !trial lets us know where each trial begins
        !and ends on this one big file
        frames_trials(n) = frames_trials(n) + 1
        write(filechannel3,FMT="(I5,1x,9(ES10.2))") Ninterpolation,rsv
    end do
    close(filechannel2)

    if (debug_flag) then
        print *, ""
        print *, allprefixes(starting_index+1:sum(alllengths(1:n)))
        do Nswitch = 1, 3
        print *, "n = ", Nswitch - 1, " neg / pos:", &
                neg_counter(Nswitch), "/", pos_counter(Nswitch), "%neg/total:", &
                neg_counter(Nswitch) *100.0 / (neg_counter(Nswitch) + pos_counter(Nswitch)), "%"
        end do
        print *, ""
        print *, "n > 2:", beyond_counter
        print *, ""
        print *, "Ninterpolation <= 9:", ten_counter
        print *, ""
    end if
end do

close(filechannel1)
close(filechannel3)

!We have to get the bounds to ALL the data
!across all experiments

!min_Ninterpolation = minval(lower_Ninterpolation)
!max_Ninterpolation = maxval(upper_Ninterpolation)
!
!do m = 1, 9
!    min_rsv(m) = minval(lower_rsv(:,m))
!    max_rsv(m) = maxval(upper_rsv(:,m))
!end do

print *, ""
print *, "comparsion_number:", comparison_number
print *, "frames_trials:", frames_trials
print *, ""

!If any of the trials has no frames, there will
!probably be trouble plotting it
if (any(frames_trials == 0)) then
    print *, "MAJOR ERROR IN TRIALS SELECTED FOR "//&
             "INTERPOLATION ANALYSIS"
    frames_trials = 0
    return
end if

!The number of bins is dynamically chosen so that
!every bin (even at the beginning and end) has an
!integer number of Ninterpolation
!This means there is some fiddling with how the
!range of the data is factored
Nbins = max_Ninterpolation - min_Ninterpolation
do
    if (Nbins < 61) then
        exit
    else if (modulo(Nbins,2) == 0) then
        Nbins = Nbins / 2
    else if (modulo(Nbins,3) == 0) then
        Nbins = Nbins / 3
    else if (modulo(Nbins,5) == 0) then
        Nbins = Nbins / 5
    else
        exit
    end if
end do

!These fallbacks exist so that this dynamic binning
!doesn't get out of hand (but then each bin may not
!have an integer number of Ninterpolation)
if (Nbins == 0) Nbins = 50
if (Nbins > 100) Nbins = 100

!The binwidths are determined from the bounds of
!the data and the number of bins designated to it
Ninterpolation_binwidth = ceiling((max_Ninterpolation -&
        min_Ninterpolation) *1.0d0/ Nbins)
rsv_binwidth = log10(max_rsv/min_rsv)*1.0d0/Nbins

!There is a multiplication by three because three
!subsets of each distribution are differentiated
allocate(Ninterpolation_binning(3*comparison_number,Nbins),&
         rsv_binning(9,3*comparison_number,Nbins),&
         zero_array(3*comparison_number))
zero_array = 0.0d0
Ninterpolation_binning = 0.0d0
rsv_binning = 0.0d0

!do n = 1, comparison_number
!    upper_Ninterpolation(n) = min_Ninterpolation
!    upper_rsv(n,:) = min_rsv
!end do

arithmetic_mean_Ninterpolation = 0.0d0
geometric_mean_Ninterpolation = 0.0d0

arithmetic_mean_rsv = 0.0d0
geometric_mean_rsv = 0.0d0

!The data for ALL trials is stored in the temporary
!interpolation file
open(filechannel1,file=gridpath5//"tmp"//interpolationfile)
do n = 1, comparison_number

    !Each trial's data set is not separated so the
    !number of lines of data read must be counted
    do m = 1, frames_trials(n)
        read(filechannel1,FMT="(I5,1x,9(ES10.2))") Ninterpolation,rsv

!        upper_Ninterpolation(n) = max(upper_Ninterpolation(n),&
!                Ninterpolation)
        arithmetic_mean_Ninterpolation(n) = &
            arithmetic_mean_Ninterpolation(n) + Ninterpolation
        geometric_mean_Ninterpolation(n) = &
            geometric_mean_Ninterpolation(n) + log10(1.0d0*Ninterpolation)

        do l = 1, 9
!            upper_rsv(n,l) = max(upper_rsv(n,l),rsv(l))
            arithmetic_mean_rsv(n,l) = &
                arithmetic_mean_rsv(n,l) + rsv(l)

            if (rsv(l) == 0.0d0) then
                geometric_mean_rsv(n,l) = &
                    geometric_mean_rsv(n,l) + log10(min_rsv(n))
            else
                geometric_mean_rsv(n,l) = &
                    geometric_mean_rsv(n,l) + log10(rsv(l))
            end if
        end do
        

        !Read each line and bin them accordingly
    
        Ninterpolation_bin = floor((Ninterpolation&
                - min_Ninterpolation)&
                / Ninterpolation_binwidth) + 1
        if (Ninterpolation_bin < 1) Ninterpolation_bin = 1
        if (Ninterpolation_bin > Nbins) Ninterpolation_bin = Nbins

        rsv_bin = floor(log10(rsv/min_rsv)/rsv_binwidth) + 1
        do l = 1, 9
            if (rsv_bin(l) < 1) rsv_bin(l) = 1
            if (rsv_bin(l) > Nbins) rsv_bin(l) = Nbins
        end do

        !The first portion of the data is ALL of
        !the data; this will be in the back of
        !all other data
        if (.true.) then
            Ninterpolation_binning(3*(n-1)+1,Ninterpolation_bin) = &
                    Ninterpolation_binning(3*(n-1)+1,Ninterpolation_bin) + &
                    1.0d0/frames_trials(n)
            do l = 1, 9
                rsv_binning(l,3*(n-1)+1,rsv_bin(l)) =&
                        rsv_binning(l,3*(n-1)+1,rsv_bin(l)) +&
                        1.0d0/frames_trials(n)
            end do
        end if

        !The second portion of the data is that
        !which has this quality (can be
        !arbitrarily picked)
        if (rsv(7) < 1.0d1) then
!       if (rsv(6) <= upper_rsv(n,5)) then
            Ninterpolation_binning(3*(n-1)+2,Ninterpolation_bin) = &
                    Ninterpolation_binning(3*(n-1)+2,Ninterpolation_bin) + &
                    1.0d0/frames_trials(n)
            do l = 1, 9
                rsv_binning(l,3*(n-1)+2,rsv_bin(l)) =&
                        rsv_binning(l,3*(n-1)+2,rsv_bin(l)) +&
                        1.0d0/frames_trials(n)
            end do
        end if

        !The third portion of the data is that
        !which has this quality (can be
        !arbitrarily picked); this will be in
        !the front of all other data
        if (rsv(7) <= 1.0d0) then
!       if (rsv(6) <= 0.5d0 * upper_rsv(n,5)) then
            Ninterpolation_binning(3*n,Ninterpolation_bin) = &
                    Ninterpolation_binning(3*n,Ninterpolation_bin) + &
                    1.0d0/frames_trials(n)
            do l = 1, 9
                rsv_binning(l,3*n,rsv_bin(l)) =&
                        rsv_binning(l,3*n,rsv_bin(l)) +&
                        1.0d0/frames_trials(n)
            end do
        end if
    end do
end do
close(filechannel1)

do n = 1, comparison_number
    arithmetic_mean_rsv(n,:) = &
        arithmetic_mean_rsv(n,:) / frames_trials(n)
    geometric_mean_rsv(n,:) = &
        10.0d0 ** (geometric_mean_rsv(n,:)  / frames_trials(n))
    
    arithmetic_mean_Ninterpolation(n) = &
        arithmetic_mean_Ninterpolation(n) * 1.0d0 / &
        frames_trials(n)
    geometric_mean_Ninterpolation = &
        10.0d0 ** (geometric_mean_Ninterpolation / &
        frames_trials(n))
end do

!Finally, record the data for
!later plotting

open(filechannel1,file=gridpath5//"InterpolationTDD_binning.dat")

if (comparison_upperlimit > comparison_lowerlimit) then
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            comparison_lowerlimit, comparison_upperlimit
else
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            min_Ninterpolation - 0.5*Ninterpolation_binwidth, &
            max_Ninterpolation + 0.5*Ninterpolation_binwidth
end if

do n = 1, comparison_number
write(filechannel1,FMT="(A1,3(E16.6,1x))") &
        "#",1.0d0*upper_Ninterpolation(n), &
        arithmetic_mean_Ninterpolation(n),&
        geometric_mean_Ninterpolation(n)
end do

write(filechannel1,FMT="(4(E16.6,1x))") min_Ninterpolation + &
        Ninterpolation_binwidth*(-0.5), &
        zero_array
do n = 1, Nbins
    write(filechannel1,FMT="(4(E16.6,1x))") min_Ninterpolation + &
            Ninterpolation_binwidth*(n-0.5),&
            Ninterpolation_binning(:,n)
end do
write(filechannel1,FMT="(4(E16.6,1x))") min_Ninterpolation + &
        Ninterpolation_binwidth*(Nbins+0.5), &
        zero_array
close(filechannel1)

do l = 1, 9

!Accept Best RMSD Distribution
if (l == 1) open(filechannel1,file=gridpath5//&
        "InterpolationARD_binning.dat")
!Interpolation RMSD Distribution
if (l == 2) open(filechannel1,file=gridpath5//&
        "InterpolationIRD_binning.dat")
!RMSD Special Variable 1 Distribution
if (l == 3) open(filechannel1,file=gridpath5//&
        "InterpolationRSV1D_binning.dat")
!RMSD Special Variable 2 Distribution
if (l == 4) open(filechannel1,file=gridpath5//&
        "InterpolationRSV2D_binning.dat")
!Accept Best Error Distribution
if (l == 5) open(filechannel1,file=gridpath5//&
        "InterpolationAED_binning.dat")
!Interpolation Error Distribution
if (l == 6) open(filechannel1,file=gridpath5//&
        "InterpolationIED_binning.dat")
!Relative Error Distribution
if (l == 7) open(filechannel1,file=gridpath5//&
        "InterpolationRED_binning.dat")
!Accept Best Coulomb Matrix Difference Distribution
if (l == 8) open(filechannel1,file=gridpath5//&
        "InterpolationACMDD_binning.dat")
!Interpolation Coulomb Matrix Difference Distribution
if (l == 9) open(filechannel1,file=gridpath5//&
        "InterpolationICMDD_binning.dat")

if (comparison_upperlimit > comparison_lowerlimit) then
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            comparison_lowerlimit, comparison_upperlimit
else
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            log10(min_rsv(l)) - 0.5*rsv_binwidth(l), &
            log10(max_rsv(l)) + 0.5*rsv_binwidth(l)
end if

do n = 1, comparison_number
write(filechannel1,FMT="(A1,3(E16.6,1x))") &
        "#",log10(upper_rsv(n,l)), &
        log10(arithmetic_mean_rsv(n,l)),&
        log10(geometric_mean_rsv(n,l))
end do

write(filechannel1,FMT="(4(E16.6,1x))")&
        log10(min_rsv(l)) +&
        rsv_binwidth(l)*(-0.5d0),&
        zero_array
do n = 1, Nbins
    write(filechannel1,FMT="(4(E16.6,1x))")&
            log10(min_rsv(l)) +&
            rsv_binwidth(l)*(n-0.5d0),&
            rsv_binning(l,:,n)
end do
write(filechannel1,FMT="(4(E16.6,1x))")&
        log10(min_rsv(l)) +&
        rsv_binwidth(l)*(Nbins+0.5d0),&
        zero_array

close(filechannel1)
end do

deallocate(Ninterpolation_binning,rsv_binning,zero_array)

end subroutine processMultipleInterpolationFiles


subroutine plotInterpolationOccurenceHeatmap()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real,dimension(Nvar) :: tmpvals,var_minvar
integer,dimension(Nvar) :: max_val_bin, val_bins
integer :: max_bin
logical :: cycle_flag

real(dp),dimension(9) :: rsv, min_rsv, max_rsv
integer :: min_Ninterpolation,max_Ninterpolation,Ninterpolation

!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
real(dp),allocatable :: heatmap_binning(:,:)

!INTEGER INCREMENTALS
integer :: n,i,j

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

var_minvar = (/ 2.0, 2.0 /)

max_val_bin = int((var_maxvar - var_minvar)/var_spacing)

allocate(heatmap_binning(max_val_bin(1),&
                         max_val_bin(2)))
heatmap_binning = 0

open(filechannel1,file=gridpath5//&
        "truncated"//interpolationfile)
read(filechannel1,FMT="(1x,3(I9,1x),2(I6,1x),"//&
        "2(9(G16.6,1x)))",iostat=iostate) &
        i, j, n,&
        min_Ninterpolation, max_Ninterpolation,&
        min_rsv,max_rsv

do
    read(filechannel1,FMT=*,iostat=iostate)&
            tmpvals, Ninterpolation, rsv
    if (iostate /= 0) exit

    !Here we define what is NOT an occurence
    if (rsv(6) <= max_rsv(5)) cycle

    val_bins = floor((tmpvals-var_minvar)/var_spacing)

    cycle_flag = .false.
    do n = 1, Nvar
        if (val_bins(n) <= 0)&
            cycle_flag = .true.
        if (val_bins(n) > max_val_bin(n))&
            cycle_flag = .true.
    end do
    
    if (cycle_flag) cycle

    if (cycle_flag) cycle

    heatmap_binning(val_bins(1),val_bins(2)) = &
        heatmap_binning(val_bins(1),val_bins(2)) + 1

end do
close(filechannel1)

max_bin = maxval(heatmap_binning)

open(filechannel1,file=gridpath5//temporaryfile1)
do i = 1, max_val_bin(1)
    tmpvals(1) = var_minvar(1) + (i-1)*var_spacing(1)
    
    do j = 1, max_val_bin(2)
        tmpvals(2) = var_minvar(2) + (j-1)*var_spacing(2)
        
        write(filechannel1,FMT=*) tmpvals, &
                heatmap_binning(i,j)
    end do
    
    write(filechannel1,FMT=*) "" 
end do
close(filechannel1)

deallocate(heatmap_binning)



open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set terminal pngcairo size 1800,1800'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//'InterpolationOccurenceHeatMap.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT='(A,4(I0.5,A))') &
        'set palette defined (',&
        0, ' "green", ',&
        1, ' "blue", ',&
        max(max_bin/2,2), ' "yellow", ',&
        max(max_bin,3), ' "red")'
write(gnuplotchannel,FMT='(A,I0.5,A)') 'set cbrange [0:',max_bin,']'
write(gnuplotchannel,FMT="(A)") 'set cblabel "Number of Occurences" font ",18" offset 1,0'

write(gnuplotchannel,FMT="(A)") 'set title "Interpolation Occurence Heatmap of an H_2 - H_2 System" font ",32" offset 0,3'
write(gnuplotchannel,*) 'set xlabel "Var1 (A)" font ",28" offset 0,-2'
write(gnuplotchannel,*) 'set xtics 1 font ",24"'
write(gnuplotchannel,*) 'set xrange [', var_minvar(1), ':', var_maxvar(1),']'
write(gnuplotchannel,*) 'set ylabel "Var2 (A)" font ",28" offset -5,0'
write(gnuplotchannel,*) 'set ytics 1 font ",24"'
write(gnuplotchannel,*) 'set yrange [', var_minvar(2), ':',var_maxvar(2),']'
write(gnuplotchannel,*) 'set cbtics'
write(gnuplotchannel,*) 'set view map'
write(gnuplotchannel,*) 'set pm3d interpolate 1,1'
write(gnuplotchannel,*) 'splot "'//gridpath5//temporaryfile1//'" u 1:2:3 w pm3d'
close(gnuplotchannel)
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine plotInterpolationOccurenceHeatmap



subroutine getRMSDinterpolation2(PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

integer, dimension(comparison_number,3) :: counters

integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

call processMultipleInterpolationFiles(counters)

if (any(counters(:,3) == 0)) then
    return
end if

open(gnuplotchannel,file=gridpath5//PNGfilename//"_"//gnuplotfile)
write(gnuplotchannel,FMT="(A,I6)") 'set term pngcairo size 2400,',1200*comparison_number
write(gnuplotchannel,*) 'set output ARG1."/../'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,FMT="(A,I0.2,A)") 'set multiplot layout ',comparison_number,&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"Interpolation Error with Varying Threshold" font ",36" offset 0,3'
!write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
!        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.3,0.950'
!write(gnuplotchannel,FMT="(A,I6,A)") &
!        'set label 2 "Ntraj = ', Ntraj_max, '" at screen 0.3,0.940'
!write(gnuplotchannel,FMT='(A,F9.4,A)') 'set label 3 "AlphaRatio = ',alpha_ratio, &
!        '" at screen 0.3,0.930'
write(gnuplotchannel,*) 'set ylabel "Frequency" font ",18"'
write(gnuplotchannel,*) 'set xtics nomirror'
write(gnuplotchannel,*) 'set grid xtics lw 2'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set format x ""'

if (comparison_number > 1) then
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '# Edit to title this plot (01)'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") &
        'plot ARG1."/'//&
        'InterpolationTDD_binning.dat" u '//'1:',2,' w boxes t \'
write(gnuplotchannel,FMT="(A)") &
        '" "\'
write(gnuplotchannel,FMT="(A,F9.6,A)") &
        ' fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

do n = 2, comparison_number-1

write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',n,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") &
        'plot ARG1."/'//&
        'InterpolationTDD_binning.dat" u '//'1:',n+1,' w boxes t \'
write(gnuplotchannel,FMT="(A)") &
        '" "\'
write(gnuplotchannel,FMT="(A,F9.6,A)") &
        ' fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

end do
end if

write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set xtics out nomirror'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Number of Points Below Threshold (Ninterpolation)" font ",24"'

write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',comparison_number,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") &
        'plot ARG1."/'//&
        'InterpolationTDD_binning.dat" u '//'1:',comparison_number+1,' w boxes t \'
write(gnuplotchannel,FMT="(A)") &
        '" "\'
write(gnuplotchannel,FMT="(A,F9.6,A)") &
        ' fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set xtics nomirror'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                         '"" log10(.00000001), '//&
                                         '"" log10(.00000005), '//&
                                          '"" log10(.0000001), '//&
                                          '"" log10(.0000005), '//&
                                           '"" log10(.000001), '//&
                                           '"" log10(.000005), '//&
                                           '""  log10(.00001), '//&
                                           '""  log10(.00005), '//&
                                           '""   log10(.0001), '//&
                                           '""   log10(.0005), '//&
                                           '""    log10(.001), '//&
                                           '""    log10(.005), '//&
                                           '""     log10(.01), '//&
                                           '""     log10(.05), '//&
                                           '""      log10(.1), '//&
                                           '""      log10(.5), '//&
                                           ' ""       log10(1), '//&
                                           ' ""       log10(5), '//&
                                           ' ""      log10(10), '//&
                                           ' ""      log10(50), '//&
                                           ' ""     log10(100), '//&
                                           ' ""     log10(500), '//&
                                           ' ""    log10(1000), '//&
                                           ' ""    log10(5000), '//&
                                           ' ""   log10(10000), '//&
                                   ')'

if (comparison_number > 1) then
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'

do n = 2, comparison_number-1

write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*(n-1)+2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*(n-1)+3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*(n-1)+4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'

end do
end if

write(gnuplotchannel,*) 'set xtics out nomirror'
write(gnuplotchannel,FMT="(A)") 'set xtics ('//&
                                         '"1e-8" log10(.00000001), '//&
                                         '"5e-8" log10(.00000005), '//&
                                          '"1e-7" log10(.0000001), '//&
                                          '"5e-7" log10(.0000005), '//&
                                           '"1e-6" log10(.000001), '//&
                                           '"5e-6" log10(.000005), '//&
                                           '"1e-5"  log10(.00001), '//&
                                           '"5e-5"  log10(.00005), '//&
                                           '"1e-4"   log10(.0001), '//&
                                           '"5e-4"   log10(.0005), '//&
                                           '"1e-3"    log10(.001), '//&
                                           '"5e-3"    log10(.005), '//&
                                           '"1e-2"     log10(.01), '//&
                                           '"5e-2"     log10(.05), '//&
                                           '"1e-1"      log10(.1), '//&
                                           '"5e-1"      log10(.5), '//&
                                           ' "1e0"       log10(1), '//&
                                           ' "5e0"       log10(5), '//&
                                           ' "1e1"      log10(10), '//&
                                           ' "5e1"      log10(50), '//&
                                           ' "1e2"     log10(100), '//&
                                           ' "5e2"     log10(500), '//&
                                           ' "1e3"    log10(1000), '//&
                                           ' "5e3"    log10(5000), '//&
                                           ' "1e4"   log10(10000), '//&
                                   ')'

write(gnuplotchannel,*) 'set xlabel "Relative Error of Accept Best to Interpolation" font ",24"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*comparison_number-1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*comparison_number,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*comparison_number+1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot -c "//&
        gridpath5//PNGfilename//"_"//gnuplotfile//" "//&
        '"'//gridpath5(1:gridpath_length+expfolder_length+5-1)//'"')

end subroutine getRMSDinterpolation2


subroutine getInterpolationplot(PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!FORMAT OF PNG FILES TO BE MADE
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

integer, dimension(comparison_number,3) :: counters
real(dp),dimension(comparison_number) :: var_max,var_Amean,var_Gmean
real(dp) :: lowerlimit, upperlimit

integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

call processMultipleInterpolationFiles(counters)

if (all(counters == 0)) then
    return
end if

open(filechannel1,file=gridpath5//&
        trim(adjustl(comparison_SATRVname))//&
        "_binning.dat")
read(filechannel1,FMT="(1x,2(E16.6,1x))") &
    lowerlimit, upperlimit

do n = 1, comparison_number
    read(filechannel1,FMT="(1x,3(E16.6,1x))") &
        var_max(n), var_Amean(n), var_Gmean(n)
end do
close(filechannel1)

open(gnuplotchannel,file=gridpath5//gnuplotfile//"_"//&
        trim(adjustl(comparison_SATRVname)))
write(gnuplotchannel,FMT="(A,I6)") 'set term pngcairo size 1200,',400*comparison_number
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output ARG1."/../'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'PL=strlen(ARG1)'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,FMT="(A,I0.2,A)") 'set multiplot layout ',comparison_number,&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"" font ",32" offset 0,-3'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
!write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set format y ""'
write(gnuplotchannel,*) 'set xtics nomirror'
write(gnuplotchannel,*) 'set grid xtics lw 2'

write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set xtics nomirror font ",16"'

if (trim(adjustl(comparison_SATRVname)) /= 'InterpolationTDD') then
    write(gnuplotchannel,*) 'set xtics ('//&
                                         '"" log10(.00000001), '//&
                                         '"" log10(.00000005), '//&
                                          '"" log10(.0000001), '//&
                                          '"" log10(.0000005), '//&
                                           '"" log10(.000001), '//&
                                           '"" log10(.000005), '//&
                                           '""  log10(.00001), '//&
                                           '""  log10(.00005), '//&
                                           '""   log10(.0001), '//&
                                           '""   log10(.0005), '//&
                                           '""    log10(.001), '//&
                                           '""    log10(.005), '//&
                                           '""     log10(.01), '//&
                                           '""     log10(.05), '//&
                                           '""      log10(.1), '//&
                                           '""      log10(.5), '//&
                                           ' ""       log10(1), '//&
                                           ' ""       log10(5), '//&
                                           ' ""      log10(10), '//&
                                           ' ""      log10(50), '//&
                                           ' ""     log10(100), '//&
                                           ' ""     log10(500), '//&
                                           ' ""    log10(1000), '//&
                                           ' ""    log10(5000), '//&
                                           ' ""   log10(10000), '//&
                                   ')'
end if

write(gnuplotchannel,FMT="(A,E16.8,A,E16.8,A)")&
        'set xrange [',lowerlimit,':',upperlimit,']'

if (comparison_number > 1) then
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '# Edit to title this plot (01)'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") 'unset arrow'
write(gnuplotchannel,FMT="(A)") 'set format x ""'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_max(1),&
        ',graph 0 to ', var_max(1), ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Amean(1),&
        ',graph 0 to ', var_Amean(1), ', graph 1 nohead front '//&
        'lw 2 lc rgb "blue"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Gmean(1),&
        ',graph 0 to ', var_Gmean(1), ', graph 1 nohead front '//&
        'lw 2 lc rgb "green"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

do n = 2, comparison_number-1

write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',n,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") 'unset arrow'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_max(n),&
        ',graph 0 to ', var_max(n), ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Amean(n),&
        ',graph 0 to ', var_Amean(n), ', graph 1 nohead front '//&
        'lw 2 lc rgb "blue"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Gmean(n),&
        ',graph 0 to ', var_Gmean(n), ', graph 1 nohead front '//&
        'lw 2 lc rgb "green"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*(n-1)+2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*(n-1)+3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*(n-1)+4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

end do
end if

write(gnuplotchannel,*) 'set xtics out nomirror font ",16"'
if (trim(adjustl(comparison_SATRVname)) /= 'InterpolationTDD') then
    write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" log10(.00000001), '//&
                                         '"5e-8" log10(.00000005), '//&
                                          '"1e-7" log10(.0000001), '//&
                                          '"5e-7" log10(.0000005), '//&
                                           '"1e-6" log10(.000001), '//&
                                           '"5e-6" log10(.000005), '//&
                                           '"1e-5"  log10(.00001), '//&
                                           '"5e-5"  log10(.00005), '//&
                                           '"1e-4"   log10(.0001), '//&
                                           '"5e-4"   log10(.0005), '//&
                                           '"1e-3"    log10(.001), '//&
                                           '"5e-3"    log10(.005), '//&
                                           '"1e-2"     log10(.01), '//&
                                           '"5e-2"     log10(.05), '//&
                                           '"1e-1"      log10(.1), '//&
                                           '"5e-1"      log10(.5), '//&
                                           ' "1e0"       log10(1), '//&
                                           ' "5e0"       log10(5), '//&
                                           ' "1e1"      log10(10), '//&
                                           ' "5e1"      log10(50), '//&
                                           ' "1e2"     log10(100), '//&
                                           ' "5e2"     log10(500), '//&
                                           ' "1e3"    log10(1000), '//&
                                           ' "5e3"    log10(5000), '//&
                                           ' "1e4"   log10(10000), '//&
                                   ')'
end if

write(gnuplotchannel,*) 'set xlabel "" font ",24"'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',comparison_number,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") 'unset arrow'
write(gnuplotchannel,FMT="(A)") 'set format x "% g"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_max(comparison_number),&
        ',graph 0 to ', var_max(comparison_number), ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Amean(comparison_number),&
        ',graph 0 to ', var_Amean(comparison_number), ', graph 1 nohead front '//&
        'lw 2 lc rgb "blue"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Gmean(comparison_number),&
        ',graph 0 to ', var_Gmean(comparison_number), ', graph 1 nohead front '//&
        'lw 2 lc rgb "green"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*comparison_number-1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*comparison_number,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*comparison_number+1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot -c "//&
        gridpath5//gnuplotfile//"_"//trim(adjustl(comparison_SATRVname))//&
        ' "'//gridpath5(1:gridpath_length+expfolder_length+5-1)//'"')

end subroutine getInterpolationplot

subroutine plotDropoff(dropoff_spacing,dropoff_spaces)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
!character(*), intent(in) :: PNGfilename

integer,intent(in) :: dropoff_spacing,dropoff_spaces
real(dp) :: max_error, min_error, max_RE, min_RE
real(dp),dimension(2*dropoff_spaces) :: dropoffErrors
integer,dimension(2*dropoff_spaces) :: dropoffErrorBins
integer :: Nbins
integer,allocatable :: dropoffErrorBinning(:,:)
real(dp),dimension(2*dropoff_spaces) :: dropoffErrorBinwidth
integer :: error_counter
integer :: iostate
integer :: i

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

max_error = 4.1198d0
min_error = 1.0712d-3

max_RE = 50.0
min_RE = 0.02

Nbins  = 30
allocate(dropoffErrorBinning(2*dropoff_spaces,Nbins))
dropoffErrorBinwidth(1:dropoff_spaces) = &
        log10(max_error/min_error)/Nbins
dropoffErrorBinwidth(dropoff_spaces+1:dropoff_spaces*2) = &
        log10(max_RE/min_RE)/Nbins

dropoffErrorBinning = 0
error_counter = 0
open(filechannel1,file=gridpath5//dropofffile)
open(filechannel2,file=gridpath5//"log"//dropofffile)
do
    read(filechannel1,iostat=iostate,FMT=*)&
        dropoffErrors
    if (iostate /= 0) exit

    write(filechannel2,FMT=*) log10(dropoffErrors)

    do i = 1, dropoff_spaces
        dropoffErrorBins(i) = floor(&
            log10(dropoffErrors(i) / min_error) / &
            dropoffErrorBinwidth(i))
    end do

    do i = dropoff_spaces+1, dropoff_spaces*2
        dropoffErrorBins(i) = floor(&
            log10(dropoffErrors(i) / min_RE) / &
            dropoffErrorBinwidth(i))
    end do

    do i = 1, dropoff_spaces * 2
        if (dropoffErrorBins(i) < 1) &
                dropoffErrorBins(i) = 1
        if (dropoffErrorBins(i) > Nbins) &
                dropoffErrorBins(i) = Nbins

        dropoffErrorBinning(i,dropoffErrorBins(i)) =&
                dropoffErrorBinning(i,dropoffErrorBins(i)) + 1
    end do

    error_counter = error_counter + 1
end do
close(filechannel1)
close(filechannel2)

open(filechannel1,file=gridpath5//temporaryfile1)
do i = 1, Nbins
    write(filechannel1,FMT=*) &
            (i-0.5) * dropoffErrorBinwidth(1) + log10(min_error), &
            dropoffErrorBinning(1:dropoff_spaces,i), &
            (i-0.5) * dropoffErrorBinwidth(dropoff_spaces+1) + log10(min_RE), &
            dropoffErrorBinning(dropoff_spaces+1:dropoff_spaces*2,i)
end do
close(filechannel1)
deallocate(dropoffErrorBinning)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo size 1200,1200'
write(gnuplotchannel,FMT="(A)") 'set encoding utf8'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//&
                                'ConsolidatedDropoff_IED.png"'
write(gnuplotchannel,FMT="(A)") "set title 'Consolidated Dropoff' font ',48'"
write(gnuplotchannel,FMT="(A)") "set xlabel 'Ninterpolation' font ',32'"
write(gnuplotchannel,FMT="(A)") "set ylabel 'Error Between Target and "//&
                                "Interpolated Energy Gradient' font ',32'"
!write(gnuplotchannel,FMT="(A)") "set logscale y"
write(gnuplotchannel,FMT="(A)") "unset key"
write(gnuplotchannel,FMT="(A,I5)") "spacing = ", dropoff_spacing
write(gnuplotchannel,FMT="(A,I5)") "spaces = ", dropoff_spaces
write(gnuplotchannel,FMT=*) "min_y = ", min_error
write(gnuplotchannel,FMT=*) "max_y = ", max_error
write(gnuplotchannel,FMT="(A)") "set xrange [0:spaces+0.5]"
write(gnuplotchannel,FMT="(A)") "set yrange [log10(min_y):log10(max_y)]"

write(gnuplotchannel,FMT="(A)") "set style fill solid 0.25 border -1"
write(gnuplotchannel,FMT="(A)") "set style boxplot outliers pointtype 7"
write(gnuplotchannel,FMT="(A)") "set style data boxplot"
write(gnuplotchannel,FMT="(A)") "unset xtics"
do i = 1, dropoff_spaces
    write(gnuplotchannel,FMT="(A,I2,A,I2,A)")&
            "set xtics add ('", i*dropoff_spacing, "' ", i, ")"
end do
write(gnuplotchannel,FMT="(A)") "plot for [i=1:spaces] '"//&
                                gridpath5//"log"//dropofffile//&
                                "' using (i):i"
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//&
        gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo size 1200,1200'
write(gnuplotchannel,FMT="(A)") 'set encoding utf8'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//&
                                'ConsolidatedDropoff_RED.png"'
write(gnuplotchannel,FMT="(A)") "set title 'Consolidated Dropoff' font ',48'"
write(gnuplotchannel,FMT="(A)") "set xlabel 'Ninterpolation' font ',32'"
write(gnuplotchannel,FMT="(A)") "set ylabel 'Relative Error Between Target and "//&
                                "Interpolated Energy Gradient' font ',32'"
!write(gnuplotchannel,FMT="(A)") "set logscale y"
write(gnuplotchannel,FMT="(A)") "unset key"
write(gnuplotchannel,FMT="(A,I5)") "spacing = ", dropoff_spacing
write(gnuplotchannel,FMT="(A,I5)") "spaces = ", dropoff_spaces
write(gnuplotchannel,FMT=*) "min_y = ", min_RE
write(gnuplotchannel,FMT=*) "max_y = ", max_RE
write(gnuplotchannel,FMT="(A)") "set xrange [0:spaces+0.5]"
write(gnuplotchannel,FMT="(A)") "set yrange [log10(min_y):log10(max_y)]"

write(gnuplotchannel,FMT="(A)") "set style fill solid 0.25 border -1"
write(gnuplotchannel,FMT="(A)") "set style boxplot outliers pointtype 7"
write(gnuplotchannel,FMT="(A)") "set style data boxplot"
write(gnuplotchannel,FMT="(A)") "unset xtics"
do i = 1, dropoff_spaces
    write(gnuplotchannel,FMT="(A,I2,A,I2,A)")&
            "set xtics add ('", i*dropoff_spacing, "' ", i, ")"
end do
write(gnuplotchannel,FMT="(A)") "plot for [i=spaces+1:2*spaces] '"//&
                                gridpath5//"log"//dropofffile//&
                                "' using (i-spaces):i"
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//&
        gridpath5//gnuplotfile)

!Now, for histograms

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo size 1200,2400'
write(gnuplotchannel,FMT="(A)") 'set encoding utf8'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//&
                                'ConsolidatedDropoff_IED_Hist.png"'
write(gnuplotchannel,FMT="(A,I5)") "spacing = ", dropoff_spacing
write(gnuplotchannel,FMT="(A,I5)") "spaces = ", dropoff_spaces
write(gnuplotchannel,FMT="(A)") "set multiplot layout spaces,1 columnsfirst "//&
                        "margins 0.1,0.95,.1,.9 spacing 0.1,0 "//&
                        "title 'Consolidated Dropoff' font ',48'"
write(gnuplotchannel,FMT="(A)") 'set style fill solid 1.0 noborder'
write(gnuplotchannel,FMT="(A)") "set ylabel 'Occurence' font ',32'"
write(gnuplotchannel,FMT="(A)") "unset xlabel"
write(gnuplotchannel,FMT="(A)") "unset xtics"
write(gnuplotchannel,FMT="(A)") "unset key"
write(gnuplotchannel,FMT=*) "min_x = ", min_error
write(gnuplotchannel,FMT=*) "max_x = ", max_error
write(gnuplotchannel,FMT="(A)") "set xrange [log10(min_x):log10(max_x)]"
write(gnuplotchannel,FMT="(A)") "set yrange [0:]"

do i = 1, dropoff_spaces
    if (i == dropoff_spaces) then
        write(gnuplotchannel,FMT="(A)") "set xtics"
        write(gnuplotchannel,FMT="(A)") &
            "set xlabel 'Error Between Target and "//&
            "Interpolated Energy Gradient' font ',32'"
    end if
    write(gnuplotchannel,FMT="(A,I0.2,A,I0.2,A)") &
            "plot '"//gridpath5//temporaryfile1//&
            "' u ($", 1, "):($", i+1, ") w boxes"
end do
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//&
        gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,FMT="(A)") 'set term pngcairo size 1200,2400'
write(gnuplotchannel,FMT="(A)") 'set encoding utf8'
write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//&
                                'ConsolidatedDropoff_RED_Hist.png"'
write(gnuplotchannel,FMT="(A,I5)") "spacing = ", dropoff_spacing
write(gnuplotchannel,FMT="(A,I5)") "spaces = ", dropoff_spaces
write(gnuplotchannel,FMT="(A)") "set multiplot layout spaces,1 columnsfirst "//&
                        "margins 0.1,0.95,.1,.9 spacing 0.1,0 "//&
                        "title 'Consolidated Dropoff' font ',48'"
write(gnuplotchannel,FMT="(A)") 'set style fill solid 1.0 noborder'
write(gnuplotchannel,FMT="(A)") "set ylabel 'Occurence' font ',32'"
write(gnuplotchannel,FMT="(A)") "unset xlabel"
write(gnuplotchannel,FMT="(A)") "unset xtics"
write(gnuplotchannel,FMT="(A)") "unset key"
write(gnuplotchannel,FMT=*) "min_x = ", min_RE
write(gnuplotchannel,FMT=*) "max_x = ", max_RE
write(gnuplotchannel,FMT="(A)") "set xrange [log10(min_x):log10(max_x)]"
write(gnuplotchannel,FMT="(A)") "set yrange [0:]"

do i = 1, dropoff_spaces
    if (i == dropoff_spaces) then
        write(gnuplotchannel,FMT="(A)") "set xtics"
        write(gnuplotchannel,FMT="(A)") &
            "set xlabel 'Relative Error Between Target and "//&
            "Interpolated Energy Gradient' font ',32'"
    end if
    write(gnuplotchannel,FMT="(A,I2,A,I2,A)") &
            "plot '"//gridpath5//temporaryfile1//&
            "' u ($", dropoff_spaces+2, "):($", dropoff_spaces+i+2, ") w boxes"
end do
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//&
        gridpath5//gnuplotfile)

return

end subroutine plotDropoff

subroutine processCheckstateFile()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

!VARIABLES IN THE INTERPOLATION FILE
real(dp),dimension(Nvar) :: vals
real(dp) :: min_rmsd, min_rmsd_prime
real(dp) :: U, KE
integer :: number_of_frames, order, neighbor_check

! For VENUS, comment out the next two lines:
!integer :: Naccept = 0
!integer :: Nalpha_tries = 1
logical :: isopened
integer :: steps_previous

integer,allocatable :: success_binning(:,:)
integer,allocatable :: failure_binning(:,:)
integer :: Nsaved,Nwasted
character(37) :: firstliner

!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

allocate(success_binning(Naccept_max,Nalpha_tries_max),&
         failure_binning(Naccept_max,Nalpha_tries_max))
success_binning = 0
failure_binning = 0

steps = 0

call system("rm "//gridpath5//"truncated"//checkstatefile)
call system("rm "//gridpath5//temporaryfile1)

open(filechannel1,file=gridpath5//checkstatefile)
open(filechannel2,file=gridpath5//"truncated"//checkstatefile)
do
    steps_previous = steps
    read(filechannel1,iostat=iostate,FMT=*)&
            number_of_frames,order,&
            neighbor_check,steps,&
            min_rmsd,min_rmsd_prime,&
            vals(1),vals(2),U,KE
    if (iostate /= 0) exit

    if (steps < steps_previous) then

        failure_binning(Naccept,Nalpha_tries) =&
            failure_binning(Naccept,Nalpha_tries) + 1

        Nalpha_tries = Nalpha_tries + 1
        Naccept = 0

        call system("rm "//gridpath5//temporaryfile1)

        if (Nalpha_tries > Nalpha_tries_max) then
            Nalpha_tries = 1

            write(filechannel2,FMT=*) number_of_frames,order,&
                    neighbor_check,steps,&
                    default_SIs(1),default_SIs(1),&
                    vals(1),vals(2),U,KE
            cycle
        end if
    end if

    if (min_rmsd_prime < outer_threshold_SI) then
        Naccept = Naccept + 1

        if (Naccept > Naccept_max) then

            success_binning(Naccept_max,Nalpha_tries) =&
                success_binning(Naccept_max,Nalpha_tries) + 1

            Nalpha_tries = 1
            Naccept = 1

            close(filechannel2)
            call system("cat "//gridpath5//temporaryfile1//&
                    " >> "//gridpath5//"truncated"//checkstatefile)
            open(filechannel2,&
                    file=gridpath5//"truncated"//checkstatefile,&
                    position="append")

            open(filechannel3,file=gridpath5//temporaryfile1)
            write(filechannel3,FMT=*) number_of_frames,order,&
                    neighbor_check,steps,min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
            close(filechannel3)
        else
            open(filechannel3,file=gridpath5//temporaryfile1,&
                 position="append")
            write(filechannel3,FMT=*) number_of_frames,order,&
                    neighbor_check,steps,min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
            close(filechannel3)
        end if
    else
        if (Naccept > 0) then

            success_binning(Naccept,Nalpha_tries) =&
                success_binning(Naccept,Nalpha_tries) + 1
    
            Nalpha_tries = 0
            Naccept = 0
    
            close(filechannel2)
            call system("cat "//gridpath5//temporaryfile1//&
                    " >> "//gridpath5//"truncated"//checkstatefile)
            open(filechannel2,&
                    file=gridpath5//"truncated"//checkstatefile,&
                    position="append")
        end if

        write(filechannel2,FMT=*) number_of_frames,order,&
                neighbor_check,steps,&
                default_SIs(1),default_SIs(1),&
                vals(1),vals(2),U,KE
    end if
end do
close(filechannel1)
close(filechannel2)

Nsaved = 0
Nwasted = 0
do Naccept = 1, Naccept_max
    Nwasted = Nwasted + Naccept * &
        sum(failure_binning(Naccept,:))
    Nsaved = Nsaved + Naccept * &
        sum(success_binning(Naccept,:))
end do

write(firstliner,FMT="(A1,4(I0.8,1x))") &
        "#",Nsaved,Nwasted,steps_previous,&
        sum(failure_binning(:,Nalpha_tries_max))

call system("sed -i '1i\"//firstliner//"' '"//&
        gridpath5//"truncated"//&
        checkstatefile//"'")

deallocate(success_binning,failure_binning)

end subroutine processCheckstateFile








subroutine processEnergyDrift()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

real(dp) :: DH, H, minDH, maxDH
integer :: metaID

integer :: totalCount

real(dp) :: DHbinwidth
integer :: Nbins, DHbin

integer,allocatable :: binnedDH(:,:)

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

Nbins = 50

totalCount = 0

mindH = 1.0d9
maxDH = -1.0d9

open(filechannel1,file=gridpath5//"energydrift.dat")
do
    read(filechannel1,iostat=iostate,FMT=*)&
        DH, H, metaID
    if (iostate /= 0) exit

    minDH = min(DH, minDH)
    maxDH = max(DH, maxDH)

    totalCount = totalCount + 1

end do
close(filechannel1)

allocate(binnedDH(3,Nbins))
binnedDH = 0
DHbinwidth = (maxDH - minDH) / Nbins

open(filechannel1,file=gridpath5//"energydrift.dat")
open(filechannel2,file=gridpath5//"binnedenergydrift.dat")
do
    read(filechannel1,iostat=iostate,FMT=*)&
        DH, H, metaID
    if (iostate /= 0) exit

    DHbin = floor((DH-minDH)/DHbinwidth) + 1
    if (DHbin < 0) DHbin = 1
    if (DHbin > Nbins) DHbin = Nbins

    binnedDH(metaID+1,DHbin) = &
        binnedDH(metaID+1,DHbin) + 1

end do
close(filechannel1)
close(filechannel2)





open(6666,file=gridpath5//"tmpDH.dat")
do DHbin = 1, Nbins
    write(6666,FMT=*) &
        (DHbin-1)*DHbinwidth + minDH, &
        binnedDH(:,DHbin)
end do
close(6666)
deallocate(binnedDH)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "binnedenergydrift.png"'
write(gnuplotchannel,*) 'set title "Energy Drifts Per Frame"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'minx = ', minDH
write(gnuplotchannel,*) 'maxx = ', maxDH
write(gnuplotchannel,*) 'set yrange [0:]'

if (.false.) then
    write(gnuplotchannel,*) 'set logscale x'
    write(gnuplotchannel,FMT="(A)") &
            'set xrange [log10(minx):log10(maxx)]'
else
    write(gnuplotchannel,FMT="(A)") &
            'set xrange [minx:maxx]'
end if


write(gnuplotchannel,*) 'set xlabel "Energy Drift (kcal/mol)"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'plot "'//gridpath5//'tmpDH.dat" u 1:2 w boxes,\'
write(gnuplotchannel,*) '     "'//gridpath5//'tmpDH.dat" u 1:3 w boxes,\'
write(gnuplotchannel,*) '     "'//gridpath5//'tmpDH.dat" u 1:4 w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)




end subroutine processEnergyDrift






subroutine readNextQP(Q,P)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

double precision :: T, V, H

double precision,dimension(3*Natoms),intent(out) :: Q, P
character(200) :: readtrajectory_aline
integer :: i,j

do
    read(filechannel6,iostat=rxnStage,&
         FMT="(A)") readtrajectory_aline
    if (rxnStage /= 0) exit
    if (readtrajectory_aline(1:5) == " XXXX") then
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline
        if (readtrajectory_aline(3:7) /= "THE C") exit
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline
        if (rxnStage /= 0) exit

        read(filechannel6,FMT=&
        "(18x,ES17.9,22x,ES17.9,A)")&
        T, V, readtrajectory_aline
        read(filechannel6,FMT=&
        "(18x,ES17.9,A)")&
        H, readtrajectory_aline
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline

        do j = 1, Natoms
            read(filechannel6,FMT=&
            "(3(F11.7,1x),3x,3(F11.7,1x),A)")&
            Q(j*3-2:j*3), P(j*3-2:j*3)
        end do
    end if
end do

return
end subroutine readNextQP


subroutine readNextQPandPDOT(Q,P,PDOT)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

double precision :: T, V, H

double precision,dimension(3*Natoms),intent(out) :: Q, P, PDOT
double precision,dimension(3*Natoms) :: Q2

real(dp),parameter :: bohr2angstrom = 0.529177d0
character(200) :: readtrajectory_aline
integer :: i,j

logical :: QP_have_been_read_in
integer :: iostate
 
QP_have_been_read_in = .false.
do
    read(filechannel6,iostat=rxnStage,&
         FMT="(A)") readtrajectory_aline
    if (rxnStage /= 0) exit
    if (readtrajectory_aline(59:66) == "gradient") then
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline
        if (rxnStage /= 0) exit

        do j = 1, Natoms
!           read(filechannel6,FMT=&
!           "(44x,3(1x,F10.6))") PDOT(j*3-2:j*3)
            read(filechannel6,iostat=iostate,FMT=&
            "(10x,3(1x,F10.6),1x,3(1x,F10.6))")&
            Q2(j*3-2:j*3), PDOT(j*3-2:j*3)
            if (iostate /= 0) cycle
        end do

    elseif (readtrajectory_aline(1:5) == " XXXX") then
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline
        if (readtrajectory_aline(3:7) /= "THE C") exit
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline
        if (rxnStage /= 0) exit

        read(filechannel6,FMT=&
        "(18x,ES17.9,22x,ES17.9,A)")&
        T, V, readtrajectory_aline
        read(filechannel6,FMT=&
        "(18x,ES17.9,A)")&
        H, readtrajectory_aline
        read(filechannel6,iostat=rxnStage,&
             FMT="(A)") readtrajectory_aline

        do j = 1, Natoms
            read(filechannel6,FMT=&
            "(3(F11.7,1x),3x,3(F11.7,1x),A)")&
            Q(j*3-2:j*3), P(j*3-2:j*3)
        end do

        QP_have_been_read_in = .true.

    end if
end do

return
end subroutine readNextQPandPDOT

subroutine readNextQandPDOTfromXYZ(Q,PDOT)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

double precision,dimension(3*Natoms),intent(out) :: Q, PDOT
character(200) :: readtrajectory_aline
integer :: i,j

integer :: iostate

read(filechannel6,FMT=*,iostat=iostate) i
if (iostate /= 0) then
  close(filechannel6)
  stop
end if

read(filechannel6,FMT=*)
do j = 1, Natoms
  read(filechannel6,FMT="(2x,A)") readtrajectory_aline
  read(readtrajectory_aline,FMT=*) &
            Q(j*3-2:j*3), PDOT(j*3-2:j*3)
  print *, j, readtrajectory_aline
end do

return
end subroutine readNextQandPDOTfromXYZ


subroutine getIMRRerrors1(Q,PDOT)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use VARIABLES
use SIMILARITY
use interactMultipleGrids
implicit none

COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC
COMMON/FORCES/NATOMS

integer :: NATOMS, NC
double precision :: NTZ,NT
double precision :: ISEED0
double precision :: TIME
double precision :: T,V,H
integer :: Ncoords

double precision,dimension(3*Natoms),intent(in) :: Q, PDOT

integer :: Ninterpolation_pool
integer, allocatable :: reorderedIndexes(:)

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)
real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp), allocatable :: inputCLS2(:,:)
real(dp), allocatable :: libUgradients(:,:,:)

real(dp),dimension(3,3) :: U
real(dp),dimension(NSIs) :: new_SIs

real(dp),dimension(3,Natoms) :: coords1,coords2
real(dp),dimension(3,Natoms) :: gradient1,gradient2
real(dp) :: min_rmsd, errorIMRR
real(dp) :: minErrorIMRR, maxErrorIMRR
real(dp) :: errorIMRRbinwidth
integer :: maxErrorIMRRoccurence

character(200) :: FMTalpha
integer :: Nbins, errorBin
integer,allocatable :: errorIMRRbinning(:,:,:)

real(dp),dimension(3*NATOMS) :: meanCoords
real(dp) :: varCoords

integer*8 :: dummyINT1, dummyINT2, dummyINT3

integer, save :: getIMRRerrors_counter = 0
integer :: iostate, step, i, j, k

Ncoords = 3 * NATOMS

coords1 = reshape(Q,(/3,Natoms/))
gradient1 = reshape(PDOT,(/3,Natoms/))

call getVarsHBrCO2(coords1,Natoms,vals,&
               Nvar,BOND_LABELLING_DATA)

!Nsort = 1
!print *, "whaaat 1"
call checkState_PCM(&
        vals,coords1,&
        gradient2,&
        dummyINT1,dummyINT2,dummyINT3)
!print *, "whaaat 2"

! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (getIMRRerrors_counter > 1) return
!if (any(vals > 6.0d0)) return
!if (Ninterpolation < 25) return
if (Ninterpolation < Ninterpolation_cap+5) return
! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ninterpolation_pool = Ninterpolation
Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)
!       20)

allocate(reorderedIndexes(Ninterpolation),&
         frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         inputCLS2(Ncoords+Ninterpolation,Ninterpolation),&
         libUgradients(3,NATOMS,Ninterpolation))

call setTarget(coords1)

minErrorIMRR = 1.0d9
maxErrorIMRR = 0.0d0

open(6666,file=gridpath5//"errorIMRRtest.dat")

do j = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1))
    end if

    do i = 1, Ninterpolation
        reorderedIndexes(i) = i
    end do

do k = 1, (Ninterpolation_pool - Ninterpolation)*2 + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part only looks at the inputs

do i = 1, Ninterpolation
    step = reorderedIndexes(i)

    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    libUgradients(:,:,i) = matmul(&
            U,gradientbuffer1(:,:,step))

    inputCLS2(1:Ncoords,i) = reshape(&
            coords2 - coords1,(/Ncoords/))
    inputCLS2(Ncoords+i,:) = 0.0d0
    inputCLS2(Ncoords+i,i) = alpha_ratio * &
            new_SIs(Nsort)**2
end do
!errorBaseline = sqrt(sum((gradient1 - &
!        libUgradients(:,:,1))**2)/Natoms)

meanCoords = 0.0d0
do i = 1, Ninterpolation
    meanCoords = meanCoords + &
                 inputCLS2(1:Ncoords,i)
end do
meanCoords = meanCoords / Ninterpolation

varCoords = 0.0d0
do i = 1, Ninterpolation
    varCoords = varCoords + &
                sum((inputCLS2(1:Ncoords,i)-meanCoords)**2)
end do
varCoords = sqrt(varCoords/Ninterpolation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outputCLS = 0.0d0
restraints = 1.0d0
restraint_values = 1.0d0

call CLS2(inputCLS2(1:Ncoords+Ninterpolation,&
                    1:Ninterpolation),&
          Ncoords+Ninterpolation,&
          Ninterpolation,&
          restraints,1,restraint_values,&
          outputCLS,frame_weights)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part looks at the interpolation

if (all(frame_weights == 0.0d0)) then
    do i = 1, Ninterpolation
        outputCLS(i) = inputCLS2(Ncoords+i,i)
    end do
    frame_weights(minloc(outputCLS(&
                  1:Ninterpolation))) = 1.0d0
end if

gradient2 = 0.0d0
do i = 1, Ninterpolation
    gradient2 = gradient2 + frame_weights(i)*&
                libUgradients(:,:,i)
end do

errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(6666,FMT="(I2,1x,I3,3(1x,ES12.7))") j,k,&
write(6666,FMT="(I2,1x,I3,3(1x,ES14.7))") j,k,&
                  sqrt(sum(meanCoords**2)/Natoms), &
                  varCoords, errorIMRR

minErrorIMRR = min(minErrorIMRR,errorIMRR)
maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (k <= (Ninterpolation_pool - Ninterpolation)) then
    reorderedIndexes(&
         maxloc(abs(frame_weights(1:Ninterpolation)))&
                    ) = Ninterpolation + k
else
    if (k == (Ninterpolation_pool - Ninterpolation) + 1) then
    do i = 1, Ninterpolation
        reorderedIndexes(i) = i
    end do
    end if

    reorderedIndexes(&
         minloc(abs(frame_weights(1:Ninterpolation)))&
                    ) = Ninterpolation + k - &
                        (Ninterpolation_pool - Ninterpolation)
end if

end do

end do

close(6666)







Nbins = 20
errorIMRRbinwidth = (maxErrorIMRR - minErrorIMRR) / Nbins
allocate(errorIMRRbinning(Nbins,3,Nalpha))

errorIMRRbinning = 0
open(6666,file=gridpath5//"errorIMRRtest.dat")
do
    read(6666,iostat=iostate,FMT=*) j,k,&
              min_rmsd, varCoords, errorIMRR
    if (iostate /= 0) exit

    errorBin = floor((errorIMRR-minErrorIMRR)/errorIMRRbinwidth) + 1
    if (errorBin < 1) errorBin = 1
    if (errorBin > Nbins) errorBin = Nbins

    if (k == 1) then
        errorIMRRbinning(errorBin,1,j) = &
            errorIMRRbinning(errorBin,1,j) + 1
    end if
    if (k <= (Ninterpolation_pool-Ninterpolation)+1) then
        errorIMRRbinning(errorBin,2,j) = &
            errorIMRRbinning(errorBin,2,j) + 1
    end if
    if (.true.) then
        errorIMRRbinning(errorBin,3,j) = &
            errorIMRRbinning(errorBin,3,j) + 1
    end if

end do
close(6666)

maxErrorIMRRoccurence = maxval(&
        errorIMRRbinning(1:Nbins,1:3,1:Nalpha))

FMTalpha = "(ES14.7"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,I6"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"errorIMRRtest_binned1.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,1,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned2.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,2,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned3.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,3,1:Nalpha)
end do
close(6666)

deallocate(errorIMRRbinning)

!getIMRRerrors_counter = getIMRRerrors_counter + 1
getIMRRerrors_counter = NC

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRR',&
        getIMRRerrors_counter,1,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "N = ', Ninterpolation,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
!write(6666,*) 'set tmargin at screen 0.8'
!write(6666,*) 'set bmargin at screen 0.2'
!write(6666,*) 'set rmargin at screen 0.8'
!write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'Nalpha = ', Nalpha
write(6666,FMT=*) 'ymax = ', maxErrorIMRRoccurence
write(6666,FMT=*) 'xmin = ', max(minErrorIMRR-errorIMRRbinwidth,0.0d0)
write(6666,FMT=*) 'xmax = ', maxErrorIMRR+errorIMRRbinwidth

write(6666,FMT="(A)") 'set multiplot layout Nalpha,1 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0.1,0'
write(6666,FMT="(A)") 'set xrange [xmin:xmax]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set style fill transparent solid 1.0'
write(6666,FMT=*) 'set boxwidth ', errorIMRRbinwidth

write(6666,FMT="(A)") 'unset xlabel'
write(6666,FMT="(A)") 'unset xtics'
do i = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1))
    end if
    write(6666,FMT="(A,F7.3,A)") &
            'set ylabel "P_{1.0e',log10(alpha_ratio),'}"'
    if (i == Nalpha) then
        write(6666,FMT="(A)") 'set xtics'
        write(6666,FMT="(A)") 'set xlabel "IMRR Error"'
    end if
    write(6666,FMT="(A,I1,A)") &
            'plot "'//gridpath5//'errorIMRRtest_binned3.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "red",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned2.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "green",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned1.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "blue"'
end do
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

deallocate(reorderedIndexes,&
           frame_weights,&
           outputCLS,restraints,&
           restraint_values,&
           inputCLS2,libUgradients)

return
end subroutine getIMRRerrors1


subroutine getIMRRerrors2(Q,PDOT)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use VARIABLES
use SIMILARITY
use interactMultipleGrids
implicit none

COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC
COMMON/FORCES/NATOMS

integer :: NATOMS, NC
double precision :: NTZ,NT
double precision :: ISEED0
double precision :: TIME
double precision :: T,V,H
integer :: Ncoords

double precision,dimension(3*Natoms),intent(in) :: Q, PDOT

integer :: Ninterpolation_pool
integer, allocatable :: reorderedIndexes(:)

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)
real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp), allocatable :: inputCLS2(:,:)
real(dp), allocatable :: libUgradients(:,:,:)

real(dp),dimension(3,3) :: U
real(dp),dimension(NSIs) :: new_SIs

real(dp),dimension(3,Natoms) :: coords1,coords2
real(dp),dimension(3,Natoms) :: gradient1,gradient2
real(dp) :: min_rmsd, errorIMRR
real(dp) :: minErrorIMRR, maxErrorIMRR
real(dp) :: errorIMRRbinwidth
integer :: maxErrorIMRRoccurence

character(200) :: FMTalpha
integer :: Nbins, errorBin
integer,allocatable :: errorIMRRbinning(:,:,:)

real(dp),allocatable :: CMvec1(:),CMvec2(:),meanCoords(:)
real(dp) :: varCoords

integer*8 :: dummyINT1, dummyINT2, dummyINT3

integer, save :: getIMRRerrors_counter = 0
integer :: iostate, step, i, j, k

Ncoords = (NATOMS*(NATOMS-1))/2

coords1 = reshape(Q,(/3,Natoms/))
gradient1 = reshape(PDOT,(/3,Natoms/))

call getVarsHBrCO2(coords1,Natoms,vals,&
               Nvar,BOND_LABELLING_DATA)

! We will test the second geometric
! representation by (1) sorting by it
! and (2) using it for the IMRR
!Nsort = 2
call checkState_PCM(&
        vals,coords1,&
        gradient2,&
        dummyINT1,dummyINT2,dummyINT3)

! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (getIMRRerrors_counter > 1) return
!if (any(vals > 6.0d0)) return
!if (Ninterpolation < 25) return
if (Ninterpolation < Ninterpolation_cap+5) return
! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ninterpolation_pool = Ninterpolation
Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)
!       20)

allocate(reorderedIndexes(Ninterpolation),&
         frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         inputCLS2(Ncoords+Ninterpolation,Ninterpolation),&
         libUgradients(3,NATOMS,Ninterpolation),&
         CMvec1(Ncoords),CMvec2(Ncoords),meanCoords(Ncoords))

call setTarget(coords1)
call getCoulombMatrixVector(&
        coords1,CMvec1)

minErrorIMRR = 1.0d9
maxErrorIMRR = 0.0d0

open(6666,file=gridpath5//"errorIMRRtest.dat")

do j = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1))
    end if

    do i = 1, Ninterpolation
        reorderedIndexes(i) = i
    end do

do k = 1, (Ninterpolation_pool - Ninterpolation)*2 + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part only looks at the inputs

do i = 1, Ninterpolation
    step = reorderedIndexes(i)

    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    libUgradients(:,:,i) = matmul(&
            U,gradientbuffer1(:,:,step))

    call getCoulombMatrixVector(&
        coords2,CMvec2)

    inputCLS2(1:Ncoords,i) = &
            CMvec2 - CMvec1
    inputCLS2(Ncoords+i,:) = 0.0d0
    inputCLS2(Ncoords+i,i) = alpha_ratio * &
            new_SIs(Nsort)**2
end do
!errorBaseline = sqrt(sum((gradient1 - &
!        libUgradients(:,:,1))**2)/Natoms)

meanCoords = 0.0d0
do i = 1, Ninterpolation
    meanCoords = meanCoords + &
                 inputCLS2(1:Ncoords,i)
end do
meanCoords = meanCoords / Ninterpolation

varCoords = 0.0d0
do i = 1, Ninterpolation
    varCoords = varCoords + &
                sum((inputCLS2(1:Ncoords,i)-meanCoords)**2)
end do
varCoords = sqrt(varCoords/Ninterpolation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outputCLS = 0.0d0
restraints = 1.0d0
restraint_values = 1.0d0

call CLS2(inputCLS2(1:Ncoords+Ninterpolation,&
                    1:Ninterpolation),&
          Ncoords+Ninterpolation,&
          Ninterpolation,&
          restraints,1,restraint_values,&
          outputCLS,frame_weights)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part looks at the interpolation

if (all(frame_weights == 0.0d0)) then
    do i = 1, Ninterpolation
        outputCLS(i) = inputCLS2(Ncoords+i,i)
    end do
    frame_weights(minloc(outputCLS(&
                  1:Ninterpolation))) = 1.0d0
end if

gradient2 = 0.0d0
do i = 1, Ninterpolation
    gradient2 = gradient2 + frame_weights(i)*&
                libUgradients(:,:,i)
end do

errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(6666,FMT="(I2,1x,I3,3(1x,ES12.7))") j,k,&
write(6666,FMT="(I2,1x,I3,3(1x,ES14.7))") j,k,&
                  sqrt(sum(meanCoords**2)/Natoms), &
                  varCoords, errorIMRR

minErrorIMRR = min(minErrorIMRR,errorIMRR)
maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (k <= (Ninterpolation_pool - Ninterpolation)) then
    reorderedIndexes(&
         maxloc(abs(frame_weights(1:Ninterpolation)))&
                    ) = Ninterpolation + k
else
    if (k == (Ninterpolation_pool - Ninterpolation) + 1) then
    do i = 1, Ninterpolation
        reorderedIndexes(i) = i
    end do
    end if

    reorderedIndexes(&
         minloc(abs(frame_weights(1:Ninterpolation)))&
                    ) = Ninterpolation + k - &
                        (Ninterpolation_pool - Ninterpolation)
end if

end do

end do

close(6666)







Nbins = 20
errorIMRRbinwidth = (maxErrorIMRR - minErrorIMRR) / Nbins
allocate(errorIMRRbinning(Nbins,3,Nalpha))

errorIMRRbinning = 0
open(6666,file=gridpath5//"errorIMRRtest.dat")
do
    read(6666,iostat=iostate,FMT=*) j,k,&
              min_rmsd, varCoords, errorIMRR
    if (iostate /= 0) exit

    errorBin = floor((errorIMRR-minErrorIMRR)/errorIMRRbinwidth) + 1
    if (errorBin < 1) errorBin = 1
    if (errorBin > Nbins) errorBin = Nbins

    if (k == 1) then
        errorIMRRbinning(errorBin,1,j) = &
            errorIMRRbinning(errorBin,1,j) + 1
    end if
    if (k <= (Ninterpolation_pool-Ninterpolation)+1) then
        errorIMRRbinning(errorBin,2,j) = &
            errorIMRRbinning(errorBin,2,j) + 1
    end if
    if (.true.) then
        errorIMRRbinning(errorBin,3,j) = &
            errorIMRRbinning(errorBin,3,j) + 1
    end if

end do
close(6666)

maxErrorIMRRoccurence = maxval(&
        errorIMRRbinning(1:Nbins,1:3,1:Nalpha))

FMTalpha = "(ES14.7"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,I6"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"errorIMRRtest_binned1.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,1,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned2.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,2,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned3.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,3,1:Nalpha)
end do
close(6666)

deallocate(errorIMRRbinning)

!getIMRRerrors_counter = getIMRRerrors_counter + 1
getIMRRerrors_counter = NC

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRR',&
        getIMRRerrors_counter,2,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "N = ', Ninterpolation,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
!write(6666,*) 'set tmargin at screen 0.8'
!write(6666,*) 'set bmargin at screen 0.2'
!write(6666,*) 'set rmargin at screen 0.8'
!write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'Nalpha = ', Nalpha
write(6666,FMT=*) 'ymax = ', maxErrorIMRRoccurence
write(6666,FMT=*) 'xmin = ', max(minErrorIMRR-errorIMRRbinwidth,0.0d0)
write(6666,FMT=*) 'xmax = ', maxErrorIMRR+errorIMRRbinwidth

write(6666,FMT="(A)") 'set multiplot layout Nalpha,1 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0.1,0'
write(6666,FMT="(A)") 'set xrange [xmin:xmax]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set style fill transparent solid 1.0'
write(6666,FMT=*) 'set boxwidth ', errorIMRRbinwidth

write(6666,FMT="(A)") 'unset xlabel'
write(6666,FMT="(A)") 'unset xtics'
do i = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1))
    end if
    write(6666,FMT="(A,F7.3,A)") &
            'set ylabel "P_{1.0e',log10(alpha_ratio),'}"'
    if (i == Nalpha) then
        write(6666,FMT="(A)") 'set xtics'
        write(6666,FMT="(A)") 'set xlabel "IMRR Error"'
    end if
    write(6666,FMT="(A,I1,A)") &
            'plot "'//gridpath5//'errorIMRRtest_binned3.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "red",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned2.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "green",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned1.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "blue"'
end do
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

deallocate(reorderedIndexes,&
           frame_weights,&
           outputCLS,restraints,&
           restraint_values,&
           inputCLS2,libUgradients,&
           CMvec1,CMvec2,meanCoords)

return
end subroutine getIMRRerrors2

subroutine getIMRRerrors3(Q,PDOT,labelID)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use VARIABLES
use SIMILARITY
use interactMultipleGrids
implicit none

COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC
COMMON/FORCES/NATOMS

integer :: NATOMS, NC
double precision :: NTZ,NT
double precision :: ISEED0
double precision :: TIME
double precision :: T,V,H

integer :: Ncoords

double precision,dimension(3*Natoms),intent(in) :: Q, PDOT
integer,intent(in) :: labelID

integer :: Ninterpolation_pool
integer, allocatable :: reorderedIndexes(:)

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)
real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp), allocatable :: inputCLS2(:,:)
real(dp), allocatable :: libUgradients(:,:,:)
real(dp), allocatable :: libenergies(:)

real(dp),dimension(3,3) :: U
real(dp),dimension(NSIs) :: new_SIs

real(dp),dimension(3,Natoms) :: coords1,coords2
real(dp),dimension(3,Natoms) :: gradient1,gradient2
real(dp) :: min_rmsd, errorIMRR, Efinal
real(dp) :: minErrorIMRR, maxErrorIMRR
real(dp) :: errorIMRRbinwidth
integer :: maxErrorIMRRoccurence

character(400) :: FMTalpha, filetext
integer :: Nbins, errorBin
integer,allocatable :: errorIMRRbinning(:,:,:)

real(dp),dimension(3*NATOMS) :: meanCoords
real(dp) :: varCoords
real(dp) :: igDistance, igDistanceMax, varCoordsMax
integer :: maxNinterpolation_current

real(dp) :: R1, R2, maxR1, maxR2, meanR1, meanR2
real(dp) :: maxErrorIMRR_cb = 0.06d0

real(dp) :: weightThreshold = 5.0d-3
integer :: errorID, Ninterpolation_current

real(dp) :: totalWeight,largestNegativeWeight
real(dp) :: l1NegativeWeight, l2NegativeWeight
real(dp),dimension((NATOMS*(NATOMS-1))/2) :: meanCM, currentCM
real(dp) :: icmDistance, varCM
real(dp),dimension((NATOMS*(NATOMS-1))/2) :: meanDM, currentDM
real(dp) :: idmDistance, varDM

character(len=*),parameter :: consolidatedErrorIMRRfile = &
                  "consolidatedErrorIMRR.dat"
real(dp) :: minR2
real(dp),dimension(Nalpha) :: chosenErrorIMRR, &
                chosenR1IMRR, chosenR2IMRR, &
                chosenLiNW, chosenL1NW, chosenL2NW, &
                chosenIGD, chosenIGV, &
                chosenIDMD, chosenIDMV, &
                chosenICMD, chosenICMV
real(dp) :: potE2

integer*8 :: dummyINT1, dummyINT2, dummyINT3

integer, save :: getIMRRerrors_counter = 0
integer :: iostate, step, i, j, k

real :: rand
integer :: randomK

logical :: output_flag
logical :: potE_flag = .false.

output_flag = (modulo(NC,100)==1)

Ncoords = 3 * NATOMS

coords1 = reshape(Q,(/3,Natoms/))
gradient1 = reshape(PDOT,(/3,Natoms/))

call getVarsHBrCO2(coords1,Natoms,vals,&
               Nvar,BOND_LABELLING_DATA)

!print *, "whaaat 1"
call checkState_PCM(&
        vals,coords1,&
        gradient2,&
        dummyINT1,dummyINT2,dummyINT3)
!print *, "whaaat 2"

! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (getIMRRerrors_counter > 1) return
!if (any(vals > 6.0d0)) return
if (Ninterpolation < 20) return !< 25) return
!if (Ninterpolation < Ninterpolation_cap+5) return
! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ninterpolation_pool = min(Ninterpolation,Ninterpolation_cap)
Ninterpolation = min(Ninterpolation,&
       Ninterpolation_cap)

!       20)
randomK = 1+5+floor(rand()*(Ninterpolation_pool-(1+5)))
randomK = min(max(randomK,1+5),Ninterpolation_pool)

allocate(reorderedIndexes(Ninterpolation),&
         frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         inputCLS2(Ncoords+Ninterpolation,Ninterpolation),&
         libUgradients(3,NATOMS,Ninterpolation),&
         libenergies(Ninterpolation))

call setTarget(coords1)

maxNinterpolation_current = 0
maxR1 = 0.0d0
maxR2 = 0.0d0
meanR1 = 0.0d0
meanR2 = 0.0d0
minErrorIMRR = 1.0d9
maxErrorIMRR = 0.0d0
igDistanceMax = 0.0d0
varCoordsMax = 0.0d0

open(6666,file=gridpath5//"errorIMRRtest.dat")

do j = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1))
    end if

    errorID = 3
    Ninterpolation_current = Ninterpolation !_pool
    k = Ninterpolation_current+1
    do i = 1, Ninterpolation_current
        reorderedIndexes(i) = i
    end do

    chosenErrorIMRR(j) = 1.0d9
    minR2 = 1.0d9

do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part only looks at the inputs

if (.false.) then
open(6667,file="/home/kazuumi/rsun_lts/kazuumi/"//&
                "dropoffs/formatted/"//&
                "10_"//"11.1674_11.2945.dat")
!               "10_"//"10.2417_10.4238.dat")
read(6667,FMT="(15(ES14.7,1x))") coords1
read(6667,FMT="(ES14.7)") V
read(6667,FMT="(15(ES14.7,1x))") gradient1

coords1 = coords1 !* 0.529177d0
call setTarget(coords1)

print *, "target"
print *, coords1(1:3,1)
print *, coords1(1:3,2)
print *, coords1(1:3,3)
print *, coords1(1:3,4)
print *, coords1(1:3,5)

do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)

    read(6667,FMT="(15(ES14.7,1x))",iostat=iostate) &
            coordsbuffer1(:,:,step)
    if (iostate /= 0) exit
    read(6667,FMT="(ES14.7)") &
            potEbuffer1(step)
    read(6667,FMT="(15(ES14.7,1x))") &
            gradientbuffer1(:,:,step)

    coordsbuffer1(:,:,step) = coordsbuffer1(:,:,step) !* 0.529177d0

    print *, "input", i
    print *, coordsbuffer1(1:3,1,step)
    print *, coordsbuffer1(1:3,2,step)
    print *, coordsbuffer1(1:3,3,step)
    print *, coordsbuffer1(1:3,4,step)
    print *, coordsbuffer1(1:3,5,step)

    call getSIs(coordsbuffer1(:,:,i),&
            coords2,U,new_SIs)
!   libUgradients(:,:,i) = &
!           gradientbuffer1(:,:,step)
    libUgradients(:,:,i) = matmul(&
            U,gradientbuffer1(:,:,step))
    libenergies(i) = &
            potEbuffer1(step)

    inputCLS2(1:Ncoords,i) = reshape(&
            coords2 - coords1,(/Ncoords/))
    inputCLS2(Ncoords+i,:) = 0.0d0
    inputCLS2(Ncoords+i,i) = alpha_ratio * &
            new_SIs(Nsort)**2
!           new_SIs(1)**2
end do
close(6667)

else
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)

    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    libUgradients(:,:,i) = matmul(&
            U,gradientbuffer1(:,:,step))
    libenergies(i) = &
            potEbuffer1(step)

    inputCLS2(1:Ncoords,i) = reshape(&
            coords2, (/Ncoords/)) !!!!!!!!!!!!!!!!!!!!!!!
!           coords2 - coords1,(/Ncoords/))
    inputCLS2(Ncoords+i,:) = 0.0d0
    inputCLS2(Ncoords+i,i) = alpha_ratio * &
            new_SIs(Nsort)**2
end do
end if

meanCoords = 0.0d0
do i = 1, Ninterpolation_current
    meanCoords = meanCoords + &
                 inputCLS2(1:Ncoords,i)
end do
meanCoords = meanCoords / Ninterpolation_current
igDistance = sqrt(sum(meanCoords**2)/Natoms)

varCoords = 0.0d0
do i = 1, Ninterpolation_current
    varCoords = varCoords + &
                sum((inputCLS2(1:Ncoords,i)-meanCoords)**2)
end do
varCoords = sqrt(varCoords/Ninterpolation_current)

igDistanceMax = max(igDistance,igDistanceMax)
varCoordsMax = max(varCoords,varCoordsMax)

meanCM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector_test(&
            coordsbuffer1(:,:,step),currentCM)
    meanCM = meanCM + currentCM
end do
call getCoulombMatrixVector_test(coords1,currentCM)
meanCM = (meanCM / Ninterpolation_current) - currentCM
icmDistance = sqrt(sum(meanCM**2)/(Natoms*(Natoms-1)*0.5))

varCM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector_test(&
            coordsbuffer1(:,:,step),currentCM)
    varCM = varCM + sum((currentCM-meanCM)**2)
end do
varCM = sqrt(varCM/Ninterpolation_current)

meanCM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector(&
            coordsbuffer1(:,:,step),currentCM)
    meanCM = meanCM + currentCM
end do
call getCoulombMatrixVector(coords1,currentCM)
meanCM = (meanCM / Ninterpolation_current) - currentCM
idmDistance = sqrt(sum(meanCM**2)/(Natoms*(Natoms-1)*0.5))

varDM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector(&
            coordsbuffer1(:,:,step),currentCM)
    varDM = varDM + sum((currentCM-meanCM)**2)
end do
varDM = sqrt(varDM/Ninterpolation_current)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outputCLS = 0.0d0
outputCLS(1:Ncoords) = reshape(coords1,(/Ncoords/)) !!!!!!!!!!!!!!!!!
restraints = 1.0d0
restraint_values = 1.0d0

!call CLS2(inputCLS2(1:Ncoords+Ninterpolation_current,&
!                    1:Ninterpolation_current),&
!          Ncoords+Ninterpolation_current,&
!          Ninterpolation_current,&
!          restraints(1:1,1:Ninterpolation_current),&
!          1,restraint_values,&
!          outputCLS(1:Ncoords+Ninterpolation_current),&
!          frame_weights(1:Ninterpolation_current))
call CLSwK(inputCLS2(1:Ncoords,&
                    1:Ninterpolation_current),&
          Ncoords,&
          Ninterpolation_current,&
          restraints(1:1,1:Ninterpolation_current),&
          1,restraint_values,&
          outputCLS(1:Ncoords),&
          frame_weights(1:Ninterpolation_current))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part looks at the interpolation

if (all(frame_weights(1:Ninterpolation_current) == 0.0d0)) then
    if (Ninterpolation_current < Ninterpolation) &
        Ninterpolation_current = &
            Ninterpolation_current + 1

    reorderedIndexes(Ninterpolation_current) = k
    k = k + 1

    if (k > Ninterpolation_pool) exit
    cycle
end if

gradient2 = 0.0d0
potE2 = 0.0d0
R2 = 0.0d0
do i = 1, Ninterpolation_current
    gradient2 = gradient2 + frame_weights(i)*&
                libUgradients(:,:,i)
    potE2 = potE2 + frame_weights(i)*&
                libenergies(i)
    R2 = R2 + abs(frame_weights(i)) * &
            inputCLS2(Ncoords+i,i)
end do
R2 = R2 / alpha_ratio
R1 = sqrt(sum(matmul(&
    inputCLS2(1:Ncoords,1:Ninterpolation_current),&
    reshape(frame_weights(1:Ninterpolation_current),&
            (/ Ninterpolation_current, 1 /)))**2))

totalWeight = sum(frame_weights(1:Ninterpolation_current))
l1NegativeWeight = abs(sum(frame_weights(1:Ninterpolation_current),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0))
l2NegativeWeight = sqrt(sum((frame_weights(1:Ninterpolation_current)**2),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0))
largestNegativeWeight = maxval(abs(frame_weights(1:Ninterpolation_current)),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0)


if (potE_flag) then
    errorIMRR = abs(V - potE2)
else
    errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (k == Ninterpolation_pool) errorID = 1

open(6667,file=gridpath5//"all"//consolidatedErrorIMRRfile,&
          position="append")
write(6667,FMT="(I8,1x,I3,3(1x,ES14.7))") &
    NC, labelID, errorIMRR, R1, R2
close(6667)

!write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2)") &
if (output_flag) &
write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2,1x,ES14.7)") &
                  j,errorID,&
                  errorIMRR,igDistance,varCoords,&
                  R1, R2, Ninterpolation_current,&
                  largestNegativeWeight
!                 totalWeight,largestNegativeWeight

maxNinterpolation_current = max(&
             Ninterpolation_current,&
          maxNinterpolation_current)
maxR1 = max(R1,maxR1)
maxR2 = max(R2,maxR2)
meanR1 = meanR1 + R1
meanR2 = meanR2 + R2
maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!
! Selection Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(ssID)
! Smallest R2
case(1)
    if (R2 < minR2) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
        chosenIGD(j) = igDistance
        chosenIGV(j) = varCoords
        chosenIDMD(j) = idmDistance
        chosenIDMV(j) = varDM
        chosenICMD(j) = icmDistance
        chosenICMV(j) = varCM
        chosenL1NW(j) = l1NegativeWeight
        chosenL2NW(j) = l2NegativeWeight
        chosenLiNW(j) = largestNegativeWeight
    end if
! Minimal Error
case(2)
    if (errorIMRR < chosenErrorIMRR(j)) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
        chosenIGD(j) = igDistance
        chosenIGV(j) = varCoords
        chosenIDMD(j) = idmDistance
        chosenIDMV(j) = varDM
        chosenICMD(j) = icmDistance
        chosenICMV(j) = varCM
        chosenL1NW(j) = l1NegativeWeight
        chosenL2NW(j) = l2NegativeWeight
        chosenLiNW(j) = largestNegativeWeight
    end if
! Random
case(3)
    if (k == randomK) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
        chosenIGD(j) = igDistance
        chosenIGV(j) = varCoords
        chosenIDMD(j) = idmDistance
        chosenIDMV(j) = varDM
        chosenICMD(j) = icmDistance
        chosenICMV(j) = varCM
        chosenL1NW(j) = l1NegativeWeight
        chosenL2NW(j) = l2NegativeWeight
        chosenLiNW(j) = largestNegativeWeight
    end if
end select

minR2 = min(minR2,R2)
minErrorIMRR = min(minErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!
! Replacement Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(rsID)
case(1)
    if (any(abs(frame_weights(1:Ninterpolation_current)) < &
        weightThreshold)) then
        do i = Ninterpolation_current, 1, -1
            if (abs(frame_weights(i)) < weightThreshold) then
                if (i < Ninterpolation_current) &
                    reorderedIndexes(i) = &
                        reorderedIndexes(i+1)
                Ninterpolation_current = &
                    Ninterpolation_current - 1
            end if
        end do
        errorID = 2
    else
        if (Ninterpolation_current < Ninterpolation) then
            Ninterpolation_current = &
                Ninterpolation_current + 1
            reorderedIndexes(Ninterpolation_current) = k
        else
            reorderedIndexes(minloc(abs(&
                  frame_weights(1:Ninterpolation_current)))) = k
        end if
    
        if (k > Ninterpolation_pool) exit
        k = k + 1
    
        errorID = 3
    end if
case(2)
    if (Ninterpolation_current < Ninterpolation) then
        Ninterpolation_current = &
            Ninterpolation_current + 1
        reorderedIndexes(Ninterpolation_current) = k
        errorID = 2
    else
        if (minval(frame_weights(1:Ninterpolation_current))<0.0d0) then
            errorID = 2
        else
            errorID = 3
        end if
        reorderedIndexes(minloc(&
              frame_weights(1:Ninterpolation_current))) = k
    end if

    if (k > Ninterpolation_pool) exit
    k = k + 1
end select

end do

end do

close(6666)

FMTalpha = "(I8,1x,I3"
do i = 1, 3*Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,ES14.7"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//consolidatedErrorIMRRfile,&
          position="append")
write(6666,FMT=trim(adjustl(FMTalpha))) &
    NC, labelID, chosenErrorIMRR, chosenR1IMRR, chosenR2IMRR
close(6666)

!FMTalpha = "(I8,1x,I3"
!do i = 1, Nalpha
!    FMTalpha = trim(adjustl(FMTalpha))//",12(1x,ES14.7)"
!end do
!FMTalpha = trim(adjustl(FMTalpha))//")"
!
!open(6666,file=gridpath5//"other"//consolidatedErrorIMRRfile,&
!          position="append")
!write(6666,FMT=trim(adjustl(FMTalpha))) &
!    NC, labelID, chosenErrorIMRR, chosenR1IMRR, chosenR2IMRR, &
!        chosenL1NW, chosenL2NW, chosenLiNW, &
!        chosenIGD, chosenIGV, &
!        chosenIDMD, chosenIDMV, &
!        chosenICMD, chosenICMV
!close(6666)







if (output_flag) then

! For the sGDML comparison

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'inputs.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
do i = 1, Ninterpolation_pool
    write(6666,FMT="(I1)") Natoms
    write(6666,FMT="(F14.9)") potEbuffer1(i) * 4.3363d-2
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "H ", coordsbuffer1(:,1,i), &
        -gradientbuffer1(:,1,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "Br", coordsbuffer1(:,2,i), &
        -gradientbuffer1(:,2,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,3,i), &
        -gradientbuffer1(:,3,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "C ", coordsbuffer1(:,4,i), &
        -gradientbuffer1(:,4,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,5,i), &
        -gradientbuffer1(:,5,i) * 51.422114d0
end do
close(6666)

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'target.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
write(6666,FMT="(I1)") Natoms
write(6666,FMT="(F14.9)") V * 4.3363d-2
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "H ", coords1(:,1), &
    -gradient1(:,1) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "Br", coords1(:,2), &
    -gradient1(:,2) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,3), &
    -gradient1(:,3) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "C ", coords1(:,4), &
    -gradient1(:,4) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,5), &
    -gradient1(:,5) * 51.422114d0
close(6666)


! To make it more readable
maxErrorIMRR = min(maxErrorIMRR, 0.06d0)

Nbins = 20
errorIMRRbinwidth = (maxErrorIMRR - minErrorIMRR) / Nbins
allocate(errorIMRRbinning(Nbins,3,Nalpha))

errorIMRRbinning = 0
open(6666,file=gridpath5//"errorIMRRtest.dat")
do
    read(6666,iostat=iostate,FMT=*) &
        j,errorID,errorIMRR,igDistance,&
        varCoords,R1,R2, Ninterpolation_current,&
        largestNegativeWeight
!       totalWeight,largestNegativeWeight
    if (iostate /= 0) exit

    errorBin = floor((errorIMRR-minErrorIMRR)/errorIMRRbinwidth) + 1
    if (errorBin < 1) errorBin = 1
    if (errorBin > Nbins) errorBin = Nbins

    if (errorID == 1) &
    errorIMRRbinning(errorBin,1,j) = &
        errorIMRRbinning(errorBin,1,j) + 1
    if (errorID == 2) &
    errorIMRRbinning(errorBin,2,j) = &
        errorIMRRbinning(errorBin,2,j) + 1

    errorIMRRbinning(errorBin,3,j) = &
        errorIMRRbinning(errorBin,3,j) + 1

end do
close(6666)

maxErrorIMRRoccurence = maxval(&
        errorIMRRbinning(1:Nbins,1:3,1:Nalpha))

meanR1 = meanR1 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))
meanR2 = meanR2 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))

FMTalpha = "(ES14.7"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,I6"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"errorIMRRtest_binned1.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,1,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned2.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,2,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned3.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,3,1:Nalpha)
end do
close(6666)

deallocate(errorIMRRbinning)

!getIMRRerrors_counter = getIMRRerrors_counter + 1
getIMRRerrors_counter = NC

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRR',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
!write(6666,FMT="(A,F10.7,A)") &
!        'set label 3 "Efinal = ', Efinal, '" at screen 0.82,0.86'
!write(6666,*) 'set tmargin at screen 0.8'
!write(6666,*) 'set bmargin at screen 0.2'
!write(6666,*) 'set rmargin at screen 0.8'
!write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'Nalpha = ', Nalpha
write(6666,FMT=*) 'ymax = ', maxErrorIMRRoccurence
write(6666,FMT=*) 'xmin = ', max(minErrorIMRR-errorIMRRbinwidth,0.0d0)
write(6666,FMT=*) 'xmax = ', maxErrorIMRR+errorIMRRbinwidth

write(6666,FMT="(A)") 'set multiplot layout Nalpha,1 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0.1,0'
write(6666,FMT="(A)") 'set xrange [xmin:xmax]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set style fill transparent solid 1.0'
write(6666,FMT=*) 'set boxwidth ', errorIMRRbinwidth

write(6666,FMT="(A)") 'unset xlabel'
write(6666,FMT="(A)") 'unset xtics'
do i = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1))
    end if
    write(6666,FMT="(A,F7.3,A)") &
            'set ylabel "P_{1.0e',log10(alpha_ratio),'}"'
    if (i == Nalpha) then
        write(6666,FMT="(A)") 'set xtics'
        write(6666,FMT="(A)") 'set xlabel "IMRR Error"'
    end if
    write(6666,FMT="(A)") 'unset arrow'
    write(6666,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', chosenErrorIMRR(i),&
        ',graph 0 to ', chosenErrorIMRR(i),&
        ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
    write(6666,FMT="(A)") 'unset label 3'
    write(6666,FMT="(A,E16.6,A)")&
        'set label 3 "', chosenErrorIMRR(i),&
        ' (E_h/a_0)" at graph 0.8, graph 0.9 '//&
        ' font ",12" front'
    write(6666,FMT="(A,I1,A)") &
            'plot "'//gridpath5//'errorIMRRtest_binned3.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "red",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned2.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "green",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned1.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "blue"'
end do
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigd',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', igDistanceMax
write(6666,FMT=*) 'ymax = ', maxErrorIMRR

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center Distance (A)"'
write(6666,FMT="(A)") 'set ylabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:3:1 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigv',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', varCoordsMax
write(6666,FMT=*) 'ymax = ', maxErrorIMRR

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center Variance (A)"'
write(6666,FMT="(A)") 'set ylabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 5:3:1 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigdN',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', igDistanceMax
write(6666,FMT=*) 'ymax = ', maxNinterpolation_current
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center distance (A)"'
write(6666,FMT="(A)") 'set ylabel "Ninterpolation"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:8:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)


open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRr1r2',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', sqrt(maxR2*meanR2)
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

!write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set xrange [:xmax*1.1]'
write(6666,FMT="(A)") 'set logscale x'

write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "R2 (A^2)"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 6:7:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRw',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', 1.5d0
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'
write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'

write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "Summed Negative Weight"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 6:9:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This part looks at the errors
! of individual frames

maxErrorIMRR = 0.0d0
open(6666,file=gridpath5//temporaryfile1)
do step = 1, Ninterpolation_pool
    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    gradient2 = matmul(&
            U,gradientbuffer1(:,:,step))

    errorIMRR = sqrt(sum((gradient2 - gradient1)**2)/Natoms)
    maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

    write(6666,FMT="(3(ES14.7,1x),ES14.7)") &
        new_SIs(1:3), errorIMRR
end do
close(6666)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,1000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorInputs',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.96'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.92'

write(6666,FMT="(A)") 'set multiplot layout 1,3 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0,0.1'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT=*) 'ymax = ', maxErrorIMRR
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set xrange [0:]'
write(6666,FMT="(A,F7.3,A)") &
        'set ylabel "Difference in Energy Gradient (E_h/a_0)"'
write(6666,FMT="(A)") 'set xlabel "RMSD (A)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//temporaryfile1//'" '//&
        'u 1:4 w p ps 2 pt 7'
write(6666,FMT="(A)") 'unset ytics'
write(6666,FMT="(A)") 'unset ylabel'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set xlabel "CMD w/ Charges"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//temporaryfile1//'" '//&
        'u 2:4 w p ps 2 pt 7'
write(6666,FMT="(A)") 'set xlabel "CMD w/o Charges"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//temporaryfile1//'" '//&
        'u 3:4 w p ps 2 pt 7'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end if

deallocate(reorderedIndexes,&
           frame_weights,&
           outputCLS,restraints,&
           restraint_values,&
           inputCLS2,libUgradients,libenergies)

return
end subroutine getIMRRerrors3


subroutine getIMRRerrors4(Q,PDOT,labelID)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use VARIABLES
use SIMILARITY
use interactMultipleGrids
implicit none

COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC
COMMON/FORCES/NATOMS

integer :: NATOMS, NC
double precision :: NTZ,NT
double precision :: ISEED0
double precision :: TIME
double precision :: T,V,H

integer :: Ncoords

double precision,dimension(3*Natoms),intent(in) :: Q, PDOT
integer,intent(in) :: labelID

integer :: Ninterpolation_pool
integer, allocatable :: reorderedIndexes(:)

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)
real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp), allocatable :: inputCLS2(:,:)
real(dp), allocatable :: libUgradients(:,:,:)
real(dp), allocatable :: libenergies(:)

real(dp),dimension(3,3) :: U
real(dp),dimension(NSIs) :: new_SIs

real(dp),dimension(3,Natoms) :: coords1,coords2
real(dp),dimension(3,Natoms) :: gradient1,gradient2
real(dp) :: min_rmsd, errorIMRR, Efinal
real(dp) :: minErrorIMRR, maxErrorIMRR
real(dp) :: errorIMRRbinwidth
integer :: maxErrorIMRRoccurence

character(400) :: FMTalpha, filetext
integer :: Nbins, errorBin
integer,allocatable :: errorIMRRbinning(:,:,:)

real(dp),allocatable :: CMvec1(:),CMvec2(:),meanCoords(:)
real(dp) :: varCoords
real(dp) :: igDistance, igDistanceMax, varCoordsMax
integer :: maxNinterpolation_current

real(dp) :: R1, R2, maxR1, maxR2, meanR1, meanR2
real(dp) :: maxErrorIMRR_cb = 0.06d0

real(dp) :: weightThreshold = 5.0d-3
integer :: errorID, Ninterpolation_current

real(dp) :: totalWeight,largestNegativeWeight

character(len=*),parameter :: consolidatedErrorIMRRfile = &
                  "consolidatedErrorIMRR.dat"
real(dp) :: minR2
real(dp),dimension(Nalpha) :: chosenErrorIMRR, &
                chosenR1IMRR, chosenR2IMRR
real(dp) :: potE2

integer*8 :: dummyINT1, dummyINT2, dummyINT3

integer, save :: getIMRRerrors_counter = 0
integer :: iostate, step, i, j, k

logical :: output_flag
logical :: potE_flag = .false.

output_flag = (modulo(NC,100)==1)

Ncoords = (NATOMS * (NATOMS-1)) / 2

coords1 = reshape(Q,(/3,Natoms/))
gradient1 = reshape(PDOT,(/3,Natoms/))

call getVarsHBrCO2(coords1,Natoms,vals,&
               Nvar,BOND_LABELLING_DATA)

!print *, "whaaat 1"
call checkState_PCM(&
        vals,coords1,&
        gradient2,&
        dummyINT1,dummyINT2,dummyINT3)
!print *, "whaaat 2"

! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (getIMRRerrors_counter > 1) return
!if (any(vals > 6.0d0)) return
 if (Ninterpolation < 20) return
!if (Ninterpolation < Ninterpolation_cap+5) return
! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Ninterpolation_pool = Ninterpolation
Ninterpolation_pool = min(Ninterpolation,Ninterpolation_cap)
Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)
!       20)

allocate(reorderedIndexes(Ninterpolation),&
         frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         inputCLS2(Ncoords+Ninterpolation,Ninterpolation),&
         libUgradients(3,NATOMS,Ninterpolation),&
         libenergies(Ninterpolation),&
         CMvec1(Ncoords),CMvec2(Ncoords),meanCoords(Ncoords))

call setTarget(coords1)
call getCoulombMatrixVector(&
        coords1,CMvec1)

maxNinterpolation_current = 0
maxR1 = 0.0d0
maxR2 = 0.0d0
meanR1 = 0.0d0
meanR2 = 0.0d0
minErrorIMRR = 1.0d9
maxErrorIMRR = 0.0d0
igDistanceMax = 0.0d0
varCoordsMax = 0.0d0

open(6666,file=gridpath5//"errorIMRRtest.dat")

do j = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1))
    end if

    errorID = 3
    Ninterpolation_current = Ninterpolation !5
    k = Ninterpolation_current + 1
    do i = 1, Ninterpolation_current
        reorderedIndexes(i) = i
    end do

    minR2 = 1.0d9
    chosenErrorIMRR(j) = 1.0d9

do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part only looks at the inputs

do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)

    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    libUgradients(:,:,i) = matmul(&
            U,gradientbuffer1(:,:,step))
    libenergies(i) = &
        potEbuffer1(step)

    call getCoulombMatrixVector(&
        coords2,CMvec2)

    inputCLS2(1:Ncoords,i) = &
            CMvec2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CMvec2 - CMvec1
    inputCLS2(Ncoords+i,:) = 0.0d0
    inputCLS2(Ncoords+i,i) = alpha_ratio * &
            new_SIs(Nsort)**2
end do

meanCoords = 0.0d0
do i = 1, Ninterpolation_current
    meanCoords = meanCoords + &
                 inputCLS2(1:Ncoords,i)
end do
meanCoords = meanCoords / Ninterpolation_current
igDistance = sqrt(sum(meanCoords**2)/Natoms)

varCoords = 0.0d0
do i = 1, Ninterpolation_current
    varCoords = varCoords + &
                sum((inputCLS2(1:Ncoords,i)-meanCoords)**2)
end do
varCoords = sqrt(varCoords/Ninterpolation_current)

igDistanceMax = max(igDistance,igDistanceMax)
varCoordsMax = max(varCoords,varCoordsMax)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outputCLS = 0.0d0
outputCLS(1:Ncoords) = CMvec1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
restraints = 1.0d0
restraint_values = 1.0d0

!call CLS2(inputCLS2(1:Ncoords+Ninterpolation_current,&
!                    1:Ninterpolation_current),&
!          Ncoords+Ninterpolation_current,&
!          Ninterpolation_current,&
!          restraints(1:1,1:Ninterpolation_current),&
!          1,restraint_values,&
!          outputCLS(1:Ncoords+Ninterpolation_current),&
!          frame_weights(1:Ninterpolation_current))
call CLSwK(inputCLS2(1:Ncoords,&
                    1:Ninterpolation_current),&
          Ncoords,&
          Ninterpolation_current,&
          restraints(1:1,1:Ninterpolation_current),&
          1,restraint_values,&
          outputCLS(1:Ncoords),&
          frame_weights(1:Ninterpolation_current))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part looks at the interpolation

if (all(frame_weights(1:Ninterpolation_current) == 0.0d0)) then
    if (Ninterpolation_current < Ninterpolation) &
        Ninterpolation_current = &
            Ninterpolation_current + 1

    reorderedIndexes(Ninterpolation_current) = k
    k = k + 1

    if (k > Ninterpolation_pool) exit
    cycle
end if

gradient2 = 0.0d0
potE2 = 0.0d0
R2 = 0.0d0
do i = 1, Ninterpolation_current
    gradient2 = gradient2 + frame_weights(i)*&
                libUgradients(:,:,i)
    potE2 = potE2 + frame_weights(i)*&
                libenergies(i)
    R2 = R2 + abs(frame_weights(i)) * &
            inputCLS2(Ncoords+i,i)
end do
R2 = R2 / alpha_ratio
R1 = sqrt(sum(matmul(&
    inputCLS2(1:Ncoords,1:Ninterpolation_current),&
    reshape(frame_weights(1:Ninterpolation_current),&
            (/ Ninterpolation_current, 1 /)))**2))

totalWeight = sum(frame_weights(1:Ninterpolation_current))
largestNegativeWeight = abs(sum(frame_weights(1:Ninterpolation_current),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0))


if (potE_flag) then
    errorIMRR = abs(V - potE2)
else
    errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (k == Ninterpolation_pool) errorID = 1

open(6667,file=gridpath5//"all"//consolidatedErrorIMRRfile,&
          position="append")
write(6667,FMT="(I8,1x,I3,3(1x,ES14.7))") &
    NC, labelID, errorIMRR, R1, R2
close(6667)

!write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2)") &
if (output_flag) &
write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2,1x,ES14.7)") &
                  j,errorID,&
                  errorIMRR,igDistance,varCoords,&
                  R1, R2, Ninterpolation_current,&
                  largestNegativeWeight
!                 totalWeight,largestNegativeWeight

maxNinterpolation_current = max(&
             Ninterpolation_current,&
          maxNinterpolation_current)
maxR1 = max(R1,maxR1)
maxR2 = max(R2,maxR2)
meanR1 = meanR1 + R1
meanR2 = meanR2 + R2
maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!
! Selection Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(ssID)
case(1)
    if (R2 < minR2) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
    end if
case(2)
    if (errorIMRR < chosenErrorIMRR(j)) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
    end if
end select

minR2 = min(minR2,R2)
minErrorIMRR = min(minErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!
! Replacement Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(rsID)
case(1)
    if (any(abs(frame_weights(1:Ninterpolation_current)) < &
        weightThreshold)) then
        do i = Ninterpolation_current, 1, -1
            if (abs(frame_weights(i)) < weightThreshold) then
                if (i < Ninterpolation_current) &
                    reorderedIndexes(i) = &
                        reorderedIndexes(i+1)
                Ninterpolation_current = &
                    Ninterpolation_current - 1
            end if
        end do
        errorID = 2
    else
        if (Ninterpolation_current < Ninterpolation) then
            Ninterpolation_current = &
                Ninterpolation_current + 1
            reorderedIndexes(Ninterpolation_current) = k
        else
            reorderedIndexes(minloc(abs(&
                  frame_weights(1:Ninterpolation_current)))) = k
        end if
    
        if (k > Ninterpolation_pool) exit
        k = k + 1
    
        errorID = 3
    end if
case(2)
    if (Ninterpolation_current < Ninterpolation) then
        Ninterpolation_current = &
            Ninterpolation_current + 1
        reorderedIndexes(Ninterpolation_current) = k
        errorID = 2
    else
        if (minval(frame_weights(1:Ninterpolation_current))<0.0d0) then
            errorID = 2
        else
            errorID = 3
        end if
        reorderedIndexes(minloc(&
              frame_weights(1:Ninterpolation_current))) = k
    end if

    if (k > Ninterpolation_pool) exit
    k = k + 1
end select

end do

end do

close(6666)


FMTalpha = "(I8,1x,I3"
do i = 1, 3*Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,ES14.7"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//consolidatedErrorIMRRfile,&
          position="append")
write(6666,FMT=trim(adjustl(FMTalpha))) &
    NC, labelID, chosenErrorIMRR, chosenR1IMRR, chosenR2IMRR
close(6666)








if (output_flag) then

! For the sGDML comparison

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'inputs.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
do i = 1, Ninterpolation_pool
    write(6666,FMT="(I1)") Natoms
    write(6666,FMT="(F14.9)") potEbuffer1(i) * 4.3363d-2
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "H ", coordsbuffer1(:,1,i), &
        -gradientbuffer1(:,1,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "Br", coordsbuffer1(:,2,i), &
        -gradientbuffer1(:,2,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,3,i), &
        -gradientbuffer1(:,3,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "C ", coordsbuffer1(:,4,i), &
        -gradientbuffer1(:,4,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,5,i), &
        -gradientbuffer1(:,5,i) * 51.422114d0
end do
close(6666)

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'target.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
write(6666,FMT="(I1)") Natoms
write(6666,FMT="(F14.9)") V * 4.3363d-2
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "H ", coords1(:,1), &
    -gradient1(:,1) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "Br", coords1(:,2), &
    -gradient1(:,2) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,3), &
    -gradient1(:,3) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "C ", coords1(:,4), &
    -gradient1(:,4) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,5), &
    -gradient1(:,5) * 51.422114d0
close(6666)

Nbins = 20
errorIMRRbinwidth = (maxErrorIMRR - minErrorIMRR) / Nbins
allocate(errorIMRRbinning(Nbins,3,Nalpha))

errorIMRRbinning = 0
open(6666,file=gridpath5//"errorIMRRtest.dat")
do
    read(6666,iostat=iostate,FMT=*) &
        j,errorID,errorIMRR,igDistance,&
        varCoords,R1,R2, Ninterpolation_current,&
        largestNegativeWeight
!       totalWeight,largestNegativeWeight
    if (iostate /= 0) exit

    errorBin = floor((errorIMRR-minErrorIMRR)/errorIMRRbinwidth) + 1
    if (errorBin < 1) errorBin = 1
    if (errorBin > Nbins) errorBin = Nbins

    if (errorID == 1) &
    errorIMRRbinning(errorBin,1,j) = &
        errorIMRRbinning(errorBin,1,j) + 1
    if (errorID == 2) &
    errorIMRRbinning(errorBin,2,j) = &
        errorIMRRbinning(errorBin,2,j) + 1

    errorIMRRbinning(errorBin,3,j) = &
        errorIMRRbinning(errorBin,3,j) + 1

end do
close(6666)

maxErrorIMRRoccurence = maxval(&
        errorIMRRbinning(1:Nbins,1:3,1:Nalpha))

meanR1 = meanR1 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))
meanR2 = meanR2 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))

FMTalpha = "(ES14.7"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,I6"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"errorIMRRtest_binned1.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,1,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned2.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,2,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned3.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,3,1:Nalpha)
end do
close(6666)

deallocate(errorIMRRbinning)

!getIMRRerrors_counter = getIMRRerrors_counter + 1
getIMRRerrors_counter = NC

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRR',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
!write(6666,FMT="(A,F10.7,A)") &
!        'set label 3 "Efinal = ', Efinal, '" at screen 0.82,0.86'
!write(6666,*) 'set tmargin at screen 0.8'
!write(6666,*) 'set bmargin at screen 0.2'
!write(6666,*) 'set rmargin at screen 0.8'
!write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'Nalpha = ', Nalpha
write(6666,FMT=*) 'ymax = ', maxErrorIMRRoccurence
write(6666,FMT=*) 'xmin = ', max(minErrorIMRR-errorIMRRbinwidth,0.0d0)
write(6666,FMT=*) 'xmax = ', maxErrorIMRR+errorIMRRbinwidth

write(6666,FMT="(A)") 'set multiplot layout Nalpha,1 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0.1,0'
write(6666,FMT="(A)") 'set xrange [xmin:xmax]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set style fill transparent solid 1.0'
write(6666,FMT=*) 'set boxwidth ', errorIMRRbinwidth

write(6666,FMT="(A)") 'unset xlabel'
write(6666,FMT="(A)") 'unset xtics'
do i = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1))
    end if
    write(6666,FMT="(A,F7.3,A)") &
            'set ylabel "P_{1.0e',log10(alpha_ratio),'}"'
    if (i == Nalpha) then
        write(6666,FMT="(A)") 'set xtics'
        write(6666,FMT="(A)") 'set xlabel "IMRR Error"'
    end if
    write(6666,FMT="(A,I1,A)") &
            'plot "'//gridpath5//'errorIMRRtest_binned3.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "red",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned2.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "green",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned1.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "blue"'
end do
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigd',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', igDistanceMax
write(6666,FMT=*) 'ymax = ', maxErrorIMRR

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center Distance (A)"'
write(6666,FMT="(A)") 'set ylabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:3:1 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigv',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', varCoordsMax
write(6666,FMT=*) 'ymax = ', maxErrorIMRR

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center Variance (A)"'
write(6666,FMT="(A)") 'set ylabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 5:3:1 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigdN',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', igDistanceMax
write(6666,FMT=*) 'ymax = ', maxNinterpolation_current
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center distance (A)"'
write(6666,FMT="(A)") 'set ylabel "Ninterpolation"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:8:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)


open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRr1r2',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', sqrt(maxR2*meanR2)
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

!write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set xrange [:xmax*1.1]'
write(6666,FMT="(A)") 'set logscale x'

write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "R2 (A^2)"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 6:7:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRw',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', 1.5d0
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'
write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'

write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "Summed Negative Weight"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 6:9:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

end if

deallocate(reorderedIndexes)
deallocate(frame_weights)
deallocate(outputCLS,restraints)
deallocate(restraint_values)
deallocate(inputCLS2,libUgradients,libenergies)
deallocate(CMvec1,CMvec2,meanCoords)

return
end subroutine getIMRRerrors4


subroutine getIMRRerrors5(Q,PDOT,labelID)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use VARIABLES
use SIMILARITY
use interactMultipleGrids
implicit none

COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC
COMMON/FORCES/NATOMS

integer :: NATOMS, NC
double precision :: NTZ,NT
double precision :: ISEED0
double precision :: TIME
double precision :: T,V,H

integer :: Ncoords

double precision,dimension(3*Natoms),intent(in) :: Q, PDOT
integer,intent(in) :: labelID

integer :: Ninterpolation_pool
integer, allocatable :: reorderedIndexes(:)

real(dp), allocatable :: Omega(:,:),outputMFN(:,:),inputMFN(:,:)
real(dp), allocatable :: libUgradients(:,:)
real(dp), allocatable :: libenergies(:)

integer :: Nfeatures
real(dp), allocatable :: inputFeatures(:,:), Feature1(:)
real(dp), allocatable :: inputJacobians(:,:), Jacobian1(:,:)

real(dp),dimension(3,3) :: U
real(dp),dimension(NSIs) :: new_SIs

real(dp),dimension(3,Natoms) :: coords1,coords2
real(dp),dimension(3,Natoms) :: gradient1,gradient2
real(dp) :: min_rmsd, errorIMRR, Efinal
real(dp) :: minErrorIMRR, maxErrorIMRR
real(dp) :: errorIMRRbinwidth
integer :: maxErrorIMRRoccurence

character(400) :: FMTalpha, filetext
integer :: Nbins, errorBin
integer,allocatable :: errorIMRRbinning(:,:,:)

integer :: maxNinterpolation_current

real(dp) :: R1, R2, maxR1, maxR2, meanR1, meanR2
real(dp) :: maxErrorIMRR_cb = 0.06d0

real(dp) :: weightThreshold = 5.0d-3
integer :: errorID, Ninterpolation_current

real(dp) :: totalWeight,largestNegativeWeight

character(len=*),parameter :: consolidatedErrorIMRRfile = &
                  "consolidatedErrorIMRR.dat"
real(dp) :: minR2
real(dp),dimension(Nalpha) :: chosenErrorIMRR, &
                chosenR1IMRR, chosenR2IMRR
real(dp) :: potE2

integer*8 :: dummyINT1, dummyINT2, dummyINT3

integer, save :: getIMRRerrors_counter = 0
integer :: iostate, step, i, j, k, l, m

logical :: output_flag

output_flag = (modulo(NC,100)==1)

Ncoords = 3 * NATOMS

coords1 = reshape(Q,(/3,Natoms/))
gradient1 = reshape(PDOT,(/3,Natoms/))

call getVarsHBrCO2(coords1,Natoms,vals,&
               Nvar,BOND_LABELLING_DATA)

!print *, "whaaat 1"
call checkState_PCM(&
        vals,coords1,&
        gradient2,&
        dummyINT1,dummyINT2,dummyINT3)
!print *, "whaaat 2"

! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (getIMRRerrors_counter > 1) return
!if (any(vals > 6.0d0)) return
if (Ninterpolation < 3) return !20) return !< 25) return
!if (Ninterpolation < Ninterpolation_cap+5) return
! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This one takes a long time so
! we cap the pool
Ninterpolation_pool = min(Ninterpolation,Ninterpolation_cap) !Ninterpolation
Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)
!       20)

allocate(reorderedIndexes(Ninterpolation),&
         Omega(Ncoords,Ncoords*Ninterpolation),&
         outputMFN(Ncoords,Ncoords*Ninterpolation+1),&
         inputMFN(Ncoords*Ninterpolation,Ncoords*Ninterpolation+1),&
         libUgradients(Ncoords,Ninterpolation),&
         libenergies(Ninterpolation))

Nfeatures = (Natoms * (Natoms-1) ) / 2
allocate(inputJacobians(Nfeatures,Ncoords*Ninterpolation),&
         Jacobian1(Nfeatures,Ncoords),&
         inputFeatures(Nfeatures,Ninterpolation),&
         Feature1(Nfeatures))

call getCoulombMatrixVector(coords1,Feature1)
call getCoulombMatrixGradient(coords1,Jacobian1)

call setTarget(coords1)

maxNinterpolation_current = 0
maxR1 = 0.0d0
maxR2 = 0.0d0
meanR1 = 0.0d0
meanR2 = 0.0d0
minErrorIMRR = 1.0d9
maxErrorIMRR = 0.0d0

open(6666,file=gridpath5//"errorIMRRtest.dat")

do j = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1))
    end if

    errorID = 3
    Ninterpolation_current = Ninterpolation_pool !5
    k = Ninterpolation_current + 1
    do i = 1, Ninterpolation_current
        reorderedIndexes(i) = i
    end do

    minR2 = 1.0d9
    chosenErrorIMRR(j) = 1.0d9

do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part only looks at the inputs

inputMFN = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)

    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    libUgradients(:,i) = reshape(matmul(&
            U,gradientbuffer1(:,:,step)),(/Ncoords/))
    libenergies(i) = &
        potEbuffer1(step)

    do m = 1, 1 !Nfeatures
        inputMFN((i-1)*Ncoords+1:i*Ncoords,m) = &
                reshape(coords2,(/Ncoords/))
    end do

    do m = (i-1)*Ncoords+1, i*Ncoords
!       inputMFN(m,m+Nfeatures) = alpha_ratio * &
        inputMFN(m,m+1) = alpha_ratio * &
                new_SIs(Nsort)**2
    end do

    call getCoulombMatrixVector(&
        coords2,inputFeatures(1:Nfeatures,i))
    call getCoulombMatrixGradient(&
        coords2,inputJacobians(1:Nfeatures,&
            (i-1)*Ncoords+1:i*Ncoords))
end do

!print *, ""
!print *, "attempting MFN"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! test
!Omega = 0.0d0
!do i = 1, 3*Natoms
!    Omega(i,i) = 1.0d0
!end do
!
!gradient2 = 0.0d0
!do i = 1, Ninterpolation_current
!    gradient2 = gradient2 + reshape(matmul(&
!            Omega(1:Ncoords,(i-1)*Ncoords+1:i*Ncoords),&
!            libUgradients(1:Ncoords,i:i)),(/3,NATOMS/))
!end do
!R1 = sqrt(sum((  matmul(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),&
!            inputMFN(1:Ncoords*Ninterpolation_current,&
!                     1:1)) - reshape(Q,(/Ncoords,1/)))**2)/Natoms)
!R2 = sum(  matmul(abs(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)),&
!           inputMFN(1:Ncoords*Ninterpolation_current,&
!                    2:Ncoords*Ninterpolation_current+1))) / alpha_ratio
!
!errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)
!
!print *, "(test) error,R1,R2:"
!print *, errorIMRR, R1, R2
!
!i = 3*Natoms
!m = 3*Natoms*Ninterpolation_current
!l = 3*Natoms*Ninterpolation_current + 1
!outputMFN(1:i,1) = Q(1:Ncoords)
!outputMFN(1:i,2:l) = 0.0d0
!
!print *, "(test) ||I_3N A - B||_F:"
!!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:l)) - &
!!               outputMFN(1:i,1:l))**2)
!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:1)) - &
!               outputMFN(1:i,1:1))**2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!i = 3*Natoms
!m = 3*Natoms*Ninterpolation_current
!l = 3*Natoms*Ninterpolation_current + 1
!!call mFN(Omega(1:i,1:m),i,m,&
!!         inputMFN(1:m,1:l),l,&
!!         outputMFN(1:i,1:l))
!
!!call sVVF(inputMFN(1:m,1),i,Ninterpolation_current,&
!!          Q(1:i),Omega(1:i,1:m))
!
!call sVVFwJ(inputFeatures(1:Nfeatures,1:Ninterpolation_current),&
!            inputJacobians(1:Nfeatures,1:Ncoords*Ninterpolation_current),&
!            Ncoords,Nfeatures,Ninterpolation_current,&
!            inputFeatures(1:Nfeatures,1),inputJacobians(1:Nfeatures,1:Ncoords),Omega(1:i,1:m))
!
!gradient2 = 0.0d0
!do i = 1, Ninterpolation_current
!    gradient2 = gradient2 + reshape(matmul(&
!            Omega(1:Ncoords,(i-1)*Ncoords+1:i*Ncoords),&
!            libUgradients(1:Ncoords,i:i)),(/3,NATOMS/))
!end do
!R1 = sqrt(sum((  matmul(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),&
!            inputMFN(1:Ncoords*Ninterpolation_current,&
!                     1:1)) - inputMFN(1:Ncoords,1:1))**2)/Natoms)
!R2 = sum(  matmul(abs(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)),&
!           inputMFN(1:Ncoords*Ninterpolation_current,&
!                    2:Ncoords*Ninterpolation_current+1))) / alpha_ratio
!
!!errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)
!errorIMRR = sqrt(sum((reshape(libUgradients(1:Ncoords,1),(/3,Natoms/))-gradient2)**2)/Natoms)
!
!print *, "(test for reproduction of 1) error,R1,R2:"
!print *, errorIMRR, R1, R2
!
!i = 3*Natoms
!m = 3*Natoms*Ninterpolation_current
!l = 3*Natoms*Ninterpolation_current + 1
!outputMFN(1:i,1) = inputMFN(1:Ncoords,1)
!outputMFN(1:i,2:l) = 0.0d0
!
!print *, "(test for reproduction of 1) ||Omega A - B||_F:"
!!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:l)) - &
!!               outputMFN(1:i,1:l))**2)
!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:1)) - &
!               outputMFN(1:i,1:1))**2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!i = 3*Natoms
!m = 3*Natoms*Ninterpolation_current
!l = 3*Natoms*Ninterpolation_current + 1
!outputMFN(1:i,1) = Q(1:Ncoords)
!outputMFN(1:i,2:l) = 0.0d0
!
!call sVVF(inputMFN(1:m,1),i,Ninterpolation_current,&
!          Q(1:i),Omega(1:i,1:m))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! This part looks at the interpolation
!
!if (all(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current) == 0.0d0)) then
!    if (Ninterpolation_current < Ninterpolation) &
!        Ninterpolation_current = &
!            Ninterpolation_current + 1
!
!    reorderedIndexes(Ninterpolation_current) = k
!    k = k + 1
!
!    print *, "what???"
!    if (k > Ninterpolation_pool) exit
!    cycle
!end if
!
!gradient2 = 0.0d0
!do i = 1, Ninterpolation_current
!    gradient2 = gradient2 + reshape(matmul(&
!            Omega(1:Ncoords,(i-1)*Ncoords+1:i*Ncoords),&
!            libUgradients(1:Ncoords,i:i)),(/3,NATOMS/))
!end do
!R1 = sqrt(sum((  matmul(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),&
!            inputMFN(1:Ncoords*Ninterpolation_current,&
!                     1:1)) - reshape(Q,(/Ncoords,1/)))**2)/Natoms)
!R2 = sum(  matmul(abs(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)),&
!           inputMFN(1:Ncoords*Ninterpolation_current,&
!                    2:Ncoords*Ninterpolation_current+1))) / alpha_ratio
!
!totalWeight = sum(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current))
!largestNegativeWeight = abs(sum(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),& 
!                       mask=Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)<0.0d0))
!
!
!errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)
!
!print *, "error,R1,R2:"
!print *, errorIMRR, R1, R2
!
!i = 3*Natoms
!m = 3*Natoms*Ninterpolation_current
!l = 3*Natoms*Ninterpolation_current + 1
!
!print *, "||Omega A - B||_F:"
!!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:l)) - &
!!               outputMFN(1:i,1:l))**2)
!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:1)) - &
!               outputMFN(1:i,1:1))**2)

!print *, ""
!print *, "Omega:"
!do i = 1, Ninterpolation_current
!    print *, "   O",i
!    do m = 1, Ncoords
!        print *, Omega(m,(i-1)*Ncoords+1:i*Ncoords)
!    end do
!    print *, ""
!end do
!print *, ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

i = 3*Natoms
m = 3*Natoms*Ninterpolation_current
l = 3*Natoms*Ninterpolation_current + 1
outputMFN(1:i,1) = Q(1:Ncoords)
outputMFN(1:i,2:l) = 0.0d0

call sVVFwJ(inputFeatures(1:Nfeatures,1:Ninterpolation_current),&
            inputJacobians(1:Nfeatures,1:Ncoords*Ninterpolation_current),&
            Ncoords,Nfeatures,Ninterpolation_current,&
            Feature1(1:Nfeatures),Jacobian1(1:Nfeatures,1:Ncoords),Omega(1:i,1:m))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part looks at the interpolation

if (all(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current) == 0.0d0)) then
    if (Ninterpolation_current < Ninterpolation) &
        Ninterpolation_current = &
            Ninterpolation_current + 1

    reorderedIndexes(Ninterpolation_current) = k
    k = k + 1

    print *, "what???"
    if (k > Ninterpolation_pool) exit
    cycle
end if

gradient2 = 0.0d0
do i = 1, Ninterpolation_current
    gradient2 = gradient2 + reshape(matmul(&
            Omega(1:Ncoords,(i-1)*Ncoords+1:i*Ncoords),&
            libUgradients(1:Ncoords,i:i)),(/3,NATOMS/))
end do
R1 = sqrt(sum((  matmul(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),&
            inputMFN(1:Ncoords*Ninterpolation_current,&
                     1:1)) - reshape(Q,(/Ncoords,1/)))**2)/Natoms)
R2 = sum(  matmul(abs(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)),&
           inputMFN(1:Ncoords*Ninterpolation_current,&
                    2:Ncoords*Ninterpolation_current+1))) / alpha_ratio

totalWeight = sum(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current))
largestNegativeWeight = abs(sum(Omega(1:Ncoords,1:Ncoords*Ninterpolation_current),& 
                       mask=Omega(1:Ncoords,1:Ncoords*Ninterpolation_current)<0.0d0))


errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)

print *, "error,R1,R2:"
print *, errorIMRR, R1, R2
!
!i = 3*Natoms
!m = 3*Natoms*Ninterpolation_current
!l = 3*Natoms*Ninterpolation_current + 1
!
!print *, "||Omega A - B||_F:"
!!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:l)) - &
!!               outputMFN(1:i,1:l))**2)
!print *, sum( (matmul(Omega(1:i,1:m),inputMFN(1:m,1:1)) - &
!               outputMFN(1:i,1:1))**2)
!
!!print *, ""
!!print *, "Omega:"
!!do i = 1, Ninterpolation_current
!!    print *, "   O",i
!!    do m = 1, Ncoords
!!        print *, Omega(m,(i-1)*Ncoords+1:i*Ncoords)
!!    end do
!!    print *, ""
!!end do
!!print *, ""
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(6667,file=gridpath5//"all"//consolidatedErrorIMRRfile,&
          position="append")
write(6667,FMT="(I8,1x,I3,3(1x,ES14.7))") &
    NC, labelID, errorIMRR, R1, R2
close(6667)

if (k == Ninterpolation_pool) errorID = 1

!write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2)") &
if (output_flag) &
write(6666,FMT="(I2,1x,I3,3(1x,ES14.7),1x,I2,1x,ES14.7)") &
                  j,errorID,&
                  errorIMRR,&
                  R1, R2, Ninterpolation_current,&
                  largestNegativeWeight

maxNinterpolation_current = max(&
             Ninterpolation_current,&
          maxNinterpolation_current)
maxR1 = max(R1,maxR1)
maxR2 = max(R2,maxR2)
meanR1 = meanR1 + R1
meanR2 = meanR2 + R2
maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!
! Selection Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(ssID)
case(1)
    if (R2 < minR2) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
    end if
case(2)
    if (errorIMRR < chosenErrorIMRR(j)) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
    end if
end select

minR2 = min(minR2,R2)
minErrorIMRR = min(minErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!
! Replacement Strategy
!!!!!!!!!!!!!!!!!!!!!!!

!select case(rsID)
!case(1)
!    if (any(abs(frame_weights(1:Ninterpolation_current)) < &
!        weightThreshold)) then
!        do i = Ninterpolation_current, 1, -1
!            if (abs(frame_weights(i)) < weightThreshold) then
!                if (i < Ninterpolation_current) &
!                    reorderedIndexes(i) = &
!                        reorderedIndexes(i+1)
!                Ninterpolation_current = &
!                    Ninterpolation_current - 1
!            end if
!        end do
!        errorID = 2
!    else
!        if (Ninterpolation_current < Ninterpolation) then
!            Ninterpolation_current = &
!                Ninterpolation_current + 1
!            reorderedIndexes(Ninterpolation_current) = k
!        else
!            reorderedIndexes(minloc(abs(&
!                  frame_weights(1:Ninterpolation_current)))) = k
!        end if
!    
!        if (k > Ninterpolation_pool) exit
!        k = k + 1
!    
!        errorID = 3
!    end if
!case(2)
!    if (Ninterpolation_current < Ninterpolation) then
!        Ninterpolation_current = &
!            Ninterpolation_current + 1
!        reorderedIndexes(Ninterpolation_current) = k
!        errorID = 2
!    else
!        if (minval(frame_weights(1:Ninterpolation_current))<0.0d0) then
!            errorID = 2
!        else
!            errorID = 3
!        end if
!        reorderedIndexes(minloc(&
!              frame_weights(1:Ninterpolation_current))) = k
!    end if
!
!    if (k > Ninterpolation_pool) exit
!    k = k + 1
!end select

!!!!!!!!!!!!!!!!!!!!!!!!!!
! (arbitrary replacement)
if (Ninterpolation_current < Ninterpolation) then
    Ninterpolation_current = &
        Ninterpolation_current + 1
    reorderedIndexes(Ninterpolation_current) = k
    errorID = 2
else
    reorderedIndexes(Ninterpolation_current) = k
    errorID = 3
end if
if (k > Ninterpolation_pool) exit
k = k + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!

end do

end do

close(6666)


FMTalpha = "(I8,1x,I3"
do i = 1, 3*Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,ES14.7"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//consolidatedErrorIMRRfile,&
          position="append")
write(6666,FMT=trim(adjustl(FMTalpha))) &
    NC, labelID, chosenErrorIMRR, chosenR1IMRR, chosenR2IMRR
close(6666)








if (output_flag) then

! For the sGDML comparison

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'inputs.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
do i = 1, Ninterpolation_pool
    write(6666,FMT="(I1)") Natoms
    write(6666,FMT="(F14.9)") potEbuffer1(i) * 4.3363d-2
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "H ", coordsbuffer1(:,1,i), &
        -gradientbuffer1(:,1,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "Br", coordsbuffer1(:,2,i), &
        -gradientbuffer1(:,2,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,3,i), &
        -gradientbuffer1(:,3,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "C ", coordsbuffer1(:,4,i), &
        -gradientbuffer1(:,4,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,5,i), &
        -gradientbuffer1(:,5,i) * 51.422114d0
end do
close(6666)

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'target.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
write(6666,FMT="(I1)") Natoms
write(6666,FMT="(F14.9)") V * 4.3363d-2
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "H ", coords1(:,1), &
    -gradient1(:,1) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "Br", coords1(:,2), &
    -gradient1(:,2) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,3), &
    -gradient1(:,3) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "C ", coords1(:,4), &
    -gradient1(:,4) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,5), &
    -gradient1(:,5) * 51.422114d0
close(6666)

Nbins = 20
errorIMRRbinwidth = (maxErrorIMRR - minErrorIMRR) / Nbins
allocate(errorIMRRbinning(Nbins,3,Nalpha))

errorIMRRbinning = 0
open(6666,file=gridpath5//"errorIMRRtest.dat")
do
    read(6666,iostat=iostate,FMT=*) &
        j,errorID,errorIMRR,&
        R1,R2, Ninterpolation_current,&
        largestNegativeWeight
!       totalWeight,largestNegativeWeight
    if (iostate /= 0) exit

    errorBin = floor((errorIMRR-minErrorIMRR)/errorIMRRbinwidth) + 1
    if (errorBin < 1) errorBin = 1
    if (errorBin > Nbins) errorBin = Nbins

    if (errorID == 1) &
    errorIMRRbinning(errorBin,1,j) = &
        errorIMRRbinning(errorBin,1,j) + 1
    if (errorID == 2) &
    errorIMRRbinning(errorBin,2,j) = &
        errorIMRRbinning(errorBin,2,j) + 1

    errorIMRRbinning(errorBin,3,j) = &
        errorIMRRbinning(errorBin,3,j) + 1

end do
close(6666)

maxErrorIMRRoccurence = maxval(&
        errorIMRRbinning(1:Nbins,1:3,1:Nalpha))

meanR1 = meanR1 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))
meanR2 = meanR2 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))

FMTalpha = "(ES14.7"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,I6"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"errorIMRRtest_binned1.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,1,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned2.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,2,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned3.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,3,1:Nalpha)
end do
close(6666)

deallocate(errorIMRRbinning)

!getIMRRerrors_counter = getIMRRerrors_counter + 1
getIMRRerrors_counter = NC

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRR',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
!write(6666,FMT="(A,F10.7,A)") &
!        'set label 3 "Efinal = ', Efinal, '" at screen 0.82,0.86'
!write(6666,*) 'set tmargin at screen 0.8'
!write(6666,*) 'set bmargin at screen 0.2'
!write(6666,*) 'set rmargin at screen 0.8'
!write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'Nalpha = ', Nalpha
write(6666,FMT=*) 'ymax = ', maxErrorIMRRoccurence
write(6666,FMT=*) 'xmin = ', max(minErrorIMRR-errorIMRRbinwidth,0.0d0)
write(6666,FMT=*) 'xmax = ', maxErrorIMRR+errorIMRRbinwidth

write(6666,FMT="(A)") 'set multiplot layout Nalpha,1 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0.1,0'
write(6666,FMT="(A)") 'set xrange [xmin:xmax]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set style fill transparent solid 1.0'
write(6666,FMT=*) 'set boxwidth ', errorIMRRbinwidth

write(6666,FMT="(A)") 'unset xlabel'
write(6666,FMT="(A)") 'unset xtics'
do i = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1))
    end if
    write(6666,FMT="(A,F7.3,A)") &
            'set ylabel "P_{1.0e',log10(alpha_ratio),'}"'
    if (i == Nalpha) then
        write(6666,FMT="(A)") 'set xtics'
        write(6666,FMT="(A)") 'set xlabel "IMRR Error"'
    end if
    write(6666,FMT="(A,I1,A)") &
            'plot "'//gridpath5//'errorIMRRtest_binned3.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "red",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned2.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "green",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned1.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "blue"'
end do
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)


open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRr1r2',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', sqrt(maxR2*meanR2)
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

!write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set xrange [:xmax*1.1]'
write(6666,FMT="(A)") 'set logscale x'

write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "R2 (A^2)"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:5:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRw',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', 1.5d0
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'
write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'

write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "Summed Negative Weight"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:7:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

end if

deallocate(reorderedIndexes)
deallocate(Omega,outputMFN,inputMFN)
deallocate(libUgradients,libenergies)

deallocate(inputJacobians,Jacobian1,&
           inputFeatures,Feature1)

return
end subroutine getIMRRerrors5



subroutine checkIMRRerrors1(Q,PDOT,labelID)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use VARIABLES
use SIMILARITY
use interactMultipleGrids
use svd1986
implicit none

COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC
COMMON/FORCES/NATOMS

integer :: NATOMS, NC
double precision :: NTZ,NT
double precision :: ISEED0
double precision :: TIME
double precision :: T,V,H

integer :: Ncoords

double precision,dimension(3*Natoms),intent(in) :: Q, PDOT
integer,intent(in) :: labelID

integer :: Ninterpolation_pool
integer, allocatable :: reorderedIndexes(:)

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)
real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp), allocatable :: inputCLS2(:,:)
real(dp), allocatable :: libUgradients(:,:,:)
real(dp), allocatable :: libenergies(:)

real(dp),dimension(3,3) :: U
real(dp),dimension(NSIs) :: new_SIs

real(dp),dimension(3,Natoms) :: coords1,coords2
real(dp),dimension(3,Natoms) :: gradient1,gradient2
real(dp) :: min_rmsd, errorIMRR, Efinal
real(dp) :: minErrorIMRR, maxErrorIMRR
real(dp) :: errorIMRRbinwidth
integer :: maxErrorIMRRoccurence

character(400) :: FMTalpha, filetext
integer :: Nbins, errorBin
integer,allocatable :: errorIMRRbinning(:,:,:)

real(dp),dimension(3*NATOMS) :: meanCoords
real(dp) :: varCoords
real(dp) :: igDistance, igDistanceMax, varCoordsMax
integer :: maxNinterpolation_current

real(dp) :: R1, R2, maxR1, maxR2, meanR1, meanR2
real(dp) :: maxErrorIMRR_cb = 0.06d0

real(dp) :: weightThreshold = 5.0d-3
integer :: errorID, Ninterpolation_current

real(dp) :: totalWeight,largestNegativeWeight
real(dp) :: l1NegativeWeight, l2NegativeWeight
real(dp),dimension((NATOMS*(NATOMS-1))/2) :: meanCM, currentCM
real(dp) :: icmDistance, varCM
real(dp),dimension((NATOMS*(NATOMS-1))/2) :: meanDM, currentDM
real(dp) :: idmDistance, varDM

character(len=*),parameter :: consolidatedErrorIMRRfile = &
                  "consolidatedErrorIMRR.dat"
real(dp) :: minR2
real(dp),dimension(Nalpha) :: chosenErrorIMRR, &
                chosenR1IMRR, chosenR2IMRR, &
                chosenLiNW, chosenL1NW, chosenL2NW, &
                chosenIGD, chosenIGV, &
                chosenIDMD, chosenIDMV, &
                chosenICMD, chosenICMV
real(dp) :: potE2

integer*8 :: dummyINT1, dummyINT2, dummyINT3

integer, save :: getIMRRerrors_counter = 0
integer :: iostate, step, i, j, k, l

real(dp),allocatable :: targetDescriptor(:),&
        inputComboDescriptors(:,:,:),&
        inputComboErrors(:,:)
real(dp),allocatable :: PD(:,:), PDinputs(:,:), G(:,:)
integer,allocatable :: inputsOrder(:)

integer :: Ncombos_max = 200
integer :: Ncombos
logical :: exit_flag

!real(dp),dimension(Nalpha,Ncombos_max)
!real(dp),allocatable :: targetDescriptor(:),&
!        inputComboDescriptors(:,:,:)
!
!real(dp),allocatable :: PDinputs(:,:)
!integer,allocatable :: inputsOrder(:)
!real(dp),dimension(Ncombos_max,Ncombos_max) :: PD
!real(dp),dimension(Ncombos_max,2) :: G

real :: rand
integer :: randomK

logical :: output_flag
logical :: potE_flag = .false.

!output_flag = (modulo(NC,100)==1)
output_flag = .true.

Ncoords = 3 * NATOMS

coords1 = reshape(Q,(/3,Natoms/))
gradient1 = reshape(PDOT,(/3,Natoms/))

call getVarsHBrCO2(coords1,Natoms,vals,&
               Nvar,BOND_LABELLING_DATA)

!print *, "whaaat 1"
call checkState_PCM(&
        vals,coords1,&
        gradient2,&
        dummyINT1,dummyINT2,dummyINT3)
!print *, "whaaat 2"

! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (getIMRRerrors_counter > 1) return
!if (any(vals > 6.0d0)) return
!if (Ninterpolation < 25) return
if (Ninterpolation < 16) return !Ninterpolation_cap+5) return
! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ninterpolation_pool = Ninterpolation !min(Ninterpolation,Ninterpolation_cap)
!Ninterpolation = min(Ninterpolation,&
!        Ninterpolation_cap)
Ninterpolation = 15
randomK = 1+5+floor(rand()*(Ninterpolation_pool-(1+5)))
randomK = min(max(randomK,1+5),Ninterpolation_pool)

allocate(reorderedIndexes(Ninterpolation),&
         frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         inputCLS2(Ncoords+Ninterpolation,Ninterpolation),&
         libUgradients(3,NATOMS,Ninterpolation),&
         libenergies(Ninterpolation))


if (output_flag) then
allocate(targetDescriptor(3),&
         inputComboErrors(Nalpha,Ncombos_max),&
         inputComboDescriptors(3,&
         Ninterpolation,Ncombos_max))
!allocate(targetDescriptor(3),&
!         inputComboErrors(Nalpha,Ninterpolation_pool),&
!         inputComboDescriptors(3,&
!         Ninterpolation,Ninterpolation_pool))

targetDescriptor(1) = &
    sqrt(sum((coords1(1:3,3)-&
              coords1(1:3,4))**2))
targetDescriptor(2) = &
    sqrt(sum((coords1(1:3,1)-&
              coords1(1:3,2))**2))
targetDescriptor(3) = &
    sqrt(sum((coords1(1:3,2)-&
              coords1(1:3,4))**2))
!    abs(dot_product(coords1(1:3,2)-coords1(1:3,3),&
!                    coords1(1:3,4)-coords1(1:3,3))) / &
!    sqrt(sum((coords1(1:3,2)-&
!              coords1(1:3,3))**2) *&
!         sum((coords1(1:3,4)-&
!              coords1(1:3,3))**2))
end if

call setTarget(coords1)

maxNinterpolation_current = 0
maxR1 = 0.0d0
maxR2 = 0.0d0
meanR1 = 0.0d0
meanR2 = 0.0d0
minErrorIMRR = 1.0d9
maxErrorIMRR = 0.0d0
igDistanceMax = 0.0d0
varCoordsMax = 0.0d0

open(6666,file=gridpath5//"errorIMRRtest.dat")

do j = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((j-1) * 1.0d0 / (Nalpha-1))
    end if

    errorID = 3
    Ninterpolation_current = Ninterpolation
    k = Ninterpolation_current+1
    Ncombos = 1
    do i = 1, Ninterpolation_current
        reorderedIndexes(i) = i
    end do

    chosenErrorIMRR(j) = 1.0d9
    minR2 = 1.0d9

do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part only looks at the inputs

do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)

    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    libUgradients(:,:,i) = matmul(&
            U,gradientbuffer1(:,:,step))
    libenergies(i) = &
            potEbuffer1(step)

    inputCLS2(1:Ncoords,i) = reshape(&
            coords2 - coords1,(/Ncoords/))
    inputCLS2(Ncoords+i,:) = 0.0d0
    inputCLS2(Ncoords+i,i) = alpha_ratio * &
            new_SIs(Nsort)**2

    if (output_flag) then
    if (j == 1) then
    inputComboDescriptors(1,i,Ncombos) = & !k-1) = &
        sqrt(sum((coordsbuffer1(1:3,3,step)-&
                  coordsbuffer1(1:3,4,step))**2)) !-&
!       targetDescriptor(1)
    inputComboDescriptors(2,i,Ncombos) = & !k-1) = &
        sqrt(sum((coordsbuffer1(1:3,1,step)-&
                  coordsbuffer1(1:3,2,step))**2)) !-&
!       targetDescriptor(2)
    inputComboDescriptors(3,i,Ncombos) = & !k-1) = &
        sqrt(sum((coordsbuffer1(1:3,2,step)-&
                  coordsbuffer1(1:3,4,step))**2)) !-&
!       targetDescriptor(3)

!       abs(dot_product(coordsbuffer1(1:3,2,step)-&
!                       coordsbuffer1(1:3,3,step),&
!                       coordsbuffer1(1:3,4,step)-&
!                       coordsbuffer1(1:3,3,step))) / &
!       sqrt(sum((coordsbuffer1(1:3,2,step)-&
!                 coordsbuffer1(1:3,3,step))**2) *&
!            sum((coordsbuffer1(1:3,4,step)-&
!                 coordsbuffer1(1:3,3,step))**2)) !-&
!       targetDescriptor(3)
    end if
    end if

end do

meanCoords = 0.0d0
do i = 1, Ninterpolation_current
    meanCoords = meanCoords + &
                 inputCLS2(1:Ncoords,i)
end do
meanCoords = meanCoords / Ninterpolation_current
igDistance = sqrt(sum(meanCoords**2)/Natoms)

varCoords = 0.0d0
do i = 1, Ninterpolation_current
    varCoords = varCoords + &
                sum((inputCLS2(1:Ncoords,i)-meanCoords)**2)
end do
varCoords = sqrt(varCoords/Ninterpolation_current)

igDistanceMax = max(igDistance,igDistanceMax)
varCoordsMax = max(varCoords,varCoordsMax)

meanCM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector_test(&
            coordsbuffer1(:,:,step),currentCM)
    meanCM = meanCM + currentCM
end do
call getCoulombMatrixVector_test(coords1,currentCM)
meanCM = (meanCM / Ninterpolation_current) - currentCM
icmDistance = sqrt(sum(meanCM**2)/(Natoms*(Natoms-1)*0.5))

varCM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector_test(&
            coordsbuffer1(:,:,step),currentCM)
    varCM = varCM + sum((currentCM-meanCM)**2)
end do
varCM = sqrt(varCM/Ninterpolation_current)

meanCM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector(&
            coordsbuffer1(:,:,step),currentCM)
    meanCM = meanCM + currentCM
end do
call getCoulombMatrixVector(coords1,currentCM)
meanCM = (meanCM / Ninterpolation_current) - currentCM
idmDistance = sqrt(sum(meanCM**2)/(Natoms*(Natoms-1)*0.5))

varDM = 0.0d0
do i = 1, Ninterpolation_current
    step = reorderedIndexes(i)
    call getCoulombMatrixVector(&
            coordsbuffer1(:,:,step),currentCM)
    varDM = varDM + sum((currentCM-meanCM)**2)
end do
varDM = sqrt(varDM/Ninterpolation_current)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outputCLS = 0.0d0
restraints = 1.0d0
restraint_values = 1.0d0

!call CLS2(inputCLS2(1:Ncoords+Ninterpolation_current,&
!                    1:Ninterpolation_current),&
!          Ncoords+Ninterpolation_current,&
!          Ninterpolation_current,&
!          restraints(1:1,1:Ninterpolation_current),&
!          1,restraint_values,&
!          outputCLS(1:Ncoords+Ninterpolation_current),&
!          frame_weights(1:Ninterpolation_current))
call CLSwK(inputCLS2(1:Ncoords,&
                    1:Ninterpolation_current),&
          Ncoords,&
          Ninterpolation_current,&
          restraints(1:1,1:Ninterpolation_current),&
          1,restraint_values,&
          outputCLS(1:Ncoords),&
          frame_weights(1:Ninterpolation_current))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This part looks at the interpolation

if (all(frame_weights(1:Ninterpolation_current) == 0.0d0)) then
!   if (Ninterpolation_current < Ninterpolation) &
!       Ninterpolation_current = &
!           Ninterpolation_current + 1

!   reorderedIndexes(Ninterpolation_current) = k
!   k = k + 1

!   if (k > Ninterpolation_pool) exit
!   cycle

    frame_weights(1) = 1.0d0
end if

gradient2 = 0.0d0
potE2 = 0.0d0
R2 = 0.0d0
do i = 1, Ninterpolation_current
    gradient2 = gradient2 + frame_weights(i)*&
                libUgradients(:,:,i)
    potE2 = potE2 + frame_weights(i)*&
                libenergies(i)
    R2 = R2 + abs(frame_weights(i)) * &
            inputCLS2(Ncoords+i,i)
end do
R2 = R2 / alpha_ratio
R1 = sqrt(sum(matmul(&
    inputCLS2(1:Ncoords,1:Ninterpolation_current),&
    reshape(frame_weights(1:Ninterpolation_current),&
            (/ Ninterpolation_current, 1 /)))**2))

totalWeight = sum(frame_weights(1:Ninterpolation_current))
l1NegativeWeight = abs(sum(frame_weights(1:Ninterpolation_current),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0))
l2NegativeWeight = sqrt(sum((frame_weights(1:Ninterpolation_current)**2),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0))
largestNegativeWeight = maxval(abs(frame_weights(1:Ninterpolation_current)),& 
                       mask=frame_weights(1:Ninterpolation_current)<0.0d0)


if (potE_flag) then
    errorIMRR = abs(V - potE2)
else
    errorIMRR = sqrt(sum((gradient1-gradient2)**2)/Natoms)
end if

if (output_flag) &
    inputComboErrors(j,Ncombos) = errorIMRR
!   inputComboErrors(j,k-1) = errorIMRR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (k == Ninterpolation_pool) errorID = 1

open(6667,file=gridpath5//"all"//consolidatedErrorIMRRfile,&
          position="append")
write(6667,FMT="(I8,1x,I3,3(1x,ES14.7))") &
    NC, labelID, errorIMRR, R1, R2
close(6667)

!write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2)") &
if (output_flag) &
write(6666,FMT="(I2,1x,I3,5(1x,ES14.7),1x,I2,1x,ES14.7)") &
                  j,errorID,&
                  errorIMRR,igDistance,varCoords,&
                  R1, R2, Ninterpolation_current,&
                  largestNegativeWeight
!                 totalWeight,largestNegativeWeight

maxNinterpolation_current = max(&
             Ninterpolation_current,&
          maxNinterpolation_current)
maxR1 = max(R1,maxR1)
maxR2 = max(R2,maxR2)
meanR1 = meanR1 + R1
meanR2 = meanR2 + R2
maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!
! Selection Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(ssID)
! Smallest R2
case(1)
    if (R2 < minR2) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
        chosenIGD(j) = igDistance
        chosenIGV(j) = varCoords
        chosenIDMD(j) = idmDistance
        chosenIDMV(j) = varDM
        chosenICMD(j) = icmDistance
        chosenICMV(j) = varCM
        chosenL1NW(j) = l1NegativeWeight
        chosenL2NW(j) = l2NegativeWeight
        chosenLiNW(j) = largestNegativeWeight
    end if
! Minimal Error
case(2)
    if (errorIMRR < chosenErrorIMRR(j)) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
        chosenIGD(j) = igDistance
        chosenIGV(j) = varCoords
        chosenIDMD(j) = idmDistance
        chosenIDMV(j) = varDM
        chosenICMD(j) = icmDistance
        chosenICMV(j) = varCM
        chosenL1NW(j) = l1NegativeWeight
        chosenL2NW(j) = l2NegativeWeight
        chosenLiNW(j) = largestNegativeWeight
    end if
! Random
case(3)
    if (k == randomK) then
        chosenErrorIMRR(j) = errorIMRR
        chosenR1IMRR(j) = R1
        chosenR2IMRR(j) = R2
        chosenIGD(j) = igDistance
        chosenIGV(j) = varCoords
        chosenIDMD(j) = idmDistance
        chosenIDMV(j) = varDM
        chosenICMD(j) = icmDistance
        chosenICMV(j) = varCM
        chosenL1NW(j) = l1NegativeWeight
        chosenL2NW(j) = l2NegativeWeight
        chosenLiNW(j) = largestNegativeWeight
    end if
end select

minR2 = min(minR2,R2)
minErrorIMRR = min(minErrorIMRR,errorIMRR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!
! Replacement Strategy
!!!!!!!!!!!!!!!!!!!!!!!

select case(rsID)
case(1)
    if (any(abs(frame_weights(1:Ninterpolation_current)) < &
        weightThreshold)) then
        do i = Ninterpolation_current, 1, -1
            if (abs(frame_weights(i)) < weightThreshold) then
                if (i < Ninterpolation_current) &
                    reorderedIndexes(i) = &
                        reorderedIndexes(i+1)
                Ninterpolation_current = &
                    Ninterpolation_current - 1
            end if
        end do
        errorID = 2
    else
        if (Ninterpolation_current < Ninterpolation) then
            Ninterpolation_current = &
                Ninterpolation_current + 1
            reorderedIndexes(Ninterpolation_current) = k
        else
            reorderedIndexes(minloc(abs(&
                  frame_weights(1:Ninterpolation_current)))) = k
        end if
    
        if (k > Ninterpolation_pool) exit
        k = k + 1
    
        errorID = 3
    end if
case(2)
    if (Ninterpolation_current < Ninterpolation) then
        Ninterpolation_current = &
            Ninterpolation_current + 1
        reorderedIndexes(Ninterpolation_current) = k
        errorID = 2
    else
        if (minval(frame_weights(1:Ninterpolation_current))<0.0d0) then
            errorID = 2
        else
            errorID = 3
        end if
        reorderedIndexes(minloc(&
              frame_weights(1:Ninterpolation_current))) = k
    end if

    if (k > Ninterpolation_pool) exit
    k = k + 1

! Just go through all the combinations until
!   (a) we run out of combinations or
!   (b) we have reached the maximum
!       number of combinations
case(3)
    call nextCombination(Ninterpolation_current,&
                         Ninterpolation_pool,&
                         reorderedIndexes(1:Ninterpolation_current),&
                         exit_flag)
    if (Ncombos == Ncombos_max) exit_flag = .true.
    if (exit_flag) then
        if (Ncombos == 1) output_flag = .false.
        exit
    else
        Ncombos = Ncombos + 1
    end if
end select

end do

end do

close(6666)

FMTalpha = "(I8,1x,I3"
do i = 1, 3*Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,ES14.7"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//consolidatedErrorIMRRfile,&
          position="append")
write(6666,FMT=trim(adjustl(FMTalpha))) &
    NC, labelID, chosenErrorIMRR, chosenR1IMRR, chosenR2IMRR
close(6666)

FMTalpha = "(I8,1x,I3"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",12(1x,ES14.7)"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"other"//consolidatedErrorIMRRfile,&
          position="append")
write(6666,FMT=trim(adjustl(FMTalpha))) &
    NC, labelID, chosenErrorIMRR, chosenR1IMRR, chosenR2IMRR, &
        chosenL1NW, chosenL2NW, chosenLiNW, &
        chosenIGD, chosenIGV, &
        chosenIDMD, chosenIDMV, &
        chosenICMD, chosenICMV
close(6666)







if (output_flag) then

! For the sGDML comparison

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'inputs.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
do i = 1, Ninterpolation_pool
    write(6666,FMT="(I1)") Natoms
    write(6666,FMT="(F14.9)") potEbuffer1(i) * 4.3363d-2
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "H ", coordsbuffer1(:,1,i), &
        -gradientbuffer1(:,1,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "Br", coordsbuffer1(:,2,i), &
        -gradientbuffer1(:,2,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,3,i), &
        -gradientbuffer1(:,3,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "C ", coordsbuffer1(:,4,i), &
        -gradientbuffer1(:,4,i) * 51.422114d0
    write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
        "O ", coordsbuffer1(:,5,i), &
        -gradientbuffer1(:,5,i) * 51.422114d0
end do
close(6666)

write(filetext,FMT="(I0.8,'_',I0.3,A)") &
        getIMRRerrors_counter,labelID,'target.xyz'
open(6666,file=gridpath4//trim(adjustl(filetext)))
write(6666,FMT="(I1)") Natoms
write(6666,FMT="(F14.9)") V * 4.3363d-2
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "H ", coords1(:,1), &
    -gradient1(:,1) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "Br", coords1(:,2), &
    -gradient1(:,2) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,3), &
    -gradient1(:,3) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "C ", coords1(:,4), &
    -gradient1(:,4) * 51.422114d0
write(6666,FMT="(A2,3(1x,F10.6),3(1x,F10.6))") &
    "O ", coords1(:,5), &
    -gradient1(:,5) * 51.422114d0
close(6666)


! To make it more readable
maxErrorIMRR = min(maxErrorIMRR, 0.06d0)

Nbins = 20
errorIMRRbinwidth = (maxErrorIMRR - minErrorIMRR) / Nbins
allocate(errorIMRRbinning(Nbins,3,Nalpha))

errorIMRRbinning = 0
open(6666,file=gridpath5//"errorIMRRtest.dat")
do
    read(6666,iostat=iostate,FMT=*) &
        j,errorID,errorIMRR,igDistance,&
        varCoords,R1,R2, Ninterpolation_current,&
        largestNegativeWeight
!       totalWeight,largestNegativeWeight
    if (iostate /= 0) exit

    errorBin = floor((errorIMRR-minErrorIMRR)/errorIMRRbinwidth) + 1
    if (errorBin < 1) errorBin = 1
    if (errorBin > Nbins) errorBin = Nbins

    if (errorID == 1) &
    errorIMRRbinning(errorBin,1,j) = &
        errorIMRRbinning(errorBin,1,j) + 1
    if (errorID == 2) &
    errorIMRRbinning(errorBin,2,j) = &
        errorIMRRbinning(errorBin,2,j) + 1

    errorIMRRbinning(errorBin,3,j) = &
        errorIMRRbinning(errorBin,3,j) + 1

end do
close(6666)

maxErrorIMRRoccurence = maxval(&
        errorIMRRbinning(1:Nbins,1:3,1:Nalpha))

meanR1 = meanR1 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))
meanR2 = meanR2 / sum(errorIMRRbinning(1:Nbins,3,1:Nalpha))

FMTalpha = "(ES14.7"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,I6"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

open(6666,file=gridpath5//"errorIMRRtest_binned1.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,1,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned2.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,2,1:Nalpha)
end do
close(6666)

open(6666,file=gridpath5//"errorIMRRtest_binned3.dat")
do i = 1, Nbins
    write(6666,FMT=trim(adjustl(FMTalpha))) &
            minErrorIMRR + &
            (i-0.5)*errorIMRRbinwidth, &
            errorIMRRbinning(i,3,1:Nalpha)
end do
close(6666)

deallocate(errorIMRRbinning)

!getIMRRerrors_counter = getIMRRerrors_counter + 1
getIMRRerrors_counter = NC

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRR',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
!write(6666,FMT="(A,F10.7,A)") &
!        'set label 3 "Efinal = ', Efinal, '" at screen 0.82,0.86'
!write(6666,*) 'set tmargin at screen 0.8'
!write(6666,*) 'set bmargin at screen 0.2'
!write(6666,*) 'set rmargin at screen 0.8'
!write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'Nalpha = ', Nalpha
write(6666,FMT=*) 'ymax = ', maxErrorIMRRoccurence
write(6666,FMT=*) 'xmin = ', max(minErrorIMRR-errorIMRRbinwidth,0.0d0)
write(6666,FMT=*) 'xmax = ', maxErrorIMRR+errorIMRRbinwidth

write(6666,FMT="(A)") 'set multiplot layout Nalpha,1 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0.1,0'
write(6666,FMT="(A)") 'set xrange [xmin:xmax]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set style fill transparent solid 1.0'
write(6666,FMT=*) 'set boxwidth ', errorIMRRbinwidth

write(6666,FMT="(A)") 'unset xlabel'
write(6666,FMT="(A)") 'unset xtics'
do i = 1, Nalpha
    if (logarithmic_alpha_flag) then
        alpha_ratio = 10.0d0**(&
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1)))
    else
        alpha_ratio = &
            alpha_start + (alpha_end-alpha_start) *&
                          ((i-1) * 1.0d0 / (Nalpha-1))
    end if
    write(6666,FMT="(A,F7.3,A)") &
            'set ylabel "P_{1.0e',log10(alpha_ratio),'}"'
    if (i == Nalpha) then
        write(6666,FMT="(A)") 'set xtics'
        write(6666,FMT="(A)") 'set xlabel "IMRR Error"'
    end if
    write(6666,FMT="(A)") 'unset arrow'
    write(6666,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', chosenErrorIMRR(i),&
        ',graph 0 to ', chosenErrorIMRR(i),&
        ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
    write(6666,FMT="(A)") 'unset label 3'
    write(6666,FMT="(A,E16.6,A)")&
        'set label 3 "', chosenErrorIMRR(i),&
        ' (E_h/a_0)" at graph 0.8, graph 0.9 '//&
        ' font ",12" front'
    write(6666,FMT="(A,I1,A)") &
            'plot "'//gridpath5//'errorIMRRtest_binned3.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "red",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned2.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "green",\'
    write(6666,FMT="(A,I1,A)") &
            '     "'//gridpath5//'errorIMRRtest_binned1.dat" '//&
            'u 1:',i+1,' w boxes lc rgb "blue"'
end do
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigd',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', igDistanceMax
write(6666,FMT=*) 'ymax = ', maxErrorIMRR

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center Distance (A)"'
write(6666,FMT="(A)") 'set ylabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:3:1 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigv',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', varCoordsMax
write(6666,FMT=*) 'ymax = ', maxErrorIMRR

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center Variance (A)"'
write(6666,FMT="(A)") 'set ylabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 5:3:1 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRigdN',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', igDistanceMax
write(6666,FMT=*) 'ymax = ', maxNinterpolation_current
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "Input Geometric Center distance (A)"'
write(6666,FMT="(A)") 'set ylabel "Ninterpolation"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 4:8:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)


open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRr1r2',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', sqrt(maxR2*meanR2)
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'

!write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set xrange [:xmax*1.1]'
write(6666,FMT="(A)") 'set logscale x'

write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "R2 (A^2)"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 6:7:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,2000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorIMRRw',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.94'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.9'
write(6666,*) 'set tmargin at screen 0.8'
write(6666,*) 'set bmargin at screen 0.2'
write(6666,*) 'set rmargin at screen 0.8'
write(6666,*) 'set lmargin at screen 0.2'
write(6666,FMT=*) 'cbmax = ', maxErrorIMRR_cb
write(6666,FMT=*) 'xmax = ', sqrt(maxR1*meanR1)
write(6666,FMT=*) 'ymax = ', 1.5d0
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'
write(6666,FMT="(A)") 'set xrange [0:xmax*1.1]'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'

write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "R1 (A)"'
write(6666,FMT="(A)") 'set ylabel "Summed Negative Weight"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (E_h/a_0)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//'errorIMRRtest.dat" '//&
        'u 6:9:3 w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This part looks at the errors
! of individual frames

maxErrorIMRR = 0.0d0
open(6666,file=gridpath5//temporaryfile1)
do step = 1, Ninterpolation_pool
    call getSIs(coordsbuffer1(:,:,step),&
            coords2,U,new_SIs)
    gradient2 = matmul(&
            U,gradientbuffer1(:,:,step))

    errorIMRR = sqrt(sum((gradient2 - gradient1)**2)/Natoms)
    maxErrorIMRR = max(maxErrorIMRR,errorIMRR)

    write(6666,FMT="(3(ES14.7,1x),ES14.7)") &
        new_SIs(1:3), errorIMRR
end do
close(6666)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 2000,1000"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'errorInputs',&
        getIMRRerrors_counter,labelID,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.82,0.96'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Nfinal = ', Ninterpolation_current,&
                             Ninterpolation_pool, '" at screen 0.82,0.92'

write(6666,FMT="(A)") 'set multiplot layout 1,3 '//&
        'margins 0.2,0.8,0.1,0.9 spacing 0,0.1'
write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT=*) 'ymax = ', maxErrorIMRR
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set xrange [0:]'
write(6666,FMT="(A,F7.3,A)") &
        'set ylabel "Difference in Energy Gradient (E_h/a_0)"'
write(6666,FMT="(A)") 'set xlabel "RMSD (A)"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//temporaryfile1//'" '//&
        'u 1:4 w p ps 2 pt 7'
write(6666,FMT="(A)") 'unset ytics'
write(6666,FMT="(A)") 'unset ylabel'
write(6666,FMT="(A)") 'set yrange [0:ymax*1.1]'
write(6666,FMT="(A)") 'set xlabel "CMD w/ Charges"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//temporaryfile1//'" '//&
        'u 2:4 w p ps 2 pt 7'
write(6666,FMT="(A)") 'set xlabel "CMD w/o Charges"'
write(6666,FMT="(A)") &
        'plot "'//gridpath5//temporaryfile1//'" '//&
        'u 3:4 w p ps 2 pt 7'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we will do the pairwise distance to
! graph creation if there are enough
! combinations

if (Ncombos > 1) then

allocate(PD(Ncombos,Ncombos),&
         PDinputs(Ninterpolation,Ninterpolation),&
         G(Ncombos,2),inputsOrder(Ninterpolation))
!allocate(PD(Ninterpolation_pool,Ninterpolation_pool),&
!         PDinputs(Ninterpolation,Ninterpolation),&
!         G(1+Ninterpolation_pool-Ninterpolation,2),&
!         inputsOrder(Ninterpolation))

do i = Ninterpolation, Ncombos !Ninterpolation_pool
    PD(i,i) = 0.0d0
    do j = i+1, Ncombos !Ninterpolation_pool

        ! The RMSD suggestion is not really
        ! doable because there will be too
        ! many permutations

        do k = 1, Ninterpolation
        do l = 1, Ninterpolation
            PDinputs(k,l) = sum((&
                inputComboDescriptors(1:3,k,i) - &
                inputComboDescriptors(1:3,l,j))**2)
        end do
        end do

        call munkres(PDinputs(1:Ninterpolation,&
                              1:Ninterpolation),&
                Ninterpolation,inputsOrder(1:Ninterpolation))

        PD(i,j) = 0.0d0
        do k = 1, Ninterpolation
            PD(i,j) = PD(i,j) + PDinputs(inputsOrder(k),k)
        end do
        PD(i,j) = sqrt(PD(i,j))
        PD(j,i) = PD(i,j)
    end do
end do

!k = 1 + Ninterpolation_pool - Ninterpolation
!call PD2G(PD(Ninterpolation:Ninterpolation_pool,&
!             Ninterpolation:Ninterpolation_pool),&
!          k,G(1:k,1:2))
!!!call PD2G(PD(1:Ncombos,1:Ncombos),&
!!!          Ncombos,G(1:Ncombos,1:2))
!         1+Ninterpolation_pool-Ninterpolation,G(1:k,1:2))

minErrorIMRR = 0.0d0
maxErrorIMRR = 0.0d0
do i = 1, Ncombos !k
do j = i+1, Ncombos !k
    minErrorIMRR = minErrorIMRR + &
        (sqrt(sum((G(i,1:2)-G(j,1:2))**2)) &
                - PD(i,j))**2
!                - PD(Ninterpolation+i-1,&
!                     Ninterpolation+j-1))**2
    maxErrorIMRR = maxErrorIMRR + PD(i,j)**2
!                 (PD(Ninterpolation+i-1,&
!                     Ninterpolation+j-1))**2
end do
end do

write(FMTalpha,FMT="(A,I0.8,'_',I0.3,A)") &
        gridpath5//'PD',&
        getIMRRerrors_counter,labelID,'.dat'
open(6666,file=trim(adjustl(FMTalpha)))

!write(FMTalpha,FMT="(I3)") &
!        1 + Ninterpolation_pool - Ninterpolation
write(FMTalpha,FMT="(I3)") Ncombos
FMTalpha = "(I3,1x,"//trim(adjustl(FMTalpha))//&
           "(1x,ES14.7))"

!do i = Ninterpolation, Ninterpolation_pool 
do i = 1, Ncombos
    write(6666,FMT=trim(adjustl(FMTalpha))) i, &
            PD(i,1:Ncombos)
!            PD(i,Ninterpolation:Ninterpolation_pool)
end do
close(6666)



write(FMTalpha,FMT="(A,I0.8,'_',I0.3,A)") &
        gridpath5//'graph',&
        getIMRRerrors_counter,labelID,'.dat'
open(6666,file=trim(adjustl(FMTalpha)))

FMTalpha = "(ES14.7,1x,ES14.7,"
do i = 1, Nalpha
    FMTalpha = trim(adjustl(FMTalpha))//",1x,ES14.7"
end do
FMTalpha = trim(adjustl(FMTalpha))//")"

do i = 1, Ncombos !k
    write(6666,FMT=FMTalpha) &
            G(i,1:2), inputComboErrors(1:Nalpha,i) * 51.422114d0
!            G(i,1:2), inputComboErrors(1:Nalpha,Ninterpolation+i-1) * 51.422114d0
end do
close(6666)

meanR1 = minval(G(1:Ncombos,1))
maxR1 = maxval(G(1:Ncombos,1))
meanR2 = minval(G(1:Ncombos,2))
maxR2 = maxval(G(1:Ncombos,2))
!meanR1 = minval(G(1:k,1))
!maxR1 = maxval(G(1:k,1))
!meanR2 = minval(G(1:k,2))
!maxR2 = maxval(G(1:k,2))

write(FMTalpha,FMT="(A,I0.8,'_',I0.3,A)") &
        gridpath5//'graph',&
        getIMRRerrors_counter,labelID,'.dat'

do j = 1, Nalpha
open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 1200,1200"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'graph',&
        getIMRRerrors_counter,j,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.02,0.97'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Ninterpolation = ', Ninterpolation,&
                 Ninterpolation_pool, '" at screen 0.02,0.93'
write(6666,FMT="(A,ES9.2,'/',ES9.2,'=',F6.2,A)") &
        'set label 3 "Stress = ', sqrt(minErrorIMRR), &
        sqrt(maxErrorIMRR), sqrt(minErrorIMRR*10000.0d0/maxErrorIMRR),&
        '%" at screen 0.02,0.89'
write(6666,FMT=*) 'i = ', j + 2
write(6666,FMT="(A)") 'set tmargin at screen 0.8'
write(6666,FMT="(A)") 'set bmargin at screen 0.2'
write(6666,FMT="(A)") 'set rmargin at screen 0.8'
write(6666,FMT="(A)") 'set lmargin at screen 0.2'
write(6666,FMT=*) 'cbmax = ', 0.15d0 !maxErrorIMRR_cb
write(6666,FMT=*) 'xmin = ', meanR1
write(6666,FMT=*) 'ymin = ', meanR2
write(6666,FMT=*) 'xmax = ', maxR1
write(6666,FMT=*) 'ymax = ', maxR2
write(6666,FMT=*) 'dx = ', (maxR1-meanR1) * 0.05
write(6666,FMT=*) 'dy = ', (maxR2-meanR2) * 0.05
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'
write(6666,FMT="(A)") 'set xrange [xmin-dx:xmax+dx]'
write(6666,FMT="(A)") 'set yrange [ymin-dy:ymax+dy]'

write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "PC1"'
write(6666,FMT="(A)") 'set ylabel "PC2"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (eV/A)"'
write(6666,FMT="(A)") &
        'plot "'//trim(adjustl(FMTalpha))//&
        '" u 1:2:(column(i)) w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)

open(6666,file=gridpath5//gnuplotfile)
write(6666,FMT="(A)") "set term pngcairo enhanced size 1200,1200"
write(6666,FMT="(A)") 'set encoding utf8'
write(6666,FMT="(A,I0.8,'_',I0.3,A)") &
        'set output "'//gridpath4//'graph',&
        getIMRRerrors_counter,100+j,'.png"'
write(6666,FMT="(A,F7.4,',',F7.4,A)") &
        'set label 1 "CV = ', vals, '" at screen 0.02,0.97'
write(6666,FMT="(A,I0.3,'/',I0.3,A)") &
        'set label 2 "Ninterpolation = ', Ninterpolation,&
                 Ninterpolation_pool, '" at screen 0.02,0.93'
write(6666,FMT="(A,ES9.2,'/',ES9.2,'=',F6.2,A)") &
        'set label 3 "Stress = ', sqrt(minErrorIMRR), &
        sqrt(maxErrorIMRR), sqrt(minErrorIMRR*10000.0d0/maxErrorIMRR),&
        '%" at screen 0.02,0.89'
write(6666,FMT=*) 'i = ', j + 2
write(6666,FMT="(A)") 'set tmargin at screen 0.8'
write(6666,FMT="(A)") 'set bmargin at screen 0.2'
write(6666,FMT="(A)") 'set rmargin at screen 0.8'
write(6666,FMT="(A)") 'set lmargin at screen 0.2'
write(6666,FMT=*) 'cbmax = ', 0.55d0 !maxErrorIMRR_cb
write(6666,FMT=*) 'xmin = ', meanR1
write(6666,FMT=*) 'ymin = ', meanR2
write(6666,FMT=*) 'xmax = ', maxR1
write(6666,FMT=*) 'ymax = ', maxR2
write(6666,FMT=*) 'dx = ', (maxR1-meanR1) * 0.05
write(6666,FMT=*) 'dy = ', (maxR2-meanR2) * 0.05
write(6666,FMT="(A)") 'set cbrange [0:cbmax]'
write(6666,FMT="(A)") 'set xrange [xmin-dx:xmax+dx]'
write(6666,FMT="(A)") 'set yrange [ymin-dy:ymax+dy]'

write(6666,FMT="(A)") 'set border lw 4'
write(6666,FMT="(A)") 'unset key'
write(6666,FMT="(A)") 'set xlabel "PC1"'
write(6666,FMT="(A)") 'set ylabel "PC2"'
write(6666,FMT="(A)") 'set cblabel "IMRR Error (eV/A)"'
write(6666,FMT="(A)") &
        'plot "'//trim(adjustl(FMTalpha))//&
        '" u 1:2:(column(i)) w p pt 7 ps 2 lt palette'
close(6666)

call system("gnuplot < "//gridpath5//gnuplotfile)
end do

deallocate(PD,PDinputs,G,inputsOrder)

end if

deallocate(targetDescriptor,&
           inputComboErrors,&
           inputComboDescriptors)

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(reorderedIndexes,&
           frame_weights,&
           outputCLS,restraints,&
           restraint_values,&
           inputCLS2,libUgradients,libenergies)

return
end subroutine checkIMRRerrors1

end module analyzeRMSDThresholdwithMultipleGrids
