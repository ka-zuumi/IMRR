module analyzeHeatMapswithMultipleGrids
implicit none



contains

subroutine analyzeTopLevelHeatMaps(PNGname)
use DOUBLE
use PARAMETERS
use ANALYSIS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
!character(gridpath_length+expfolder_length) :: gridpath4
!character(gridpath_length+expfolder_length+5) :: gridpath5
character(*),intent(in) :: PNGname

real,dimension(Nvar) :: var_minvar = (/ 1.0, 1.5 /)
integer,dimension(Nvar) :: var_index
integer,dimension(Nvar) :: Nbin_max, Nbin_offset

integer :: population, population_max
integer :: iostate
logical :: exists
real,dimension(Nvar) :: vals
real(dp) :: potE
real(dp),dimension(3,Natoms) :: coords, gradient
character(50) :: var_filename

integer,allocatable :: POPheatmap(:,:)

integer :: i, j, n

!gridpath4 = gridpath0//expfolder
!gridpath5 = gridpath4//intermediatefolder

population_max = 20

do Ngrid = 1, Ngrid_max
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            Ngrid
    gridpath2 = gridpath0//Ngrid_text//"/grid/"
    
    Nbin_max = (var_maxvar - var_minvar) / var_spacing
    Nbin_offset = var_minvar / var_spacing
    
    allocate(POPheatmap(Nbin_max(1),Nbin_max(2)))
    POPheatmap = 0
    population_max = 0
    
    open(filechannel2,file=gridpath5//temporaryfile1)
    do i = 1, Nbin_max(1)
    do j = 1, Nbin_max(2)
    
        var_index = Nbin_offset + (/i, j/)
    
    !   write(var_filename,FMT=var_multipleFMT&
    !     (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
    !       var_index * multiplier(:,Norder+1)
        write(var_filename,FMT=var_multipleFMT&
          (1+0*multipleFMT_length:(0+1)*multipleFMT_length) )&
            var_index * multiplier(:,0+1)

        inquire(file=gridpath2//trim(var_filename),exist=exists)
print *,"checking "//trim(var_filename)
        population = 0

        if (exists) then
            if (unreadable_flag) then
            open(filechannel1,file=gridpath2//trim(var_filename),&
                    form="unformatted")
            do
                read(filechannel1,iostat=iostate) vals
print *,vals
                if (iostate /= 0) exit
                read(filechannel1) n
print *,n
                read(filechannel1) potE
print *,potE
                read(filechannel1) coords
print *,coords
                read(filechannel1) gradient
print *,gradient
            
                population = population + 1
            end do
            else
            open(filechannel1,file=gridpath2//trim(var_filename))
            do
                read(filechannel1,FMT=FMT1,advance="no",iostat=iostate) vals
print *,vals
                if (iostate /= 0) exit
                read(filechannel1,FMT="(I5)",advance="no") n
print *,n
                read(filechannel1,FMT="(ES10.2)",advance="no") potE
print *,potE
                read(filechannel1,FMT=FMT3,advance="no") coords
print *,coords
                read(filechannel1,FMT=FMT3) gradient
print *,gradient

                population = population + 1
            end do
            end if
            close(filechannel1)
        end if
    
        POPheatmap(i,j) = population
        population_max = max(population,population_max)
!       write(filechannel2,FMT="(F5.2,1x,F5.2,1x,I8)")&
        write(filechannel2,FMT="("//var_singleFMT(1:singleFMT_length)//&
                            ",1x,"//var_singleFMT(1:singleFMT_length)//",1x,I8)")&
                (var_index - (/ 0.5, 0.5 /))*multiplier(:,0+1), population
    end do
    
    write(filechannel2,FMT=*) ""
    end do
    
    close(filechannel2)
    deallocate(POPheatmap)
    
    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,FMT="(A)") 'set terminal pngcairo size 2400,1800'
    write(gnuplotchannel,FMT="(A)") 'set output "'//gridpath4//PNGname//'_'//Ngrid_text//'.png"'
    write(gnuplotchannel,FMT="(A)") 'unset key'
    !write(gnuplotchannel,FMT="(A)") 'set palette defined ( 0 0 1 0, 0.3333 0 0 1, 0.6667 1 0 0,\'
    !write(gnuplotchannel,FMT="(A)") '     1 1 0.6471 0 )'
    write(gnuplotchannel,FMT=*) 'cmax = ', population_max
    write(gnuplotchannel,FMT=*) 'xmin = ', var_minvar(1)
    write(gnuplotchannel,FMT=*) 'xmax = ', var_maxvar(1)
    write(gnuplotchannel,FMT=*) 'ymin = ', var_minvar(2)
    write(gnuplotchannel,FMT=*) 'ymax = ', var_maxvar(2)
    write(gnuplotchannel,FMT="(A)") 'set palette defined ('//&
             '0 "white",'//&
             '0.5 "white",'//&
             '0.5 "blue",'//&
             'cmax "red")'
    write(gnuplotchannel,FMT="(A)") 'set cbrange [0:cmax]'
    write(gnuplotchannel,FMT="(A)") 'set cblabel "Number of Frames in the Cell" font ",18" offset 1,0'
    write(gnuplotchannel,FMT="(A)") 'set title "Configurational Heatmap of an HBr - CO_2 System" font ",32" offset 0,3'
!   write(gnuplotchannel,FMT="(A)") 'set title "Configurational Heatmap of an H_2 - H_2 System" font ",32" offset 0,3'
    write(gnuplotchannel,FMT="(A)") 'set xlabel "Var1 (A)" font ",28" offset 0,-2'
    write(gnuplotchannel,FMT="(A)") 'set xtics 1 font ",24"'
    write(gnuplotchannel,FMT="(A)") 'set xrange [xmin:xmax]'
    write(gnuplotchannel,FMT="(A)") 'set ylabel "Var2 (A)" font ",28" offset -5,0'
    write(gnuplotchannel,FMT="(A)") 'set ytics 1 font ",24"'
    write(gnuplotchannel,FMT="(A)") 'set yrange [ymin:ymax]'
    write(gnuplotchannel,FMT="(A)") 'set cbtics'
    write(gnuplotchannel,FMT="(A)") 'set view map'
    write(gnuplotchannel,FMT="(A)") 'set pm3d interpolate 1,1'
    write(gnuplotchannel,FMT="(A)") 'splot "'//gridpath5//temporaryfile1//'" u 1:2:3 w pm3d'
    close(gnuplotchannel)
    
    call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)
    
end do

end subroutine analyzeTopLevelHeatMaps

!subroutine analyzeHeatMaps1()
!use PARAMETERS
!use mapCellData
!use ANALYSIS
!implicit none
!
!!FORMATTING OF JPG FILES
!character(5) :: variable_length_text
!character(Ngrid_text_length) :: Ngrid_text
!character(Ngrid_text_length+1) :: folder_text
!character(6) :: Ntraj_text
!
!!I/O HANDLING
!integer :: iostate
!
!!INTEGER INCREMENTALS
!integer :: n
!
!
!write(variable_length_text,FMT="(I5)") Ngrid_text_length
!do Ngrid = 1, Ngrid_total
!
!        !The folders are named starting from 001 by increments of 1
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!
!        !We will produce a basic parent-level heat map for the grid
!        open(filechannel1,file=gridpath0//Ngrid_text//"/"//counter0file)
!        do n = 1, counter0_max
!        read(filechannel1,FMT="(I8)") counter0(n)
!        end do
!        close(filechannel1)
!
!        call mapCell(0,counter0,counter0_max,overcrowd0,bounds1,bounds2,multiplier1_0,multiplier2_0,gridpath0//Ngrid_text//"/")
!end do
!
!end subroutine analyzeHeatMaps1
!
!
!


end module analyzeHeatMapswithMultipleGrids
