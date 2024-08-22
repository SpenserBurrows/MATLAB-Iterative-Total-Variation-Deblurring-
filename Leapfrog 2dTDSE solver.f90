      program FinalBurrows
      implicit none 
      integer :: x,y,t,nx,ny,nt,yu,yd,xu,xd,slit,potus40,n ! define integer variables
      real*8 :: dx,dy,dt,d,w,xi,yi,sigma,x0,y0,kx,c,v0,A ! define double variables 
      real*8, ALLOCATABLE :: Rarr(:,:,:) ! define array for R(x,y,t)
      real*8, ALLOCATABLE :: Iarr(:,:,:) ! define array for I(x,y,t)
      real*8, ALLOCATABLE :: Varr(:,:) ! define array for V(x,y)
      
      762 format (a,ES14.2,x,a,ES14.2,x,a,ES14.2,x,a,ES14.2) !format statements
      556 format (ES14.8,a) !format statements
      9 format (a,i8.4,x,a,ES14.2,x,a,ES14.2) !format statements

open(78, file='t7.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(77, file='CS.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(76, file='t6.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(75, file='t5.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(74, file='t4.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(73, file='t3.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(72, file='t2.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(71, file='t1.csv',form='formatted',status='unknown')!open csv file for painless excel graphing
open(70, file='v.csv',form='formatted',status='unknown')!open csv file for painless excel graphing

potus40=2! 0=wall, 1=Mr. Gorbachev, tear down this wall!


v0=dble(10**3)!set absurdly high potential barrier
nx=200 !set number of position steps
ny=200 !nx=ny
nt=7000 !set number time steps

dx=dble(.005) !set size of position steps 
dy=dx
dt=dble(0.0000005) !set size of time steps

A=dble(10.0)!set amplitude
sigma=dble(.1)!set sigma 
kx=dble(40.0)!set k
x0=xi(50,dx)!set inital x postion for wavepacket
y0=yi(100,dy)!set inital y position for wavepacket
c=dt/(dble(2.0)*dx**2)!set constant 

d=10!slit spacing
w=10!slit width
slit=110!slit x position

print 762, 'Wavepacket starting at x0=',x0,'and y0=',y0,' with k=',kx,'and sigma=',sigma

print 9, 'Running for ',nt,'loops with dx=dy=',dx,' and dt=',dt

ALLOCATE (Varr(nx,ny)) !allocate Potential array 

!intialize potential array
n=1
do x=1,nx !x loop
do y=1,ny !y loop

    if (abs(x-slit)<5) then ! at x = slit
        if (potus40==0) then ! double slit

        if (abs(y-100) < d .or. abs(y-100) > d+w) then
            Varr(x,y)=v0
        else
        Varr(x,y)=0
        endif
        endif !for double slit

    if (potus40==2) then ! diffraction grating

            if (abs(y-n*d) < d) then
                Varr(x,y)=v0
            else
            Varr(x,y)=0
            endif

        if (abs(y-n*d)== d+w) then
        n=n+2
        endif
        else    
            Varr(x,y)=0
        endif !for grating
    else    
        Varr(x,y)=0
    endif !for x=slit



enddo !enddo x
enddo !enddo y

    if (potus40==1) then ! v=0 option for testing
        Varr=0
    endif

        do y=1,ny
            do x=1,nx
            write(70,556,advance='no') varr(x,y),','            
            enddo
            write(70,*) 
            enddo 



        t=1 !initialize real at t=1

ALLOCATE (Rarr(nx,ny,nt)) !allocate real array 

do x = 2, nx-1 ! x loop
    do y = 2, ny-1! y loop
    
        Rarr(x,y,1)=A*dexp(-((xi(x,dx)-x0)/sigma)**2)*cos(xi(x,dx)*kx)*dexp(-((yi(y,dy)-y0)/sigma)**2)

    
    enddo !y loop
    enddo !x loop
      
!initialize imaginary at t=1+1/2

ALLOCATE (Iarr(nx,ny,nt)) !allocate Imaginary array 

do x = 2,nx-1 ! x loop
    xu=x+1
    xd=x-1
    do y = 2, ny-1! y loop
        yu=y+1
        yd=y-1

        Iarr(x,y,1)=A*dexp(-((xi(x,dx)-x0)/sigma)**2)*sin(xi(x,dx)*kx)*dexp(-((yi(y,dy)-y0)/sigma)**2) !I(x,y,t=1)
        Iarr(x,y,1)=Iarr(x,y,1)+c/dble(2.0)*(Rarr(xu,y,1)+Rarr(xd,y,1)+(Rarr(x,yu,1)+Rarr(x,yd,1)-dble(4)*Rarr(x,y,1)))!I(x,y,t=1.5)

    enddo !y loop
    enddo !x loop
        

!leapfrog time!

do t = 1, nt-1 !time loop
     
    do x = 2, nx-1 ! x loop
        xu=x+1
        xd=x-1
    do y = 2, ny-1! y loop

        yu=y+1
        yd=y-1

Rarr(x,y,t+1)=Rarr(x,y,t)-c*(Iarr(xu,y,t)+Iarr(xd,y,t)+Iarr(x,yu,t)+Iarr(x,yd,t)-dble(4)*Iarr(x,y,t))-dt*varr(x,y)*Iarr(x,y,t)!R(x,y,t+dt)

Iarr(x,y,t+1)=Iarr(x,y,t)+c*(Rarr(xu,y,t)+Rarr(xd,y,t)+Rarr(x,yu,t)+Rarr(x,yd,t)-dble(4)*Rarr(x,y,t))+dt*varr(x,y)*Rarr(x,y,t) !I(x,y,t+1/2dt)


    enddo !y loop
    enddo !x loop

enddo !time loop

print*,'Writing files...'!be patient!



    t=int(5)!write to file for first time step
    print*,'t=',t!print time position

    do y=1,ny
    do x=1,nx
        write(71,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array             
    enddo
        write(71,*) 
    enddo

    t=int(nt*.25)!write to file for second time step
    print*,'t=',t!print time position

    do y=1,ny
    do x=1,nx
        write(72,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array        
    enddo
        write(72,*) 
    enddo

    t=int(nt*.8)!write to file for 3rd time step
    print*,'t=',t!print time position

    do y=1,ny
    do x=1,nx
        write(73,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array       
    enddo
        write(73,*) 
    enddo
    
    t=int(nt*.85)!write to file for 4th time step

    print*,'t=',t!print time position
    do y=1,ny
    do x=1,nx
        write(74,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array        
    enddo
        write(74,*) 
    enddo
        
  t=int(nt*.90)!write to file for 5th time step

    print*,'t=',t!print time position
  do y=1,ny
    do x=1,nx
        write(75,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array           
    enddo
        write(75,*) 
    enddo

    t=int(nt*.95)!write to file for 6th time step
    print*,'t=',t!print time position

        do y=1,ny
        do x=1,nx
            write(76,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array       
        enddo
            write(76,*) 
    enddo
           
    t=nt-1!write to file for last time step
    print*,'t=',t!print time position

    do y=1,ny
    do x=1,nx
        write(78,556,advance='no') iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2,','!write as 2d array            
    enddo
        write(78,*) 
    enddo

        x=slit+45!write to file for intesnity cross section at last time step
        write(77,*) 'y,Intensity'
        do y=1,ny
            write(77,*) yi(y,dy),',', iarr(x,y,t)*iarr(x,y,t+1)+rarr(x,y,t)**2!write for graphing          
        enddo
        
    
    DEALLOCATE(Iarr) !dealocate the imaginary array 
    DEALLOCATE(Rarr) !dealocate the real array 
    DEALLOCATE(Varr) !dealocate the potential array 

close(70) !close file
close(71) !close file
close(72) !close file
close(73) !close file
close(74) !close file
close(75) !close file
close(76) !close file
close(77) !close file
close(78) !close file

print*,'Done!'!yay!

end program FinalBurrows
 
    Double Precision FUNCTION xi(x,dx)!scale x 
    implicit none
    integer ::x ! define integer variables
    real*8 :: dx!define local variables for use in the function
    xi=dble(x)*dx         !result to return 
    return
    end

    Double Precision FUNCTION yi(y,dy)!scale y
    implicit none
    integer ::y ! define integer variables
    real*8 :: dy!define local variables for use in the function
    yi=dble(y)*dy         !result to return 
    return
    end



