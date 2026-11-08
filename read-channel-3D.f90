          implicit doubleprecision (a-h,o-z) 
          real*8,allocatable,dimension(:,:,:):: d,u,v,w,T,d1,u1,v1,w1,T1
          real*8,allocatable:: xx(:),yy(:),zz(:)
          real*8 tt,us1,us2,uy,Tw,ut,h1,h2,hs
	      integer istep
	      integer IT
	      integer Cstart,Cend,Cinterval
	      character(len=120) filename
	      character(len=120) savename
	      character(len=8) Cnum
	      logical file_exist
!-----------------------------------------
          gamma=1.4d0
          open(33,file="../opencfd.in")
          do k=1,4
          read(33,*)
          enddo
          read(33,*) nx,ny,nz ! 第5行开始读
          do k=1,7
          read(33,*)
          enddo
	  read(33,*) SLx,SLy,SLz  ! 第13行开始读
	  read(33,*)
	 
          read(33,*) Re,Ama,tmp,tmp,tmp,tmp,T_inf  ! 第15行开始读
          do k=1,13
		  read(33,*)
		  enddo
		  read(33,*) tmp,tmp,tmp,tpm,Tw    ! 第29行开始读
		  close(33)  
          print*, "nx,ny=",nx,ny,nz
		  print*, "Re,Ama,T_inf=",Re,Ama,T_inf
          print*, "Tw=",Tw

          Tsb=110.4/T_inf
          Amu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
          p00=1.d0/(gamma*Ama*Ama)

       allocate(xx(nx),yy(ny),zz(nz))
       allocate(d(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),T(nx,ny,nz))
!------------------------------------------------------
	   PI=3.14159265358979d0
	   hx = (SLx) / (nx-1)
	   hy = (SLy) / (ny-1)
       hz = (SLz) / (nz-1)
!---------read grid message --------------------------
	  open(66,file="../ocfd-grid.dat",form="unformatted")
	  !read(66) xx   ! X和Z方向为均匀网格，不需要读入
	  read(66) yy    ! 读入Y方向的非均匀网格
	  close(66)

          
      Do i=1,nx
          xx(i)=(i-1.d0)/(nx-1.d0)*SLx
      enddo  
    
      Do i=1,nz
          zz(i)=(i-1.d0)/(nz-1.d0)*SLz
      enddo 

       !   DO i= 1,nx
	      !xx(i)= -1*SLx/2 + (i-1)*hx
	      !
       !   ENDDO
	  
	      !DO j= 1,ny
	      !yy(j)= (j-1)*hy
	      
          !ENDDO
      Istep=0
	  tt=0.d0



!---------read Series.txt--------------------------
      
      
      open(19,file="series.txt",form="formatted")
      read(19,*)    
      read(19,*) Cstart,Cend,Cinterval
      close(19)  
        
       
!---------read data file------------------------
      do IT=Cstart,Cend,Cinterval
          
         write(Cnum,"(I8.8)") mod(IT,100000000)
         filename='../OCFD'//Cnum//'.dat'
         inquire(file=filename,exist=file_exist)          
         if(file_exist) then
         print*,"ok,",filename
         open(33,file=filename,form="unformatted")
         read(33) Istep,tt
         print*,Istep,tt
		 do k=1,nz
	     read(33) d(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) u(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) v(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) w(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) T(:,:,k)
	     enddo
         close(33)
      
        
         print*, "read data OK ..."
         print*, "Istep,tt=",Istep,tt
        
        savename='flow3d'//Cnum//'.plt'
        open(44,file=savename,form="formatted")
        write(44,*) "variables=x,y,z,<greek>r</greek>,u,v,w,T,p"
     	write(44,*) "zone i= ",nx, " j= ",ny, " k= ",nz
	    
        do k=1,nz
	do j=1,ny
        do i=1,nx
	       write(44,"(9f16.8)") xx(i),yy(j),zz(k),d(i,j,k),u(i,j,k),v(i,j,k),&
               w(i,j,k),T(i,j,k), d(i,j,k)*T(i,j,k)*p00
	    enddo
        enddo
        enddo
	
	    close(44)
	    
	    endif
	    
	   !i0=1000    !1545,1470,1395,1245,1095,945,  
    !   print*,"x=",xx(i0)
    !   
    !   hs=yy(2)-yy(1)
    !   print*,u(10,2)
	   !us1=u(i0,2) ;  us2=u(i0,3) ; h1=yy(2) ; h2=yy(3)
	   !uy=(h2*h2*us1-h1*h1*us2)/(h2*h2*h1-h1*h1*h2)  ! 2nd order
	   !tw=Amu_w*uy           !壁面切应力
	   !ut=sqrt(tw/d(i0,1))    !壁面摩擦系数
    !   yp=d(i0,1)*ut/Amu_w
    !   print*,yp
	   ! 
    !   open(55,file='us.dat')     
    !   do j=2,ny-1
    !      yp=yy(j)*d(i0,1)*ut/Amu_w       !为什么还多乘以一个密度
	   !   write(55,"(2f15.6)") yp,u(i0,j)/ut
	   !enddo
	   !close(55) 
	  
	  enddo
        
 
 !----------------------------------------     
      end
