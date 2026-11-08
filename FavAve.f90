     implicit doubleprecision (a-h,o-z) 
          real*8,allocatable,dimension(:,:,:):: d,u,v,w,T             ! 瞬时流场数据
          real*8,allocatable,dimension(:):: d_a,u_a,v_a,w_a,T_a   ! Reynolds平均
          real*8,allocatable,dimension(:,:,:):: d_ar,u_ar,v_ar,w_ar,T_ar   ! Reynolds平均读入
	  real*8,allocatable,dimension(:,:,:):: dd,du,dv,dw,dT   
	  real*8,allocatable,dimension(:,:,:):: duu,dvv,dww
          real*8,allocatable,dimension(:):: df,uf,vf,wf,Tf  !单个变量的Favre空间平均值
          real*8,allocatable,dimension(:):: uf_rms,vf_rms,wf_rms  ! Favre平均的RMS   
          real*8,allocatable:: xx(:),yy(:),zz(:)
          real*8 tt,us1,us2,uy,Tw,ut,h1,h2,hs
	  real*8 Tsb,T_inf, miu_w,p00
	  real*8 tau_w, Re_tau, u_tau
          real*8,allocatable,dimension(:)::yplus,uplus
	 
	      integer Istep
	      integer IT
	      integer K_average
	      integer Cstart,Cend,Cinterval
	      character(len=120) filename
	      character(len=120) savename
	      character(len=8) Cnum
	      logical file_exist
!-----------------------------------------
          gamma=1.4d0
          open(32,file="../opencfd.in")
          do k=1,4
          read(32,*)
          enddo
          read(32,*) nx,ny,nz ! 第5行开始读
          do k=1,7
          read(32,*)
          enddo
	  read(32,*) SLx,SLy,SLz  ! 第13行开始读
	  read(32,*)
	 
          read(32,*) Re,Ama,tmp,tmp,tmp,tmp,T_inf  ! 第15行开始读
          do k=1,13
		  read(32,*)
		  enddo
		  read(32,*) tmp,tmp,tmp,tpm,Tw    ! 第29行开始读
		  close(32)  
          print*, "nx,ny=",nx,ny,nz
		  print*, "Re,Ama,T_inf=",Re,Ama,T_inf
          print*, "Tw=",Tw

          Tsb=110.4/T_inf
	  miu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
          !Amu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
          p00=1.d0/(gamma*Ama*Ama)

       allocate(xx(nx),yy(ny),zz(nz))
       allocate(d(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),T(nx,ny,nz))
       allocate(d_a(ny),u_a(ny),v_a(ny),w_a(ny),T_a(ny))
       allocate(d_ar(nx,ny,nz),u_ar(nx,ny,nz),v_ar(nx,ny,nz),w_ar(nx,ny,nz),T_ar(nx,ny,nz))
       allocate(dd(nx,ny,nz),du(nx,ny,nz),dv(nx,ny,nz),dw(nx,ny,nz),dT(nx,ny,nz))
       allocate( duu(nx,ny,nz),dvv(nx,ny,nz),dww(nx,ny,nz) )
       allocate(yplus(ny),uplus(ny))
       allocate(df(ny),uf(ny),vf(ny),wf(ny),Tf(ny))
       allocate(uf_rms(ny),vf_rms(ny),wf_rms(ny))
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

!--------- read Series.txt --------------------------
      
      
      open(19,file="series.txt",form="formatted")
      read(19,*)    
      read(19,*) Cstart,Cend,Cinterval
      close(19)


!---------read .average data------------------------
          
         write(Cnum,"(I8.8)") mod(Cend,100000000)
         filename='../OCFD'//Cnum//'.dat.average'
         inquire(file=filename,exist=file_exist)          
         if(file_exist) then
         print*,"ok,",filename
         open(33,file=filename,form="unformatted")
         read(33) Istep,tt
         print*,Istep,tt
		 do k=1,nz
	     read(33) d_ar(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) u_ar(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) v_ar(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) w_ar(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) T_ar(:,:,k)
	     enddo
         close(33)

         d_a=0.d0;u_a=0.d0;v_a=0.d0;w_a=0.d0;T_a=0.d0
	 
	 ! 做对.average数据做XZ平面空间平均
	do j=1,ny  
	  do k=1,nz
          do i=1,nx
	    d_a(j)=d_a(j)+d_ar(i,j,k)
	    u_a(j)=u_a(j)+u_ar(i,j,k)
	    v_a(j)=v_a(j)+v_ar(i,j,k)
	    w_a(j)=w_a(j)+w_ar(i,j,k)
	    T_a(j)=T_a(j)+T_ar(i,j,k)
	  enddo 
	  enddo
	    d_a(j)=d_a(j)/(nx*nz)
	    u_a(j)=u_a(j)/(nx*nz)
	    v_a(j)=v_a(j)/(nx*nz)
	    w_a(j)=w_a(j)/(nx*nz)
	    T_a(j)=T_a(j)/(nx*nz)
	enddo
      
        
         print*, "read data OK ..."
	 endif


!--------- do Favre average ------------------------

      K_average = 0
      dd= 0.d0;du= 0.d0;dv= 0.d0;dw= 0.d0;dT= 0.d0
      duu= 0.d0;dvv= 0.d0;dww= 0.d0


    do IT=Cstart,Cend,Cinterval
          
         write(Cnum,"(I8.8)") mod(IT,100000000)
         filename='../OCFD'//Cnum//'.dat'
         inquire(file=filename,exist=file_exist)          
          if(file_exist) then
          print*,"ok,",filename
          open(34,file=filename,form="unformatted")
          read(34) Istep,tt
          print*,Istep,tt
		 do k=1,nz
	     read(34) d(:,:,k)
	     enddo
		 do k=1,nz
	     read(34) u(:,:,k)
	     enddo
		 do k=1,nz
	     read(34) v(:,:,k)
	     enddo
		 do k=1,nz
	     read(34) w(:,:,k)
	     enddo
		 do k=1,nz
	     read(34) T(:,:,k)
	     enddo
          close(34)
          endif
        
         print*, "read data OK ..."
         print*, "Istep,tt=",Istep,tt

      
	     do k=1,nz
             do j=1,ny
             do i=1,nx     
             
	         dd(i,j,k) = (k_average*dd(i,j,k)+d(i,j,k)*d(i,j,k))/(K_average+1.d0)
		 du(i,j,k) = (k_average*du(i,j,k)+d(i,j,k)*u(i,j,k))/(K_average+1.d0)
		 dv(i,j,k) = (k_average*dv(i,j,k)+d(i,j,k)*v(i,j,k))/(K_average+1.d0)
		 dw(i,j,k) = (k_average*dw(i,j,k)+d(i,j,k)*w(i,j,k))/(K_average+1.d0)
		 dT(i,j,k) = (k_average*dT(i,j,k)+d(i,j,k)*T(i,j,k))/(K_average+1.d0)
		 ! 求<ρρ>、<ρu>、<ρv>、<ρw>、<ρT>
		 !IF (i==nx/2 .and. j==ny/2 .and. k==nz/2) THEN
                      !print*, "du,dv,dw=",du(i,j,k),dv(i,j,k),dw(i,j,k)
		 !ENDIF
		 

                 
		 duu(i,j,k) = (k_average*duu(i,j,k)+d(i,j,k)*u(i,j,k)*u(i,j,k))/(K_average+1.d0)
		 dvv(i,j,k) = (k_average*dvv(i,j,k)+d(i,j,k)*v(i,j,k)*v(i,j,k))/(K_average+1.d0)
		 dww(i,j,k) = (k_average*dww(i,j,k)+d(i,j,k)*w(i,j,k)*w(i,j,k))/(K_average+1.d0)
		 ! 求<ρuu>、<ρvv>、<ρww>

             enddo
	     enddo
             enddo 
	     k_average = k_average + 1
	    
    enddo    


!---------make spacial average------------------------

        df=0.d0;uf=0.d0;vf=0.d0;wf=0.d0;Tf=0.d0
	! Favre average
	do j=1,ny  
	  do k=1,nz
          do i=1,nx
	    df(j)=df(j)+dd(i,j,k)/d_a(j)
	    uf(j)=uf(j)+du(i,j,k)/d_a(j)
	    vf(j)=vf(j)+dv(i,j,k)/d_a(j)
	    wf(j)=wf(j)+dw(i,j,k)/d_a(j)
	    Tf(j)=Tf(j)+dT(i,j,k)/d_a(j)
	  enddo 
	  enddo
	    df(j)=df(j)/(nx*nz)
	    uf(j)=uf(j)/(nx*nz)
	    vf(j)=vf(j)/(nx*nz)
	    wf(j)=wf(j)/(nx*nz)
	    Tf(j)=Tf(j)/(nx*nz)
	enddo
        !print*, "uf,vf,wf=",uf(ny/2),vf(ny/2),wf(ny/2)



        uf_rms=0.d0;vf_rms=0.d0;wf_rms=0.d0
	! Favre RMS
	do j=1,ny  
	  do k=1,nz
          do i=1,nx
	    uf_rms(j)=uf_rms(j)+sqrt( abs(duu(i,j,k)/d_a(j)-du(i,j,k)/d_a(j)* &
	                du(i,j,k)/d_a(j)) )
	    vf_rms(j)=vf_rms(j)+sqrt( abs(dvv(i,j,k)/d_a(j)-dv(i,j,k)/d_a(j)* &
	                dv(i,j,k)/d_a(j)) )
	    wf_rms(j)=wf_rms(j)+sqrt( abs(dww(i,j,k)/d_a(j)-dw(i,j,k)/d_a(j)* &
	                dw(i,j,k)/d_a(j)) )
	  enddo
	  enddo
	    uf_rms(j)=uf_rms(j)/(nx*nz)
	    vf_rms(j)=vf_rms(j)/(nx*nz)
	    wf_rms(j)=wf_rms(j)/(nx*nz)
	enddo




!---------output u+ , y+ and normal Reynolds stress ------------------------

        tau_w = miu_w * (u_a(2)-u_a(1))/(yy(2)-yy(1)) 
	u_tau = sqrt(tau_w/d_a(1))              ! 密度取壁面处的密度
	Re_tau = u_tau*SLy*d_a(1)/(2.d0*miu_w)  ! 变量均取壁面处的值
        print*, "rho_w=",d_a(1)
	print*, "tau_w, u_tau, Re_tau=",tau_w, u_tau, Re_tau

	do j=1,ny  
	  yplus(j) = Re_tau *2 / SLy * yy(j)
	  uplus(j) = u_a(j) / u_tau
	enddo

	do j=1,ny  
	  df(j)=df(j)
	  uf(j)=uf(j)/u_tau
      vf(j)=vf(j)/u_tau
	  wf(j)=wf(j)/u_tau
	  Tf(j)=Tf(j)
	enddo
	! 注意！此处密度ρ和温度T没有做归一化处理

        do j=1,ny  
	  uf_rms(j)=uf_rms(j)/u_tau
          vf_rms(j)=vf_rms(j)/u_tau
	  wf_rms(j)=wf_rms(j)/u_tau
	enddo



!--------- save files ------------------------
	     
	     ! 输出Favre平均的二进制文件
	     savename='OCFDFavreAverage.dat'
	     	     
	     open(43,file=savename,form="unformatted")
	     write(43) Istep,tt
		 do k=1,nz
	     write(43) ( (dd(i,j,k)/d_ar(i,j,k),i=1,nx), j=1,ny)
	     enddo
		 do k=1,nz
	     write(43) ( (du(i,j,k)/d_ar(i,j,k),i=1,nx), j=1,ny)
	     enddo
		 do k=1,nz
	     write(43) ( (dv(i,j,k)/d_ar(i,j,k),i=1,nx), j=1,ny)
	     enddo
		 do k=1,nz
	     write(43) ( (dw(i,j,k)/d_ar(i,j,k),i=1,nx), j=1,ny)
	     enddo
		 do k=1,nz
	     write(43) ( (dT(i,j,k)/d_ar(i,j,k),i=1,nx), j=1,ny)
	     enddo
       
	    close(43)
	     
	     
	     savename='Favre-Average.dat'
	     	     
	     open(44,file=savename,form="formatted")
	write(44,*) 'Title = "Comput." '
        write(44,'(A160)',advance='no') "variables=y<sup>+</sup>,{<greek>r</greek>},{u<sup>+</sup>''u<sup>+</sup>''}, &
	    {v<sup>+</sup>''v<sup>+</sup>''},{w<sup>+</sup>''w<sup>+</sup>'',{T}}"
	write(44,*) 'ZONE T="Comput." '
     	!write(44,*) "zone i= ",nx, " j= ",ny
	    
        do j=2,ny-1
       
	   write(44,"(6f16.8)") yplus(j), df(j),uf(j),vf(j),wf(j),Tf(j)
	   
        enddo
       
	    close(44)

             




             savename='Favre-RMS.dat'

	     open(45,file=savename,form="formatted")
	write(45,*) 'Title = "Comput." '
        write(45,'(A130)',advance='no') "variables=y<sup>+</sup>,{u<sup>+</sup>''<sub>rms</sub>}, &
	    {v<sup>+</sup>''<sub>rms</sub>},{w<sup>+</sup>''<sub>rms</sub>'}"
	write(45,*) 'ZONE T="Comput." '
     	!write(45,*) "zone i= ",nx, " j= ",ny
	    
        do j=2,ny-1
       
	   write(45,"(4f16.8)") yplus(j), uf_rms(j),vf_rms(j),wf_rms(j)
	   
        enddo
       
	    close(45)



            
             savename='Reynolds-average.dat'

	     open(46,file=savename,form="formatted")
	write(46,*) 'Title = "Comput." '
        write(46,'(A100)',advance='no') "variables=y<sup>+</sup>,<<greek>r</greek>>,<u<sup>+</sup>> &
	    <v<sup>+</sup>>,<w<sup>+</sup>>,<T>"
	write(46,*) 'ZONE T="Comput." '
     	!write(46,*) "zone i= ",nx, " j= ",ny
	    
        do j=2,ny-1
       
	   write(46,"(6f16.8)") yplus(j), d_a(j),u_a(j)/u_tau,v_a(j)/u_tau,w_a(j)/u_tau,T_a(j)
	   
        enddo
       
	    close(46)

	    
	    

 !----------------------------------------     
      end
