    implicit doubleprecision (a-h,o-z) 
    !real*8,allocatable,dimension(:,:,:):: d,u,v,w,T
    real*8,allocatable,dimension(:,:,:):: d_ar,u_ar,v_ar,w_ar,T_ar  ! 每个点的时均变量
	real*8,allocatable,dimension(:,:,:):: dd,du,dv,dw,dT
    real*8,allocatable,dimension(:)::da,ua,va,wa,Ta

    real*8,allocatable,dimension(:,:):: tau_wb,tau_ws
    real*8,allocatable,dimension(:,:):: u_taub,u_taus
    real*8,allocatable,dimension(:,:):: Re_taub,Re_taus
    real*8,allocatable,dimension(:,:):: Cfb,Cfs
    real*8,allocatable,dimension(:)::Cfba,Cfsa

    real*8,allocatable:: xx(:),yy(:),zz(:)
    real*8 tt,us1,us2,uy,Tw,ut,h1,h2,hs
	real*8 miu_w, tau_w, u_tau, Re_tau
    real*8,allocatable,dimension(:)::yplus, uplus,Uvd
	integer istep
	integer IT
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
    print*, "nx,ny,nz=",nx,ny,nz
    print*, "Re,Ama,T_inf=",Re,Ama,T_inf 
    print*, "Tw=",Tw
    Tsb=110.4/T_inf
    miu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)   
    !Amu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)  ! 2024.5.16替换     
    p00=1.d0/(gamma*Ama*Ama)
    print*, "Tsb, miu_w, p00=",Tsb, miu_w, p00

 
    allocate(xx(nx),yy(ny),zz(nz))
    !allocate(d(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),T(nx,ny,nz)) 
    allocate(d_ar(nx,ny,nz),u_ar(nx,ny,nz),v_ar(nx,ny,nz),w_ar(nx,ny,nz),T_ar(nx,ny,nz))  ! 每个点的时均变量
    allocate(dd(nx,ny,nz),du(nx,ny,nz),dv(nx,ny,nz),dw(nx,ny,nz),dT(nx,ny,nz))
    allocate(da(ny),ua(ny),va(ny),wa(ny),Ta(ny))
    allocate(yplus(ny),uplus(ny),Uvd(ny)) 
    allocate(tau_wb(nx,nz),tau_ws(nx,nz))
    allocate(u_taub(nx,nz),u_taus(nx,nz))
    allocate(Re_taub(nx,nz),Re_taus(nx,nz))
    allocate(Cfb(nx,nz),Cfs(nx,nz))
    allocate(Cfba(nx),Cfsa(nx))
!------------------------------------------------------
	PI=3.14159265358979d0
    hx = (SLx) / (nx-1)
    hy = (SLy) / (ny-1)   
    hz = (SLz) / (nz-1)
!---------read grid message --------------------------
	  
    open(66,file="../ocfd-grid.dat",form="unformatted")
    !read(66) xx   ! X和Z方向为均匀网格，不需要读入
    read(66) yy    ! 读入Y方向的均匀网格
    close(66)
    
    yy=yy+1    ! 根据y区域选择是否加这句      
    Do i=1,nx 
        xx(i)=(i-1.d0)/(nx-1.d0)*SLx
    enddo  

    Do i=1,nz 
        zz(i)=(i-1.d0)/(nz-1.d0)*SLz
    enddo 
  
    Istep=0
    tt=0.d0

!---------read series.txt--------------------------
      
 
    open(19,file="series.txt",form="formatted")
    read(19,*)    
    read(19,*) temp,IT,temp
    close(19)  

   
!---------read data file------------------------
          
         
    write(Cnum,"(I8.8)") mod(IT,100000000)     
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

        print*, "read data OK ..."
    else
        print*, ".dat.average data is missing !!"
        stop
    endif

!---------read Favre-average data------------------------

         
    filename='OCFDFavreAverage.dat'
    inquire(file=filename,exist=file_exist)          
    if(file_exist) then
        print*,"ok,",filename
        open(33,file=filename,form="unformatted")
        read(33) Istep,tt
        print*,Istep,tt
        do k=1,nz
            read(33) dd(:,:,k)
        enddo
        do k=1,nz
            read(33) du(:,:,k)
        enddo
        do k=1,nz
            read(33) dv(:,:,k)
        enddo
        do k=1,nz
            read(33) dw(:,:,k)
        enddo
        do k=1,nz
            read(33) dT(:,:,k)
        enddo
        close(33)
        
        print*, "read Favre-average data OK ..."
    else
        print*, "Favre-averaged data is missing !!"
        stop
    endif

!---------make spacial average------------------------

        
    da=0.d0;ua=0.d0;va=0.d0;wa=0.d0;Ta=0.d0
		
	do j=1,ny  
        do k=1,nz
            do i=1,nx
                da(j)=da(j)+d_ar(i,j,k)
                ua(j)=ua(j)+u_ar(i,j,k)
                va(j)=va(j)+v_ar(i,j,k)
                wa(j)=wa(j)+w_ar(i,j,k)
                Ta(j)=Ta(j)+T_ar(i,j,k)
            enddo 	 
        enddo
	    da(j)=da(j)/(nx*nz)
	    ua(j)=ua(j)/(nx*nz)
	    va(j)=va(j)/(nx*nz)
	    wa(j)=wa(j)/(nx*nz)
	    Ta(j)=Ta(j)/(nx*nz)
    enddo
	
	db=0.d0;Ub=0.d0
        
    do j=1,ny    
        db=db+da(j)   
    enddo
        
    db=db/ny                      ! bulk density

	do k=1,nz  
        do j=1,ny
            do i=1,nx
                Ub=Ub+du(i,j,k)   ! <ρU>={U}*<ρ>
            enddo
        enddo 
    enddo
	
        
    Ub=Ub/(nx*ny*nz)/db           ! bulk mean velocity
	print*, "rho_bulk,U_bulk=",db,Ub



!---------compute u+ , y+ and Re_τ ------------------------

    Re_taub_ave = 0.d0; Re_taus_ave = 0.d0
	
    do k=1,nz
        do i=1,nx
            us1=u_ar(i,2,k) ;  us2=u_ar(i,3,k) ; h1=yy(2)-yy(1) ; h2=yy(3)-yy(1)	  
            tau_wb(i,k)=miu_w*(h2*h2*us1-h1*h1*us2)/(h2*h2*h1-h1*h1*h2)  ! 2nd order,grad(u)	
            !tau_wb = miu_w * (ua(2)-ua(1))/(yy(2)-yy(1))      ! 1st order,grad(u)
            u_taub(i,k) = sqrt(tau_wb(i,k)/d_ar(i,1,k))        ! 密度取吹气壁面处的密度
            Re_taub(i,k) = u_taub(i,k)*SLy*d_ar(i,1,k)/(2.d0*miu_w)    ! 变量均取吹气壁面处的值
            Re_taub_ave = Re_taub_ave + Re_taub(i,k)

	  
            us1=u_ar(i,ny-1,k) ;  us2=u_ar(i,ny-2,k) ; h1=yy(ny-1)-yy(ny) ; h2=yy(ny-2)-yy(ny)	  
            tau_ws(i,k)=-miu_w*(h2*h2*us1-h1*h1*us2)/(h2*h2*h1-h1*h1*h2)  ! 2nd order,grad(u)
            !tau_ws = miu_w * (ua(ny-1)-ua(ny))/(yy(ny)-yy(ny-1))      ! 1st order,grad(u) 
            u_taus(i,k) = sqrt(tau_ws(i,k)/d_ar(i,ny,k))              ! 密度取吸气壁面处的密度
            Re_taus(i,k) = u_taus(i,k)*SLy*d_ar(i,ny,k)/(2.d0*miu_w)  ! 变量均取吸气壁面处的值
            Re_taus_ave = Re_taus_ave + Re_taus(i,k)

            
            !Re_tau = (Re_taub+Re_taus)/2.d0
        enddo
	enddo

    Re_taub_ave=Re_taub_ave/(nx*nz)
	Re_taus_ave=Re_taus_ave/(nx*nz)

	print*, "Re_taub_ave, Re_taus_ave=", Re_taub_ave, Re_taus_ave



!--------------------- compute Cf ------------------------


    ! Cf摩擦阻力系数 Cf=2×τw/(ρw×U^2) 
    Cfb_ave=0.d0; Cfs_ave=0.d0
    Cfba=0.d0;Cfsa=0.d0
        
    do i=1,nx
        do k=1,nz
            Cfb(i,k) = 2*tau_wb(i,k)/(d_ar(i,1,k)*Ub**2)
            Cfs(i,k) = 2*tau_ws(i,k)/(d_ar(i,ny,k)*Ub**2)
 
            Cfb_ave = Cfb_ave + Cfb(i,k)
            Cfs_ave = Cfs_ave + Cfs(i,k)
            !Cf  = (Cfb+Cfs)/2.d0
            
            Cfba(i)=Cfba(i)+Cfb(i,k)
            Cfsa(i)=Cfsa(i)+Cfs(i,k)
        enddo
    enddo

    Cfb_ave=Cfb_ave/(nx*nz)
	Cfs_ave=Cfs_ave/(nx*nz)
	Cf_ave = (Cfb_ave+Cfs_ave)/2.d0
    
    Cfba=Cfba/nz
    Cfsa=Cfsa/nz

    print*, "Cfb_ave, Cfs_ave=",Cfb_ave,Cfs_ave
	print*, "Cf_average=",Cf_ave
    

!--------------------- output Cf ------------------------    
    open(44,file="Cf_channel_XY.plt",form="formatted")
    write(44,*) "variables=x,z,C<sub>fb</sub>,C<sub>fs</sub>"
    write(44,*) "zone i= ",nx, " j= ",nz
    do k=1,nz
        do i=1,nx
            write(44,*) xx(i),zz(k),Cfb(i,k),Cfs(i,k) 
        enddo
    enddo  
    close(44)
    
    
    open(45,file="Cf_channel_line.dat",form="formatted")
    write(45,*) "variables=x,C<sub>fb</sub>,C<sub>fs</sub>"  
    do i=1,nx    
        write(45,*) xx(i),Cfba(i),Cfsa(i) 
    enddo
    close(45)
    
!---------compute R , S and G ------------------------

! R 减阻率       R=(Cf0-Cfb/2-Cfs/2)/Cf0

    print*, "Input Cf0:"
	read(*,*) Cf0
    R  = (Cf0-(Cfb_ave+Cfs_ave)/2.d0)/Cf0
	
	print*, "Drag reduction rate: R=",R

! Win 吹气功率   Win=0.5 × ρw × vw^3
        
	print*, "Input blowing velocity vw: "
    read(*,*) vw
    Win=0.5*da(1)*vw**3

    print*, "Blowing power: Win=",Win

! S 净节能率     S=(Cf0-Cfb-Win/U)/Cf0

    S = (Cf0-Cfb_ave-Win/Ub)/Cf0

	print*, "Net-energy saving rate: S=",S


! G 控制增益     G=U×(Cf0-Cfb)/Win

    G=Ub*(Cf0-Cfb_ave)/Win

	print*, "Control Gain: G=",G




!---------output control effect parameters ------------------------

      open(49,file="ControlEffect.txt",form="formatted")
      write(49,*)"rho_bulk, U_bulk =",db,Ub
      
      write(49,*)

      !write(49,*) "rho_wb, rho_ws =",da(1),da(ny)
      !write(49,*) "tau_wb, tau_ws =",tau_wb, tau_ws
      !write(49,*) "u_taub, u_taus =",u_taub, u_taus
      !write(49,*) "Re_taub, Re_taus =",Re_taub,Re_taus
      !write(49,*) "Re_tau =",Re_tau
      write(49,*)"Re_taub_ave, Re_taus_ave=", Re_taub_ave, Re_taus_ave
      

      
      write(49,*)

      write(49,*) "Cfb_ave, Cfs_ave =",Cfb_ave,Cfs_ave 
      write(49,*) "Cf_average =",Cf_ave
      write(49,*) "Drag reduction rate: R =",R
      write(49,*) "Blowing power: Win =",Win
      write(49,*) "Net-energy saving rate: S =",S
      write(49,*) "Control Gain: G =",G
 

      close(49)


        



end






!!---------output u+ , y+ and Uvd+ ------------------------
!
!        
!        savename='UYplus'//Cnum//'.dat'
!        open(44,file=savename,form="formatted")
!	write(44,*) 'Title = "Comput." '
!        write(44,*) "variables=y<sup>+</sup>,<greek>`</greek>u<sup>+</sup>, &
!	    <greek>`</greek>u<sub>VD</sub><sup>+</sup>"
!	write(44,*) 'ZONE T="Comput." '
!     	!write(44,*) "zone i= ",nx, " j= ",ny
!	    
!        do j=2,ny-1
!       
!	   write(44,"(3f16.8)") yplus(j), uplus(j), Uvd(j)
!	   
!        enddo
!       
!	    close(44)
!	    
