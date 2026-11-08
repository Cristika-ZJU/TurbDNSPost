!========================================================================================================================================= 
!                              本程序利用OpenCFD数据将湍动能分解，得到湍动能方程的各个量
!            
!                              **归一化过程在输出时进行，所有项/(ρ*U_τ*U_τ*U_τ)，其中ρ=1
!
!=========================================================================================================================================     
     implicit doubleprecision (a-h,o-z)         
	  real*8,allocatable,dimension(:,:,:):: d,u,v,w,T,p                   ! 瞬时流场数据
         !real*8,allocatable,dimension(:,:,:):: d_fr,u_fr,v_fr,w_fr,T_fr    ! 读入Favre平均数据
          real*8,allocatable,dimension(:,:,:):: da,ua,va,wa,Ta,pa    !单个变量的Reynolds平均值
	  real*8,allocatable,dimension(:):: d_a,u_a,v_a,w_a,T_a               ! 时均流场空间平均（只和法向高度y有关）
          real*8,allocatable,dimension(:,:,:):: df,uf,vf,wf,Tf    !单个变量的Favre平均值
          real*8,allocatable,dimension(:,:,:):: uf_x,uf_y,uf_z    
	  real*8,allocatable,dimension(:,:,:):: vf_x,vf_y,vf_z   
	  real*8,allocatable,dimension(:,:,:):: wf_x,wf_y,wf_z    !Favre平均速度的空间导数 ∂{U}/∂X

	  real*8,allocatable,dimension(:,:,:):: duu,dvv,dww,duv,dvw,duw                        ! <ρU''U''>
          real*8,allocatable,dimension(:,:,:):: duuu,duuv,duuw,dvvu,dvvv,dvvw,dwwu,dwwv,dwww   ! <ρU''U''U''>
          real*8,allocatable,dimension(:,:,:):: duuu_x,duuv_y,duuw_z
	  real*8,allocatable,dimension(:,:,:):: dvvu_x,dvvv_y,dvvw_z
	  real*8,allocatable,dimension(:,:,:):: dwwu_x,dwwv_y,dwww_z                           ! ∂<ρU''U''U''>/∂X

	  real*8,allocatable,dimension(:,:,:):: u2,v2,w2           ! U''
	  real*8,allocatable,dimension(:,:,:):: ufr,vfr,wfr        ! <U''>
	  real*8,allocatable,dimension(:,:,:):: pa_x,pa_y,pa_z     ! ∂<p>/∂X
	  real*8,allocatable,dimension(:,:,:):: p1u2,p1v2,p1w2     ! <p'U''>
	  real*8,allocatable,dimension(:,:,:):: pu_x,pv_y,pw_z     ! ∂<p'U''>/∂X
	  !real*8,allocatable,dimension(:,:,:):: u2_x,v2_y,w2_y     ! ∂U''/∂X
	  real*8,allocatable,dimension(:,:,:):: p1ux,p1vy,p1wz     ! <p' ∂U''/∂X>

	  real*8,allocatable,dimension(:,:,:):: miu,miua               ! μ,<μ>
          real*8,allocatable,dimension(:,:,:):: ua_x,ua_y,ua_z     
	  real*8,allocatable,dimension(:,:,:):: va_x,va_y,va_z
	  real*8,allocatable,dimension(:,:,:):: wa_x,wa_y,wa_z     ! ∂<U>/∂X
	  real*8,allocatable,dimension(:,:,:):: sigxx,sigyy,sigzz,sigxy,sigyz,sigxz          ! <σ>
	  real*8,allocatable,dimension(:,:,:):: sig1xx,sig1yy,sig1zz,sig1xy,sig1yz,sig1xz    ! σ'
	  real*8,allocatable,dimension(:,:,:):: sigxx_x,sigyy_y,sigzz_z     
	  real*8,allocatable,dimension(:,:,:):: sigxy_x,sigxy_y,sigyz_y,sigyz_z,sigxz_x,sigxz_z   ! ∂<σ>/∂X
	  real*8,allocatable,dimension(:,:,:):: u2_x,u2_y,u2_z,v2_x,v2_y,v2_z,w2_x,w2_y,w2_z  ! ∂U''/∂X
	  real*8,allocatable,dimension(:,:,:):: siaxx,siayy,siazz,siaxy,siayz,siaxz           ! σ'表达式第三项(Favre脉动的Reynolds平均)

	  real*8,allocatable,dimension(:,:,:):: u2sxx1,u2sxy1,u2sxz1
	  real*8,allocatable,dimension(:,:,:):: v2sxy1,v2syy1,v2syz1
	  real*8,allocatable,dimension(:,:,:):: w2sxz1,w2syz1,w2szz1         ! <U''σ'>
	  real*8,allocatable,dimension(:,:,:):: u2sxx1_x,u2sxy1_y,u2sxz1_z
	  real*8,allocatable,dimension(:,:,:):: v2sxy1_x,v2syy1_y,v2syz1_z
	  real*8,allocatable,dimension(:,:,:):: w2sxz1_x,w2syz1_y,w2szz1_z   ! ∂<U''σ'>/∂X
	  real*8,allocatable,dimension(:,:,:):: sxx1u2x,sxy1u2y,sxz1u2z
	  real*8,allocatable,dimension(:,:,:):: sxy1v2x,syy1v2y,syz1v2z
	  real*8,allocatable,dimension(:,:,:):: sxz1w2x,syz1w2y,szz1w2z      !<σ' ∂U''/∂X>
!	  real*8,allocatable,dimension(:,:,:):: u2sxxx,u2sxyy,u2sxzz
!	  real*8,allocatable,dimension(:,:,:):: v2sxyx,v2syyy,v2syzz
!	  real*8,allocatable,dimension(:,:,:):: w2sxzx,w2syzy,w2szzz         !<U> ∂<σ>/∂X

          real*8,allocatable,dimension(:,:,:):: Pk,Tk,Pik1,Pik2,Pik3,Dk,Ek,Mk
	  real*8,allocatable,dimension(:):: Pka,Tka,Pik1a,Pik2a,Pik3a,Dka,Eka,Mka
	  
          real*8,allocatable:: xx(:),yy(:),zz(:)
          real*8 tt,us1,us2,uy,Tw,ut,h1,h2,hs
	  real*8 Tsb,T_inf, miu_w,p00
	  real*8 tau_w, Re_tau, u_tau
          real*8,allocatable,dimension(:)::yplus,uplus
	 
	      integer Istep, IT
	      integer i,j,k
	      integer K_average
	      integer Cstart,Cend,Cinterval
	      character(len=120) filename
	      character(len=120) savename
	      character(len=8) Cnum
	      logical file_exist
!--------------------------------------------------------------------------------------------------
          
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
	  miu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)    ! 以最少变量形式进行无量纲化
          !miu_w=(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)          ! 以来流值进行无量纲化
          !Amu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
          p00=1.d0/(gamma*Ama*Ama)

!--------------------------------------------------------------------------------------------------

       allocate(xx(nx),yy(ny),zz(nz))
       allocate(d(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),T(nx,ny,nz),p(nx,ny,nz))
      !allocate(d_fr(nx,ny,nz),u_fr(nx,ny,nz),v_fr(nx,ny,nz),w_fr(nx,ny,nz),T_fr(nx,ny,nz))
       allocate(da(nx,ny,nz),ua(nx,ny,nz),va(nx,ny,nz),wa(nx,ny,nz),Ta(nx,ny,nz),pa(nx,ny,nz))
       allocate(df(nx,ny,nz),uf(nx,ny,nz),vf(nx,ny,nz),wf(nx,ny,nz),Tf(nx,ny,nz))
 
       allocate( duu(nx,ny,nz),dvv(nx,ny,nz),dww(nx,ny,nz) )
       allocate( duv(nx,ny,nz),dvw(nx,ny,nz),duw(nx,ny,nz) )
       allocate( uf_x(nx,ny,nz),uf_y(nx,ny,nz),uf_z(nx,ny,nz) )
       allocate( vf_x(nx,ny,nz),vf_y(nx,ny,nz),vf_z(nx,ny,nz) )
       allocate( wf_x(nx,ny,nz),wf_y(nx,ny,nz),wf_z(nx,ny,nz) )

       allocate( duuu(nx,ny,nz),duuv(nx,ny,nz),duuw(nx,ny,nz) )
       allocate( dvvu(nx,ny,nz),dvvv(nx,ny,nz),dvvw(nx,ny,nz) )
       allocate( dwwu(nx,ny,nz),dwwv(nx,ny,nz),dwww(nx,ny,nz) )
       allocate( duuu_x(nx,ny,nz),duuv_y(nx,ny,nz),duuw_z(nx,ny,nz) )
       allocate( dvvu_x(nx,ny,nz),dvvv_y(nx,ny,nz),dvvw_z(nx,ny,nz) )
       allocate( dwwu_x(nx,ny,nz),dwwv_y(nx,ny,nz),dwww_z(nx,ny,nz) ) 

       allocate( u2(nx,ny,nz),v2(nx,ny,nz),w2(nx,ny,nz) )          ! U''
       allocate( ufr(nx,ny,nz),vfr(nx,ny,nz),wfr(nx,ny,nz))        ! <U''>
       allocate( pa_x(nx,ny,nz),pa_y(nx,ny,nz),pa_z(nx,ny,nz) )    ! ∂<p>/∂X
       allocate( p1u2(nx,ny,nz),p1v2(nx,ny,nz),p1w2(nx,ny,nz) )    ! <p'U''>
       allocate( pu_x(nx,ny,nz),pv_y(nx,ny,nz),pw_z(nx,ny,nz) )    ! ∂<p'U''>/∂X
       !allocate( u2_x(nx,ny,nz),v2_y(nx,ny,nz),w2_y(nx,ny,nz) )    ! ∂U''/∂X
       allocate( p1ux(nx,ny,nz),p1vy(nx,ny,nz),p1wz(nx,ny,nz) )    ! <p' ∂U''/∂X>

       allocate( miu(nx,ny,nz),miua(nx,ny,nz) )
       allocate( ua_x(nx,ny,nz),ua_y(nx,ny,nz),ua_z(nx,ny,nz) )
       allocate( va_x(nx,ny,nz),va_y(nx,ny,nz),va_z(nx,ny,nz) )
       allocate( wa_x(nx,ny,nz),wa_y(nx,ny,nz),wa_z(nx,ny,nz) )    ! ∂<U>/∂X
       allocate(sigxx(nx,ny,nz),sigyy(nx,ny,nz),sigzz(nx,ny,nz),sigxy(nx,ny,nz),sigyz(nx,ny,nz),sigxz(nx,ny,nz))
       allocate(sig1xx(nx,ny,nz),sig1yy(nx,ny,nz),sig1zz(nx,ny,nz),sig1xy(nx,ny,nz),sig1yz(nx,ny,nz),sig1xz(nx,ny,nz))
       allocate( sigxx_x(nx,ny,nz),sigyy_y(nx,ny,nz),sigzz_z(nx,ny,nz) )
       allocate( sigxy_x(nx,ny,nz),sigxy_y(nx,ny,nz),sigyz_y(nx,ny,nz) )
       allocate( sigyz_z(nx,ny,nz),sigxz_x(nx,ny,nz),sigxz_z(nx,ny,nz) )   ! ∂<σ>/∂X
       allocate( u2_x(nx,ny,nz),u2_y(nx,ny,nz),u2_z(nx,ny,nz) )
       allocate( v2_x(nx,ny,nz),v2_y(nx,ny,nz),v2_z(nx,ny,nz) )
       allocate( w2_x(nx,ny,nz),w2_y(nx,ny,nz),w2_z(nx,ny,nz) )            ! ∂U''/∂X
       allocate( siaxx(nx,ny,nz),siayy(nx,ny,nz),siazz(nx,ny,nz) )
       allocate( siaxy(nx,ny,nz),siayz(nx,ny,nz),siaxz(nx,ny,nz) )         ! σ'表达式第三项(Favre脉动的Reynolds平均)

       allocate( u2sxx1(nx,ny,nz),u2sxy1(nx,ny,nz),u2sxz1(nx,ny,nz) )
       allocate( v2sxy1(nx,ny,nz),v2syy1(nx,ny,nz),v2syz1(nx,ny,nz) )
       allocate( w2sxz1(nx,ny,nz),w2syz1(nx,ny,nz),w2szz1(nx,ny,nz) )        ! <U''σ'>
       allocate( u2sxx1_x(nx,ny,nz),u2sxy1_y(nx,ny,nz),u2sxz1_z(nx,ny,nz) )
       allocate( v2sxy1_x(nx,ny,nz),v2syy1_y(nx,ny,nz),v2syz1_z(nx,ny,nz) )
       allocate( w2sxz1_x(nx,ny,nz),w2syz1_y(nx,ny,nz),w2szz1_z(nx,ny,nz) )  ! ∂<U''σ'>/∂X
       allocate( sxx1u2x(nx,ny,nz),sxy1u2y(nx,ny,nz),sxz1u2z(nx,ny,nz) )
       allocate( sxy1v2x(nx,ny,nz),syy1v2y(nx,ny,nz),syz1v2z(nx,ny,nz) )
       allocate( sxz1w2x(nx,ny,nz),syz1w2y(nx,ny,nz),szz1w2z(nx,ny,nz) )      ! <σ' ∂U''/∂X>

       allocate( Pk(nx,ny,nz),Tk(nx,ny,nz),Pik1(nx,ny,nz),Pik2(nx,ny,nz))
       allocate( Pik3(nx,ny,nz),Dk(nx,ny,nz),Ek(nx,ny,nz),Mk(nx,ny,nz))
       allocate( Pka(ny),Tka(ny),Pik1a(ny),Pik2a(ny),Pik3a(ny),Dka(ny),Eka(ny),Mka(ny) )
 

       allocate(yplus(ny),uplus(ny))
       allocate(d_a(ny),u_a(ny),v_a(ny),w_a(ny),T_a(ny))

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
    
      Do k=1,nz
          zz(i)=(k-1.d0)/(nz-1.d0)*SLz
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
         print*,"read",filename,"ok"
         open(33,file=filename,form="unformatted")
         read(33) Istep,tt
         print*,Istep,tt
		 do k=1,nz
	     read(33) da(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) ua(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) va(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) wa(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) Ta(:,:,k)
	     enddo
         close(33)


	 
       ! make spacial average for time-averaged data
       do j=1,ny  
	  do k=1,nz
          do i=1,nx
	    d_a(j)=d_a(j)+da(i,j,k)
	    u_a(j)=u_a(j)+ua(i,j,k)
	    v_a(j)=v_a(j)+va(i,j,k)
	    w_a(j)=w_a(j)+wa(i,j,k)
	    T_a(j)=T_a(j)+Ta(i,j,k)
	  enddo 
	  enddo
	    d_a(j)=d_a(j)/(nx*nz)
	    u_a(j)=u_a(j)/(nx*nz)
	    v_a(j)=v_a(j)/(nx*nz)
	    w_a(j)=w_a(j)/(nx*nz)
	    T_a(j)=T_a(j)/(nx*nz)
	enddo

!---------output u+ , y+ and normal Reynolds stress ------------------------

        tau_w = miu_w * (u_a(2)-u_a(1))/(yy(2)-yy(1)) 
	u_tau = sqrt(tau_w/d_a(1))              ! 密度取壁面处的密度
	Re_tau = u_tau*SLy*d_a(1)/(2.d0*miu_w)  ! 变量均取壁面处的值
	delta_niu = miu_w/(d_a(1)*u_tau)
        print*, "rho_w=",d_a(1)
	print*, "tau_w, u_tau, Re_tau=",tau_w, u_tau, Re_tau
	print*, "delta_niu=",delta_niu

	do j=1,ny  
	  yplus(j) = Re_tau *2 / SLy * (yy(j)-yy(1))
	  uplus(j) = u_a(j) / u_tau
	enddo

!	do j=1,ny  
!	  df(j)=df(j)
!	  uf(j)=uf(j)/u_tau
!          vf(j)=vf(j)/u_tau
!	  wf(j)=wf(j)/u_tau
!	  Tf(j)=Tf(j)
!	enddo
!	! 注意！此处密度ρ和温度T没有做归一化处理

       ! Normalization for X、Y、Z
!	do k=1,nz
!	     do j=1,ny  
!             do i=1,nx
!	         xx(i)=xx(i)/delta_niu
!                 yy(j)=yy(j)/delta_niu
!	         zz(k)=zz(k)/delta_niu 
!        enddo 
!	enddo
!	enddo
!	hx = hx / delta_niu
!	hy = hy / delta_niu
!	hz = hz / delta_niu

	   


	   ! Normalization for time-averaged data
!             do k=1,nz
!	     do j=1,ny  
!             do i=1,nx
!	         da(i,j,k)=da(i,j,k)/d_a(1)
!	         ua(i,j,k)=ua(i,j,k)/u_tau
!                 va(i,j,k)=va(i,j,k)/u_tau
!	         wa(i,j,k)=wa(i,j,k)/u_tau 
!		 Ta(i,j,k)=Ta(i,j,k)/T_a(1)
!             enddo 
!	     enddo
!	     enddo





	 !print*," I am HERE !!!"	


!---------read FavreAverage.dat data------------------------
          
	     filename='OCFDFavreAverage.dat'
	     	     
	     open(33,file=filename,form="unformatted")
	     read(33) Istep,tt
		 do k=1,nz
	     read(33) df(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) uf(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) vf(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) wf(:,:,k)
	     enddo
		 do k=1,nz
	     read(33) Tf(:,:,k)
	     enddo
       
	    close(33)

	   ! Normalization for Favre-averaged data
!             do k=1,nz
!	     do j=1,ny  
!             do i=1,nx
!	         df(i,j,k)=df(i,j,k)/d_a(1)
!	         uf(i,j,k)=uf(i,j,k)/u_tau
!                 vf(i,j,k)=vf(i,j,k)/u_tau
!	         wf(i,j,k)=wf(i,j,k)/u_tau 
!		 Tf(i,j,k)=Tf(i,j,k)/T_a(1)
!             enddo 
!	     enddo
!	     enddo


!	 df=0.d0;uf=0.d0;vf=0.d0;wf=0.d0;Tf=0.d0;
!	 ! 做对Favre平均数据做XZ平面空间平均
!	do j=1,ny  
!	  do k=1,nz
!          do i=1,nx
!	    df(j)=df(j)+d_fr(i,j,k)
!	    uf(j)=uf(j)+u_fr(i,j,k)
!	    vf(j)=vf(j)+v_fr(i,j,k)
!	    wf(j)=wf(j)+w_fr(i,j,k)
!	    Tf(j)=Tf(j)+T_fr(i,j,k)
!	  enddo 
!	  enddo
!	    df(j)=df(j)/(nx*nz)
!	    uf(j)=uf(j)/(nx*nz)
!	    vf(j)=vf(j)/(nx*nz)
!	    wf(j)=wf(j)/(nx*nz)
!	    Tf(j)=Tf(j)/(nx*nz)
!	enddo
      
        
         print*, "read Favre-averaged data OK ..."

	 endif


!--------- Kinetic Energy Decomposition ------------------------

      	 call dx3d(ua,ua_x,nx,ny,nz,hx)
	 !call dy3d(ua,ua_y,nx,ny,nz,hy)
	 call dy3d(ua,ua_y,nx,ny,nz,yy)
	 call dz3d(ua,ua_z,nx,ny,nz,hz)
	 call dx3d(va,va_x,nx,ny,nz,hx)
	 !call dy3d(va,va_y,nx,ny,nz,hy)
	 call dy3d(va,va_y,nx,ny,nz,yy)
	 call dz3d(va,va_z,nx,ny,nz,hz)
	 call dx3d(wa,wa_x,nx,ny,nz,hx)
	 !call dy3d(wa,wa_y,nx,ny,nz,hy)
	 call dy3d(wa,wa_y,nx,ny,nz,yy)
	 call dz3d(wa,wa_z,nx,ny,nz,hz)


	 ! 求∂<U>/∂X


	 do k=1,nz
	 do j=1,ny
	 do i=1,nx
	    pa(i,j,k)=p00*da(i,j,k)*Ta(i,j,k)           
	    miua(i,j,k)=1.0/Re*(1.0+Tsb)*sqrt(Ta(i,j,k)**3)/(Tsb+Ta(i,j,k))
	    
	    sigxx(i,j,k)=miua(i,j,k)*(4/3*ua_x(i,j,k)-2/3*(va_y(i,j,k)+wa_z(i,j,k)))
	    sigyy(i,j,k)=miua(i,j,k)*(4/3*va_y(i,j,k)-2/3*(ua_x(i,j,k)+wa_z(i,j,k)))
	    sigzz(i,j,k)=miua(i,j,k)*(4/3*wa_z(i,j,k)-2/3*(ua_x(i,j,k)+va_y(i,j,k)))
	    sigxy(i,j,k)=miua(i,j,k)*(ua_y(i,j,k)+va_x(i,j,k))
        sigyz(i,j,k)=miua(i,j,k)*(va_z(i,j,k)+wa_y(i,j,k))
	    sigxz(i,j,k)=miua(i,j,k)*(ua_z(i,j,k)+wa_x(i,j,k))

	 enddo
	 enddo
	 enddo
	 ! 求<σ>
 

	 call dx3d(sigxx,sigxx_x,nx,ny,nz,hx)
	 !call dy3d(sigxy,sigxy_y,nx,ny,nz,hy)
	 call dy3d(sigxy,sigxy_y,nx,ny,nz,yy)
	 call dz3d(sigxz,sigxz_z,nx,ny,nz,hz)
	 call dx3d(sigxy,sigxy_x,nx,ny,nz,hx)
	 !call dy3d(sigyy,sigyy_y,nx,ny,nz,hy)
	 call dy3d(sigyy,sigyy_y,nx,ny,nz,yy)
	 call dz3d(sigyz,sigyz_z,nx,ny,nz,hz)
	 call dx3d(sigxz,sigxz_x,nx,ny,nz,hx)
	 !call dy3d(sigyz,sigyz_y,nx,ny,nz,hy)
	 call dy3d(sigyz,sigyz_y,nx,ny,nz,yy)
	 call dz3d(sigzz,sigyz_z,nx,ny,nz,hz)
	 ! 求∂<σ>/∂X

	 call dx3d(uf,uf_x,nx,ny,nz,hx)
	 !call dy3d(uf,uf_y,nx,ny,nz,hy)
	 call dy3d(uf,uf_y,nx,ny,nz,yy)
	 call dz3d(uf,uf_z,nx,ny,nz,hz)
	 call dx3d(vf,vf_x,nx,ny,nz,hx)
	 !call dy3d(vf,vf_y,nx,ny,nz,hy)
	 call dy3d(vf,vf_y,nx,ny,nz,yy)
	 call dz3d(vf,vf_z,nx,ny,nz,hz)
	 call dx3d(wf,wf_x,nx,ny,nz,hx)
	 !call dy3d(wf,wf_y,nx,ny,nz,hy)
	 call dy3d(wf,wf_y,nx,ny,nz,yy)
         call dz3d(wf,wf_z,nx,ny,nz,hz)
	 ! 求∂{U}/∂X

	 call dx3d(pa,pa_x,nx,ny,nz,hx)
	 !call dy3d(pa,pa_y,nx,ny,nz,hy)
	 call dy3d(pa,pa_y,nx,ny,nz,yy)
	 call dz3d(pa,pa_z,nx,ny,nz,hz)
	 ! 求∂<p>/∂X



      K_average = 0
      duu=0.d0;dvv=0.d0;dww=0.d0;duv=0.d0;dvw=0.d0;dvw=0.d0;
      duuu=0.d0;duuv=0.d0;duuw=0.d0;dvvu=0.d0;dvvv=0.d0;dvvw=0.d0;dwwu=0.d0;dwwv=0.d0;dwww=0.d0;
      ufr=0.d0;vfr=0.d0;wfr=0.d0;  p1u2=0.d0;p1v2=0.d0;p1w2=0.d0
      p1ux=0.d0;p1vy=0.d0;p1wz=0.d0
      siaxx=0.d0;siayy=0.d0;siazz=0.d0;siaxy=0.d0;siayz=0.d0;siaxz=0.d0
      u2sxx1=0.d0;u2sxy1=0.d0;u2szx1=0.d0
      v2sxy1=0.d0;v2syy1=0.d0;v2syz1=0.d0
      w2sxz1=0.d0;w2syz1=0.d0;w2szz1=0.d0
      sxx1u2x=0.d0;sxy1u2y=0.d0;sxz1u2z=0.d0  
      sxy1v2x=0.d0;syy1v2y=0.d0;syz1v2z=0.d0
      sxz1w2x=0.d0;syz1w2y=0.d0;szz1w2z=0.d0

! ********************* 循环开始 ************************
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

         print*, "Read instantaneous data OK ..."
         print*, "Istep,tt=",Istep,tt

       ! Normalization for instantaneous data
         do k=1,nz
	 do j=1,ny  
         do i=1,nx
!	    d(i,j,k)=d(i,j,k)/d_a(1)
!	    u(i,j,k)=u(i,j,k)/u_tau
!           v(i,j,k)=v(i,j,k)/u_tau
!	    w(i,j,k)=w(i,j,k)/u_tau
!	    T(i,j,k)=T(i,j,k)/T_a(1)

            p(i,j,k)=p00*d(i,j,k)*T(i,j,k)
	    miu(i,j,k)=1.0/Re*(1.0+Tsb)*sqrt(T(i,j,k)**3)/(Tsb+T(i,j,k))
         enddo 
	 enddo
	 enddo

        
         

         
	     do k=1,nz
	     do j=1,ny
	     do i=1,nx
	         u2(i,j,k)=u(i,j,k)-uf(i,j,k)
		 v2(i,j,k)=v(i,j,k)-vf(i,j,k)
		 w2(i,j,k)=w(i,j,k)-wf(i,j,k)
             enddo
	     enddo
	     enddo
	         ! 求U''

	         call dx3d(u2,u2_x,nx,ny,nz,hx)
	         !call dy3d(u2,u2_y,nx,ny,nz,hy)
		 call dy3d(u2,u2_y,nx,ny,nz,yy)
	         call dz3d(u2,u2_z,nx,ny,nz,hz)
	         call dx3d(v2,v2_x,nx,ny,nz,hx)
	         !call dy3d(v2,v2_y,nx,ny,nz,hy)
		 call dy3d(v2,v2_y,nx,ny,nz,yy)
	         call dz3d(v2,v2_z,nx,ny,nz,hz)
	         call dx3d(w2,w2_x,nx,ny,nz,hx)
	         !call dy3d(w2,w2_y,nx,ny,nz,hy)
		 call dy3d(w2,w2_y,nx,ny,nz,yy)
	         call dz3d(w2,w2_z,nx,ny,nz,hz)
	         ! 求∂U''/∂X
	     
	     do k=1,nz
             do j=1,ny
             do i=1,nx                 
	         duu(i,j,k) = (k_average*duu(i,j,k)+d(i,j,k)* u2(i,j,k) &
		                       *u2(i,j,k))/(K_average+1.d0)
	         dvv(i,j,k) = (k_average*dvv(i,j,k)+d(i,j,k)* v2(i,j,k) &
		                       *v2(i,j,k))/(K_average+1.d0)
	         dww(i,j,k) = (k_average*dww(i,j,k)+d(i,j,k)* w2(i,j,k) &
		                       *w2(i,j,k))/(K_average+1.d0)
		 
		 duv(i,j,k) = (k_average*duv(i,j,k)+d(i,j,k)* u2(i,j,k) &
		                       *v2(i,j,k))/(K_average+1.d0)
		 dvw(i,j,k) = (k_average*dvw(i,j,k)+d(i,j,k)* v2(i,j,k) &
		                       *w2(i,j,k))/(K_average+1.d0)
		 duw(i,j,k) = (k_average*duw(i,j,k)+d(i,j,k)* u2(i,j,k) &
		                       *w2(i,j,k))/(K_average+1.d0)
             enddo
	     enddo
             enddo 
	         ! 求<ρU''U''>



	     ! 求∂{U}/∂X，见循环之前，此处略



	     do k=1,nz
             do j=1,ny
             do i=1,nx                 
		 duuu(i,j,k) = (k_average*duuu(i,j,k)+0.5*d(i,j,k)* (u2(i,j,k))**2  &
		                       *(u2(i,j,k)))/(K_average+1.d0)
		 duuv(i,j,k) = (k_average*duuv(i,j,k)+0.5*d(i,j,k)* (u2(i,j,k))**2  &
		                       *(v2(i,j,k)))/(K_average+1.d0)
	         duuw(i,j,k) = (k_average*duuw(i,j,k)+0.5*d(i,j,k)* (u2(i,j,k))**2  &
		                       *(w2(i,j,k)))/(K_average+1.d0)
	         dvvu(i,j,k) = (k_average*dvvu(i,j,k)+0.5*d(i,j,k)* (v2(i,j,k))**2  &
		                       *(u2(i,j,k)))/(K_average+1.d0)
		 dvvv(i,j,k) = (k_average*dvvv(i,j,k)+0.5*d(i,j,k)* (v2(i,j,k))**2  &
		                       *(v2(i,j,k)))/(K_average+1.d0)
	         dvvw(i,j,k) = (k_average*dvvw(i,j,k)+0.5*d(i,j,k)* (v2(i,j,k))**2  &
		                       *(w2(i,j,k)))/(K_average+1.d0)
		 dwwu(i,j,k) = (k_average*dwwu(i,j,k)+0.5*d(i,j,k)* (w2(i,j,k))**2  &
		                       *(u2(i,j,k)))/(K_average+1.d0)
		 dwwv(i,j,k) = (k_average*dwwv(i,j,k)+0.5*d(i,j,k)* (w2(i,j,k))**2  &
		                       *(v2(i,j,k)))/(K_average+1.d0)
	         dwww(i,j,k) = (k_average*dwww(i,j,k)+0.5*d(i,j,k)* (w2(i,j,k))**2  &
		                       *(w2(i,j,k)))/(K_average+1.d0)
             enddo
	     enddo
             enddo 
		 ! 求<0.5ρU''U''U''>


		 call dx3d(duuu,duuu_x,nx,ny,nz,hx)
		 !call dy3d(duuv,duuv_y,nx,ny,nz,hy)
		 call dy3d(duuv,duuv_y,nx,ny,nz,yy)
		 call dz3d(duuw,duuw_z,nx,ny,nz,hz)
		 call dx3d(dvvu,dvvu_x,nx,ny,nz,hx)
		 !call dy3d(dvvv,dvvv_y,nx,ny,nz,hy)
                 call dy3d(dvvv,dvvv_y,nx,ny,nz,yy)
		 call dz3d(dvvw,dvvw_z,nx,ny,nz,hz)
		 call dx3d(dwwu,dwwu_x,nx,ny,nz,hx)
		 !call dy3d(dwwv,dwwv_y,nx,ny,nz,hy)
		 call dy3d(dwwv,dwwv_y,nx,ny,nz,yy)
		 call dz3d(dwww,dwww_z,nx,ny,nz,hz)
		 ! 求∂<0.5ρU''U''U''> / ∂X


             do k=1,nz
             do j=1,ny
             do i=1,nx
		 ufr(i,j,k) = (k_average*ufr(i,j,k)+ u2(i,j,k))/(K_average+1.d0)
		 vfr(i,j,k) = (k_average*vfr(i,j,k)+ v2(i,j,k))/(K_average+1.d0)
		 wfr(i,j,k) = (k_average*wfr(i,j,k)+ w2(i,j,k))/(K_average+1.d0)
	     enddo
	     enddo
             enddo
		 ! 求<U''>


	     ! 求∂<p>/∂X，见循环之前，此处略


             do k=1,nz
             do j=1,ny
             do i=1,nx
		 p1u2(i,j,k) = (k_average*p1u2(i,j,k)+(p(i,j,k)-pa(i,j,k))*u2(i,j,k))/(K_average+1.d0)
		 p1v2(i,j,k) = (k_average*p1v2(i,j,k)+(p(i,j,k)-pa(i,j,k))*v2(i,j,k))/(K_average+1.d0)
		 p1w2(i,j,k) = (k_average*p1w2(i,j,k)+(p(i,j,k)-pa(i,j,k))*w2(i,j,k))/(K_average+1.d0)
	     enddo
	     enddo
             enddo
		 ! 求∂<p'U''>

		 call dx3d(p1u2,pu_x,nx,ny,nz,hx)
		 !call dy3d(p1v2,pv_y,nx,ny,nz,hy)
		 call dy3d(p1v2,pv_y,nx,ny,nz,yy)
		 call dz3d(p1w2,pw_z,nx,ny,nz,hz)
		 ! 求∂<p'U''>/∂X



		 call dx3d(u2,u2_x,nx,ny,nz,hx)
		 !call dy3d(u2,u2_y,nx,ny,nz,hy)
		 call dy3d(u2,u2_y,nx,ny,nz,yy)
		 call dz3d(u2,u2_z,nx,ny,nz,hz)
		 call dx3d(v2,v2_x,nx,ny,nz,hx)
		 !call dy3d(v2,v2_y,nx,ny,nz,hy)
		 call dy3d(v2,v2_y,nx,ny,nz,yy)
		 call dz3d(v2,v2_z,nx,ny,nz,hz)
		 call dx3d(w2,w2_x,nx,ny,nz,hx)
		 !call dy3d(w2,w2_y,nx,ny,nz,hy)
		 call dy3d(w2,w2_y,nx,ny,nz,yy)
		 call dz3d(w2,w2_z,nx,ny,nz,hz)
		 ! 求∂U''/∂X


             do k=1,nz
             do j=1,ny
             do i=1,nx
		 p1ux(i,j,k) = (k_average*p1ux(i,j,k)+(p(i,j,k)-pa(i,j,k))*u2_x(i,j,k))/(K_average+1.d0)
		 p1vy(i,j,k) = (k_average*p1vy(i,j,k)+(p(i,j,k)-pa(i,j,k))*v2_y(i,j,k))/(K_average+1.d0)
		 p1wz(i,j,k) = (k_average*p1wz(i,j,k)+(p(i,j,k)-pa(i,j,k))*w2_z(i,j,k))/(K_average+1.d0)
	     enddo
	     enddo
             enddo
		 ! 求<p' ∂U''/∂X>




             do k=1,nz
             do j=1,ny
             do i=1,nx
		 siaxx(i,j,k)=(k_average*siaxx(i,j,k)+miu(i,j,k)*(4/3*u2_x(i,j,k) &
		                       -2/3*(v2_y(i,j,k)+w2_z(i,j,k))))/(K_average+1.d0)
		 siayy(i,j,k)=(k_average*siayy(i,j,k)+miu(i,j,k)*(4/3*v2_y(i,j,k) &
		                       -2/3*(u2_x(i,j,k)+w2_z(i,j,k))))/(K_average+1.d0)
		 siazz(i,j,k)=(k_average*siazz(i,j,k)+miu(i,j,k)*(4/3*w2_z(i,j,k) &
		                       -2/3*(u2_x(i,j,k)+v2_y(i,j,k))))/(K_average+1.d0)
		 siaxy(i,j,k)=(k_average*siaxy(i,j,k)+miu(i,j,k)*(u2_y(i,j,k)+v2_x(i,j,k)))/(K_average+1.d0)
		 siayz(i,j,k)=(k_average*siayz(i,j,k)+miu(i,j,k)*(v2_z(i,j,k)+w2_y(i,j,k)))/(K_average+1.d0)
		 siaxz(i,j,k)=(k_average*siaxz(i,j,k)+miu(i,j,k)*(u2_z(i,j,k)+w2_x(i,j,k)))/(K_average+1.d0)
                 ! 求σ'第三部分——带Reynolds平均的那部分
		 sig1xx(i,j,k)=miu(i,j,k)*(4/3*(u2_x(i,j,k)+uf_x(i,j,k)) &
		                      -2/3*(v2_y(i,j,k)+w2_z(i,j,k)+vf_y(i,j,k)+wf_z(i,j,k)))-siaxx(i,j,k)
		 sig1yy(i,j,k)=miu(i,j,k)*(4/3*(v2_y(i,j,k)+vf_y(i,j,k)) &
		                      -2/3*(u2_x(i,j,k)+w2_z(i,j,k)+uf_x(i,j,k)+wf_z(i,j,k)))-siayy(i,j,k)
		 sig1zz(i,j,k)=miu(i,j,k)*(4/3*(w2_z(i,j,k)+wf_z(i,j,k)) &
		                      -2/3*(u2_x(i,j,k)+v2_y(i,j,k)+uf_x(i,j,k)+vf_y(i,j,k)))-siazz(i,j,k)
		 sig1xy(i,j,k)=miu(i,j,k)*(u2_y(i,j,k)+v2_x(i,j,k)+uf_y(i,j,k)+vf_x(i,j,k))-siaxy(i,j,k)
		 sig1yz(i,j,k)=miu(i,j,k)*(v2_z(i,j,k)+w2_y(i,j,k)+vf_z(i,j,k)+wf_y(i,j,k))-siayz(i,j,k)
		 sig1xz(i,j,k)=miu(i,j,k)*(u2_z(i,j,k)+w2_x(i,j,k)+uf_z(i,j,k)+wf_x(i,j,k))-siaxz(i,j,k)
		 ! 求σ'



		 u2sxx1(i,j,k)=(k_average*u2sxx1(i,j,k)+u2(i,j,k)*sig1xx(i,j,k))/(K_average+1.d0)
		 u2sxy1(i,j,k)=(k_average*u2sxy1(i,j,k)+u2(i,j,k)*sig1xy(i,j,k))/(K_average+1.d0)
		 u2sxz1(i,j,k)=(k_average*u2sxz1(i,j,k)+u2(i,j,k)*sig1xz(i,j,k))/(K_average+1.d0)
		 v2sxy1(i,j,k)=(k_average*v2sxy1(i,j,k)+v2(i,j,k)*sig1xy(i,j,k))/(K_average+1.d0)
                 v2syy1(i,j,k)=(k_average*v2syy1(i,j,k)+v2(i,j,k)*sig1yy(i,j,k))/(K_average+1.d0)
		 v2syz1(i,j,k)=(k_average*v2syz1(i,j,k)+v2(i,j,k)*sig1yz(i,j,k))/(K_average+1.d0)
		 w2sxz1(i,j,k)=(k_average*w2sxz1(i,j,k)+w2(i,j,k)*sig1xz(i,j,k))/(K_average+1.d0)
		 w2syz1(i,j,k)=(k_average*w2syz1(i,j,k)+w2(i,j,k)*sig1yz(i,j,k))/(K_average+1.d0)
		 w2szz1(i,j,k)=(k_average*w2szz1(i,j,k)+w2(i,j,k)*sig1zz(i,j,k))/(K_average+1.d0)
		 ! 求<U''σ'>
	     enddo
	     enddo
             enddo


		 call dx3d(u2sxx1,u2sxx1_x,nx,ny,nz,hx)		 
                 !call dy3d(u2sxy1,u2sxy1_y,nx,ny,nz,hy)
		 call dy3d(u2sxy1,u2sxy1_y,nx,ny,nz,yy)		 
		 call dz3d(u2sxz1,u2sxz1_z,nx,ny,nz,hz)
		 call dx3d(v2sxy1,v2sxy1_x,nx,ny,nz,hx)
		 !call dy3d(v2syy1,v2syy1_y,nx,ny,nz,hy)
		 call dy3d(v2syy1,v2syy1_y,nx,ny,nz,yy)
		 call dz3d(v2syz1,v2syz1_z,nx,ny,nz,hz)
		 call dx3d(w2sxz1,w2sxz1_x,nx,ny,nz,hx)
		 !call dy3d(w2syz1,w2syz1_y,nx,ny,nz,hy)
		 call dy3d(w2syz1,w2syz1_y,nx,ny,nz,yy)
		 call dz3d(w2szz1,w2szz1_z,nx,ny,nz,hz)
		 ! 求∂<U''σ'>/∂X

 
	     do k=1,nz
             do j=1,ny
             do i=1,nx
		 sxx1u2x(i,j,k)=(k_average*sxx1u2x(i,j,k)+sig1xx(i,j,k)*u2_x(i,j,k))/(K_average+1.d0)
		 sxy1u2y(i,j,k)=(k_average*sxy1u2y(i,j,k)+sig1xy(i,j,k)*u2_y(i,j,k))/(K_average+1.d0)
		 sxz1u2z(i,j,k)=(k_average*sxz1u2z(i,j,k)+sig1xz(i,j,k)*u2_z(i,j,k))/(K_average+1.d0)
	         sxy1v2x(i,j,k)=(k_average*sxy1v2x(i,j,k)+sig1xy(i,j,k)*v2_x(i,j,k))/(K_average+1.d0)
	         syy1v2y(i,j,k)=(k_average*syy1v2y(i,j,k)+sig1yy(i,j,k)*v2_y(i,j,k))/(K_average+1.d0)
		 syz1v2z(i,j,k)=(k_average*syz1v2z(i,j,k)+sig1yz(i,j,k)*v2_z(i,j,k))/(K_average+1.d0)
		 sxz1w2x(i,j,k)=(k_average*sxz1w2x(i,j,k)+sig1xz(i,j,k)*w2_x(i,j,k))/(K_average+1.d0)
		 syz1w2y(i,j,k)=(k_average*syz1w2y(i,j,k)+sig1yz(i,j,k)*w2_y(i,j,k))/(K_average+1.d0)
		 szz1w2z(i,j,k)=(k_average*szz1w2z(i,j,k)+sig1zz(i,j,k)*w2_z(i,j,k))/(K_average+1.d0)
		  
	     enddo
	     enddo
             enddo
                 ! 求<σ' ∂U''/∂X>


 

	     k_average = k_average + 1


	    
    enddo  

! ********************* 循环结束 ************************

   ! 开始组装湍动能方程分解的各个量 
     do k=1,nz
     do j=1,ny
     do i=1,nx

             Pk(i,j,k)=-duu(i,j,k)*uf_x(i,j,k)-duv(i,j,k)*uf_y(i,j,k)-duw(i,j,k)*uf_z(i,j,k) &
		       -duv(i,j,k)*vf_x(i,j,k)-dvv(i,j,k)*vf_y(i,j,k)-dvw(i,j,k)*vf_z(i,j,k) &
	               -duw(i,j,k)*wf_x(i,j,k)-dvw(i,j,k)*wf_y(i,j,k)-dww(i,j,k)*wf_z(i,j,k) 

             Tk(i,j,k)=-duuu_x(i,j,k)-duuv_y(i,j,k)-duuw_z(i,j,k) &
		       -dvvu_x(i,j,k)-dvvv_y(i,j,k)-dvvw_z(i,j,k) &
		       -dwwu_x(i,j,k)-dwwv_y(i,j,k)-dwww_z(i,j,k)

	     Pik1(i,j,k)=-ufr(i,j,k)*pa_x(i,j,k)-vfr(i,j,k)*pa_y(i,j,k)-wfr(i,j,k)*pa_z(i,j,k)
	     Pik2(i,j,k)=-pu_x(i,j,k)-pv_y(i,j,k)-pw_z(i,j,k)
             Pik3(i,j,k)=-p1ux(i,j,k)-p1vy(i,j,k)-p1wz(i,j,k)

             Dk(i,j,k)=u2sxx1_x(i,j,k)+u2sxy1_y(i,j,k)+u2sxz1_z(i,j,k)+  &
		       v2sxy1_x(i,j,k)+v2syy1_y(i,j,k)+v2syz1_z(i,j,k)+  & 
                       w2sxz1_x(i,j,k)+w2syz1_y(i,j,k)+w2szz1_z(i,j,k)
		 
	     Ek(i,j,k)=-sxx1u2x(i,j,k)-sxy1u2y(i,j,k)-sxz1u2z(i,j,k)  &
                       -sxy1v2x(i,j,k)-syy1v2y(i,j,k)-syz1v2z(i,j,k)  &
                       -sxz1w2x(i,j,k)-syz1w2y(i,j,k)-szz1w2z(i,j,k)
		 
	     Mk(i,j,k)=ufr(i,j,k)*sigxx_x(i,j,k) + ufr(i,j,k)*sigxy_y(i,j,k) + ufr(i,j,k)*sigxz_z(i,j,k) + &
		       vfr(i,j,k)*sigxy_x(i,j,k) + vfr(i,j,k)*sigyy_y(i,j,k) + vfr(i,j,k)*sigyz_z(i,j,k) + &
		       wfr(i,j,k)*sigxz_x(i,j,k) + wfr(i,j,k)*sigyz_y(i,j,k) + wfr(i,j,k)*sigzz_z(i,j,k)



     enddo
     enddo
     enddo




!---------make spacial average for TKE------------------------

        Pka=0.d0;Tka=0.d0;Pik1a=0.d0;Pik2a=0.d0;Pik3a=0.d0;Dka=0.d0;Eka=0.d0;Mka=0.d0
	! Favre average
	do j=1,ny  
	  do k=1,nz
          do i=1,nx
	    Pka(j)=Pka(j)+Pk(i,j,k)
	    Tka(j)=Tka(j)+Tk(i,j,k)
	    Pik1a(j)=Pik1a(j)+Pik1(i,j,k)
	    Pik2a(j)=Pik2a(j)+Pik2(i,j,k)
	    Pik3a(j)=Pik3a(j)+Pik3(i,j,k)
	    Dka(j)=Dka(j)+Dk(i,j,k)
	    Eka(j)=Eka(j)+Ek(i,j,k)
	    Mka(j)=Mka(j)+Mk(i,j,k)


	  enddo 
	  enddo

            Pka(j)=Pka(j)/(nx*nz)
	    Tka(j)=Tka(j)/(nx*nz)
	    Pik1a(j)=Pik1a(j)/(nx*nz)
	    Pik2a(j)=Pik2a(j)/(nx*nz)
	    Pik3a(j)=Pik3a(j)/(nx*nz)
	    Dka(j)=Dka(j)/(nx*nz)
	    Eka(j)=Eka(j)/(nx*nz)
	    Mka(j)=Mka(j)/(nx*nz)
	enddo
        !print*, "uf,vf,wf=",uf(ny/2),vf(ny/2),wf(ny/2)

        



!--------- save files ------------------------
	     
	     
	savename='TKE.dat'
	     	     
	open(44,file=savename,form="formatted")
	write(44,*) 'Title = "Comput." '
        write(44,'(A)',advance='no') "variables=<I>y<sup>+</sup></I>,<I>P</I><sub>K</sub><sup>+</sup>,<I>T</I><sub>K</sub><sup>+</sup>,&
	   <I><greek>P</greek></I><sub>K1</sub><sup>+</sup>,<I><greek>P</greek></I><sub>K2</sub><sup>+</sup>,&
	   <I><greek>P</greek></I><sub>K3</sub><sup>+</sup>,<I><greek>P</greek></I><sub>K</sub><sup>+</sup>,<I>D</I><sub>K</sub><sup>+</sup>,&
	   <I><greek>e</greek></I><sub>K</sub><sup>+</sup>,<I>M</I><sub>K</sub><sup>+</sup>"
	write(44,*) 'ZONE T="Comput." '
     	!write(44,*) "zone i= ",nx, " j= ",ny
	    
        do j=2,ny-1
       
	   write(44,"(10f16.8)") yplus(j), Pka(j)/(u_tau**3/delta_niu),Tka(j)/(u_tau**3/delta_niu), &
	   Pik1a(j)/(u_tau**3/delta_niu),Pik2a(j)/(u_tau**3/delta_niu),Pik3a(j)/(u_tau**3/delta_niu), &
	   (Pik1a(j)+Pik2a(j)+Pik3a(j))/(u_tau**3/delta_niu),&
	   Dka(j)/(u_tau**3/delta_niu),Eka(j)/(u_tau**3/delta_niu),Mka(j)/(u_tau**3/delta_niu)
	   
        enddo
       
	    close(44)

             


	    
	    
end
 !=================================================================================================================  
 !                                                Subroutine
 !=================================================================================================================  
      


    ! 求一维空间导数
    subroutine dx3d(f,fx,nx,ny,nz,hx)      ! (待求导 数组,导数数组,网格节点个数*3,平均网格尺寸)
	implicit doubleprecision (a-h,o-z)
    real*8 f(nx,ny,nz),fx(nx,ny,nz)
	    b1=8.d0/(12.d0*hx)
	    b2=1.d0/(12.d0*hx)
        a1=1.d0/(60.d0*hx)
        a2=-3.d0/(20.d0*hx)
        a3=3.d0/(4.d0*hx)

        do k=1,nz
	do j=1,ny

	  do i=4,nx-3
              fx(i,j,k)=a1*(f(i+3,j,k)-f(i-3,j,k)) +a2*(f(i+2,j,k)-f(i-2,j,k)) +a3*(f(i+1,j,k)-f(i-1,j,k))      
	  enddo

              fx(1,j,k)=(-3.d0*f(1,j,k)+4.d0*f(2,j,k)-f(3,j,k))/(2.d0*hx)            !2nd order 
              fx(2,j,k)=(-2.d0*f(1,j,k)-3.d0*f(2,j,k)+6.d0*f(3,j,k)-f(4,j,k))/(6.d0*hx)  !3rd order
	      fx(3,j,k)=b1*(f(4,j,k)-f(2,j,k))-b2*(f(5,j,k)-f(1,j,k))
	      fx(nx-2,j,k)=b1*(f(nx-1,j,k)-f(nx-3,j,k))-b2*(f(nx,j,k)-f(nx-4,j,k))
              fx(nx-1,j,k)=(f(nx-3,j,k)-6.d0*f(nx-2,j,k)+3.d0*f(nx-1,j,k) +2.d0*f(nx,j,k))/(6.d0*hx)
              fx(nx,j,k)=(f(nx-2,j,k)-4.d0*f(nx-1,j,k)+3.d0*f(nx,j,k))/(2.d0*hx)

	enddo
	enddo
        

    end



!    subroutine dy3d(f,fy,nx,ny,nz,hy)      ! (待求导 数组,导数数组,网格节点个数*3,平均网格尺寸)
!	implicit doubleprecision (a-h,o-z)
!    real*8 f(nx,ny,nz),fy(nx,ny,nz)
!	    b1=8.d0/(12.d0*hy)
!	    b2=1.d0/(12.d0*hy)
!        a1=1.d0/(60.d0*hy)
!        a2=-3.d0/(20.d0*hy)
!        a3=3.d0/(4.d0*hy)
!
!        do k=1,nz
!	do i=1,nx
!
!	  do j=4,ny-3
!            fy(i,j,k)=a1*(f(i,j+3,k)-f(i,j-3,k)) +a2*(f(i,j+2,k)-f(i,j-2,k)) +a3*(f(i,j+1,k)-f(i,j-1,k))
!          enddo
!
!            fy(i,1,k)=(-3.d0*f(i,1,k)+4.d0*f(i,2,k)-f(i,3,k))/(2.d0*hy)            !2nd order 
!            fy(i,2,k)=(-2.d0*f(i,1,k)-3.d0*f(i,2,k)+6.d0*f(i,3,k)-f(i,4,k))/(6.d0*hy)  !3rd order
!	    fy(i,3,k)=b1*(f(i,4,k)-f(i,2,k))-b2*(f(i,5,k)-f(i,1,k))
!	    fy(i,ny-2,k)=b1*(f(i,ny-1,k)-f(i,ny-3,k))-b2*(f(i,ny,k)-f(i,ny-4,k))
!            fy(i,ny-1,k)=(f(i,ny-3,k)-6.d0*f(i,ny-2,k)+3.d0*f(i,ny-1,k) +2.d0*f(i,ny,k))/(6.d0*hy)
!            fy(i,ny,k)=(f(i,ny-2,k)-4.d0*f(i,ny-1,k)+3.d0*f(i,ny,k))/(2.d0*hy)
!
!	enddo
!	enddo
!      
!
!    end

    subroutine dy3d(f,fy,nx,ny,nz,yy)      ! (待求导 数组,导数数组,网格节点个数*3,y的坐标)
	implicit doubleprecision (a-h,o-z)
    real*8 f(nx,ny,nz),fy(nx,ny,nz),yy(ny)
    do k=1,nz
    do i=1,nx
       do j=1,ny-1  
          fy(i,j,k)=(f(i,j+1,k)-f(i,j,k))/(yy(j+1)-yy(j))
       enddo
       fy(i,ny,k)=fy(i,ny-1,k)
    enddo
    enddo

    end





    subroutine dz3d(f,fz,nx,ny,nz,hz)      ! (待求导 数组,导数数组,网格节点个数*3,平均网格尺寸)
	implicit doubleprecision (a-h,o-z)
    real*8 f(nx,ny,nz),fz(nx,ny,nz)
	    b1=8.d0/(12.d0*hz)
	    b2=1.d0/(12.d0*hz)
        a1=1.d0/(60.d0*hz)
        a2=-3.d0/(20.d0*hz)
        a3=3.d0/(4.d0*hz)

        do j=1,ny
	do i=1,nx

	  do k=4,nz-3
            fz(i,j,k)=a1*(f(i,j,k+3)-f(i,j,k-3)) +a2*(f(i,j,k+2)-f(i,j,k-2)) +a3*(f(i,j,k+1)-f(i,j,k-1))
          enddo

            fz(i,j,1)=(-3.d0*f(i,j,1)+4.d0*f(i,j,2)-f(i,j,3))/(2.d0*hz)            !2nd order 
            fz(i,j,2)=(-2.d0*f(i,j,1)-3.d0*f(i,j,2)+6.d0*f(i,j,3)-f(i,j,4))/(6.d0*hz)  !3rd order
	    fz(i,j,3)=b1*(f(i,j,4)-f(i,j,2))-b2*(f(i,j,5)-f(i,j,1))
	    fz(i,j,nz-2)=b1*(f(i,j,nz-1)-f(i,j,nz-3))-b2*(f(i,j,nz)-f(i,j,nz-4))
            fz(i,j,nz-1)=(f(i,j,nz-3)-6.d0*f(i,j,nz-2)+3.d0*f(i,j,nz-1) +2.d0*f(i,j,nz))/(6.d0*hz)
            fz(i,j,nz)=(f(i,j,nz-2)-4.d0*f(i,j,nz-1)+3.d0*f(i,j,nz))/(2.d0*hz)

	enddo
	enddo
      

    end
