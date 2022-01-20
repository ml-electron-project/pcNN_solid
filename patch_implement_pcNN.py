############################
# Patch to implement pcNN for VASP ver. 5.4.4. written by Ryo Nagai
# email: r-nag@issp.u-tokyo.ac.jp
#

# This add-on is only available to owners of a valid VASP license.
# This add-on is distributed under an agreement from VASP Software GmbH.

# - This add-on comes without any waranty.
# - Commercial use and patent use are prohibited.
# - If the add-on is modified by inclusion of additional source code lines of the Add-on to the VASP software code,
#   the developer must send the revised version of the add-on to VASP Software GmbH for further review.


# 1. Backup vasp/src/metagga.F.
# 2. Run this patch. You are asked the path to metagga.F.
# 3. metagga_patched.F is generated. Replace metagga.F with it.
# 4. Compile VASP.

# ## Citing this XC functinoal
# @misc{nagai2021machinelearningbased,
#       title={Machine-Learning-Based Exchange-Correlation Functional with Physical Asymptotic Constraints},
#       author={Ryo Nagai and Ryosuke Akashi and Osamu Sugino},
#       year={2021},
#       eprint={2111.15593},
#       archivePrefix={arXiv},
#       primaryClass={cond-mat.mtrl-sci}
# }
############################

def patch():
    print("input path to metagga.F:")
    target_file = input()
    output = "metagga_patched.F"
    with open(target_file, "r") as f:
        l_strip = [s.strip("\n") for s in f.readlines()]

        p1s = [p1]
        p2s = [p2]
        p3s = [p3]
        p4s = [p4]
        p5s = [p5]
        p6s = [p6]
        for i, li in enumerate(l_strip):
            if "MODULE setxcmeta" in li:
                l_strip = l_strip[:i]+p1s+l_strip[i:]
                break
        for i, li in enumerate(l_strip):
            if "=='PBE') THEN" in li:
                l_strip = l_strip[:i+4]+p2s+l_strip[i+4:]
                break
        for i, li in enumerate(l_strip):
            if "! PBE, for testing mainly" in li:
                l_strip = l_strip[:i+3]+p3s+l_strip[i+3:]
                break

        for i, li in enumerate(l_strip):
            if "END MODULE metalib" in li:
                l_strip = l_strip[:i-3]+p4s+l_strip[i-3:]
                break
        for i, li in enumerate(l_strip):
            if "SUBROUTINE XC_META_" in li:
                l_strip = l_strip[:i+1]+p5s+l_strip[i+1:]
                break
        for i, li in enumerate(l_strip):
            if "SUBROUTINE METAGGASPIN(&" in li:
                l_strip = l_strip[:i+22]+p6s+l_strip[i+22:]
                break

        with open(output, "w") as fo:
            fo.write('\n'.join(l_strip))


p1 = """  module readNN
      integer, PARAMETER :: hidden=100
      real*8,save,  dimension(:,:) :: w1(hidden,2), w2(hidden,hidden), w3(hidden,hidden), w4(1,hidden)
      real*8,save,  dimension(:,:) :: w1c(hidden,2), w2c(hidden,hidden), w3c(hidden,hidden*2), w4c(1,hidden)
      real*8,save,  dimension(:):: b1(hidden), b2(hidden), b3(hidden),b4(1)
      real*8,save,  dimension(:):: b1c(hidden), b2c(hidden), b3c(hidden),b4c(1)
      contains
      subroutine loadNN(w1, w2, w3, w4, b1, b2, b3, b4, &
        & w1c, w2c, w3c, w4c, b1c, b2c, b3c, b4c)
        real*8,  dimension(:,:) :: w1, w2, w3, w4
        real*8,  dimension(:,:) :: w1c, w2c, w3c, w4c
        real*8  :: b1(:), b2(:), b3(:), b4(:)
        real*8  :: b1c(:), b2c(:), b3c(:), b4c(:)
        integer :: n, m, i
        
        n=hidden
        m=2
        open(80, file='nnparams/w1.txt')
        do i =1,n
          read(80,*) (w1(i,j), j=1,m)
        end do
        close(80)
    
        n=hidden
        m=hidden
        open(81, file='nnparams/w2.txt')
        do i =1,n
          read(81,*) (w2(i,j), j=1,m)
        end do
        close(81)
    
        n=hidden
        m=hidden
        open(82, file='nnparams/w3.txt')
        do i =1,n
          read(82,*) (w3(i,j), j=1,m)
        end do
        close(82)
    
        n=1
        m=hidden
        open(83, file='nnparams/w4.txt')
        do i =1,n
          read(83,*) (w4(i,j), j=1,m)
        end do
        close(83)
    
        n=hidden
        open(84, file='nnparams/b1.txt')
        do i =1,n
          read(84,*) b1(i)
        end do
        close(84)
    
        n=hidden
        open(85, file='nnparams/b2.txt')
        do i =1,n
          read(85,*) b2(i)
        end do
        close(85)
    
        n=hidden
        open(86, file='nnparams/b3.txt')
        do i =1,n
          read(86,*) b3(i)
        end do
        close(86)
    
        n=1
        open(87, file='nnparams/b4.txt')
        do i =1,n
          read(87,*) b4(i)
        end do
        close(87)

        n=hidden
        m=2
        open(88, file='nnparams/w1c.txt')
        do i =1,n
          read(88,*) (w1c(i,j), j=1,m)
        end do
        close(88)
    
        n=hidden
        m=hidden
        open(89, file='nnparams/w2c.txt')
        do i =1,n
          read(89,*) (w2c(i,j), j=1,m)
        end do
        close(89)
    
        n=hidden
        m=hidden*2
        open(90, file='nnparams/w3c.txt')
        do i =1,n
          read(90,*) (w3c(i,j), j=1,m)
        end do
        close(90)
    
        n=1
        m=hidden
        open(91, file='nnparams/w4c.txt')
        do i =1,n
          read(91,*) (w4c(i,j), j=1,m)
        end do
        close(91)
    
        n=hidden
        open(92, file='nnparams/b1c.txt')
        do i =1,n
          read(92,*) b1c(i)
        end do
        close(92)
    
        n=hidden
        open(93, file='nnparams/b2c.txt')
        do i =1,n
          read(93,*) b2c(i)
        end do
        close(93)
    
        n=hidden
        open(94, file='nnparams/b3c.txt')
        do i =1,n
          read(94,*) b3c(i)
        end do
        close(94)
    
        n=1
        open(95, file='nnparams/b4c.txt')
        do i =1,n
          read(95,*) b4c(i)
        end do
        close(95)
      end subroutine loadNN

  end module readNN
"""

p2 = """      ELSEIF (SZNAM(1:6)=='NNMGGA') THEN
          ID_METAGGA=777
          LMETA_NEEDS_POT=.TRUE.
          LMETA_NEEDS_MU=.TRUE.
          call loadNN(w1, w2, w3, w4, b1, b2, b3, b4, &
          &w1c, w2c, w3c, w4c, b1c, b2c, b3c, b4c)
          """

p3 = """      ELSEIF (ID_METAGGA==777) THEN
        ! NN-metaGGA
          CALL NNMGGA_XC(&
           RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW,&
           Exc_NNM,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)
           ! Sum everything
           EXC=Exc_NNM/(RHOUP+RHODW)
           dEXCdRHOup=VXD1
           dEXCdRHOdw=VXD2
           dEXCdABSNABup=VXDD1
           dEXCdABSNABdw=VXDD2
           dEXCdTAUup=AMUXD1
           dEXCdTAUdw=AMUXD2

         !  write(*,*) RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
         !& EXC,dEXCdRHOup,dEXCdRHOdw,dEXCdABSNABup,dEXCdABSNABdw,dEXCdTAUup,dEXCdTAUdw
           ! from Hartree to Rydberg
           EXC=EXC*2
           dEXCdRHOup=dEXCdRHOup*2
           dEXCdRHOdw=dEXCdRHOdw*2
           dEXCdABSNABup=dEXCdABSNABup*2
           dEXCdABSNABdw=dEXCdABSNABdw*2
           dEXCdTAUup=dEXCdTAUup*2
           dEXCdTAUdw=dEXCdTAUdw*2
"""
p4 = """
      subroutine NNMGGA_XC(&
        RU,RD,DRU,DRD,DRT,TAUU,TAUD,&
        Exc_NNM,VXCD1,VXCDD1,VXCD2,VXCDD2,AMUXCD1,AMUXCD2)
        use :: readNN
        IMPLICIT None
        integer :: n,m,i,j
        REAL(q) Exc_NNM, RU,RD,DRU,DRD,DRT,TAUU,TAUD,RHO
        REAL(q) exc, fxc, fx_u, fx_d, fc, DRUDcos
        REAL(q) VXCD1,VXCD2,VXCDD1,VXCDD2,AMUXCD1,AMUXCD2
        REAL(q) Ex_SCAN1,VXD11,VXDD11,VXD12,VXDD12,AMUXD11,AMUXD12
        REAL(q) Ex_SCAN2,VXD21,VXDD21,VXD22,VXDD22,AMUXD21,AMUXD22
        REAL(q) Ec_SCAN,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2
        real*8, dimension(:) :: x(7)
        real*8, dimension(:,:) :: dfx_dnu(1,4), dfx_dnd(1,4), dfc_dn(1,4)
        real*8, dimension(:) :: t(4), logt(4), g1(hidden), h1(hidden), g2(hidden), h2(hidden), g3(hidden), h3(hidden), g4(1), h4(1), eunif_x(7)
        real*8, dimension(:) :: fxc_g4(1), fxc_g3(hidden), fxc_g2(hidden), fxc_g1(hidden), fxc_h3(hidden),fxc_h2(hidden), fxc_h1(hidden)
        real*8, dimension(:) :: fxc_logt(4), fxc_t(4), fxc_logt_t(4), t0_x(7), t1_x(7), t2_x(7), t3_x(7), fxc_t_x(7), fxc_x(7)
        real*8, dimension(:,:):: t_x(4,7)
        real*8, PARAMETER :: THRD=1./3.,THRD4=4./3.,THRD5=5./3.
        real*8, PARAMETER :: a = 4.0, sles=0.2
        real*8, dimension(:) :: backlogt(4)
      ! SCAN
         ! Exchange
        CALL VSCANx(&
        &   RU,RU,DRU,DRU,DRU*2,TAUU,TAUU, &
        &   Ex_SCAN1,VXD11,VXDD11,VXD12,VXDD12,AMUXD11,AMUXD12)
        Ex_SCAN1=Ex_SCAN1/(RU+RU)
        CALL VSCANx(&
        &   RD,RD,DRD,DRD,DRD*2,TAUD,TAUD, &
        &   Ex_SCAN2,VXD21,VXDD21,VXD22,VXDD22,AMUXD21,AMUXD22)
        Ex_SCAN2=Ex_SCAN2/(RD+RD)
          ! Correlation
          CALL VSCANc(&
        &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
        &   Ec_SCAN,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2)
        Ec_SCAN=Ec_SCAN/(RU+RD)
          ! Sum everything
          ! EXC_SCAN=(Ex_SCAN+Ec_SCAN)/(RU+RD)
          ! dEXCdRHOup_SCAN=VXD1+VCD1
          ! dEXCdRHOdw_SCAN=VXD2+VCD2
          ! dEXCdABSNABup_SCAN=VXDD1+VCDD1
          ! dEXCdABSNABdw_SCAN=VXDD2+VCDD2
          ! dEXCdTAUup_SCAN=AMUXD1+AMUCD1
          ! dEXCdTAUdw_SCAN=AMUXD2+AMUCD2
          !write(*,*) RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
          !&Exc_SCAN,dEXCdRHOup_SCAN,dEXCdRHOdw_SCAN,dEXCdABSNABup_SCAN,dEXCdABSNABdw_SCAN,dEXCdTAUup_SCAN,dEXCdTAUdw_SCAN

        call calc_x(RU*2,DRU*2,TAUU*2, fx_u, dfx_dnu)
        call calc_x(RD*2,DRD*2,TAUD*2, fx_d, dfx_dnd)
        
        call calc_c(RU, RD, DRT, TAUU+TAUD, fc, dfc_dn)
        


        exc=0.5*(2.*RU*Ex_SCAN1*fx_u+2.*RD*Ex_SCAN2*fx_d)/(RU+RD)
        exc=exc+fc*Ec_SCAN
        Exc_NNM=exc*(RU+RD)

        
        VXCD1=0.5*(VXD11+VXD12)*fx_u+2*RU*Ex_SCAN1*dfx_dnu(1,1)
        VXCD1=VXCD1+VCD1*fc+(RU+RD)*Ec_SCAN*dfc_dn(1,1)

        VXCD2=0.5*(VXD21+VXD22)*fx_d+2*RD*Ex_SCAN2*dfx_dnd(1,1)
        VXCD2=VXCD2+VCD2*fc+(RU+RD)*Ec_SCAN*dfc_dn(1,2)

        DRUDcos=0.5*(DRT**2.-DRU**2.-DRD**2.)

        !TODO 下の式で正しいがなぜ合うのかわからない(exchange)
        VXCDD1=0.5*(VXDD11*2.)*fx_u+RU*Ex_SCAN1*(dfx_dnu(1,3)*2.)
        VXCDD1=VXCDD1+VCDD1*fc+(RU+RD)*Ec_SCAN*(dfc_dn(1,3)*0.5/DRT*(2.*DRU+2*DRUDcos/DRU))

        VXCDD2=0.5*(VXDD21*2.)*fx_d+RD*Ex_SCAN2*(dfx_dnd(1,3)*2.)
        VXCDD2=VXCDD2+VCDD2*fc+(RU+RD)*Ec_SCAN*(dfc_dn(1,3)*0.5/DRT*(2.*DRD+2*DRUDcos/DRD))

        AMUXCD1=0.5*(AMUXD11+AMUXD12)*fx_u+2*RU*Ex_SCAN1*dfx_dnu(1,4)
        AMUXCD1=AMUXCD1+AMUCD1*fc+(RU+RD)*Ec_SCAN*dfc_dn(1,4)

        AMUXCD2=0.5*(AMUXD21+AMUXD22)*fx_d+2*RD*Ex_SCAN2*dfx_dnd(1,4)
        AMUXCD2=AMUXCD2+AMUCD2*fc+(RU+RD)*Ec_SCAN*dfc_dn(1,4)
            
      end subroutine NNMGGA_XC
    

      subroutine calc_x(n, dn, tau, f, df_dn)
        use :: readNN
        implicit none
        real*8, PARAMETER :: PI =3.141592653589793238
        integer i,j,k,l
        integer, PARAMETER :: nconds=2, nt=2, nsize=4 
        real*8, PARAMETER :: THRD=1./3.,THRD2=2./3.,THRD4=4./3.,THRD5=5./3., delta_x=1.
    
        real*8 :: n, dn, tau
        real*8 :: s, tmpval, csum, unif, iunif, tauunif, itauunif, dunif_dn, ft
        real*8, dimension(:) :: t(nt), f0(nconds), cs(nconds)
        real*8, dimension(:,:):: t0(nconds, nt), f_t0(nconds)
        real*8, dimension(:,:):: dc_dt(nconds,nt)
        real*8, dimension(:,:):: dft0_dt0(nconds,nt), dft_dt(1,nt)
        
        real*8, dimension(:):: df0_dt(nconds,nt)
        real*8, dimension(:)::dt_dn(nt, nsize)
        real*8, dimension(:,:,:):: dc_dt0(nconds, nconds, nt), dt0_dt(nconds, nt, nt)
        
        real*8:: f
        real*8, dimension(:,:):: df_dn(1, nsize)
        
        df0_dt=0.
        unif=(n+1.e-7)**THRD
        iunif=1.0/unif
        s=(dn**2.0+0.1**(56.0/3.0))**0.5/unif**4.
        t(1)=tanh(s/1.0)
        tauunif=2.871234000188191*unif**5.
        itauunif=1.0/tauunif
        t(2)=tanh((tau*itauunif-1.0)/1.0)
    
        dunif_dn=THRD*(n+1.e-7)**(-THRD2)!correct
    
        tmpval=(1.-t(1))*(1.+t(1))
        dt_dn(1,1)=tmpval*(-4.)*s*iunif*dunif_dn
        dt_dn(1,2)=0
        dt_dn(1,3)=tmpval*0.5*s/(dn**2.0+0.1**(56.0/3.0))*2.0*dn
        dt_dn(1,4)=0
    
        tmpval=(1.-t(2))*(1.+t(2))
        dt_dn(2,1)=tmpval*(-5.)*tau*itauunif*iunif*dunif_dn
        dt_dn(2,2)=0
        dt_dn(2,3)=0
        dt_dn(2,4)=tmpval*1.0/tauunif
    
    
        call NN_x(t, ft, dft_dt)
    
        !uniform electron gas
        t0(1,1)=0
        t0(1,2)=0
        f0(1)=1
    
        !f vanish like s^(-1/2) as s goes to inf
        t0(2,1)=1
        t0(2,2)=t(2)
        f0(2)=1
    
        dt0_dt=0.
        dt0_dt(2,2,2)=1.
        
        do i =1,2
          call NN_x(t0(i,:), f_t0(i), dft0_dt0(i,:))
        end do 
        
        call construct_c(t, t0, cs, dc_dt, dc_dt0, delta_x)
        
        csum=0.
        f=0.
        do i =1,2
          f=f+(ft-f_t0(i)+f0(i))*cs(i)
          csum=csum+cs(i)
        end do 
        f=f/csum
    
        !!!!!!backprop!!!!!!
        call backpropagation(nconds, nt, f, ft, f_t0, f0, cs, csum, dft_dt, dft0_dt0, &
        &df0_dt, dc_dt, dc_dt0, dt0_dt, dt_dn, df_dn)
        call shifted_softplus1(f, df_dn)
        !correct: dft_dn, dt0_dn
        ! write(*,*) f
        ! write(*,*) df_dn(1,1),df_dn(1,2),df_dn(1,3),df_dn(1,4) !TODO 値が合わないので修正
        !memo: df_dt0_dt と df_dft_dt がほぼ同じ値なため桁落ちが起きる
        ! write(*,*) (df_dt0_dt(1,i), i=1,2)
        ! write(*,*) (dft0_dt(2,i), i=1,2)#これはあう
        ! write(*,*) dft_dt(1,1)*csum/cs(2),dft_dn(1,2)*csum/cs(2),dft_dn(1,3)*csum/cs(2),dft_dn(1,4)*csum/cs(2)
        !write(*,*) dft0_dt0(1,1),dft0_dt0(1,2),dft0_dt0(2,1),dft0_dt0(2,2)!合格
        ! write(*,*) f_t0(1,1)
    
        !write(*,*) dt0_dt(2,2,2)*dt_dn(2,1), dt0_dt(2,2,2)*dt_dn(2,2), dt0_dt(2,2,2)*dt_dn(2,3), dt0_dt(2,2,2)*dt_dn(2,4)
        !write(*,*) dt_dn(1,1),dt_dn(1,2),dt_dn(1,3),dt_dn(1,4)
        ! write(*,*) dt_dn(2,1),dt_dn(2,2),dt_dn(2,3),dt_dn(2,4)
      end subroutine calc_x
    
      subroutine calc_c(nup, ndw, dn, tau, fc, dfc_dn)
        use :: readNN
        implicit none
        real*8, PARAMETER :: PI =3.141592653589793238
        integer i,j,k,l
        integer, PARAMETER :: nconds=3, nt=4, nsize=4 
        real*8, PARAMETER :: THRD=1./3.,THRD2=2./3.,THRD4=4./3.,THRD5=5./3., delta_c=1.
    
        real*8 :: nup, ndw, n, dn, tau, div, ft
        real*8 :: s, tmpval, tmpval2, csum, unif, iunif, tauunif, itauunif, dunif_dn, tmp1, tmp2
        real*8, dimension(:) :: t(nt), f0(nconds), cs(nconds), t_low(nt), t_low0(nt), tmp3(1,nt), tmp4(1,nt)
        real*8 :: f_low, f_low0
        real*8, dimension(:,:):: t0(nconds, nt), f_t0(nconds)
        real*8, dimension(:,:):: dc_dt(nconds,nt), df_dft(1,1)
        real*8, dimension(:,:):: dft0_dt0(nconds,nt), dft_dt(1,nt)
        
        real*8, dimension(:):: df0_dt(nconds,nt)
        real*8, dimension(:):: dfdt_low(1,nt), dfdt_low0(1,nt)
        real*8, dimension(:)::dt_dn(nt, nsize), dt_low_dt(nt, nt), dt_low0_dt(nt, nt)
        real*8, dimension(:,:,:):: dc_dt0(nconds, nconds, nt), dt0_dt(nconds, nt, nt)
        
        real*8:: fc
        real*8, dimension(:,:):: dfc_dn(1, nsize)
        real*8, dimension(:,:):: dc_dn(nconds,nt)
        df0_dt=0.
        dt_low_dt=0.
        dt_low0_dt=0.
        dt0_dt=0.
    
        n=nup+ndw
        unif=(n+1.e-7)**THRD
        iunif=1.0/unif
        t(1)=tanh(unif/1.0)
    
        div=1./(n+1.e-7)
        t(2)=(1+(nup-ndw)*div)**THRD4+(1-(nup-ndw)*div)**THRD4
        t(2)=tanh(t(2)*0.5/1.0)
    
        s=(dn**2.0+0.1**(56.0/3.0))**0.5/unif**4.
        t(3)=tanh(s/1.0)
    
        tauunif=5.921762640653615*unif**5.
        itauunif=1.0/tauunif
        t(4)=tanh((tau*itauunif-1.0)/1.0)
    
    
    
        dunif_dn=THRD*(n+1.e-7)**(-THRD2)!correct
    
        tmpval=(1.-t(1))*(1.+t(1))
        dt_dn(1,1)=tmpval*dunif_dn
        dt_dn(1,2)=tmpval*dunif_dn
        dt_dn(1,3)=0
        dt_dn(1,4)=0
        
        tmpval=(1.-t(2))*(1.+t(2))
        tmpval2=div*div
        tmp1=THRD4*(1+(nup-ndw)*div)**THRD
        tmp2=THRD4*(1-(nup-ndw)*div)**THRD
        dt_dn(2,1)=tmp1*((n+1.e-7)*1.-1.*(nup-ndw))*tmpval2+tmp2*((n+1.e-7)*(-1.)-1.*(ndw-nup))*tmpval2
        dt_dn(2,1)=tmpval*dt_dn(2,1)*0.5
        dt_dn(2,2)=tmp1*((n+1.e-7)*(-1.)-1.*(nup-ndw))*tmpval2+tmp2*((n+1.e-7)*1.-1.*(ndw-nup))*tmpval2
        dt_dn(2,2)=tmpval*dt_dn(2,2)*0.5
        dt_dn(2,3)=0
        dt_dn(2,4)=0
        
        tmpval=(1.-t(3))*(1.+t(3))
        tmpval2=tmpval*(-4.)*s*iunif*dunif_dn
        dt_dn(3,1)=tmpval2
        dt_dn(3,2)=tmpval2
        dt_dn(3,3)=tmpval*0.5*s/(dn**2.0+0.1**(56.0/3.0))*2.0*dn
        dt_dn(3,4)=0
    
        tmpval=(1.-t(4))*(1.+t(4))
        tmpval2=tmpval*(-5.)*tau*itauunif*iunif*dunif_dn
        dt_dn(4,1)=tmpval2
        dt_dn(4,2)=tmpval2
        dt_dn(4,3)=0
        dt_dn(4,4)=tmpval*1.0*itauunif
    
        call NN_c(t, ft, dft_dt)
    
        !uniform electron gas
        t0(1,1)=t(1)
        t0(1,2)=t(2)
        t0(1,3)=0
        t0(1,4)=0
        f0(1)=1.
    
        !f vanish like s^(-1/2) as s goes to inf
        t0(2,1)=0
        t0(2,2)=t(2)
        t0(2,3)=t(3)
        t0(2,4)=t(4)
    
        t_low(1)=t(1)
        t_low(2)=0.
        t_low(3)=t(3)
        t_low(4)=t(4)
    
        dt_low_dt(1,1)=1.
        dt_low_dt(3,3)=1.
        dt_low_dt(4,4)=1.
    
        t_low0(1)=0.
        t_low0(2)=0.
        t_low0(3)=t(3)
        t_low0(4)=t(4)
        dt_low0_dt(3,3)=1.
        dt_low0_dt(4,4)=1.
        call NN_c(t_low, f_low, dfdt_low(1,:))
        call NN_c(t_low0, f_low0, dfdt_low0(1,:))
        
        f0(2)=f_low-f_low0+1.
        tmp3=matmul(dfdt_low, dt_low_dt)
        tmp4=matmul(dfdt_low0, dt_low0_dt)
        do i=1,4
          df0_dt(2,i)=tmp3(1,i)-tmp4(1,i)
        end do
    
        t0(3,1)=1.
        t0(3,2)=t(2)
        t0(3,3)=t(3)
        t0(3,4)=t(4)
        f0(3)=1.
        
        dt0_dt(1,1,1)=1.
        dt0_dt(1,2,2)=1.
    
        dt0_dt(2,2,2)=1.
        dt0_dt(2,3,3)=1.
        dt0_dt(2,4,4)=1.
    
        dt0_dt(3,2,2)=1.
        dt0_dt(3,3,3)=1.
        dt0_dt(3,4,4)=1.
    
        do i =1,nconds
          call NN_c(t0(i,:), f_t0(i), dft0_dt0(i,:))
        end do 
        
        call construct_c(t, t0, cs, dc_dt, dc_dt0, delta_c)
    
        fc=0.
        csum=0.
        do i =1,nconds
            fc=fc+(ft-f_t0(i)+f0(i))*cs(i)
            csum=csum+cs(i)
        end do 
        
        fc=fc/csum

        call backpropagation(nconds, nt, fc, ft, f_t0, f0, cs, csum, dft_dt, dft0_dt0, &
        &df0_dt, dc_dt, dc_dt0, dt0_dt, dt_dn, dfc_dn)
        call shifted_softplus1(fc, dfc_dn)
    
      end subroutine calc_c
    
      function vmul(x,y)
        real*8 :: x(:),y(:)
        real*8 :: vmul(size(x))
        do i=1, size(y)
          vmul(i)=x(i)*y(i)
        end do
      end function vmul
    
      function vadd(x,y)
        real*8 :: x(:),y(:)
        real*8 :: vadd(size(x))
        do i=1, size(y)
          vadd(i)=x(i)+y(i)
        end do
      end function vadd 
    
      function shifted_softplus0(x)
        real*8, PARAMETER :: l2 =0.6931471805599453
        real*8 :: x(:)
        real*8 :: shifted_softplus0(size(x))
        do i=1, size(x)
          shifted_softplus0(i)=1.0/l2*log(1.0+exp(2.0*l2*x(i)))
        end do
      end function shifted_softplus0
    
      function back_shifted_softplus0(x)
        real*8, PARAMETER :: l2 =0.6931471805599453
        real*8 :: x(:), tmp
        real*8 :: back_shifted_softplus0(size(x))
        do i=1, size(x)
          tmp=exp(2.0*l2*x(i))
          back_shifted_softplus0(i)=2.0*tmp/(1.0+tmp)
        end do
      end function back_shifted_softplus0
    
      subroutine shifted_softplus1(f, df_dn)
        implicit none
        real*8, PARAMETER :: l2 =0.6931471805599453
        real*8 :: f, df_dn(:,:), tmp, tmp2
        integer :: i
    
        tmp=exp(2.0*l2*(f-1.0))
        tmp2=2.0*tmp/(1.0+tmp)
    
        f=1.0/l2*log(1.0+tmp)
        do i=1, size(df_dn(1,:))
          df_dn(1,i)=tmp2*df_dn(1,i)
        end do
      end subroutine shifted_softplus1
    
      
    
      subroutine NN_x(t, ft, dft_dt)
        use :: readNN
        implicit none
        real*8, dimension(:) :: t(2)
        real*8, dimension(:,:) :: dft_dt(1,2)
        real*8 ::ft
        real*8, dimension(:) :: g1(hidden), h1(hidden), g2(hidden), h2(hidden), g3(hidden), h3(hidden), g4(1)
        real*8, dimension(:):: dft_dh3(hidden), dft_dg3(hidden), dft_dh2(hidden),dft_dg2(hidden),dft_dh1(hidden),dft_dg1(hidden)
    
        g1=vadd(matmul(t, transpose(w1)),b1)
        h1=shifted_softplus0(g1)
        g2=vadd(matmul(h1, transpose(w2)),b2)
        h2=shifted_softplus0(g2)
        g3=vadd(matmul(h2, transpose(w3)),b3)
        h3=shifted_softplus0(g3)
        g4=vadd(matmul(h3, transpose(w4)),b4)
        ft=g4(1)
    
        dft_dh3=w4(1,:) !matmul(dft_dg4,w4)
        dft_dg3=vmul(back_shifted_softplus0(g3),dft_dh3)
        dft_dh2=matmul(dft_dg3,w3)
        dft_dg2=vmul(back_shifted_softplus0(g2),dft_dh2)
        dft_dh1=matmul(dft_dg2,w2)
        dft_dg1=vmul(back_shifted_softplus0(g1),dft_dh1)
        dft_dt(1,:)=matmul(dft_dg1,w1)
    
      end subroutine NN_x
    
      subroutine NN_c(t, ft, dft_dt)
        use :: readNN
        implicit none
        real*8, dimension(:) :: t(4)
        real*8, dimension(:,:) :: dft_dt(1,4)
        real*8 ::ft
        real*8, dimension(:) :: g11(hidden), h11(hidden), g21(hidden), g12(hidden), h12(hidden)
        real*8, dimension(:) :: g22(hidden), h21(hidden), h22(hidden), h2(hidden*2), g3(hidden), h3(hidden), g4(1)
        real*8, dimension(:):: dft_dh3(hidden), dft_dg3(hidden), dft_dh2(hidden*2),dft_dg2(hidden),dft_dh1(hidden),dft_dg1(hidden)
        
        g11=vadd(matmul(t(1:2), transpose(w1c)),b1c)
        h11=shifted_softplus0(g11)
        g21=vadd(matmul(h11, transpose(w2c)),b2c)
        h21=shifted_softplus0(g21)
    
        g12=vadd(matmul(t(3:4), transpose(w1)),b1)
        h12=shifted_softplus0(g12)
        g22=vadd(matmul(h12, transpose(w2)),b2)
        h22=shifted_softplus0(g22)
    
        h2(1:100)=h21
        h2(101:200)=h22
    
        g3=vadd(matmul(h2, transpose(w3c)),b3c)
        h3=shifted_softplus0(g3)
        g4=vadd(matmul(h3, transpose(w4c)),b4c)
        ft=g4(1)
    
        ! if (back==1) then 
        dft_dh3=w4c(1,:) !matmul(dft_dg4,w4)
        dft_dg3=vmul(back_shifted_softplus0(g3),dft_dh3)
        dft_dh2=matmul(dft_dg3,w3c)
        
        dft_dg2=vmul(back_shifted_softplus0(g21),dft_dh2(1:100))
        dft_dh1=matmul(dft_dg2,w2c)
        dft_dg1=vmul(back_shifted_softplus0(g11),dft_dh1)
        dft_dt(1,1:2)=matmul(dft_dg1,w1c)
    
        dft_dg2=vmul(back_shifted_softplus0(g22),dft_dh2(101:200))
        dft_dh1=matmul(dft_dg2,w2)
        dft_dg1=vmul(back_shifted_softplus0(g12),dft_dh1)
        dft_dt(1,3:4)=matmul(dft_dg1,w1)
        ! end if 
      end subroutine NN_c
    
      subroutine construct_c(t, t0, cs, dc_dt, dc_dt0, delta)
        implicit none
        real*8 :: t(:), t0(:,:), cs(:), dc_dt(:,:), dc_dt0(:,:,:), delta
        real*8 :: delta2, tmpval, ci_dj, denomi2
        integer :: i,j,k,l, st, st0, i2, j2
        real*8, dimension(:) :: dis(size(t0(:,1))), numer(size(t0(:,1))), denomi(size(t0(:,1)))
        real*8, dimension(:,:) :: dis_ij(size(t0(:,1)),size(t0(:,1)))
        
        st=size(t)
        st0=size(t0(:,1))
        
        delta2=delta**2.
        do i =1,st0
          dis(i)=0.
          do j =1,st
            dis(i)=dis(i)+(t(j)-t0(i,j))**2.
          end do
          dis(i)=tanh(dis(i)/delta2)
        end do 
        
        do i =1,st0
          do j =i+1,st0
            dis_ij(i,j)=0.
            do k =1,st
              dis_ij(i,j)=dis_ij(i,j)+(t0(i,k)-t0(j,k))**2.
            end do 
            dis_ij(i,j)=tanh(dis_ij(i,j)/delta2)
            dis_ij(j,i)=dis_ij(i,j)
          end do 
        end do 
        
        do i =1,st0
          denomi(i)=1.
          numer(i)=1.
          do j =1,st0
            if (j==i) cycle
            denomi(i)=denomi(i)*dis(j)
            numer(i)=numer(i)*dis_ij(i,j)
          end do 
          cs(i)=denomi(i)/numer(i)
        end do 
        
        do i =1,st0
          do j =1,st0
    
            if (i==j) then 
              do k=1,st
                tmpval=0.0
                do l=1,st0
                  if (l==i) cycle
                  tmpval=tmpval-2.*(t0(i,k)-t0(l,k))*(1.-dis_ij(i,l))*(1.+dis_ij(i,l))/dis_ij(i,l)
                end do
                dc_dt0(i,j,k)=cs(i)*tmpval/delta2
              end do
            else 
              do k=1,st
                denomi2=1.
                do j2 =1,st0
                  if (j2==i .or. j2==j) cycle
                  denomi2=denomi2*dis(j2)
                end do 
                ci_dj=denomi2/numer(i)
                !ci_dj=c(i)/dis(j) implemented to avoid zero division
    
                tmpval=-2.*(t0(j,k)-t0(i,k))*(1.-dis_ij(i,j))*(1.+dis_ij(i,j))*cs(i)/dis_ij(i,j)&
                &+2.*(t0(j,k)-t(k))*(1.-dis(j))*(1.+dis(j))*ci_dj
    
                dc_dt0(i,j,k)=tmpval/delta2
              end do
            end if
          end do 
        end do
    
        do i =1,st0
          do j =1,st
            tmpval=0
            do k=1,st0
              if (i==k) cycle
    
              denomi2=1.
              do j2 =1,st0
                if (j2==i .or. j2==k) cycle
                denomi2=denomi2*dis(j2)
              end do 
              ci_dj=denomi2/numer(i)
              !ci_dj=c(i)/dis(k) implemented to avoid zero division
    
              tmpval=tmpval+2.*ci_dj*(1.-dis(k))*&
              &(1.+dis(k))*(t(j)-t0(k,j))
            end do
            dc_dt(i,j)=tmpval/delta2
          end do 
        end do
      end subroutine construct_c
    
      subroutine backpropagation(nconds, nt, f, ft, f_t0, f0, cs, csum, dft_dt, dft0_dt0, &
        &df0_dt, dc_dt, dc_dt0, dt0_dt, dt_dn, df_dn)
        implicit none
        integer i,j,k,l
        integer:: nconds, nt
        real*8 :: f, ft, f_t0(:), f0(:), csum, cs(:), dft_dt(:,:), dft0_dt0(:,:), df0_dt(:,:), dc_dt(:,:), dc_dt0(:,:,:)
        real*8 :: dt_dn(:,:), df_dn(:,:), dt0_dt(:,:,:)
    
        real*8, PARAMETER :: THRD=1./3.,THRD2=2./3.,THRD4=4./3.,THRD5=5./3., delta_c=1.
        real*8 :: df_dft(1,1), df_dft0(1,nconds),df_df0(1,nconds),df_dc(1,nconds)
        real*8 :: df_df0_dt(1,nt), df_dt0(1,nconds,nt), df_dft_dt(1,nt), df_dc_dt(1,nt)
        real*8 :: df_dt0_dt(1,nt), df_dt(1,nt)
    
        df_dft(1,1)=1.
        
        do i =1,nconds
          df_dft0(1,i)=-cs(i)/csum
          df_df0(1,i)=cs(i)/csum 
          df_dc(1,i)=(ft-f_t0(i)+f0(i))/csum-f/csum
        end do 
    
        df_df0_dt=matmul(df_df0, df0_dt)
    
        do j=1,nt
          do i=1,nconds
            df_dt0(1,i,j)=df_dft0(1,i)*dft0_dt0(i,j)
          end do
        end do
    
        df_dft_dt=matmul(df_dft, dft_dt)
    
    
        do i=1,nconds
          do k=1,nconds
            do j=1,nt
              do l=1,nt
                dc_dt(i,j)=dc_dt(i,j)+dc_dt0(i,k,l)*dt0_dt(k,l,j)
              end do
            end do
          end do
        end do
    
        df_dc_dt=matmul(df_dc,dc_dt)
    
        df_dt0_dt=0.
        do i=1,nt
          do k=1,nt
            do j=1,nconds
              df_dt0_dt(1,i)=df_dt0_dt(1,i)+df_dt0(1,j,k)*dt0_dt(j,k,i)
            end do
          end do
        end do
        
        do i=1,nt
            df_dt(1,i)=df_dt0_dt(1,i)+df_dft_dt(1,i)+df_dc_dt(1,i)+ df_df0_dt(1,i)
        end do
    
        df_dn=matmul(df_dt,dt_dn)
    
      end subroutine backpropagation
"""
p5 = """      USE readNN"""

p6 = """      REAL(q) Exc_NNM"""
if __name__ == '__main__':
    patch()
