c********************************************************************c
c                                                                    c
c                            VARIA2.F                                c
c                                                                    c
c********************************************************************c
c
c	function varia(x,itype,z1,icall,omega0,lambda0)
c	real lambda0
c
c This routine provides linear spectrums of dark matters in several 
c cosmological models. The spectrums are normalized in one of the two ways. 
c One is by \sigma_8=1 at the present time for pure CDM models (1-4); 
c and the other is by the COBE quadrupole which is 16uK for 
c all MDM models and the HDM model. Scaling with each other is easy.
c The code must be initialized with the same parameters (except
c icall) as the application required. This routine is different from
c the routine varia1.f in ~/topsimu.
c
c Written by Yipeng Jing 26/05/1993; revised on 28/10/1995.
c
c Input variable:
c       x-------the wavenumber in h/Mpc at the present time;
c       itype---the notation for a model as described in following;
c       z1------1+z, z is the redshift of the simulation beginning.
c       icall---0 for initialization; 1 for later use.
c       omega0--the density parameter at z=0;
c       lambda0-the lambda parameter at z=0;
c Output variable:
c       varia()-the power spectrum;
c
c Other routines:
c       growf--the growth factor in the cosmological models of pressure zero.
c
c If itype =1: the power spectrum of cold dark matter with \Gamma=0.5;
c              this is the spectrrum of the standard CDM model:
c              i.e., Omega=1 of cold dark matter, h=0.5, Zeldovich primordial 
c              spectrum, negligible amount of baryon; From Bardeen etal.(1986);
c If itype =2: the power spectrum of cold dark matter with \Gamma=0.25;
c If itype =3: the power spectrum of cold dark matter with \Gamma=0.20;
c If itype =4: the power spectrum of cold dark matter with \Gamma=0.15;
c if itype =5: the CDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.6 of cold dark matter,
c              Omega_b=0.1 of baryons and Omega_h=0.3 of 1 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Holtzman (1989)
c if itype =6: the CDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.69 of cold dark matter,
c              Omega_b=0.01 of baryons and Omega_h=0.3 of 3 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Holtzman (1989)
c if itype =7: the CDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.69 of cold dark matter,
c              Omega_b=0.01 of baryons and Omega_h=0.3 of 1 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Holtzman (1989)
c if itype =8: the CDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.6 of cold dark matter,
c              Omega_b=0.1 of baryons and Omega_h=0.3 of 1 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Klypin et al. (1993).
c if itype =9: the HDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.6 of cold dark matter,
c              Omega_b=0.1 of baryons and Omega_h=0.3 of 1 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Klypin et al. (1993).
c if itype=10: the CDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.6 of cold dark matter,
c              Omega_b=0.1 of baryons and Omega_h=0.3 of 1 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Klypin et al. (1993). (revised version)
c if itype=11: the HDM spectrum in a hybrid model of cold dark matter, 
c              neutrinos and baryon: i.e. Omega_c=0.6 of cold dark matter,
c              Omega_b=0.1 of baryons and Omega_h=0.3 of 1 flavor 
c              of neutrinos; h=0.5 and Zeldovich primordial spectrum.
c              From Klypin et al. (1993). (revised version)
c if itype=12: the HDM spectrum from Bardeen et al. (1986);
c if itype=20: the CDM spectrum of \Gamma=0.2 with a bump;
c
c

      function varia(x,itype,z1,icall,omega0,lambda0)
      real lambda0
      parameter (hubb=0.71,Theta=2.728/2.7,x0=0.0704)
      save const
      if(icall.eq.0)then
         const=4./9.*(2.7e-5)**2*(2.9979e3)**4
     &        *(growf(1./z1,omega0,lambda0)/growf(1.,omega0,lambda0))**2
         if(itype.eq.1)then
            const=const/1.541517973
         else if(itype.eq.2)then
            const=const/0.3229024708
         else if(itype.eq.3.or.itype.eq.20)then
            const=const/0.1854164004
         else if(itype.eq.4)then
            const=const/0.8778662980E-01
         else if(itype.eq.12)then
            const=const/2.081929
         else if(itype.eq.22)then
            omegab0=0.166
            pk=tk(omega0,omegab0,hubb,Theta,0,x)
            const=23367.3125*1.062760
     &        *(growf(1./z1,omega0,lambda0)/growf(1.,omega0,lambda0))**2
         else if(itype.eq.23)then
            omegab0=0.166
            pk=tk(omega0,omegab0,hubb,Theta,0,x)
            const=21768.69*1.063837
     &        *(growf(1./z1,omega0,lambda0)/growf(1.,omega0,lambda0))**2
         else if(itype.eq.24)then
            omegab0=0.166
            pk=tkcf(x,omegab0,0)
            const=21768.69*1.079418
     &        *(growf(1./z1,omega0,lambda0)/growf(1.,omega0,lambda0))**2
         else if(itype.eq.25)then
            omegab0=0.166
            pk=tkcf(x,omegab0,0)
            const=24096.75
     &        *(growf(1./z1,omega0,lambda0)/growf(1.,omega0,lambda0))**2

         else if(itype.eq.30)then
            omegab0=0.155
            pk=tkcf_planck(x,omegab0,0)
            const=17705.33
     &        *(growf(1./z1,omega0,lambda0)/growf(1.,omega0,lambda0))**2
         endif
         return
      endif

      if(x.le.1.e-8.and.itype.ne.22)then
         varia=const*x
         
      else if(itype.eq.1)then   ! CDM spectrum
         omehi=1./0.5
         x1 = x*omehi
         varia=const/x*(log(1.+2.34*x1)/(2.34*omehi*(1.+ 3.89*x1 + 
     &        (16.1*x1)**2 + (5.46*x1)**3 + (6.71*x1)**4 )**0.25))**2

      else if(itype.eq.2)then   ! low-Omega-CDM spectrum
         omehi=1./(0.25)
         x1 = x*omehi
         varia=const/x*(log(1.+2.34*x1)/(2.34*omehi*(1.+ 3.89*x1 + 
     &        (16.1*x1)**2 + (5.46*x1)**3 + (6.71*x1)**4 )**0.25))**2
         
      else if(itype.eq.3)then   ! low-Omega-CDM spectrum
         omehi=1./(0.20)
         x1 = x*omehi
         varia=const/x*(log(1.+2.34*x1)/(2.34*omehi*(1.+ 3.89*x1 + 
     &        (16.1*x1)**2 + (5.46*x1)**3 + (6.71*x1)**4 )**0.25))**2

      else if(itype.eq.4)then   ! low-Omega-CDM spectrum
         omehi=1./(0.15)
         x1 = x*omehi
         varia=const/x*(log(1.+2.34*x1)/(2.34*omehi*(1.+ 3.89*x1 + 
     &        (16.1*x1)**2 + (5.46*x1)**3 + (6.71*x1)**4 )**0.25))**2

      else if(itype.eq.5)then   ! 1 flavor 0.1 baryon
         x1= x*0.5
         varia=const*x/(1.-0.4678*sqrt(x1)+16.95*x1
     &        +12.03*x1**1.5+652.1*x1**2)**2
         
      else if(itype.eq.6)then   ! CDM spectrum in 3 flavors 0.01 baryon
         x1= x*0.5
         varia=const*x/(1.-0.1158*sqrt(x1)+0.2238*x1
     &        +164.6*x1**1.5+458.3*x1**2)**2

      else if (itype.eq.7)then  ! 1 flavor 0.01 baryon
         x1= x*0.5
         varia=const*x/(1.-0.3988*sqrt(x1)+18.43*x1
     &        -3.340*x1**1.5+503.0*x1**2)**2

      else if (itype.eq.8)then  ! the spetrum of CDM at z in mixing model
c of 60% CDM, 10% baryon and 30% neutrinos of one flavor.
         a=1./z1
         x1=x*0.5
         varia=const*log(1.+18.*x1*sqrt(a))/((9.*sqrt(a))*(1.+1.2
     &        *sqrt(x1)-27.*x1+347.*(1.-sqrt(a)/5.)*x1**1.5-18.
     &        *(1.-0.32*a*a)*x1*x1)**2)
         
      else if (itype.eq.9)then  ! the spetrum of HDM at z in mixing model
c of 60% CDM, 10% baryon and 30% neutrinos of one flavor.
         a=1./z1
         x1=x*0.5
         q=x1/sqrt(a)
         varia=const*log(1.+18.*x1*sqrt(a))/((9.*sqrt(a))*(1.+1.2
     &        *sqrt(x1)-27.*x1+347.*(1.-sqrt(a)/5.)*x1**1.5-18.
     &   *(1.-0.32*a*a)*x1*x1)**2)*exp(-q/16.)/(1.+0.03*q+0.67*q**2)
           
      else if (itype.eq.10)then ! the spetrum of CDM at z in mixing model
c     of 60% CDM, 10% baryon and 30% neutrinos of one flavor (revised).
         a=1./z1
         a07=a**0.7
         x1=x*0.5
         varia=const*2.*log(1.+12.5*x1*a07)/((12.5*a07)*(1.+3.93*(1.
     &        -a07/3.2)*sqrt(x1)-60.2*(1.-a07/3.5)*x1+409.*(1.
     &        -a07/5.)*x1**1.5-17.1*(1.-a07*a07/15.)*x1*x1)**2)
         
      else if (itype.eq.11)then ! the spetrum of HDM at z in mixing model
c     of 60% CDM, 10% baryon and 30% neutrinos of one flavor. (revised)
         a=1./z1
         x1=x*0.5
         q=x1/sqrt(a)
         varia=const*2.*log(1.+12.5*x1*a07)/((12.5*a07)*(1.+3.93*(1.
     &        -a07/3.2)*sqrt(x1)-60.2*(1.-a07/3.5)*x1+409.*(1.
     &        -a07/5.)*x1**1.5-17.1*(1.-a07*a07/15.)*x1*x1)**2)
     &        *exp(-q/23.)/(1.+0.2*q+0.75*q**2)

      else if(itype.eq.12)then  ! HDM spectrum
         omehi=1./0.765
         x1 = x*omehi
         varia=const*x*exp(-0.32*(2.6*x1)-(2.6*x1)**2)/
     &        (1+1.6*x1+(4.*x1)**1.5+(0.92*x1)**2)**2
         
      else if(itype.eq.20)then  ! low-Omega-CDM spectrum
         omehi=1./(0.20)
         x1 = x*omehi
         varia=const/x*(log(1.+2.34*x1)/(2.34*omehi*(1.+ 3.89*x1 + 
     &        (16.1*x1)**2 + (5.46*x1)**3 + (6.71*x1)**4 )**0.25))**2
     &        *(1.+exp(-((x-0.06)/0.015)**2*0.5))
      else if(itype.eq.22)then  !WMAP best-fitting model and Eisenstein & Hu P(k)
         if(x.le.7.65e-4)then
            varia=const*(x/x0)
         else
            omegab0=0.166
            x1=x*hubb
            tkx=tk(omega0,omegab0,hubb,Theta,1,x1)
            varia=const*(x/x0)**(0.93-0.0155*log(x/x0))*tkx*tkx 
         endif

      else if(itype.eq.23)then  !WMAP best-fitting model and Eisenstein & Hu P(k)
         omegab0=0.166
         x1=x*hubb
         tkx=tk(omega0,omegab0,hubb,Theta,1,x1)
         varia=const*(x/x0) *tkx*tkx 
      else if(itype.eq.24)then  !WMAP CMBFAST
         omegab0=0.166
         tkx=tkcf(x,omegab0,1)
         varia=const*(x/x0) *tkx*tkx 
      else if(itype.eq.25)then  !WMAP CMBFAST wmap ns=0.968
         omegab0=0.166
         tkx=tkcf(x,omegab0,1)
         varia=const*(x/x0)**0.968 *tkx*tkx 
      else if(itype.eq.30)then  ! CMBFAST Planck ns=0.9603
         omegab0=0.155
         tkx=tkcf_planck(x,omegab0,1)
         varia=const*(x/x0)**0.9603 *tkx*tkx 
      endif
        
      end

      function tk(Omega0,Omegab0,h,Theta,set,kk)
      parameter (e=2.7183)
      integer set
      real kk,f,Tc
      real sv,j0,Tb,Th,Ti
      real Tj
      real s,zeq,zd,Req,Rd,b1,b2,keq,alphac,a1,a2,betac
      real alphab,betanode,betab,ksilk,bb1,bb2,G
      save Omega0h2,zeq,keq,b2,b1,zd,Rd,Req,s,ksilk,a1,a2,
     &    alphac,bb1,bb2,betac,alphab,betab,betanode


      if (set.eq.0) then
c     define some constants	
         Omega0h2=Omega0*h**2
         zeq=2.5*10**4*Omega0h2*Theta**(-4)
         keq=7.46*0.01*Omega0h2*Theta**(-2)
         b2=0.238*Omega0h2**(0.223)
         b1=0.313*Omega0h2**(-0.419)*(1+0.607*Omega0h2**(0.674))
         zd=1291.*Omega0h2**(0.251)/(1+0.659*Omega0h2**0.828)
     $        *(1+b1*(Omega0h2*Omegab0)**b2)
         Rd=31.5*Omega0h2*Omegab0*Theta**(-4)/(zd/(10**3))
         Req=31.5*Omega0h2*Omegab0*Theta**(-4)/(zeq/(10**3))
         s=2./(3.*keq)*sqrt(6./Req)*(log(sqrt(1.+Rd)+sqrt(Rd+Req))
     $        -log(1+sqrt(Req)))
         ksilk=1.6*(Omega0h2*Omegab0)**(0.52)*Omega0h2**(0.73)
     $        *(1+(10.4*Omega0h2)**(-0.95))
         a1=(46.9*Omega0h2)**(0.67)*(1+(32.1*Omega0h2)**(-0.532))
         a2=(12.*Omega0h2)**(0.424)*(1+(45.*Omega0h2)**(-0.582))
         alphac=a1**(-Omegab0)*a2**(-Omegab0**3)
         bb1=0.944/(1+(458.*Omega0h2)**(-0.708))
         bb2=(0.395*Omega0h2)**(-0.0266)
         betac=1./(1+bb1*((1-Omegab0)**bb2-1))
         alphab=2.07*keq*s*(1+Rd)**(-0.75)*G((1+zeq)/(1+zd))
         betab=0.5+Omegab0+(3.-2.*Omegab0)*sqrt((17.2*Omega0h2)**2+1)
         betanode=8.41*Omega0h2**(0.435)
	endif
        call T0(kk,alphac,betac,Th,num,keq,e)
        call T0(kk,1.,betac,Ti,num,keq,e)
        call T0(kk,1.,1.,Tj,num,keq,e)
        f=1./(1.+(kk*s/5.4)**4)
        Tb=Tj/(1.+(kk*s/5.2)**2)+alphab/
     $       (1.+(betab/kk/s)**3)
     &       *e**(-(kk/ksilk)**1.4)
        sv=s/((1.+(betanode/kk/s)**3)**(1./3.))
        j0=(sin(kk*sv))/(kk*sv)
        Tc=f*Ti+(1.-f)*Th
        tk=Omegab0*Tb*j0+(1.-Omegab0)*Tc
	end
	
	subroutine T0(kk,x2,x3,Tv,num,keq,e)
	real kk,Tv,keq,e,q
	real x2,x3,c
	q=kk/(13.41*keq)
        c=14.2/x2+386./(1.+69.9*q**1.08)
        Tv=log(e+1.8*x3*q)/(log(e+1.8*x3*q)+c*q**2)
	end

	function G(y)
	G=y*(-6.*sqrt(1.+y)+(2.+3.*y)*(log(sqrt(1.+y)+1.)-
     &    log(sqrt(1.+y)-1.)))
	 end   
	    
	
c-----read the transfer function given by the cmbfast, always be careful
c-----that the data in tarans0.dat is what you want
c
        function tkcf(k,Omegab0,set)
        parameter(nbin=1000)
        real k,Omegab0,tk(nbin),kk(nbin),tkc,tkb
        integer set
        save tk,kk,kbin,dlogk

        if(set.eq.0)then
           kbin=0
           print*,"WARNING: transf0.dat consistent with parameters?"
           open(1,file="/home/ypjing/adaptive/src/transf0.dat",
     &          status="old")
 2         read(1,*,end=199)k,tkc,tkb
           kbin=kbin+1
           kk(kbin)=k
           tk(kbin)=tkc*(1.-Omegab0)+tkb*Omegab0
           goto 2
 199       close(1)
           do i=2,kbin
              tk(i)=tk(i)/tk(1)
           enddo
           tk(1)=1.
           dlogk=log(kk(2))-log(kk(1))
           return
        endif

        if(k.le.kk(1))then
           tkcf=1.
        else 
           ibin=(log(k)-log(kk(1)))/dlogk+1.
           if(ibin.le.kbin-1)then
              tkcf=(tk(ibin+1)-tk(ibin))/(kk(ibin+1)-kk(ibin))*
     &             (k-kk(ibin))+tk(ibin)
           else
              print*,"k too large in tkcf"
              stop
           endif
        endif
        end
	

	
c-----read the transfer function given by the cmbfast, always be careful
c-----that the data in planck_transf.dat is what you want
c
        function tkcf_planck(k,Omegab0,set)
        parameter(nbin=10000)
        real k,Omegab0,tk(nbin),kk(nbin)
        integer set
        save tk,kk,kbin,dlogk
        if(set.eq.0)then
           kbin=0
           print*,"WARNING: planck_transf.dat consistent with 
     & parameters?"
           open(1,file="/home/ypjing/adaptive/src/planck_transf.dat",
     &          status="old")
c           open(1,file="planck_transf.dat",
c     &          status="old")
 2         read(1,*,end=199)k,tkc,tkb
           kbin=kbin+1
           kk(kbin)=k
           tk(kbin)=tkc*(1.-Omegab0)+tkb*Omegab0
           goto 2
 199       close(1)
           do i=2,kbin
              tk(i)=tk(i)/tk(1)
           enddo
           tk(1)=1.
           dlogk=log(kk(2))-log(kk(1))
           return
        endif

        if(k.le.kk(1))then
           tkcf_planck=1.
        else 
           ibin=(log(k)-log(kk(1)))/dlogk+1.
           if(ibin.le.kbin-1)then
              tkcf_planck=(tk(ibin+1)-tk(ibin))/(kk(ibin+1)-kk(ibin))*
     &             (k-kk(ibin))+tk(ibin)
           else
              print*,"k too large in tkcf_planck"
              stop
           endif
        endif

        end
	
	
	
	







