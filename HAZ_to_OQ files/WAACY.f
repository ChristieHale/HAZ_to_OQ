
      subroutine calc_occurRates_WAACY ( fault_area, Mchar, SR, bvalue, mpdf_param,
     1                                   ifault, ifltwidth, Mminhaz, coef_area, magStep,
     2                                   occurRate, nMag1 )
      
c     this subroutine calculates the occurance rates for the WAACY magnitude pdf
c     reference: "Alternative Magnitude-Frequency Distribution for Evaluating Seismic
c                 Hazard for Multi-Segment Ruptures" Wooddell, Abrahamson, 
c                 Acevedo-Cabrera, Youngs. Draft provided by Wooddell on 10/1/2020

      implicit none
      include 'max_dims.h'

      integer ifault, ifltwidth, nMag, iMag, nMag1, jMag
      real fault_area, Mchar, SR, bvalue, mpdf_param(MAX_FLT,MAX_WIDTH,6), 
     1     pi, sigmaM, Mmaxspec, bvalue_tail, F, Mcharadd, Mcap, Mmaxcheck,
     2     Mmax, beta, Nsig, beta_tail, M1, M2, Mmin, c1, c2, c3, M0exp1, 
     3     M0char, M0exp2, d1, d2, a, a1, dM, Mag, fWaacy, pmag(MAX_DM),
     4     sumtop, sumbot, M0(MAX_DM), rup_area, MA_b, MA_a, Aratio(MAX_DM),
     5     mult, F1, F_prime, pmag1(MAX_DM), sumbot1, Mminhaz, mu, N_Mmin0,
     6     N_Mminhaz, coef_area(2,MAX_FLT), phi10, phi20, magStep, pmag2(MAX_DM),
     7     occurRate(MAX_DM), Mag1, Magl, Magh, pmag3(MAX_DM), sumpmag, sumpmag1,
     8     sumpmag2, sumpmag3 
      real*8 phi1, phi2 

c     calculate constant
      pi = 4. * atan(1.)

c     set WAACY parameters from input file
      sigmaM = mpdf_param(ifault,ifltwidth,1)
      Mmaxspec = mpdf_param(ifault,ifltwidth,2)
      bvalue_tail = mpdf_param(ifault,ifltwidth,3)
      F = mpdf_param(ifault,ifltwidth,4)
      Mcharadd = mpdf_param(ifault,ifltwidth,5)
      Mcap = mpdf_param(ifault,ifltwidth,6)
      
c     set magnitude area parameters from input file      
      MA_a = coef_area(1,ifault)
      MA_b = coef_area(2,ifault)
      
c     set maximum magnitude (confusing - go over this with Norm)
      if (Mcharadd .eq. 0) then
        Mmaxcheck = Mmaxspec
      else
        Mmaxcheck = Mchar + Mcharadd
      endif
      if (Mmaxcheck .gt. Mcap) then
        Mmax = Mcap
      else
        Mmax = Mmaxcheck
      endif 

c     calculate additional WAACY parameters
      beta = alog(10.0)*bvalue
      Nsig = 1.5 !fixed per documentation
      beta_tail = alog(10.0)*bvalue_tail
      M1 = Mchar - 0.25
      M2 = Mchar + Nsig*sigmaM

c     set minimum magnitude
      Mmin = 0.
      
c     calculate WAACY equations (not dependent on F or Mmin)
      c1 = (1/(sqrt(2*pi)*sigmaM))*exp(((-1)*(1.5*sigmaM)**2)/(2*(sigmaM**2)))  !eq. 2
      call ndtr((-0.25/sigmaM), phi1) 
      call ndtr(1.5, phi2)
      phi10 = real(phi1)
      phi20 = real(phi2)
      c2 = 1/(phi1-(phi2))  !eq. 3
      c3 = (1-exp(-beta_tail*(Mmax-M2)))/beta_tail  !eq. 4
      M0exp1 = (beta*(10**16.05)*(exp((-beta+3.45)*M1)-1))/((1-exp(-beta*m1))*(-beta+3.45))  !eq. 9 
      M0char = (10**(1.5*Mchar+16.05))*(2.63*(sigmaM-0.2)+1.19)  !eq. 10
      M0exp2 = (beta_tail*(10**(16.05))*exp(beta_tail*M2)*(exp((-beta_tail+3.45)*Mmax)-exp((-beta_tail+3.45)*M2)))/
     1         ((1-exp(-beta_tail*(Mmax-M2)))*(-beta_tail+3.45))  !eq. 11
      d2 = ((1/(1+c1*c2*c3))*M0char)+(((c1*c2*c3)/(1+c1*c2*c3))*M0exp2)  !eq. 8 
      
c     calculate WAACY equations (dependent on F or Mmin)  
      d1 = M0exp1/(exp(-beta*Mmin)-exp(-beta*M1))  !eq. 7    
      a = (d2*F)/(d1*(1-F)+d2*F)  ! eq. 5
      a1 = (d1*(1-F))/(d1*(1-F)+(d2*F))  ! eq. 6
      
c     set magnitude discretization (small for moment balancing)
      dM = 0.01
      
c     calculate magnitude pdf (using F)
c     loop over magnitudes
      nMag = (Mmax - Mmin)/dM      
      do iMag=1,nMag  
        Mag = Mmin + (iMag-0.5)*dM  
        if (Mag .gt. Mmin .and. Mag .lt. M1) then
          fWaacy = a*(beta*exp(-beta*(Mag-Mmin)))/(1-exp(-beta*(M1-Mmin)))  !eq. 1        
        else if (Mag .ge. M1 .and. Mag .le. M2) then
          fWaacy = a1*(c2/(1+c1*c2*c3))*(1/(sqrt(2*pi)*sigmaM))*exp(((-1)*((Mag-Mchar)**2))/(2*sigmaM**2))  !eq. 1
        else if (Mag .gt. M2 .and. Mag .lt. Mmax) then
          fWaacy = a1*((c1*c2)/(1+c1*c2*c3))*exp(-beta_tail*(Mag-M2))  !eq. 1
        else
          write (*,*) 'magnitude increment is outside of Mmin, Mmax range'
          pause
        endif 
        pmag(iMag) = fWaacy*dM
      enddo
c     check sum of pdf
      sumpmag = 0.
      do iMag=1,nMag
        sumpmag = sumpmag + pmag(iMag)
      enddo 

c     calculate F'
      sumtop = 0.
      sumbot = 0.
      do iMag=1,nMag
        Mag = Mmin + (iMag-0.5)*dM 
        M0(iMag) = 10**((1.5*Mag)+16.05)
        rup_area = 10**(MA_b*Mag+MA_a)
        Aratio(iMag) = min(fault_area/rup_area,1.)
        mult = pmag(iMag)*M0(iMag)*Aratio(iMag)
        if (Mag .le. M1) then
          sumtop = sumtop + mult
        endif
        if (Mag .le. Mmax) then
          sumbot = sumbot + mult
        endif
      enddo
      F1 = sumtop / sumbot !eq. 20
      F_prime = (F/F1)*F  !eq. 19

c     Redo - calculate WAACY equations (dependent on F')      
      a = (d2*F_prime)/(d1*(1-F_prime)+d2*F_prime)  ! eq. 5
      a1 = (d1*(1-F_prime))/(d1*(1-F_prime)+(d2*F_prime))  ! eq. 6
        
c     Redo - calculate magnitude pdf (using F')
c     loop over magnitudes
      do iMag=1,nMag  
        Mag = Mmin + (iMag-0.5)*dM  
        if (Mag .gt. Mmin .and. Mag .lt. M1) then
          fWaacy = a*(beta*exp(-beta*(Mag-Mmin)))/(1-exp(-beta*(M1-Mmin)))  !eq. 1        
        else if (Mag .ge. M1 .and. Mag .le. M2) then
          fWaacy = a1*(c2/(1+c1*c2*c3))*(1/(sqrt(2*pi)*sigmaM))*exp(((-1)*((Mag-Mchar)**2))/(2*sigmaM**2))  !eq. 1
        else if (Mag .gt. M2 .and. Mag .lt. Mmax) then
          fWaacy = a1*((c1*c2)/(1+c1*c2*c3))*exp(-beta_tail*(Mag-M2))  !eq. 1
        else
          write (*,*) 'magnitude increment is outside of Mmin, Mmax range'
          pause
        endif 
        pmag1(iMag) = fWaacy*dM
      enddo    
c     check sum of pdf
      sumpmag1 = 0.
      do iMag=1,nMag
        sumpmag1 = sumpmag1 + pmag1(iMag)
      enddo 

c     calculate N(M>Mminhaz)
      sumbot = 0.
      sumbot1 = 0.
      do iMag=1,nMag
        Mag = Mmin + (iMag-0.5)*dM 
        mult = pmag1(iMag)*M0(iMag)*Aratio(iMag)
        sumbot = sumbot + mult
        if (Mag .gt. Mminhaz) then
          sumbot1 = sumbot1 + pmag1(iMag)
        endif
      enddo      
c     set mu = 3E+11 dyne/cm2      
      mu = 3*10.**11.      
c     calculate rate (convert SR from mm to cm, convert fault area from km2 to cm2)               
      N_Mmin0 = fault_area*SR*mu*0.1*100000*100000/sumbot  !eq. 15
      N_Mminhaz = N_Mmin0*sumbot1  !eq. 16

c     Redo - calculate WAACY equations (dependent on F' and Mminhaz)
      d1 = M0exp1/(exp(-beta*Mminhaz)-exp(-beta*M1))  !eq. 7      
      a = (d2*F_prime)/(d1*(1-F_prime)+d2*F_prime)  ! eq. 5
      a1 = (d1*(1-F_prime))/(d1*(1-F_prime)+(d2*F_prime))  ! eq. 6

c     calculate occurence rates for OpenQuake
c     Redo - calculate magnitude pdf (using F' and Mminhaz)
c     set magnitude discretization (from input file) but sum using smaller dM
      nMag = (Mmax - Mminhaz)/dM
      do iMag=1,nMag  
        Mag = Mminhaz + (iMag-0.5)*dM   
        if (Mag .gt. Mminhaz .and. Mag .lt. M1) then
          fWaacy = a*(beta*exp(-beta*(Mag-Mminhaz)))/(1-exp(-beta*(M1-Mminhaz)))  !eq. 1        
        else if (Mag .ge. M1 .and. Mag .le. M2) then
          fWaacy = a1*(c2/(1+c1*c2*c3))*(1/(sqrt(2*pi)*sigmaM))*exp(((-1)*((Mag-Mchar)**2))/(2*sigmaM**2))  !eq. 1
        else if (Mag .gt. M2 .and. Mag .lt. Mmax) then
          fWaacy = a1*((c1*c2)/(1+c1*c2*c3))*exp(-beta_tail*(Mag-M2))  !eq. 1
        else
          write (*,*) 'magnitude increment is outside of Mmin, Mmax range'
          pause
        endif 
        pmag2(iMag) = fWaacy*dM
      enddo 
c     check sum of pdf
      sumpmag2 = 0.
      do iMag=1,nMag
        sumpmag2 = sumpmag2 + pmag2(iMag)
      enddo 
 
      nMag1 = (Mmax - Mminhaz)/magstep
      do jMag=1,nMag1
        Mag1 = Mminhaz + (jMag-0.5)*magstep
        Magl = Mag1 - (0.5*magstep)
        Magh = Mag1 + (0.5*magstep)
        pmag3(jMag) = 0.
        do iMag=1,nMag
          Mag = Mminhaz + (iMag-0.5)*dM
          if (Mag .lt. Magh .and. Mag .gt. Magl) then
            pmag3(jMag) = pmag3(jMag) + pmag2(iMag)
          endif
        enddo
        occurRate(jMag) = N_Mminhaz*pmag3(jMag)
      enddo
c     check sum of pdf
      sumpmag3 = 0.
      do jMag=1,nMag1
        sumpmag3 = sumpmag3 + pmag3(jMag)
      enddo 

c     check that probability density functions sum to 1
      write (42,*) sumpmag, sumpmag1, sumpmag2, sumpmag3

      return 
      end

c ------------------------------------------------------------------

      subroutine ndtr( X, P)    
C          Reference: Abramowitz and Stegan equation 7.1.26                                             
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                               

      implicit none

      real x
      real*8 p, x1, x2, p1, a1, a2, a3, a4, a5, t
      data p1, a1, a2, a3, a4, a5 / 0.3275911, 0.254829592, 
     1     -0.284496736, 1.421413741, -1.453152027, 1.061405429 /
     
      if ( x .lt. 0. ) then
	    x1 = abs(x)
	  else
	    x1 = x
      endif
	  
        x2 = x1/(sqrt(2.))
        t = 1/(1+(p1*x2))
        p = 1-0.5*((a1*t)+(a2*(t**2))+(a3*(t**3))+(a4*(t**4))+(a5*(t**
     1      5)))*(exp(-(x2**2)))
     
      if ( x .gt. 0. ) then
	    p = 1. - p
      endif

      return 
      end
      
c ------------------------------------------------------------------
      
