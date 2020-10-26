
      subroutine calc_avalue ( fault_area, Mmax, SR, bvalue, avalue )

c     this subroutine calculates the a-value for a truncated gutenberg richter
c     magnitude pdf       

      implicit none
      include 'max_dims.h'

      integer nMag, iMag
      real fault_area, Mmax, bvalue, beta, dM, Mmin, Mag, 
     1     sum_inc_prob, sum_inc_Moprob, a, b, inc_prob, inc_Mo, inc_Mprob,
     2     inc_Moprob, Mean_Mo_eqk, mu, NMmin, SR, avalue

c     calculate beta
      beta = alog(10.0)*bvalue

c     set magnitude discretization (small for moment balancing)
      dM = 0.01

c     set Mmin (for moment balancing)
      Mmin = 0.0
      
c     loop over magnitudes
      nMag = (Mmax - Mmin)/dM
      sum_inc_prob = 0.0
      sum_inc_Moprob = 0.0
      do iMag=1,nMag  
        Mag = Mmin + (iMag-0.5)*dM
        a = beta*exp(-beta*(Mag-Mmin))
        b = 1 - exp(-beta*(Mmax-Mmin))
        inc_prob = a/b*dM
        inc_Mo = 10**(1.5*Mag+16.05)
        inc_Mprob = Mag*inc_prob
        inc_Moprob = inc_Mo*inc_prob
        
        sum_inc_prob = sum_inc_prob + inc_prob
        sum_inc_Moprob = sum_inc_Moprob + inc_Moprob 
      enddo

c     check that probability density function sums to 1
      write (42,*) sum_inc_prob

      Mean_Mo_eqk = sum_inc_Moprob
      
c     set mu = 3E+11 dyne/cm2      
      mu = 3*10.**11.
      
c     calculate rate (convert SR from mm to cm, convert fault area from km2 to cm2)      
      NMmin = fault_area*SR*mu*0.1*100000*100000/Mean_Mo_eqk
      
      avalue = log10(NMmin)

      return 
      end
