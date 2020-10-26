      program HAZ_to_OQ

c     This program reads in HAZ45 input files and calculates requisite inputs for OpenQuake input files
c     Limited to specific applications
c     Current capabilities: 
c     calculate a-values for Truncated Gutenberg Richter Mag Distribution
c     calculate occurrence rates for WAACY Mag Distribution

      implicit none
      include 'max_dims.h' 
      include 'declare.h'
   
c     Write Program information to the screen.
      write (*,*) '*********************************'
      write (*,*) '*        HAZ_to_OQ Code         *'
      write (*,*) '*         October, 2020          *'
      write (*,*) '*********************************'
      write (*,*)

      write (*,*) 'Enter 1 to calculate a-values for Truncated Gutenberg Richter Mag PDF'
      write (*,*) 'Enter 2 to calculate occurrence rates for WAACY Mag PDF'
      read (*,*) runFlag

c     Read in run file
      call read_run (a, b, siteX, siteY)

c     Read in info from fault file
      call Read_fault_data ( b, nflt, nfp, nDip, DeltaDip, Dip, DipWt, nSR, SR, SRWt, 
     1     nbvalue, bvalue, bvalueWt, nThick, Thick, ThickWt, nMagChar, MagChar, 
     2     MagCharWt, nwidth, fault_width, fZ, flat, flong, nDD, xstep, mpdf_param,
     3     Mminhaz, coef_area, magStep )              

c     Write headers for output files
      write (40,*) "Thick ", "Dip ", "fault_area ", "Mchar "
      if (runFlag .eq. 1) then
        write (41,*) "Slip Rate ", "avalue "
      else if (runFlag .eq. 2) then
        write (41,*) "Slip Rate ", "Occurrence Rates "
      endif
      if (runFlag .eq. 1) then
        write (42,*) "sum_inc_prob "
      else if (runFlag .eq. 2) then
        write (42,*) "sumpmag ", "sumpmag1 ", "sumpmag2 ", "sumpmag3 "
      endif

c     Loop over sources
      do ifault=1,nflt
      
c     Loop over fault widths (combo crustal thickness and dip)
      do ifltwidth=1,nwidth(ifault)
      
c         Set bottom of fault
          call Flt_bottom(nfp(ifault), ifault, Dip(ifault,ifltwidth,1),  
     1         fault_width(ifault,ifltwidth), fZ, flat, flong, nDD(ifault))

c         Convert lon,lat to x,y in km 
          call Convert_coordinates (nfp(ifault), ifault, nDD(ifault), 
     1         siteX, siteY, fLat, fLong, fZ, nPts, xFlt, yFlt, zFlt, x0, y0, z0) 

c         Turn fault into a grid and calculate fault area
          call Calc_fault_area ( xFlt, yFlt, zFlt, nfp(ifault), nDD(ifault), x0, y0, z0,
     1         xstep(ifault), fault_area)  

c         Loop over characteristic magnitude
          do imag=1,nMagChar(ifault,ifltwidth)
          
c           Loop over slip rate
            do isr=1,nSR(ifault)
            
c             Loop over bvalue
              do ib=1,nbvalue(ifault)

                if (runFlag .eq. 1) then            
c               Calculate a-values for truncated exponential magnitude distribution
                call calc_avalue (fault_area, MagChar(ifault,ifltwidth,imag), SR(ifault,isr), 
     1               bvalue(ifault,ib), avalue)             

                else if (runFlag .eq. 2) then
c               Calculate occurrence rates for WAACY magnitude distribution
                call calc_occurRates_WAACY (fault_area, MagChar(ifault,ifltwidth,imag), SR(ifault,isr),
     1               bvalue(ifault,ib), mpdf_param, ifault, ifltwidth, Mminhaz(ifault), coef_area, 
     2               magStep(ifault), occurRate, nMag)    
                endif
        

                write (40,*) fault_width(ifault,ifltwidth), Dip(ifault,ifltwidth,1), fault_area, MagChar(ifault,ifltwidth,imag)
                if (runFlag .eq. 1) then
                  write (41,*) SR(ifault,isr), avalue
                else if (runFlag .eq. 2) then
                  write (41,*) SR(ifault,isr), (occurRate(idM),idM=1,nMag)
                endif
                
              enddo
            enddo          
          enddo
        enddo        
      enddo       
      close(40)
      close(41)
      close(42)

      write (*,*) 'Finished'
      pause
      
      end
            
c -----------------------------------------------------------------

      subroutine read_run (a, b, siteX, siteY)
      
      implicit none

      integer a, b
      real siteX, siteY
      character*100 runfile, faultfile    
      
        write (*,*) 'Enter the input file name to run the model.'
        read (*,'( a100)') runfile
        open (a,file=runfile,status='old')             
        read (a,'( a80)') faultfile
        open (b,file=faultfile,status='old')
        read (a,*) siteX
        read (a,*) siteY

      end        
