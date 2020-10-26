c    Declarations for Main file

      integer a, b
      parameter (a=9, b=10)

      integer nFlt, nDip, nSR(MAX_FLT), nThick(MAX_FLT), imag, isr, 
     1        nfp(MAX_FLT), nwidth(MAX_FLT), nDD(MAX_FLT), ifault, ib,
     2        ifltwidth, nPts, nbvalue(MAX_FLT), nMagChar(MAX_FLT,MAX_N1),
     3        nMag, idM, runFlag  
      
      real DeltaDip(MAX_FLT,MAX_N1), DipWt(MAX_FLT,MAX_N1), SR(MAX_FLT,MAXPARAM),
     1     SRWt(MAX_FLT,MAXPARAM), Thick(MAX_FLT,MAX_N1), ThickWt(MAX_FLT,MAX_N1), 
     2     MagChar(MAX_FLT,MAX_WIDTH,MAX_N1), MagCharWt(MAX_FLT,MAX_WIDTH,MAX_N1),
     3     dip(MAX_FLT,MAX_WIDTH,MAX_SEG), fault_width(MAX_FLT,MAX_WIDTH),
     4     fZ(MAX_FLT,MAX_DD,MAX_SEG), flat(MAX_FLT,MAX_DD,MAX_SEG), 
     5     flong(MAX_FLT,MAX_DD,MAX_SEG), siteX, siteY, xFlt(MAX_DD,MAX_SEG),
     6     yFlt(MAX_DD,MAX_SEG), zFlt(MAX_DD,MAX_SEG), x0, y0, z0, fault_area, 
     7     xstep(MAX_FLT), bvalue(MAX_FLT,MAX_N1), bvalueWt(MAX_FLT,MAX_N1), avalue,
     8     mpdf_param(MAX_FLT,MAX_WIDTH,6), Mminhaz(MAX_FLT), coef_area(2,MAX_FLT),
     9     magStep(MAX_FLT), occurRate(MAX_DM)   
