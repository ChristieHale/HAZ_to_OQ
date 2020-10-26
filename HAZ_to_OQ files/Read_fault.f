c     compatible with HAZ45.3   

      subroutine Read_fault_data ( b, nFlt, nfp, n_Dip, deltaDip1, dip, dipWt1,
     1     nSR1, sr1, wt_sr1, n_bValue3, bValue3, bValueWt3, nThick1, faultThick1, 
     2     faultThickWt1, nRefMag2, refMag2, refMagWt2, nWidth, faultWidth, fZ, 
     3     fLat, fLong, nDownDip, hxStep, mpdf_param1, minMag, coef_area, magStep) 
     
      implicit none
      include 'max_dims.h'

      integer b, synHWFlag(MAX_FLT,MAX_SYN), nsyn_Case(MAX_FLT), 
     1        synjcalc(MAX_FLT), fltDirect(MAX_FLT), synchron(MAX_FLT),      
     2        fIndex(3,MAX_FLT), nWidth(MAX_FLT), sourceType(MAX_FLT), 
     3        attenType(MAX_FLT), grid_n(MAX_FLT), n_bValue3(MAX_FLT),
     4        nfp(MAX_FLT), nDownDip(MAX_FLT), HWsource7, iRate, nb1 
      integer nMag(MAX_FLT), nRupArea(MAX_FLT), nRupWidth(MAX_FLT),
     1        nRate, n_bValue, nRefMag(MAX_N1), ifsystem, isyn, 
     2        nParamVar(MAX_FLT,MAX_WIDTH), nfsystem,  
     3        fsys(MAX_FLT), iDepthModel(MAX_FLT),
     4        nFtype(MAX_FLT), faultFlag(MAX_FLT,100,MAX_FLT) 
      integer nFtype1(MAX_FLT), nRefMag2(MAX_FLT,MAX_N1), iDip, 
     1        iWidth, nThick1, nSR, nSR1(MAX_FLT), nMoRate, nRecInt, 
     2        ii, ipt, nFlt, iCoor, iFlt0, k, nFlt2, i, iflt2, igrid, 
     3        n_Dip, nActRate, iRecur, iThickDip, iThick1, nRefMag0, iFM 
      integer iflt, nSegModel, nMagRecur, nFtypeModels, 
     1        nFM, iRefMag, i_bValue, segModelFlag(MAX_FLT,100), 
     2        nSegModel0(MAX_FLT), ncountS7(MAX_FLT)
      real synmag(MAX_FLT,MAX_SYN), syndistRup(MAX_FLT,MAX_SYN),
     1     syndistJB(MAX_FLT,MAX_SYN), syndistSeismo(MAX_FLT,MAX_SYN), 
     2     synftype(MAX_FLT,MAX_SYN), synhypo(MAX_FLT,MAX_SYN), 
     3     synwt(MAX_FLT,MAX_SYN), syn_dip(MAX_FLT,MAX_SYN), 
     4     syn_zTOR(MAX_FLT,MAX_SYN), syn_RupWidth(MAX_FLT,MAX_SYN)
      real syn_RX(MAX_FLT,MAX_SYN), syn_Ry0(MAX_FLT,MAX_SYN),
     1     rateParam1(MAX_N1), RateWt1(MAX_N1), al_segWt(MAX_FLT), 
     2     dipWt1(MAX_N1), deltaDip1(MAX_N1), bValue1(MAX_N1), 
     3     bValueWt1(MAX_N1), bValue2(MAX_N1), faultThick1(MAX_N1), 
     4     faultThickWt1(MAX_N1), refMag1(MAX_N1,MAX_N1)
      real refMagWt1(MAX_N1,MAX_N1), refMag0(MAX_N1), refMagWt0(MAX_N1),
     1     magRecurWt1(MAX_N1), magRecur1(MAX_N1), probAct(MAX_FLT),
     2     minlat, maxlat, minlong, maxlong, grid_a(MAX_FLT,MAX_GRID),
     3     grid_lat(MAX_FLT,MAX_GRID), grid_long(MAX_FLT,MAX_GRID), 
     4     grid_top(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT) 
      real grid_dlat(MAX_FLT), minMag(MAX_FLT), magStep(MAX_FLT), 
     1     hxStep(MAX_FLT), hyStep(MAX_FLT), minDepth(MAX_FLT), 
     2     sampleStep(MAX_FLT), rP1(MAXPARAM), 
     3     rateParam(MAX_FLT,MAXPARAM,MAX_WIDTH), rP2(MAXPARAM),
     4     rateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH), rp3(MAXPARAM)
      real beta(MAX_FLT,MAXPARAM,MAX_WIDTH), rp4(MAXPARAM), mtest,
     1     magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH), rp5(MAXPARAM),
     2     magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH), sigWidth(MAX_FLT),
     3     faultWidth(MAX_FLT,MAX_WIDTH), ftype_wt(MAX_FLT,MAX_N1),
     4     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6), rp6(MAXPARAM)
      real maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_width(2,MAX_FLT),
     1     maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_area(2,MAX_FLT), 
     2     sigArea(MAX_FLT), fLong(MAX_FLT,MAX_DD,MAX_SEG), 
     3     fLat(MAX_FLT,MAX_DD,MAX_SEG), fZ(MAX_FLT,MAX_DD,MAX_SEG), 
     4     ftype(MAX_FLT,MAX_N1), faultWidthWt(MAX_FLT,MAX_WIDTH)
      real ftype1(MAX_FLT,MAX_N1), ftype_wt1(MAX_FLT,MAX_N1), 
     1     ftmodelwt(MAX_N1), wt_srBranch, wt_ActrateBranch,
     2     mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM), dip1, top, 
     3     scaleRate(MAX_FLT), segWt(MAX_FLT,MAX_FLT)    
      real mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM), wt_recIntBranch,
     1     sr(MAXPARAM), wt_sr(MAXPARAM), actRate(MAXPARAM), 
     2     actRateWt1(MAXPARAM), MoRate(MAXPARAM), wt_MoRate(MAXPARAM),
     3     MoRateDepth(MAXPARAM), MoRDepth(MAXPARAM), segWt1(MAX_FLT),
     4     rec_Int(MAXPARAM), wt_recInt(MAXPARAM), dip2, testMaxMag    
      real dip(MAX_FLT,MAX_WIDTH, MAX_SEG), segModelWt1(MAX_FLT,100), 
     1     wt_MoRateBranch, sum, ProbAct0, depthParam(MAX_FLT,5), 
     2     rateType1(MAXPARAM), RateType(MAX_FLT,MAXPARAM,MAX_WIDTH), 
     3     magS7(MAX_FLT,MAX_S7), rateS7(MAX_FLT,MAX_S7), 
     4     distS7(MAX_FLT,MAX_S7), DipS7(MAX_FLT,MAX_S7) 
      real mechS7(MAX_FLT,MAX_S7), refMag2(MAX_FLT,MAX_WIDTH,MAX_N1), 
     1     refMagWt2(MAX_FLT,MAX_WIDTH,MAX_N1), sr1(MAX_FLT,MAXPARAM),
     2     wt_sr1(MAX_FLT,MAXPARAM), bValue3(MAX_FLT,MAX_N1), bValueWt3(MAX_FLT,MAX_N1),
     3     mpdf_param1(MAX_FLT,MAX_WIDTH,6)
      character*80 fName(MAX_FLT), fName1

c     Input Fault Parameters
      read (b,*,err=3001) iCoor
      read (b,*,err=3002) NFLT

      iflt = 0
      ifsystem = 1
      
      DO iFlt0=1,NFLT
       read (b,'( a80)',err=3003) fName1
       read (b,*,err=3004) probAct0

c      Read number of segmentation models for this fault system
       read (b,*,err=3005) nSegModel
       read (b,*,err=3006) (segWt(iFlt0,k),k=1,nSegModel)

c      Read total number of fault segments
       read (b,*,err=3007) nFlt2
      
       do i=1,nSegModel
         read (b,*,err=3008) (faultFlag(iFlt0,i,k),k=1,nFlt2)
       enddo

       do iflt2=1,nflt2
        iFlt = iFlt + 1
        call S21_CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c       Set fault indexes
        fIndex(1,iflt) = iflt0
        fIndex(2,iflt) = iflt2

c       Find total weight for this segment
        sum = 0.
        do i=1,nSegModel
          sum = sum + segWt(iFlt0,i) * faultFlag(iFlt0,i,iFlt2)
        enddo
        segWt1(iFlt) = sum

c       Set segment flags for this fault system
        sum = 0.
        do i=1,nSegModel
          segModelFlag(iFlt,i) = faultFlag(iFlt0,i,iFlt2)
          segModelWt1(iFlt,i) = segWt(iFlt0,i)
        enddo
        nSegModel0(iFlt) = nSegModel

c       Read name of this segment
        read(b,'( a80)',err=3009) fName(iFlt)
        write (*,'( i5,2x,a80)') iFlt, fName(iFlt)

c       Read type of source (planar, areal, grid1, grid2, irregular)
        read(b,*,err=3010) sourceType(iFlt),attenType(iFlt),sampleStep(iFlt),
     1                 fltDirect(iFlt), synchron(iFlT)

c       Now read in the synchronous rupture case parameters if needed.
        if (synchron(iFlt) .gt. 0) then
           read (b,*,err=3011) nsyn_Case(iFlt), synjcalc(iFlt)

c          Now read in the magnitude, distance and weigths for each synchronous rupture case.
           do isyn=1,nsyn_Case(iFlt)
             read (b,*,err=3012) synmag(iFlt,isyn),syndistRup(iFlt,isyn),
     1                 syndistJB(iFlt,isyn),syndistSeismo(iFlt,isyn),
     2                 synHWFlag(iFlt,isyn),synhypo(iFlt,isyn),
     3                 synftype(iFlt,isyn), 
     4                 syn_dip(iFlt,isyn), syn_zTOR(iFlt,isyn), 
     5                 syn_RupWidth(iFlt,isyn), syn_RX(iFlt,isyn), syn_Ry0(iFlt,isyn),
     5                 synwt(iFlt,isyn)
           enddo
        endif

c       Read aleatory segmentation wts
        read (b,*,err=3013) al_segWt(iFlt)

c       Check for standard fault source or areal source
        if ( sourceType(iFlt) .eq. 1 .or. sourceType(iFLt) .eq. 2) then
          read (b,*,err=3014) dip1, top
          read(b,*,err=3015) nfp(iFlt)

          call S21_CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
          do ipt=1,nfp(iFlt)
            read (b,*,err=3016) fLong(iFlt,1,ipt), fLat(iFlt,1,ipt)
            fZ(iFlt,1,ipt) = top
            if (sourceType(iFlt) .eq. 2) then
               grid_top(iFlt,ipt) = top
            endif
          enddo
          nDownDip(iFlt) = 1
        endif

c        Check for grid source (w/o depth)
         if ( sourceType(iFlt) .eq. 3 ) then
           read (b,*,err=3014)  dip1, top

           call S30_RdGrid1 (iFlt, grid_a, grid_dlong, grid_dlat, grid_n,
     1             grid_long, grid_lat, minLat, minLong,  maxLat, maxLong, scaleRate(iFlt) )

           do igrid=1,grid_n(iFlt)
             grid_top(iFlt,igrid) = top
           enddo
           nDownDip(iFlt) = 1
         endif

c        Check for grid source (w/ depth)
         if ( sourceType(iFlt) .eq. 4 ) then
              read (b,*,err=3014)  dip1

              call S30_RdGrid2 (iFlt,grid_a,grid_dlong,grid_dlat,grid_n,
     1             grid_long, grid_lat, minLat, minLong, maxLat, maxLong, scaleRate(iFlt),
     2             grid_top)

              nDownDip(iFlt) = 1
         endif

c        Check for custom fault source
         if ( sourceType(iFlt) .eq. 5 .or. sourceType(iFlt) .eq. 6) then
           read(b,*,err=3017) nDownDip(iFLt), nfp(iFlt)  
   
           call S21_CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
           do ipt=1,nfp(iFlt)
              read (b,*,err=3018) (fLong(iFlt,k,ipt), fLat(iFlt,k,ipt), fZ(iflt,k,ipt), k=1,nDownDip(iFlt) ) 

          enddo
         endif

c        Check for Sourcetype 7 (UCERF)
         if ( sourceType(iFlt) .eq. 7 ) then
            read (b,*)  top, HWsource7
            call S30_RdSource7 (iFlt, mags7, rates7, dists7, dips7, mechs7, ncountS7 )
            write (*,*) 'NcountS7 =', ncountS7(iFlt)
         endif

c        Read dip Variation
         if ( sourceType(iFlt) .lt. 5 ) then
           read (b,*,err=3019) n_Dip
           call S21_CheckDim ( n_Dip, MAX_N1, 'MAX_N1    ' )
           read (b,*,err=3020) (deltaDip1(i),i=1,n_Dip)
           read (b,*,err=3021) (dipWt1(i),i=1,n_Dip)
           call S21_CheckWt ( dipWt1, n_Dip, fName(iFlt), 'Dip                 ' )
         else
           n_Dip = 1
           deltaDip1(1) = 0.
           dipWt1(1) = 1.
         endif

c        Read b-values (not for activity rate cases)
         read (b,*,err=3022) n_bValue
         n_bValue3(iflt)=n_bValue
         call S21_CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
         if ( n_bValue .gt. 0 ) then
           read (b,*,err=3023) (bValue1(i),i=1,n_bValue)
           read (b,*,err=3024) (bValueWt1(i),i=1,n_bValue)
           do i=1,n_bValue
             bValue3(iflt,i) = bValue1(i)
             bValueWt3(iflt,i) = bValueWt1(i)
           enddo
         endif

c        Read activity rate - b-value pairs
         read (b,*,err=3025) nActRate 

         if ( nActRate .ne. 0 ) then
           do ii=1,nActRate
             read (b,*,err=3026) bValue2(ii), actRate(ii), actRateWt1(ii)
C     Scale activity rate for gridded source (i.e., sourcetype 3 or 4) based on limited lat and long values.
             if (sourcetype(iFlt) .eq. 3 .or. sourcetype(iFlt) .eq. 4) then
                if (scalerate(iFlt) .ne. 1.0) then
                   actRate(ii) = actRate(ii)*Scalerate(iFlt)
                endif
             endif
           enddo

           call S21_CheckWt ( actRateWt1, nActRate, fName(iFlt), 'ActRate             ' )
         endif

c        Read weights for rate methods
         read (b,*,err=3027) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
         sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
         if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFlt)
              stop 99
         endif

c        Read slip-rates
         read (b,*,err=3028) nSR
         nSR1(iflt)=nSR
         if ( nSR .gt. 0 ) then
           read (b,*,err=3029) (sr(k),k=1,nSR)
           read (b,*,err=3030) (wt_sr(k),k=1,nSR)
           do k=1,nSR
             sr1(iflt,k) = sr(k)
             wt_sr1(iflt,k) = wt_sr(k)
           enddo
           call S21_CheckWt (wt_sr, nSR, fName(iFlt), 'Slip Rates          ')
         endif

c        Read recurrence intervals
         read (b,*,err=3031) nRecInt
         if ( nRecInt .gt. 0 ) then
           read (b,*,err=3032) (rec_Int(k),k=1,nRecInt)
           read (b,*,err=3033) (wt_recInt(k),k=1,nRecInt)
           call S21_CheckWt (wt_recInt, nRecInt, fName(iFlt), 'Recurrence Intervals')
         endif

c        Read moment-rates
         read (b,*,err=3034) nMoRate
         if ( nMoRate .gt. 0 ) then
           read (b,*,err=3035) (MoRate(k),k=1,nMoRate)
           read (b,*,err=3036) (MoRateDepth(k),k=1,nMoRate)
           read (b,*,err=3037) (wt_MoRate(k),k=1,nMoRate)
           call S21_CheckWt (wt_MoRate, nMoRate, fName(iFlt), 'Moment Rates        ')
         endif

c        Check that slip-rates are not used for areal sources
         if ( sourceType(iFlt) .ge. 2 .and. sourceType(iFlt) .le. 4 ) then               
           if (nSR .ne. 0) then
              write (*,'( 2x,''Error: slip-rates not allowed for areal sources'')')
              write (*,'( 2x,''Flt: '',a60)') fname(iFlt)
              stop 99
           endif
         endif

c        Load into single array called "rate_param"
c        Set MoRDepth=1.0 or the inverse for latter scaling of MoRates
         nRate = nSR + nActRate + nRecInt + nMoRate
         call S21_CheckDim ( nRate, MAX_N1, 'MAX_N1    ' )
         do k=1,nSR
            rateParam1(k) = sr(k)
            rateWt1(k) = wt_sr(k)*wt_srbranch              
            rateType1(k) = 1
            MoRDepth(k) = 1.0
         enddo
         do k=1,nActRate
            rateParam1(k+nSR) = actRate(k)
            rateWt1(k+nSR) = actRateWt1(k) * wt_actRateBranch             
            rateType1(k+nSR) = 2
            MoRDepth(k+nSR) = 1.0
         enddo
         do k=1,nRecInt
            rateParam1(k+nSR+nActRate) = rec_Int(k)
            rateWt1(k+nSR+nActRate) = wt_recInt(k) *wt_recIntBranch              
            rateType1(k+nSR+nActRate) = 3
            MoRDepth(k+nSR+nActRate) = 1.0
         enddo
         do k=1,nMoRate
            rateParam1(k+nSR+nActRate+nRecInt) = MoRate(k)
            rateWt1(k+nSR+nActRate+nRecInt) = wt_MoRate(k) *wt_MoRateBranch              
            rateType1(k+nSR+nActRate+nRecInt) = 4
            MoRDepth(k+nSR+nActRate+nRecInt) = 1.0/MoRateDepth(nMoRate)
         enddo

c        Read Mag recurrence weights (char, exp, etc.)
         read (b,*,err=3038) nMagRecur
         call S21_CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
         read (b,*,err=3039) (magRecur1(i),i=1,nMagRecur)
         read (b,*,err=3040) (magRecurWt1(i),i=1,nMagRecur)
         call S21_CheckWt ( magRecurWt1,nMagRecur, fName(iFlt), 'Mag Recur           ' )

c        Read in corresponding magnitude parameters. 
         do iRecur=1,nMagRecur
           if (magRecur1(iRecur) .eq. 4) then
             read (b,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
             rp6(iRecur) = 0.0
c          Read in necessary parameters for bi-exponential distribution.
           elseif (magRecur1(iRecur) .eq. 5) then
             read (b,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
             rp6(iRecur) = 0.0
           elseif (magRecur1(iRecur) .eq. 10) then
             read (b,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur),
     1       rp6(iRecur)
           else              
             read (b,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur)
             rp4(iRecur) = 0.0
             rp5(iRecur) = 0.0
             rp6(iRecur) = 0.0
           endif
         enddo

c        Read seismogenic thickness
         if ( sourceType(iFlt) .lt. 5) then
           read (b,*,err=3042) nThick1
           call S21_CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
           read (b,*,err=3043) (faultThick1(i),i=1,nThick1)
           read (b,*,err=3044) (faultThickWt1(i),i=1,nThick1)         
           call S21_CheckWt (faultThickWt1, nThick1, fName(iFlt), 'Seismo Thick        ')
         else
           nThick1 = 1
           faultTHick1(1) = -99.
           faultThickWt1(1) = 1.
         endif
         
c        Read depth pdf
         read (b,*,err=3045) iDepthModel(iFlt), (depthParam(iflt,k),k=1,3)     

c        Read reference mags for each fault thickness
         iThickDip = 1
         do iThick1=1,nThick1
           read (b,*,err=3047) nRefMag0
           read (b,*,err=3048) (refMag0(i),i=1,nRefMag0)
           read (b,*,err=3049) (refMagWt0(i),i=1,nRefMag0)

           if ( nRefMag0 .ne. 0 ) then
             call S21_CheckWt ( refMagWt0, nrefMag0, fName(iFlt), 'Ref Mag                 ' )
           endif       
              
c          Copy these ref magnitudes for each addtional dip (no correlation allowed)
           do iDip=1,n_Dip
             nRefMag(iThickDip) = nRefMag0
             nRefMag2(iFlt,iThickDip) = nRefMag(iThickDip)
             do i=1,nRefMag0
               refMag1(iThickDip,i) = refMag0(i)
               refMag2(iFlt,iThickDip,i) = refMag1(iThickDip,i)
               refMagWt1(iThickDip,i) = refMagWt0(i)
               refMagWt2(iFlt,iThickDip,i) = refMagWt1(iThickDip,i) 
             enddo
             iThickDip = iThickDip + 1
           enddo
         enddo 

c        Read min mag, step sizes, and rupture variability info
         read (b,*,err=3050) minMag(iFlt), magStep(iFlt), hxStep(iFlt), 
     1              hyStep(iFlt), nRupArea(iFlt), nRupWidth(iFlt), minDepth(iFlt)

C       Allow hxstep and hystep to be different for sourceType 2
        if (sourceType(iFlt).ne.2) then
C         Check that Hxstep = Hystep
          if (hxstep(iFlt) .ne. hystep(iFlt) ) then
            write (*,*) 'Hxstep not equal to Hystep for fault: '
            write (*,*) fname1
            write (*,*) 'Hxstep = ', hxstep(iFlt) 
            write (*,*) 'Hystep = ', hystep(iFlt) 
            write (*,*) 'These values must be equal or errors can occur!'
            write (*,*) 'Check input fault file.'
            stop 99
          endif
        endif

         read (b,*,err=3051) (coef_area(k,iFlt),k=1,2), sigArea(iFlt)
         read (b,*,err=3052) (coef_width(k,iFlt),k=1,2), sigWidth(iFlt)

c        Read ftype Models
         read (b,*,err=3053) nFtypeModels
         do iFM=1,nFtypeModels
           read (b,*,err=3054) ftmodelwt(iFM)
           read (b,*,err=3055) nFtype1(iFM)
           read (b,*,err=3056) (ftype1(iFM,k),k=1,nFtype1(iFM))
           read (b,*,err=3057) (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
           call S21_CheckWt1a (Ftype_Wt1, nFtype1(iFM), iFM, 
     1          MAX_N1, fName(iFlt), 'Fault Mech          ' )
         enddo

c        Load up Ftype Models and weights.
         nFM = 1
         do iFM=1,nFtypeModels
           do k=1,nFtype1(iFM)
             ftype(iFlt,nFM) = ftype1(iFM,k)
             ftype_wt(iFlt,nFM) = ftype_wt1(iFM,k)*ftmodelwt(iFM)
             nFm = nFm + 1
           enddo
         enddo
         nFm = nFm - 1
         nFtype(iFlt) = nFm
         call S21_CheckDim ( nFtype(iFlt), MAX_FTYPE, 'MAX_FTYPE ' )

c        Load up parameter variations into large single dimension arrays
         testMaxMag = 0.
         iWidth = 0
         do iThick1=1,nThick1
        
           do iDip=1,n_Dip
           iWidth = iWidth + 1
           call S21_CheckDim ( iWidth, MAX_WIDTH, 'MAX_WIDTH' )   
           
             dip2 = dip1 + deltaDip1(iDip)
             dip(iFlt,iWidth,1) = dip2
             faultWidth(iFlt,iWidth) = faultThick1(iThick1)
             faultWidthWt(iFlt,iWidth) = faultThickWt1(iThick1) * dipWt1(iDip)

           i = 0
           mtest = 0.0

           do iRecur=1,nMagRecur    
           
            do iRate=1,nRate
              if ( rateType1(iRate) .eq. 2 ) then
                nb1 = 1
              else
                nb1 = n_bvalue
              endif
              do i_bValue=1,nb1
               do iRefMag=1,nRefMag(iWidth)
                 i = i + 1
                 call S21_CheckDim ( i, MAXPARAM, 'MAXPARAM  ' )
                 magRecur(iFlt,i,iWidth) = magRecur1(iRecur)
                 magRecurWt(iFlt,i,iWidth) = magRecurWt1(iRecur)

c                Scale moment rate by reference thickness.
                 if (rateType1(iRate) .eq. 4) then
                    rateParam(iFlt,i,iWidth) = rateParam1(iRate)*MoRDepth(iRate)*faultWidth(iFlt,iWidth)
                 else
                    rateParam(iFlt,i,iWidth) = rateParam1(iRate)
                 endif

                 rateType(iFlt,i,iWidth) = rateType1(iRate)
                 if ( rateType1(iRate) .eq. 2 ) then
                   beta(iFlt,i,iWidth) = bValue2(iRate - nSR)*alog(10.)
                 else
                   beta(iFlt,i,iWidth) = bValue1(i_bValue)*alog(10.)
                 endif
                 if ( rateType1(iRate) .eq. 2 ) then
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate)
                 else
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate) * bValueWt1(i_bValue) 
                 endif
                 maxMagWt(iFlt,i,iWidth) = refMagWt1(iWidth,iRefMag)

c                Set max mag
                 if (magRecur1(iRecur) .eq. 0 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 1 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP1(iRecur)
                 elseif (magRecur1(iRecur) .eq. 3 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 6 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 elseif (magRecur1(iRecur) .eq. 10 ) then
                   if ( rp5(iRecur) .eq. 0. ) then
                     maxMag(iFlt,i,iWidth) = rp2(iRecur)
                   else
                     maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rp5(iRecur)
                   endif
                 endif
                 if ( maxMag(iFlt,i,iWidth) .gt. testMaxMag ) then
                   nMag(iFlt) = ceiling((maxMag(iFlt,i,iWidth) - minMag(iFLt) ) / magStep(iFlt))
                   if (sourceType(iFlt) .eq. 7) then
                     nMag(iFlt) = ncountS7(iFlt)
                   endif
                     if (nMag(iFlt) .eq. 0 ) then
                       nMag(iFlt) = 1
                     endif    
                   testMaxMag = maxMag(iFlt,i,iWidth)
                 endif               

c  temp fix - the maxmag, refmag problem for waacy
                 if (magRecur1(iRecur) .eq. 10 ) then
                     maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 endif
                          
                 mpdf_param(iFlt,i,iWidth,1) = rP1(iRecur)
                 mpdf_param(iFlt,i,iWidth,2) = rP2(iRecur)
                 mpdf_param(iFlt,i,iWidth,3) = rP3(iRecur)
                 mpdf_param(iFlt,i,iWidth,4) = rP4(iRecur)
                 mpdf_param(iFlt,i,iWidth,5) = rP5(iRecur)
                 mpdf_param(iFlt,i,iWidth,6) = rP6(iRecur)
                 
                 mpdf_param1(iFlt,iWidth,1) = rP1(iRecur)
                 mpdf_param1(iFlt,iWidth,2) = rP2(iRecur)
                 mpdf_param1(iFlt,iWidth,3) = rP3(iRecur)
                 mpdf_param1(iFlt,iWidth,4) = rP4(iRecur)
                 mpdf_param1(iFlt,iWidth,5) = rP5(iRecur)
                 mpdf_param1(iFlt,iWidth,6) = rP6(iRecur)

c      Test for the maximum magnitude value for all realizations of this fault.
                 if (maxMag(iFlt,i,iWidth) .ge. mtest) then
                   mtest = maxMag(iFlt,i,iWidth)
                 endif

                 mmagout(iFlt,iWidth,iRefMag) = refMag1(iWidth,iRefMag)
                 mmagoutwt(iFlt,iWidth,iRefMag) = refMagWt1(iWidth,iRefMag)
     
c     End loop over iRefMag
                enddo
c     End loop over ib_value
               enddo
c     End loop over iRate
              enddo
c     End loop over iRecur
             enddo
             nParamVar(iFlt,iWidth) = i
c     End loop over iDip
            enddo
c     End loop over iThick1
           enddo
           probAct(iFlt) = probAct0
           nWidth(iFlt) = iWidth
           fsys(iFlt) = ifsystem
c     End loop over iFlt2 (# segments)
          enddo
          ifsystem = ifsystem + 1          
c     End loop over iFlt
      enddo
      nfsystem = ifsystem - 1
      nFlt = iflt
      
      close (b)
      
      return
 3001 write (*,'( 2x,''Flt file error:  iCorr'')')
      stop 99
 3002 write (*,'( 2x,''Flt file error:  nFlt'')')
      stop 99
 3003 write (*,'( 2x,''Flt file error:  flt sys name'')')
      stop 99
 3004 write (*,'( 2x,''Flt file error:  prob Act'')')
      stop 99
 3005 write (*,'( 2x,''Flt file error:  nSegModel'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3006 write (*,'( 2x,''Flt file error:  Seg wts'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3007 write (*,'( 2x,''Flt file error:  number of segments'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3008  write (*,'( 2x,''Flt file error: seg flags for seg model'', i5)') i
       write (*,'( 2x,''From fault: '',a80)') fName1
       stop 99
 3009 write (*,'( 2x,''Flt file error: seg name for seg model '', i5)') i
      write (*,'( 2x,''From fault: '',a80)') fName1
      stop 99
 3010 write (*,'( 2x,''Flt file error: source type line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3011 write (*,'( 2x,''Flt file error: nSynRup line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3012 write (*,'( 2x,''Flt file error: syn Rup param line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3013 write (*,'( 2x,''Flt file error: Aleatory seg wt '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3014 write (*,'( 2x,''Flt file error: dip and top '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3015 write (*,'( 2x,''Flt file error: number flt point '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3016 write (*,'( 2x,''Flt file error: long lat values '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3017 write (*,'( 2x,''Flt file error: nDowndip, nFP '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3018 write (*,'( 2x,''Flt file error: long lat values '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3019 write (*,'( 2x,''Flt file error: nDip '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3020 write (*,'( 2x,''Flt file error: delta Dip '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3021 write (*,'( 2x,''Flt file error: Dip wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3022 write (*,'( 2x,''Flt file error: n b-value '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3023 write (*,'( 2x,''Flt file error: delta b-value '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3024 write (*,'( 2x,''Flt file error: b-value wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3025 write (*,'( 2x,''Flt file error: nActRate'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3026 write (*,'( 2x,''Flt file error: b, activity, wt'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3027 write (*,'( 2x,''Flt file error: wts for activity rate approach'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3028 write (*,'( 2x,''Flt file error: number Slip rates '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3029 write (*,'( 2x,''Flt file error: slip rates '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3030 write (*,'( 2x,''Flt file error: slip-rate wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3031 write (*,'( 2x,''Flt file error: number Rec Intervals '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3032 write (*,'( 2x,''Flt file error: Rec Intervals '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3033 write (*,'( 2x,''Flt file error: Rec Interval wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3034 write (*,'( 2x,''Flt file error: Number Moment rate'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3035 write (*,'( 2x,''Flt file error: Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3036 write (*,'( 2x,''Flt file error: depths for Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3037 write (*,'( 2x,''Flt file error: wts for Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3038 write (*,'( 2x,''Flt file error: number mag pdf'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3039 write (*,'( 2x,''Flt file error: mag pdf flag'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3040 write (*,'( 2x,''Flt file error: mag pdf wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3041 write (*,'( 2x,''Flt file error: param for mag pdf'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3042 write (*,'( 2x,''Flt file error: number crustal thickness '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3043 write (*,'( 2x,''Flt file error: crustal thicknesses '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3044 write (*,'( 2x,''Flt file error: crustal thickness wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3045 write (*,'( 2x,''Flt file error: depth pdf and param'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3047 write (*,'( 2x,''Flt file error: depth pdf and param'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3048 write (*,'( 2x,''Flt file error: number Ref Mag'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3049 write (*,'( 2x,''Flt file error: Ref mags'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3050 write (*,'( 2x,''Flt file error: Ref mag wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3051 write (*,'( 2x,''Flt file error: coeff A(M) model'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3052 write (*,'( 2x,''Flt file error: coeff W(M) model'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3053 write (*,'( 2x,''Flt file error: number Ftype models'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3054 write (*,'( 2x,''Flt file error: Ftype mode wt'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3055 write (*,'( 2x,''Flt file error: number Ftype'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3056 write (*,'( 2x,''Flt file error: Ftypes'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3057 write (*,'( 2x,''Flt file error: Ftype Aleatory wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99

      end

c ----------------------------

      subroutine S21_CheckDim ( n, nMax, name )

      implicit none
      
      character(len=*) :: name
      integer n, nMax
      
      if ( n .gt. nMax ) then
        write (*,*)
        write (*,'( 2x,''Array Dimension Too Small'')')
        write (*,'( 2x,''Increase '',a10,''to '',i5)') name, n
        write (*,*)
        stop 99
      endif
      return
      end

c  --------------------------

      subroutine S21_CheckWt ( x, n, fName, name )

      implicit none
      include 'max_dims.h'
      
      integer i, n
      real x(MAXPARAM), sum
      character*20 name, fName
      
      sum = 0.
      do i=1,n
        sum = sum + x(i)
      enddo
      sum = sum - 1.
      if ( abs(sum) .gt. 0.001 ) then
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 2x,a11)') name
        write (*,'( i5, 10f10.4)') n, (x(i),i=1,n)
        write (*,'( 2x,''sum ='',f12.7)') sum
        write (*,'( 2x,a80)') fName
        stop 99
      endif
      return
      end

c  --------------------------

      subroutine S21_CheckWt1a ( x, n, j, n1, fName, name  )

      implicit none
      include 'max_dims.h'
      
      integer n, j, n1, i
      real x(MAX_FLT,MAX_N1), sum
      character*20 fName, name

      sum = 0.
      do i=1,n
        sum = sum + x(j,i)
      enddo
      if ( abs(sum-1.)  .gt. 0.001 ) then
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 10f10.4)') (x(j,i),i=1,n),sum
        write (*,'( 2x,a20)') name
        write (*,'( 2x,a80)') fName
        stop 99
      endif
      return
      end

c  -------------------------------------------------------------------

       subroutine S30_RdGrid1 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale )

      implicit none
      include 'max_dims.h'
         
      integer grid_n(MAX_FLT), iFlt, nHead, i, n, j
      real grid_long(MAX_FLT,MAX_GRID), grid_lat(MAX_FLT,MAX_GRID),
     1     grid_a(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT), 
     2     grid_dlat(MAX_FLT), minLat, minLong, maxLat, maxLong, rate_scale     
      real*8 sum, sum1  
      character*80 filein, dummy  
                
      read (10,'( a80)') filein
      open (11,file=filein,status='old',err=2100)
          
c     Read header
      read (11,*,err=2001) nHead
      write (*,'( i5)') nHead
      
      do i=1,nHead
        read (11,'( a80)',err=2002) dummy
        write (*,'( a80)') dummy
      enddo
      read (11,*,err=2003) n, grid_dlong(iFlt), grid_dlat(iFlt)
      if ( n .gt. MAX_GRID ) then
        write (*,'( 2x,''Increase MAX_GRID to '',i7)') n
        stop 99
      endif
      write (*,'( i5)') n
          
c     read grid data keeping only grid points with non-zero rates
c     and between min and max long and lat
      j = 1
      sum = 0.
      sum1 = 0.
      do i=1,n
        read (11,*,err=2004) grid_long(iFlt,j), grid_lat(iFlt,j),grid_a(iFlt,j)
        sum = sum + grid_a(iFlt,j)
        if ( grid_a(iFlt,j) .gt. 0. 
     1       .and. grid_long(iFlt,j) .ge. minLong  
     1       .and. grid_long(iFlt,j) .le. maxLong 
     2       .and. grid_lat(iFlt,j) .ge. minLat 
     2       .and. grid_lat(iFlt,j).le. maxLat ) then
          sum1 = sum1 + grid_a(iFlt,j)
          j = j + 1
        endif
      enddo
      close (11)
      grid_n(iFlt) = j - 1
      rate_scale = sum1/sum

      
c     Note: Rate_scale keeps track of the activity rate that is removed because
c     it is at too large a distance (but is still in the activity rate of the input)            

      return
      
 2100 write (*,'( 2x,''bad gridded seismicity file'')')
      write (*,'( a80)') filein
      stop 99
 2001 write (*,'( 2x,''gridded seismicity file error: nHeader'')')
      stop 99
 2002 write (*,'( 2x,''gridded seismicity file error: header'')')
      stop 99
 2003 write (*,'( 2x,''gridded seismicity file error: n, grid_dLong, grid_dlat'')')
      stop 99
 2004 write (*,'( 2x,''gridded seismicity file error: long, lat, a-value'')')
      write (*,'( 2x,''entry:'',i5)') i
      stop 99
      end
      
c ----------------------------------------------------------------------
     
      subroutine S30_RdGrid2 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale,
     2           grid_top )

      implicit none
      include 'max_dims.h'

      integer grid_n(1), iFlt, nHead, i, n, j
      real grid_long(MAX_FLT,1), grid_lat(MAX_FLT,1), grid_a(MAX_FLT,1),
     1     grid_dlong(1), grid_dlat(1), grid_top(MAX_FLT,1), minLat, 
     2     minLong, maxLat, maxLong, dummy, rate_scale
      real*8 sum, sum1    
      character*80 filein                     

      read (10,'( a80)') filein
      open (11,file=filein,status='old')
          
c     Read header
      read (11,*) nHead
      do i=1,nHead
        read (11,'( a1)') dummy
      enddo
      read (11,*) n, grid_dlong(iFlt), grid_dlat(iFlt)
      if ( n .gt. MAX_GRID ) then
        write (*,'( 2x,''Increase MAX_GRID to '',i7)') n
        stop 99
      endif
          
c     read grid data keeping only grid points with non-zero rates
c     and between min and max long and lat
      j = 1
      sum = 0.
      sum1 = 0.
      do i=1,n
        read (11,*) grid_long(iFlt,j), grid_lat(iFlt,j),grid_top(iFlt,j), grid_a(iFlt,j)
        sum = sum + grid_a(iFlt,j)
        if ( grid_a(iFlt,j) .gt. 0. 
     1       .and. grid_long(iFlt,j) .ge. minLong  
     1       .and. grid_long(iFlt,j) .le. maxLong 
     2       .and. grid_lat(iFlt,j) .ge. minLat 
     2       .and. grid_lat(iFlt,j).le. maxLat ) then
          sum1 = sum1 + grid_a(iFlt,j)
          j = j + 1
        endif
      enddo
      close (11)
      grid_n(iFlt) = j - 1
      rate_scale = sum1/sum
      
c     Note: Rate_scale keeps track of the activity rate that is removed because
c     it is at too large a distance (but is still in the activity rate of the input)            

      return
      end

c  --------------------------------------------------------------------

      subroutine S30_RdSource7 ( iFlt, mag, rate, dist, dip, mech, ncount )
 
      implicit none
      include 'max_dims.h'
      
      integer ncount(MAX_FLT), srflag, rupid, nearID, i, iFlt         
      real rate(MAX_FLT,MAX_S7), mag(MAX_FLT,MAX_S7), dist(MAX_FLT,MAX_S7), 
     1     lat, long, strike, Dip(MAX_FLT,MAX_S7), rake, mech(MAX_FLT,MAX_S7)
      character*80 filein
                
      read (10,'( a80)') filein
      open (11,file=filein,status='old')
          
      ncount(iflt) = 0
      do i=1,MAX_S7
         read (11,*,END=777) srflag,rupID,mag(iFlt,i),rate(iFlt,i),nearID,dist(iFlt,i),
     1             lat,long,strike,Dip(iFlt,i),rake
     
         if (rake .ge. 30.0 .and. rake .le. 150.0) then    
           mech(iFlt,i) = 1.0
         elseif (rake .ge. -120.0 .and. rake .le. -60.0) then
           mech(iFlt,i) = -1.0
         else
           mech(iFLt,i) = 0.0
         endif
      
         ncount(iFlt) = ncount(iFlt) + 1
      enddo

      write (*,*) 'More than 70,000 source in data file.'
      write (*,*) 'Need to separate into two files.'
      stop 99
      
 777  continue

C     Classify the Rake angle with Mechanims choices of SS(0), RV(1), or NM(-1)
C     The following Rake Angles are used to classify Fault Mechanims:
C     Reverse: 30<=Rake<=150  (includes RV/OB as RV)
C     Normal:  -120<=Rake<=-60
C     Strike-Slip: -180<=Rake<-120
C                  -60<Rake<=30
C                  150<Rake<=180
C        (includes NM/OB as SS)

      return
      end

     

