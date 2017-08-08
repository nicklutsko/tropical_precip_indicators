load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin 


file_names   =new (20,string)
file_names_A =new (20,string)

file_names(0)   ="Amon_ACCESS1-0_historical_r1i1p1_197901-200512"
file_names_A(0) ="Amon_ACCESS1-0_amip_r1i1p1_197901-200512"

file_names(1)   ="Amon_ACCESS1-3_historical_r1i1p1_197901-200512"
file_names_A(1) ="Amon_ACCESS1-3_amip_r1i1p1_197901-200512"

file_names(2)   ="Amon_bcc-csm1-1_historical_r1i1p1_197901-200512"
file_names_A(2) ="Amon_bcc-csm1-1_amip_r1i1p1_197901-200512"

file_names(3)   ="Amon_BNU-ESM_historical_r1i1p1_197901-200512"
file_names_A(3) ="Amon_BNU-ESM_amip_r1i1p1_197901-200512"

file_names(4)   ="Amon_CanESM2_historical_r1i1p1_197901-200512"
file_names_A(4) ="Amon_CanAM4_amip_r1i1p1_197901-200512"

file_names(5)   ="Amon_CCSM4_historical_r1i1p1_197901-200512"
file_names_A(5) ="Amon_CCSM4_amip_r1i1p1_197901-200512"

file_names(6)   ="Amon_CNRM-CM5_historical_r1i1p1_197901-200512"
file_names_A(6) ="Amon_CNRM-CM5_amip_r1i1p1_197901-200512"

file_names(7)   ="Amon_FGOALS-s2_historical_r1i1p1_197901-200512"
file_names_A(7) ="Amon_FGOALS-s2_amip_r1i1p1_197901-200512"

file_names(8)   ="Amon_GFDL-CM3_form_historical_r1i1p1_197901-200512"
file_names_A(8) ="Amon_GFDL-CM3_form_amip_r1i1p1_197901-200512"

file_names(9)   ="Amon_GISS-E2-R_historical_r1i1p1_197901-200512"
file_names_A(9) ="Amon_GISS-E2-R_amip_r1i1p1_197901-200512"

file_names(10)   ="Amon_HadGEM2-ES_historical_r1i1p1_197901-200512"
file_names_A(10) ="Amon_HadGEM2-A_amip_r1i1p1_197901-200512"

file_names(11)   ="Amon_inmcm4_historical_r1i1p1_197901-200512"
file_names_A(11) ="Amon_inmcm4_amip_r1i1p1_197901-200512"

file_names(12)   ="Amon_IPSL-CM5A-LR_historical_r1i1p1_197901-200512"
file_names_A(12) ="Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200512"

file_names(13)   ="Amon_IPSL-CM5A-MR_historical_r1i1p1_197901-200512"
file_names_A(13) ="Amon_IPSL-CM5A-MR_amip_r1i1p1_197901-200512"

file_names(14)   ="Amon_IPSL-CM5B-LR_historical_r1i1p1_197901-200512"
file_names_A(14) ="Amon_IPSL-CM5B-LR_amip_r1i1p1_197901-200512"

file_names(15)   ="Amon_MIROC5_historical_r1i1p1_197901-200512"
file_names_A(15) ="Amon_MIROC5_amip_r1i1p1_197901-200512"

file_names(16)   ="Amon_MPI-ESM-LR_historical_r1i1p1_197901-200512"
file_names_A(16) ="Amon_MPI-ESM-LR_amip_r1i1p1_197901-200512"

file_names(17)   ="Amon_MPI-ESM-MR_historical_r1i1p1_197901-200512"
file_names_A(17) ="Amon_MPI-ESM-MR_amip_r1i1p1_197901-200512"

file_names(18)   ="Amon_MRI-CGCM3_historical_r1i1p1_197901-200512"
file_names_A(18) ="Amon_MRI-CGCM3_amip_r1i1p1_197901-200512"

file_names(19)   ="Amon_NorESM1-M_historical_r1i1p1_197901-200512"
file_names_A(19) ="Amon_NorESM1-M_amip_r1i1p1_197901-200512"

print("AMIP:")
print(file_names_A)

print("CMIP:")
print(file_names)

m_number=dimsizes(file_names)
v_number=dimsizes(file_names)


  pr_W        = new((/2,m_number/),float)
  pr_E        = new((/2,m_number/),float)
  pr_A        = new((/2,m_number/),float)
  pr_S        = new((/2,m_number/),float)
  pr_I        = new((/2,m_number/),float)
  pr_D        = new((/2,m_number/),float)
  pr_H        = new((/2,m_number/),float)
  pr_L        = new((/2,m_number/),float)
  pr_R        = new((/2,m_number/),float)

  pr_Ai        = new((/2,m_number/),float)
  pr_Ei        = new((/2,m_number/),float)
  pr_Mi        = new((/2,m_number/),float)

 PI=acos(-1.)

trop_lat = 30.
crit_lat=20.

index_lat=20.
eq_lat = 2.


;-----------------------------------------------------------------------------------------------------------------------------------------------
; AMIP:

do j_f=1,m_number

 file1=addfile("../data_1979_2005/merged_files/zonmean/"+file_names_A(j_f-1)+"_tm_ZM.nc","r")

  lat=file1->lat
  lat_number = dimsizes(lat)
  precip=86400.*file1->pr

  l_0 = 0

  do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -trop_lat .AND. lat(j_lat-2) .LT. -trop_lat) then 
    l_ST = j_lat-1
   end if

   if(lat(j_lat-1) .GE. -crit_lat .AND. lat(j_lat-2) .LT. -crit_lat) then 
    c_ST = j_lat-1
   end if

   if(lat(j_lat-1) .GE. -index_lat .AND. lat(j_lat-2) .LT. -index_lat) then 
    i_ST = j_lat-1
   end if

   if(lat(j_lat-1) .GE. -eq_lat .AND. lat(j_lat-2) .LT. -eq_lat) then 
    i_Seq = j_lat-1
   end if

   if(lat(j_lat-1) .GE. 0. .AND. lat(j_lat-2) .LT. 0.) then
     if (lat(j_lat-1) .EQ. 0.) then
       l_0 = j_lat-1
       l_0N = j_lat-1
       l_0S = j_lat-2 
     else
       l_0N = j_lat-1
       l_0S = j_lat-2 
     end if
   end if

   if(lat(j_lat) .GE. trop_lat .AND. lat(j_lat-1) .LT. trop_lat) then 
     l_NT = j_lat-1
   end if

   if(lat(j_lat) .GE. crit_lat .AND. lat(j_lat-1) .LT. crit_lat) then 
    c_NT = j_lat-1
   end if

   if(lat(j_lat) .GE. index_lat .AND. lat(j_lat-1) .LT. index_lat) then 
    i_NT = j_lat-1
   end if

   if(lat(j_lat) .GE. eq_lat .AND. lat(j_lat-1) .LT. eq_lat) then 
    i_Neq = j_lat-1
   end if
 
   end do

;--------------------------------------------------------------------------------
; Find absolute precip. maximum
  
  
  ITCZ=c_ST + maxind(precip(0,c_ST:c_NT,0))
  pr_I(0,j_f-1) = abs(doubletofloat(lat(ITCZ)))

;--------------------------------------------------------------------------------
; Find second ITCZ (if exists)

  s_max_ind = 0

    if (lat(ITCZ) .GT. 0. ) then
    do j_lat= 1,(l_0N-c_ST)
     if ( precip(0,l_0N-j_lat,0) .gt. precip(0,l_0N-j_lat+1,0)) then
       if (abs(lat(l_0N) - lat(l_0N-j_lat)) .GT. 1.) then
         if (s_max_ind .eq. 0) then   
           lat_max = l_0N-j_lat
           s_max_ind = s_max_ind + 1
         else
           if (precip(0,l_0N-j_lat,0) .gt. precip(0,lat_max,0)) then
            lat_max = l_0N-j_lat
           end if
         end if
       end if
     end if   
    end do
  else
    do j_lat = 1,(c_NT-l_0S)
      if ( precip(0,l_0S+j_lat,0) .gt. precip(0,l_0S+j_lat-1,0)) then
        if (abs(lat(l_0S) - lat(l_0S+j_lat)) .GT. 1.) then    
          if (s_max_ind .eq. 0) then
            lat_max = l_0S+j_lat
            s_max_ind = s_max_ind + 1
          else
            if (precip(0,l_0S+j_lat,0) .gt. precip(0,lat_max,0)) then
              lat_max = l_0S+j_lat
            end if
          end if
        end if
      end if  
    end do
  end if

  if (s_max_ind .gt. 0) then
    pr_S(0,j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(lat_max)))
  else
    pr_S(0,j_f-1) = 0. 
  end if
;--------------------------------------------------------------------------------
; Calculate precip at the equator

  if (l_0 .gt. 0)
    pr_E(0,j_f-1) = precip(0,l_0,0)
  else
    pr_E(0,j_f-1) = (precip(0,l_0N,0)+precip(0,l_0S,0))/2.
  end if
;--------------------------------------------------------------------------------
; Calculate difference in precip between the ITCZ and the equator

  pr_D(0,j_f-1) = precip(0,ITCZ,0) - pr_E(0,j_f-1)
;--------------------------------------------------------------------------------
; Find difference in precip between ITCZ and other hemisphere
  if (lat(ITCZ) .GT. 0. ) then
    precip_o = max(precip(0,l_ST:l_0S,0))
  else
    precip_o = max(precip(0,l_0N:l_NT,0))
  end if

  ;if (s_max_ind .gt. 0) then
  ;  pr_A(0,j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
  ;else
  ;  pr_A(0,j_f-1) = pr_D(0,j_f-1) 
  ;end if 

  if (s_max_ind .gt. 0) then
    if (lat(ITCZ) .gt. 0) then
      pr_A(0,j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
    else
      pr_A(0,j_f-1) = -(precip(0,ITCZ,0) - precip(0,lat_max,0))
    end if
  else
    if (lat(ITCZ) .gt. 0) then
      pr_A(0,j_f-1) = pr_D(0,j_f-1) 
    else
      pr_A(0,j_f-1) = -pr_D(0,j_f-1)
    end if
  end if


;--------------------------------------------------------------------------------
; Find difference in latitude between subtropical minima

  N_minimum =l_0N+ minind(precip(0,l_0N:l_NT,0))
  S_minimum =l_ST+ minind(precip(0,l_ST:l_0S,0))
  pr_W(0,j_f-1) = abs(doubletofloat(lat(N_minimum)-lat(S_minimum)))

  if (lat(ITCZ) .GT. 0.) then
    pr_H(0,j_f-1) = precip(0,ITCZ,0) - precip(0,N_minimum,0)
    pr_L(0,j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(N_minimum))) 
  else
    pr_H(0,j_f-1) = precip(0,ITCZ,0) - precip(0,S_minimum,0) 
    pr_L(0,j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(S_minimum)))
  end if

  pr_R(0,j_f-1) = precip(0,N_minimum,0) - precip(0,S_minimum,0)

;-------------------------------------------------------------------------------------  
; Old Indicators:
; Tropical Mean 20S-20N:

 do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -index_lat .AND. lat(j_lat-1) .LE. index_lat) then 
    if(lat(j_lat-2) .LT. -index_lat) then
      DN=2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI))+ \
            precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)))/2.     \
           *doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.+index_lat)/180.*PI
    else
      if(lat(j_lat) .GT. index_lat) then
        DN=DN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                            doubletofloat(index_lat-(lat(j_lat-1)+ lat(j_lat-2))/2.)/180.*PI
       else
        DN=DN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                    doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2.)/180.*PI
      end if
    end if
   end if
 end do

 DN = DN/doubletofloat(2.*PI*(sin(index_lat/180.*PI)-sin(-index_lat/180.*PI)))
 pr_Mi(0,j_f-1) = DN

;------------------------------------------------------------------------------------------
; Asymmetry:

 do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -index_lat .AND. lat(j_lat-1) .LE. 0.) then 
    if(lat(j_lat-2) .LT. -index_lat) then
      DN_SH=2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI))+ \
            precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)))/2.     \
           *doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.+index_lat)/180.*PI
    else
      if(lat(j_lat) .GT. 0.) then
        DN_SH=DN_SH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                            doubletofloat((lat(j_lat-1)+ lat(j_lat-2))/2.)/180.*PI
       else
        DN_SH=DN_SH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                    doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2.)/180.*PI
     end if
    end if
   end if

   if(lat(j_lat-1) .GE. 0. .AND. lat(j_lat-1) .LE. index_lat) then
    if(lat(j_lat-2) .LT. 0.) then
      DN_NH=2*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                  precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \         
                 ((lat(j_lat-1)+ lat(j_lat))/2.)/180.*PI
    else
      if(lat(j_lat) .GT. index_lat) then
        DN_NH=DN_NH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                        precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2.   * \
                      (index_lat-(lat(j_lat-1)+ lat(j_lat)-2)/2.)/180.*PI
      else
        DN_NH=DN_NH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                                  precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                     ((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2. )/180.*PI
      end if
    end if 
   end if
  end do

  DN_SH=DN_SH/doubletofloat(2.*PI*(-sin(-index_lat/180.*PI)))
  DN_NH=DN_NH/doubletofloat(2.*PI*(sin(index_lat/180.*PI))) 

   pr_Ai(0,j_f-1) = doubletofloat((DN_NH-DN_SH)/DN)


;------------------------------------------------------------------------------------------------------------
; Relative Eq. Prec.:

 do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -eq_lat .AND. lat(j_lat-1) .LE. eq_lat) then 
    if(lat(j_lat-2) .LT. -eq_lat) then
      EN=2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI))+ \
            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2.     \
           *doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.+eq_lat)/180.*PI
    else
      if(lat(j_lat) .GT. eq_lat) then
        EN=EN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                            doubletofloat(eq_lat-(lat(j_lat-1)+ lat(j_lat-2))/2.)/180.*PI
       else
        EN=EN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                    doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2.)/180.*PI
      end if
    end if
   end if
 end do

 EN = EN/doubletofloat(2.*PI*(sin(eq_lat/180.*PI)-sin(-eq_lat/180.*PI)))
 pr_Ei(0,j_f-1) = doubletofloat(EN/DN) -1.
 ;pr_Ei(j_f-1) = doubletofloat(EN)



  delete(DN)
  delete(EN)
  delete(DN_SH)
  delete(DN_NH)

  delete(l_0S)
  delete(l_0N)
  delete(l_0)
  delete(l_ST)
  delete(l_NT)
  delete(c_ST)
  delete(c_NT)

  delete(file1)
  delete(precip)
  delete(lat)
  delete(lat_number)

  delete(ITCZ)
  if (s_max_ind .gt. 0) then
    delete(lat_max)
  end if
  delete(s_max_ind)
  delete(N_minimum)
  delete(S_minimum)
 
 end do



;-----------------------------------------------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------------------------------------------
; CMIP:


do j_f=1,m_number

 file1=addfile("../../../HISTORICAL/plots/data_1979_2005/merged_files/zonmean/"+file_names(j_f-1)+"_tm_ZM.nc","r")

  lat=file1->lat
  lat_number = dimsizes(lat)
  precip=86400.*file1->pr

  l_0 = 0

  do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -trop_lat .AND. lat(j_lat-2) .LT. -trop_lat) then 
    l_ST = j_lat-1
   end if

   if(lat(j_lat-1) .GE. -crit_lat .AND. lat(j_lat-2) .LT. -crit_lat) then 
    c_ST = j_lat-1
   end if

   if(lat(j_lat-1) .GE. -index_lat .AND. lat(j_lat-2) .LT. -index_lat) then 
    i_ST = j_lat-1
   end if

   if(lat(j_lat-1) .GE. -eq_lat .AND. lat(j_lat-2) .LT. -eq_lat) then 
    i_Seq = j_lat-1
   end if

   if(lat(j_lat-1) .GE. 0. .AND. lat(j_lat-2) .LT. 0.) then
     if (lat(j_lat-1) .EQ. 0.) then
       l_0 = j_lat-1
       l_0N = j_lat-1
       l_0S = j_lat-2 
     else
       l_0N = j_lat-1
       l_0S = j_lat-2 
     end if
   end if

   if(lat(j_lat) .GE. trop_lat .AND. lat(j_lat-1) .LT. trop_lat) then 
     l_NT = j_lat-1
   end if

   if(lat(j_lat) .GE. crit_lat .AND. lat(j_lat-1) .LT. crit_lat) then 
    c_NT = j_lat-1
   end if

   if(lat(j_lat) .GE. index_lat .AND. lat(j_lat-1) .LT. index_lat) then 
    i_NT = j_lat-1
   end if

   if(lat(j_lat) .GE. eq_lat .AND. lat(j_lat-1) .LT. eq_lat) then 
    i_Neq = j_lat-1
   end if
 
   end do

;--------------------------------------------------------------------------------
; Find absolute precip. maximum
  
  
  ITCZ=c_ST + maxind(precip(0,c_ST:c_NT,0))
  pr_I(1,j_f-1) = abs(doubletofloat(lat(ITCZ)))

;--------------------------------------------------------------------------------
; Find second ITCZ (if exists)

  s_max_ind = 0

    if (lat(ITCZ) .GT. 0. ) then
    do j_lat= 1,(l_0N-c_ST)
     if ( precip(0,l_0N-j_lat,0) .gt. precip(0,l_0N-j_lat+1,0)) then
       if (abs(lat(l_0N) - lat(l_0N-j_lat)) .GT. 1.) then
         if (s_max_ind .eq. 0) then   
           lat_max = l_0N-j_lat
           s_max_ind = s_max_ind + 1
         else
           if (precip(0,l_0N-j_lat,0) .gt. precip(0,lat_max,0)) then
            lat_max = l_0N-j_lat
           end if
         end if
       end if
     end if   
    end do
  else
    do j_lat = 1,(c_NT-l_0S)
      if ( precip(0,l_0S+j_lat,0) .gt. precip(0,l_0S+j_lat-1,0)) then
        if (abs(lat(l_0S) - lat(l_0S+j_lat)) .GT. 1.) then    
          if (s_max_ind .eq. 0) then
            lat_max = l_0S+j_lat
            s_max_ind = s_max_ind + 1
          else
            if (precip(0,l_0S+j_lat,0) .gt. precip(0,lat_max,0)) then
              lat_max = l_0S+j_lat
            end if
          end if
        end if
      end if  
    end do
  end if

  if (s_max_ind .gt. 0) then
    pr_S(1,j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(lat_max)))
  else
    pr_S(1,j_f-1) = 0. 
  end if
;--------------------------------------------------------------------------------
; Calculate precip at the equator

  if (l_0 .gt. 0)
    pr_E(1,j_f-1) = precip(0,l_0,0)
  else
    pr_E(1,j_f-1) = (precip(0,l_0N,0)+precip(0,l_0S,0))/2.
  end if
;--------------------------------------------------------------------------------
; Calculate difference in precip between the ITCZ and the equator

  pr_D(1,j_f-1) = precip(0,ITCZ,0) - pr_E(1,j_f-1)
;--------------------------------------------------------------------------------
; Find difference in precip between ITCZ and other hemisphere
  if (lat(ITCZ) .GT. 0. ) then
    precip_o = max(precip(0,l_ST:l_0S,0))
  else
    precip_o = max(precip(0,l_0N:l_NT,0))
  end if

  ;if (s_max_ind .gt. 0) then
  ;  pr_A(1,j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
  ;else
  ;  pr_A(1,j_f-1) = pr_D(1,j_f-1) 
  ;end if 

  if (s_max_ind .gt. 0) then
    if (lat(ITCZ) .gt. 0) then
      pr_A(1,j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
    else
      pr_A(1,j_f-1) = -(precip(0,ITCZ,0) - precip(0,lat_max,0))
    end if
  else
    if (lat(ITCZ) .gt. 0) then
      pr_A(1,j_f-1) = pr_D(1,j_f-1) 
    else
      pr_A(1,j_f-1) = -pr_D(1,j_f-1)
    end if
  end if


;--------------------------------------------------------------------------------
; Find difference in latitude between subtropical minima

  N_minimum =l_0N+ minind(precip(0,l_0N:l_NT,0))
  S_minimum =l_ST+ minind(precip(0,l_ST:l_0S,0))
  pr_W(1,j_f-1) = abs(doubletofloat(lat(N_minimum)-lat(S_minimum)))

  if (lat(ITCZ) .GT. 0.) then
    pr_H(1,j_f-1) = precip(0,ITCZ,0) - precip(0,N_minimum,0)
    pr_L(1,j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(N_minimum))) 
  else
    pr_H(1,j_f-1) = precip(0,ITCZ,0) - precip(0,S_minimum,0) 
    pr_L(1,j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(S_minimum)))
  end if

  pr_R(1,j_f-1) = precip(0,N_minimum,0) - precip(0,S_minimum,0)

;-------------------------------------------------------------------------------------  
; Old Indicators:
; Tropical Mean 20S-20N:

 do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -index_lat .AND. lat(j_lat-1) .LE. index_lat) then 
    if(lat(j_lat-2) .LT. -index_lat) then
      DN=2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI))+ \
            precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)))/2.     \
           *doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.+index_lat)/180.*PI
    else
      if(lat(j_lat) .GT. index_lat) then
        DN=DN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                            doubletofloat(index_lat-(lat(j_lat-1)+ lat(j_lat-2))/2.)/180.*PI
       else
        DN=DN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                    doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2.)/180.*PI
      end if
    end if
   end if
 end do

 DN = DN/doubletofloat(2.*PI*(sin(index_lat/180.*PI)-sin(-index_lat/180.*PI)))
 pr_Mi(1,j_f-1) = DN

;------------------------------------------------------------------------------------------
; Asymmetry:

 do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -index_lat .AND. lat(j_lat-1) .LE. 0.) then 
    if(lat(j_lat-2) .LT. -index_lat) then
      DN_SH=2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI))+ \
            precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)))/2.     \
           *doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.+index_lat)/180.*PI
    else
      if(lat(j_lat) .GT. 0.) then
        DN_SH=DN_SH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                            doubletofloat((lat(j_lat-1)+ lat(j_lat-2))/2.)/180.*PI
       else
        DN_SH=DN_SH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                    doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2.)/180.*PI
     end if
    end if
   end if

   if(lat(j_lat-1) .GE. 0. .AND. lat(j_lat-1) .LE. index_lat) then
    if(lat(j_lat-2) .LT. 0.) then
      DN_NH=2*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                  precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \         
                 ((lat(j_lat-1)+ lat(j_lat))/2.)/180.*PI
    else
      if(lat(j_lat) .GT. index_lat) then
        DN_NH=DN_NH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                        precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2.   * \
                      (index_lat-(lat(j_lat-1)+ lat(j_lat)-2)/2.)/180.*PI
      else
        DN_NH=DN_NH+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                                  precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                     ((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2. )/180.*PI
      end if
    end if 
   end if
  end do

  DN_SH=DN_SH/doubletofloat(2.*PI*(-sin(-index_lat/180.*PI)))
  DN_NH=DN_NH/doubletofloat(2.*PI*(sin(index_lat/180.*PI))) 

   pr_Ai(1,j_f-1) = doubletofloat((DN_NH-DN_SH)/DN)


;------------------------------------------------------------------------------------------------------------
; Relative Eq. Prec.:

 do j_lat=2,lat_number-1
   if(lat(j_lat-1) .GE. -eq_lat .AND. lat(j_lat-1) .LE. eq_lat) then 
    if(lat(j_lat-2) .LT. -eq_lat) then
      EN=2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI))+ \
            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2.     \
           *doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.+eq_lat)/180.*PI
    else
      if(lat(j_lat) .GT. eq_lat) then
        EN=EN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                            doubletofloat(eq_lat-(lat(j_lat-1)+ lat(j_lat-2))/2.)/180.*PI
       else
        EN=EN+2.*PI*(precip(0,j_lat-1,0)*doubletofloat(cos(lat(j_lat-1)/180.*PI)) + \
                            precip(0,j_lat,0)*doubletofloat(cos(lat(j_lat)/180.*PI)))/2. * \
                    doubletofloat((lat(j_lat-1)+ lat(j_lat))/2.-(lat(j_lat-2)+ lat(j_lat-1))/2.)/180.*PI
      end if
    end if
   end if
 end do

 EN = EN/doubletofloat(2.*PI*(sin(eq_lat/180.*PI)-sin(-eq_lat/180.*PI)))
 pr_Ei(1,j_f-1) = doubletofloat(EN/DN) -1.
 ;pr_Ei(j_f-1) = doubletofloat(EN)



  delete(DN)
  delete(EN)
  delete(DN_SH)
  delete(DN_NH)

  delete(l_0S)
  delete(l_0N)
  delete(l_0)
  delete(l_ST)
  delete(l_NT)
  delete(c_ST)
  delete(c_NT)

  delete(file1)
  delete(precip)
  delete(lat)
  delete(lat_number)

  delete(ITCZ)
  if (s_max_ind .gt. 0) then
    delete(lat_max)
  end if
  delete(s_max_ind)
  delete(N_minimum)
  delete(S_minimum)
 
 end do



    cor_matrix_EE = covcorm(pr_E,(/1,0/))

    cor_matrix_DD = covcorm(pr_D,(/1,0/))

    cor_matrix_AA = covcorm(pr_A,(/1,0/))

    cor_matrix_HH = covcorm(pr_H,(/1,0/))

    cor_matrix_RR = covcorm(pr_R,(/1,0/))

    cor_matrix_II = covcorm(pr_I,(/1,0/))
 
    cor_matrix_SS = covcorm(pr_S,(/1,0/))

    cor_matrix_LL = covcorm(pr_L,(/1,0/))

    cor_matrix_WW = covcorm(pr_W,(/1,0/)) 


    cor_matrix_EiEi = covcorm(pr_Ei,(/1,0/))

    cor_matrix_MiMi = covcorm(pr_Mi,(/1,0/))

    cor_matrix_AiAi = covcorm(pr_Ai,(/1,0/))



    print("EE:")
    print(cor_matrix_EE(1)^2)

    print("DD:")
    print(cor_matrix_DD(1)^2)

    print("AA:")
    print(cor_matrix_AA(1)^2)

    print("HH:")
    print(cor_matrix_HH(1)^2)

    print("RR:")
    print(cor_matrix_RR(1)^2)

    print("II:")
    print(cor_matrix_II(1)^2)

    print("SS:")
    print(cor_matrix_SS(1)^2)

    print("LL:")
    print(cor_matrix_LL(1)^2)

    print("WW:")
    print(cor_matrix_WW(1)^2)


    print("AiAi:")
    print(cor_matrix_AiAi(1)^2)

    print("EiEi:")
    print(cor_matrix_EiEi(1)^2)

    print("MiMi:")
    print(cor_matrix_MiMi(1)^2)


end