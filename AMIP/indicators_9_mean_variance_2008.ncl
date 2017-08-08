load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin 

file_names=new (24,string)
;file_names=new (12,string)
file_names(0)="Amon_ACCESS1-0_amip_r1i1p1_197901-200812"
file_names(1)="Amon_ACCESS1-3_amip_r1i1p1_197901-200812"
file_names(2)="Amon_bcc-csm1-1_amip_r1i1p1_197901-200812"
file_names(3)="Amon_BNU-ESM_amip_r1i1p1_197901-200812"

file_names(4)="Amon_CanAM4_amip_r1i1p1_197901-200812"
file_names(5)="Amon_CCSM4_amip_r1i1p1_197901-200812"
file_names(6)="Amon_CNRM-CM5_amip_r1i1p1_197901-200812"
file_names(7)="Amon_FGOALS-s2_amip_r1i1p1_197901-200812"

file_names(8)="Amon_GFDL-CM3_form_amip_r1i1p1_197901-200812"
file_names(9)="Amon_GFDL-HIRAM-C180_form_amip_r1i1p1_197901-200812"
file_names(10)="Amon_GFDL-HIRAM-C360_form_amip_r1i1p1_197901-200812"
file_names(11)="Amon_GISS-E2-R_amip_r1i1p1_197901-200812"

file_names(12)="Amon_HadGEM2-A_amip_r1i1p1_197901-200812"
file_names(13)="Amon_inmcm4_amip_r1i1p1_197901-200812"
file_names(14)="Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200812"
file_names(15)="Amon_IPSL-CM5A-MR_amip_r1i1p1_197901-200812"

file_names(16)="Amon_IPSL-CM5B-LR_amip_r1i1p1_197901-200812"
file_names(17)="Amon_MIROC5_amip_r1i1p1_197901-200812"
file_names(18)="Amon_MPI-ESM-LR_amip_r1i1p1_197901-200812"
file_names(19)="Amon_MPI-ESM-MR_amip_r1i1p1_197901-200812"

file_names(20)="Amon_MRI-AGCM3-2H_amip_r1i1p1_197901-200812"
file_names(21)="Amon_MRI-AGCM3-2S_amip_r1i1p1_197901-200812"
file_names(22)="Amon_MRI-CGCM3_amip_r1i1p1_197901-200812"
file_names(23)="Amon_NorESM1-M_amip_r1i1p1_197901-200812"


m_number=dimsizes(file_names)
v_number=dimsizes(file_names)


  pr_W        = new(m_number,float)
  pr_E        = new(m_number,float)
  pr_A        = new(m_number,float)
  pr_S        = new(m_number,float)
  pr_I        = new(m_number,float)
  pr_D        = new(m_number,float)
  pr_H        = new(m_number,float)
  pr_L        = new(m_number,float)
  pr_R        = new(m_number,float)

  pr_W_s        = new(m_number,float)
  pr_E_s        = new(m_number,float)
  pr_A_s        = new(m_number,float)
  pr_S_s        = new(m_number,float)
  pr_I_s        = new(m_number,float)
  pr_D_s        = new(m_number,float)
  pr_H_s        = new(m_number,float)
  pr_L_s        = new(m_number,float)
  pr_R_s        = new(m_number,float)

  pr_W_se       = new(m_number,float)
  pr_E_se       = new(m_number,float)
  pr_A_se       = new(m_number,float)
  pr_S_se       = new(m_number,float)
  pr_I_se       = new(m_number,float)
  pr_D_se       = new(m_number,float)
  pr_H_se       = new(m_number,float)
  pr_L_se       = new(m_number,float)
  pr_R_se       = new(m_number,float)

  ind_vec     = new(m_number,integer)

 ind_vec = -1

 PI=acos(-1.)

trop_lat = 30.
crit_lat=20.

ct=0

do j_f=1,m_number

 file1=addfile("../data_1979_2008/merged_files/zonmean/"+file_names(j_f-1)+"_tm_ZM.nc","r")

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

   if(lat(j_lat-1) .GE. crit_lat .AND. lat(j_lat-2) .LT. crit_lat) then 
    c_NT = j_lat-1
   end if
 
   end do
;--------------------------------------------------------------------------------
; Find absolute precip. maximum
  
  
  ITCZ=c_ST + maxind(precip(0,c_ST:c_NT,0))
  pr_I(j_f-1) = abs(doubletofloat(lat(ITCZ)))

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
    pr_S(j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(lat_max)))
    ct=ct+1
    ind_vec(ct-1) = j_f-1
  else
    pr_S(j_f-1) = 0. 
  end if
;--------------------------------------------------------------------------------
; Calculate precip at the equator

  if (l_0 .gt. 0)
    pr_E(j_f-1) = precip(0,l_0,0)
  else
    pr_E(j_f-1) = (precip(0,l_0N,0)+precip(0,l_0S,0))/2.
  end if
;--------------------------------------------------------------------------------
; Calculate difference in precip between the ITCZ and the equator

  pr_D(j_f-1) = precip(0,ITCZ,0) - pr_E(j_f-1)
;--------------------------------------------------------------------------------
; Find difference in precip between ITCZ and other hemisphere
  if (lat(ITCZ) .GT. 0. ) then
    precip_o = max(precip(0,l_ST:l_0S,0))
  else
    precip_o = max(precip(0,l_0N:l_NT,0))
  end if

  ;if (s_max_ind .gt. 0) then
  ;  pr_A(j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
  ;else
  ;  pr_A(j_f-1) = pr_D(j_f-1) 
  ;end if 

  if (s_max_ind .gt. 0) then
    if (lat(ITCZ) .gt. 0) then
      pr_A(j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
    else
      pr_A(j_f-1) = -(precip(0,ITCZ,0) - precip(0,lat_max,0))
    end if
  else
    if (lat(ITCZ) .gt. 0) then
      pr_A(j_f-1) = pr_D(j_f-1) 
    else
      pr_A(j_f-1) = -pr_D(j_f-1)
    end if
  end if  

;--------------------------------------------------------------------------------
; Find difference in latitude between subtropical minima

  N_minimum =l_0N+ minind(precip(0,l_0N:l_NT,0))
  S_minimum =l_ST+ minind(precip(0,l_ST:l_0S,0))
  pr_W(j_f-1) = abs(doubletofloat(lat(N_minimum)-lat(S_minimum)))

  if (lat(ITCZ) .GT. 0.) then
    pr_H(j_f-1) = precip(0,ITCZ,0) - precip(0,N_minimum,0)
    pr_L(j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(N_minimum))) 
  else
    pr_H(j_f-1) = precip(0,ITCZ,0) - precip(0,S_minimum,0) 
    pr_L(j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(S_minimum)))
  end if

  pr_R(j_f-1) = precip(0,N_minimum,0) - precip(0,S_minimum,0)

;-------------------------------------------------------------------------------------  
  

  

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

;-------------------------------------------------------------------------------------------------------------------------
;GPCP

 file1=addfile("../data_1979_2008/merged_files/zonmean/pr_GPCP-SG_L3_v2.2_197901-200812_tm_ZM.nc","r")

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

   if(lat(j_lat-1) .GE. crit_lat .AND. lat(j_lat-2) .LT. crit_lat) then 
    c_NT = j_lat-1
   end if
 
   end do
;--------------------------------------------------------------------------------
; Find absolute precip. maximum
  
  
  ITCZ=c_ST + maxind(precip(0,c_ST:c_NT,0))
  pr_Ig = abs(doubletofloat(lat(ITCZ)))

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
    pr_Sg = abs(doubletofloat(lat(ITCZ) - lat(lat_max)))
  else
    pr_Sg = 0. 
  end if

;--------------------------------------------------------------------------------
; Calculate precip at the equator

  if (l_0 .gt. 0)
    pr_Eg = precip(0,l_0,0)
  else
    pr_Eg = (precip(0,l_0N,0)+precip(0,l_0S,0))/2.
  end if
;--------------------------------------------------------------------------------
; Calculate difference in precip between the ITCZ and the equator

  pr_Dg = precip(0,ITCZ,0) - pr_Eg
;--------------------------------------------------------------------------------
; Find difference in precip between ITCZ and other hemisphere
  if (lat(ITCZ) .GT. 0. ) then
    precip_o = max(precip(0,l_ST:l_0S,0))
  else
    precip_o = max(precip(0,l_0N:l_NT,0))
  end if

  ;if (s_max_ind .gt. 0) then
  ;  pr_Ag = precip(0,ITCZ,0) - precip(0,lat_max,0)
  ;else
  ;  pr_Ag = pr_Dg 
  ;end if 

  if (s_max_ind .gt. 0) then
    if (lat(ITCZ) .gt. 0) then
      pr_Ag = precip(0,ITCZ,0) - precip(0,lat_max,0)
    else
      pr_Ag = -(precip(0,ITCZ,0) - precip(0,lat_max,0))
    end if
  else
    if (lat(ITCZ) .gt. 0) then
      pr_Ag = pr_Dg 
    else
      pr_Ag = -pr_Dg
    end if
  end if

;--------------------------------------------------------------------------------
; Find difference in latitude between subtropical minima

  N_minimum =l_0N+ minind(precip(0,l_0N:l_NT,0))
  S_minimum =l_ST+ minind(precip(0,l_ST:l_0S,0))
  pr_Wg = abs(doubletofloat(lat(N_minimum)-lat(S_minimum)))

  if (lat(ITCZ) .GT. 0.) then
    pr_Hg = precip(0,ITCZ,0) - precip(0,N_minimum,0)
    pr_Lg = abs(doubletofloat(lat(ITCZ) - lat(N_minimum))) 
  else
    pr_Hg = precip(0,ITCZ,0) - precip(0,S_minimum,0) 
    pr_Lg = abs(doubletofloat(lat(ITCZ) - lat(S_minimum)))
  end if


  pr_Rg = precip(0,N_minimum,0) - precip(0,S_minimum,0)

;-------------------------------------------------------------------------------------  

do j_f =1,m_number  
 pr_W_s(j_f-1) = (pr_W(j_f-1)-pr_Wg)^2
 pr_E_s(j_f-1) = (pr_E(j_f-1)-pr_Eg)^2  
 pr_A_s(j_f-1) = (pr_A(j_f-1)-pr_Ag)^2
 pr_S_s(j_f-1) = (pr_S(j_f-1)-pr_Sg)^2
 pr_I_s(j_f-1) = (pr_I(j_f-1)-pr_Ig)^2  
 pr_D_s(j_f-1) = (pr_D(j_f-1)-pr_Dg)^2
 pr_H_s(j_f-1) = (pr_H(j_f-1)-pr_Hg)^2
 pr_L_s(j_f-1) = (pr_L(j_f-1)-pr_Lg)^2  
 pr_R_s(j_f-1) = (pr_R(j_f-1)-pr_Rg)^2
end do

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
 

;---------------------------------------------------------------------------------------------------------------------------

;ERA-interim

 file1=addfile("../data_1979_2008/merged_files/zonmean/pr_Amon_reanalysis_IFS-Cy31r2_197901-200812_tm_ZM.nc","r")

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

   if(lat(j_lat-1) .GE. crit_lat .AND. lat(j_lat-2) .LT. crit_lat) then 
    c_NT = j_lat-1
   end if
 
   end do
;--------------------------------------------------------------------------------
; Find absolute precip. maximum
  
  
  ITCZ=c_ST + maxind(precip(0,c_ST:c_NT,0))
  pr_Ie = abs(doubletofloat(lat(ITCZ)))

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
    pr_Se = abs(doubletofloat(lat(ITCZ) - lat(lat_max)))
  else
    pr_Se = 0. 
  end if

;--------------------------------------------------------------------------------
; Calculate precip at the equator

  if (l_0 .gt. 0)
    pr_Ee = precip(0,l_0,0)
  else
    pr_Ee = (precip(0,l_0N,0)+precip(0,l_0S,0))/2.
  end if
;--------------------------------------------------------------------------------
; Calculate difference in precip between the ITCZ and the equator

  pr_De = precip(0,ITCZ,0) - pr_Ee
;--------------------------------------------------------------------------------
; Find difference in precip between ITCZ and other hemisphere
  if (lat(ITCZ) .GT. 0. ) then
    precip_o = max(precip(0,l_ST:l_0S,0))
  else
    precip_o = max(precip(0,l_0N:l_NT,0))
  end if

  ;if (s_max_ind .gt. 0) then
  ;  pr_Ae = precip(0,ITCZ,0) - precip(0,lat_max,0)
  ;else
  ;  pr_Ae = pr_De
  ;end if

  if (s_max_ind .gt. 0) then
    if (lat(ITCZ) .gt. 0) then
      pr_Ae = precip(0,ITCZ,0) - precip(0,lat_max,0)
    else
      pr_Ae = -(precip(0,ITCZ,0) - precip(0,lat_max,0))
    end if
  else
    if (lat(ITCZ) .gt. 0) then
      pr_Ae = pr_De 
    else
      pr_Ae = -pr_De
    end if
  end if
 

;--------------------------------------------------------------------------------
; Find difference in latitude between subtropical minima

  N_minimum =l_0N+ minind(precip(0,l_0N:l_NT,0))
  S_minimum =l_ST+ minind(precip(0,l_ST:l_0S,0))
  pr_We = abs(doubletofloat(lat(N_minimum)-lat(S_minimum)))

  if (lat(ITCZ) .GT. 0.) then
    pr_He = precip(0,ITCZ,0) - precip(0,N_minimum,0)
    pr_Le = abs(doubletofloat(lat(ITCZ) - lat(N_minimum))) 
  else
    pr_He = precip(0,ITCZ,0) - precip(0,S_minimum,0) 
    pr_Le = abs(doubletofloat(lat(ITCZ) - lat(S_minimum)))
  end if

  pr_Re = precip(0,N_minimum,0) - precip(0,S_minimum,0)

;-------------------------------------------------------------------------------------  

do j_f =1,m_number  
 pr_W_se(j_f-1) = (pr_W(j_f-1)-pr_We)^2
 pr_E_se(j_f-1) = (pr_E(j_f-1)-pr_Ee)^2  
 pr_A_se(j_f-1) = (pr_A(j_f-1)-pr_Ae)^2
 pr_S_se(j_f-1) = (pr_S(j_f-1)-pr_Se)^2
 pr_I_se(j_f-1) = (pr_I(j_f-1)-pr_Ie)^2  
 pr_D_se(j_f-1) = (pr_D(j_f-1)-pr_De)^2
 pr_H_se(j_f-1) = (pr_H(j_f-1)-pr_He)^2
 pr_L_se(j_f-1) = (pr_L(j_f-1)-pr_Le)^2  
 pr_R_se(j_f-1) = (pr_R(j_f-1)-pr_Re)^2
end do

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

 ;-------------------------------------------------------------------------------------------------






pr_Sn = new(ct,float)

vt = 0

print(pr_S(:))

do j_f=1,m_number
 if (ind_vec(j_f-1) .ne. -1) then
   vt = vt + 1
   pr_Sn(vt-1) = pr_S(ind_vec(j_f-1))
 end if
end do


 print("W gpcp:")
 print(pr_Wg-0.)

 print("E gpcp:")
 print(pr_Eg-0.)

 print("A gpcp:")
 print(pr_Ag-0.)
 
 print("S gpcp:")
 print(pr_Sg-0.)

 print("I gpcp:")
 print(pr_Ig-0.)

 print("D gpcp:")
 print(pr_Dg-0.)

 print("H gpcp:")
 print(pr_Hg-0.)

 print("L gpcp:")
 print(pr_Lg-0.)

 print("R gpcp:")
 print(pr_Rg-0.)


 print("W ERA:")
 print(pr_We-0.)

 print("E ERA:")
 print(pr_Ee-0.)

 print("A ERA:")
 print(pr_Ae-0.)
 
 print("S ERA:")
 print(pr_Se-0.)

 print("I ERA:")
 print(pr_Ie-0.)

 print("D ERA:")
 print(pr_De-0.)

 print("H ERA:")
 print(pr_He-0.)

 print("L ERA:")
 print(pr_Le-0.)

 print("R ERA:")
 print(pr_Re-0.)




; print("Sn mean:")
; print(sum(pr_Sn(:))/ct)

 print("W mean:")
 print(sum(pr_W(:))/m_number)

 print("E mean:")
 print(sum(pr_E(:))/m_number)

 print("A mean:")
 print(sum(pr_A(:))/m_number)
 
 print("S mean:")
 print(sum(pr_S(:))/m_number)

 print("I mean:")
 print(sum(pr_I(:))/m_number)

 print("D mean:")
 print(sum(pr_D(:))/m_number)

 print("H mean:")
 print(sum(pr_H(:))/m_number)

 print("L mean:")
 print(sum(pr_L(:))/m_number)

 print("R mean:")
 print(sum(pr_R(:))/m_number)


 print("W dev gpcp:")
 print((sum(pr_W_s(:))/m_number)^0.5)

 print("E dev gpcp:")
 print((sum(pr_E_s(:))/m_number)^0.5)

 print("A dev gpcp:")
 print((sum(pr_A_s(:))/m_number)^0.5)
 
 print("S dev gpcp:")
 print((sum(pr_S_s(:))/m_number)^0.5)

 print("I dev gpcp:")
 print((sum(pr_I_s(:))/m_number)^0.5)

 print("D dev gpcp:")
 print((sum(pr_D_s(:))/m_number)^0.5)

 print("H dev gpcp:")
 print((sum(pr_H_s(:))/m_number)^0.5)

 print("L dev gpcp:")
 print((sum(pr_L_s(:))/m_number)^0.5)

 print("R dev gpcp:")
 print((sum(pr_R_s(:))/m_number)^0.5)



 print("W rel dev gpcp:")
 print((sum(pr_W_s(:))/m_number)^0.5/pr_Wg)

 print("E rel dev gpcp:")
 print((sum(pr_E_s(:))/m_number)^0.5/pr_Eg)

 print("A rel dev gpcp:")
 print((sum(pr_A_s(:))/m_number)^0.5/pr_Ag)
 
 print("S rel dev gpcp:")
 print((sum(pr_S_s(:))/m_number)^0.5/pr_Sg)

 print("I rel dev gpcp:")
 print((sum(pr_I_s(:))/m_number)^0.5/pr_Ig)

 print("D rel dev gpcp:")
 print((sum(pr_D_s(:))/m_number)^0.5/pr_Dg)

 print("H rel dev gpcp:")
 print((sum(pr_H_s(:))/m_number)^0.5/pr_Hg)

 print("L rel dev gpcp:")
 print((sum(pr_L_s(:))/m_number)^0.5/pr_Lg)

 print("R rel dev gpcp:")
 print((sum(pr_R_s(:))/m_number)^0.5/pr_Rg)



 print("W dev era:")
 print((sum(pr_W_se(:))/m_number)^0.5)

 print("E dev era:")
 print((sum(pr_E_se(:))/m_number)^0.5)

 print("A dev era:")
 print((sum(pr_A_se(:))/m_number)^0.5)
 
 print("S dev era:")
 print((sum(pr_S_se(:))/m_number)^0.5)

 print("I dev era:")
 print((sum(pr_I_se(:))/m_number)^0.5)

 print("D dev era:")
 print((sum(pr_D_se(:))/m_number)^0.5)

 print("H dev era:")
 print((sum(pr_H_se(:))/m_number)^0.5)

 print("L dev era:")
 print((sum(pr_L_se(:))/m_number)^0.5)

 print("R dev era:")
 print((sum(pr_R_se(:))/m_number)^0.5)



 print("W rel dev era:")
 print((sum(pr_W_se(:))/m_number)^0.5/pr_We)

 print("E rel dev era:")
 print((sum(pr_E_se(:))/m_number)^0.5/pr_Ee)

 print("A rel dev era:")
 print((sum(pr_A_se(:))/m_number)^0.5/pr_Ae)
 
 print("S rel dev era:")
 print((sum(pr_S_se(:))/m_number)^0.5/pr_Se)

 print("I rel dev era:")
 print((sum(pr_I_se(:))/m_number)^0.5/pr_Ie)

 print("D rel dev era:")
 print((sum(pr_D_se(:))/m_number)^0.5/pr_De)

 print("H rel dev era:")
 print((sum(pr_H_se(:))/m_number)^0.5/pr_He)

 print("L rel dev era:")
 print((sum(pr_L_se(:))/m_number)^0.5/pr_Le)

 print("R rel dev era:")
 print((sum(pr_R_se(:))/m_number)^0.5/pr_Re)



   value_matrix99 = new((/2,ct/),float)
   value_matrix99(0,:) = pr_Sn(:)
   value_matrix99(1,:) = pr_Sn(:)


   value_matrix00 = new((/2,m_number/),float)
   value_matrix00(0,:) = pr_W(:)
   value_matrix00(1,:) = pr_W(:)

   value_matrix11 = new((/2,m_number/),float)
   value_matrix11(0,:) = pr_E(:)
   value_matrix11(1,:) = pr_E(:)

   value_matrix22 = new((/2,m_number/),float)
   value_matrix22(0,:) = pr_A(:)
   value_matrix22(1,:) = pr_A(:)   

   value_matrix33 = new((/2,m_number/),float)
   value_matrix33(0,:) = pr_S(:)
   value_matrix33(1,:) = pr_S(:)

   value_matrix44 = new((/2,m_number/),float)
   value_matrix44(0,:) = pr_I(:)
   value_matrix44(1,:) = pr_I(:)

   value_matrix55 = new((/2,m_number/),float)
   value_matrix55(0,:) = pr_D(:)
   value_matrix55(1,:) = pr_D(:)

   value_matrix66 = new((/2,m_number/),float)
   value_matrix66(0,:) = pr_H(:)
   value_matrix66(1,:) = pr_H(:)

   value_matrix77 = new((/2,m_number/),float)
   value_matrix77(0,:) = pr_L(:)
   value_matrix77(1,:) = pr_L(:)

   value_matrix88 = new((/2,m_number/),float)
   value_matrix88(0,:) = pr_R(:)
   value_matrix88(1,:) = pr_R(:)



    cor_matrix_SnSn = covcorm(value_matrix99,(/0,1/))
  
    cor_matrix_WW = covcorm(value_matrix00,(/0,1/))
    cor_matrix_EE = covcorm(value_matrix11,(/0,1/))
    cor_matrix_AA = covcorm(value_matrix22,(/0,1/))
    cor_matrix_SS = covcorm(value_matrix33,(/0,1/))
    cor_matrix_II = covcorm(value_matrix44,(/0,1/))
    cor_matrix_DD = covcorm(value_matrix55,(/0,1/))
    cor_matrix_HH = covcorm(value_matrix66,(/0,1/))
    cor_matrix_LL = covcorm(value_matrix77,(/0,1/))
    cor_matrix_RR = covcorm(value_matrix88,(/0,1/))
 
 ;print(pr_Sn)


;    print("SnSn:")
;    print(cor_matrix_SnSn(0,0)^0.5)
  

    print("WW:")
    print(cor_matrix_WW(0,0)^0.5)

    print("EE:")
    print(cor_matrix_EE(0,0)^0.5)

    print("AA:")
    print(cor_matrix_AA(0,0)^0.5)

    print("SS:")
    print(cor_matrix_SS(0,0)^0.5)

    print("II:")
    print(cor_matrix_II(0,0)^0.5)

    print("DD:")
    print(cor_matrix_DD(0,0)^0.5)  

    print("HH:")
    print(cor_matrix_HH(0,0)^0.5)

    print("LL:")
    print(cor_matrix_LL(0,0)^0.5)

    print("RR:")
    print(cor_matrix_RR(0,0)^0.5)  



;    print("SnSn rel:")
;    print(cor_matrix_SnSn(0,0)^0.5/sum(pr_Sn(:))*ct)


    print("WW rel:")
    print(cor_matrix_WW(0,0)^0.5/sum(pr_W(:))*m_number)

    print("EE rel:")
    print(cor_matrix_EE(0,0)^0.5/sum(pr_E(:))*m_number)

    print("AA rel:")
    print(cor_matrix_AA(0,0)^0.5/sum(pr_A(:))*m_number)

    print("SS rel:")
    print(cor_matrix_SS(0,0)^0.5/sum(pr_S(:))*m_number)

    print("II rel:")
    print(cor_matrix_II(0,0)^0.5/sum(pr_I(:))*m_number)

    print("DD rel:")
    print(cor_matrix_DD(0,0)^0.5/sum(pr_D(:))*m_number)  

    print("HH rel:")
    print(cor_matrix_HH(0,0)^0.5/sum(pr_H(:))*m_number)

    print("LL rel:")
    print(cor_matrix_LL(0,0)^0.5/sum(pr_L(:))*m_number)

    print("RR rel:")
    print(cor_matrix_RR(0,0)^0.5/sum(pr_R(:))*m_number) 


   

end