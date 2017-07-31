load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin 

file_names=new (8,string)

file_names(0)="Amon_CCSM4_aquaControl_r1i1p1_197901-198312"
file_names(1)="Amon_CNRM-CM5_aquaControl_r1i1p1_197901-198312"
file_names(2)="Amon_FGOALS-s2_aquaControl_r1i1p1_197901-198312"
file_names(3)="Amon_IPSL-CM5A-LR_aquaControl_r1i1p1_197901-198312"
file_names(4)="Amon_MIROC5_aquaControl_r1i1p1_197901-198312"
file_names(5)="Amon_MPI-ESM-LR_aquaControl_r1i1p1_197901-198312"
file_names(6)="Amon_MPI-ESM-MR_aquaControl_r1i1p1_197901-198312"
file_names(7)="Amon_MRI-CGCM3_aquaControl_r1i1p1_197901-198312"


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


  pr_Ai        = new(m_number,float)
  pr_Ai_NH     = new(m_number,float)
  pr_Ai_SH     = new(m_number,float)
  pr_Ei        = new(m_number,float)
  pr_Mi        = new(m_number,float)

  pr_Ae        = new(m_number,float)
  pr_Ae_NH     = new(m_number,float)
  pr_Ae_SH     = new(m_number,float)
  pr_Ee        = new(m_number,float)
  pr_Me        = new(m_number,float)

 PI=acos(-1.)

trop_lat = 30.
crit_lat=20.

index_lat=20.
eq_lat = 2.

do j_f=1,m_number

 file1=addfile("../data_1979_1983/zonmean/"+file_names(j_f-1)+"_tm_ZM.nc","r")

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
  pr_I(j_f-1) = abs(doubletofloat(lat(ITCZ)))

;--------------------------------------------------------------------------------
; Find second ITCZ (if exists)

  N_maximum =l_0N+ maxind(precip(0,l_0N:l_NT,0))
  S_maximum =l_ST+ maxind(precip(0,l_ST:l_0S,0))

  N_minimum =l_0N+ minind(precip(0,l_0N:l_NT,0))
  S_minimum =l_ST+ minind(precip(0,l_ST:l_0S,0))

  s_max_ind = 0


   if (lat(ITCZ) .GT. 0. ) then
    do j_lat= 1,(l_0N-S_minimum)
     if ( precip(0,l_0N-j_lat,0) .gt. precip(0,l_0N-j_lat+1,0)) then
       if (abs(lat(l_0N) - lat(l_0N-j_lat)) .GT. 1.) then
         if (s_max_ind .eq. 0) then  
           print(precip(0,l_0N-j_lat,0)-0.)
           print(precip(0,l_0N-j_lat+1,0)-0.) 
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
    do j_lat = 1,(N_minimum-l_0S)
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

  ;pr_A(j_f-1) = precip(0,N_maximum,0) - precip(0,S_maximum,0)

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
 pr_Mi(j_f-1) = DN

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

  pr_Ai_NH(j_f-1) = doubletofloat(DN_NH)
  pr_Ai_SH(j_f-1) = doubletofloat(DN_SH)

   pr_Ai(j_f-1) = doubletofloat((DN_NH-DN_SH)/DN)


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
 pr_Ei(j_f-1) = doubletofloat(EN/DN) -1.
 ;pr_Ei(j_f-1) = doubletofloat(EN)



;-------------------------------------------------------------------------------------  
;calculate estimated tropical mean:

  pr_Me(j_f-1) = doubletofloat( (cos(lat(S_minimum)*PI/180.)*precip(0,S_minimum,0) + cos(lat(S_maximum)*PI/180.)*precip(0,S_maximum,0)) * (lat(S_maximum)- lat(S_minimum))*PI/180. - \  
               (cos(lat(S_maximum)*PI/180.)*precip(0,S_maximum,0) + pr_E(j_f-1)) * lat(S_maximum)*PI/180. + \
               (cos(lat(N_maximum)*PI/180.)*precip(0,N_maximum,0) + pr_E(j_f-1)) * lat(N_maximum)*PI/180. +  \
               (cos(lat(N_minimum)*PI/180.)*precip(0,N_minimum,0) + cos(lat(N_maximum)*PI/180.)*precip(0,N_maximum,0)) * (lat(N_minimum)- lat(N_maximum))*PI/180.)

  pr_Me(j_f-1) = pr_Me(j_f-1) / doubletofloat(2.*(sin(lat(N_minimum)*PI/180.) - sin(lat(S_minimum)*PI/180.)))


;calculate estimated equatorial precipitation index:

  pr_Ee(j_f-1) = pr_E(j_f-1) / pr_Me(j_f-1)


;calculate estimated asymm. precipitation index:

  ;pr_Ae(j_f-1) = doubletofloat( (cos(lat(N_maximum)*PI/180.)*precip(0,N_maximum,0) + pr_E(j_f-1)) * lat(N_maximum)*PI/180. +  \
  ;             (cos(lat(N_minimum)*PI/180.)*precip(0,N_minimum,0) + cos(lat(N_maximum)*PI/180.)*precip(0,N_maximum,0)) * (lat(N_minimum)- lat(N_maximum))*PI/180.) / \
  ;             doubletofloat(2.*sin(lat(N_minimum)*PI/180.)) - \
  ;             doubletofloat( (cos(lat(S_minimum)*PI/180.)*precip(0,S_minimum,0) + cos(lat(S_maximum)*PI/180.)*precip(0,S_maximum,0)) * (lat(S_maximum)- lat(S_minimum))*PI/180. - \  
  ;             (cos(lat(S_maximum)*PI/180.)*precip(0,S_maximum,0) + pr_E(j_f-1)) * lat(S_maximum)*PI/180. ) / \  
  ;             doubletofloat(-2.*sin(lat(S_minimum)*PI/180.))
  ;
  ; pr_Ae(j_f-1) = pr_Ae(j_f-1) / pr_Me(j_f-1)

  pr_Ae_NH(j_f-1) = doubletofloat((cos(lat(N_maximum)*PI/180.)*precip(0,N_maximum,0) + pr_E(j_f-1)) * lat(N_maximum)*PI/180. +  \
               (cos(lat(N_minimum)*PI/180.)*precip(0,N_minimum,0) + cos(lat(N_maximum)*PI/180.)*precip(0,N_maximum,0)) * (lat(N_minimum)- lat(N_maximum))*PI/180.)/2.
  pr_Ae_NH(j_f-1) = pr_Ae_NH(j_f-1)/doubletofloat(sin(lat(N_minimum)*PI/180.))

  pr_Ae_SH(j_f-1) = doubletofloat( (cos(lat(S_minimum)*PI/180.)*precip(0,S_minimum,0) + cos(lat(S_maximum)*PI/180.)*precip(0,S_maximum,0)) * (lat(S_maximum)- lat(S_minimum))*PI/180. - \  
               (cos(lat(S_maximum)*PI/180.)*precip(0,S_maximum,0) + pr_E(j_f-1)) * lat(S_maximum)*PI/180.) / 2. 

  pr_Ae_SH(j_f-1) = pr_Ae_SH(j_f-1) / doubletofloat(-sin(lat(S_minimum)*PI/180.))

  pr_Ae(j_f-1) = (pr_Ae_NH(j_f-1) - pr_Ae_SH(j_f-1)) / pr_Me(j_f-1) 



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


  print(pr_Ei(:))
  print(pr_Ai(:))
  print(pr_Mi(:))
  

   value_matrixWEi = new((/2,m_number/),float)
   value_matrixWEi(0,:) = pr_W(:)
   value_matrixWEi(1,:) = pr_Ei(:)

   value_matrixWAi = new((/2,m_number/),float)
   value_matrixWAi(0,:) = pr_W(:)
   value_matrixWAi(1,:) = pr_Ai(:)

   value_matrixWMi = new((/2,m_number/),float)
   value_matrixWMi(0,:) = pr_W(:)
   value_matrixWMi(1,:) = pr_Mi(:)   

 
   value_matrixEEi = new((/2,m_number/),float)
   value_matrixEEi(0,:) = pr_E(:)
   value_matrixEEi(1,:) = pr_Ei(:)

   value_matrixEAi = new((/2,m_number/),float)
   value_matrixEAi(0,:) = pr_E(:)
   value_matrixEAi(1,:) = pr_Ai(:)

   value_matrixEMi = new((/2,m_number/),float)
   value_matrixEMi(0,:) = pr_E(:)
   value_matrixEMi(1,:) = pr_Mi(:)  


   value_matrixAEi = new((/2,m_number/),float)
   value_matrixAEi(0,:) = pr_A(:)
   value_matrixAEi(1,:) = pr_Ei(:)

   value_matrixAAi = new((/2,m_number/),float)
   value_matrixAAi(0,:) = pr_A(:)
   value_matrixAAi(1,:) = pr_Ai(:)

   value_matrixAMi = new((/2,m_number/),float)
   value_matrixAMi(0,:) = pr_A(:)
   value_matrixAMi(1,:) = pr_Mi(:)  


   value_matrixSEi = new((/2,m_number/),float)
   value_matrixSEi(0,:) = pr_S(:)
   value_matrixSEi(1,:) = pr_Ei(:)

   value_matrixSAi = new((/2,m_number/),float)
   value_matrixSAi(0,:) = pr_S(:)
   value_matrixSAi(1,:) = pr_Ai(:)

   value_matrixSMi = new((/2,m_number/),float)
   value_matrixSMi(0,:) = pr_S(:)
   value_matrixSMi(1,:) = pr_Mi(:)  


   value_matrixIEi = new((/2,m_number/),float)
   value_matrixIEi(0,:) = pr_I(:)
   value_matrixIEi(1,:) = pr_Ei(:)

   value_matrixIAi = new((/2,m_number/),float)
   value_matrixIAi(0,:) = pr_I(:)
   value_matrixIAi(1,:) = pr_Ai(:)

   value_matrixIMi = new((/2,m_number/),float)
   value_matrixIMi(0,:) = pr_I(:)
   value_matrixIMi(1,:) = pr_Mi(:)  


   value_matrixDEi = new((/2,m_number/),float)
   value_matrixDEi(0,:) = pr_D(:)
   value_matrixDEi(1,:) = pr_Ei(:)

   value_matrixDAi = new((/2,m_number/),float)
   value_matrixDAi(0,:) = pr_D(:)
   value_matrixDAi(1,:) = pr_Ai(:)

   value_matrixDMi = new((/2,m_number/),float)
   value_matrixDMi(0,:) = pr_D(:)
   value_matrixDMi(1,:) = pr_Mi(:)  


   value_matrixHEi = new((/2,m_number/),float)
   value_matrixHEi(0,:) = pr_H(:)
   value_matrixHEi(1,:) = pr_Ei(:)

   value_matrixHAi = new((/2,m_number/),float)
   value_matrixHAi(0,:) = pr_H(:)
   value_matrixHAi(1,:) = pr_Ai(:)

   value_matrixHMi = new((/2,m_number/),float)
   value_matrixHMi(0,:) = pr_H(:)
   value_matrixHMi(1,:) = pr_Mi(:)  


   value_matrixLEi = new((/2,m_number/),float)
   value_matrixLEi(0,:) = pr_L(:)
   value_matrixLEi(1,:) = pr_Ei(:)

   value_matrixLAi = new((/2,m_number/),float)
   value_matrixLAi(0,:) = pr_L(:)
   value_matrixLAi(1,:) = pr_Ai(:)

   value_matrixLMi = new((/2,m_number/),float)
   value_matrixLMi(0,:) = pr_L(:)
   value_matrixLMi(1,:) = pr_Mi(:)


   value_matrixREi = new((/2,m_number/),float)
   value_matrixREi(0,:) = pr_R(:)
   value_matrixREi(1,:) = pr_Ei(:)

   value_matrixRAi = new((/2,m_number/),float)
   value_matrixRAi(0,:) = pr_R(:)
   value_matrixRAi(1,:) = pr_Ai(:)

   value_matrixRMi = new((/2,m_number/),float)
   value_matrixRMi(0,:) = pr_R(:)
   value_matrixRMi(1,:) = pr_Mi(:)




   value_matrixEeEi = new((/2,m_number/),float)
   value_matrixEeEi(0,:) = pr_Ee(:)
   value_matrixEeEi(1,:) = pr_Ei(:)

   value_matrixEeAi = new((/2,m_number/),float)
   value_matrixEeAi(0,:) = pr_Ee(:)
   value_matrixEeAi(1,:) = pr_Ai(:)

   value_matrixEeMi = new((/2,m_number/),float)
   value_matrixEeMi(0,:) = pr_Ee(:)
   value_matrixEeMi(1,:) = pr_Mi(:)


   value_matrixAeEi = new((/2,m_number/),float)
   value_matrixAeEi(0,:) = pr_Ae(:)
   value_matrixAeEi(1,:) = pr_Ei(:)

   value_matrixAeAi = new((/2,m_number/),float)
   value_matrixAeAi(0,:) = pr_Ae(:)
   value_matrixAeAi(1,:) = pr_Ai(:)

   value_matrixAeMi = new((/2,m_number/),float)
   value_matrixAeMi(0,:) = pr_Ae(:)
   value_matrixAeMi(1,:) = pr_Mi(:)


   value_matrixMeEi = new((/2,m_number/),float)
   value_matrixMeEi(0,:) = pr_Me(:)
   value_matrixMeEi(1,:) = pr_Ei(:)

   value_matrixMeAi = new((/2,m_number/),float)
   value_matrixMeAi(0,:) = pr_Me(:)
   value_matrixMeAi(1,:) = pr_Ai(:)

   value_matrixMeMi = new((/2,m_number/),float)
   value_matrixMeMi(0,:) = pr_Me(:)
   value_matrixMeMi(1,:) = pr_Mi(:)
  

 
    cor_matrix_WEi = covcorm(value_matrixWEi,(/1,0/))
    cor_matrix_WAi = covcorm(value_matrixWAi,(/1,0/))
    cor_matrix_WMi = covcorm(value_matrixWMi,(/1,0/))

    cor_matrix_EEi = covcorm(value_matrixEEi,(/1,0/))
    cor_matrix_EAi = covcorm(value_matrixEAi,(/1,0/))
    cor_matrix_EMi = covcorm(value_matrixEMi,(/1,0/))

    cor_matrix_AEi = covcorm(value_matrixAEi,(/1,0/))
    cor_matrix_AAi = covcorm(value_matrixAAi,(/1,0/))
    cor_matrix_AMi = covcorm(value_matrixAMi,(/1,0/))

    cor_matrix_SEi = covcorm(value_matrixSEi,(/1,0/))
    cor_matrix_SAi = covcorm(value_matrixSAi,(/1,0/))
    cor_matrix_SMi = covcorm(value_matrixSMi,(/1,0/))

    cor_matrix_IEi = covcorm(value_matrixIEi,(/1,0/))
    cor_matrix_IAi = covcorm(value_matrixIAi,(/1,0/))
    cor_matrix_IMi = covcorm(value_matrixIMi,(/1,0/))

    cor_matrix_DEi = covcorm(value_matrixDEi,(/1,0/))
    cor_matrix_DAi = covcorm(value_matrixDAi,(/1,0/))
    cor_matrix_DMi = covcorm(value_matrixDMi,(/1,0/))

    cor_matrix_HEi = covcorm(value_matrixHEi,(/1,0/))
    cor_matrix_HAi = covcorm(value_matrixHAi,(/1,0/))
    cor_matrix_HMi = covcorm(value_matrixHMi,(/1,0/))

    cor_matrix_LEi = covcorm(value_matrixLEi,(/1,0/))
    cor_matrix_LAi = covcorm(value_matrixLAi,(/1,0/))
    cor_matrix_LMi = covcorm(value_matrixLMi,(/1,0/))

    cor_matrix_REi = covcorm(value_matrixREi,(/1,0/))
    cor_matrix_RAi = covcorm(value_matrixRAi,(/1,0/))
    cor_matrix_RMi = covcorm(value_matrixRMi,(/1,0/))


    cor_matrix_EeEi = covcorm(value_matrixEeEi,(/1,0/))
    cor_matrix_EeAi = covcorm(value_matrixEeAi,(/1,0/))
    cor_matrix_EeMi = covcorm(value_matrixEeMi,(/1,0/))

    cor_matrix_AeEi = covcorm(value_matrixAeEi,(/1,0/))
    cor_matrix_AeAi = covcorm(value_matrixAeAi,(/1,0/))
    cor_matrix_AeMi = covcorm(value_matrixAeMi,(/1,0/))

    cor_matrix_MeEi = covcorm(value_matrixMeEi,(/1,0/))
    cor_matrix_MeAi = covcorm(value_matrixMeAi,(/1,0/))
    cor_matrix_MeMi = covcorm(value_matrixMeMi,(/1,0/))


   
    print("WEi:")
    print(cor_matrix_WEi(1)^2)

    print("WAi:")
    print(cor_matrix_WAi(1)^2)

    print("WMi:")
    print(cor_matrix_WMi(1)^2)


    print("EEi:")
    print(cor_matrix_EEi(1)^2)

    print("EAi:")
    print(cor_matrix_EAi(1)^2)

    print("EMi:")
    print(cor_matrix_EMi(1)^2)


    print("AEi:")
    print(cor_matrix_AEi(1)^2)

    print("AAi:")
    print(cor_matrix_AAi(1)^2)

    print("AMi:")
    print(cor_matrix_AMi(1)^2)


    print("SEi:")
    print(cor_matrix_SEi(1)^2)

    print("SAi:")
    print(cor_matrix_SAi(1)^2)

    print("SMi:")
    print(cor_matrix_SMi(1)^2)


    print("IEi:")
    print(cor_matrix_IEi(1)^2)

    print("IAi:")
    print(cor_matrix_WAi(1)^2)

    print("IMi:")
    print(cor_matrix_IMi(1)^2)


    print("DEi:")
    print(cor_matrix_DEi(1)^2)

    print("DAi:")
    print(cor_matrix_DAi(1)^2)

    print("DMi:")
    print(cor_matrix_DMi(1)^2)


    print("HEi:")
    print(cor_matrix_HEi(1)^2)

    print("HAi:")
    print(cor_matrix_HAi(1)^2)

    print("HMi:")
    print(cor_matrix_HMi(1)^2)


    print("LEi:")
    print(cor_matrix_LEi(1)^2)

    print("LAi:")
    print(cor_matrix_LAi(1)^2)

    print("LMi:")
    print(cor_matrix_LMi(1)^2)


    print("REi:")
    print(cor_matrix_REi(1)^2)

    print("RAi:")
    print(cor_matrix_RAi(1)^2)

    print("RMi:")
    print(cor_matrix_RMi(1)^2)



    print("EeEi:")
    print(cor_matrix_EeEi(1)^2)

    print("EeAi:")
    print(cor_matrix_EeAi(1)^2)

    print("EeMi:")
    print(cor_matrix_EeMi(1)^2)


    print("AeEi:")
    print(cor_matrix_AeEi(1)^2)

    print("AeAi:")
    print(cor_matrix_AeAi(1)^2)

    print("AeMi:")
    print(cor_matrix_AeMi(1)^2)


    print("MeEi:")
    print(cor_matrix_MeEi(1)^2)

    print("MeAi:")
    print(cor_matrix_MeAi(1)^2)

    print("MeMi:")
    print(cor_matrix_MeMi(1)^2)
    

    ;print(pr_Mi)
    ;print(pr_Me)

    ;print(pr_Ae_SH)
    ;print(pr_Ai_SH)
    ;print(pr_Ae_NH)
    ;print(pr_Ai_NH)

    ;print(pr_Ae)
    ;print(pr_Ai)
    

 
end