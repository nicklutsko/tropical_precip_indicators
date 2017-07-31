load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin 

file_names=new (26,string)

file_names(0)="Amon_ACCESS1-0_historical_r1i1p1_197901-200512"
file_names(1)="Amon_ACCESS1-3_historical_r1i1p1_197901-200512"
file_names(2)="Amon_bcc-csm1-1_historical_r1i1p1_197901-200512"
file_names(3)="Amon_BNU-ESM_historical_r1i1p1_197901-200512"

file_names(4)="Amon_CanESM2_historical_r1i1p1_197901-200512"
file_names(5)="Amon_CCSM4_historical_r1i1p1_197901-200512"
file_names(6)="Amon_CESM1-BGC_historical_r1i1p1_197901-200512"
file_names(7)="Amon_CNRM-CM5_historical_r1i1p1_197901-200512"

file_names(8)="Amon_CSIRO-Mk3-6-0_historical_r1i1p1_197901-200512"
file_names(9)="Amon_FGOALS-s2_historical_r1i1p1_197901-200512"
file_names(10)="Amon_GFDL-CM3_form_historical_r1i1p1_197901-200512"
file_names(11)="Amon_GFDL-ESM2G_form_historical_r1i1p1_197901-200512"

file_names(12)="Amon_GFDL-ESM2M_form_historical_r1i1p1_197901-200512"
file_names(13)="Amon_GISS-E2-R_historical_r1i1p1_197901-200512"
file_names(14)="Amon_GISS-E2-H_historical_r1i1p1_197901-200512"
file_names(15)="Amon_HadGEM2-ES_historical_r1i1p1_197901-200512"

file_names(16)="Amon_inmcm4_historical_r1i1p1_197901-200512"
file_names(17)="Amon_IPSL-CM5A-LR_historical_r1i1p1_197901-200512"
file_names(18)="Amon_IPSL-CM5A-MR_historical_r1i1p1_197901-200512"
file_names(19)="Amon_IPSL-CM5B-LR_historical_r1i1p1_197901-200512"

file_names(20)="Amon_MIROC5_historical_r1i1p1_197901-200512"
file_names(21)="Amon_MIROC-ESM_historical_r1i1p1_197901-200512"
file_names(22)="Amon_MPI-ESM-LR_historical_r1i1p1_197901-200512"
file_names(23)="Amon_MPI-ESM-MR_historical_r1i1p1_197901-200512"

file_names(24)="Amon_MRI-CGCM3_historical_r1i1p1_197901-200512"
file_names(25)="Amon_NorESM1-M_historical_r1i1p1_197901-200512"

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
  pr_Ai_pred   = new(m_number,float)
  pr_Ei        = new(m_number,float)
  pr_Mi        = new(m_number,float)

 PI=acos(-1.)

trop_lat = 30.
crit_lat=20.

index_lat=20.
eq_lat = 2.

do j_f=1,m_number

 file1=addfile("../data_1979_2005/merged_files/zonmean/"+file_names(j_f-1)+"_tm_ZM.nc","r")

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















    ;ind_numb = 10
    ;X_matrix = new((/ind_numb,m_number/),float)
    ;X_matrix(0,:) = 1. 
    ;X_matrix(1,:) = pr_E(:)
    ;X_matrix(2,:) = pr_D(:)
    ;X_matrix(3,:) = pr_A(:)
    ;X_matrix(4,:) = pr_H(:)
    ;X_matrix(5,:) = pr_R(:)
    ;X_matrix(6,:) = pr_I(:)
    ;X_matrix(7,:) = pr_S(:)
    ;X_matrix(8,:) = pr_L(:)
    ;X_matrix(9,:) = pr_W(:)


   ;ind_numb = 8
   ; X_matrix = new((/ind_numb,m_number/),float)
   ; X_matrix(0,:) = 1. 
   ; X_matrix(1,:) = pr_E(:)
   ; X_matrix(2,:) = pr_D(:)
   ; X_matrix(3,:) = pr_A(:)
   ; X_matrix(4,:) = pr_H(:)
   ; X_matrix(5,:) = pr_R(:)
   ; X_matrix(6,:) = pr_I(:)
   ; X_matrix(7,:) = pr_L(:)


   ;ind_numb = 7
   ;X_matrix = new((/ind_numb,m_number/),float)
   ;X_matrix(0,:) = 1. 
   ;X_matrix(1,:) = pr_E(:)
   ;X_matrix(2,:) = pr_A(:)
   ;X_matrix(3,:) = pr_H(:)
   ;X_matrix(4,:) = pr_R(:)
   ;X_matrix(5,:) = pr_I(:)
   ;X_matrix(6,:) = pr_L(:)


   ;ind_numb = 6
   ;X_matrix = new((/ind_numb,m_number/),float)
   ;X_matrix(0,:) = 1. 
   ;X_matrix(1,:) = pr_E(:)
   ;X_matrix(2,:) = pr_A(:)
   ;X_matrix(3,:) = pr_R(:)
   ;X_matrix(4,:) = pr_I(:)
   ;X_matrix(5,:) = pr_L(:)

   ;ind_numb = 5
   ;X_matrix = new((/ind_numb,m_number/),float)
   ;X_matrix(0,:) = 1. 
   ;X_matrix(1,:) = pr_E(:)
   ;X_matrix(2,:) = pr_A(:)
   ;X_matrix(3,:) = pr_I(:)
   ;X_matrix(4,:) = pr_L(:)

   ;ind_numb = 4
   ;X_matrix = new((/ind_numb,m_number/),float)
   ;X_matrix(0,:) = 1. 
   ;X_matrix(1,:) = pr_A(:)
   ;X_matrix(2,:) = pr_I(:)
   ;X_matrix(3,:) = pr_L(:)

   ind_numb = 3
   X_matrix = new((/ind_numb,m_number/),float)
   X_matrix(0,:) = 1. 
   X_matrix(1,:) = pr_A(:)
   X_matrix(2,:) = pr_I(:)

   
 
 
   
 

    ;cor_matrix = covcorm(X_matrix,(/1,0/))
  
  
    ;print("E:")
    ;print(cor_matrix(1)-0.)
  
    beta=reg_multlin(pr_Ai(:),X_matrix,False)
    print("beta:")
    print(beta)


    ;pr_Ai_pred(:)  = beta(0)   +beta(1)*pr_E(:)+ beta(2)*pr_D(:)+ beta(3)*pr_A(:)+ beta(4)*pr_M(:) + \
    ;                    beta(5)*pr_R(:)+ beta(6)*pr_I(:)+ beta(7)*pr_S(:)+ beta(8)*pr_L(:)+ beta(9)*pr_W(:) 

    pr_Ai_pred(:)  = beta(0)
    do j_i=2,ind_numb
     pr_Ai_pred(:) =  pr_Ai_pred(:) + beta(j_i-1)*X_matrix(j_i-1,:)
    end do
    
    ;delete(cor_matrix) 

    pred_matrix = new((/2,m_number/),float)
    pred_matrix(0,:) = pr_Ai_pred(:)
    pred_matrix(1,:) = pr_Ai(:)
 

    cor_matrix = covcorm(pred_matrix,(/1,0/))
  
  
    print("Asymmetry:")
    print((cor_matrix(1)-0.)^2)     

end