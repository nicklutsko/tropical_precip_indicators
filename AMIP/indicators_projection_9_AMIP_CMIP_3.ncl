load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin 

;--------------------------------------------------------------------------------------------------------------------------
;AMIP
;--------------------------------------------------------------------------------------------------------------------------

file_names=new (24,string)
;file_names=new (12,string)
file_names(0)="Amon_ACCESS1-0_amip_r1i1p1_197901-200512"
file_names(1)="Amon_ACCESS1-3_amip_r1i1p1_197901-200512"
file_names(2)="Amon_bcc-csm1-1_amip_r1i1p1_197901-200512"
file_names(3)="Amon_BNU-ESM_amip_r1i1p1_197901-200512"

file_names(4)="Amon_CanAM4_amip_r1i1p1_197901-200512"
file_names(5)="Amon_CCSM4_amip_r1i1p1_197901-200512"
file_names(6)="Amon_CNRM-CM5_amip_r1i1p1_197901-200512"
file_names(7)="Amon_FGOALS-s2_amip_r1i1p1_197901-200512"

file_names(8)="Amon_GFDL-CM3_form_amip_r1i1p1_197901-200512"
file_names(9)="Amon_GFDL-HIRAM-C180_form_amip_r1i1p1_197901-200512"
file_names(10)="Amon_GFDL-HIRAM-C360_form_amip_r1i1p1_197901-200512"
file_names(11)="Amon_GISS-E2-R_amip_r1i1p1_197901-200512"

file_names(12)="Amon_HadGEM2-A_amip_r1i1p1_197901-200512"
file_names(13)="Amon_inmcm4_amip_r1i1p1_197901-200512"
file_names(14)="Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200512"
file_names(15)="Amon_IPSL-CM5A-MR_amip_r1i1p1_197901-200512"

file_names(16)="Amon_IPSL-CM5B-LR_amip_r1i1p1_197901-200512"
file_names(17)="Amon_MIROC5_amip_r1i1p1_197901-200512"
file_names(18)="Amon_MPI-ESM-LR_amip_r1i1p1_197901-200512"
file_names(19)="Amon_MPI-ESM-MR_amip_r1i1p1_197901-200512"

file_names(20)="Amon_MRI-AGCM3-2H_amip_r1i1p1_197901-200512"
file_names(21)="Amon_MRI-AGCM3-2S_amip_r1i1p1_197901-200512"
file_names(22)="Amon_MRI-CGCM3_amip_r1i1p1_197901-200512"
file_names(23)="Amon_NorESM1-M_amip_r1i1p1_197901-200512"

print("AMIP names")
print(file_names)

m_number=dimsizes(file_names)
v_number=dimsizes(file_names)

precip_eq = new(m_number,float)

  pr_W        = new(m_number,float)
  pr_H        = new(m_number,float)
  pr_E        = new(m_number,float)
  pr_A        = new(m_number,float)
  pr_S        = new(m_number,float)
  pr_I        = new(m_number,float)
  pr_D        = new(m_number,float)
  pr_M        = new(m_number,float)
  pr_L        = new(m_number,float)
  pr_R        = new(m_number,float)

 PI=acos(-1.)

trop_lat = 30.
crit_lat=20.

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
  pr_I(j_f-1) = doubletofloat(lat(ITCZ))

;--------------------------------------------------------------------------------
; Find second ITCZ (if exists)

 
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

;  if (s_max_ind .gt. 0) then
;    pr_A(j_f-1) = precip(0,ITCZ,0) - precip(0,lat_max,0)
;  else
;    pr_A(j_f-1) = pr_D(j_f-1) 
;  end if 

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

 



;---------------------------------------------------------------------------------------
;PREDICTION WITH 9 PARAMETERS



 index_num = 9
 index_vec = new((/m_number,index_num/),float)


 index_vec(:,0) =   pr_W(:)
 index_vec(:,1) =   pr_E(:)
 index_vec(:,2) =   pr_A(:)
 index_vec(:,3) =   pr_S(:)
 index_vec(:,4) =   pr_I(:)
 index_vec(:,5) =   pr_D(:)
 index_vec(:,6) =   pr_H(:)
 index_vec(:,7) =   pr_L(:)
 index_vec(:,8) =   pr_R(:)


;-----------------------------------------------------------------------------------------------
 res_number=180
 pr_hr = new((/m_number,res_number/),float)
 lat_vec = new(res_number,float)
 do j_res=1,res_number
  lat_vec(j_res-1) = -90.+180./res_number/2. + 180./res_number*(j_res-1)
 end do 

  limit_lat =80.
  low_number= floattoint((-limit_lat-(-90.+180./res_number/2.))/180.*res_number)+1 
  high_number = res_number - low_number-1

 do j_f=1,m_number

  file1=addfile("../data_1979_2005/merged_files/zonmean/"+file_names(j_f-1)+"_tm_ZM.nc","r")
  
  lat=file1->lat
  lat_number = dimsizes(lat)
    
  do j_res=low_number,high_number
  
      pr_lat= lat_vec(j_res-1)
      do j_lat=2,lat_number-1
        if(lat(j_lat) .GE. pr_lat .AND. lat(j_lat-1) .LT. pr_lat) then 
          l_pr = j_lat-1   
        end if
      end do
    
      pr_hr(j_f-1,j_res-1)  = 86400.*(abs( ((pr_lat-doubletofloat(lat(l_pr)))*file1->pr(0,l_pr+1,0) + (doubletofloat(lat(l_pr+1))-pr_lat)*file1->pr(0,l_pr,0)) \
                        /(doubletofloat(lat(l_pr+1))-doubletofloat(lat(l_pr)))))
  
  end do

  delete(file1)
  delete(lat)
  delete(lat_number)
  delete(l_pr)
  delete(pr_lat)
 end do

;-------------------------------------------------------------------------------------------------------


 min_matrix = new((/index_num,index_num/),float)
 min_vec    = new(index_num,float)
 par_vec    = new((/res_number,index_num/),float)
 contrib_vec    = new((/index_num,res_number/),float)
 imp_precip = new((/m_number,res_number/),float)
 prec_error = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix = new((/index_num,index_num/),float)
  min_vec    =  new(index_num,float) 
  ind_mean =  new(index_num,float) 

  do l_ind = 1,index_num
    do k_ind = l_ind,index_num
       min_matrix(l_ind-1,k_ind-1) = sum(index_vec(:,l_ind-1)*index_vec(:,k_ind-1))
       min_matrix(k_ind-1,l_ind-1) = min_matrix(l_ind-1,k_ind-1)
    end do
    min_vec(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec(:,l_ind-1))
    ind_mean(l_ind-1) = sum(index_vec(:,l_ind-1))/m_number    
  end do  

  parameters = solve_linsys(min_matrix,min_vec)
  par_vec(j_res-1,:) = parameters(:)  
  contrib_vec(:,j_res-1) =  parameters(:)*ind_mean(:)

  do j_f = 1,m_number
    imp_precip(j_f-1,j_res-1) =sum(parameters(:)*index_vec(j_f-1,:))
    prec_error(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters)
  delete(min_vec)
  delete(min_matrix)
  delete(ind_mean)    

 end do

;---------------------------------------------------------------------------------------------------------------------------
;PREDICTION WITH 6 PARAMETERS

 index_num6 = 6
 index_vec6 = new((/m_number,index_num6/),float)

 index_vec6(:,0) =   pr_L(:)
 index_vec6(:,1) =   pr_E(:)
 index_vec6(:,2) =   pr_A(:)
 index_vec6(:,3) =   pr_S(:)
 index_vec6(:,4) =   pr_H(:)
 index_vec6(:,5) =   pr_D(:)


 min_matrix6 = new((/index_num6,index_num6/),float)
 min_vec6    = new(index_num6,float)
 par_vec6    = new((/res_number,index_num6/),float)
 contrib_vec6    = new((/index_num6,res_number/),float)
 imp_precip6 = new((/m_number,res_number/),float)
 prec_error6 = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix6 = new((/index_num6,index_num6/),float)
  min_vec6    =  new(index_num6,float) 
  ind_mean6 =  new(index_num6,float) 

  do l_ind = 1,index_num6
    do k_ind = l_ind,index_num6
       min_matrix6(l_ind-1,k_ind-1) = sum(index_vec6(:,l_ind-1)*index_vec6(:,k_ind-1))
       min_matrix6(k_ind-1,l_ind-1) = min_matrix6(l_ind-1,k_ind-1)
    end do
    min_vec6(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec6(:,l_ind-1))
    ind_mean6(l_ind-1) = sum(index_vec6(:,l_ind-1))/m_number    
  end do  

  parameters6 = solve_linsys(min_matrix6,min_vec6)
  par_vec6(j_res-1,:) = parameters6(:)  
  contrib_vec6(:,j_res-1) =  parameters6(:)*ind_mean6(:)

  do j_f = 1,m_number
    imp_precip6(j_f-1,j_res-1) =sum(parameters6(:)*index_vec6(j_f-1,:))
    prec_error6(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip6(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters6)
  delete(min_vec6)
  delete(min_matrix6)
  delete(ind_mean6)    

 end do



;----------------------------------------------------------------------------------------------------------
;PREDICTION WITH 4 PARAMETERS

 index_num4 = 4
 index_vec4 = new((/m_number,index_num4/),float)

 index_vec4(:,0) =   pr_W(:)
 index_vec4(:,1) =   pr_E(:)
 index_vec4(:,2) =   pr_D(:)
 index_vec4(:,3) =   pr_H(:)


 min_matrix4 = new((/index_num4,index_num4/),float)
 min_vec4    = new(index_num4,float)
 par_vec4    = new((/res_number,index_num4/),float)
 contrib_vec4    = new((/index_num4,res_number/),float)
 imp_precip4 = new((/m_number,res_number/),float)
 prec_error4 = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix4 = new((/index_num4,index_num4/),float)
  min_vec4    =  new(index_num4,float) 
  ind_mean4 =  new(index_num4,float) 

  do l_ind = 1,index_num4
    do k_ind = l_ind,index_num4
       min_matrix4(l_ind-1,k_ind-1) = sum(index_vec4(:,l_ind-1)*index_vec4(:,k_ind-1))
       min_matrix4(k_ind-1,l_ind-1) = min_matrix4(l_ind-1,k_ind-1)
    end do
    min_vec4(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec4(:,l_ind-1))
    ind_mean4(l_ind-1) = sum(index_vec4(:,l_ind-1))/m_number    
  end do  

  parameters4 = solve_linsys(min_matrix4,min_vec4)
  par_vec4(j_res-1,:) = parameters4(:)  
  contrib_vec4(:,j_res-1) =  parameters4(:)*ind_mean4(:)

  do j_f = 1,m_number
    imp_precip4(j_f-1,j_res-1) =sum(parameters4(:)*index_vec4(j_f-1,:))
    prec_error4(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip4(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters4)
  delete(min_vec4)
  delete(min_matrix4)
  delete(ind_mean4)    

 end do




;----------------------------------------------------------------------------------------------------------------------------
; PREDICTION 3 PARAMETERS


 index_num3 = 3
 index_vec3 = new((/m_number,index_num3/),float)

 index_vec3(:,0) =   pr_W(:)
 index_vec3(:,1) =   pr_E(:)
 index_vec3(:,2) =   pr_D(:)



 min_matrix3 = new((/index_num3,index_num3/),float)
 min_vec3    = new(index_num3,float)
 par_vec3    = new((/res_number,index_num3/),float)
 contrib_vec3    = new((/index_num3,res_number/),float)
 imp_precip3 = new((/m_number,res_number/),float)
 prec_error3 = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix3 = new((/index_num3,index_num3/),float)
  min_vec3    =  new(index_num3,float) 
  ind_mean3 =  new(index_num3,float) 

  do l_ind = 1,index_num3
    do k_ind = l_ind,index_num3
       min_matrix3(l_ind-1,k_ind-1) = sum(index_vec3(:,l_ind-1)*index_vec3(:,k_ind-1))
       min_matrix3(k_ind-1,l_ind-1) = min_matrix3(l_ind-1,k_ind-1)
    end do
    min_vec3(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec3(:,l_ind-1))
    ind_mean3(l_ind-1) = sum(index_vec3(:,l_ind-1))/m_number    
  end do  

  parameters3 = solve_linsys(min_matrix3,min_vec3)
  par_vec3(j_res-1,:) = parameters3(:)  
  contrib_vec3(:,j_res-1) =  parameters3(:)*ind_mean3(:)

  do j_f = 1,m_number
    imp_precip3(j_f-1,j_res-1) =sum(parameters3(:)*index_vec3(j_f-1,:))
    prec_error3(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip3(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters3)
  delete(min_vec3)
  delete(min_matrix3)
  delete(ind_mean3)    

 end do




;--------------------------------------------------------------------------------------------------------------------------------
 ;print(imp_precip(0,:))


 
   wks_type = "eps"
  wks_type@wkPSResolution = 100000

  wks   = gsn_open_wks (wks_type,"predicted_precip_AMIP_CMIP")
  plot = new(6,graphic)



;---------------------------------------------------------------------------------------------------------------
 
  variance_array = new ((/5,res_number/),float)
  correlation_array = new ((/4,res_number/),float)
 
  do j_res = low_number,high_number

    value_matrix_pr = new ((/2,m_number/),float)
    value_matrix_pr(0,:) = pr_hr(:,j_res-1) 
    value_matrix_pr(1,:) = pr_hr(:,j_res-1)
    cov_matrix_pr = covcorm(value_matrix_pr,(/0,1/))
    variance_array(4,j_res-1) = cov_matrix_pr(0,0)
 
    value_matrix_imp = new ((/2,m_number/),float)
    value_matrix_imp(0,:) = imp_precip(:,j_res-1)
    value_matrix_imp(1,:) = imp_precip(:,j_res-1)
    cov_matrix_imp = covcorm(value_matrix_imp,(/0,1/))
    variance_array(0,j_res-1) = cov_matrix_imp(0,0)

    value_matrix_imp3 = new ((/2,m_number/),float)
    value_matrix_imp3(0,:) = imp_precip3(:,j_res-1)
    value_matrix_imp3(1,:) = imp_precip3(:,j_res-1)
    cov_matrix_imp3 = covcorm(value_matrix_imp3,(/0,1/))
    variance_array(3,j_res-1) = cov_matrix_imp3(0,0)

    value_matrix_imp4 = new ((/2,m_number/),float)
    value_matrix_imp4(0,:) = imp_precip4(:,j_res-1)
    value_matrix_imp4(1,:) = imp_precip4(:,j_res-1)
    cov_matrix_imp4 = covcorm(value_matrix_imp4,(/0,1/))
    variance_array(2,j_res-1) = cov_matrix_imp4(0,0)

    value_matrix_imp6 = new ((/2,m_number/),float)
    value_matrix_imp6(0,:) = imp_precip6(:,j_res-1)
    value_matrix_imp6(1,:) = imp_precip6(:,j_res-1)
    cov_matrix_imp6 = covcorm(value_matrix_imp6,(/0,1/))
    variance_array(1,j_res-1) = cov_matrix_imp6(0,0)


    value_matrix_cor      = new ((/2,m_number/),float)
    value_matrix_cor(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor(1,:) = imp_precip(:,j_res-1)
    cov_matrix_cor = covcorm(value_matrix_cor,(/1,1/))
    correlation_array(0,j_res-1) = cov_matrix_cor(0,1)

    value_matrix_cor3      = new ((/2,m_number/),float)
    value_matrix_cor3(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor3(1,:) = imp_precip3(:,j_res-1)
    cov_matrix_cor3 = covcorm(value_matrix_cor3,(/1,1/))
    correlation_array(3,j_res-1) = cov_matrix_cor3(0,1)

    value_matrix_cor4      = new ((/2,m_number/),float)
    value_matrix_cor4(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor4(1,:) = imp_precip4(:,j_res-1)
    cov_matrix_cor4 = covcorm(value_matrix_cor4,(/1,1/))
    correlation_array(2,j_res-1) = cov_matrix_cor4(0,1)

    value_matrix_cor6      = new ((/2,m_number/),float)
    value_matrix_cor6(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor6(1,:) = imp_precip6(:,j_res-1)
    cov_matrix_cor6 = covcorm(value_matrix_cor6,(/1,1/))
    correlation_array(1,j_res-1) = cov_matrix_cor6(0,1)   
 
 
    delete(value_matrix_pr)
    delete(value_matrix_imp)
    delete(value_matrix_imp3)
    delete(value_matrix_imp4)
    delete(value_matrix_imp6)
    delete(value_matrix_cor)
    delete(value_matrix_cor3)
    delete(value_matrix_cor4)
    delete(value_matrix_cor6)

    delete(cov_matrix_pr)
    delete(cov_matrix_imp)
    delete(cov_matrix_imp3)
    delete(cov_matrix_imp4) 
    delete(cov_matrix_imp6)

    delete(cov_matrix_cor)
    delete(cov_matrix_cor3)
    delete(cov_matrix_cor4)
    delete(cov_matrix_cor6)
    
  end do

  res2   = True
  res2@tiMainString = "c"
  res2@tiMainFont = 21
  res2@tiMainFontHeightF = 0.03
  res2@tiYAxisString  = "StD. [mm day~S~-1~N~]"
  res2@tiYAxisFont    = 21
  res2@tiYAxisFontHeightF = 0.025 
  res2@tiXAxisFont    = 21
  res2@tiXAxisFontHeightF = 0.025
  res2@tmXBLabelsOn = False

  res2@tiMainOffsetYF = -0.018
  res2@tiMainOffsetXF = -0.13
  res2@vpWidthF          = 0.5
  res2@vpHeightF         = 0.5
  res2@tmYLLabelFontHeightF = 0.025 
  res2@tmXBLabelFontHeightF = 0.025

  res2@tmXBLabelFont= 21
  res2@tmXBMajorOutwardLengthF= 0.0
  res2@tmXBMinorOutwardLengthF= 0.0
  res2@tmXBMajorThicknessF= 2.0
  res2@tmYLLabelFont= 21
  res2@tmYLMajorOutwardLengthF= 0.0
  res2@tmYLMinorOutwardLengthF= 0.0
  res2@tmYLMajorThicknessF= 2.0

  res2@tmBorderThicknessF=  2.

  res2@trXMaxF     = 0.5
  res2@trXMinF     = -0.5

 ;res2@trYMaxF     = 200000.
 ;res2@trYMinF     = -200000.

  res2@trYMaxF     = 1.0
  res2@trYMinF     = 0.

  res2@gsnDraw         = False
  res2@gsnFrame         = False

  res2@xyLineColors   = (/"firebrick3","blue","orange","black","gray34"/)
  res2@xyDashPatterns = (/0,1,3,2,0/)
  res2@xyLineThicknesses =(/2.,2.,2.,2.,2./)

  res2@tmXBMode = "explicit"
  res2@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  res2@tmXBLabels = (/-60,-30,0,30,60/)
  res2@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)
 
  plot(2)=gsn_csm_xy(wks,sin(lat_vec(:)/180.*PI),variance_array(:,:),res2)
;---------------------------------------------------------------------------------------------------------------

 RMSE_vec = new ((/4,res_number/),float) 

 do j_res = low_number,high_number
  RMSE_vec(0,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip(:,j_res-1))^2)/m_number)
  RMSE_vec(3,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip3(:,j_res-1))^2)/m_number)
  RMSE_vec(2,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip4(:,j_res-1))^2)/m_number)
  RMSE_vec(1,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip6(:,j_res-1))^2)/m_number)
 end do


   RMSE_vec_int = new (4,float) 
   RMSE_vec_int(:) = 0.
   do j_res = 1,res_number
     if ( lat_vec(j_res-1) .gt. -20. .and. lat_vec(j_res-1) .lt. 20.) then
         RMSE_vec_int(:) = RMSE_vec_int(:) + \
                           2.*PI*(RMSE_vec(:,j_res-1)*cos(lat_vec(j_res-1)/180.*PI) + RMSE_vec(:,j_res)*doubletofloat(cos(lat_vec(j_res)/180.*PI)))/2. * \
                           ((lat_vec(j_res-1)+ lat_vec(j_res))/2.-(lat_vec(j_res-2)+ lat_vec(j_res-1))/2. )/180.*PI                                             
     end if 
   end do
   RMSE_vec_int(:) = RMSE_vec_int(:)/doubletofloat(2.*PI*(sin(20.*PI/180.) - sin(-20.*PI/180.)))
   print("RMSE AMIP:")
   print(RMSE_vec_int)
   

  res3   = True
  res3@tiMainString = "a"
  res3@tiMainFont = 21
  res3@tiMainFontHeightF = 0.03
  res3@tiYAxisString  = "RMSE [mm day~S~-1~N~]"
  res3@tiYAxisFont    = 21
  res3@tiYAxisFontHeightF = 0.025 
  res3@tiXAxisFont    = 21
  res3@tiXAxisFontHeightF = 0.025
  res3@tmXBLabelsOn = False

  res3@tiMainOffsetYF = -0.018
  res3@tiMainOffsetXF = -0.13
  res3@vpWidthF          = 0.5
  res3@vpHeightF         = 0.5
  res3@tmYLLabelFontHeightF = 0.025 
  res3@tmXBLabelFontHeightF = 0.025

  res3@tmXBLabelFont= 21
  res3@tmXBMajorOutwardLengthF= 0.0
  res3@tmXBMinorOutwardLengthF= 0.0
  res3@tmXBMajorThicknessF= 2.0
  res3@tmYLLabelFont= 21
  res3@tmYLMajorOutwardLengthF= 0.0
  res3@tmYLMinorOutwardLengthF= 0.0
  res3@tmYLMajorThicknessF= 2.0

  res3@tmBorderThicknessF=  2.

  res3@trXMaxF     = 0.5
  res3@trXMinF     = -0.5

 ;res3@trYMaxF     = 200000.
 ;res3@trYMinF     = -200000.

  res3@trYMaxF     = 0.7
  res3@trYMinF     = 0.

  res3@gsnDraw         = False
  res3@gsnFrame         = False


  res3@xyLineColors   = (/"firebrick3","blue","orange","black"/)
  res3@xyDashPatterns = (/0,1,3,2/)
  res3@xyLineThicknesses =(/2.,2.,2.,2./)

  res3@tmXBMode = "explicit"
  res3@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  res3@tmXBLabels = (/-60,-30,0,30,60/)
  res3@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)
 
  plot(0)=gsn_csm_xy(wks,sin(lat_vec(:)/180.*PI),RMSE_vec(:,:),res3)

;-----------------------------------------------------------------------------------------------------------
  res4   = True
  res4@tiMainFont = 21
  res4@tiMainString = "e"
  res4@tiMainFontHeightF = 0.03
  res4@tiYAxisString  = "Correlation"
  res4@tiYAxisFont    = 21
  res4@tiYAxisFontHeightF = 0.025 
  res4@tiXAxisString  = "Latitude"
  res4@tiXAxisFont    = 21
  res4@tiXAxisFontHeightF = 0.025

  res4@tiMainOffsetYF = -0.018
  res4@tiMainOffsetXF = -0.13
  res4@vpWidthF          = 0.5
  res4@vpHeightF         = 0.5
  res4@tmYLLabelFontHeightF = 0.025 
  res4@tmXBLabelFontHeightF = 0.025

  res4@tmXBLabelFont= 21
  res4@tmXBMajorOutwardLengthF= 0.0
  res4@tmXBMinorOutwardLengthF= 0.0
  res4@tmXBMajorThicknessF= 2.0
  res4@tmYLLabelFont= 21
  res4@tmYLMajorOutwardLengthF= 0.0
  res4@tmYLMinorOutwardLengthF= 0.0
  res4@tmYLMajorThicknessF= 2.0

  res4@tmBorderThicknessF=  2.

  res4@trXMaxF     = 0.5
  res4@trXMinF     = -0.5

 ;res4@trYMaxF     = 200000.
 ;res4@trYMinF     = -200000.

  res4@trYMaxF     = 1.
  res4@trYMinF     = 0.

  res4@gsnDraw         = False
  res4@gsnFrame         = False


  res4@xyLineColors   = (/"firebrick3","blue","orange","black","gray34"/)
  res4@xyDashPatterns = (/0,1,3,2,0/)
  res4@xyLineThicknesses =(/2.,2.,2.,2.,2./)

  res4@tmXBMode = "explicit"
  res4@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  res4@tmXBLabels = (/-60,-30,0,30,60/)
  res4@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)
 
  plot(4)=gsn_csm_xy(wks,sin(lat_vec(:)/180.*PI),correlation_array(:,:),res4)
  





;-----------------------------------------------------------------------------------------------------------------

  wks2   = gsn_open_wks (wks_type,"fitting_func_9")
  plot2 = new(4,graphic)

  prec_func = new ((/5,res_number/),float)
  lat_func  = new ((/4,res_number/),float)

; index_vec(:,0) =   pr_W(:)
; index_vec(:,1) =   pr_E(:)
; index_vec(:,2) =   pr_A(:)
; index_vec(:,3) =   pr_S(:)
; index_vec(:,4) =   pr_I(:)
; index_vec(:,5) =   pr_D(:)
; index_vec(:,6) =   pr_H(:)
; index_vec(:,7) =   pr_L(:)
; index_vec(:,8) =   pr_R(:)


;print("A")



  prec_func(0,:) = par_vec(:,1)
  prec_func(1,:) = par_vec(:,2)
  prec_func(2,:) = par_vec(:,5)
  prec_func(3,:) = par_vec(:,6)
  prec_func(4,:) = par_vec(:,8)

  lat_func(0,:) = par_vec(:,0)
  lat_func(1,:) = par_vec(:,3)
  lat_func(2,:) = par_vec(:,4)
  lat_func(3,:) = par_vec(:,7)
  
  ;print(lat_func(0,:))

  resp   = True
  resp@tiMainString = "Zonal-mean P"
  resp@tiMainFont = 21
  resp@tiMainFontHeightF = 0.02
  resp@tiYAxisString  = "P [mm day~S~-1~N~]"
  resp@tiYAxisFont    = 21
  resp@tiYAxisFontHeightF = 0.02 
  resp@tiXAxisString  = "Latitude"
  resp@tiXAxisFont    = 21
  resp@tiXAxisFontHeightF = 0.02

  resp@tiMainOffsetYF = -0.01
  resp@tiMainOffsetXF = -0.07
  resp@vpWidthF          = 0.5
  resp@vpHeightF         = 0.5
  resp@tmYLLabelFontHeightF = 0.02 
  resp@tmXBLabelFontHeightF = 0.02

  resp@tmXBLabelFont= 21
  resp@tmXBMajorOutwardLengthF= 0.0
  resp@tmXBMinorOutwardLengthF= 0.0
  resp@tmXBMajorThicknessF= 2.0
  resp@tmYLLabelFont= 21
  resp@tmYLMajorOutwardLengthF= 0.0
  resp@tmYLMinorOutwardLengthF= 0.0
  resp@tmYLMajorThicknessF= 2.0

  resp@tmBorderThicknessF=  2.

  resp@trXMaxF     = 0.5
  resp@trXMinF     = -0.5

 resp@trYMaxF     = 2.
 resp@trYMinF     = -2.


  resp@gsnDraw         = False
  resp@gsnFrame         = False

  resp@xyLineColors   = (/"red","blue","grey34","purple4","yellow"/)
  resp@xyDashPatterns = (/1,1,1,1/)

  resp@tmXBMode = "explicit"
  resp@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  resp@tmXBLabels = (/-60,-30,0,30,60/)
  resp@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)


 ;print(sin(lat_vec(:)/180.*PI))

  plot2(0)=gsn_csm_xy(wks2,sin(lat_vec(:)/180.*PI),prec_func(:,:),resp)


;----------------------------------------------------------------------------------------------------
;print("B")

  resl   = True
  resl@tiMainString = "Zonal-mean P"
  resl@tiMainFont = 21
  resl@tiMainFontHeightF = 0.02
  resl@tiYAxisString  = "P [mm day~S~-1~N~]"
  resl@tiYAxisFont    = 21
  resl@tiYAxisFontHeightF = 0.02 
  resl@tiXAxisString  = "Latitude"
  resl@tiXAxisFont    = 21
  resl@tiXAxisFontHeightF = 0.02

  resl@tiMainOffsetYF = -0.01
  resl@tiMainOffsetXF = -0.07
  resl@vpWidthF          = 0.5
  resl@vpHeightF         = 0.5
  resl@tmYLLabelFontHeightF = 0.02 
  resl@tmXBLabelFontHeightF = 0.02

  resl@tmXBLabelFont= 21
  resl@tmXBMajorOutwardLengthF= 0.0
  resl@tmXBMinorOutwardLengthF= 0.0
  resl@tmXBMajorThicknessF= 2.0
  resl@tmYLLabelFont= 21
  resl@tmYLMajorOutwardLengthF= 0.0
  resl@tmYLMinorOutwardLengthF= 0.0
  resl@tmYLMajorThicknessF= 2.0

  resl@tmBorderThicknessF=  2.

  resl@trXMaxF     = 0.5
  resl@trXMinF     = -0.5

 resl@trYMaxF     = 2.
 resl@trYMinF     = -2.


  resl@gsnDraw         = False
  resl@gsnFrame         = False


  resl@xyLineColors   = (/"seagreen","black","orange","black"/)
  resl@xyDashPatterns = (/0,0,0,0/)

  resl@tmXBMode = "explicit"
  resl@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  resl@tmXBLabels = (/-60,-30,0,30,60/)
  resl@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)

  plot2(1)=gsn_csm_xy(wks2,sin(lat_vec(:)/180.*PI),lat_func(:,:),resl)
;-----------------------------------------------------------------------------------------------------
 ;print("C")
  
  resc   = True
  resc@tiMainString = "Zonal-mean P"
  resc@tiMainFont = 21
  resc@tiMainFontHeightF = 0.02
  resc@tiYAxisString  = "P [mm day~S~-1~N~]"
  resc@tiYAxisFont    = 21
  resc@tiYAxisFontHeightF = 0.02 
  resc@tiXAxisString  = "Latitude"
  resc@tiXAxisFont    = 21
  resc@tiXAxisFontHeightF = 0.02

  resc@tiMainOffsetYF = -0.01
  resc@tiMainOffsetXF = -0.07
  resc@vpWidthF          = 0.5
  resc@vpHeightF         = 0.5
  resc@tmYLLabelFontHeightF = 0.02 
  resc@tmXBLabelFontHeightF = 0.02

  resc@tmXBLabelFont= 21
  resc@tmXBMajorOutwardLengthF= 0.0
  resc@tmXBMinorOutwardLengthF= 0.0
  resc@tmXBMajorThicknessF= 2.0
  resc@tmYLLabelFont= 21
  resc@tmYLMajorOutwardLengthF= 0.0
  resc@tmYLMinorOutwardLengthF= 0.0
  resc@tmYLMajorThicknessF= 2.0

  resc@tmBorderThicknessF=  2.

  resc@trXMaxF     = 0.5
  resc@trXMinF     = -0.5

  resc@trYMaxF     = 6.
  resc@trYMinF     = -2.


  resc@gsnDraw         = False
  resc@gsnFrame         = False


  resc@xyLineColors   = (/"seagreen","red","blue","black","orange","grey34","purple4","black","yellow"/)
  resc@xyDashPatterns = (/0,1,1,0,0,1,1,0,1/)

  resc@tmXBMode = "explicit"
  resc@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  resc@tmXBLabels = (/-60,-30,0,30,60/)
  resc@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)

  ;print(res_number)

 

;print("Ca")
  plot2(2)=gsn_csm_xy(wks2,sin(lat_vec(14:(res_number-15))/180.*PI),contrib_vec(:,14:(res_number-15)),resc)
;print("D")




delete(m_number)
delete(v_number)
delete(precip_eq)

delete(pr_W)
delete(pr_E)
delete(pr_A)
delete(pr_S)
delete(pr_I)
delete(pr_D)
delete(pr_H)
delete(pr_L)
delete(pr_R)

delete(index_num)
delete(index_vec)

delete(res_number)
delete(pr_hr)
delete(lat_vec)


delete(par_vec)
delete(contrib_vec)
delete(imp_precip)
delete(prec_error)


delete(index_vec3)
delete(par_vec3)
delete(contrib_vec3)
delete(imp_precip3)
delete(prec_error3)

delete(index_vec4)
delete(par_vec4)
delete(contrib_vec4)
delete(imp_precip4)
delete(prec_error4)

delete(index_vec6)
delete(par_vec6)
delete(contrib_vec6)
delete(imp_precip6)
delete(prec_error6)


delete(variance_array)
delete(correlation_array)
delete(RMSE_vec)
delete(RMSE_vec_int)

delete(res2)
delete(res3)
delete(res4)











;--------------------------------------------------------------------------------------------------------------------------
;CMIP
;--------------------------------------------------------------------------------------------------------------------------


file_names_C=new (26,string)

file_names_C(0)="Amon_ACCESS1-0_historical_r1i1p1_197901-200512"
file_names_C(1)="Amon_ACCESS1-3_historical_r1i1p1_197901-200512"
file_names_C(2)="Amon_bcc-csm1-1_historical_r1i1p1_197901-200512"
file_names_C(3)="Amon_BNU-ESM_historical_r1i1p1_197901-200512"

file_names_C(4)="Amon_CanESM2_historical_r1i1p1_197901-200512"
file_names_C(5)="Amon_CCSM4_historical_r1i1p1_197901-200512"
file_names_C(6)="Amon_CESM1-BGC_historical_r1i1p1_197901-200512"
file_names_C(7)="Amon_CNRM-CM5_historical_r1i1p1_197901-200512"

file_names_C(8)="Amon_CSIRO-Mk3-6-0_historical_r1i1p1_197901-200512"
file_names_C(9)="Amon_FGOALS-s2_historical_r1i1p1_197901-200512"
file_names_C(10)="Amon_GFDL-CM3_form_historical_r1i1p1_197901-200512"
file_names_C(11)="Amon_GFDL-ESM2G_form_historical_r1i1p1_197901-200512"

file_names_C(12)="Amon_GFDL-ESM2M_form_historical_r1i1p1_197901-200512"
file_names_C(13)="Amon_GISS-E2-R_historical_r1i1p1_197901-200512"
file_names_C(14)="Amon_GISS-E2-H_historical_r1i1p1_197901-200512"
file_names_C(15)="Amon_HadGEM2-ES_historical_r1i1p1_197901-200512"

file_names_C(16)="Amon_inmcm4_historical_r1i1p1_197901-200512"
file_names_C(17)="Amon_IPSL-CM5A-LR_historical_r1i1p1_197901-200512"
file_names_C(18)="Amon_IPSL-CM5A-MR_historical_r1i1p1_197901-200512"
file_names_C(19)="Amon_IPSL-CM5B-LR_historical_r1i1p1_197901-200512"

file_names_C(20)="Amon_MIROC5_historical_r1i1p1_197901-200512"
file_names_C(21)="Amon_MIROC-ESM_historical_r1i1p1_197901-200512"
file_names_C(22)="Amon_MPI-ESM-LR_historical_r1i1p1_197901-200512"
file_names_C(23)="Amon_MPI-ESM-MR_historical_r1i1p1_197901-200512"

file_names_C(24)="Amon_MRI-CGCM3_historical_r1i1p1_197901-200512"
file_names_C(25)="Amon_NorESM1-M_historical_r1i1p1_197901-200512"

print("CMIP Names")
print(file_names_C)

m_number=dimsizes(file_names_C)
v_number=dimsizes(file_names_C)

precip_eq = new(m_number,float)

  pr_W        = new(m_number,float)
  pr_H        = new(m_number,float)
  pr_E        = new(m_number,float)
  pr_A        = new(m_number,float)
  pr_S        = new(m_number,float)
  pr_I        = new(m_number,float)
  pr_D        = new(m_number,float)
  pr_L        = new(m_number,float)
  pr_R        = new(m_number,float)

 PI=acos(-1.)

trop_lat = 30.
crit_lat=20.

do j_f=1,m_number

 file1=addfile("../../../HISTORICAL/plots/data_1979_2005/merged_files/zonmean/"+file_names_C(j_f-1)+"_tm_ZM.nc","r")

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
  pr_I(j_f-1) = doubletofloat(lat(ITCZ))

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

 



;---------------------------------------------------------------------------------------
;PREDICTION WITH 9 PARAMETERS



 index_num = 9
 index_vec = new((/m_number,index_num/),float)


 index_vec(:,0) =   pr_W(:)
 index_vec(:,1) =   pr_E(:)
 index_vec(:,2) =   pr_A(:)
 index_vec(:,3) =   pr_S(:)
 index_vec(:,4) =   pr_I(:)
 index_vec(:,5) =   pr_D(:)
 index_vec(:,6) =   pr_H(:)
 index_vec(:,7) =   pr_L(:)
 index_vec(:,8) =   pr_R(:)


;-----------------------------------------------------------------------------------------------
 res_number=180
 pr_hr = new((/m_number,res_number/),float)
 lat_vec = new(res_number,float)
 do j_res=1,res_number
  lat_vec(j_res-1) = -90.+180./res_number/2. + 180./res_number*(j_res-1)
 end do 

  limit_lat =80.
  low_number= floattoint((-limit_lat-(-90.+180./res_number/2.))/180.*res_number)+1 
  high_number = res_number - low_number-1

 do j_f=1,m_number

  file1=addfile("../../../HISTORICAL/plots/data_1979_2005/merged_files/zonmean/"+file_names_C(j_f-1)+"_tm_ZM.nc","r")
  
  lat=file1->lat
  lat_number = dimsizes(lat)
    
  do j_res=low_number,high_number
  
      pr_lat= lat_vec(j_res-1)
      do j_lat=2,lat_number-1
        if(lat(j_lat) .GE. pr_lat .AND. lat(j_lat-1) .LT. pr_lat) then 
          l_pr = j_lat-1   
        end if
      end do
    
      pr_hr(j_f-1,j_res-1)  = 86400.*(abs( ((pr_lat-doubletofloat(lat(l_pr)))*file1->pr(0,l_pr+1,0) + (doubletofloat(lat(l_pr+1))-pr_lat)*file1->pr(0,l_pr,0)) \
                        /(doubletofloat(lat(l_pr+1))-doubletofloat(lat(l_pr)))))
  
  end do

  delete(file1)
  delete(lat)
  delete(lat_number)
  delete(l_pr)
  delete(pr_lat)
 end do

;-------------------------------------------------------------------------------------------------------


 min_matrix = new((/index_num,index_num/),float)
 min_vec    = new(index_num,float)
 par_vec    = new((/res_number,index_num/),float)
 contrib_vec    = new((/index_num,res_number/),float)
 imp_precip = new((/m_number,res_number/),float)
 prec_error = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix = new((/index_num,index_num/),float)
  min_vec    =  new(index_num,float) 
  ind_mean =  new(index_num,float) 

  do l_ind = 1,index_num
    do k_ind = l_ind,index_num
       min_matrix(l_ind-1,k_ind-1) = sum(index_vec(:,l_ind-1)*index_vec(:,k_ind-1))
       min_matrix(k_ind-1,l_ind-1) = min_matrix(l_ind-1,k_ind-1)
    end do
    min_vec(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec(:,l_ind-1))
    ind_mean(l_ind-1) = sum(index_vec(:,l_ind-1))/m_number    
  end do  

  parameters = solve_linsys(min_matrix,min_vec)
  par_vec(j_res-1,:) = parameters(:)  
  contrib_vec(:,j_res-1) =  parameters(:)*ind_mean(:)

  do j_f = 1,m_number
    imp_precip(j_f-1,j_res-1) =sum(parameters(:)*index_vec(j_f-1,:))
    prec_error(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters)
  delete(min_vec)
  delete(min_matrix)
  delete(ind_mean)    

 end do

;---------------------------------------------------------------------------------------------------------------------------
;PREDICTION WITH 6 PARAMETERS

 index_num6 = 6
 index_vec6 = new((/m_number,index_num6/),float)

 index_vec6(:,0) =   pr_W(:)
 index_vec6(:,1) =   pr_E(:)
 index_vec6(:,2) =   pr_A(:)
 index_vec6(:,3) =   pr_S(:)
 index_vec6(:,4) =   pr_H(:)
 index_vec6(:,5) =   pr_D(:)


 min_matrix6 = new((/index_num6,index_num6/),float)
 min_vec6    = new(index_num6,float)
 par_vec6    = new((/res_number,index_num6/),float)
 contrib_vec6    = new((/index_num6,res_number/),float)
 imp_precip6 = new((/m_number,res_number/),float)
 prec_error6 = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix6 = new((/index_num6,index_num6/),float)
  min_vec6    =  new(index_num6,float) 
  ind_mean6 =  new(index_num6,float) 

  do l_ind = 1,index_num6
    do k_ind = l_ind,index_num6
       min_matrix6(l_ind-1,k_ind-1) = sum(index_vec6(:,l_ind-1)*index_vec6(:,k_ind-1))
       min_matrix6(k_ind-1,l_ind-1) = min_matrix6(l_ind-1,k_ind-1)
    end do
    min_vec6(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec6(:,l_ind-1))
    ind_mean6(l_ind-1) = sum(index_vec6(:,l_ind-1))/m_number    
  end do  

  parameters6 = solve_linsys(min_matrix6,min_vec6)
  par_vec6(j_res-1,:) = parameters6(:)  
  contrib_vec6(:,j_res-1) =  parameters6(:)*ind_mean6(:)

  do j_f = 1,m_number
    imp_precip6(j_f-1,j_res-1) =sum(parameters6(:)*index_vec6(j_f-1,:))
    prec_error6(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip6(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters6)
  delete(min_vec6)
  delete(min_matrix6)
  delete(ind_mean6)    

 end do


;---------------------------------------------------------------------------------------------------------------------------
;PREDICTION WITH 4 PARAMETERS

 index_num4 = 4
 index_vec4 = new((/m_number,index_num4/),float)

 index_vec4(:,0) =   pr_S(:)
 index_vec4(:,1) =   pr_E(:)
 index_vec4(:,2) =   pr_A(:)
 index_vec4(:,3) =   pr_D(:)


 min_matrix4 = new((/index_num4,index_num4/),float)
 min_vec4    = new(index_num4,float)
 par_vec4    = new((/res_number,index_num4/),float)
 contrib_vec4    = new((/index_num4,res_number/),float)
 imp_precip4 = new((/m_number,res_number/),float)
 prec_error4 = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix4 = new((/index_num4,index_num4/),float)
  min_vec4    =  new(index_num4,float) 
  ind_mean4 =  new(index_num4,float) 

  do l_ind = 1,index_num4
    do k_ind = l_ind,index_num4
       min_matrix4(l_ind-1,k_ind-1) = sum(index_vec4(:,l_ind-1)*index_vec4(:,k_ind-1))
       min_matrix4(k_ind-1,l_ind-1) = min_matrix4(l_ind-1,k_ind-1)
    end do
    min_vec4(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec4(:,l_ind-1))
    ind_mean4(l_ind-1) = sum(index_vec4(:,l_ind-1))/m_number    
  end do  

  parameters4 = solve_linsys(min_matrix4,min_vec4)
  par_vec4(j_res-1,:) = parameters4(:)  
  contrib_vec4(:,j_res-1) =  parameters4(:)*ind_mean4(:)

  do j_f = 1,m_number
    imp_precip4(j_f-1,j_res-1) =sum(parameters4(:)*index_vec4(j_f-1,:))
    prec_error4(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip4(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters4)
  delete(min_vec4)
  delete(min_matrix4)
  delete(ind_mean4)    

 end do










;----------------------------------------------------------------------------------------------------------------------------
; PREDICTION 3 PARAMETERS


 index_num3 = 3
 index_vec3 = new((/m_number,index_num3/),float)

 index_vec3(:,0) =   pr_S(:)
 index_vec3(:,1) =   pr_E(:)
 index_vec3(:,2) =   pr_D(:)



 min_matrix3 = new((/index_num3,index_num3/),float)
 min_vec3    = new(index_num3,float)
 par_vec3    = new((/res_number,index_num3/),float)
 contrib_vec3    = new((/index_num3,res_number/),float)
 imp_precip3 = new((/m_number,res_number/),float)
 prec_error3 = new((/m_number,res_number/),float)

 do j_res = 1,res_number

  min_matrix3 = new((/index_num3,index_num3/),float)
  min_vec3    =  new(index_num3,float) 
  ind_mean3 =  new(index_num3,float) 

  do l_ind = 1,index_num3
    do k_ind = l_ind,index_num3
       min_matrix3(l_ind-1,k_ind-1) = sum(index_vec3(:,l_ind-1)*index_vec3(:,k_ind-1))
       min_matrix3(k_ind-1,l_ind-1) = min_matrix3(l_ind-1,k_ind-1)
    end do
    min_vec3(l_ind-1) = sum(pr_hr(:,j_res-1)*index_vec3(:,l_ind-1))
    ind_mean3(l_ind-1) = sum(index_vec3(:,l_ind-1))/m_number    
  end do  

  parameters3 = solve_linsys(min_matrix3,min_vec3)
  par_vec3(j_res-1,:) = parameters3(:)  
  contrib_vec3(:,j_res-1) =  parameters3(:)*ind_mean3(:)

  do j_f = 1,m_number
    imp_precip3(j_f-1,j_res-1) =sum(parameters3(:)*index_vec3(j_f-1,:))
    prec_error3(j_f-1,j_res-1) = (pr_hr(j_f-1,j_res-1)-imp_precip3(j_f-1,j_res-1))/pr_hr(j_f-1,j_res-1)
  end do

  delete(parameters3)
  delete(min_vec3)
  delete(min_matrix3)
  delete(ind_mean3)    

 end do




;---------------------------------------------------------------------------------------------------------------
 
  variance_array = new ((/5,res_number/),float)
  correlation_array = new ((/4,res_number/),float)
 
  do j_res = low_number,high_number

    value_matrix_pr = new ((/2,m_number/),float)
    value_matrix_pr(0,:) = pr_hr(:,j_res-1) 
    value_matrix_pr(1,:) = pr_hr(:,j_res-1)
    cov_matrix_pr = covcorm(value_matrix_pr,(/0,1/))
    variance_array(4,j_res-1) = cov_matrix_pr(0,0)
 
    value_matrix_imp = new ((/2,m_number/),float)
    value_matrix_imp(0,:) = imp_precip(:,j_res-1)
    value_matrix_imp(1,:) = imp_precip(:,j_res-1)
    cov_matrix_imp = covcorm(value_matrix_imp,(/0,1/))
    variance_array(0,j_res-1) = cov_matrix_imp(0,0)

    value_matrix_imp3 = new ((/2,m_number/),float)
    value_matrix_imp3(0,:) = imp_precip3(:,j_res-1)
    value_matrix_imp3(1,:) = imp_precip3(:,j_res-1)
    cov_matrix_imp3 = covcorm(value_matrix_imp3,(/0,1/))
    variance_array(3,j_res-1) = cov_matrix_imp3(0,0)

    value_matrix_imp4 = new ((/2,m_number/),float)
    value_matrix_imp4(0,:) = imp_precip4(:,j_res-1)
    value_matrix_imp4(1,:) = imp_precip4(:,j_res-1)
    cov_matrix_imp4 = covcorm(value_matrix_imp4,(/0,1/))
    variance_array(2,j_res-1) = cov_matrix_imp4(0,0)

    value_matrix_imp6 = new ((/2,m_number/),float)
    value_matrix_imp6(0,:) = imp_precip6(:,j_res-1)
    value_matrix_imp6(1,:) = imp_precip6(:,j_res-1)
    cov_matrix_imp6 = covcorm(value_matrix_imp6,(/0,1/))
    variance_array(1,j_res-1) = cov_matrix_imp6(0,0)


    value_matrix_cor      = new ((/2,m_number/),float)
    value_matrix_cor(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor(1,:) = imp_precip(:,j_res-1)
    cov_matrix_cor = covcorm(value_matrix_cor,(/1,1/))
    correlation_array(0,j_res-1) = cov_matrix_cor(0,1)

    value_matrix_cor3      = new ((/2,m_number/),float)
    value_matrix_cor3(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor3(1,:) = imp_precip3(:,j_res-1)
    cov_matrix_cor3 = covcorm(value_matrix_cor3,(/1,1/))
    correlation_array(3,j_res-1) = cov_matrix_cor3(0,1)

    value_matrix_cor4      = new ((/2,m_number/),float)
    value_matrix_cor4(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor4(1,:) = imp_precip4(:,j_res-1)
    cov_matrix_cor4 = covcorm(value_matrix_cor4,(/1,1/))
    correlation_array(2,j_res-1) = cov_matrix_cor4(0,1)

    value_matrix_cor6      = new ((/2,m_number/),float)
    value_matrix_cor6(0,:) = pr_hr(:,j_res-1)
    value_matrix_cor6(1,:) = imp_precip6(:,j_res-1)
    cov_matrix_cor6 = covcorm(value_matrix_cor6,(/1,1/))
    correlation_array(1,j_res-1) = cov_matrix_cor6(0,1)   
 
 
    delete(value_matrix_pr)
    delete(value_matrix_imp)
    delete(value_matrix_imp3)
    delete(value_matrix_imp4)
    delete(value_matrix_imp6)
    delete(value_matrix_cor)
    delete(value_matrix_cor3)
    delete(value_matrix_cor4)
    delete(value_matrix_cor6)

    delete(cov_matrix_pr)
    delete(cov_matrix_imp)
    delete(cov_matrix_imp3) 
    delete(cov_matrix_imp4) 
    delete(cov_matrix_imp6)

    delete(cov_matrix_cor)
    delete(cov_matrix_cor3) 
    delete(cov_matrix_cor4) 
    delete(cov_matrix_cor6)
    
  end do

  res2   = True
  res2@tiMainString = "d"
  res2@tiMainFont = 21
  res2@tiMainFontHeightF = 0.03
  res2@tiYAxisFont    = 21
  res2@tiYAxisFontHeightF = 0.025
  res2@tmYLLabelsOn = False 
  res2@tiXAxisFont    = 21
  res2@tiXAxisFontHeightF = 0.025
  res2@tmXBLabelsOn = False

  res2@tiMainOffsetYF = -0.01
  res2@tiMainOffsetXF = -0.13
  res2@vpWidthF          = 0.5
  res2@vpHeightF         = 0.5
  res2@tmYLLabelFontHeightF = 0.025 
  res2@tmXBLabelFontHeightF = 0.025

  res2@tmXBLabelFont= 21
  res2@tmXBMajorOutwardLengthF= 0.0
  res2@tmXBMinorOutwardLengthF= 0.0
  res2@tmXBMajorThicknessF= 2.0
  res2@tmYLLabelFont= 21
  res2@tmYLMajorOutwardLengthF= 0.0
  res2@tmYLMinorOutwardLengthF= 0.0
  res2@tmYLMajorThicknessF= 2.0

  res2@tmBorderThicknessF=  2.

  res2@trXMaxF     = 0.5
  res2@trXMinF     = -0.5

 ;res2@trYMaxF     = 200000.
 ;res2@trYMinF     = -200000.

  res2@trYMaxF     = 1.0
  res2@trYMinF     = 0.

  res2@gsnDraw         = False
  res2@gsnFrame         = False

  res2@xyLineColors   = (/"firebrick3","blue","orange","black","gray34"/)
  res2@xyDashPatterns = (/0,1,3,2,0/)
  res2@xyLineThicknesses =(/2.,2.,2.,2.,2./)


  res2@xyExplicitLegendLabels = (/" 9 indices"," 6 indices"," 4 indices"," 3 indices"," multi-model ensemble"/)

  res2@pmLegendDisplayMode    = "Always" 
  res2@lgJustification        = "TopRight"
 ;res2@trYMinF                = -3.5      ; Leave space at bottom for legend
  res2@pmLegendOrthogonalPosF = -0.95     ;  Move into plot
  res2@pmLegendParallelPosF   = -0.3       ; Move to right
  res2@pmLegendWidthF         = 0.3       ; Decrease width
  res2@pmLegendHeightF        = 0.3       ; Decrease height
  res2@lgBoxMinorExtentF      = 0.1       ; Shorten the legend lines
  res2@lgLabelFontHeightF    = 0.023     ; Change the font size
  res2@lgLabelFont           = 21
  res2@lgPerimOn              = False
  res2@lgLineThicknessF      = 5. 


  res2@tmXBMode = "explicit"
  res2@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  res2@tmXBLabels = (/-60,-30,0,30,60/)
  res2@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)
 
  plot(3)=gsn_csm_xy(wks,sin(lat_vec(:)/180.*PI),variance_array(:,:),res2)
;---------------------------------------------------------------------------------------------------------------

 RMSE_vec = new ((/4,res_number/),float) 

 do j_res = low_number,high_number
  RMSE_vec(0,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip(:,j_res-1))^2)/m_number)
  RMSE_vec(3,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip3(:,j_res-1))^2)/m_number)
  RMSE_vec(2,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip4(:,j_res-1))^2)/m_number)
  RMSE_vec(1,j_res-1) = sqrt(sum((pr_hr(:,j_res-1) - imp_precip6(:,j_res-1))^2)/m_number)
 end do


   RMSE_vec_int = new (4,float) 
   RMSE_vec_int(:) = 0.
   do j_res = 1,res_number
     if ( lat_vec(j_res-1) .gt. -20. .and. lat_vec(j_res-1) .lt. 20.) then
         RMSE_vec_int(:) = RMSE_vec_int(:) + \
                           2.*PI*(RMSE_vec(:,j_res-1)*cos(lat_vec(j_res-1)/180.*PI) + RMSE_vec(:,j_res)*doubletofloat(cos(lat_vec(j_res)/180.*PI)))/2. * \
                           ((lat_vec(j_res-1)+ lat_vec(j_res))/2.-(lat_vec(j_res-2)+ lat_vec(j_res-1))/2. )/180.*PI                                             
     end if 
   end do
   RMSE_vec_int(:) = RMSE_vec_int(:)/doubletofloat(2.*PI*(sin(20.*PI/180.) - sin(-20.*PI/180.)))
   print("RMSE CMIP:")
   print(RMSE_vec_int)



  res3   = True
  res3@tiMainString = "b"
  res3@tiMainFont = 21
  res3@tiMainFontHeightF = 0.03
  res3@tiYAxisFont    = 21
  res3@tmYLLabelsOn = False 
  res3@tiYAxisFontHeightF = 0.025 
  res3@tiXAxisFont    = 21
  res3@tiXAxisFontHeightF = 0.025
  res3@tmXBLabelsOn = False

  res3@tiMainOffsetYF = -0.01
  res3@tiMainOffsetXF = -0.13
  res3@vpWidthF          = 0.5
  res3@vpHeightF         = 0.5
  res3@tmYLLabelFontHeightF = 0.025 
  res3@tmXBLabelFontHeightF = 0.025

  res3@tmXBLabelFont= 21
  res3@tmXBMajorOutwardLengthF= 0.0
  res3@tmXBMinorOutwardLengthF= 0.0
  res3@tmXBMajorThicknessF= 2.0
  res3@tmYLLabelFont= 21
  res3@tmYLMajorOutwardLengthF= 0.0
  res3@tmYLMinorOutwardLengthF= 0.0
  res3@tmYLMajorThicknessF= 2.0

  res3@tmBorderThicknessF=  2.

  res3@trXMaxF     = 0.5
  res3@trXMinF     = -0.5

 ;res3@trYMaxF     = 200000.
 ;res3@trYMinF     = -200000.

  res3@trYMaxF     = 0.7
  res3@trYMinF     = 0.

  res3@gsnDraw         = False
  res3@gsnFrame         = False

  res3@xyLineColors   = (/"firebrick3","blue","orange","black"/)
  res3@xyDashPatterns = (/0,1,3,2/)
  res3@xyLineThicknesses =(/2.,2.,2.,2./)

  res3@tmXBMode = "explicit"
  res3@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  res3@tmXBLabels = (/-60,-30,0,30,60/)
  res3@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)
 
  plot(1)=gsn_csm_xy(wks,sin(lat_vec(:)/180.*PI),RMSE_vec(:,:),res3)

;-----------------------------------------------------------------------------------------------------------
  res4   = True
  res4@tiMainString = "f"
  res4@tiMainFont = 21
  res4@tiMainFontHeightF = 0.03
  res4@tiYAxisFont    = 21
  res4@tiYAxisFontHeightF = 0.025 
  res4@tmYLLabelsOn = False
  res4@tiXAxisString  = "Latitude"
  res4@tiXAxisFont    = 21
  res4@tiXAxisFontHeightF = 0.025

  res4@tiMainOffsetYF = -0.01
  res4@tiMainOffsetXF = -0.13
  res4@vpWidthF          = 0.5
  res4@vpHeightF         = 0.5
  res4@tmYLLabelFontHeightF = 0.025 
  res4@tmXBLabelFontHeightF = 0.025

  res4@tmXBLabelFont= 21
  res4@tmXBMajorOutwardLengthF= 0.0
  res4@tmXBMinorOutwardLengthF= 0.0
  res4@tmXBMajorThicknessF= 2.0
  res4@tmYLLabelFont= 21
  res4@tmYLMajorOutwardLengthF= 0.0
  res4@tmYLMinorOutwardLengthF= 0.0
  res4@tmYLMajorThicknessF= 2.0

  res4@tmBorderThicknessF=  2.

  res4@trXMaxF     = 0.5
  res4@trXMinF     = -0.5

 ;res4@trYMaxF     = 200000.
 ;res4@trYMinF     = -200000.

  res4@trYMaxF     = 1.
  res4@trYMinF     = 0.

  res4@gsnDraw         = False
  res4@gsnFrame         = False

  res4@xyLineColors   = (/"firebrick3","blue","orange","black","gray34"/)
  res4@xyDashPatterns = (/0,1,3,2,0/)
  res4@xyLineThicknesses =(/2.,2.,2.,2.,2./)

  res4@tmXBMode = "explicit"
  res4@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  res4@tmXBLabels = (/-60,-30,0,30,60/)
  res4@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)
 
  plot(5)=gsn_csm_xy(wks,sin(lat_vec(:)/180.*PI),correlation_array(:,:),res4)
  
;-----------------------------------------------------------------------------------------------------------
; draw panel with white space added
 resP                 = True
 resP@gsnPanelYWhiteSpacePercent = 0
 resP@gsnPanelXWhiteSpacePercent = 0
 resP@gsnPanelXF = (/0.15,0.46,0.15,0.46,0.15,0.46/)
 resP@gsnPanelYF = (/0.97,0.97,0.66,0.66,0.35,0.35/)
 ;resP@gsnPanelYF = (/-0.9,-0.9,-0.35,-0.35,-0.2,-0.2/)
 gsn_panel(wks,plot,(/3,2/),resP)





;-----------------------------------------------------------------------------------------------------------------

  wks2   = gsn_open_wks (wks_type,"fitting_func_9")
  plot2 = new(4,graphic)

  prec_func = new ((/5,res_number/),float)
  lat_func  = new ((/4,res_number/),float)

; index_vec(:,0) =   pr_W(:)
; index_vec(:,1) =   pr_E(:)
; index_vec(:,2) =   pr_A(:)
; index_vec(:,3) =   pr_S(:)
; index_vec(:,4) =   pr_I(:)
; index_vec(:,5) =   pr_D(:)
; index_vec(:,6) =   pr_M(:)
; index_vec(:,7) =   pr_L(:)
; index_vec(:,8) =   pr_R(:)


;print("A")



  prec_func(0,:) = par_vec(:,1)
  prec_func(1,:) = par_vec(:,2)
  prec_func(2,:) = par_vec(:,5)
  prec_func(3,:) = par_vec(:,6)
  prec_func(4,:) = par_vec(:,8)

  lat_func(0,:) = par_vec(:,0)
  lat_func(1,:) = par_vec(:,3)
  lat_func(2,:) = par_vec(:,4)
  lat_func(3,:) = par_vec(:,7)
  
  ;print(lat_func(0,:))

  resp   = True
  resp@tiMainString = "Zonal-mean P"
  resp@tiMainFont = 21
  resp@tiMainFontHeightF = 0.02
  resp@tiYAxisString  = "P [mm day~S~-1~N~]"
  resp@tiYAxisFont    = 21
  resp@tiYAxisFontHeightF = 0.02 
  resp@tiXAxisString  = "Latitude"
  resp@tiXAxisFont    = 21
  resp@tiXAxisFontHeightF = 0.02

  resp@tiMainOffsetYF = -0.01
  resp@tiMainOffsetXF = -0.07
  resp@vpWidthF          = 0.5
  resp@vpHeightF         = 0.5
  resp@tmYLLabelFontHeightF = 0.02 
  resp@tmXBLabelFontHeightF = 0.02

  resp@tmXBLabelFont= 21
  resp@tmXBMajorOutwardLengthF= 0.0
  resp@tmXBMinorOutwardLengthF= 0.0
  resp@tmXBMajorThicknessF= 2.0
  resp@tmYLLabelFont= 21
  resp@tmYLMajorOutwardLengthF= 0.0
  resp@tmYLMinorOutwardLengthF= 0.0
  resp@tmYLMajorThicknessF= 2.0

  resp@tmBorderThicknessF=  2.

  resp@trXMaxF     = 0.5
  resp@trXMinF     = -0.5

 resp@trYMaxF     = 2.
 resp@trYMinF     = -2.


  resp@gsnDraw         = False
  resp@gsnFrame         = False

  resp@xyLineColors   = (/"red","blue","grey34","purple4","yellow"/)
  resp@xyDashPatterns = (/1,1,1,1/)

  resp@tmXBMode = "explicit"
  resp@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  resp@tmXBLabels = (/-60,-30,0,30,60/)
  resp@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)


 ;print(sin(lat_vec(:)/180.*PI))

  plot2(0)=gsn_csm_xy(wks2,sin(lat_vec(:)/180.*PI),prec_func(:,:),resp)


;----------------------------------------------------------------------------------------------------
;print("B")

  resl   = True
  resl@tiMainString = "Zonal-mean P"
  resl@tiMainFont = 21
  resl@tiMainFontHeightF = 0.02
  resl@tiYAxisString  = "P [mm day~S~-1~N~]"
  resl@tiYAxisFont    = 21
  resl@tiYAxisFontHeightF = 0.02 
  resl@tiXAxisString  = "Latitude"
  resl@tiXAxisFont    = 21
  resl@tiXAxisFontHeightF = 0.02

  resl@tiMainOffsetYF = -0.01
  resl@tiMainOffsetXF = -0.07
  resl@vpWidthF          = 0.5
  resl@vpHeightF         = 0.5
  resl@tmYLLabelFontHeightF = 0.02 
  resl@tmXBLabelFontHeightF = 0.02

  resl@tmXBLabelFont= 21
  resl@tmXBMajorOutwardLengthF= 0.0
  resl@tmXBMinorOutwardLengthF= 0.0
  resl@tmXBMajorThicknessF= 2.0
  resl@tmYLLabelFont= 21
  resl@tmYLMajorOutwardLengthF= 0.0
  resl@tmYLMinorOutwardLengthF= 0.0
  resl@tmYLMajorThicknessF= 2.0

  resl@tmBorderThicknessF=  2.

  resl@trXMaxF     = 0.5
  resl@trXMinF     = -0.5

 resl@trYMaxF     = 2.
 resl@trYMinF     = -2.


  resl@gsnDraw         = False
  resl@gsnFrame         = False


  resl@xyLineColors   = (/"seagreen","black","orange","black"/)
  resl@xyDashPatterns = (/0,0,0,0/)

  resl@tmXBMode = "explicit"
  resl@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  resl@tmXBLabels = (/-60,-30,0,30,60/)
  resl@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)

  plot2(1)=gsn_csm_xy(wks2,sin(lat_vec(:)/180.*PI),lat_func(:,:),resl)
;-----------------------------------------------------------------------------------------------------
 ;print("C")
  
  resc   = True
  resc@tiMainString = "Zonal-mean P"
  resc@tiMainFont = 21
  resc@tiMainFontHeightF = 0.02
  resc@tiYAxisString  = "P [mm day~S~-1~N~]"
  resc@tiYAxisFont    = 21
  resc@tiYAxisFontHeightF = 0.02 
  resc@tiXAxisString  = "Latitude"
  resc@tiXAxisFont    = 21
  resc@tiXAxisFontHeightF = 0.02

  resc@tiMainOffsetYF = -0.01
  resc@tiMainOffsetXF = -0.07
  resc@vpWidthF          = 0.5
  resc@vpHeightF         = 0.5
  resc@tmYLLabelFontHeightF = 0.02 
  resc@tmXBLabelFontHeightF = 0.02

  resc@tmXBLabelFont= 21
  resc@tmXBMajorOutwardLengthF= 0.0
  resc@tmXBMinorOutwardLengthF= 0.0
  resc@tmXBMajorThicknessF= 2.0
  resc@tmYLLabelFont= 21
  resc@tmYLMajorOutwardLengthF= 0.0
  resc@tmYLMinorOutwardLengthF= 0.0
  resc@tmYLMajorThicknessF= 2.0

  resc@tmBorderThicknessF=  2.

  resc@trXMaxF     = 0.5
  resc@trXMinF     = -0.5

  resc@trYMaxF     = 6.
  resc@trYMinF     = -2.


  resc@gsnDraw         = False
  resc@gsnFrame         = False


  resc@xyLineColors   = (/"seagreen","red","blue","black","orange","grey34","purple4","black","yellow"/)
  resc@xyDashPatterns = (/0,1,1,0,0,1,1,0,1/)

  resc@tmXBMode = "explicit"
  resc@tmXBValues = (/sin(-60./180*PI),sin(-30./180*PI),sin(0.),sin(30./180*PI),sin(60./180*PI)/)
  resc@tmXBLabels = (/-60,-30,0,30,60/)
  resc@tmXBMinorValues =sin((/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)/180.*PI)

  ;print(res_number)

 

;print("Ca")
  plot2(2)=gsn_csm_xy(wks2,sin(lat_vec(14:(res_number-15))/180.*PI),contrib_vec(:,14:(res_number-15)),resc)
;print("D")


;----------------------------------------------------------------------------------------------------
; draw panel with white space added
 resP2                 = True
 resP2@gsnPanelYWhiteSpacePercent = 0
 resP2@gsnPanelXWhiteSpacePercent = 0

 resP2@gsnPanelXF = (/0.08,0.41,0.08,0.41,0.08,0.41/)
 ;resP2@gsnPanelYF = (/-0.6,-0.6,-0.35,-0.35,-0.1,-0.1/)

 gsn_panel(wks2,plot2,(/3,2/),resP2)




;print("E")



end