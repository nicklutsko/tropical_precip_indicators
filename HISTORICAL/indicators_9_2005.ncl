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

  pr_H        = new(m_number,float)
  pr_W        = new(m_number,float)
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

  if (precip(0,N_minimum,0) .le. precip(0,S_minimum,0)) then
    ;pr_L(j_f-1) = abs(doubletofloat(lat(N_minimum)))
    pr_M(j_f-1) = pr_E(j_f-1) - precip(0,N_minimum,0)
    
  else   
    ;pr_L(j_f-1) = abs(doubletofloat(lat(S_minimum)))
    pr_M(j_f-1) = pr_E(j_f-1) - precip(0,S_minimum,0)
  end if

  pr_R(j_f-1) = precip(0,N_minimum,0) - precip(0,S_minimum,0)

  if (lat(ITCZ) .GT. 0.) then
    pr_H(j_f-1) = precip(0,ITCZ,0) - precip(0,N_minimum,0)
    pr_L(j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(N_minimum))) 
  else
    pr_H(j_f-1) = precip(0,ITCZ,0) - precip(0,S_minimum,0) 
    pr_L(j_f-1) = abs(doubletofloat(lat(ITCZ) - lat(S_minimum)))
  end if

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


  print(pr_L(:))
  print(pr_W(:))
  

   value_matrixHW = new((/2,m_number/),float)
   value_matrixHW(0,:) = pr_H(:)
   value_matrixHW(1,:) = pr_W(:)

   value_matrixHE = new((/2,m_number/),float)
   value_matrixHE(0,:) = pr_H(:)
   value_matrixHE(1,:) = pr_E(:)

   value_matrixHA = new((/2,m_number/),float)
   value_matrixHA(0,:) = pr_H(:)
   value_matrixHA(1,:) = pr_A(:)

   value_matrixHS = new((/2,m_number/),float)
   value_matrixHS(0,:) = pr_H(:)
   value_matrixHS(1,:) = pr_S(:)   

   value_matrixHI = new((/2,m_number/),float)
   value_matrixHI(0,:) = pr_H(:)
   value_matrixHI(1,:) = pr_I(:)

   value_matrixHD = new((/2,m_number/),float)
   value_matrixHD(0,:) = pr_H(:)
   value_matrixHD(1,:) = pr_D(:)

   value_matrixHL = new((/2,m_number/),float)
   value_matrixHL(0,:) = pr_H(:)
   value_matrixHL(1,:) = pr_L(:)

   value_matrixHR = new((/2,m_number/),float)
   value_matrixHR(0,:) = pr_H(:)
   value_matrixHR(1,:) = pr_R(:)   


   value_matrixWE = new((/2,m_number/),float)
   value_matrixWE(0,:) = pr_W(:)
   value_matrixWE(1,:) = pr_E(:)

   value_matrixWA = new((/2,m_number/),float)
   value_matrixWA(0,:) = pr_W(:)
   value_matrixWA(1,:) = pr_A(:)

   value_matrixWS = new((/2,m_number/),float)
   value_matrixWS(0,:) = pr_W(:)
   value_matrixWS(1,:) = pr_S(:)   

   value_matrixWI = new((/2,m_number/),float)
   value_matrixWI(0,:) = pr_W(:)
   value_matrixWI(1,:) = pr_I(:)

   value_matrixWD = new((/2,m_number/),float)
   value_matrixWD(0,:) = pr_W(:)
   value_matrixWD(1,:) = pr_D(:)

   value_matrixWL = new((/2,m_number/),float)
   value_matrixWL(0,:) = pr_W(:)
   value_matrixWL(1,:) = pr_L(:)

   value_matrixWR = new((/2,m_number/),float)
   value_matrixWR(0,:) = pr_W(:)
   value_matrixWR(1,:) = pr_R(:)   
   

   value_matrixEA = new((/2,m_number/),float)
   value_matrixEA(0,:) = pr_E(:)
   value_matrixEA(1,:) = pr_A(:)

   value_matrixES = new((/2,m_number/),float)
   value_matrixES(0,:) = pr_E(:)
   value_matrixES(1,:) = pr_S(:)   

   value_matrixEI = new((/2,m_number/),float)
   value_matrixEI(0,:) = pr_E(:)
   value_matrixEI(1,:) = pr_I(:)

   value_matrixED = new((/2,m_number/),float)
   value_matrixED(0,:) = pr_E(:)
   value_matrixED(1,:) = pr_D(:)

   value_matrixEL = new((/2,m_number/),float)
   value_matrixEL(0,:) = pr_E(:)
   value_matrixEL(1,:) = pr_L(:)

   value_matrixER = new((/2,m_number/),float)
   value_matrixER(0,:) = pr_E(:)
   value_matrixER(1,:) = pr_R(:)   


   value_matrixAS = new((/2,m_number/),float)
   value_matrixAS(0,:) = pr_A(:)
   value_matrixAS(1,:) = pr_S(:)   

   value_matrixAI = new((/2,m_number/),float)
   value_matrixAI(0,:) = pr_A(:)
   value_matrixAI(1,:) = pr_I(:)

   value_matrixAD = new((/2,m_number/),float)
   value_matrixAD(0,:) = pr_A(:)
   value_matrixAD(1,:) = pr_D(:)

   value_matrixAL = new((/2,m_number/),float)
   value_matrixAL(0,:) = pr_A(:)
   value_matrixAL(1,:) = pr_L(:)

   value_matrixAR = new((/2,m_number/),float)
   value_matrixAR(0,:) = pr_A(:)
   value_matrixAR(1,:) = pr_R(:)   


   value_matrixSI = new((/2,m_number/),float)
   value_matrixSI(0,:) = pr_S(:)
   value_matrixSI(1,:) = pr_I(:)

   value_matrixSD = new((/2,m_number/),float)
   value_matrixSD(0,:) = pr_S(:)
   value_matrixSD(1,:) = pr_D(:)

   value_matrixSL = new((/2,m_number/),float)
   value_matrixSL(0,:) = pr_S(:)
   value_matrixSL(1,:) = pr_L(:)

   value_matrixSR = new((/2,m_number/),float)
   value_matrixSR(0,:) = pr_S(:)
   value_matrixSR(1,:) = pr_R(:)   


   value_matrixID = new((/2,m_number/),float)
   value_matrixID(0,:) = pr_I(:)
   value_matrixID(1,:) = pr_D(:)  

   value_matrixIL = new((/2,m_number/),float)
   value_matrixIL(0,:) = pr_I(:)
   value_matrixIL(1,:) = pr_L(:)

   value_matrixIR = new((/2,m_number/),float)
   value_matrixIR(0,:) = pr_I(:)
   value_matrixIR(1,:) = pr_R(:)


   value_matrixDL = new((/2,m_number/),float)
   value_matrixDL(0,:) = pr_D(:)
   value_matrixDL(1,:) = pr_L(:)

   value_matrixDR = new((/2,m_number/),float)
   value_matrixDR(0,:) = pr_D(:)
   value_matrixDR(1,:) = pr_R(:)


   value_matrixLR = new((/2,m_number/),float)
   value_matrixLR(0,:) = pr_L(:)
   value_matrixLR(1,:) = pr_R(:)
   
    cor_matrix_HW = covcorm(value_matrixHW,(/1,0/)) 
    cor_matrix_HE = covcorm(value_matrixHE,(/1,0/))
    cor_matrix_HA = covcorm(value_matrixHA,(/1,0/))
    cor_matrix_HS = covcorm(value_matrixHS,(/1,0/))
    cor_matrix_HI = covcorm(value_matrixHI,(/1,0/))
    cor_matrix_HD = covcorm(value_matrixHD,(/1,0/))
    cor_matrix_HL = covcorm(value_matrixHL,(/1,0/))
    cor_matrix_HR = covcorm(value_matrixHR,(/1,0/))

    cor_matrix_WE = covcorm(value_matrixWE,(/1,0/))
    cor_matrix_WA = covcorm(value_matrixWA,(/1,0/))
    cor_matrix_WS = covcorm(value_matrixWS,(/1,0/))
    cor_matrix_WI = covcorm(value_matrixWI,(/1,0/))
    cor_matrix_WD = covcorm(value_matrixWD,(/1,0/))
    cor_matrix_WL = covcorm(value_matrixWL,(/1,0/))
    cor_matrix_WR = covcorm(value_matrixWR,(/1,0/))

    cor_matrix_EA = covcorm(value_matrixEA,(/1,0/))
    cor_matrix_ES = covcorm(value_matrixES,(/1,0/))
    cor_matrix_EI = covcorm(value_matrixEI,(/1,0/))
    cor_matrix_ED = covcorm(value_matrixED,(/1,0/))
    cor_matrix_EL = covcorm(value_matrixEL,(/1,0/))
    cor_matrix_ER = covcorm(value_matrixER,(/1,0/))

    cor_matrix_AS = covcorm(value_matrixAS,(/1,0/))
    cor_matrix_AI = covcorm(value_matrixAI,(/1,0/))
    cor_matrix_AD = covcorm(value_matrixAD,(/1,0/))
    cor_matrix_AL = covcorm(value_matrixAL,(/1,0/))
    cor_matrix_AR = covcorm(value_matrixAR,(/1,0/))

    cor_matrix_SI = covcorm(value_matrixSI,(/1,0/))
    cor_matrix_SD = covcorm(value_matrixSD,(/1,0/))
    cor_matrix_SL = covcorm(value_matrixSL,(/1,0/))
    cor_matrix_SR = covcorm(value_matrixSR,(/1,0/))

    cor_matrix_ID = covcorm(value_matrixID,(/1,0/))
    cor_matrix_IL = covcorm(value_matrixIL,(/1,0/))
    cor_matrix_IR = covcorm(value_matrixIR,(/1,0/))

    cor_matrix_DL = covcorm(value_matrixDL,(/1,0/))
    cor_matrix_DR = covcorm(value_matrixDR,(/1,0/))

    cor_matrix_LR = covcorm(value_matrixLR,(/1,0/))
 

    print("HW:")
    print(cor_matrix_HW(1)^2)

    print("HE:")
    print(cor_matrix_HE(1)^2)

    print("HA:")
    print(cor_matrix_HA(1)^2)

    print("HS:")
    print(cor_matrix_HS(1)^2)

    print("HI:")
    print(cor_matrix_HI(1)^2)

    print("HD:")
    print(cor_matrix_HD(1)^2)

    print("HL:")
    print(cor_matrix_HL(1)^2)

    print("HR:")
    print(cor_matrix_HR(1)^2)

    print("   -------   ")

    print("WE:")
    print(cor_matrix_WE(1)^2)

    print("WA:")
    print(cor_matrix_WA(1)^2)

    print("WS:")
    print(cor_matrix_WS(1)^2)

    print("WI:")
    print(cor_matrix_WI(1)^2)

    print("WD:")
    print(cor_matrix_WD(1)^2)

    print("WL:")
    print(cor_matrix_WL(1)^2)

    print("WR:")
    print(cor_matrix_WR(1)^2)

    print("   -------   ")

    print("EA:")
    print(cor_matrix_EA(1)^2)

    print("ES:")
    print(cor_matrix_ES(1)^2)

    print("EI:")
    print(cor_matrix_EI(1)^2)

    print("ED:")
    print(cor_matrix_ED(1)^2)

    print("EL:")
    print(cor_matrix_EL(1)^2)

    print("ER:")
    print(cor_matrix_ER(1)^2)

    print("   -------   ")

    print("AS:")
    print(cor_matrix_AS(1)^2)

    print("AI:")
    print(cor_matrix_AI(1)^2)

    print("AD:")
    print(cor_matrix_AD(1)^2)

    print("AL:")
    print(cor_matrix_AL(1)^2)

    print("AR:")
    print(cor_matrix_AR(1)^2)

    print("   -------   ")

    print("SI:")
    print(cor_matrix_SI(1)^2)

    print("SD:")
    print(cor_matrix_SD(1)^2)

    print("SL:")
    print(cor_matrix_SL(1)^2)

    print("SR:")
    print(cor_matrix_SR(1)^2)

    print("   -------   ")

    print("ID:")
    print(cor_matrix_ID(1)^2)

    print("IL:")
    print(cor_matrix_IL(1)^2)

    print("IR:")
    print(cor_matrix_IR(1)^2)

    print("   -------   ")

    print("DL:")
    print(cor_matrix_DL(1)^2)

    print("DR:")
    print(cor_matrix_DR(1)^2)

    print("   -------   ")
    
    print("LR:")
    print(cor_matrix_LR(1)^2)  

    
 

   

end