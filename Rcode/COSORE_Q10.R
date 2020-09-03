
## Day and night Q10 comparison based on COSORE
### prepare COSORE data
get_cosore_data <- function (){
  csr_database() %>% 
    dplyr::filter(grepl("Rh",CSR_MSMT_VAR)) %>% 
    dplyr::select(CSR_DATASET, CSR_LONGITUDE, CSR_LATITUDE, CSR_DATE_BEGIN, CSR_DATE_END) %>% 
    mutate(n_yr = year(CSR_DATE_END) - year(CSR_DATE_BEGIN)) ->
    csr_rh_dset
  
  csr_rh_dset$CSR_DATASET
  csr_dataset("d20190424_ZHANG_maple")$ports
  csr_dataset("d20190424_ZHANG_maple")$data %>% 
    filter(CSR_PORT %in% c(1,7)) %>% 
    mutate(depth = 5, dset = "dset1") %>% 
    rename('CSR_SM' = 'CSR_SM30',
           'CSR_TS' = 'CSR_T5') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set1
  
  csr_dataset("d20190424_ZHANG_oak")$ports
  csr_dataset("d20190424_ZHANG_oak")$data %>% 
    filter(CSR_PORT %in% c(1,2,7,8)) %>% 
    mutate(depth = 5, dset = "dset2") %>% 
    rename('CSR_SM' = 'CSR_SM30',
           'CSR_TS' = 'CSR_T5') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set2
  
  csr_dataset("d20190430_DESAI")$ports
  csr_dataset("d20190430_DESAI")$data %>% 
    filter(CSR_PORT %in% c(4)) ->
    d_set3
  colnames(d_set3)
  d_set3 %>% 
    mutate(CSR_TS = coalesce(CSR_T7, CSR_T8, CSR_T9, CSR_T14, CSR_T11, CSR_T15, CSR_T16.5, CSR_T20, CSR_T20.5),
           CSR_SM = coalesce(CSR_SM4, CSR_SM18)) %>% 
    mutate(depth = 14, dset = "dset3") %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set3
  
  csr_dataset("d20190504_SAVAGE_hf006-05")$ports
  csr_dataset("d20190504_SAVAGE_hf006-05")$data %>% 
    filter(CSR_PORT %in% c(1,3,5,7)) %>% 
    mutate(depth = 10, dset = "dset4") %>% 
    rename('CSR_SM' = 'CSR_SM10',
           'CSR_TS' = 'CSR_T10') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set4
  
  csr_dataset("d20190517_MAURITZ")$ports
  csr_dataset("d20190517_MAURITZ")$data %>% 
    filter(CSR_PORT %in% c(1,3,5,7,9)) %>% 
    mutate(depth = 5, dset = "dset5") %>% 
    rename('CSR_SM' = 'CSR_SM5',
           'CSR_TS' = 'CSR_T5') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set5
  
  csr_dataset("d20190610_SIHI_H2")$ports
  csr_dataset("d20190610_SIHI_H2")$data %>% 
    filter(CSR_PORT %in% c(20:22)) %>% 
    mutate(depth = 10, dset = "dset6") %>% 
    rename('CSR_SM' = 'CSR_SM10',
           'CSR_TS' = 'CSR_T10') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set6
  
  csr_dataset("d20190617_SCOTT_WKG")$ports
  csr_dataset("d20190617_SCOTT_WKG")$data %>% 
    filter(CSR_PORT %in% c(5:7)) %>% 
    mutate(depth = 5, dset = "dset7") %>% 
    rename('CSR_SM' = 'CSR_SM5',
           'CSR_TS' = 'CSR_T5') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set7
  
  csr_dataset("d20190830_LIANG")$ports
  csr_dataset("d20190830_LIANG")$data %>% 
    filter(CSR_PORT %in% c(1,3,4,5,8,9,10,11,13,15)) %>% 
    rename('CSR_TS' = 'CSR_T5') %>% 
    mutate(depth = 5, dset = "dset8", CSR_SM = NA) %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set8
  
  csr_dataset("d20200109_HIRANO_PDB")$ports
  csr_dataset("d20200109_HIRANO_PDB")$data %>% 
    filter(CSR_PORT %in% c(1)) %>% 
    mutate(depth = 5, dset = "dset9") %>% 
    rename('CSR_SM' = 'CSR_SM30',
           'CSR_TS' = 'CSR_T5') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set9
  
  csr_dataset("d20200122_BLACK")$ports
  csr_dataset("d20200122_BLACK")$data %>% 
    filter(CSR_PORT %in% c(3,9,12)) %>% 
    mutate(depth = 2, dset = "dset10") %>% 
    rename('CSR_SM' = 'CSR_SM7.5',
           'CSR_TS' = 'CSR_T2') %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set10
  d_set10 %>% 
    filter(is.na(CSR_SM))
  
  csr_dataset("d20200331_PEICHL")$ports
  csr_dataset("d20200331_PEICHL")$data %>% 
    filter(CSR_PORT %in% c(3,6,9,12)) %>% 
    rename('CSR_TS' = 'CSR_T10') %>% 
    mutate(depth = 10, dset = "dset11", CSR_SM = NA) %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set11
  
  csr_dataset("d20200417_ARAIN_TP39")$ports
  csr_dataset("d20200417_ARAIN_TP39")$data %>% 
    filter(CSR_PORT %in% c(3,8)) %>% 
    rename('CSR_TS' = 'CSR_T5',
           'CSR_SM' = 'CSR_SM5') %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM) %>%
    mutate(depth = 5, dset = "dset12") %>% 
    filter(!is.na(CSR_TS)) %>% 
    dplyr::select(CSR_PORT, CSR_TIMESTAMP_BEGIN, CSR_TIMESTAMP_END, CSR_FLUX_CO2, CSR_TS, CSR_SM, depth, dset) ->
    d_set12
  
  bind_rows(
    d_set1, d_set2, d_set3, d_set4, d_set5, d_set6, 
    d_set7, d_set8, d_set9, d_set10, d_set11, d_set12) ->
    d_set_all
  
  # Prepare time information and ID
  d_set_all %>% 
    mutate(YYYMMDDD = as.Date(CSR_TIMESTAMP_END, "%Y-%m"),
           end_hour = hour(CSR_TIMESTAMP_END),
           DN_label = case_when(
             end_hour %in% c(7:18) ~ "Day", #need to discuss when is day time
             TRUE ~ "Night")) %>% 
    filter(CSR_FLUX_CO2 > 0) %>% 
    mutate(Rh_log = log(CSR_FLUX_CO2),
           ID = paste0(dset, "-", CSR_PORT, "-", year(CSR_TIMESTAMP_END),"-", month(CSR_TIMESTAMP_END), "-", week(CSR_TIMESTAMP_END))) ->
    d_set_all
  
  return(d_set_all)
}

# calculate Q10 for cosore data
get_cosore_q10 <- function (sdata) {
  sdata %>% 
    dplyr::select(ID) %>% 
    unique() ->
    d_set_id
  
  out <- data.frame() # create a data frame to hold the results
  for(i in 1:nrow(d_set_id)){
    sdata %>% 
      filter(ID == d_set_id$ID[i] & DN_label == "Day") ->
      sub_day
    
    sdata %>% 
      filter(ID == d_set_id$ID[i] & DN_label == "Night") ->
      sub_night
    
    n_day <- nrow(sub_day)
    n_night <- nrow(sub_night)
    
    m_day <- try(lm(Rh_log ~ CSR_TS, data = sub_day))
    m_night <- try(lm(Rh_log ~ CSR_TS, data = sub_night))
    
    R2_day <- ifelse(is(m_day, "try-error"), NA, summary(m_day)$r.squared)
    R2_night <- ifelse(is(m_night, "try-error"), NA, summary(m_night)$r.squared)
    
    # calculate Q10
    Q10_day <- ifelse(nrow(sub_day) < 4, NA, csr_rh_Q10(sub_day)) # only calculate Q10 when have more than 6 observations
    Q10_night <- ifelse(nrow(sub_night) < 4, NA, csr_rh_Q10(sub_night)) # only calculate Q10 when have more than 6 observations
    
    # calculate mean ST and SM
    m_ts_day <- mean(sub_day$CSR_TS, na.rm = T)
    m_sm_day <- mean(sub_day$CSR_SM, na.rm = T)
    m_ts_night <- mean(sub_night$CSR_TS, na.rm = T)
    m_sm_night <- mean(sub_night$CSR_SM, na.rm = T)
    
    print(paste0("*****", i))
    
    # combine results
    out <- rbind(out,
                 data.frame(i, d_set_id$ID[i], Q10_day, Q10_night, m_ts_day, m_sm_day, m_ts_night, m_sm_night,
                            n_day, n_night, R2_day, R2_night))
    
  }
  return(out)
}

