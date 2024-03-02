#MANUSCRIPT: Selective social tolerance drives differentiated relationships among wild female chimpanzees
#AUTH: Stephanie Fox, Martin N. Muller, Natalia Camargo Peña, Nicole Thompson González, Zarin Machanda, Emily Otali, Richard Wrangham, Melissa Emery Thompson
# corresponding author: safox@ucsb.edu 
# Department of Anthropology
# University of California Santa Barbara


# Table of Contents ####
# DATASETS
# 1. Summary of social diff from obs and null data
# 2. Summary of PPI from obs and null data 
# 3. Dyadic PPI for logistic regression
# 4. Targeted behavioural data: 5m index, joint movement, hinde index, cofeeding, grooming
# 5. Long term behavioural data: 5m index, aggression, grooming
# 6. Vigilance dataset
# 7. Scratching dataset

# ANALYSIS 
# 1. Preparation of permuted social networks
# 2. Calculation of social differentiation (obs and null)
# 3. Calculation of dyadic PPI values (obs and null)
# 4. Stability analysis 1: PPI vs null distribution
# 5. Stability analysis 2: PPI logistic regression
# 6. Example joint movements model
# 7. Example mutual association model
# 8. Example aggression model
# 9. Example cofeeding model
# 10. Example vigilance model
# 11. Example scratching model
# 12. Example model checks
# 13. Example of permuted models and summarize results
# 14. Example model comparison

# NOTES ####
# - Instead of including the 1000 permuted datasets used to form the null distributions, 
# I am including the code on how I made them. I am not including the raw data for constructing them. 
# - For models where I permuted the model 1000 times while switching ID1 and ID2, 
# example models are provided below, and then below that is example code for how I 
# permuted them and summarized results 



#:####
# GET DATA ####
load("publicdata_fc_tolerance.Rdata")
#rename these datasets so that they suit all the code below
# 1. Summary of social diff from obs and null data
sdiff_all1yr <- public_sdiff_all1yr
# 2. Summary of PPI from obs and null data 
ppi_summary_top2 <- public_ppi_summary_top2 
# 3. Dyadic PPI for logistic regression
pp_log_2.data <- public_pp_log_2.data 
# 4. Targeted behavioural data: 5m index, joint movement, hinde index, cofeeding, grooming
behav_undirected <- public_behav_undirected 
# 5. Long term behavioural data: 5m index, aggression, grooming
dyad_w5_agg2 <- public_dyad_w5_agg2 
# 6. Vigilance dataset
AVS_logdata <- public_AVS_logdata 
# 7. Scratching dataset
CH_logdata <- public_CH_logdata 


#:####
# ANALYSIS ####
# 1. Preparation of permuted social networks ####
# NOTE: this code is to show HOW I did this, but I have not provided the raw 5m prox scans

## method:  
# make a list of periods: 
# filter to period 
# peel off the w5 
# shuffle it
# reattach it
# make a new column "fix" for where if "focal = w5" = "Y"
# make a new column for where group by --> more then 2 of same individual 
# use the new column to identify rows to fix
# for loop: 
# go through every line
# if "fix" column = y then
# make object called "l2check" = FALSE
# while(isFALSE(l2check))
# ID1 = the id from the w5 column 
# ID2 = the id from the focal column 
# [j] <- randomly choose another line number
# ID3 = the id from the random w5 column 
# ID4 = the id from the random focal column 
# check if ID1 != ID4 AND ID2 != ID3, then:
# new df <- filter by the date time and the focal ID and the new ID (scan A)
# new df <- filter by the date time and the focal ID and the new ID (scan B)
# check #if row length is 0 for new df A and new df B. if true, then:
# check if the date of the scan is within adult hood dates for swapped ID, then:
#set l2check = TRUE
# end of the while loop
# set the new IDs w5[i] = ID# and w5[j] = ID#
# end of if statement
# end of for loop 
# end of period 
# save this dataset to the slot in the list 

list.f5m1yr <- vector("list", length = 1000)
list.period <- vector("list", length = 10)
#pb <- progress_bar$new(format = " running [:bar] :percent eta: :eta", total = 10, clear = FALSE, width= 10)
set.seed(123) #must set seed outside of the loop so that it doesn't reset each time the loop runs 
p <- c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")
t <- Sys.time()
for(i in seq(1000)){
  #pb$tick()
  #Sys.sleep(1/10)
  for(k in seq(10)){
    per <- p[k]
    f5m_period = f5m1yr %>% filter(period == per)
    w5data <- f5m_period %>% select(w5)
    w5shuffle <- sample(w5data$w5, replace = FALSE) %>% data.frame() %>% rename(w5 = ".")
    rand1 <- f5m_period %>% # add the shuffled w5 back on 
      select(-w5) %>%
      bind_cols(w5shuffle)%>%
      left_join(adulthood %>% select(chimp_id, adulthood_start, adulthood_end), by = c(w5 = "chimp_id")) %>%
      mutate(fix = ifelse(Focal == w5, "y", "n")) %>% # if focal = w5 then this line needs to be fixed
      mutate(fix = ifelse(Date_Time < adulthood_start | Date_Time > adulthood_end, "y", fix))
    doubles <- rand1 %>% #find cases where same w5 occurs twice in the same scan
      group_by(Date_Time, w5) %>% 
      summarise(n = n()) %>% 
      filter(n >1) %>%
      mutate(fix2 = "y") #tag that these need fixing
    rand2 <- rand1 %>%
      left_join(doubles %>% select(Date_Time, w5, fix2), by = c("Date_Time", "w5")) %>% #add the fix2 tags onto the main dataframe
      mutate(fix = ifelse(fix2 == "y" & !is.na(fix2), "y", fix)) %>%
      select(-fix2) 
    #fix the problem rows by randomly choosing a new row, checking criteria, then flipping the w5 IDs if its ok
    for(j in 1:length(rand2$Focal)){
      if(rand2$fix[j] == "y"){
        IDA <- rand2$Focal[j]
        IDB <- rand2$w5[j]
        ticker <- 0 
        l2check = FALSE
        while(isFALSE(l2check)){
          ticker <- ticker + 1
          l2 <- sample(nrow(rand2), 1) #choose the new row with which to try flipping
          IDa <- rand2$Focal[l2]
          IDb <- rand2$w5[l2]
          #check if IDs are ok
          if(IDA != IDb & IDa != IDB){ #will not create a focal = w5 problem
            #check if there is more than one line per scan
            dt1 <- rand2$Date_Time[j]
            dt2 <- rand2$Date_Time[l2]
            #filter to minidataset of any lines where the new ID is already in the scan
            dt_rows1 <- rand2 %>% filter(Date_Time == dt1 & Focal == IDA & w5 == IDb) %>% nrow()
            dt_rows2 <- rand2 %>% filter(Date_Time == dt2 & Focal == IDa & w5 == IDB) %>% nrow()
            if(dt_rows1 == 0){ #if there is no rows then no other of the same ID in that scan already
              if(dt_rows2 == 0){ #same with the other ID and other scan
                #check if the new female ID was an adult at the time
                if(rand2$adulthood_start[l2] <= rand2$Date[j] & rand2$adulthood_end[l2] >= rand2$Date[j]){
                  if(rand2$adulthood_start[j] <= rand2$Date[l2] & rand2$adulthood_end[j] >= rand2$Date[l2]){
                    l2check = TRUE
                    rand2$w5[j] <- IDb
                    rand2$w5[l2] <- IDB
                  }}}}}}}}
    
    list.period[[k]] <- rand2
  }
  rand_all <- do.call("rbind", list.period)  
  list.f5m1yr[[i]] <- rand_all
}
Sys.time() - t
# 3 min per two runs 


save(list.f5m1yr, file = "listf5m1yr.Rdata")


# 2. Calculation of social differentiation (obs and null)####

# social_differentiation function is from MNWeiss/aninet (https://github.com/MNWeiss/aninet/blob/master/R/social_differentiation.R)
#need to get w5 count (numerator) and scan count (denominator) into vectors
t <- Sys.time()
list.sdiff.real <- vector("list", length = 5)
p <- c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")
for(i in seq(10)){
  per <- p[i]
  num <- w5_rates1yr %>% filter(period == per) %>% pull(w5_count)
  den <- w5_rates1yr %>% filter(period == per) %>% pull(scans_dyad) 
  sdiff <- social_differentiation(num, den, initial.params = c(0.1, 0.1), nsim = 10000) 
  list.sdiff.real[i] <- sdiff[2:2]
}
sdiff_real1yr <- do.call("rbind", list.sdiff.real) %>% 
  data.frame() %>%
  rename(sdiff_period = ".") %>%
  mutate(period = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")) %>%
  select(period, sdiff_period)
Sys.time() - t

mean.sdiff.real <- mean(sdiff_real1yr$sdiff_period)
sd.sdiff.real <- sd(sdiff_real1yr$sdiff_period)

#load in the sdiff_summary1yr which is the 1000 permuted S values and sd turned to means and CI
sdiff_all1yr <- sdiff_real1yr %>%
  left_join(sdiff_summary1yr)
#graph observed vs null 
graph.socdiff <- ggplot(sdiff_all1yr, aes(x = period, y = sdiff_period, group = 1)) + 
  geom_point() +
  #geom_line(col='red') +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.3) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), alpha = 0.4) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_blank(), 
        text = element_text(size = 10)) +
  xlab("Period") + ylab("Estimate of Social \n Differentiation (S)") 
graph.socdiff

# Turn into permuted w5 counts
list.w5rand1yr <- vector("list", length = 1000)
#make a blank dataframe each loop will fill in
w5_blanks <- w5_rates1yr %>% 
  mutate(w5_count = NA) %>% 
  select(-p_start, -p_end, -scans_ID1, -scans_ID2, -w5_rate)
t <- Sys.time()
for(i in seq(1000)){
  dat <- list.f5m1yr[[i]] %>% data.frame
  randcounts <- w5_blanks 
  #count all of the w5 for each dyad
  for(j in 1:length(undir_dyads1yr$ID1)){
    IDA <- undir_dyads1yr$ID1[j] 
    IDB <- undir_dyads1yr$ID2[j]
    per <- undir_dyads1yr$period[j]
    count <- dat %>% filter(Focal == IDA & w5 == IDB & period == per |
                              Focal == IDB & w5 == IDA & period == per) %>%
      tally()
    randcounts$w5_count[j] = count
  }
  randcounts <- randcounts %>% data.frame %>% 
    mutate(w5_count = as.numeric(w5_count), 
           w5rate = w5_count/scans_dyad)
  list.w5rand1yr[[i]] = randcounts
}
Sys.time() - t 
#3 min per two runs
#save(list.w5rand1yr, file = "listw5rand1yr.Rdata")

# Null network differentiation scores
list.socdiff.final.1yr <- vector("list", length = 1000)
p <- c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")
t <- Sys.time()
for(i in seq(1000)){
  randcounts <- list.w5rand1yr[[i]] %>% data.frame
  #take the randcounts dataset now and use it to calc soc diff per period
  #each socdiff will be stored in a list, which we compile into a list of 5 soc.diff (one for each period)
  randcounts <- randcounts %>% unnest(cols = c(w5_count))
  list.socdiff.rand <- vector("list", length = 10)
  for(k in seq(10)){
    per <- p[k]
    num <- randcounts %>% filter(period == per) %>% pull(w5_count)
    #the denom is still coming from the total scans per female, which is counted above
    den <- w5_rates1yr %>% filter(period == per) %>% pull(scans_dyad) 
    sdiff <- social_differentiation(num, den, initial.params = c(0.1, 0.1), nsim = 1) 
    list.socdiff.rand[k] <- sdiff[2:2]
  }
  #take the mini list of
  socdiff.rand <- do.call("rbind", list.socdiff.rand) %>% 
    data.frame() %>%
    rename(sdiff_period = ".") %>%
    mutate(period = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")) %>%
    select(period, sdiff_period)
  list.socdiff.final.1yr[[i]] <- socdiff.rand
}
Sys.time() - t
# takes about 5 minutes

sdiff_perm1yr <- do.call("rbind", list.socdiff.final.1yr) %>% data.frame
save(sdiff_perm.1yr, file = "sdiff_perm.Rdata")

#check distributions and get 95%CI per period
sd2010 <- sdiff_perm1yr %>% filter(period == "2010")
hist(sd2010$sdiff_period)
sd2010.95 <- get_confidence_interval(sd2010, level = .95) #this is from package infer
sd2010.80 <- get_confidence_interval(sd2010, level = .80)

sd2011 <- sdiff_perm1yr %>% filter(period == "2011")
hist(sd2011$sdiff_period)
sd2011.95 <- get_confidence_interval(sd2011, level = .95) #this is from package infer
sd2011.80 <- get_confidence_interval(sd2011, level = .80)

sd2012 <- sdiff_perm1yr %>% filter(period == "2012")
hist(sd2012$sdiff_period)
sd2012.95 <- get_confidence_interval(sd2012, level = .95) #this is from package infer
sd2012.80 <- get_confidence_interval(sd2012, level = .80)

sd2013 <- sdiff_perm1yr %>% filter(period == "2013")
hist(sd2013$sdiff_period)
sd2013.95 <- get_confidence_interval(sd2013, level = .95) #this is from package infer
sd2013.80 <- get_confidence_interval(sd2013, level = .80)

sd2014 <- sdiff_perm1yr %>% filter(period == "2014")
hist(sd2014$sdiff_period)
sd2014.95 <- get_confidence_interval(sd2014, level = .95) #this is from package infer
sd2014.80 <- get_confidence_interval(sd2014, level = .80)

sd2015 <- sdiff_perm1yr %>% filter(period == "2015")
hist(sd2015$sdiff_period)
sd2015.95 <- get_confidence_interval(sd2015, level = .95) #this is from package infer
sd2015.80 <- get_confidence_interval(sd2015, level = .80)

sd2016 <- sdiff_perm1yr %>% filter(period == "2016")
hist(sd2016$sdiff_period)
sd2016.95 <- get_confidence_interval(sd2016, level = .95) #this is from package infer
sd2016.80 <- get_confidence_interval(sd2016, level = .80)

sd2017 <- sdiff_perm1yr %>% filter(period == "2017")
hist(sd2017$sdiff_period)
sd2017.95 <- get_confidence_interval(sd2017, level = .95) #this is from package infer
sd2017.80 <- get_confidence_interval(sd2017, level = .80)

sd2018 <- sdiff_perm1yr %>% filter(period == "2018")
hist(sd2018$sdiff_period)
sd2018.95 <- get_confidence_interval(sd2018, level = .95) #this is from package infer
sd2018.80 <- get_confidence_interval(sd2018, level = .80)

sd2019 <- sdiff_perm1yr %>% filter(period == "2019")
hist(sd2019$sdiff_period)
sd2019.95 <- get_confidence_interval(sd2019, level = .95) #this is from package infer
sd2019.80 <- get_confidence_interval(sd2019, level = .80)


CI95 <- sd2010.95 %>%
  rbind(sd2011.95, sd2012.95, sd2013.95, sd2014.95, sd2015.95, sd2016.95, sd2017.95, sd2018.95, sd2019.95) %>%
  mutate(period = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")) %>%
  rename(lower95 = lower_ci, 
         upper95 = upper_ci)

CI80 <- sd2010.80 %>%
  rbind(sd2011.80, sd2012.80, sd2013.80, sd2014.80, sd2015.80, sd2016.80, sd2017.80, sd2018.80, sd2019.80) %>%
  mutate(period = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")) %>%
  rename(lower80 = lower_ci, 
         upper80 = upper_ci)

sdiff_summary1yr <-sdiff_perm1yr %>% 
  group_by(period) %>%
  summarise(mean.sdiff = mean(sdiff_period), 
            sd.sdiff = sd(sdiff_period), 
            n.sdiff = n(), 
            max.sdiff = max(sdiff_period), 
            min.sdiff = min(sdiff_period), 
            se.sdiff = sd.sdiff/sqrt(n.sdiff)) %>%
  left_join(CI80) %>%
  left_join(CI95)

save(sdiff_summary1yr, file = "sdiff_summary1yr.Rdata")

#make a dataframe with the sdiff values and the perm confidence set
sdiff_all1yr <- sdiff_real1yr %>%
  left_join(sdiff_summary1yr)

graph.socdiff <- ggplot(sdiff_all1yr, aes(x = period, y = sdiff_period, group = 1)) + 
  geom_point() +
  #geom_line(col='red') +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.3) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), alpha = 0.4) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_blank(), 
        text = element_text(size = 10)) +
  xlab("Period") + ylab("Estimate of Social \n Differentiation (S)") 
graph.socdiff
ggsave(graph.socdiff, file = "graph.socdiff.jpeg")


# 3. Calculation of dyadic PPI values (obs and null) ####

#first make a df with each dyad repeated twice - ID1 ID2 then ID2 ID1 - bc we need to rank from each female's perspective
ppi_rawdat <- w5_rates1yr %>% 
  select(ID1, ID2, period, w5_rate, kinship) %>%
  full_join(w5_rates1yr %>% select(ID1, ID2, period, w5_rate, kinship), 
            by = c(ID1 = "ID2", ID2 = "ID1", period = "period", w5_rate = "w5_rate", kinship = "kinship"))

#for each female need to order by top partner through last within each period
#will do that inside the same loop where we identify top partners

ppi_template <- ppi_rawdat %>% distinct(ID1, period) %>%
  mutate(pp1 = NA, 
         pp2 = NA)
list.ppiranks <- vector("list", length = nrow(ppi_template))
for(i in 1:length(ppi_template$pp1)){
  ID <- ppi_template$ID1[i]
  pd <- ppi_template$period[i]
  minidat <- ppi_rawdat %>% filter(ID1 == ID & period == pd) %>%
    arrange(desc(w5_rate)) %>%
    mutate(pprank = rank(desc(w5_rate)))
  list.ppiranks[[i]] <- minidat
  ppi_template$pp1[i] <- minidat$ID2[1]
  ppi_template$pp2[i] <- minidat$ID2[2]
}

ppiranks1yr <- do.call("rbind", list.ppiranks) %>% data.frame() %>%
  mutate(pprank = round(pprank, digits = 0))  %>%
  mutate(pprank_factor = as.factor(pprank)) %>%
  mutate(period = as.numeric(as.character(period)))

#make this for the null networks too:
list.ppirand_top2 <- vector("list", length = 1000)
t <- Sys.time()
for(i in seq(1000)){
  dat <- list.w5rand1yr[[i]]
  #need to make a new dataset with each dyad repeated twice in each order
  minidata <- dat %>% 
    full_join(dat %>% select(ID1, ID2, period, w5rate), 
              by = c(ID1 = "ID2", ID2 = "ID1", period = "period", w5rate = "w5rate")) %>%
    select(-w5_count, -scans_dyad)
  ppi_template <- minidata %>% distinct(ID1, period)
  
  list.ppiranks <- vector("list", length = nrow(ppi_template))
  for(j in 1:length(ppi_template$ID1)){
    ID <- ppi_template$ID1[j]
    pd <- ppi_template$period[j]
    minidat2 <- minidata %>% filter(ID1 == ID & period == pd) %>%
      arrange(desc(w5rate)) %>%
      mutate(pprank = rank(desc(w5rate)))
    list.ppiranks[[j]] <- minidat2
  }
  ppiranks1yr <- do.call("rbind", list.ppiranks) %>% data.frame() %>%
    mutate(pprank = round(pprank, digits = 0))  %>%
    mutate(period = as.numeric(as.character(period)))
  
  #new dataset with ppi
  ppi <- ppiranks1yr %>% distinct(ID1, period) %>%
    mutate(years_comp = NA,
           n = NA,
           u = NA, 
           x = NA, 
           mean_w5 = NA)
  for(k in 1:length(ppi$ID1)){
    who <- ppi$ID1[k]
    t1 <- ppi$period[k]
    t2 <- ppi$period[k] + 1
    minidata3 <- ppiranks1yr %>% 
      filter(ID1 == who & period == t1 | ID1 == who & period == t2) %>%
      filter(pprank <= 2)
    #number of partners we are comparing (bc for HH she only has 1 top partner because all others tie):
    ppi$n[k] <- minidata3 %>% filter(period == t1) %>% nrow() 
    #how many years being compared (if no t2 then years_comp = 1 and can be elim below):
    ppi$years_comp[k] <- minidata3 %>% distinct(period) %>% nrow()
    #sum of distinct top partners:
    ppi$u[k] <- minidata3 %>% distinct(ID2) %>% nrow()
    #number of partners available in t1 that were not also available as potential partner at t2:
    minidata4 <- minidata3 %>% 
      filter(period == t1) %>%
      left_join(ppiranks1yr %>% mutate(present = "Y") %>% 
                  filter(ID1 == who & period == t2) %>% select(ID2, present), 
                by = c("ID2" = "ID2")) %>%
      filter(!is.na(present)) 
    #max 2 females would be missing in the next year, so minus number who were available the next year
    ppi$x[k] <- 2-nrow(minidata4)
    #what is meanw5 rate for that female for her top partners during those two years? 
    ppi$mean_w5[k] <- mean(minidata4$w5rate)
  }
  
  #we only want to calculate ppi for females where 2 years are available to compare, 
  #and are comparing 3 top partners (sometimes there is a tie)
  #and are comparing max 6 distinct top partners (sometimes there is a tie)
  #and are comparing only females where at least 1 of their top partners was still available the next year
  ppi2 <- ppi %>% filter(years_comp == 2 & n == 2 & u <= 6 & x <=2) %>%
    mutate(t1 = period, 
           t2 = period + 1) %>% 
    select(ID1, t1, t2, n, u, x, mean_w5) %>%
    mutate(ppi_calc = (2*n - u)/(2*n - n - x)) %>%
    group_by(t1) %>%
    summarise(yr_ppi = mean(ppi_calc))
  
  list.ppirand_top2[[i]] <- ppi2
}
Sys.time() - t
save(list.ppirand_top2, file = "list_ppirand_top2.Rdata")

ppi_perm1yr_top2 <- do.call("rbind", list.ppirand_top2) %>% data.frame()
save(ppi_perm1yr_top2, file = "ppi_perm1yr_top2.Rdata")

ppi_2010 <- ppi_perm1yr_top2 %>% filter(t1 == 2010) %>% select(yr_ppi)
hist(ppi_2010$yr_ppi)
ppi2010.95 <- get_confidence_interval(ppi_2010, level = .95) #this is from package infer
ppi2010.80 <- get_confidence_interval(ppi_2010, level = .80)

ppi_2011 <- ppi_perm1yr_top2 %>% filter(t1 == 2011) %>% select(yr_ppi)
hist(ppi_2011$yr_ppi)
ppi2011.95 <- get_confidence_interval(ppi_2011, level = .95) #this is from package infer
ppi2011.80 <- get_confidence_interval(ppi_2011, level = .80)

ppi_2012 <- ppi_perm1yr_top2 %>% filter(t1 == 2012) %>% select(yr_ppi)
hist(ppi_2012$yr_ppi)
ppi2012.95 <- get_confidence_interval(ppi_2012, level = .95) #this is from package infer
ppi2012.80 <- get_confidence_interval(ppi_2012, level = .80)

ppi_2013 <- ppi_perm1yr_top2 %>% filter(t1 == 2013) %>% select(yr_ppi)
hist(ppi_2013$yr_ppi)
ppi2013.95 <- get_confidence_interval(ppi_2013, level = .95) #this is from package infer
ppi2013.80 <- get_confidence_interval(ppi_2013, level = .80)

ppi_2014 <- ppi_perm1yr_top2 %>% filter(t1 == 2014) %>% select(yr_ppi)
hist(ppi_2014$yr_ppi)
ppi2014.95 <- get_confidence_interval(ppi_2014, level = .95) #this is from package infer
ppi2014.80 <- get_confidence_interval(ppi_2014, level = .80)

ppi_2015 <- ppi_perm1yr_top2 %>% filter(t1 == 2015) %>% select(yr_ppi)
hist(ppi_2015$yr_ppi)
ppi2015.95 <- get_confidence_interval(ppi_2015, level = .95) #this is from package infer
ppi2015.80 <- get_confidence_interval(ppi_2015, level = .80)

ppi_2016 <- ppi_perm1yr_top2 %>% filter(t1 == 2016) %>% select(yr_ppi) %>% filter(!is.na(yr_ppi))
hist(ppi_2016$yr_ppi) #sometimes we get ties that results in NaN or inf calculations (I think) so I'm taking them out
ppi2016.95 <- get_confidence_interval(ppi_2016, level = .95) #this is from package infer
ppi2016.80 <- get_confidence_interval(ppi_2016, level = .80)

ppi_2017 <- ppi_perm1yr_top2 %>% filter(t1 == 2017) %>% select(yr_ppi) %>% filter(!is.na(yr_ppi))
hist(ppi_2017$yr_ppi) #sometimes we get ties that results in NaN or inf calculations (I think) so I'm taking them out
ppi2017.95 <- get_confidence_interval(ppi_2017, level = .95) #this is from package infer
ppi2017.80 <- get_confidence_interval(ppi_2017, level = .80)

ppi_2018 <- ppi_perm1yr_top2 %>% filter(t1 == 2018) %>% select(yr_ppi)
hist(ppi_2018$yr_ppi)
ppi2018.95 <- get_confidence_interval(ppi_2018, level = .95) #this is from package infer
ppi2018.80 <- get_confidence_interval(ppi_2018, level = .80)

#there is no 2019 data bc no 2020 data to compare with 
# ppi_2019 <- ppi_perm1yr_top2 %>% filter(t1 == 2019) %>% select(yr_ppi)
# hist(ppi_2019$yr_ppi)
# ppi2019.95 <- get_confidence_interval(ppi_2019, level = .95) #this is from package infer
# ppi2019.80 <- get_confidence_interval(ppi_2019, level = .80)


CI95_ppi <- ppi2010.95 %>%
  rbind(ppi2011.95, ppi2012.95, ppi2013.95, ppi2014.95, ppi2015.95, ppi2016.95, ppi2017.95, ppi2018.95) %>%
  mutate(period = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")) %>%
  mutate(period = as.numeric(period)) %>%
  rename(lower95 = lower_ci, 
         upper95 = upper_ci)

CI80_ppi <- ppi2010.80 %>%
  rbind(ppi2011.80, ppi2012.80, ppi2013.80, ppi2014.80, ppi2015.80, ppi2016.80, ppi2017.80, ppi2018.80) %>%
  mutate(period = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")) %>%
  mutate(period = as.numeric(period)) %>%
  rename(lower80 = lower_ci, 
         upper80 = upper_ci)

ppi_perm_summary1yr_top2 <- ppi_perm1yr_top2 %>% 
  rename(period = "t1") %>%
  group_by(period) %>%
  summarise(mean.ppi = mean(yr_ppi), 
            sd.ppi = sd(yr_ppi), 
            n.ppi = n(), 
            max.ppi = max(yr_ppi), 
            min.ppi = min(yr_ppi), 
            se.ppi = sd.ppi/sqrt(n.ppi)) %>%
  left_join(CI80_ppi) %>%
  left_join(CI95_ppi)

save(ppi_perm_summary1yr_top2, file = "ppi_perm_summary1yr_top2.Rdata")


# 4. Stability analysis 1: PPI vs null distribution ####
# calculate PPI for top 2 partners
#new dataset with ppi
ppi_top2 <- pp_who %>% select(ID1, period) %>%
  mutate(years_comp = NA,
         n = NA,
         u = NA, 
         x = NA, 
         mean_w5 = NA)
for(i in 1:length(ppi_top2$ID1)){
  who <- ppi_top2$ID1[i]
  t1 <- ppi_top2$period[i]
  t2 <- ppi_top2$period[i] + 1
  minidata <- ppiranks1yr %>% 
    filter(ID1 == who & period == t1 | ID1 == who & period == t2) %>%
    filter(pprank <= 2)
  #number of partners we are comparing (bc for HH she only has 1 top partner because all others tie):
  ppi_top2$n[i] <- minidata %>% filter(period == t1) %>% nrow() 
  #how many years being compared (if no t2 then years_comp = 1 and can be elim below):
  ppi_top2$years_comp[i] <- minidata %>% distinct(period) %>% nrow()
  #sum of distinct top partners:
  ppi_top2$u[i] <- minidata %>% distinct(ID2) %>% nrow()
  #number of partners available in t1 that were not also available as potential partner at t2:
  minidata2 <- minidata %>% 
    filter(period == t1) %>%
    left_join(ppiranks1yr %>% mutate(present = "Y") %>% 
                filter(ID1 == who & period == t2) %>% select(ID2, present), 
              by = c("ID2" = "ID2")) %>%
    filter(!is.na(present)) 
  #max 2 females would be missing in the next year, so minus number who were available the next year
  ppi_top2$x[i] <- 2-nrow(minidata2)
  #what is meanw5 rate for that female for her top partners during those two years? 
  ppi_top2$mean_w5[i] <- mean(minidata$w5_rate)
}

ppi_top2_2 <- ppi_top2 %>% filter(years_comp == 2 & n == 2) %>%
  mutate(t1 = period, 
         t2 = period + 1) %>% 
  select(ID1, t1, t2, n, u, x, mean_w5) %>%
  mutate(ppi_calc = (2*n - u)/(2*n - n -x))


# Graph PPI top 2 vs null
#make a summary table with avg PPI per year (all females)
# "ppi_perm_summary1yr_top2.Rdata" is the data summarized from the 1000 null datasets

ppi_summary_top2 <- ppi_top2_2 %>% group_by(t1) %>%
  summarise(yr_ppi = mean(ppi_calc)) %>%
  left_join(ppi_perm_summary1yr_top2, by = c("t1" = "period")) %>%
  rename(period = t1)

graph.ppitop2 <- ggplot(ppi_summary_top2, aes(x = period, y = yr_ppi, group = 1)) + 
  geom_point() +
  #geom_line(col='red') +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.3) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), alpha = 0.4) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_blank(), 
        text = element_text(size = 10)) +
  xlab("Period") + ylab("Partner Preference Index (PPI) \n (top two partners)") 
graph.ppitop2


# 5. Stability analysis 2: PPI logistic regression ####

#does assoc strength one year (w5) predict whether you will be a top p next year? 
# this permuted 1000 times, example code for how to run permutations is below 
pp_log_2.data <- pp_log_2 %>% filter(pprank <= 2) %>% filter(period != 2019) %>% filter(!is.na(nextyear))
log.top2p <- glmmTMB(nextyear ~ 1 
                     + scale(w5_rate)
                     + kinship
                     + (1|ID1)
                     + (1|ID2),
                     family = binomial, 
                     data = pp_log_2)


# 6. Example joint movements model ####
# this permuted 1000 times, example code for how to run permutations is below 
wt_ttf.data <- behav_undirected %>% filter(dpc>0) 
lmer.wt_ttf <- glmmTMB(wt_ttf_count ~  
                         scale(P5dyad_rate) +
                         kinship +
                         #kinship*scale(P5dyad_rate) +
                         #log(dpc) +
                         offset(log(dpc)) +
                         (1|ID1) +
                         (1|ID2),
                       data = wt_ttf.data, 
                       family = nbinom1) 

# 7. Example mutual association model ####
# this permuted 1000 times, example code for how to run permutations is below 
hinde5.data <- behav_undirected %>% filter(!is.na(hinde_index_5)) %>% mutate(hinde_index_5 = hinde_index_5/100)
lmer.hinde5 <- glmmTMB(hinde_index_5 ~
                         scale(P5dyad_rate) +
                         #I(scale(P5dyad_rate)^2) +
                         kinship +
                         #kinship*scale(P5dyad_rate)+
                         #offset(log(dpc)) + (no relationship between obs time and hinde so don't need to put in)
                         (1|ID1) +
                         (1|ID2),
                       data = hinde5.dat,
                       #family = gaussian)
                       #ziformula = ~1,
                       family = tweedie(link = "log"))


# 8. Example aggression model ####
# this permuted 1000 times, example code for how to run permutations is below 
agg.data <- dyad_w5_agg2
lmer.agg <- glmmTMB(agg_counts1yr ~  
                      scale(w5_rate) +
                      kinship +
                      #kinship*scale(P5dyad_rate) +
                      #log(dpc) +
                      offset(log(partycount_d)) +
                      (1|ID1) +
                      (1|ID2),
                    data = agg.dat, 
                    family = nbinom2) 

# 9. Example cofeeding model ####
# this permuted 1000 times, example code for how to run permutations is below 
mcdur.data <- behav_undirected %>% filter(focal_period_dpc_mcdur > 0)
lmer.mcdur <- glmmTMB(mc_dur ~ 
                        scale(P5nf_dyad_rate) +
                        kinship +
                        offset(log(focal_period_dpc_mcdur)) +
                        (1|ID1) +
                        (1|ID2),
                      data = mcdur.data,
                      family = tweedie(link ="log"))

# 10. Example vigilance model ####

lmer.scans <- glmmTMB(scan_post ~  
                        scale(P5dyad_rate) +
                        kinship +
                        #kinship*scale(P5dyad_rate) +
                        #offset(log(A5forAV)) +
                        (1|Actor) +
                        (1|Receiver),
                      data = AVS_logdata,
                      family = binomial)

# 11. Example scratching model ####

lmer.scratches <- glmmTMB(scratch_post ~1 +  
                            scale(P5dyad_rate) +
                            kinship +
                            #kinship*scale(P5dyad_rate) +
                            #offset(log(A5forAV)) +
                            (1|Actor) +
                            (1|Receiver),
                          data = CH_logdata,
                          family = binomial)

# 12. Example model checks ####

summary(lmer.mcdur)
testDispersion(lmer.mcdur)
plotResiduals(simulateResiduals(lmer.mcdur))
testResiduals(lmer.mcdur)
check_collinearity(lmer.mcdur)


# 13. Example of permuting models and summarize results ####

set.seed(5)
list.pp_log_2.coef <- vector("list", length = 1000)
list.pp_log_2.confint <- vector("list", length = 1000)
pp_log_2.data <- pp_log_2 %>% filter(pprank <= 2) %>% filter(period != 2019) %>% filter(!is.na(nextyear))
pp_log_2.datnoID <- pp_log_2.data %>% select(-ID1, -ID2)

t <- Sys.time()
for(i in seq(1000)){
  converge_check <- FALSE
  while(isFALSE(converge_check)) {
    
    IDs <- pp_log_2.data %>% select(ID1, ID2) %>%
      apply(., 1, sample) %>% t() %>% data.frame() %>%
      rename(ID1 = X1, ID2 = X2)
    pp_log_2.modeldat <- pp_log_2.datnoID %>% cbind(., IDs) 
    
    log.top2p <- glmmTMB(nextyear ~ 1 
                         + scale(w5_rate)
                         + kinship
                         #+ (1|period)
                         + (1|ID1)
                         + (1|ID2),
                         family = binomial, 
                         data = pp_log_2.modeldat)
    
    if(!is.na(summary(log.top2p)[["AICtab"]][["AIC"]])){
      if(!is.na(summary(log.top2p)$coefficients$cond[[1, 2]])){
        
        list.pp_log_2.coef[[i]]<- summary(log.top2p)$coefficients$cond %>%
          data.frame() %>%
          mutate(run = i) %>%
          rownames_to_column(var = "parameter")
        
        list.pp_log_2.confint[[i]] <- confint(log.top2p, full = T) %>%
          data.frame() %>%
          mutate(run = i) %>%
          rownames_to_column(var = "parameter")
        
        converge_check <- TRUE
      }}}}
Sys.time() - t 


summary(log.top2p)
testDispersion(log.top2p)
plot(simulateResiduals(log.top2p))
diagnose(log.top2p)

# Summarized results log top 2  

confint.pp_log_2.perm <- do.call("rbind", list.pp_log_2.confint)
confint.pp_log_2.avg <- confint.pp_log_2.perm %>%
  mutate(parameter = gsub("cond.", "", parameter)) %>% 
  group_by(parameter) %>%
  summarise(avg.ci.low = mean(X2.5..), 
            avg.ci.high = mean(X97.5..))

coefs.pp_log_2.perm <- do.call("rbind", list.pp_log_2.coef) #PICK UP HERE
pp_log_2.avg <- coefs.pp_log_2.perm %>%
  group_by(parameter) %>%
  summarise(n = n(), #just want to check n to make sure we are getting all the data
            Est_mean = mean(Estimate), 
            SE_mean = mean(Std..Error), 
            z_mean = mean(z.value),
            pval_mean = mean(Pr...z..),
            pval_prop = length(which(Pr...z.. < 0.05))/n*100) %>% 
  ungroup() %>%
  left_join(confint.pp_log_2.avg) %>%
  column_to_rownames(., var = "parameter") %>%
  select(-z_mean, -n) %>%
  round(., digits = 3)  %>%
  mutate(Model = ifelse(row_number() == 1, "Top Two Partner Consistency by Strength", "")) %>%
  rownames_to_column()

# 14. Example model comparison ####
set.seed(5)
list.agg3grm.coef <- vector("list", length = 1000)
list.agg3grm.confint <- vector("list", length = 1000)
list.agg3grm.dredge <- vector("list", length = 1000)
agg3grm.data <- dyad_w5_agg2
agg3grm.datnoID <- agg3grm.data %>% select(-ID1, -ID2)

options(na.action = "na.fail")
t <- Sys.time()
for(i in seq(1000)){
  converge_check <- FALSE
  while(isFALSE(converge_check)) {
    IDs <- agg3grm.data %>% select(ID1, ID2) %>%
      apply(., 1, sample) %>% t() %>% data.frame() %>%
      rename(ID1 = X1, ID2 = X2)
    agg3grm.modeldat <- agg3grm.datnoID %>% cbind(., IDs) 
    lmer.agg3grm <- glmmTMB(agg_counts1yr ~  
                              scale(w5_rate) + #how often in 5 given how often each is scanned (during a focal)
                              scale(grm_rate) + #rate of grooming given how often each is scanned (during a focal)
                              kinship +
                              offset(log(partycount_d)) +
                              (1|ID1) +
                              (1|ID2),
                            data = agg3grm.modeldat, 
                            family = nbinom2) 
    if(!is.na(summary(lmer.agg3grm)[["AICtab"]][["AIC"]])){ #just checking that it did converge before proceeding
      if(!is.na(summary(lmer.agg3grm)$coefficients$cond[[1, 2]])){ #just checking that it did converge before proceeding
        
        dredge.agg3grm <- dredge(lmer.agg3grm, fixed = c("cond(kinship)", "cond(offset(log(partycount_d)))"))
        list.agg3grm.dredge[[i]] <- dredge.agg3grm %>%
          data.frame() %>%
          mutate(run = i) %>%
          rownames_to_column(var = "model_num") %>%
          mutate(rank = c(1, 2, 3, 4))
        
        #getting confidence intervals for each model w/n dredge is annoying to do while keeping track of which model the CIs belong to
        models <- (get.models(dredge.agg3grm, TRUE)) #get details on each model
        ci <- lapply(models, confint) #pull CIs from each model
        ci1 <- as.data.frame(ci$'1') %>% mutate(model = "Null Model")
        ci2 <- as.data.frame(ci$'2') %>% mutate(model = "Grooming")
        ci3 <- as.data.frame(ci$'3') %>% mutate(model = "Proximity")
        ci4 <- as.data.frame(ci$'4') %>% mutate(model = "Proximity and Grooming")
        list.agg3grm.confint[[i]] <- rbind(ci1, ci2, ci3, ci4) %>% #compile all CIs into table
          data.frame() %>%
          mutate(run = i) %>%
          rownames_to_column(var = "parameter")
        
        converge_check <- TRUE
      }}}}
Sys.time() - t #16 seconds for two


