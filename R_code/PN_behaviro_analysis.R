## For anti behavior tracking data analysis
## behavior videos are analyzed with deeplabcut
## created on 03192020
## last updated 08042020

cc_anti_behavior <- function(path_ID, frame_rate, length_pix) {
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   ## only select the position of head
    as.data.frame() %>% 
    mutate(V2 = V2- V8, V3 = V3- V9, V5 = V5-V8, V6 = V6 -V9) %>%  ##normalize from FIJI
    select(str_c("V", 1:6)) 
  
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")
  
  ## remove the frame for the first few frame
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    subset(., Prob >0.5 | Frame> 60)
  
  ## change the dim from pix to mm (the dim of one chammer is 165mm)
  # value 387 for anti experiment 20200309
  # value 447 for miniscope 20191209
  r_pix_mm <- 165 / length_pix
  dat_anti <- dat_anti %>% 
    mutate(Head_x = abs(Head_x * r_pix_mm), Head_y = abs(Head_y * r_pix_mm)) %>% 
    mutate(Tail_x = abs(Tail_x * r_pix_mm), Tail_y = abs(Tail_y * r_pix_mm))
    
  chamber_div <- 165.5 ## the lenght of two chambels and the gap between
  
  ## assign each frame in hot as 1 (48) and nor plate(30) as 2
  dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
  ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[2]/nrow(dat_anti))
  
  p_anti <- ggplot(dat_anti, aes(Head_y, Head_x, colour=chamber))+
    geom_point(shape=1)+
    theme_void()+
    theme(legend.position = "none")
  
  ## calculate how many times the mouse cross the border, and latency
  digit <- ifelse(dat_anti$Head_y< chamber_div, 1,2)
  
  for (i in c(2: (length(digit)-20))) {
    digit[i]<- ifelse(digit[i]!= digit[i-1] & digit[i+20] == digit[i], digit[i], digit[i-1] )
  }
  
  dat_anti$digit<- digit
  cross_digit <- diff(digit)
  cross_hot_nor <- which(cross_digit==1)+1
  
  cross_nor_hot <- which(cross_digit==-1)+1
  
  total_cross <- length(c(cross_hot_nor, cross_nor_hot))
  # calculate the latency of fist moving 
  # frame_rate <- 10 ##10 Hz
  cross_latency_1st  <- ifelse(digit[1]==2,2, cross_hot_nor[1]/10)
  cross_latency_1st <- ifelse(cross_latency_1st <2, 2, cross_latency_1st)
  
  # latency of cross back from nor to hot
  cross_latency_2nd <- ifelse(length(cross_nor_hot)==0, nrow(dat_anti)/10 - cross_latency_1st, 
                              cross_nor_hot[1]/10 - cross_latency_1st)
  
  
  ## calculate the distance and speed of movement
  x_diff <- c(0, diff(dat_anti$Head_x))
  y_diff <- c(0, diff(dat_anti$Head_y))
  dis_frame <- sqrt(x_diff^2 + y_diff^2)
  move_speed <- dis_frame*10
  
  move_speed_t <- zoo(move_speed) %>% 
    rollapply(., width = 10, by = 10, FUN = mean, na.rm = TRUE, align = "left") %>% 
    as.vector()
  peak_threshold <- mean(move_speed_t) +  sd(move_speed_t)
  peak_speed <- findpeaks(move_speed_t, nups = 2,minpeakdistance = 2, minpeakheight = peak_threshold) %>% 
    .[,2]*10 %>% 
    sort()
  latency_accerlation <- (peak_speed[which(peak_speed> cross_hot_nor[1])][1] - cross_hot_nor[1])/10
  # total distance, normalization
  total_distance <- round(sum(dis_frame), digits = 2)
  
  # speed in different chambers
  speed_hot <- mean(move_speed[which(dat_anti$chamber == "Hot")])
  speed_nor <- mean(move_speed[which(dat_anti$chamber == "Nor")])
  speed_compare <- speed_nor/speed_hot
  # average speed (3s) before crossing the chambers
  speed_to_nor_1st <- ifelse(cross_hot_nor[1]< 20, mean(move_speed[0:cross_hot_nor[1]]) ,mean(move_speed[(cross_hot_nor[1]-20):cross_hot_nor[1]]))
  speed_to_hot_1st <- ifelse(length(cross_nor_hot)==0, NA, mean(move_speed[(cross_nor_hot[1]-20):cross_nor_hot[1]]))
  
  dat_anti_result <- tibble(ratio=ratio_time, latency1= cross_latency_1st, latency2 = cross_latency_2nd,
                            total_distance = total_distance, speed1 = speed_to_nor_1st, speed2= speed_to_hot_1st, total_cross= total_cross,
                            speed_compare=speed_compare, latency_accerlation= latency_accerlation)
  return(dat_anti_result)
}

## for opto anti behavior------
## for mice m20
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m20")

path_trace_m20 <- list.files( pattern = "*.csv" )
dat_behavior_m20 <- mapply(cc_anti_behavior, path_trace_m20, frame_rate = 10,length_pix=387, SIMPLIFY = F)

## for mice m19
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m19")

path_trace_m19 <- list.files( pattern = "*.csv" )
dat_behavior_m19 <- mapply(cc_anti_behavior, path_trace_m19, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m22
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m22")

path_trace_m22 <- list.files( pattern = "*.csv" )
dat_behavior_m22 <- mapply(cc_anti_behavior, path_trace_m22, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m23
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m23")

path_trace_m23 <- list.files( pattern = "*.csv" )
dat_behavior_m23 <- mapply(cc_anti_behavior, path_trace_m23, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m63
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m63")

path_trace_m63 <- list.files( pattern = "*.csv" )
dat_behavior_m63 <- mapply(cc_anti_behavior, path_trace_m63, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m851
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m851")

path_trace_m851 <- list.files( pattern = "*.csv" )
dat_behavior_m851 <- mapply(cc_anti_behavior, path_trace_m851, frame_rate = 10, length_pix=387,SIMPLIFY = F)


## for miniscope implanted behavior----
# for m3
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m3")

path_trace_m3 <- list.files( pattern = "*.csv" )[c(3,7)]
dat_behavior_m3 <- mapply(cc_anti_behavior, path_trace_m3, frame_rate = 20, length_pix=478,SIMPLIFY = F) %>% 
  do.call(rbind,.)

# for m7
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m7")

path_trace_m7 <- list.files( pattern = "*.csv" )[c(3,7)]
dat_behavior_m7 <- mapply(cc_anti_behavior, path_trace_m7, frame_rate = 20, length_pix=478,SIMPLIFY = F) %>% 
  do.call(rbind,.)


# for m17
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m17")

path_trace_m17 <- list.files( pattern = "*.csv" )[c(3,7)]
dat_behavior_m17 <- mapply(cc_anti_behavior, path_trace_m17, frame_rate = 20, length_pix=387,SIMPLIFY = F) %>% 
  do.call(rbind,.)

# for m18
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m18")

path_trace_m18 <- list.files( pattern = "*.csv" )[c(3,7)]
dat_behavior_m18 <- mapply(cc_anti_behavior, path_trace_m18, frame_rate = 20, length_pix=387,SIMPLIFY = F) %>% 
  do.call(rbind,.)

# for m855
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m855")

path_trace_m855 <- list.files( pattern = "*.csv" )
dat_behavior_m855 <- mapply(cc_anti_behavior, path_trace_m855, frame_rate = 10, length_pix=364,SIMPLIFY = F) %>% 
  do.call(rbind,.)

# for m857
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m857")

path_trace_m857 <- list.files( pattern = "*.csv" )
dat_behavior_m857 <- mapply(cc_anti_behavior, path_trace_m857, frame_rate = 10, length_pix=364,SIMPLIFY = F) %>% 
  do.call(rbind,.)

## combine behavior and miniscope data
dat_mini_area <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  arrange(match(ID, mouse_ID))

dat_mini_ratio <- dat_cell_area %>% 
  select(ID,Day, ratio) %>% 
  arrange(match(ID, mouse_ID))


dat_behavior <- rbind(dat_behavior_m3, dat_behavior_m7, dat_behavior_m17, dat_behavior_m18, dat_behavior_m855, dat_behavior_m857) %>% 
  mutate(ID = rep(mouse_ID, each=2), Day = rep(c("Pre", "Test"), length(mouse_ID))) %>% 
  mutate(Sum = dat_mini_area$sum, Ratio = dat_mini_ratio$ratio) 


p_sum_latency_test <- dat_behavior %>% 
  filter(Day =="Test") %>% 
  ggplot(., aes(Sum, latency2))+
  geom_point()+
  geom_smooth(method = "lm", colour="indianred" )+
  labs(x="AUG (z-score*s)", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_sum_latency_pre <- dat_behavior %>% 
  filter(Day =="Pre") %>% 
  ggplot(., aes(Sum, latency2))+
  geom_point()+
  geom_smooth(method = "lm", colour="deepskyblue4" )+
  labs(x="AUG (z-score*s)", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_sum_latency_pre_cor <- dat_behavior %>% 
  filter(Day =="Pre") %>% 
  cor.test(~latency2 + Sum,.)

p_sum_latency_test_cor <- dat_behavior %>% 
  filter(Day =="Test") %>% 
  cor.test(~latency2 + Sum,.)

p_sum_latency_test <- dat_behavior %>% 
  filter(Day =="Test") %>% 
  ggplot(., aes(Sum, ratio*100))+
  geom_point()+
  geom_smooth(method = "lm", colour="indianred" )+
  labs(x="AUG (z-score*s)", y="Time spend in chammber 2 (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_sum_latency_pre <- dat_behavior %>% 
  filter(Day =="Pre") %>% 
  ggplot(., aes(Sum, ratio*100))+
  geom_point()+
  geom_smooth(method = "lm", colour="deepskyblue4" )+
  labs(x="AUG (z-score*s)", y="Time spend in chammber 2 (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


p_anti_cor <- plot_grid(p_sum_latency_pre, p_sum_latency_test, nrow = 1)


# correlation between ca2+ activity and pain behavior
dat_anti_pain <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/pn_anti_pain_behavior.xlsx") %>% 
  mutate(Latency_licking = (First_licking - First_crossing)/Frame_rate,
         Latency_rearing = (First_rearing - First_crossing)/Frame_rate) 

linear_scale <- function(x){
  y = (x-min(x))/(max(x)- min(x))
  return(y)
}
dat_anti_pain_cor <- dat_anti_area %>% 
  filter(Day =="Test") %>% 
  arrange(match(ID, mouse_ID)) %>% 
  select(sum) %>% 
  cbind(., dat_anti_pain[,c("ID", "Latency_licking", "Latency_rearing")]) %>% 
  mutate(., Latency_licking = linear_scale(Latency_licking),Latency_rearing = linear_scale (Latency_rearing) ) %>% 
  pivot_longer(.,-c(sum,ID), names_to = "Variable" )

p_anti_behavior <- ggplot(dat_anti_pain_cor, aes(sum, value, colour = Variable))+
  geom_point()+
  geom_smooth(method = "lm" , se = F)+
  labs(x="AUC (s.d.*s)", y="Norm. latency")+
  scale_colour_manual(values=c( "black", "gray"))+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = c(0.2, 0.8))


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_behavior.pdf", width = 80/25.6, height =60/25.6, family = "Arial")
p_anti_behavior
dev.off()   

# p_anti_rearing <- ggplot(dat_anti_pain_cor, aes(sum, Latency_rearing))+
#   geom_point()+
#   geom_smooth(method = "lm",colour="indianred" )+
#   labs(x="AUC (z-score*s)", y="Latency of 1st rearing (s)")+
#   theme(axis.line.x = element_line(),
#         axis.line.y = element_line(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
#   theme(legend.position = 'none')
# 
# p_anti_pain_cor <- plot_grid(p_anti_licking, p_anti_rearing, nrow = 1)

## correlation test
t_cor_licking <- dat_anti_pain_cor %>% 
  cor.test(~Latency_licking + sum,.)

t_cor_rearing <- dat_anti_pain_cor %>% 
  cor.test(~Latency_rearing + sum,.)

p_cor <- plot_grid(p_anti_cor, p_anti_pain_cor, nrow = 2)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor.pdf", width = 80/25.6, height =60/25.6, family = "Arial")
p_anti_licking
dev.off()   
## for HM4Di behavior------
# for m10
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m10")

path_trace_m10 <- list.files( pattern = "*.csv" )
dat_behavior_m10 <- mapply(cc_anti_behavior, path_trace_m10, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m11
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m11")

path_trace_m11 <- list.files( pattern = "*.csv" )
dat_behavior_m11 <- mapply(cc_anti_behavior, path_trace_m11, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m12
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m12")

path_trace_m12 <- list.files( pattern = "*.csv" )
dat_behavior_m12 <- mapply(cc_anti_behavior, path_trace_m12, frame_rate = 10, length_pix=455, SIMPLIFY = F)

# for m13
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m13")

path_trace_m13 <- list.files( pattern = "*.csv" )
dat_behavior_m13 <- mapply(cc_anti_behavior, path_trace_m13, frame_rate = 10, length_pix=455, SIMPLIFY = F)


## combine all anti data and plot----
mouse_anti_list <- list(dat_behavior_m17, dat_behavior_m18, dat_behavior_m3, dat_behavior_m7, dat_behavior_m855)

## for 1st latency, day2-3 as Ctrl, day5-6 as anti and day7 in pain
latency_1st_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[2]])
latency_1st_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[2]])
latency_1st_ctrl <- rowMeans(cbind(latency_1st_ctrl_d2, latency_1st_ctrl_d3))

latency_1st_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[2]])
latency_1st_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[2]])
latency_1st_anti <- rowMeans(cbind(latency_1st_ctrl_d5, latency_1st_ctrl_d6))

latency_1st_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[2]])

anti_latency <- c(latency_1st_ctrl, latency_1st_anti, latency_1st_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(latency_1st_ctrl))
dat_anti_1st_latency <- data.frame(Group = anti_group, Latency = anti_latency)
dat_anti_1st_latency$Group <- factor(dat_anti_1st_latency$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_latency_sta <- ddply(dat_anti_1st_latency, .(Group), summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))

p_anti_1st_latency <- ggplot(dat_anti_1st_latency, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Latency of 1st border crossing (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 2nd latencyv(crossing back), day2-3 as Ctrl, day5-6 as anti and day7 in pain
latency_2nd_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[3]])
latency_2nd_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[3]])
latency_2nd_ctrl <- rowMeans(cbind(latency_2nd_ctrl_d2, latency_2nd_ctrl_d3))

latency_2nd_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[3]])
latency_2nd_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[3]])
latency_2nd_anti <- rowMeans(cbind(latency_2nd_ctrl_d5, latency_2nd_ctrl_d6))

latency_2nd_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[3]])

anti_2nd_latency <- c(latency_2nd_ctrl, latency_2nd_anti, latency_2nd_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(latency_2nd_ctrl))
dat_anti_2nd_latency <- data.frame(Group = anti_group, Latency = anti_2nd_latency)
dat_anti_2nd_latency$Group <- factor(dat_anti_2nd_latency$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_2nd_latency_sta <- ddply(dat_anti_2nd_latency, .(Group), summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))

p_anti_2nd_latency <- ggplot(dat_anti_2nd_latency, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 1st speed, day2-3 as Ctrl, day5-6 as anti and day7 in pain
speed_1st_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[5]])
speed_1st_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[5]])
speed_1st_ctrl <- rowMeans(cbind(speed_1st_ctrl_d2, speed_1st_ctrl_d3))

speed_1st_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[5]])
speed_1st_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[5]])
speed_1st_anti <- rowMeans(cbind(speed_1st_ctrl_d5, speed_1st_ctrl_d6))

speed_1st_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[5]])

anti_speed <- c(speed_1st_ctrl, speed_1st_anti, speed_1st_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(speed_1st_ctrl))
dat_anti_1st_speed <- data.frame(Group = anti_group, speed = anti_speed)
dat_anti_1st_speed$Group <- factor(dat_anti_1st_speed$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_speed_sta <- ddply(dat_anti_1st_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_1st_speed <- ggplot(dat_anti_1st_speed, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Speed of crossing (mm/s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 2nd speed (crossing back), day2-3 as Ctrl, day5-6 as anti and day7 in pain
speed_2nd_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[6]])
speed_2nd_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[6]])
speed_2nd_ctrl <- rowMeans(cbind(speed_2nd_ctrl_d2, speed_2nd_ctrl_d3), na.rm = T)

speed_2nd_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[6]])
speed_2nd_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[6]])[-4]
speed_2nd_anti <- rowMeans(cbind(speed_2nd_ctrl_d5, speed_2nd_ctrl_d6), na.rm = T)
speed_2nd_anti <- speed_2nd_anti[!is.na(speed_2nd_anti)]
speed_2nd_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[6]])

anti_2nd_speed <- c(speed_2nd_ctrl, speed_2nd_anti, speed_2nd_anti_pain)
anti_group_2nd <- rep(c("Ctrl", "Anti", "Anti_pain"), c(length(speed_2nd_ctrl), length(speed_2nd_anti), length(speed_2nd_anti_pain)))
dat_anti_2nd_speed <- data.frame(Group = anti_group_2nd, speed = anti_2nd_speed)
dat_anti_2nd_speed$Group <- factor(dat_anti_2nd_speed$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_2nd_speed_sta <- ddply(dat_anti_2nd_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_2nd_speed <- ggplot(dat_anti_2nd_speed, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Speed of crossing (mm/s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## for total number of crossing, day2-3 as Ctrl, day5-6 as anti and day7 in pain
num_cross_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[7]])
num_cross_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[7]])
num_cross_ctrl <- rowMeans(cbind(num_cross_ctrl_d2, num_cross_ctrl_d3))

num_cross_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[7]])
num_cross_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[7]])
num_cross_anti <- rowMeans(cbind(num_cross_ctrl_d5, num_cross_ctrl_d6))

num_cross_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[7]])

anti_num_cross <- c(num_cross_ctrl, num_cross_anti, num_cross_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(num_cross_ctrl))
dat_anti_num_cross <- data.frame(Group = anti_group, num = anti_num_cross)
dat_anti_num_cross$Group <- factor(dat_anti_num_cross$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_speed_sta <- ddply(dat_anti_1st_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_num_cross <- ggplot(dat_anti_num_cross, aes(Group, num, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Total crossing #")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for ratio of time in Nor chammber, day2-3 as Ctrl, day5-6 as anti and day7 in pain
ratio_time_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[1]])
ratio_time_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[1]])
ratio_time_ctrl <- rowMeans(cbind(ratio_time_ctrl_d2, ratio_time_ctrl_d3))

ratio_time_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[1]])
ratio_time_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[1]])
ratio_time_anti <- rowMeans(cbind(ratio_time_ctrl_d5, ratio_time_ctrl_d6))

ratio_time_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[1]])

anti_ratio_time <- c(ratio_time_ctrl, ratio_time_anti, ratio_time_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(ratio_time_ctrl))
dat_anti_ratio_time <- data.frame(Group = anti_group, ratio = anti_ratio_time)
dat_anti_ratio_time$Group <- factor(dat_anti_ratio_time$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_speed_sta <- ddply(dat_anti_1st_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_ratio_time <- ggplot(dat_anti_ratio_time, aes(Group, ratio, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## combine all figures
p_anti_com <- plot_grid(p_anti_1st_latency, p_anti_num_cross, p_anti_2nd_latency, 
                        p_anti_1st_speed, p_anti_ratio_time, p_anti_2nd_speed, nrow = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_com.pdf", width = 140/25.6, height = 120/25.6, family = "Arial")
p_anti_com
dev.off()



## compare ctrl and opto, combine all data and plot-----
mouse_opto <- list(dat_behavior_m19, dat_behavior_m20, dat_behavior_m23, dat_behavior_m63, dat_behavior_m851)
mouse_wt <- list(dat_behavior_m16, dat_behavior_m17, dat_behavior_m18, dat_behavior_m22, dat_behavior_m3, dat_behavior_m7, dat_behavior_m27)

## for time spend in nor plate
prop_wt <- sapply(mouse_wt, function(x) x[[7]][[1]])
prop_opto <- sapply(mouse_opto, function(x) x[[7]][[1]])
prop_comb <- c(prop_wt, prop_opto)
prop_group <- rep(c("Ctrl", "Opto"), c(length(prop_wt), length(prop_opto)))

dat_anti_prop <- data.frame(Group = prop_group, Value = prop_comb)
dat_anti_prop$Group <- factor(dat_anti_prop$Group, levels = c("Ctrl", "Opto"))

p_anti_prop <- ggplot(dat_anti_prop, aes(Group, Value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## for 1st latency to nor chamber
latency_1st_wt <- sapply(mouse_wt, function(x) x[[7]][[2]])
latency_1st_opto <- sapply(mouse_opto, function(x) x[[7]][[2]])
latency_1st_comb <- c(latency_1st_wt, latency_1st_opto)
latency_1st_group <- rep(c("Ctrl", "Opto"), c(length(latency_1st_wt), length(latency_1st_opto)))

dat_anti_latency_1st <- data.frame(Group = latency_1st_group, Latency = latency_1st_comb)
dat_anti_latency_1st$Group <- factor(dat_anti_latency_1st$Group, levels = c("Ctrl", "Opto"))

p_anti_latency_1st <- ggplot(dat_anti_latency_1st, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Latency of 1st border crossing")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## For the first latency to nor
latency_2nd_wt <- sapply(mouse_wt, function(x) x[[7]][[3]])
latency_2nd_opto <- sapply(mouse_opto, function(x) x[[7]][[3]])
latency_2nd_comb <- c(latency_2nd_wt, latency_2nd_opto)
latency_2nd_group <- rep(c("Ctrl", "Opto"), c(length(latency_2nd_wt), length(latency_2nd_opto)))

dat_anti_latency_2nd <- data.frame(Group = latency_2nd_group, Latency = latency_2nd_comb)
dat_anti_latency_2nd$Group <- factor(dat_anti_latency_2nd$Group, levels = c("Ctrl", "Opto"))

p_anti_latency_2nd <- ggplot(dat_anti_latency_2nd, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Latency of crossing back")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 1st speed to nor chamber
speed_1st_wt <- sapply(mouse_wt, function(x) x[[7]][[5]])
speed_1st_opto <- sapply(mouse_opto, function(x) x[[7]][[5]])
speed_1st_comb <- c(speed_1st_wt, speed_1st_opto)
speed_1st_group <- rep(c("Ctrl", "Opto"), c(length(speed_1st_wt), length(speed_1st_opto)))

dat_anti_speed_1st <- data.frame(Group = speed_1st_group, speed = speed_1st_comb)
dat_anti_speed_1st$Group <- factor(dat_anti_speed_1st$Group, levels = c("Ctrl", "Opto"))

p_anti_speed_1st <- ggplot(dat_anti_speed_1st, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="speed of 1st border crossing")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## Speed of crossing back
speed_2nd_wt <- sapply(mouse_wt, function(x) x[[7]][[6]])
speed_2nd_opto <- sapply(mouse_opto, function(x) x[[7]][[6]])
speed_2nd_comb <- c(speed_2nd_wt, speed_2nd_opto)
speed_2nd_group <- rep(c("Ctrl", "Opto"), c(length(speed_2nd_wt), length(speed_2nd_opto)))

dat_anti_speed_2nd <- data.frame(Group = speed_2nd_group, speed = speed_2nd_comb)
dat_anti_speed_2nd$Group <- factor(dat_anti_speed_2nd$Group, levels = c("Ctrl", "Opto"))

p_anti_speed_2nd <- ggplot(dat_anti_speed_2nd, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="speed of crossing back")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## total crossing time
n_cross_wt <- sapply(mouse_wt, function(x) x[[7]][[7]])
n_cross_opto <- sapply(mouse_opto, function(x) x[[7]][[7]])
n_cross_comb <- c(n_cross_wt, n_cross_opto)
n_cross_group <- rep(c("Ctrl", "Opto"), c(length(n_cross_wt), length(n_cross_opto)))

dat_anti_n_cross <- data.frame(Group = n_cross_group, Value = n_cross_comb)
dat_anti_n_cross$Group <- factor(dat_anti_n_cross$Group, levels = c("Ctrl", "Opto"))

p_anti_n_cross <- ggplot(dat_anti_n_cross, aes(Group, Value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_opto_com <- plot_grid(p_anti_latency_1st, p_anti_prop, p_anti_latency_2nd,p_anti_speed_2nd, nrow = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_opto_com.pdf", width = 80/25.6, height = 120/25.6, family = "Arial")
p_opto_com
dev.off()

## for HM4Di behavior------
# for m10
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m10")

path_trace_m10 <- list.files( pattern = "*.csv" )
dat_behavior_m10 <- mapply(cc_anti_behavior, path_trace_m10, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m11
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m11")

path_trace_m11 <- list.files( pattern = "*.csv" )
dat_behavior_m11 <- mapply(cc_anti_behavior, path_trace_m11, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m12
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m12")

path_trace_m12 <- list.files( pattern = "*.csv" )
dat_behavior_m12 <- mapply(cc_anti_behavior, path_trace_m12, frame_rate = 10, length_pix=455, SIMPLIFY = F)

# for m13
setwd("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m13")

path_trace_m13 <- list.files( pattern = "*.csv" )
dat_behavior_m13 <- mapply(cc_anti_behavior, path_trace_m13, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## combine and compare
mouse_HM4di <- list(dat_behavior_m10, dat_behavior_m11)
mouse_ctrl <- list(dat_behavior_m12, dat_behavior_m13)

## for time spend in nor plate
prop_ctrl <- sapply(mouse_ctrl, function(x) x[[7]][[1]])
prop_HM4di <- sapply(mouse_HM4di, function(x) x[[7]][[1]])
prop_comb <- c(prop_ctrl, prop_HM4di)
prop_group <- rep(c("Ctrl", "HM4di"), c(length(prop_ctrl), length(prop_HM4di)))

dat_anti_prop <- data.frame(Group = prop_group, Value = prop_comb)
dat_anti_prop$Group <- factor(dat_anti_prop$Group, levels = c("Ctrl", "HM4di"))

p_anti_prop <- ggplot(dat_anti_prop, aes(Group, Value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## for WT ctrl to prove the placebo analgesia-----
dat_anti_ctrl <- vector(mode = "list", 7)
dat_anti_con <- vector(mode = "list", 7)
length_pix <- c(363,362,363,364,362,364,365)
day_group <- c("Rm","Rm", "Pre", "Rm","Rm","Rm", "Test")
for(i in c(1:7)){
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",i)
  path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",i)
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="Ctrl", ID=str_c("m",1:9), .before = 'ratio')
  
  dat_anti_con[[i]]<- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="Cond.", ID=str_c("m",1:10), .before = 'ratio')
}

dat_anti_wt <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  dplyr::filter(Day!="Rm") %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  mutate(Day=factor(Day, levels = c("Pre",  "Test")), Group=factor(Group, levels = c("Ctrl", "Cond.")))

dat_anti_wt_sta <- ddply(dat_anti_wt,.(variable,Day, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  filter(Day=="Test")

p_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")

p_ratio_test <- dat_anti_wt %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_anti_wt %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
 scale_y_continuous(limits = c(0, 150), expand = c(0, 0))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(100,102,102,100),"path")+
  annotate("text",x=2,y=102, label="*", size=5)

p_latency1_test <- dat_anti_wt %>% 
  filter(variable=="latency1"& Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 90))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(70,72,72,70),"path")+
  annotate("text",x=2,y=72, label="**", size=5)

p_latency2_test <- dat_anti_wt %>% 
  filter(variable=="latency2"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_anti_wt %>% 
  filter(variable=="speed_compare") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.05,2.05,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_speed_compare_test <- dat_anti_wt %>% 
  filter(variable=="speed_compare"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_total_cross <- dat_anti_wt %>% 
  filter(variable=="total_cross"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.02,2.02,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_total_cross_test <- dat_anti_wt %>% 
  filter(variable=="total_cross"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_pain_compare1 <- plot_grid(p_ratio, p_latency1, p_latency2, p_speed_compare, ncol = 2)

## only plot the latency of crossing back and the ratio of staying
p_pain_compare1 <- plot_grid(p_latency2, p_ratio, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_compare1.pdf", width = 100/25.6, height = 70/25.6, family = "Arial")
p_pain_compare1
dev.off()

## for pain behavir, like licking, and wearing-----
# mannualy analyzed by chong (07062020)
dat_wt_manual <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl_07042020_manualy.xlsx") %>% 
  mutate('Rearing'= Total_rearing/Total_frame*600, 'First_rearing'=First_rearing*0.1, 
         'Guarding'= (Guarding - Frame_1st_crossing)*0.1, 'Acceleration'= (Acceleration -Frame_1st_crossing )*0.1) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  mutate(Group=factor(Group, levels=c("Ctrl", "Cond.")))

dat_wt_manual_sta <- dat_wt_manual %>% 
  ddply(., .( variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 50), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(40,42,42,40),"path")+
  annotate("text",x=1.5,y=42, label="**", size=5)

p_licking_test <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_wt_manual %>% 
  filter(variable=="First_rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(60,62,62,60),"path")+
  annotate("text",x=1.5,y=62, label="*", size=5)

p_1st_rearing_test <- dat_wt_manual %>% 
  filter(variable=="First_rearing") %>% 
  wilcox.test(value~Group,.)

p_rearing <- dat_wt_manual %>% 
  filter(variable=="Rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="# rearing per min ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(1, 10))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(9,9.1,9.1,9),"path")+
  annotate("text",x=1.5,y=9.1, label="*", size=5)

p_rearing_test<- dat_wt_manual %>% 
  filter(variable=="Rearing") %>% 
  wilcox.test(value~Group,.)

p_jump <- dat_wt_manual %>% 
  filter(variable=="Jump") %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(300, 2000))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(1850,1900,1900,1850),"path")+
  annotate("text",x=1.5,y=1900, label="*", size=5)

p_jump_test<- dat_wt_manual %>% 
  filter(variable=="Jump") %>% 
  wilcox.test(value~Group,.)

p_Acceleration <- dat_wt_manual %>% 
  filter(variable=="Acceleration") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency to 1st acceleration (s) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(75,77,77,75),"path")+
  annotate("text",x=1.5,y=77, label="***", size=5)

p_Acceleration_test<- dat_wt_manual %>% 
  filter(variable=="Acceleration") %>% 
  wilcox.test(value~Group,.)
p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_rearing, p_Acceleration, nrow = 1)
p_con <- plot_grid(p_pain_compare1, p_pain_compare2, nrow = 2, rel_heights = c(2, 1))


p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_rearing, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_compare2.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_pain_compare2
dev.off()

## correlation analysis of the latency2 and pain behavior
dat_anti_cor1 <- dat_anti_wt %>% 
  filter(variable=="latency2"& Day == "Test") %>% 
  arrange(ID, Group)

dat_anti_cor2 <- dat_wt_manual %>% 
  filter(variable == "First_licking") %>% 
  arrange(ID, Group)

dat_anti_cor3 <- dat_wt_manual %>% 
  filter(variable == "First_rearing") %>% 
  arrange(ID, Group)



dat_anti_cor <- tibble(Latency = dat_anti_cor1$value, Licking = dat_anti_cor2$value, Rearing = dat_anti_cor3$value)
p_cor_latency_licking <- dat_anti_cor %>% 
  ggplot(., aes(Latency, Licking))+
  geom_point()+
  geom_smooth(method = "lm") +
  labs(x="Latency of crossing back (s)", y="Latency of first licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

p_cor_latency_rearing <- dat_anti_cor %>% 
  ggplot(., aes(Latency, Rearing))+
  geom_point()+
  geom_smooth(method = "lm") +
  labs(x="Latency of crossing back (s)", y="Latency of first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

p_cor <- plot_grid(p_cor_latency_licking, p_cor_latency_rearing, ncol = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_cor
dev.off()

## correlation test of
p_cor_test_likcing <- dat_anti_cor %>% 
  cor.test(~ Latency + Licking, .)

p_cor_test_likcing <- dat_anti_cor %>% 
  cor.test(~ Latency + Rearing, .)
## heat map plot to show the path----

cc_behavior_reshape <- function(path_ID, frame_rate, length_pix) {
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F)
  
  dat_anti <- dat_anti[-c(1:3), ]
  ## only slect the position of head
  dat_anti <- as.data.frame(apply(dat_anti, 2, as.numeric))
  
  ## normalize the data, zero point value from fiji
  dat_anti$V2<- dat_anti$V2 - dat_anti$V8
  dat_anti$V3<- dat_anti$V3 - dat_anti$V9
  dat_anti$V5 <- dat_anti$V5 - dat_anti$V8
  dat_anti$V6 <- dat_anti$V6 - dat_anti$V9
  dat_anti <- dat_anti[,1:6]
  
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")
  
  ## remove the frame for the first few frame
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ]
  dat_anti <- subset(dat_anti, Prob >0.5 | Frame> 60)
  
  r_pix_mm <- 165 / length_pix
  dat_anti$Head_x <- abs(dat_anti$Head_x * r_pix_mm)
  dat_anti$Head_y <- abs(dat_anti$Head_y * r_pix_mm)
  
  dat_anti$Tail_x <- abs(dat_anti$Tail_x * r_pix_mm)
  dat_anti$Tail_y <- abs(dat_anti$Tail_y * r_pix_mm)
  
  chamber_div <- 165.5 ## the lenght of two chambels and the gap between
  
  ## assign each frame in hot as 1 (48) and nor plate(30) as 2
  dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
  return(dat_anti)
 
}

## for m1 day1
path_ctrl <- list(str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",1),
                  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",5),
                  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",7))
cc_plot <-function(x){
  p_heat <- ggplot(x, aes(x=Head_y, y=Head_x) ) +
    geom_path(color="indianred") + # darkcyan and indianred
    scale_fill_viridis()+ 
    theme_void()
  
  return(p_heat)
}
  
dat_heat_ctrl <- path_ctrl %>% 
  mapply(function(x) list.files(x, pattern = ".csv", full.names = T),., SIMPLIFY = F) %>% 
  mapply(function(x) x[1], ., SIMPLIFY = F) %>% 
  mapply(cc_behavior_reshape, ., frame_rate = 10, length_pix=364, SIMPLIFY = F) %>% 
  mapply(cc_plot,., SIMPLIFY = F) 

p_heat_ctrl <- plot_grid(dat_heat_ctrl[[1]], dat_heat_ctrl[[2]], dat_heat_ctrl[[3]], ncol = 1)
  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_ctrl.pdf", width = 90/25.6, height = 120/25.6, family = "Arial")
p_heat_ctrl
dev.off()

path_con <- list(str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",1),
                  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",5),
                  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",7))


dat_heat_con <- path_con %>% 
  mapply(function(x) list.files(x, pattern = ".csv", full.names = T),., SIMPLIFY = F) %>% 
  mapply(function(x) x[1], ., SIMPLIFY = F) %>% 
  mapply(cc_behavior_reshape, ., frame_rate = 10, length_pix=364, SIMPLIFY = F) %>% 
  mapply(cc_plot,., SIMPLIFY = F) 

p_heat_con <- plot_grid(dat_heat_con[[1]], dat_heat_con[[2]], dat_heat_con[[3]], ncol = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_con.pdf", width = 90/25.6, height = 120/25.6, family = "Arial")
p_heat_con
dev.off()

path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",i)

## for naloxone injected mice-----

dat_anti_ctrl <- vector(mode = "list", 2)
dat_anti_con <- vector(mode = "list", 2)
length_pix <- c(362,363)
day_group <- c("Pre",  "Test")
for(i in c(1:2)){
  d <- c(3,7)
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/naloxone_08042020/naloxone_ctrl/d",d[i])
  path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/naloxone_08042020/naloxone_mice/d",d[i])
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="Ctrl", ID=str_c("m",1:10), .before = 'ratio')
  
  dat_anti_con[[i]]<- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="Nalox.", ID=str_c("m",1:10), .before = 'ratio')
}

dat_anti_nalox <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  mutate(Day=factor(Day, levels = c("Pre",  "Test")), Group=factor(Group, levels = c("Ctrl", "Nalox.")))

dat_anti_nalox_sta <- ddply(dat_anti_nalox,.(variable,Day, Group), summarise, 'value'= mean(value))


p_ratio <- dat_anti_nalox %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(18, 17))+
  scale_color_manual(values=c("indianred", "mediumpurple4"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none', legend.title = element_blank())

p_ratio_test <- dat_anti_nalox %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

## for aov test
p_ratio_test <- dat_anti_nalox %>% 
  filter(variable=="ratio" ) %>% 
  aov(value ~ Group + Day, .)

p_latency_test <- dat_anti_nalox %>% 
  filter(variable=="latency2" ) %>% 
  aov(value ~ Group + Day, .)
summary(p_latency_test)

p_latency1 <- dat_anti_nalox %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(18, 17))+
  scale_color_manual(values=c("indianred", "mediumpurple4"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())

p_latency1_test <- dat_anti_nalox %>% 
  filter(variable=="latency1"& Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_nalox %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(18, 17))+
  scale_color_manual(values=c("indianred", "mediumpurple4"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150))+
  theme(legend.position = "none")

p_latency2_test <- dat_anti_nalox %>% 
  filter(variable=="latency2"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_anti_nalox %>% 
  filter(variable=="speed_compare") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "mediumpurple4"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_speed_compare_test <- dat_anti_nalox %>% 
  filter(variable=="speed_compare"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_total_cross <- dat_anti_nalox %>% 
  filter(variable=="total_cross"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.02,2.02,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_total_cross_test <- dat_anti_nalox %>% 
  filter(variable=="total_cross"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_pain_nalox <- plot_grid(p_ratio, p_latency2, ncol = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_nalox.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_pain_nalox
dev.off()

## manually sorted pain behavior
dat_nalox_pain <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/naloxone_08042020/wt_nalox_08042020_manualy.xlsx") %>% 
  mutate(First_licking = (First_licking - Frame_1st_crossing)/10, First_rearing = (First_rearing-Frame_1st_crossing)/10, 
         Guarding = (Guarding - Frame_1st_crossing)/10, Acceleration = (Acceleration- Frame_1st_crossing)/10, Jump = (Jump- Frame_1st_crossing)/10) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Nalox.")))

dat_nalox_pain_sta <- dat_nalox_pain %>% 
  ddply(.,.(variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_nalox_pain %>% 
  filter(variable=="First_licking") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "mediumpurple4"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 70), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(62,63,63,62),"path")+
  annotate("text",x=1.5,y=63, label="**", size=5)

p_licking_test <- dat_nalox_pain %>% 
  filter(variable=="First_licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_nalox_pain %>% 
  filter(variable=="First_rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "mediumpurple4"))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_1st_rearing_test <- dat_nalox_pain %>% 
  filter(variable=="First_rearing") %>% 
  t.test(value~Group,.)


p_jump <- dat_nalox_pain %>% 
  filter(variable=="Jump1") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "mediumpurple4"))+
  labs(x="", y="Latency to 1st jump (s) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(190,195,195,190),"path")+
  annotate("text",x=1.5,y=195, label="**", size=5)

p_jump_test<- dat_nalox_pain %>% 
  filter(variable=="Jump1") %>% 
  wilcox.test(value~Group,.)

p_rearing_freq <- dat_nalox_pain %>% 
  filter(variable=="Freq_rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "mediumpurple4"))+
  labs(x="", y="# rearings per mins ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(110,115,115,110),"path")+
  annotate("text",x=1.5,y=115, label="***", size=5)

p_Acceleration_test<- dat_nalox_pain %>% 
  filter(variable=="Acceleration") %>% 
  wilcox.test(value~Group,.)

p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_Acceleration, p_jump,nrow = 1)
p_con_nalox <- plot_grid(p_pain_compare1, p_pain_compare2, nrow = 2, rel_heights = c(2, 1))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_con_nalox.pdf", width = 140/25.6, height = 180/25.6, family = "Arial")
p_con_nalox
dev.off()


## for optical activation-----
dat_anti_ctrl <- vector(mode = "list", 2)
dat_anti_con <- vector(mode = "list", 2)
length_pix <- c(362,363)
day_group <- c("Pre",  "Test")
for(i in c(1:2)){
  d <- c(3,8)
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08042020/Ctrl/d",d[i])
  path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08042020/ChR2/d",d[i])
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="eYFP", ID=str_c("m1_",1:2), .before = 'ratio')
  
  dat_anti_con[[i]]<- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="ChR2", ID=str_c("m1_",3:5), .before = 'ratio')
}

dat_anti_optic <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  mutate(Day=factor(Day, levels = c("Pre",  "Test")), Group=factor(Group, levels = c("eYFP", "ChR2")))

dat_anti_optic_sta <- ddply(dat_anti_optic,.(variable,Day, Group), summarise, 'value'= mean(value))


p_ratio <- dat_anti_optic %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(Day, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none')

p_ratio_test <- dat_anti_optic %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_anti_optic %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150), expand = c(0, 0))+
  theme(legend.position = c(0.8, 0.8))

p_latency1_test <- dat_anti_optic %>% 
  filter(variable=="latency1"& Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_optic %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 90))+
  theme(legend.position = "none")

p_latency2_test <- dat_anti_optic %>% 
  filter(variable=="latency2"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_anti_optic %>% 
  filter(variable=="speed_compare") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_speed_compare_test <- dat_anti_optic %>% 
  filter(variable=="speed_compare"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_total_cross <- dat_anti_optic %>% 
  filter(variable=="total_cross"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.02,2.02,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_total_cross_test <- dat_anti_optic %>% 
  filter(variable=="total_cross"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_pain_compare1 <- plot_grid(p_ratio, p_latency1, p_latency2, p_speed_compare, ncol = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_opto.pdf", width = 140/25.6, height = 120/25.6, family = "Arial")
p_pain_compare1
dev.off()

## for optical inhibition-----
dat_anti_ctrl <- vector(mode = "list", 2)
dat_anti_con <- vector(mode = "list", 2)
length_pix <- c(362,363)
day_group <- c("Pre",  "Post")
for(i in c(1:2)){
  d <- c(3,7)
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/ctrl/d",d[i])
  path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/NpHR/d",d[i])
  
  
  dat_anti_ctrl1<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_extract(list.files(path_ctrl,pattern = ".csv", full.names = T ), regex("m\\d+")) )
  
  
  
  dat_anti_ctrl[[i]] <- dat_anti_ctrl1 %>% 
    add_column(Day = day_group[i], Group="eYFP", .before = 'ratio')
  
  dat_anti_con1 <- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_extract(list.files(path_con,pattern = ".csv", full.names = T ), regex("m\\d+")))
  
  
  dat_anti_con[[i]] <- dat_anti_con1 %>% 
    add_column(Day = day_group[i], Group="NpHR",  .before = 'ratio')
}

## combine data
dat_anti_optic <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  mutate( Day = factor(Day, levels = c("Pre", "Post")),Group=factor(Group, levels = c("eYFP", "NpHR")))

## only choose the mice which show the placebo anagesia
# ID_mouse_placebo <- dat_anti_optic %>% 
#   filter(variable == 'latency1') %>% 
#   filter(change <=0)
# ## filter the data by mouse ID selected
# dat_anti_optic_filter <- dat_anti_optic %>% 
#   filter(ID %in% ID_mouse_placebo$ID) %>% 
#   mutate( Group=factor(Group, levels = c("eYFP", "NpHR")))

dat_anti_optic_sta <- ddply(dat_anti_optic,.(variable,Day, Group), summarise, 'value'= mean(value))


p_ratio <- dat_anti_optic %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(19, 17))+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none')

p_ratio_test <- dat_anti_optic %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_anti_optic %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 180), expand = c(0, 0))+
  theme(legend.position = c(0.8, 0.8))

p_latency1_test <- dat_anti_optic %>% 
  filter(variable=="latency1"& Day=="Post") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_optic %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(19, 17))+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 120))+
  theme(legend.position = "none")

p_latency2_test <- dat_anti_optic %>% 
  filter(variable=="latency2"& Day=="Post") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_anti_optic %>% 
  filter(variable=="speed_compare") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_speed_compare_test <- dat_anti_optic %>% 
  filter(variable=="speed_compare"& Day=="Test") %>% 
  wilcox.test(value~Group,.)



p_pain_compare1 <- plot_grid(p_latency2,p_ratio, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_opto_inhibit.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_pain_compare1
dev.off()

## compare mannulay sorted behavior
dat_inhibit_manual <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/wt_optic_08252020_manualy.xlsx") %>% 
  mutate('Rearing'= (First_rearing - Frame_1st_crossing)/10, 'Licking'= (First_licking -Frame_1st_crossing)/10 , 'Jumping'= (Jump -Frame_1st_crossing)/10) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  filter(Group != "ChR2") %>% 
  mutate(Group=factor(Group, levels=c("eYFP", "NpHR")))

dat_inhibit_manual_sta <- dat_inhibit_manual %>% 
  ddply(., .( variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_inhibit_manual %>% 
  filter(variable=="Licking") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17))+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 60), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(55,56,56,55),"path")+
  annotate("text",x=1.5,y=56, label="**", size=5)

p_licking_test <- dat_inhibit_manual %>% 
  filter(variable=="Licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_inhibit_manual %>% 
  filter(variable=="Rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17))+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_1st_rearing_test <- dat_inhibit_manual %>% 
  filter(variable=="Rearing") %>% 
  wilcox.test(value~Group,.)

p_Jump <- dat_inhibit_manual %>% 
  filter(variable=="Jump1") %>% 
  ggplot(., aes(Group, value/10, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17))+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to 1st jumping ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(185,187,187,185),"path")+
  annotate("text",x=1.5,y=187, label="***", size=5)

p_jumping_test<- dat_inhibit_manual %>% 
  filter(variable=="Jump1") %>% 
  wilcox.test(value~Group,.)


p_num_rearing <- dat_inhibit_manual %>% 
  filter(variable=="num_rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17))+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="# rearing / min ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_num_rearing_test<- dat_inhibit_manual %>% 
  filter(variable=="num_rearing") %>% 
  wilcox.test(value~Group,.)

p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_num_rearing,p_Jump, nrow = 1, rel_widths = c(1.05, 1.1, 1, 1.1))
p_con <- plot_grid(p_pain_compare1, p_pain_compare2, nrow = 2, rel_heights = c(2, 1))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_optic_inhibit.pdf", width = 132/25.6, height = 60/25.6, family = "Arial")
p_pain_compare2
dev.off()

## for optical excitation anti behavior------
dat_anti_ctrl <- vector(mode = "list", 2)
dat_anti_ChR2 <- vector(mode = "list", 2)
length_pix <- c(362,363)
day_group <- c("Pre",  "Post")
for(i in c(1:2)){
  d <- c(3,7)
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/ctrl/d",d[i])
  path_ChR2 <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/ChR2/d",d[i])
  
  
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_c('C',str_extract(list.files(path_ctrl,pattern = ".csv", full.names = T ), regex("m\\d+")) )) %>% 
    add_column(Day = day_group[i], Group="eYFP", .before = 'ratio')
  
  
  
  dat_anti_ChR2[[i]]<- list.files(path_ChR2,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_c('C',str_extract(list.files(path_ChR2,pattern = ".csv", full.names = T ), regex("m\\d+")) )) %>% 
    add_column(Day = day_group[i], Group="ChR2",  .before = 'ratio')
  
}

## combine data
dat_anti_optic <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_ChR2)) %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  mutate( Day = factor(Day, levels = c("Pre", "Post")),Group=factor(Group, levels = c("eYFP", "ChR2")))

## only choose the mice which show the placebo anagesia
# ID_mouse_placebo <- dat_anti_optic %>% 
#   filter(variable == 'latency1') %>% 
#   filter(change <=0)
# ## filter the data by mouse ID selected
# dat_anti_optic_filter <- dat_anti_optic %>% 
#   filter(ID %in% ID_mouse_placebo$ID) %>% 
#   mutate( Group=factor(Group, levels = c("eYFP", "NpHR")))

dat_anti_optic_sta <- ddply(dat_anti_optic,.(variable,Day, Group), summarise, 'value'= mean(value))


p_ratio <- dat_anti_optic %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 15))+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none')

p_ratio_test <- dat_anti_optic %>% 
  filter(variable=="ratio" & Day=="Post") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_anti_optic %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 15))+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 180), expand = c(0, 0))+
  theme(legend.position = 'none')

p_latency1_test <- dat_anti_optic %>% 
  filter(variable=="latency1"& Day=="Post") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_optic %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 15))+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

p_latency2_test <- dat_anti_optic %>% 
  filter(variable=="latency2"& Day=="Post") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_anti_optic %>% 
  filter(variable=="speed_compare") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 15))+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_speed_compare_test <- dat_anti_optic %>% 
  filter(variable=="speed_compare"& Day=="Post") %>% 
  wilcox.test(value~Group,.)



p_pain_compare1 <- plot_grid( p_latency1, p_latency2,p_ratio, p_speed_compare, ncol = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_opto_inhibit.pdf", width = 140/25.6, height = 120/25.6, family = "Arial")
p_pain_compare1
dev.off()

## compare mannulay sorted behavior
dat_excite_manual <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/wt_optic_08252020_manualy.xlsx") %>% 
  mutate('Rearing'= (First_rearing - Frame_1st_crossing)/10, 'Licking'= (First_licking -Frame_1st_crossing)/10 , 'Jumping'= (Jump -Frame_1st_crossing)/10) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  filter(Group != "NpHR") %>% 
  mutate(Group=factor(Group, levels=c("eYFP", "ChR2")))

dat_excite_manual_sta <- dat_excite_manual %>% 
  ddply(., .( variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_excite_manual %>% 
  filter(variable=="Licking") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  scale_shape_manual(values=c(20, 15))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")
  # annotate(x=c(1,1,2,2), y=c(55,56,56,55),"path")+
  # annotate("text",x=1.5,y=56, label="**", size=5)

p_licking_test <- dat_excite_manual %>% 
  filter(variable=="Licking") %>% 
  wilcox.test(value~Group, .) 

p_1st_rearing <- dat_excite_manual %>% 
  filter(variable=="Rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  scale_shape_manual(values=c(20, 15))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_1st_rearing_test <- dat_excite_manual %>% 
  filter(variable=="Rearing") %>% 
  wilcox.test(value~Group,.)

p_Jump <- dat_excite_manual %>% 
  filter(variable=="Jump1") %>% 
  ggplot(., aes(Group, value/10, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  scale_shape_manual(values=c(20, 15))+
  labs(x="", y="Latency to 1st jumping ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_jumping_test<- dat_excite_manual %>% 
  filter(variable=="Jump1") %>% 
  wilcox.test(value~Group,.)


p_num_rearing <- dat_excite_manual %>% 
  filter(variable=="num_rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("springgreen4", "royalblue4"))+
  scale_shape_manual(values=c(20, 15))+
  labs(x="", y="# rearing / min ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_num_rearing_test<- dat_excite_manual %>% 
  filter(variable=="num_rearing") %>% 
  wilcox.test(value~Group,.)

p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_num_rearing,p_Jump, nrow = 1)
p_con <- plot_grid(p_pain_compare1, p_pain_compare2, nrow = 2, rel_heights = c(2, 1))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_optic_excite.pdf", width = 140/25.6, height = 180/25.6, family = "Arial")
p_con
dev.off()
## for pn cpp experiment data analysis (inhibiiton)-----
cc_anti_cpp <- function(path_ID, frame_rate){
  
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   ## only select the position of head
    as.data.frame() %>% 
    select(str_c("V", c(1:4,8:9))) 
  
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob",  "reference_x", "reference_y")
  
  ## remove the frame for the first few frame
  ## the length and width of the CPP box (426.6,416 )
  limit_x <- mean(dat_anti$reference_x)
  limit_y <- mean(dat_anti$reference_y)

  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    subset(., Prob >0.5 | Frame> 60) %>% 
    mutate(Head_x = ifelse(Head_x < limit_x , limit_x, Head_x), Head_y = ifelse(Head_y < (limit_y - 416), (limit_y-416), Head_y)) %>% 
    mutate(Head_x = ifelse(Head_x > limit_x + 426.6 , limit_x +426.6, Head_x), Head_y = ifelse(Head_y > (limit_y + 430), (limit_y+430), Head_y)) %>% 
    mutate(reference_y = mean(reference_y)) %>% 
    mutate(chamber = ifelse(Head_y > (reference_y +10), "right", "left"))
  
  ## for the stimulation side
  test_d <- str_extract(path_ID, regex("d\\d+"))
  m_ID <- str_extract(path_ID, regex("m\\d+"))
  
  if (test_d == "d1"){

    ratio_time <- dat_anti %>% 
      count(chamber) %>% 
      mutate(ratio = n/sum(n)*100, ID = m_ID) %>% 
      select(!n)
    
  } else{
    stim_side <- names(sort(table(dat_anti$chamber[1:10]),decreasing=TRUE)[1])
    
    ratio_time <- dat_anti %>% 
      count(chamber) %>% 
      mutate(ratio1 = n/sum(n)*100) %>% 
      mutate(stim = stim_side, ID1 = m_ID) %>% 
      rename(chamber1 = chamber)
  }


  p_anti <- ggplot(dat_anti, aes(Head_y, Head_x,colour=chamber))+
    geom_point(shape=1, size=1)+
    theme_void()+
    theme(legend.position = "none")
  
  return(ratio_time)
}



## For inhibition
path_cpp_inhibit_pre <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_cpp/inhibition/d1",pattern = ".csv", full.names = T )
cpp_inhibit_pre <- mapply(cc_anti_cpp, path_cpp_inhibit_pre, frame_rate = 10, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  remove_rownames() %>% 
  mutate(Group = "NpHR")

path_cpp_inhibit_test <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_cpp/inhibition/d4",pattern = ".csv", full.names = T )
cpp_inhibit_test <- mapply(cc_anti_cpp, path_cpp_inhibit_test, frame_rate = 10, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  remove_rownames() 

dat_cpp_inhibt <- cbind(cpp_inhibit_pre, cpp_inhibit_test) %>% 
  filter(chamber1 == stim) %>% 
  select(ID, ratio, ratio1, Group)
  


## for eYFP
path_cpp_eYFP_pre <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_cpp/EGFP/d1",pattern = ".csv", full.names = T )
cpp_eYFP_pre <- mapply(cc_anti_cpp, path_cpp_eYFP_pre, frame_rate = 10, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  remove_rownames() %>% 
  mutate(Group = "eYFP")

path_cpp_eYFP_test <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_cpp/EGFP/d4",pattern = ".csv", full.names = T )
cpp_eYFP_test <- mapply(cc_anti_cpp, path_cpp_eYFP_test, frame_rate = 10, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  remove_rownames() 

dat_cpp_eYFP <- cbind(cpp_eYFP_pre, cpp_eYFP_test) %>% 
  filter(chamber1 == stim) %>% 
  select(ID, ratio, ratio1, Group)


## For excitation
path_cpp_excite_pre <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_cpp/excitation/d1",pattern = ".csv", full.names = T )
cpp_excite_pre <- mapply(cc_anti_cpp, path_cpp_excite_pre, frame_rate = 10, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  remove_rownames() %>% 
  mutate(Group = "ChR2") 

path_cpp_excite_test <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_cpp/excitation/d4",pattern = ".csv", full.names = T )
cpp_excite_test <- mapply(cc_anti_cpp, path_cpp_excite_test, frame_rate = 10, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  remove_rownames()

dat_cpp_excite <- cbind(cpp_excite_pre, cpp_excite_test) %>% 
  filter(chamber1 == stim) %>% 
  select(ID, ratio, ratio1, Group) %>% 
  mutate(ID = str_c(ID,"e"))


dat_cpp <- rbind(dat_cpp_eYFP, dat_cpp_inhibt, dat_cpp_excite) %>% 
  mutate(change = (ratio1-ratio)/ratio*100) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2")))

## change the value based on stimulation side
dat_cpp$change[dat_cpp$ID=="m41"] <- -(dat_cpp$change[dat_cpp$ID=="m41"])
dat_cpp$change[dat_cpp$ID=="m44"] <- -(dat_cpp$change[dat_cpp$ID=="m44"])
dat_cpp$change[dat_cpp$ID=="m41e"] <- -(dat_cpp$change[dat_cpp$ID=="m41e"])
dat_cpp$change[dat_cpp$ID=="m43e"] <- -(dat_cpp$change[dat_cpp$ID=="m43e"])

dat_cpp_sta <- dat_cpp %>% 
  ddply(., .(Group), summarise,n=length(change),mean=mean(change),sd=sd(change),se=sd(change)/sqrt(length(change)))
  

p_cpp <- dat_cpp %>% 
  ggplot(., aes(Group, change, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_color_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="CPP score (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = 'none')

t_cpp <- dat_cpp %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', 'ChR2'))) %>% 
  aov(change~Group,.)
summary(t_cpp)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cpp_combin.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_cpp
dev.off()

## measuring the thermal and mechanical pain with optic inhibition------
p_von_inhibit <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  group_by(Group, ID) %>% 
  # summarize(d=Threshold[Light == "on"]-Threshold[Light == "off"]) %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', 'ChR2'))) %>% 
  mutate(Light = factor(Light, levels = c("off", "on"))) %>% 
  ggplot(., aes(interaction(Light, Group),Threshold, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="pain threshold (g)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 6))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_von_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  group_by(Group, ID) %>% 
  summarize(d=Threshold[Light == "on"]-Threshold[Light == "off"]) %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', 'ChR2'))) %>% 
  ddply(., .(Group), summarise,n=length(d),mean=mean(d),sd=sd(d),se=sd(d)/sqrt(length(d)))
  
  
  
t_von <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  group_by(Group, ID) %>% 
  #summarize(d=Threshold[Light == "on"]-Threshold[Light == "off"]) %>% 
  #mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', 'ChR2'))) %>% 
  aov(Threshold~Group*Light,.)
summary(t_von)
TukeyHSD(t_von, ordered = T)  

p_har_inhibit <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = "Har_optic_inhibition") %>% 
  select("Group", "ID","Light", "Latency") %>%
  group_by(Group, ID) %>% 
  summarize(d=Latency[Light == "on"]-Latency[Light == "off"]) %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR'))) %>% 
  ggplot(., aes(Group, d, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.25,shape = 1)+
  scale_colour_manual(values=c("springgreen4", "yellow3" ))+
  labs(x="", y="in withdraw latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_inhibit_combin <- plot_grid(p_von_inhibit, p_har_inhibit, ncol = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_von_inhibit.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_von_inhibit
dev.off()

p_hot_inhibit1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = 4) %>% 
  filter(Temp == 48) %>% 
  select("Group", "ID", "Latency_withdraw") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  ggplot(., aes(Group, Latency_withdraw, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Withdrawal latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_hot1_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Pn_anti_10072020.xlsx", sheet = 4) %>% 
  filter(Temp == 48) %>% 
  select("Group", "ID", "Latency_withdraw") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  ddply(., .(Group), summarise,n=length((Latency_withdraw)),mean=mean((Latency_withdraw)),sd=sd((Latency_withdraw)),se=sd((Latency_withdraw))/sqrt(length((Latency_withdraw))))

t_hot1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Pn_anti_10072020.xlsx", sheet = 4) %>% 
  filter(Temp == 48) %>% 
  select("Group", "ID", "Latency_withdraw") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  aov(Latency_withdraw~Group, .)
summary(t_hot1)

p_hot_inhibit2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Pn_anti_10072020.xlsx", sheet = 4) %>% 
  filter(Temp == 52) %>% 
  select("Group", "ID", "Latency_withdraw") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  ggplot(., aes(Group, Latency_withdraw, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Withdrawal latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_hot2_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Pn_anti_10072020.xlsx", sheet = 4) %>% 
  filter(Temp == 52) %>% 
  select("Group", "ID", "Latency_withdraw") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  ddply(., .(Group), summarise,n=length((Latency_withdraw)),mean=mean((Latency_withdraw)),sd=sd((Latency_withdraw)),se=sd((Latency_withdraw))/sqrt(length((Latency_withdraw))))

t_hot2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Pn_anti_10072020.xlsx", sheet = 4) %>% 
  filter(Temp == 52) %>% 
  select("Group", "ID", "Latency_withdraw") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  aov(Latency_withdraw~Group, .)
summary(t_hot2)

p_hot_inhibit <- plot_grid(p_hot_inhibit1, p_hot_inhibit2, nrow = 1)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_hot_inhibit.pdf", width = 82/25.6, height = 60/25.6, family = "Arial")
p_hot_inhibit
dev.off()

## for von frey with multiple filaments

p_von_mutiple <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = "von_frey") %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  ddply(., .(Group, Num_von), summarise,n = length(Change),mean=mean(Change),sd=sd(Change),se=sd(Change)/sqrt(length(Change))) %>% 
  ggplot(., aes(Num_von, mean, color=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  geom_point(aes(colour = Group, shape = Group),size=2)+
  geom_line()+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Withdraw latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_continuous(trans='log10')+
  scale_y_continuous(limits = c(-0.8, 1))+
  geom_hline(yintercept = 0, linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.8))


p_von_mutiple1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = "von_frey") %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR', "ChR2"))) %>% 
  ddply(., .(Group, Num_von), summarise,n = length(Withdraw_on),mean=mean(Withdraw_on),sd=sd(Withdraw_on),se=sd(Withdraw_on)/sqrt(length(Withdraw_on))) %>% 
  ggplot(., aes(Num_von, mean, color=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  geom_point(aes(colour = Group, shape = Group),size=2)+
  geom_line()+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Withdraw latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_continuous(trans='log10')+
  scale_y_continuous(limits = c(0, 1.2))+
  geom_hline(yintercept = 0, linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

library(tidyverse)
library(ggpubr)
library(rstatix)

t_von_mutiple <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = "von_frey") %>% 
  mutate(Num_von = as.factor(Num_von)) %>% 
  aov(Change ~ Group + Num_von + Group:Num_von, .)
summary(t_von_mutiple )


TukeyHSD(t_von_mutiple)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_von_multiple1.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_von_mutiple1
dev.off()

t_von_mutiple <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = "von_frey") %>% 
  mutate(Num_von = as.factor(Num_von)) %>% 
  aov(Change ~ Group + Num_von + Group:Num_von, .)
summary(t_von_mutiple )

TukeyHSD(t_von_mutiple, which = "Group:Num_von")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_von_multiple.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_von_mutiple
dev.off()
## for the rotarod experiment -----

p_rotarod <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Rotarod_11092020.xlsx", sheet = 1) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2"))) %>% 
  ggplot(., aes(Group, Diff, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_color_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="D latency of fall (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = 'none')

t_rotarod <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Rotarod_11092020.xlsx", sheet = 1) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2"))) %>% 
  aov(Diff~Group, .)
summary(t_rotarod)  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_rotarod.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_rotarod
dev.off()


## for penk opto experiment-----
dat_anti_ctrl <- vector(mode = "list", 2)
dat_anti_ChR2 <- vector(mode = "list", 2)
dat_anti_NpHR <- vector(mode = "list", 2)
length_pix <- c(362,363)
day_group <- c("Pre",  "Test")
for(i in c(1:2)){
  d <- c(3,7)
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/Penk_opto/ctrl/d",d[i])
  path_ChR2 <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/Penk_opto/ChR2/d",d[i])
  path_NpHR <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/Penk_opto/NpHR/d",d[i])
  
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_extract(list.files(path_ctrl, pattern = ".csv", full.names = T ), regex("m\\d+"))) %>% 
    add_column(Day = day_group[i], Group="eYFP", .before = 'ratio')
  
  dat_anti_ChR2[[i]]<- list.files(path_ChR2,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_extract(list.files(path_ChR2, pattern = ".csv", full.names = T ), regex("m\\d+"))) %>% 
    add_column(Day = day_group[i], Group="ChR2", .before = 'ratio')
  
  dat_anti_NpHR[[i]]<- list.files(path_NpHR,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    mutate(ID = str_extract(list.files(path_NpHR, pattern = ".csv", full.names = T ), regex("m\\d+"))) %>% 
    add_column(Day = day_group[i], Group="NpHR", .before = 'ratio')
}

dat_penk_optic <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_ChR2)) %>% 
  rbind(., do.call(rbind, dat_anti_NpHR)) %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  mutate(Day=factor(Day, levels = c("Pre",  "Test")), Group=factor(Group, levels = c("eYFP", "NpHR","ChR2")))

dat_penk_optic_sta <- ddply(dat_penk_optic,.(variable,Day, Group), summarise, 'value'= mean(value))


p_ratio <- dat_penk_optic %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(Day, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none')

p_ratio_test <- dat_penk_optic %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_penk_optic %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150), expand = c(0, 0))+
  theme(legend.position = c(0.8, 0.8))

p_latency1_test <- dat_penk_optic %>% 
  filter(variable=="latency1"& Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_penk_optic %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 90))+
  theme(legend.position = "none")

p_latency2_test <- dat_penk_optic %>% 
  filter(variable=="latency2"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_penk_optic %>% 
  filter(variable=="speed_compare") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_speed_compare_test <- dat_anti_optic %>% 
  filter(variable=="speed_compare"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_total_cross <- dat_penk_optic %>% 
  filter(variable=="total_cross"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.02,2.02,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_total_cross_test <- dat_anti_optic %>% 
  filter(variable=="total_cross"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_pain_compare1 <- plot_grid(p_ratio, p_latency1, p_latency2, p_speed_compare, ncol = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_opto.pdf", width = 140/25.6, height = 120/25.6, family = "Arial")
p_pain_compare1
dev.off()
## mannual sorted pain behavior--
dat_penk_manual <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/wt_optic_08252020_manualy.xlsx", sheet = 2) %>% 
  mutate('Rearing'= (First_rearing - Frame_1st_crossing)/10, 'Licking'= (First_licking -Frame_1st_crossing)/10 , 'Jumping'= (Jump -Frame_1st_crossing)/10) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  mutate(Group=factor(Group, levels=c("EYFP", "Nphr", "Chr2")))

dat_penk_manual_sta <- dat_penk_manual %>% 
  ddply(., .( variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_penk_manual %>% 
  filter(variable=="Licking") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 60), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(55,56,56,55),"path")+
  annotate("text",x=1.5,y=56, label="**", size=5)

p_licking_test <- dat_penk_manual %>% 
  filter(variable=="Licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_penk_manual %>% 
  filter(variable=="Rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_1st_rearing_test <- dat_penk_manual %>% 
  filter(variable=="Rearing") %>% 
  wilcox.test(value~Group,.)

p_Jump <- dat_penk_manual %>% 
  filter(variable=="Jumping") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to 1st jumping ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(185,187,187,185),"path")+
  annotate("text",x=1.5,y=187, label="**", size=5)

p_jumping_test<- dat_inhibit_manual %>% 
  filter(variable=="Jumping") %>% 
  wilcox.test(value~Group,.)


p_num_rearing <- dat_penk_manual %>% 
  filter(variable=="num_rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="# rearing / min ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_num_rearing_test<- dat_inhibit_manual %>% 
  filter(variable=="num_rearing") %>% 
  wilcox.test(value~Group,.)

p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_num_rearing,p_Jump, nrow = 1)
p_con <- plot_grid(p_pain_compare1, p_pain_compare2, nrow = 2, rel_heights = c(2, 1))

## for penk opto ethovision------
p_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 2) %>% 
  select(Group, Day, ID,  Ratio) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2")), Day = factor(Day, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(interaction(Day, Group), Ratio, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")

p_latency1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 2) %>% 
  select(Group, Day, ID,  Latency1) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2")), Day = factor(Day, levels = c("Pre", "Post"))) %>%
  ggplot(., aes(interaction(Day, Group), Latency1, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Latency of crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")


p_latency2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 2) %>% 
  select(Group, Day, ID,  Latency2) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2")), Day = factor(Day, levels = c("Pre", "Post"))) %>%
  ggplot(., aes(interaction(Day, Group), Latency2, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_pain_compare1 <- plot_grid(p_latency2, p_ratio, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_compare_penk.pdf", width = 150/25.6, height = 70/25.6, family = "Arial")
p_pain_compare1
dev.off()

p_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 1) %>% 
  select(Group, ID,  First_licking) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2"))) %>%
  ggplot(., aes(Group, First_licking, group= Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")


p_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 1) %>% 
  select(Group, ID,  First_rearing) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2"))) %>%
  ggplot(., aes(Group, First_rearing, group= Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_num_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 1) %>% 
  select(Group, ID,  Num_rearing) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2"))) %>%
  ggplot(., aes(Group, Num_rearing, group= Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="# of rearing")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_jump <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pn_penk_behavior_ethovision/PN_penk_pn_behavior.xlsx", sheet = 1) %>% 
  select(Group, ID,  First_jumping) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR", "ChR2"))) %>%
  ggplot(., aes(Group, First_jumping, group= Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("springgreen4", "orange2", "dodgerblue4"))+
  labs(x="", y="Latency of jumping")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_pain_compare2 <- plot_grid(p_licking, p_rearing, p_num_rearing,p_jump, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_compare2_penk.pdf", width = 200/25.6, height = 60/25.6, family = "Arial")
p_pain_compare2
dev.off()

## for Dor-cre mice with optical inhibition-------
dat_dor_opto <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_anti_10072020.xlsx", sheet = "Dor_optic_anti") %>% 
  as_tibble() %>% 
  mutate(First_licking = First_licking - Latency_1st_crossing, First_rearing = First_rearing -Latency_1st_crossing, Jumping = Jumping - Latency_1st_crossing, Crossing_back = Crossing_back - Latency_1st_crossing ) %>%
  mutate(cross_diff = Crossing_back - Crossing_back_d3) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "NpHR")))

dat_dor_opto_sta <- dat_dor_opto %>% 
  pivot_longer(-c(Group, ID), names_to = "Variable") %>% 
  ddply(., .( Variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) 

p_dor_licking <-  dat_dor_opto %>% 
  ggplot(., aes(Group, First_licking, Group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ 
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(19, 17))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(10, 50))+
  theme(legend.position = "none")


p_dor_rearing <-  dat_dor_opto %>% 
  ggplot(., aes(Group, First_rearing, Group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ 
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(19, 17))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(20, 120))+
  theme(legend.position = "none")

p_dor_jumping <-  dat_dor_opto %>% 
  ggplot(., aes(Group, Jumping, Group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ 
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(19, 17))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(80, 200))+
  theme(legend.position = "none")

p_dor_crossback <-  dat_dor_opto %>% 
  ggplot(., aes(Group, cross_diff, Group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+ 
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(19, 17))+
  scale_colour_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")



wilcox.test(First_licking~Group, dat_dor_opto)
wilcox.test(First_rearing~Group, dat_dor_opto)
wilcox.test(Jumping~Group, dat_dor_opto)
wilcox.test(Crossing_back~Group, dat_dor_opto)

p_dor_anti <- plot_grid(p_dor_licking, p_dor_rearing, p_dor_jumping, ncol = 3)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dor_anti.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_dor_anti
dev.off()

## open field and rota rod for Dor-cre mice------
dat_dor_openfield <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/open_field_dor.xlsx") %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR")))

p_distance <- dat_dor_openfield %>% 
  filter(variable == "Distance") %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line()+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="Time (min)", y="Distance (cm)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  annotate("rect", xmin = 2, xmax = 4, ymin = 100, ymax = 350 ,alpha = .1,fill = "yellow")+
  annotate("rect", xmin = 6, xmax = 8, ymin = 100, ymax = 350 ,alpha = .1,fill = "yellow")+
  theme(legend.title = element_blank())

 

p_velocity <- dat_dor_openfield %>% 
  filter(variable == "Velocity") %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line()+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="Time (min)", y="Velocity (cm/s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  annotate("rect", xmin = 2, xmax = 4, ymin =3, ymax = 12 ,alpha = .1,fill = "yellow")+
  annotate("rect", xmin = 6, xmax = 8, ymin = 3, ymax = 12 ,alpha = .1,fill = "yellow")+
  theme(legend.title = element_blank())

t_dat_dor_openfield <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/open_field_dor.xlsx") %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  filter(variable == "Velocity") %>% 
  mutate(Light = rep(rep(c("Off", "On"), each = 8), 2)) %>% 
  aov(mean~ Light * Group ,.)

summary(t_dat_dor_openfield)
  

p_freq <- dat_dor_openfield %>% 
  filter(variable == "Freq") %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line()+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="Time (min)", y="Freq")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  annotate("rect", xmin = 2, xmax = 4, ymin =2, ymax = 12 ,alpha = .1,fill = "yellow")+
  annotate("rect", xmin = 6, xmax = 8, ymin = 2, ymax = 12 ,alpha = .1,fill = "yellow")+
  theme(legend.title = element_blank())

p_duration <- dat_dor_openfield %>% 
  filter(variable == "Duration") %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line()+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="Time (min)", y="Duration in center (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  annotate("rect", xmin = 2, xmax = 4, ymin =3, ymax = 16 ,alpha = .1,fill = "yellow")+
  annotate("rect", xmin = 6, xmax = 8, ymin = 3, ymax = 16 ,alpha = .1,fill = "yellow")+
  theme(legend.title = element_blank())

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_velocity.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
p_velocity
dev.off()

## for Rota rod
p_rota_inhibit <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 3) %>% 
  group_by(Group, ID) %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR'))) %>% 
  mutate(Light = factor(Light, levels = c("Off", "On"))) %>% 
  ggplot(., aes(interaction(Light, Group),Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="Latency to fall (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

t_rota_inhibit <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 3) %>% 
  group_by(Group, ID) %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR'))) %>% 
  mutate(Light = factor(Light, levels = c("Off", "On"))) %>% 
  aov(Latency~ Group *Light,.)

summary(t_rota_inhibit)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dor_rota.pdf", width = 65/25.6, height = 60/25.6, family = "Arial")
p_rota_inhibit
dev.off()

## Von frey for Dor cells in Pn------
p_von_inhibit <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  group_by(Group, ID) %>% 
  # summarize(d=Threshold[Light == "on"]-Threshold[Light == "off"]) %>% 
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR'))) %>% 
  mutate(Light = factor(Light, levels = c("off", "on"))) %>% 
  ggplot(., aes(interaction(Light, Group),Threshold, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="pain threshold (g)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_von_inhibit_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  group_by(Group, ID) %>% 
  aov(Threshold~Group*Light,.) %>% 
  summary
  
t_von_dor <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  filter(Light == "on") %>% 
  wilcox.test(Threshold~Group,., paired = F)

t_von_dor1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 1) %>% 
  select("Group", "ID","Light", "Threshold") %>%
  filter(Group == "NpHR") %>% 
  wilcox.test(Threshold~Light,., paired = T)

p_hot_inhibit1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 2) %>% 
  filter(Temp == 48) %>% 
  select("Group", "ID", "Latency") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR'))) %>% 
  ggplot(., aes(Group, Latency, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="Withdrawal latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

t_hot_inhibit1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 2) %>% 
  filter(Temp == 48) %>% 
  wilcox.test(Latency~Group,.)
  

p_hot_inhibit2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 2) %>% 
  filter(Temp == 52) %>% 
  select("Group", "ID", "Latency") %>%
  mutate(Group = factor(Group, levels = c('eYFP', 'NpHR'))) %>% 
  ggplot(., aes(Group, Latency, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_shape_manual(values=c(19, 17, 15))+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="Withdrawal latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_hot <- plot_grid(p_hot_inhibit1, p_hot_inhibit2, ncol = 2)
p_dor_pain <- plot_grid(p_von_inhibit, p_hot, ncol = 2)

dat_hot_inhibit1_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 2) %>% 
  ddply(., .(Temp, Group), summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency))) 
  

t_hot_inhibit1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pn_dor/Pn_dor_vonfrey.xlsx", sheet = 2) %>% 
  filter(Temp == 52) %>% 
  wilcox.test(Latency~Group,.)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dor_pain_von.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_von_inhibit
dev.off()

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dor_pain_hot.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_hot
dev.off()

## Dor cre hm4di-----
p_von_inhibit <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Von") %>% 
  ddply(., .(Group, Time), summarise,n=length(Threshold),mean=mean(Threshold),sd=sd(Threshold),se=sd(Threshold)/sqrt(length(Threshold))) %>%
  ggplot(., aes(Time, mean, colour=Group))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  scale_colour_manual(values=c("turquoise4", "yellow4"))+
  labs(x="", y="pain threshold (g)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")


