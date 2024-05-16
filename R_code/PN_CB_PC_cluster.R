## CLustering analysis for PC data------
setwd('/Users/cchen/Documents/neuroscience/Pn project/Pn_project/')
source("Functions/c_miniscope_matlab_d7.R")

back_crossing_frame <- tibble(ID = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame = c(1974, 2588, 1612, 1294, 2210, 2698), 
                              Last_frame= c(3016, 3168, 3106, 3180, 3726, 3640))

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_d7, mouse_file, SIMPLIFY = F)

## do the clustering analysis

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)

fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")

dat_cell_trace_cluster1 <- kmeans(t(dat_cb_combine), centers = 2, iter.max = 10)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cb_cluster.pdf", width = 80/25.6, height = 65/25.6, family = "Arial")
fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")
dev.off()
fviz_cluster(dat_cell_trace_cluster1, t(dat_cb_combine), geom = "point", ellipse.type = "convex", ggtheme = theme_bw())

set.seed(42)
dat_cell_trace_cluster <- kmeans(t(dat_cb_combine), centers = 2, iter.max = 10)[[1]]


# Heatmap plot
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Crossing", "Crossing_back", "Last_crossing")
score_range <- range(dat_cell_trace)



for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
    group_by(name, cluster) %>%
    summarise(mean = mean(value), .groups = "drop") %>%
    arrange(cluster, mean)
  
  dat_rect <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
    ddply(., .(cluster), summarise, n = length(name)/40) %>% 
    add_column(xmin = 2, xmax = 2.1) %>% 
    mutate(ymax = cumsum(n) ) %>% 
    mutate(ymin = c(1, lag(ymax)[-1]))
  
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(name = factor(name, levels = cell_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value))+ 
    geom_tile(height=2)+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="", y="")+
    geom_vline(xintercept = 0, col= 'red', linetype=2)+
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
    theme(legend.position = "none")+
    annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
  
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Crossing, p_heat_Crossing_back, p_heat_Last_crossing, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_com_cluster.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()    


# box plot
p_cell_firing_corssing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Crossing")

p_cell_firing_corssing_back <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Crossing_back")

p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  ggplot(., aes(Group, mean, group = Group, color = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=name), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 6))+
  theme(legend.position = 'none')

t_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  filter(cluster ==2) %>% 
  aov(value~Group,.)

dat_cell_firing_sta <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  filter(cluster ==1) %>% 
  ddply(.,.(Group), summarise, mean= mean(value))

summary(t_cell_firing)
TukeyHSD(t_cell_firing)

## for individual mice

p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ddply(.,.(name, Group, cluster, ID), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.( Group, cluster, ID), summarise, mean = mean(mean_acti)) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.5, 3))+
  theme(legend.position = 'none')
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_firing.pdf", width = 60/25.6, height = 65/25.6, family = "Arial" )
p_cell_firing
dev.off()

## compare PC in two clusters during Pre, Cond, and Post-----
rm(list = setdiff(ls(), "dat_cell_trace_cluster"))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    
    if (t1_p < 40){
      column_means = colMeans(dat_trace[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_trace))
      
      # Combine the new matrix with the original one
      dat_trace = rbind(new_rows, dat_trace)
      
    } else {
      dat_trace <- dat_trace
    }
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_trace[(t1_p-39):(t1_p+40),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)

comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")
score_range <- range(dat_cell_trace)

for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
    group_by(name, cluster) %>%
    summarise(mean = mean(value), .groups = "drop") %>%
    arrange(cluster, mean)
  
  dat_rect <- dat_cell_cond_trace %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
    ddply(., .(cluster), summarise, n = length(name)/40) %>% 
    add_column(xmin = 2, xmax = 2.1) %>% 
    mutate(ymax = cumsum(n) ) %>% 
    mutate(ymin = c(1, lag(ymax)[-1]))
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(name = factor(name, levels = cell_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value))+ 
    geom_tile(height=2)+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="", y="")+
    geom_vline(xintercept = 0, col= 'red', linetype=2)+
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
    theme(legend.position = "none")+
    annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
  
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_pc.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()


## boxplot 

p_combine_cluster <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=name), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="Time (s)", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 8))+
  theme(legend.position = 'none')

t_combine_cluster <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==2) %>% 
  #filter(Group != "Post") %>% 
  #wilcox.test(mean~Group,., paired = T)
  aov(mean~Group, .)

dat_cell_trace_cluster1 <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==1) %>% 
  ddply(., .(Group), summarize, mean_acti = mean(mean), n = length(mean),se=sd(mean)/sqrt(length(mean))) 



summary(t_combine_cluster)
TukeyHSD(t_combine_cluster)
## plot for each mouse
p_combine_cluster_mouse <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  ddply(., .(ID,Group, cluster ), summarise, mean_acti = mean(mean)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, mean_acti, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 4))+
  theme(legend.position = 'none')


t_combine_cluster_mouse <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  ddply(., .(ID,Group, cluster ), summarise, mean_acti = mean(mean)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==2) %>% 
  # filter(Group != "Post") %>% 
  # wilcox.test(mean_acti~Group,., paired = T)
  aov(mean_acti~Group, .)
 
dat_combine_cluster_mouse_sta <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  ddply(., .(ID,Group, cluster ), summarise, mean_acti = mean(mean)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==1) %>% 
  ddply(., .(Group), summarize, mean_acti1 = mean(mean_acti), n = length(mean_acti),se=sd(mean_acti)/sqrt(length(mean_acti))) 

summary(t_combine_cluster_mouse)
TukeyHSD(t_combine_cluster_mouse)



dat_for_pair <- dat_combine_cluster_mouse %>% 
  filter(Group != "Pre") 

pairwise.wilcox.test(dat_for_pair$mean_acti, dat_for_pair$Group, p.adjust.methods="bonf")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_combine_cluster.pdf", width = 65/25.6, height = 65/25.6, family = "Arial" )
p_combine_cluster_mouse
dev.off()



## amplitude of each neurons in this cluster
rm(list = setdiff(ls(), "dat_cell_trace_cluster"))

c_miniscope_matlab_ft_time <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_spike_time <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Event_time/", pattern = ID, full.names = T)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale) 
    
    dat_event_time <- read.xlsx(dat_spike_time[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      as.matrix()
    
    # Create a new matrix to hold the result
    dat_event_peak <- matrix(0, nrow = nrow(dat_event_time), ncol = ncol(dat_event_time))
    
    # Identify positions where event_time_matrix is 1
    event_positions <- dat_event_time == 1
    
    # Replace "1"s in result_matrix with corresponding values from data_trace_matrix
    dat_event_peak[event_positions] <- dat_trace[event_positions]
    
    # Optionally, convert result_matrix back to a data frame
    dat_event_peak <- as.data.frame(dat_event_peak)
    
    # result_df now contains the raw values from data_trace where event_time had "1"s, and "0" otherwise
    
    
    ## number of rows to be binned
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    
    if (t1_p < 40){
      column_means = colMeans(dat_event_peak[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_event_peak))
      
      # Combine the new matrix with the original one
      dat_event_peak = rbind(new_rows, dat_event_peak)
      
    } else {
      dat_event_peak <- dat_event_peak
    }
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim <- dat_event_peak[(t1_p-39):(t1_p+40),] 
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace_event <- mapply(c_miniscope_matlab_ft_time, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace_event <- lapply(dat_cell_trace_event, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) 


dat_cell_cond_trace_event <- lapply(dat_cell_trace_event, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group) 


dat_cb_combine_event <- lapply(dat_cell_trace_event, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(dat_cell_pre_trace_event, dat_cell_cond_trace_event,.)

h_cell_peak <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(value != 0) %>% 
  filter(cluster ==1) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(x = value, color = Group)) +
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="Peak amplitude (s.d.)", y="Cummulative probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme(legend.position = c(0.1, 0.8), legend.title = element_blank())

t_cell_peak <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(value != 0) %>% 
  filter(cluster ==1) %>% 
  filter(Group != "Post") %>% 
  ks.test(value~Group,.)


h_cell_peak_cluster2 <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(value != 0) %>% 
  filter(cluster ==2) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(x = value, color = Group)) +
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="Peak amplitude (s.d.)", y="Cummulative probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme(legend.position = c(0.1, 0.8), legend.title = element_blank())


p_amp_ecdf <- plot_grid(h_cell_peak, h_cell_peak_cluster2, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_amp_ecdf.pdf", width = 110/25.6, height = 60/25.6, family = "Arial" )
p_amp_ecdf
dev.off()

## calculate the mean and sd of spikes in Pre group
dat_cell_pre_trace_event_sta <- lapply(dat_cell_trace_event, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na() %>% 
  ddply(., .(Group), summarise, mean = mean(value), sd = sd(value), se=sd(value)/sqrt(length(value))) 




p_hist_pre <- lapply(dat_cell_trace_event, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na() %>% 
  ggplot(.,aes(x=value))+
  stat_ecdf(geom = "step", col = "#8491B4FF")+
  labs(x="Peak amplitude (s.d.)", y="Cummulative probability")+
   theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")
  


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_hist_pre.pdf", width = 60/25.6, height = 60/25.6, family = "Arial" )
p_hist_pre
dev.off()

## the mean + 1sd is 3 
p_peak_amp <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(value != 0) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
  filter(value > 3) %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Peak amplitude (s.d)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(3, 15))+
  theme(legend.title = element_blank(), legend.position = "none")

t_peak_amp <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(value != 0) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
  filter(value > 3) %>% 
  filter(cluster==2) %>% 
  aov(value~Group,.)

dat_peak_amp_sta <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
  filter(value > 3) %>% 
  filter(cluster==2) %>% 
  ddply(.,.(Group), summarise, mean= mean(value), n = length(value), se=sd(value)/sqrt(length(value)))



summary(t_peak_amp)
TukeyHSD(t_peak_amp)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_peak_amp.pdf", width = 70/25.6, height = 65/25.6, family = "Arial" )
p_peak_amp
dev.off()


## calculate the firing frequency of these events with larger amplitude 

p_peak_amp_frequency <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(., .(name, Group, cluster), summarise, n = sum(value != 0)) %>% 
  mutate(Freq = n/4) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
  ggplot(., aes(Group, Freq, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group), width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Freq. of Ca2+ events (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 10))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_firing_frequency_sta <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  #filter(value < 1) %>% 
  ddply(., .(name, Group, cluster), summarise, n = sum(value != 0)) %>% 
  mutate(Freq = n/4) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
  ddply(.,.(Group, cluster), summarise, mean= mean(Freq), n = length(Freq), se=sd(Freq)/sqrt(length(Freq)))

t_peak_amp_frequency <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(., .(name, Group, cluster), summarise, n = sum(value != 0)) %>% 
  mutate(Freq = n/4) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==2) %>% 
  filter(Group != "Post") %>% 
  wilcox.test(Freq~Group,., paired = T)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_peak_freq.pdf", width = 70/25.6, height = 65/25.6, family = "Arial" )
p_peak_amp_frequency
dev.off()





## d^2, cross-day aligned neurons-----
back_crossing_frame <- tibble(ID_mouse = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame_d3= c(436, 938,1670, 808, 1562, 774),
                              Frame_d6 = c(3570, 2116, 2794, 4100, 2912,2002),
                              Frame_d7 = c(1974, 2588, 1612, 1294, 1562, 2698))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1, 5)
  t_crossing <- c(1,  5)
  
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    select(V1, V5)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    t1_p_back <- t_crossing_back[i]
    
    ## for corssing less than 2s
    if (t1_p < 40){
      column_means = colMeans(dat_trace[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_trace))
      
      # Combine the new matrix with the original one
      dat_trace_add = rbind(new_rows, dat_trace)
      
    } else {
      dat_trace_add <- dat_trace
    }
    ## extract cell activity when they cross the border
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_trace_add[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim_cross <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    colnames(dat_stim_cross) <- str_c(ID,"Cell", 1: ncol(dat_stim_cross))
    
    ## extract cell activity when cross back
    dat_stim2 <- dat_trace[(t1_p_back-40):(t1_p_back+40-1),] 
    
    dat_stim_back <- aggregate(dat_stim2,list(rep(1:(nrow(dat_stim2)%/%n+1),each=n,len=nrow(dat_stim2))),mean)[-1]
    colnames(dat_stim_back) <- str_c(ID,"Cell", 1: ncol(dat_stim_back))
    
    
    ## calculate the d' of each cells
    comp_time <- seq(-2, 1.9, by=0.1)
    dat_stim_cross_sta <- dat_stim_cross %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_cross=mean(value),sd_cross=sd(value))
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_back=mean(value),sd_back=sd(value))
    
    dat_stim_combine <- full_join(dat_stim_cross_sta, dat_stim_back_sta, by = 'name') %>% 
      mutate(sd_pool = sqrt(sd_cross^2 + sd_back^2)/sqrt(2)) %>% 
      mutate(sd_pool = ifelse(sd_pool < 1e-06, 1e-06, sd_pool)) %>% 
      mutate(d = (mean_cross - mean_back)/sd_pool) %>% 
      mutate(d2 = d^2) %>% 
      dplyr::select(name, d2) %>% 
      mutate(ID = ID)
    
    dat_stim_trace[[i]] <- dat_stim_combine
    
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre") 
  
# 
# dat_cell_d_cond <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
#   do.call(rbind,.) %>% 
#   mutate(Group = "Cond")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre,.) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(., .(ID, cluster, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


dat_cell_d_combine_sta <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(., .(ID, Group, cluster), summarise,n=length(d2),value=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ddply(., .(Group, cluster),summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_d_combine <- dat_cell_d_combine %>% 
  filter(cluster ==1) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF",  "#3C5488FF"))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 60))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine <- dat_cell_d_combine %>% 
  filter(cluster==1) %>% 
  wilcox.test(mean~Group, ., paired = T)

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre,.) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  filter(cluster==1) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF",  "#3C5488FF"))+
  labs(x="", y="(d')^2" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 500))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre" ,"Post"))) %>% 
  filter(cluster ==1) %>% 
  wilcox.test(d2~ Group,., paired = T)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_combine.pdf", width = 45/25.6, height = 60/25.6, family = "Arial")
p_d_combine
dev.off()


## compare the neuron activity and their time of crossing------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,3,4, 5)
  t_crossing <- c(1,3, 4, 5)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V3, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  group_day <- c("Pre","D5", "Cond", "Post")
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_t_crossing <- tibble(Crossing = dat_trace1$crossing[t_crossing], Group = group_day, ID = ID)
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      #apply(., 2, function(x) runmed(x, k = 41)) %>% 
      apply(., 2, scale) 
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    
    if (t1_p < 40){
      column_means = colMeans(dat_trace[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_trace))
      
      # Combine the new matrix with the original one
      dat_trace = rbind(new_rows, dat_trace)
      
    } else {
      dat_trace <- dat_trace
    }
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_trace[(t1_p-39):(t1_p+40),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(list(dat_stim_trace, dat_t_crossing))
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_trace1 <- dat_cell_trace %>% 
  lapply(., function(x) x[[1]])

comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre","D5", "Cond", "Post")
score_range <- range(dat_cell_trace1)

for(i in c(1:4)){
  cell_order <- lapply(dat_cell_trace1, function(x)  x[[3]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
  p_heat <- lapply(dat_cell_trace1, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(name = factor(name, levels = cell_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value))+ 
    geom_tile(height=2)+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="", y="")+
    geom_vline(xintercept = 0, col= 'red', linetype=2)+
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
    theme(legend.position = "none")
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Pre, p_heat_D5,p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_pc.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## for individual mouse
dat_cell_pre_trace <- lapply(dat_cell_trace1, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_d5_trace <- lapply(dat_cell_trace1, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "D5") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_cond_trace <- lapply(dat_cell_trace1, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_ca_trace <- lapply(dat_cell_trace1, function(x) x[[4]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) %>% 
  rbind(dat_cell_pre_trace, dat_cell_d5_trace, dat_cell_cond_trace,.) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre","D5", "Cond", "Post"))) %>% 
  filter(cluster ==2) %>% 
  ddply(., .(name, Group, ID), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(Group, ID), summarise, mean_acti1 = mean(mean_acti), n= length(mean_acti)) %>% 
  as_tibble()


dat_t_crossing <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Crossing = Crossing/10) %>% 
  mutate(Group = factor(Group, levels = c("Pre","D5", "Cond", "Post")))

p_Ca_cross_correlation <- left_join(dat_cell_ca_trace, dat_t_crossing) %>% 
  filter(Crossing < 200) %>% 
  ggplot(., aes(mean_acti1, Crossing))+
  geom_point(aes(colour = Group),  size=2)+
  geom_smooth(method = "lm", se=F)+
  labs(x="Mean activity of PCs in cluster 1 (s.d.)", y="Latency of 1st crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme(legend.position = 'none')


cor.test(p_Ca_cross_correlation$mean_acti1, p_Ca_cross_correlation$Crossing)  

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_Ca_cross_correlation.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_Ca_cross_correlation
dev.off()


## random crossing of PC neurons, placebo state---------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  
  length_pre <- dat_trace1$traces[[1]] %>% 
    nrow()
  t_crossing_pre <- sample(c(40:(length_pre-40)), 100)
  
  length_cond <- dat_trace1$traces[[4]] %>% 
    nrow()
  t_crossing_cond <- sample(c(40:(length_cond-40)), 100)
  
  length_post <- dat_trace1$traces[[5]] %>% 
    nrow()
  t_crossing_post <- sample(c(40:(length_post-40)), 100)
  
  t_crossing <- list(t_crossing_pre, t_crossing_cond, t_crossing_post)
  
  global_ID_mouse <- dat_global_cell_combine %>% 
    filter(mouse_ID == ID)
  
  
  
  dat_stim_trace <- vector(mode = 'list', length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(global_ID_mouse[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      #select(all_of(global_cell)) %>% 
      #apply(., 2, function(x) runmed(x, k = 41)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    
    t_crossing_day <- t_crossing[[i]]
    t_crossing_mean <- c()
    
    for (j in seq_along(t_crossing_day)){
      t1_p <- t_crossing_day[j]
      
      dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
      
      dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1] %>% 
        apply(., 2, mean) 
      
      t_crossing_mean <- cbind(t_crossing_mean,dat_stim)   
      
      
    }
    
    dat_stim_trace[[i]] <- t_crossing_mean %>% 
      as_tibble() %>% 
      apply(., 1, mean) %>% 
      as_tibble()
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  unlist() 

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  unlist()

dat_cell_post_trace <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  unlist()
p_cell_acti <- c(dat_cell_pre_trace, dat_cell_cond_trace, dat_cell_post_trace) %>% 
  as_tibble() %>% 
  add_column(Group = c(rep("Pre", length(dat_cell_pre_trace)), rep("Cond", length(dat_cell_cond_trace)), rep("Post", length(dat_cell_post_trace)))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1.5, 8))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_anti_pc.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_cell_acti
dev.off()