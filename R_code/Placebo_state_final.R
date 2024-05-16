## final figures for Ca2+ imaging data

## get the cell ID of cells become more active
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

cell_order <- c_miniscope_matlab_ft(mouse_file[[3]])[[3]] %>% 
  apply(., 2, mean) 

cell_order_ID <- cell_order[order(cell_order)] %>% 
  names() %>% 
  tail(.,10) %>% 
  readr::parse_number()



  


## plot the pre and post of cell traces for main figure---

dat_trace1 <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))[[3]] %>% 
  raveio::read_mat(.)

cross_ID <- dat_trace1$global_map %>% 
  as_tibble() %>% 
  select(V1, V4, V5) %>% 
  mutate_all(na_if, 0) %>% 
  drop_na() %>% 
  filter(V5 %in% cell_order_ID)



plot_range <- range(-1.1, 7.8)
## for pre
global_cell_pre <- pull(cross_ID[,1])

t_crossing_pre <- c(dat_trace1[2][[1]][1], 2293, 2517)

p_pre_trace <- dat_trace1[[4]] %>% 
  as_tibble() %>% 
  select(all_of(global_cell_pre)) %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  rowid_to_column(., "Time") %>% 
  pivot_longer(-Time) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  scale_y_continuous(limits = plot_range)+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank()) +
  annotate("rect", xmin = t_crossing_pre-40, xmax = t_crossing_pre + 40, ymin = -Inf, ymax = Inf,
           alpha = .2,fill = "blue")

## for cond
global_cell_cond <- pull(cross_ID[,2])

t_crossing_cond <- dat_trace1[2][[1]][4]
p_cond_trace <- dat_trace1[[4]] %>% 
  as_tibble() %>% 
  select(all_of(global_cell_cond)) %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  rowid_to_column(., "Time") %>% 
  pivot_longer(-Time) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  facet_grid(rows = vars(name))+
  scale_y_continuous(limits = plot_range)+
  theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank()) +
  annotate("rect", xmin = t_crossing_cond-40, xmax = t_crossing_cond + 40, ymin = -Inf, ymax = Inf,
           alpha = .2,fill = "blue")

## for traces after conditoning
global_cell_post <- pull(cross_ID[,3])

t_crossing_post <- c(dat_trace1[2][[1]][5], 1003, 3942)
p_post_trace <- dat_trace1[[8]] %>% 
  as_tibble() %>% 
  select(all_of(global_cell_pre)) %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  rowid_to_column(., "Time") %>% 
  pivot_longer(-Time) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  facet_grid(rows = vars(name))+
  theme_void()+
  scale_y_continuous(limits = plot_range)+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank())+
  annotate("rect", xmin = t_crossing_post-40, xmax = t_crossing_post + 40, ymin = -Inf, ymax = Inf,
           alpha = .2,fill = "blue")

p_traces <- plot_grid(p_pre_trace, p_cond_trace,p_post_trace, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_traces
dev.off()


## plot the traces, with random activity--------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## heatmap plot the activity, ordered by cell activity
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")
score_range <- range(dat_cell_trace)

for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
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
    theme(legend.position = "none")
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## calculate cell activity change mouse by mouse (cell becomes more active)
cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 3)) %>% 
    add_column(Group = rep(c("Pre", "Cond", "Post"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) 
  
  dat_trace_acti <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = m_id)
  
  # calculate percentage of cell numer increased activity
  ratio_pre_con <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Cond - Pre) %>% 
    dplyr::filter(diff > 0)
  
  ratio_pre_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Post - Pre) %>% 
    dplyr::filter(diff > 0)
  
  ratio_con_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Post - Cond) %>% 
    dplyr::filter(diff > 0)
  
  
  dat_ratio <- c(nrow(ratio_pre_con), nrow(ratio_pre_post), nrow(ratio_con_post))/ncol(dat_trace[[1]])
  
  name_con_pre <- ratio_pre_con$name
  name_post_pre <- ratio_pre_post$name
  name_post_con <- ratio_con_post$name
  return(list(dat_trace_acti, dat_ratio, name_con_pre, name_post_pre, name_post_con))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  rename(name = Group, value= cell_acti ) %>% 
  bind_rows(., p_random) %>% 
  add_column(Group = rep(c("Real", "Random"), each = nrow(p_random))) %>% 
  mutate(Group = factor(Group, levels = c("Real", "Random"))) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ddply(., .(Group, name), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(name, y = mean, group=Group, colour= Group))+
  geom_point(aes(shape = Group), size = 2 )+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line()+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 0.6))+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))



p_test <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  aov(cell_acti~ Group,.)

summary(p_test)



dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "greater")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_dat_anti
dev.off()


## d^2, cross-day aligned neurons-----
back_crossing_frame <- tibble(ID_mouse = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame_d3= c(2292, 498, 610, 535,1044,1770),
                              Frame_d7 = c(1003, 1632, 2113,1923,1572,1892))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4, 8)
  t_crossing <- c(1, 5)
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    t1_p_back <- t_crossing_back[i]
    
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


dat_cell_d_combine_sta <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),value=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ddply(., .(Group),summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 60))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 600))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1,rel_widths = c(.95, 1) )

dat_combine_cell_sta <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ddply(., .(Group),summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2)))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 70/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()

## plot the trace and show event detection-------
dat_trace1 <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))[[3]] %>% 
  raveio::read_mat(.)

cross_ID <- dat_trace1$global_map %>% 
  as_tibble() %>% 
  select(V1, V4, V5) %>% 
  mutate_all(na_if, 0) %>% 
  drop_na() %>% 
  filter(V5 %in% cell_order_ID)



plot_range <- range(-1.1, 7.8)
## for pre
global_cell_pre <- pull(cross_ID[,1])


p_pre_trace <- dat_trace1[[4]] %>% 
  as_tibble() %>% 
  select(all_of(global_cell_pre)) %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  rowid_to_column(., "Time") %>% 
  pivot_longer(-Time) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  scale_y_continuous(limits = plot_range)+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank()) 

p_pre_event <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/m3_d3_spike_time.xlsx", colNames = F, rowNames = F) %>%
  as_tibble() %>% 
  select(all_of(global_cell_pre)) %>% 
  rowid_to_column(., "Time") %>% 
  pivot_longer(-Time) %>% 
  filter(value == 1) %>% 
  ggplot(., aes(Time, value))+
  geom_tile()+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank()) 

p_event <- plot_grid(p_pre_trace, p_pre_event, ncol = 1)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_event.pdf", width = 90/25.6, height = 90/25.6, family = "Arial")
p_event
dev.off()

## calculate the spiking rate-----

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_event_day <- rep(0, 3)
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    dat_event <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    dat_event_day[i] <- dat_event1
    
  }
  return(list(dat_stim_trace, dat_event_day))
}



## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## plot the firing rate

p_dat_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_dat_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  aov(value~name,.)

summary(t_dat_firing)
TukeyHSD(t_dat_firing)

dat_cell_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post")))

pairwise.t.test(dat_cell_firing$value, dat_cell_firing$name, paired = T, alternative = "greater")


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_dat_firing
dev.off()


## for comparing crossing, back-crossing and last crossing------
back_crossing_frame <- tibble(ID = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame = c(1003, 1632, 2113,1923,1572,1892), 
                              Last_frame= c(1756, 2777, 2536, 2238, 3062, 2306))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    filter(ID== ID_mouse )
  t_crossing <- c(dat_trace1[[2]][5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID_mouse, full.names = T)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  dat_event_day <- rep(0, 3)
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1[[8]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    dat_event <- read.xlsx(dat_spike[3], colNames = F) %>% 
      as_tibble() 
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- t_crossing[i]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    dat_event_day[i] <- dat_event1
  }
  return(list(dat_stim_trace, dat_event_day))
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


p_dat_firing_d7 <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Crossing =V1, Crossing_back= V2, Last_crossing = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_dat_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Crossing =V1, Crossing_back= V2, Last_crossing = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(value~name,.)

summary(t_dat_firing)
TukeyHSD(t_dat_firing)

dat_firing_sta <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Crossing =V1, Crossing_back= V2, Last_crossing = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ddply(., .(name), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  


## compare each cells
p_cell_firing_corssing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Crossing")
  
p_cell_firing_corssing_back <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Crossing_back")

p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(mean~Group,.)


cell_firing_sta <-  lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, value1 = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ddply(., .(Group), summarise,n=length(value1),mean=mean(value1),sd=sd(value1),se=sd(value1)/sqrt(length(value1)))
  



summary(t_cell_firing)
TukeyHSD(t_cell_firing)

p_d7_combin <- plot_grid(p_cell_firing, p_dat_firing_d7, nrow = 1, rel_widths = c(1, 0.8))
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing_d7.pdf", width = 90/25.6, height = 60/25.6, family = "Arial" )
p_d7_combin
dev.off()


## time profile of ca2+ activity-----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_event_day <- data.frame(Pre= numeric(80), Cond= integer(80), Post = integer(80))
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    dat_event <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>%
      rowMeans() 
    
    
    dat_event_day[,i]<- dat_event1
    
  }
  dat_event_day$ID <- ID
  return(list(dat_stim_trace, dat_event_day))
}



## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## plot the firing rate

p_dat_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = rep(seq(-2, 2, length.out = 80), 6)) %>% 
  pivot_longer(-c(ID, Time)) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ddply(.,.(Time, name), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, color = name))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=name), alpha=0.2, linetype=0)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="Time (s)", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 0.4))+
  theme(legend.title = element_blank())

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing_time.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_dat_firing
dev.off()


## d^2 cross-day aligned neurons-----
back_crossing_frame <- tibble(ID_mouse = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame_d3= c(2292, 498, 610, 535,1044,1770),
                              Frame_d7 = c(1003, 1632, 2113,1923,1572,1892))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4, 8)
  t_crossing <- c(1, 5)
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)[c(1,3)]
  
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    t1_p_back <- t_crossing_back[i]
    
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 8))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.5, 250))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1, rel_widths = c(0.8,1))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()


## ## for d3, d5, d6, d7 and latency to first crossing -------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  
  num_compare <- c(4,6:8)
  t_crossing <- c(1,3:5)
  
  
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V3, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  
  comp_time <- seq(-2, 1.9, by=0.1)
  
  dat_trace_combine <- dat_stim_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 4)) %>% 
    add_column(Group = rep(c("D3", "D5", "D6", "D7"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = ID) %>% 
    add_column(Latency_cross = dat_trace1$crossing[t_crossing]/20)
  
  return(dat_trace_combine)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

p_ca_crossing <- dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  #filter(Latency_cross < 150) %>% 
  ggplot(., aes(cell_acti, Latency_cross))+
  geom_point(aes(colour = ID),  size=2)+
  geom_smooth(method = "lm", se=F)+
  labs(x="Mean activity of rACC-Pn neurons (s.d.)", y="Latency of 1st crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme(legend.title = element_blank())

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ca_crossing.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_ca_crossing
dev.off()

p_cor <-  dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  filter(Latency_cross < 150) 

cor.test(p_cor$cell_acti, p_cor$Latency_cross)  


## for IT neurons-----------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))

  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = 3)
  
  
  for (i in c(1:3)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[i]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    t1_p <- dat_trace1$crossing[i]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim <- dat_trace[(t1_p-20):(t1_p+20-1),] 
    
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_EXTRACT", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## heatmap plot the activity, ordered by cell activity
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")
score_range <- range(dat_cell_trace)

for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    as_tibble() %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    as_tibble() %>% 
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

p_heat <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## calculate cell activity change mouse by mouse (cell becomes more active)
cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 3)) %>% 
    add_column(Group = rep(c("Pre", "Cond", "Post"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) 
  
  dat_trace_acti <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = m_id)
  
  # calculate percentage of cell numer increased activity
  ratio_pre_con <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Cond - Pre) %>% 
    dplyr::filter(diff > 0)
  
  ratio_pre_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Post - Pre) %>% 
    dplyr::filter(diff > 0)
  
  ratio_con_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Post - Cond) %>% 
    dplyr::filter(diff > 0)
  
  
  dat_ratio <- c(nrow(ratio_pre_con), nrow(ratio_pre_post), nrow(ratio_con_post))/ncol(dat_trace[[1]])
  
  name_con_pre <- ratio_pre_con$name
  name_post_pre <- ratio_pre_post$name
  name_post_con <- ratio_con_post$name
  return(list(dat_trace_acti, dat_ratio, name_con_pre, name_post_pre, name_post_con))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, cell_acti, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.6, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

p_test <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  aov(cell_acti~ Group,.)

summary(p_test)

dat_anti_sta <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ddply(., .(Group),summarise,n=length(cell_acti),mean=mean(cell_acti),sd=sd(cell_acti),se=sd(cell_acti)/sqrt(length(cell_acti)))

dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "greater")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_dat_anti
dev.off()

## compare all activity by all cells
dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Pre")

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Cond")

p_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  ggplot(., aes(Group, mean_acti, group = Group))+
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
  scale_y_continuous(limits = c(-1.5, 5))+
  theme(legend.title = element_blank(), legend.position = "none")

p_acti_combine <- plot_grid( p_cell_acti,p_dat_anti, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acti_combine.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_acti_combine
dev.off()

t_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  aov(mean_acti~Group,.)

summary(t_cell_acti)
TukeyHSD(t_cell_acti, ordered = T)

pairwise.t.test(t_cell_acti$mean_acti, t_cell_acti$Group, paired = T, alternative = "greater")

cell_acti_sta <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  ddply(., .(Group), summarise,n=length(mean_acti),mean=mean(mean_acti),sd=sd(mean_acti),se=sd(mean_acti)/sqrt(length(mean_acti)))


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_portion.pdf", width = 60/25.6, height = 65/25.6, family = "Arial")
p_cell_portion
dev.off()

## 2. For d^2
back_crossing_frame <- tibble(ID_mouse = c( "m18", "m19", "m20"), Frame_d3= c(2277, 303, 674),
                              Frame_d7 = c(1051, 1853, 546))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))

  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  num_compare <- c(1,3)
  t_crossing <- c(1, 3)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
 
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,num_compare[i]])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned

    t1_p <- dat_trace1$crossing[t_crossing[i]]
    t1_p_back <- t_crossing_back[i]
    
    ## extract cell activity when they cross the border
    dat_stim_cross <- dat_trace[(t1_p-20):(t1_p+20-1),] 
    
    colnames(dat_stim_cross) <- str_c(ID,"Cell", 1: ncol(dat_stim_cross))
    
    ## extract cell activity when cross back
    dat_stim_back <- dat_trace[(t1_p_back-20):(t1_p_back+20-1),] 
    
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_EXTRACT", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.5, 600))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1, rel_widths = c(0.8,1))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()


# 3. For 1st crossing, crossing back and last crossing

back_crossing_frame <- tibble(ID = c( "m18", "m19", "m20"), Frame = c(1050, 1853, 546), 
                              Last_frame= c(1886, 1975, 2225)) 

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    filter(ID== ID_mouse )
  t_crossing <- c(dat_trace1$crossing[3], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1$traces[[3]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned

    t1_p <- t_crossing[i]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim <- dat_trace[(t1_p-20):(t1_p+20-1),] 
    
    
    
    colnames(dat_stim) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_EXTRACT", full.names = T))


dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Crossing", "Crossing_back", "Last_crossing")
score_range <- range(dat_cell_trace)


for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    as_tibble() %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    as_tibble() %>% 
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

p_heat <- plot_grid(p_heat_Crossing, p_heat_Crossing_back, p_heat_Last_crossing, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## calculate cell activity change mouse by mouse
cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 3)) %>% 
    add_column(Group = rep(c("Crossing", "Crossing_back", "Last_crossing"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) 
  
  dat_trace_acti <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = m_id)
  
  # calculate percentage of cell numer increased activity
  ratio_pre_con <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Crossing_back - Crossing) %>% 
    dplyr::filter(diff > 0)
  
  ratio_pre_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Last_crossing - Crossing) %>% 
    dplyr::filter(diff > 0)
  
  ratio_con_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Last_crossing - Crossing_back) %>% 
    dplyr::filter(diff > 0)
  
  dat_ratio <- c(nrow(ratio_pre_con), nrow(ratio_pre_post), nrow(ratio_con_post))/ncol(dat_trace[[1]])
  
  return(list(dat_trace_acti, dat_ratio))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti_cross <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, cell_acti, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")



p_test <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(cell_acti~ Group,.)

summary(p_test)
TukeyHSD(p_test)

dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "less")

## compare for each cells
p_cell_firing_corssing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Crossing")

p_cell_firing_corssing_back <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Crossing_back")

p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1.2, 5))+
  theme(legend.title = element_blank(), legend.position = "none")

t_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>%
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(mean~Group,.)

summary(t_cell_firing)
TukeyHSD(t_cell_firing)


p_d7_combin <- plot_grid(p_cell_firing, p_dat_anti_cross, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing_d7.pdf", width = 80/25.6, height = 60/25.6, family = "Arial" )
p_d7_combin
dev.off()



## 4. compare Pre, cond and Post with spiking events
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_spikes/", pattern = ID, full.names = T)
  
  dat_stim_trace <- rep(0, 3)
  
  
  for (i in c(1:3)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
    ## number of rows to be binned
    t1_p <- dat_trace1$crossing[i]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim <- dat_trace[(t1_p-20):(t1_p+20-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    
    
    dat_stim_trace[i] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_EXTRACT", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


p_random <- dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c( "m18", "m19", "m20"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 0.5))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_random.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_random
dev.off()


# 5. For 1st crossing, crossing back and last crossing

back_crossing_frame <- tibble(ID = c( "m18", "m19", "m20"), Frame = c(1050, 1853, 546), 
                              Last_frame= c(1886, 1975, 2225)) 

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    filter(ID== ID_mouse )
  t_crossing <- c(dat_trace1$crossing[3], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_spikes/", pattern = ID_mouse, full.names = T)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- read.xlsx(dat_spike[3], colNames = F) %>% 
      as_tibble() 
    
    ## number of rows to be binned
    
    t1_p <- t_crossing[i]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim <- dat_trace[(t1_p-20):(t1_p+20-1),] 
    
    
    
    colnames(dat_stim) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_EXTRACT", full.names = T))


dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Crossing", "Crossing_back", "Last_crossing")
score_range <- range(dat_cell_trace)


for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    as_tibble() %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    as_tibble() %>% 
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

p_heat <- plot_grid(p_heat_Crossing, p_heat_Crossing_back, p_heat_Last_crossing, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 3)) %>% 
    add_column(Group = rep(c("Crossing", "Crossing_back", "Last_crossing"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) 
  
  dat_trace_acti <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = m_id)
  
  # calculate percentage of cell numer increased activity
  ratio_pre_con <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Crossing_back - Crossing) %>% 
    dplyr::filter(diff > 0)
  
  ratio_pre_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Last_crossing - Crossing) %>% 
    dplyr::filter(diff > 0)
  
  ratio_con_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Last_crossing - Crossing_back) %>% 
    dplyr::filter(diff > 0)
  
  dat_ratio <- c(nrow(ratio_pre_con), nrow(ratio_pre_post), nrow(ratio_con_post))/ncol(dat_trace[[1]])
  
  return(list(dat_trace_acti, dat_ratio))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti_cross <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, cell_acti, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 0.8))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing_d7.pdf", width = 40/25.6, height = 60/25.6, family = "Arial" )
p_dat_anti_cross
dev.off()

p_test <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(cell_acti~ Group,.)

summary(p_test)
TukeyHSD(p_test)

dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "less")

## compare for each cells
p_cell_firing_corssing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Crossing")

p_cell_firing_corssing_back <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Crossing_back")

p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1.2, 5))+
  theme(legend.title = element_blank(), legend.position = "none")

t_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>%
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(mean~Group,.)

summary(t_cell_firing)
TukeyHSD(t_cell_firing)


p_d7_combin <- plot_grid(p_cell_firing, p_dat_anti_cross, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing_d7.pdf", width = 80/25.6, height = 60/25.6, family = "Arial" )
p_d7_combin
dev.off()


## plot the label cells

num_m18 <- c(37, 36, 28, 18)
num_m19 <- c(58, 55, 52, 35)
num_m20 <- c(140, 139, 138, 80)

dat_num <- tibble(m18 = num_m18, m19 = num_m19, m20 = num_m20) %>% 
  add_column(Group = c("Pre", "Cond", "Post", "Align")) %>%
  pivot_longer(-Group) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post", "Align"))) %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF", "violetred3"))+
  labs(x="", y="No. of cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(10, 150))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_number.pdf", width = 55/25.6, height = 65/25.6, family = "Arial")
dat_num
dev.off()

## plot the trace to show

dat_trace1 <- raveio::read_mat("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/IT_EXTRACT/m18_IT.mat")$traces[[1]] %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  select(sample(1:37, 15)) %>% 
  add_column(Time = 1: 2428) %>% 
  pivot_longer(-Time) %>% 
  ggplot(., aes(Time, value, Group = name, col= name))+
  geom_line() +
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank())


## plot the noise and signal baseline level with data from Fatih------

p_noise_cor_all <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Noise correaltion (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 1))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_all <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  aov(value~name,.)

dat_noise_cor_all_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  ddply(., .(name), summarise, mean = mean(value))
summary(t_noise_cor_all)
TukeyHSD(t_noise_cor_all, ordered = T)


## averaged by mice
p_noise_cor_mice <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Noise correlation (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.02, 0.1))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_mice <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  aov(value~name,.)

summary(t_noise_cor_mice)
TukeyHSD(t_noise_cor_mice, ordered = T)


p_noise_all <- plot_grid(p_noise_cor_all, p_noise_cor_mice, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_noise_all.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_noise_all
dev.off()
## for signal

p_signal_cor_all <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_signal_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Signal correlation (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 1))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_all <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_signal_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  aov(value~name,.)

summary(t_signal_cor_all)
TukeyHSD(t_signal_cor_all, ordered = T)

## averaged by mice
p_signal_cor_mice <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_signal_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Signal correlation (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.1, 0.3))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_mice <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/all_signal_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  aov(value~name,.)

summary(t_signal_cor_mice)
TukeyHSD(t_signal_cor_mice, ordered = T)

p_signal_all <- plot_grid(p_signal_cor_all, p_signal_cor_mice, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_signal_all.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_signal_all
dev.off()


## For back-crossing noise

p_noise_cor_all_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_noise_cor.xlsx", rowNames = F, cols = c(1:3)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha = 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="Noise correlation during crossing back (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 1))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_all_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_noise_cor.xlsx", rowNames = F, cols = c(1:3)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Post) %>% 
  pivot_longer(-X1) %>% 
  wilcox.test(value~name,.)



## averaged by mice
p_noise_cor_mice_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_noise_cor.xlsx", rowNames = F, cols = c(1, 4,5)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Post = 3) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="Noise correlation during crossing back (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.1, 0.2))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_mice <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_noise_cor.xlsx", rowNames = F, cols = c(1, 4,5)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Post = 3) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  wilcox.test(value~name,.)


p_signal_all_back <- plot_grid(p_noise_cor_all_back, p_noise_cor_mice_back, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_signal_all_back.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_signal_all_back
dev.off()
## for signal

p_signal_cor_all_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_signal_cor.xlsx", rowNames = F, cols = c(1:3)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre",  "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="Signal correlation during crossing back (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_all_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_signal_cor.xlsx", rowNames = F, cols = c(1:3)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Post) %>% 
  pivot_longer(-X1) %>% 
  wilcox.test(value~name,.)



## averaged by mice
p_signal_cor_mice_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_signal_cor.xlsx", rowNames = F, cols = c(1, 4, 5)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Post = 3) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="Signal correlation (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.05, 0.35))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_mice_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_signal_cor.xlsx", rowNames = F, cols = c(1, 4, 5)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Post = 3) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  wilcox.test(value~name,.)

p_signal_mice_back <- plot_grid(p_signal_cor_all_back, p_signal_cor_mice_back, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_signal_mice_back.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_signal_mice_back
dev.off()
## for first crossing


p_noise_cor_all_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Noise correaltion druing first crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_all_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  aov(value~name,.)

summary(t_noise_cor_all_d7)
TukeyHSD(t_noise_cor_all_d7, ordered = T)

## averaged by mice
p_noise_cor_mice_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Noise correlation during first crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.02, 0.3))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_mice_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  aov(value~name,.)

summary(t_noise_cor_mice_d7)
TukeyHSD(t_noise_cor_mice_d7, ordered = T)

p_noise_crossing_d7 <- plot_grid(p_noise_cor_all_d7, p_noise_cor_mice_d7, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_noise_crossing_d7.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_noise_crossing_d7
dev.off()

## for signal

p_signal_cor_all_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_signal_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Signal correlation during first crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_all_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_signal_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  aov(value~name,.)

summary(t_signal_cor_all_d7)
TukeyHSD(t_signal_cor_all_d7, ordered = T)

## averaged by mice
p_signal_cor_mice_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_signal_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Signal correlation during first crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.08, 0.32))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_mice_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_signal_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  aov(value~name,.)

summary(t_signal_cor_mice_d7)
TukeyHSD(t_signal_cor_mice_d7, ordered = T)

p_signal_crossing_d7 <- plot_grid(p_signal_cor_all_d7, p_signal_cor_mice_d7, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_signal_crossing_d7.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_signal_crossing_d7
dev.off()

## for random crossing

p_noise_cor_all_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Noise correaltion \n with random crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_all_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_noise_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  aov(value~name,.)

summary(t_noise_cor_all_random)
TukeyHSD(t_noise_cor_all_random, ordered = T)

## averaged by mice
p_noise_cor_mice_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Noise correlation \n with random crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.1, 0.2))+
  theme(legend.title = element_blank(), legend.position = "none")


t_noise_cor_mice_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  aov(value~name,.)

dat_noise_cor_mice_random_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  ddply(., .(name), summarise, mean = mean(value))


summary(t_noise_cor_mice_random)
TukeyHSD(t_noise_cor_mice_random, ordered = T)


p_noise_cor_random <- plot_grid(p_noise_cor_all_random, p_noise_cor_mice_random, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
tiff("p_noise_cor_random.tiff", width= 90/25.6, height= 60/25.6, units="in", res=300)
p_noise_cor_random
dev.off()
## for signal

p_signal_cor_all_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_signal_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_violin()+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Signal correlation with \n random crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_all_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_signal_cor.xlsx", rowNames = F, cols = c(1:4)) %>% 
  as_tibble() %>% 
  select(X1, Pre, Cond, Post) %>% 
  pivot_longer(-X1) %>% 
  aov(value~name,.)

summary(t_signal_cor_all_random)
TukeyHSD(t_signal_cor_all_random)

## averaged by mice
p_signal_cor_mice_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_signal_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour =name, shape =name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Signal correlation with \n random crossing (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.1, 0.3))+
  theme(legend.title = element_blank(), legend.position = "none")


t_signal_cor_mice_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_signal_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  aov(value~name,.)

summary(t_signal_cor_mice_random)
TukeyHSD(t_signal_cor_mice_random)

p_signal_cor_random <- plot_grid(p_signal_cor_all_random, p_signal_cor_mice_random, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")

tiff("p_signal_cor_random.tiff", width= 90/25.6, height= 60/25.6, units="in", res=300)
p_signal_cor_random

dev.off()


# cairo_pdf("p_signal_cor_random.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
# p_signal_cor_random
# dev.off()


## combine first crossing, back crossing and random
p_noise_cor_mice_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/back_crossing_noise_cor.xlsx", rowNames = F, cols = c(1, 4,5)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Post = 3) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  add_column(Group = "Back_crossing")

p_noise_cor_mice_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/first_crossing_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  select(-Cond) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  add_column(Group = "First_crossing")


p_noise_cor_mice_random <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Singal_random/random_noise_cor.xlsx", rowNames = F, cols = c(1, 5:7)) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  rename(ID = 1, Pre = 2, Cond = 3, Post = 4) %>% 
  select(-Cond) %>% 
  drop_na() %>% 
  pivot_longer(-ID) %>% 
  add_column(Group = "Random_crossing")

dat_noise_com <- rbind(p_noise_cor_mice_d7, p_noise_cor_mice_back, p_noise_cor_mice_random) %>% 
  as_tibble() %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  mutate(Group = factor(Group, levels = c("First_crossing", "Back_crossing", "Random_crossing"))) %>% 
  ggplot(., aes(x=Group, y=value, fill=name)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  labs(x="", y="Noise correlation (r)")+
  theme(axis.line.x = element_line(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.05, 0.25))+
  theme(legend.title = element_blank())

t_dat_noise_com <- rbind(p_noise_cor_mice_d7, p_noise_cor_mice_back, p_noise_cor_mice_random) %>% 
  as_tibble() %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  mutate(Group = factor(Group, levels = c("First_crossing", "Back_crossing", "Random_crossing"))) %>% 
  aov(value ~Group *name, .)

summary(t_dat_noise_com)  
TukeyHSD(t_dat_noise_com)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("dat_noise_com.pdf", width = 92/25.6, height = 70/25.6, family = "Arial")
dat_noise_com
dev.off()

## compare two 
dat_noise_com <- rbind(p_noise_cor_mice_d7, p_noise_cor_mice_back, p_noise_cor_mice_random) %>% 
  as_tibble() %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  mutate(Group = factor(Group, levels = c("First_crossing", "Back_crossing", "Random_crossing"))) %>% 
  filter(name == "Pre") %>% 
  ggplot(., aes(x=Group, y=value, fill=Group)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  labs(x="", y="Noise correlation (r)", title = "Pre")+
  theme(axis.line.x = element_line(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.05, 0.25))+
  theme(legend.title = element_blank(), legend.position = "none")


dat_noise_com_post <- rbind(p_noise_cor_mice_d7, p_noise_cor_mice_back, p_noise_cor_mice_random) %>% 
  as_tibble() %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  mutate(Group = factor(Group, levels = c("First_crossing", "Back_crossing", "Random_crossing"))) %>% 
  filter(name == "Post") %>% 
  ggplot(., aes(x=Group, y=value, fill=Group)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  labs(x="", y="Noise correlation (r)", title = "Post")+
  theme(axis.line.x = element_line(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.05, 0.25))+
  theme(legend.title = element_blank(), legend.position = "none")

p_combine <- plot_grid(dat_noise_com, dat_noise_com_post, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("dat_noise_com.pdf", width = 120/25.6, height = 70/25.6, family = "Arial")
p_combine
dev.off()

t_aov <- rbind(p_noise_cor_mice_d7, p_noise_cor_mice_back, p_noise_cor_mice_random) %>% 
  as_tibble() %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  mutate(Group = factor(Group, levels = c("First_crossing", "Back_crossing", "Random_crossing"))) %>% 
  filter(name == "Post") %>% 
  aov(value ~ Group, .)

summary(t_aov)
TukeyHSD(t_aov)
