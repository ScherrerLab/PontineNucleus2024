
library(readABF)
## for the ephys data
##for the intrincit properties------
dat_ap_current <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = 1) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond")))

dat_ap_current_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = 1) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  ddply(.,.(Group, name), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  
# for RMP
p_rmp <- dat_ap_current %>% 
  select(Group, RMP) %>% 
  ggplot(., aes(Group, RMP, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Resting membrane potential (mV)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')




t_rmp <- dat_ap_current %>% 
  select(Group, RMP) %>% 
  wilcox.test(RMP~Group,.)

# for Rin
p_Rin <- dat_ap_current %>% 
  select(Group, Rin) %>% 
  ggplot(., aes(Group, Rin, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Input resistance (MΩ)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

t_Rin <- dat_ap_current %>% 
  select(Group, Rin) %>% 
  wilcox.test(Rin~Group,.)

## for Ap_peak

p_Ap_peak <- dat_ap_current %>% 
  select(Group, Ap_peak) %>% 
  ggplot(., aes(Group, Ap_peak, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Peak amplitude (mV)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


t_Ap_peak <- dat_ap_current %>% 
  select(Group, Ap_peak) %>% 
  wilcox.test(Ap_peak~Group,.)

## for sag
p_Sag <- dat_ap_current %>% 
  select(Group, Sag) %>% 
  ggplot(., aes(Group, Sag, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Sag amplitude (mV)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

t_Sag <- dat_ap_current %>% 
  select(Group, Sag) %>% 
  wilcox.test(Sag~Group,.)

## for half width of the ap
p_AP_t50 <- dat_ap_current %>% 
  select(Group, AP_t50) %>% 
  ggplot(., aes(Group, AP_t50, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Peak amplitude (mV)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

t_AP_t50 <- dat_ap_current %>% 
  select(Group, AP_t50) %>% 
  wilcox.test(AP_t50~Group,.)

## for num_ap_1st
p_Num_ap_1st <- dat_ap_current %>% 
  select(Group, Num_ap_1st) %>% 
  group_by(Group, Num_ap_1st) %>% 
  summarise(count1 = n()) %>% 
  mutate(perc = count1/sum(count1)) %>% 
  mutate(Num_ap_1st = factor(Num_ap_1st, levels = c(4, 3, 2, 1))) %>% 
  ggplot(., aes(x = Group, y = perc*100, fill = Num_ap_1st)) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "", y = "Percent", fill = "Num_ap") +
  scale_fill_manual(values = c("indianred","darkgray","gray", "darkcyan"))+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ap_1st.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_Num_ap_1st
dev.off()

Num_ap_1st_sta <- dat_ap_current %>% 
  select(Group, Num_ap_1st) %>% 
  ddply(., .(Group), summarise,n=length(Num_ap_1st),mean=mean(Num_ap_1st),sd=sd(Num_ap_1st),se=sd(Num_ap_1st)/sqrt(length(Num_ap_1st)))
  
t_Num_ap_1st <- dat_ap_current %>% 
  select(Group, Num_ap_1st) %>% 
  wilcox.test(Num_ap_1st~Group,.)

p_ap_property <- plot_grid(p_rmp,  p_Num_ap_1st, nrow = 1)
p_ap_passive <- plot_grid(p_rmp, p_Rin,p_Ap_peak,p_AP_t50, p_Sag, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ap_passive.pdf", width = 200/25.6, height = 60/25.6, family = "Arial")
p_ap_passive
dev.off()


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ap_property.pdf", width = 70/25.6, height = 60/25.6, family = "Arial")
p_ap_property
dev.off()
## for the ap firing frequecny------
p_freq_current <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "current injection") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Current), summarise,n=length(Num_ap),mean=mean(Num_ap),sd=sd(Num_ap),se=sd(Num_ap)/sqrt(length(Num_ap))) %>% 
  ggplot(., aes(x =as.factor(Current), y = mean,group=Group,colour=Group))+
  geom_point(aes(shape = Group), size = 2 )+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line( )+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x= "Input Current (pA)", y= "APs Frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = c(0.2, 0.8), legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_freq_current.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_freq_current
dev.off()

## for spontaneous release-----
p_spon_freq <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, IEI) %>% 
  filter(IEI !=0) %>% 
  ggplot(., aes(IEI/10, colour= Group))+
  stat_ecdf(geom = "step") +
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x= "Inter-event interval (ms)", y= "Cummulative probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


  
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_freq.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_spon_freq
dev.off()


dat_spon_freq_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, ID,IEI) %>% 
  filter(IEI !=0) %>% 
  ddply(., .(Group, ID), summarise,IEI=mean(IEI)) %>% 
  ddply(., .(Group), summarise,n=length(IEI),mean=mean(IEI),sd=sd(IEI),se=sd(IEI)/sqrt(length(IEI)))

t_spon_freq_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, ID,IEI) %>% 
  filter(IEI !=0) %>% 
  ddply(., .(Group, ID), summarise,IEI=mean(IEI)) %>% 
  wilcox.test(IEI~Group, .)

p_spon_freq_bar <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, ID,IEI) %>% 
  filter(IEI !=0) %>% 
  ddply(., .(Group, ID), summarise,IEI=mean(IEI)) %>% 
  mutate(Freq = round(1000/IEI, digits = 2)) %>% 
  ggplot(., aes(Group, Freq, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_freq_bar.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_spon_freq_bar
dev.off()
  

## for amp
p_spon_amp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <500) %>% 
  ggplot(., aes(Amp, colour= Group))+
  stat_ecdf(geom = "step") +
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x= "Peak amplitude (pA)", y= "Cummulative probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_amp.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_spon_amp
dev.off()


dat_spon_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, ID,Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <500) %>% 
  ddply(., .(Group, ID), summarise,Amp=mean(Amp)) %>% 
  ddply(., .(Group), summarise,n=length(Amp),mean=mean(Amp),sd=sd(Amp),se=sd(Amp)/sqrt(length(Amp)))

  
t_spon_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, ID,Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <500) %>% 
  ddply(., .(Group, ID), summarise,Amp=mean(Amp)) %>% 
  wilcox.test(Amp~Group, .)


p_spon_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  select(Group, ID,Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <500) %>% 
  ddply(., .(Group, ID), summarise,Amp=mean(Amp)) %>% 
  ggplot(., aes(Group, Amp, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Peak amplitude (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_amp_sta.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_spon_amp_sta
dev.off()

## For the PPR analysis------
p_PPR <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PPR") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Interval), summarise,n=length(Ratio),mean=mean(Ratio),sd=sd(Ratio),se=sd(Ratio)/sqrt(length(Ratio))) %>% 
  mutate(Interval = factor(Interval, levels = c("20","50","100","200","500"))) %>% 
  ggplot(., aes(Interval, mean,colour=Group))+
  geom_point(aes(colour=Group))+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
  labs(x= "Time interval (ms)", y= "Paired pulse ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))
 
## plot the trace of ACC-Pn trace-----
library(readABF)

# plot the trace without drug perfusion--
trace_acc_pn <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/04012021/21401000.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  slice(501:1000) %>% 
  mutate(Time = seq(0, 100, length.out = 500)) %>% 
  gather(Sweep, value, -Time)

trace_acc_pn_sta <- trace_acc_pn %>% 
  ddply(., .(Time), summarise, mean= mean(value))

p_acc_pn_trace <- ggplot() + 
  geom_line(data = trace_acc_pn, aes(x = Time, y = value), color = "gray90") +
  geom_line(data = trace_acc_pn_sta, aes(x = Time, y = mean), color = "black") +
  theme_void()+
  scale_y_continuous(limits = c(-400, 50))+
  annotate(x=c(80,80,90,90), y=c(-300),"path")+
  annotate(x=c(90), y=c(-300, -300, -200, -200),"path")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acc_pn_trace.pdf", width = 60/25.6, height = 45/25.6, family = "Arial")
p_acc_pn_trace
dev.off() 
## plot with CNQX (10µM)  
trace_acc_pn_cnqx <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/04012021/21401004.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  slice(501:1000) %>% 
  mutate(Time = seq(0, 100, length.out = 500)) %>% 
  gather(Sweep, value, -Time)

trace_acc_pn_cnqx_sta <- trace_acc_pn_cnqx %>% 
  ddply(., .(Time), summarise, mean= mean(value))

p_acc_pn_cnqx_trace <- ggplot() + 
  geom_line(data = trace_acc_pn_cnqx, aes(x = Time, y = value), color = "gray90") +
  geom_line(data = trace_acc_pn_cnqx_sta, aes(x = Time, y = mean), color = "black") +
  theme_void()+
  scale_y_continuous(limits = c(-400, 50))+
  annotate(x=c(80,80,90,90), y=c(-300),"path")+
  annotate(x=c(90), y=c(-300, -300, -200, -200),"path")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acc_pn_trace_cnqx.pdf", width = 60/25.6, height = 45/25.6, family = "Arial")
p_acc_pn_cnqx_trace
dev.off() 
## plot trace with TTX
trace_acc_pn_TTX <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/04012021/21401002.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  slice(501:1000) %>% 
  mutate(Time = seq(0, 100, length.out = 500)) %>% 
  gather(Sweep, value, -Time)

trace_acc_pn_TTX_sta <- trace_acc_pn_TTX %>% 
  ddply(., .(Time), summarise, mean= mean(value))

p_acc_pn_TTX_trace <- ggplot() + 
  geom_line(data = trace_acc_pn_TTX, aes(x = Time, y = value), color = "gray90") +
  geom_line(data = trace_acc_pn_TTX_sta, aes(x = Time, y = mean), color = "black") +
  theme_void()+
  scale_y_continuous(limits = c(-400, 50))+
  annotate(x=c(80,80,90,90), y=c(-300),"path")+
  annotate(x=c(90), y=c(-300, -300, -200, -200),"path")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acc_pn_trace_TTX.pdf", width = 60/25.6, height = 45/25.6, family = "Arial")
p_acc_pn_TTX_trace
dev.off()
## plot trace with TTX and 4-AP
trace_acc_pn_4ap <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/04012021/21401003.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  slice(501:1000) %>% 
  mutate(Time = seq(0, 100, length.out = 500)) %>% 
  gather(Sweep, value, -Time)

trace_acc_pn_4ap_sta <- trace_acc_pn_4ap %>% 
  ddply(., .(Time), summarise, mean= mean(value))

p_acc_pn_4ap_trace <- ggplot() + 
  geom_line(data = trace_acc_pn_4ap, aes(x = Time, y = value), color = "gray90") +
  geom_line(data = trace_acc_pn_4ap_sta, aes(x = Time, y = mean), color = "black") +
  theme_void()+
  scale_y_continuous(limits = c(-400, 50))+
  annotate(x=c(80,80,90,90), y=c(-300),"path")+
  annotate(x=c(90), y=c(-300, -300, -200, -200),"path")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acc_pn_trace_4ap.pdf", width = 60/25.6, height = 45/25.6, family = "Arial")
p_acc_pn_4ap_trace
dev.off()


## plot the current injection trace-----
trace_acc_cond <- readABF("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Current\ injection/21406007.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 15000, length.out = 15000)) %>% 
  ggplot(., aes(Time, V4))+
  geom_line(color = "indianred")+
  theme_void()+
  scale_y_continuous(limits = c(-80, 40))+
  annotate(x=c(13000,13000,15000,15000), y=c(-50),"path")+ # 200ms
  annotate(x=c(15000), y=c(-50, -50, -30, -30),"path") # 20 mv

trace_acc_ctrl <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Current\ injection/21407017.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 15000, length.out = 15000)) %>% 
  ggplot(., aes(Time, V4))+
  geom_line(color = "darkcyan")+
  theme_void()+
  scale_y_continuous(limits = c(-80, 40))+
  annotate(x=c(13000,13000,15000,15000), y=c(-50),"path")+ # 200ms
  annotate(x=c(15000), y=c(-50, -50, -30, -30),"path") # 20 mv

p_current_trace <- plot_grid(trace_acc_ctrl, trace_acc_cond, ncol = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_current_trace.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_current_trace
dev.off()


## plot the trace to show spontaneous release-----
p_spon_trace_cond <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/spon_trace_examp/cond_filtered_trace_21325009.abf")[['data']][[1]] %>%
  as.vector() %>% 
  .[1800001 : 2300000] %>% 
  plot(., type="l", axes=FALSE, ylim=c(-80, 0), xlab="", ylab="", col = "indianred" )

p_spon_trace_ctrl <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/spon_trace_examp/ctrl_filtered_trace_21324001.abf")[['data']][[1]] %>%
  as.vector() %>% 
  .[1000001 : 1500000] %>% 
  plot(., type="l", axes=FALSE, ylim=c(-80, 0), xlab="", ylab="", col= "darkcyan" )

## For the PPR_electrical analysis------
p_PPR <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PPR_electrical") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group, Interval), summarise,n=length(Ratio),mean=mean(Ratio),sd=sd(Ratio),se=sd(Ratio)/sqrt(length(Ratio))) %>% 
  ggplot(., aes(Interval, mean,colour=Group))+
  geom_line()+
  geom_point(aes(shape = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x= "Time interval (ms)", y= "Paired pulse ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_continuous(trans = log10_trans(), breaks = c(20,50,100,200,500),
                     labels = c("20","50","100","200","500"))+
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())

t_ppr <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PPR_electrical") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  aov(Ratio ~ Group + Interval, .)
summary(t_ppr)  
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ppr.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_PPR
dev.off()

## plot the trace to show ppr, 50 ms
dat_ppr_trace_ctrl <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PPR/Ctrl/PPR_50ms_ctrl.abf")[['data']] %>%
  do.call(cbind,.) %>% 
  rowMeans() %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 1000, length.out = 50000)) %>% 
  ggplot(., aes(Time, value)) +
  geom_line(colour = "darkcyan") +
  theme_void()+
  scale_x_continuous(limits = c(50, 170))


dat_ppr_trace_cond <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PPR/Ctrl/PPR_50ms_cond.abf")[['data']] %>%
  do.call(cbind,.) %>% 
  rowMeans() %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 1000, length.out = 10000)) %>% 
  ggplot(., aes(Time, value)) +
  geom_line(colour = "indianred") +
  theme_void()+
  scale_x_continuous(limits = c(50, 170))

p_ppr_trace <- plot_grid(dat_ppr_trace_ctrl, dat_ppr_trace_cond, ncol = 1)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ppr_trace.pdf", width = 100/25.6, height = 70/25.6, family = "Arial")
p_ppr_trace
dev.off()
## For AMPA/NMDA ratio-----

p_ampa_nmda <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "AMPA_NMDA_Eletric") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Ratio, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="AMPA/NMDA ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_ampa_nmda_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "AMPA_NMDA_Eletric") %>% 
  ddply(., .(Group),  summarise, n=length(Ratio),mean=mean(Ratio),sd=sd(Ratio),se=sd(Ratio)/sqrt(length(Ratio)))

p_ratio <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "AMPA_NMDA_Eletric") %>% 
  wilcox.test(Ratio~Group, .)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ampa_nmda.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_ampa_nmda
dev.off()


## plot the trace to show
dat_trace_ampa <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/example_trace/", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[501:1500], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:50]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable') %>% 
  ggplot(., aes(Time, value, group =Variable))+
  geom_line()+
  theme_void()+
  scale_y_continuous(limits = c(-200, 200))+
  annotate(x=c(80,80,100,100), y=c(-150),"path")+ # 20ms
  annotate(x=c(100), y=c(-150, -150, -100, -100),"path") # 50 mv
  



## function to get the amplitude of AMPA and NMDA current from multiple traces
cc_ampa_nmda <- function (path_cell) {
  dat_trace <- as.list(list.files(path = path_cell, pattern = "*.abf", full.names = T)) %>% 
    mapply(function(x) readABF(x)[['data']], ., SIMPLIFY = F) %>% 
    mapply(function(x) do.call(cbind, x), ., SIMPLIFY = F) %>% 
    mapply(function (x) rowMeans(x), ., SIMPLIFY = F) %>% 
    mapply(function(x) x[2500:7500], ., SIMPLIFY = F) %>% 
    mapply(function(x) x - mean(x[1:350]), ., SIMPLIFY = F) 
    
  ## find the location of Ampa peak 
  dat_trace_80 <- dat_trace[[1]][500:2000]*(-1)
  ampa_position <- findpeaks(dat_trace_80, nups = 1, minpeakdistance = 200, npeaks = 3, threshold = 0) %>% 
    as_tibble() %>% 
    subset(V2 < 600) %>% 
    arrange(., desc(V1) )
  
  
  
  ## get the peak of ampa current, average 1 ms
  Ampa_peak_range <- ampa_position[1,2] %>% 
    as.numeric() %>% 
    c((.-25): (.+24)) + 500
  
  Ampa_peak_baseline <- c(1:350)
  
  Ampa_peak <- dat_trace %>% 
    mapply(function (x) mean(x[Ampa_peak_range]) - mean(x[Ampa_peak_baseline]), .) 
  
  Ampa_peak_ratio <- Ampa_peak/abs(Ampa_peak[1])
  
  
  ## get the NMDA peak, 50ms after the stimulation at 58 ms
  point_nmda <- 400+50*55
  Nmda_peak_range <- c((point_nmda-25): (point_nmda+24))
  Nmda_peak <- dat_trace %>% 
    mapply(function (x) mean(x[Nmda_peak_range]) - mean(x[Ampa_peak_baseline]), .)
  
  Nmda_peak_ratio <- Nmda_peak/Nmda_peak[8]
  ## plot the trace to show the Ampa and NMDA peaks
  
  v_hold <- seq(-80, 60, by = 20)
  ID_cell <- str_extract(path_cell, regex("Cell\\d+"))
  
  #plot(-dat_trace_80, type = "l", main = ID_cell)
  #points(ampa_position[1,2], -ampa_position[1,1], col="red")

  dat_ratio <- tibble(ID = ID_cell,Hold = v_hold, ratio_ampa = Ampa_peak_ratio, ratio_nmda = Nmda_peak_ratio)
  
  ratio_ampa_nmda <- abs(Ampa_peak[2]/Nmda_peak[8])
  
  ratio_ampa_nmda <- ifelse(ratio_ampa_nmda > 6, 6, ratio_ampa_nmda)
  
  
  return(list(dat_ratio, ratio_ampa_nmda))

}

## for the cond. group
r_ampa_nmda_cond <- as.list(list.files('~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/AMPA_NMDA/Cond/', full.names = T)) %>% 
  mapply(cc_ampa_nmda, ., SIMPLIFY = F) %>% 
  mapply(function(x) x[[2]], .)

r_ampa_nmda_ctrl <- as.list(list.files('~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/AMPA_NMDA/Ctrl/', full.names = T)) %>% 
  mapply(cc_ampa_nmda, ., SIMPLIFY = F) %>% 
  mapply(function(x) x[[2]], .)

r_num <- c( r_ampa_nmda_ctrl, r_ampa_nmda_cond)
Group <- c(rep("Ctrl", length(r_ampa_nmda_ctrl)), rep("Cond.", length(r_ampa_nmda_cond)))

p_ratio <- tibble(Group = Group, ratio = r_num) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, ratio, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="AMPA / NMDA ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')
  
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ratio.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_ratio
dev.off()

## plot the ampla and nmda point
dat_ampa_nmda_cond <- as.list(list.files('~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/AMPA_NMDA/Cond/', full.names = T)) %>% 
  mapply(cc_ampa_nmda, ., SIMPLIFY = F) %>% 
  mapply(function(x) x[[1]], ., SIMPLIFY = F) %>% 
  do.call(rbind,. ) %>% 
  as_tibble() %>% 
  mutate(Group = "Cond.")

dat_ampa_nmda_ctrl <- as.list(list.files('~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/AMPA_NMDA/Ctrl/', full.names = T)) %>% 
  mapply(cc_ampa_nmda, ., SIMPLIFY = F) %>% 
  mapply(function(x) x[[1]], ., SIMPLIFY = F) %>% 
  do.call(rbind,. ) %>% 
  as_tibble() %>% 
  mutate(Group = "Ctrl")

p_ampa <- rbind(dat_ampa_nmda_cond, dat_ampa_nmda_ctrl) %>% 
  select(!ratio_nmda) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group, Hold),  summarise, n=length(ratio_ampa),mean=mean(ratio_ampa),sd=sd(ratio_ampa),se=sd(ratio_ampa)/sqrt(length(ratio_ampa))) %>% 
  ggplot(., aes(Hold, mean, colour= Group))+
  geom_point(size=2)+
  geom_smooth(method = "lm", se = F, linetype = "dashed", size=0.5)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ampa.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_ampa
dev.off()

p_nmda <- rbind(dat_ampa_nmda_cond, dat_ampa_nmda_ctrl) %>% 
  select(!ratio_ampa) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group, Hold),  summarise, n=length(ratio_nmda),mean=mean(ratio_nmda),sd=sd(ratio_nmda),se=sd(ratio_nmda)/sqrt(length(ratio_nmda))) %>% 
  ggplot(., aes(Hold, mean, colour= Group))+
  geom_point(size=2)+
  geom_line(linetype = "dashed", size=0.5)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_nmda.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_nmda
dev.off()

## kinetic properties of the evoked EPSC------

p_amp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "Evoke current") %>% 
  select(Group, ID, Amp) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Amp, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Amplitude (pA) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_latency <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "Evoke current") %>% 
  select(Group, ID, latency) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, latency, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_risetime <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "Evoke current") %>% 
  select(Group, ID, risetime) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, risetime, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Risetime (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_t50 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "Evoke current") %>% 
  select(Group, ID, t50) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, t50, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Half duration (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_epsc_prop <- plot_grid(p_amp, p_latency, p_risetime, p_t50, ncol = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_epsc_prop.pdf", width = 90/25.6, height = 120/25.6, family = "Arial")
p_epsc_prop
dev.off()
## plot trace to show EPSC kinetic

trace_epsc_ctrl <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/EPSC/EPSC_ctrl_trace.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  unlist() %>% 
  .[2500:10000]
trace_epsc_ctrl <- trace_epsc_ctrl - mean(trace_epsc_ctrl[1:50])
  

trace_epsc_cond <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/EPSC/EPSC_cond_trace.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  unlist() %>% 
  .[2500:10000] 
trace_epsc_cond <- trace_epsc_cond - mean(trace_epsc_cond[1:50])

trace_epsc <- c(trace_epsc_ctrl, trace_epsc_cond)
trace_time <- rep(seq(0, 150, length.out = length(trace_epsc_cond)), 2)

p_trace <- tibble(Time = trace_time,Trace = trace_epsc ) %>% 
  mutate(Group = c(rep("Ctrl", length(trace_epsc_ctrl)), rep("Cond.", length(trace_epsc_cond)))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Time, Trace, color = Group))+
  geom_line()+
  theme_void()+
  theme(legend.position = "none")+
  annotate(x=c(130,130,150,150), y=c(-200),"path")+
  annotate(x=c(150), y=c(-200, -200, -150, -150),"path")


## kinetic properties of the evoked EPSC without SR95531------

p_amp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "EPSC_property_nosr") %>% 
  select(Group, ID, Amp) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Amp, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Amplitude (pA) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_latency <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "EPSC_property_nosr") %>% 
  select(Group, ID, Latency) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Latency, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_risetime <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "EPSC_property_nosr") %>% 
  select(Group, ID, Risetime) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Risetime, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Risetime (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_t50 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "EPSC_property_nosr") %>% 
  select(Group, ID, t50) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, t50, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Half duration (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_decay <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "EPSC_property_nosr") %>% 
  select(Group, ID, Decay) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Decay, group = Group)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Half duration (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")
p_epsc_prop <- plot_grid(p_amp, p_latency, p_risetime, p_t50, ncol = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_epsc_prop.pdf", width = 100/25.6, height = 140/25.6, family = "Arial")
p_epsc_prop
dev.off()
## Compare the E-I ratio----

p_EI_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Ratio, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="EPSC/IPSC ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_EI_latency <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Diff_latency1, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="EPSC-IPSC delay (ms) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

dat_EI_latency_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group), summarise,n=length(Diff_latency1),mean=mean(Diff_latency1),sd=sd(Diff_latency1),se=sd(Diff_latency1)/sqrt(length(Diff_latency1)))

dat_EI_ratio_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group), summarise,n=length(Ratio),mean=mean(Ratio),sd=sd(Ratio),se=sd(Ratio)/sqrt(length(Ratio)))

t_EI_latency <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  wilcox.test(Diff_latency1~Group, .)

t_EI_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  wilcox.test(Ratio~Group, .)

p_EI <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  mutate(c_e = E_peak/70, c_i = I_peak/80) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(c_e, c_i, group = Group))+
  geom_point(aes(colour= Group,shape = Group), size = 2)+
  geom_smooth(method = "lm",colour = "black", se= F, size = 1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="EPSC peak conductance (nS)", y="IPSC peak conductance (nS) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_latency.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_EI_latency
dev.off()

p_ratio <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  wilcox.test(Ratio~Group, .)


p_EI <- plot_grid(p_EI_ratio, p_EI_latency, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_EI
dev.off()

## correlation analysis
p_cor_test_ei_ctrl<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  filter(Group =="Ctrl") %>% 
  cor.test(~ E_peak + I_peak, .)

p_cor_test_ei_cond<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "E_I") %>% 
  filter(Group =="Cond.") %>% 
  cor.test(~ E_peak + I_peak, .)
## plot the E_I trace
EI_ctrl <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/E_I_trace/ctrl/", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[501:1500], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:50]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble()

EI_ctrl_trace <- EI_ctrl %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')

EI_ctrl_mean <- cbind(rowMeans(EI_ctrl[,1:10]), rowMeans(EI_ctrl[,11:20])) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')
  

p_EI_ctrl <- ggplot()+
  geom_line(data = EI_ctrl_trace, aes(x= Time, y = value, group = Variable), colour = "gray90")+
  geom_line(data = EI_ctrl_mean, aes(x= Time, y = value, group = Variable), colour = "darkcyan", size = 1)+
  theme_void()+
  scale_y_continuous(limits = c(-300, 550))+
  annotate(x=c(80,80,100,100), y=c(-200),"path")+ # 20ms
  annotate(x=c(100), y=c(-200, -200, -100, -100),"path") # 50 mv

## for the cond. trace
EI_cond <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/E_I_trace/cond/", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']][1:10], .) %>% 
  mapply(function(x) x[501:1500], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:50]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble()

EI_cond_trace <- EI_cond %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')

EI_cond_mean <- cbind(rowMeans(EI_cond[,1:10]), rowMeans(EI_cond[,11:20])) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')


p_EI_cond <- ggplot()+
  geom_line(data = EI_cond_trace, aes(x= Time, y = value, group = Variable), colour = "gray90")+
  geom_line(data = EI_cond_mean, aes(x= Time, y = value, group = Variable), colour = "indianred", size = 1)+
  theme_void()+
  scale_y_continuous(limits = c(-300, 550))+
  annotate(x=c(80,80,100,100), y=c(-200),"path")+ # 20ms
  annotate(x=c(100), y=c(-200, -200, -100, -100),"path") # 100 mv

p_EI <- plot_grid(p_EI_ctrl, p_EI_cond, nrow = 2)
  

## plot trace to show the feedforward inhibition follwing excitation
EI_ctrl <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/E_I_trace/E_I/Ctrl", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[501:1500], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:50]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble()

EI_ctrl_trace <- EI_ctrl %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')

EI_ctrl_mean <- rowMeans(EI_ctrl) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')


p_EI_ctrl <- ggplot()+
  geom_line(data = EI_ctrl_trace, aes(x= Time, y = value, group = Variable), colour = "gray90")+
  geom_line(data = EI_ctrl_mean, aes(x= Time, y = value, group = Variable), colour = "darkcyan", size = 1)+
  theme_void()+
  scale_y_continuous(limits = c(-150, 200))+
  annotate(x=c(80,80,100,100), y=c(-100),"path")+ # 20ms
  annotate(x=c(100), y=c(-100, -100, -50, -50),"path") # 50 mv

## for the cond trace
EI_cond <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/E_I_trace/E_I/Cond", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[501:1500], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:50]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble()

EI_cond_trace <- EI_cond %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')

EI_cond_mean <- rowMeans(EI_cond) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 1000)) %>% 
  pivot_longer(!Time, names_to = 'Variable')


p_EI_cond <- ggplot()+
  geom_line(data = EI_cond_trace, aes(x= Time, y = value, group = Variable), colour = "gray90")+
  geom_line(data = EI_cond_mean, aes(x= Time, y = value, group = Variable), colour = "indianred", size = 1)+
  theme_void()+
  scale_y_continuous(limits = c(-150, 200))+
  annotate(x=c(80,80,100,100), y=c(-100),"path")+ # 20ms
  annotate(x=c(100), y=c(-100, -100, -50, -50),"path") # 50 mv


## for the LTP analysis-------
dat_ltp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "LTP") %>% 
  group_by(ID) %>% 
  mutate_at(4, function(x) x/mean(x[1:60])*100) %>% 
  group_by(group_sweep = cut(Sweep, breaks = seq(0, 460, 10))) %>%
  ddply(., .(Group,ID, group_sweep), summarise,mean_amp=mean(Amplitude)) %>% 
  ddply(., .(Group, group_sweep),summarise,mean=mean(mean_amp),sd=sd(mean_amp),se=sd(mean_amp)/sqrt(length(mean_amp))) %>% 
  mutate(Time = rep(seq(0.5, 46, 1),2)) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Time, mean, group = Group, colour= Group))+
  geom_point(aes(colour = Group, shape = Group), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="Time (min)", y="Amplitude (norm)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(50, 280), expand = c(0, 0))+
  geom_hline(yintercept=100, linetype="dashed", color = "gray")+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.9))
  
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ltp.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
dat_ltp
dev.off()
  
dat_ltp1 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "LTP") %>% 
  group_by(ID) %>% 
  mutate_at(4, function(x) x/mean(x[1:60])*100) %>% 
  group_by(group_sweep = cut(Sweep, breaks = seq(0, 460, 10))) %>%
  ddply(., .(Group,ID, group_sweep), summarise,mean_amp=mean(Amplitude)) %>% 
  mutate(Time = rep(seq(0.5, 46, 1),10)) %>% 
  ggplot(., aes(Time, mean_amp, colour= Group))+
  geom_point()+
  facet_grid(cols = vars(ID))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="Time (min)", y="Amplitude (norm)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_hline(yintercept=100, linetype="dashed", color = "gray")+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.9))

## plot the trace to show LTP of Ctrl and cond
p_ltp_ctrl <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/LTP/Ctrl/", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[1001:3000], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:100]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 2000)) %>% 
  pivot_longer(!Time, names_to = 'Variable') %>%
  mutate(Variable = factor(Variable, levels = c("V1", "V2"))) %>% 
  ggplot(., aes(Time, value, group =Variable, colour = Variable))+
  geom_line()+
  scale_color_manual(values=c("darkcyan", alpha("darkcyan", 0.5)))+
  theme_void()+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(-200, 100))+
  annotate(x=c(80,80,100,100), y=c(-100),"path")+ # 20ms
  annotate(x=c(100), y=c(-100, -100, -50, -50),"path") # 50 mv

p_ltp_cond <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/LTP/Cond/", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[1001:3000], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:100]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 2000)) %>% 
  pivot_longer(!Time, names_to = 'Variable') %>%
  mutate(Variable = factor(Variable, levels = c("V1", "V2"))) %>% 
  ggplot(., aes(Time, value, group =Variable, colour = Variable))+
  geom_line()+
  scale_color_manual(values=c("indianred", alpha("indianred", 0.5)))+
  theme_void()+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(-200, 100))+
  annotate(x=c(80,80,100,100), y=c(-100),"path")+ # 20ms
  annotate(x=c(100), y=c(-100, -100, -50, -50),"path") # 50 mv


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ltp_ctrl.pdf", width = 40/25.6, height = 40/25.6, family = "Arial")
p_ltp_ctrl
dev.off()

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ltp_cond.pdf", width = 40/25.6, height = 40/25.6, family = "Arial")
p_ltp_cond
dev.off()


p_ltp_example <- read.xlsx("~cchen2/Downloads/LTP_example.xlsx", sheet = 1) %>% 
  mutate(Time = seq(0, 46, length.out = 460)) %>% 
  ggplot(., aes(Time, -Peak )) +
  geom_point(size = 1)+
  labs(x = "")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 500))
  
p_ltp_example_rs <- read.xlsx("~cchen2/Downloads/LTP_example.xlsx", sheet = 1) %>% 
  mutate(Time = seq(0, 46, length.out = 460)) %>% 
  ggplot(., aes(Time, Rs )) +
  geom_line()+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 8))

p_ltp_combine <- plot_grid(p_ltp_example, p_ltp_example_rs,align = "v", ncol = 1,rel_heights = c(0.7, 0.3) )

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ltp_combine.pdf", width = 75/25.6, height = 70/25.6, family = "Arial")
p_ltp_combine
dev.off()

## for the PV interneruons----
p_amp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Amp) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Amp, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Amplitude (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Amp) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group),summarise,n=length(Amp),mean=mean(Amp),sd=sd(Amp),se=sd(Amp)/sqrt(length(Amp)))

t_test_amp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Amp) %>% 
  wilcox.test(Amp~Group, .)

t_test_risetime <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Risetime) %>% 
  wilcox.test(Risetime~Group, .)

t_test_t50 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, t50) %>% 
  wilcox.test(t50~Group, .)

p_risetime <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Risetime) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Risetime, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Risetime (ms)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')
  
p_t50 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, t50) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, t50, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Half-duration (ms)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_latency <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Latency) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Latency, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Latency (ms)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_latency_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Latency) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group),summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))


t_test_latency <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Latency) %>% 
  wilcox.test(Latency~Group, .)


p_pv_ipsc <- plot_grid(p_amp, p_risetime, p_t50, p_latency, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pv_ipsc.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_pv_ipsc
dev.off()

p_amp2_amp1 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Amp2_Amp1) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Amp2_Amp1, colour = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Amp2/Amp1")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_amp10_amp1 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV") %>% 
  select(Group, Amp10_Amp1) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, Amp10_Amp1, colour = Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Amplitude (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## calculate the SD of latency and peak
dat_pv_sd <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV_variance") %>% 
  pivot_longer(., -c(ID, Group), names_to = "Variable", values_to = "value" ) %>% 
  ddply(., .(ID, Group, Variable), summarise, sd=sd(value))
  

p_sd_latency <- dat_pv_sd %>% 
  filter(Variable == "Latency") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, sd, colour = Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Amplitude (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_sd_peak <- dat_pv_sd %>% 
  filter(Variable == "Peak") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ggplot(., aes(Group, sd, colour = Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Amplitude (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')
  

## plot the PPR  
  
p_PPR <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV_PPR") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  ddply(., .(Group, Inter_stim), summarise,n=length(Ratio),mean=mean(Ratio),sd=sd(Ratio),se=sd(Ratio)/sqrt(length(Ratio))) %>% 
  ggplot(., aes(Inter_stim, mean,colour=Group))+
  geom_line()+
  geom_point(aes(shape = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x= "Time interval (ms)", y= "Paired pulse ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_continuous(trans = log10_trans(), breaks = c(20,50,100,200,500),
                     labels = c("20","50","100","200","500"))+
  theme(legend.position = c(0.2, 0.8), legend.title = element_blank())

t_ppr <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV_PPR") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  aov(Ratio ~ Group + Inter_stim, .)
summary(t_ppr)
p_PPR_100 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV_PPR") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  filter(Inter_stim == 100) %>% 
  ggplot(., aes(Group, Ratio, colour = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()

p_PPR_100 <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PV_PPR") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  filter(Inter_stim == 100) %>%
  wilcox.test(Ratio ~ Group, .)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ppr.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_PPR
dev.off()

## plot the trace to show GABA current---
p_pv_sr <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PV/", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[1001:3000], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:100]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 100, length.out = 2000)) %>% 
  pivot_longer(!Time, names_to = 'Variable') %>%
  mutate(Variable = factor(Variable, levels = c("V1", "V2"))) %>% 
  ggplot(., aes(Time, value, group =Variable, colour = Variable))+
  geom_line()+
  scale_color_manual(values=c("gray", "black"))+
  theme_void()+
  theme(legend.position = "none")+
  annotate(x=c(80,80,100,100), y=c(-600),"path")+ # 20ms
  annotate(x=c(100), y=c(-600, -600, -400, -400),"path") # 200pa 


## plot the trace to show PV interneuron firing
trace_pv_ap<- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/pv_trace/pv_ap_firing.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 15000, length.out = 15000)) %>% 
  ggplot(., aes(Time, V1))+
  geom_line(color = "black")+
  theme_void()+
  scale_y_continuous(limits = c(-80, 40))+
  annotate(x=c(13000,13000,15000,15000), y=c(-50),"path")+ # 200ms
  annotate(x=c(15000), y=c(-50, -50, -30, -30),"path") # 20 mv

## plot trace to show light evoked inhibition
PV_trace_ctrl <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PV_traces/Ctrl", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[2401:5000], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:200]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble()

PV_ctrl_trace <- PV_trace_ctrl %>% 
  mutate(Time = seq(from = 0,  to = 130, length.out =2600 )) %>% 
  pivot_longer(!Time, names_to = 'Variable')

PV_ctrl_mean <- rowMeans(PV_trace_ctrl) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 130, length.out = 2600)) %>% 
  pivot_longer(!Time, names_to = 'Variable')


p_PV_ctrl <- ggplot()+
  geom_line(data = PV_ctrl_trace, aes(x= Time, y = value, group = Variable), colour = "gray90")+
  geom_line(data = PV_ctrl_mean, aes(x= Time, y = value, group = Variable), colour = "darkcyan", size = 1)+
  theme_void()+
  scale_y_continuous(limits = c(-1500, 100))+
  annotate(x=c(100,100,120,120), y=c(-1200),"path")+ # 20ms
  annotate(x=c(120), y=c(-1200, -1200, -1000, -1000),"path") # 50 mv

## for the cond trace
PV_trace_cond <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PV_traces/Cond", pattern = "*.abf", full.names = T )) %>% 
  mapply(function(x) readABF(x)[['data']], .) %>% 
  mapply(function(x) x[2401:5000], ., SIMPLIFY = F) %>% 
  mapply(function(x) x - mean(x[1:200]), ., SIMPLIFY = F) %>% 
  do.call(cbind, .) %>% 
  as_tibble()

PV_cond_trace <- PV_trace_cond %>% 
  mutate(Time = seq(from = 0,  to = 130, length.out =2600 )) %>% 
  pivot_longer(!Time, names_to = 'Variable')

PV_cond_mean <- rowMeans(PV_trace_cond) %>% 
  as_tibble() %>% 
  mutate(Time = seq(from = 0,  to = 130, length.out = 2600)) %>% 
  pivot_longer(!Time, names_to = 'Variable')


p_PV_cond <- ggplot()+
  geom_line(data = PV_cond_trace, aes(x= Time, y = value, group = Variable), colour = "gray90")+
  geom_line(data = PV_cond_mean, aes(x= Time, y = value, group = Variable), colour = "indianred", size = 1)+
  theme_void()+
  scale_y_continuous(limits = c(-1500, 100))+
  annotate(x=c(100,100,120,120), y=c(-1200),"path")+ # 20ms
  annotate(x=c(120), y=c(-1200, -1200, -1000, -1000),"path") # 50 mv

p_pv_trace <- plot_grid(p_PV_ctrl, p_PV_cond, ncol = 1)

## plot the trace to show latency diff
p_latency_dif <- rbind(PV_ctrl_mean[251:450, ], PV_cond_mean[251:450,]) %>% 
  mutate(Group = rep(c("Ctrl", "Cond"), each = 200)) %>% 
  ggplot(., aes(Time, value, colour = Group))+
  geom_line()+
  theme_void()+
  annotate(x=c(15,15,17,17), y=c(-1000),"path")+ # 2ms
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_latency_diff.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_latency_dif
dev.off() 
  

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pv_trace.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_pv_trace
dev.off()


## behavior for PAC ephys-----
p_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PAC_ephys_behavior") %>% 
  filter(variable=="Ratio" ) %>% 
  ggplot(., aes(Day, mean, group=Group, colour = Group))+
  geom_line()+
  geom_point(aes(shape = Group, colour = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
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
  scale_y_continuous(limits = c(20, 100))+
  theme(legend.position = "none")

t_ratio_ephy <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "PAC_ephys_behavior") %>% 
  filter(variable=="Ratio" )
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pac_ephys.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_ratio
dev.off()


## for td-negative cell properties------
p_freq_current <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "TD_ctrl_current_injection") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Current), summarise,n=length(Num_ap),mean=mean(Num_ap),sd=sd(Num_ap),se=sd(Num_ap)/sqrt(length(Num_ap))) %>% 
  ggplot(., aes(x =as.factor(Current), y = mean,group=Group,colour=Group))+
  geom_point(aes(shape = Group), size = 2 )+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line( )+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x= "Input Current (pA)", y= "APs Frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = c(0.2, 0.8), legend.title = element_blank())


t_freq_current <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "TD_ctrl_current_injection") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  mutate(Current = as.factor(Current)) %>% 
  aov(Num_ap ~ Current + Group,. )

summary(t_freq_current)




setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_freq_current.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_freq_current
dev.off()

## plot the trace
# for cond.


trace_IT_current_cond <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/example_trace/IT/21o27014.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  select(1, 2, 4, 6,8) %>% 
  mutate(Time = seq(0, 1500, length.out = 15000)) %>% 
  slice(500:13000) %>% 
  gather(Sweep, value, -Time) %>% 
  ggplot(., aes(Time, value, group = Sweep))+
  geom_line(color = "indianred")+
  scale_colour_manual(values=cc)+
  theme_void()+
  scale_y_continuous(limits = c(-100, 80))+
  theme(legend.title = element_blank())

## for control
trace_IT_current_ctrl <- readABF("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/example_trace/IT/21o20000.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  select(1, 2, 4, 6,8) %>% 
  mutate(Time = seq(0, 1500, length.out = 15000)) %>% 
  slice(500:13000) %>% 
  gather(Sweep, value, -Time) %>% 
  ggplot(., aes(Time, value, group = Sweep))+
  geom_line(color = "darkcyan")+
  scale_colour_manual(values=cc)+
  theme_void()+
  scale_y_continuous(limits = c(-100, 80))+
  theme(legend.title = element_blank())


## for spontenous release
p_spon_freq <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, IEI) %>% 
  filter(IEI !=0) %>% 
  ggplot(., aes(IEI/10, colour= Group))+
  stat_ecdf(geom = "step") +
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x= "Inter-event interval (ms)", y= "Cumulative probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(),legend.position = c(0.8, 0.4))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_freq.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_spon_freq
dev.off()


dat_spon_freq_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, ID,IEI) %>% 
  filter(IEI !=0) %>% 
  ddply(., .(Group, ID), summarise,IEI=mean(IEI)) %>% 
  ddply(., .(Group), summarise,n=length(IEI),mean=mean(IEI),sd=sd(IEI),se=sd(IEI)/sqrt(length(IEI)))

t_spon_freq_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, ID,IEI) %>% 
  filter(IEI !=0) %>% 
  ddply(., .(Group, ID), summarise,IEI=mean(IEI)) %>% 
  wilcox.test(IEI~Group, .)

p_spon_freq_bar <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, ID,IEI) %>% 
  filter(IEI !=0) %>% 
  ddply(., .(Group, ID), summarise,IEI=mean(IEI)) %>% 
  mutate(Freq = round(1000/IEI, digits = 2)) %>% 
  ggplot(., aes(Group, Freq, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_freq_bar.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_spon_freq_bar
dev.off()

## for amp
p_spon_amp <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <100) %>% 
  ggplot(., aes(Amp, colour= Group))+
  stat_ecdf(geom = "step") +
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x= "Peak amplitude (pA)", y= "Cummulative probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_amp.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_spon_amp
dev.off()


dat_spon_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, ID,Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <500) %>% 
  ddply(., .(Group, ID), summarise,Amp=mean(Amp)) %>% 
  ddply(., .(Group), summarise,n=length(Amp),mean=mean(Amp),sd=sd(Amp),se=sd(Amp)/sqrt(length(Amp)))


t_spon_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, ID,Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <500) %>% 
  ddply(., .(Group, ID), summarise,Amp=mean(Amp)) %>% 
  wilcox.test(Amp~Group, .)


p_spon_amp_sta <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/ephys/PN_ephys.xlsx", sheet = "td_ctrl_mini") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond."))) %>% 
  select(Group, ID,Amp) %>% 
  mutate(Amp = -Amp) %>% 
  filter(Amp > 0 & Amp <100) %>% 
  ddply(., .(Group, ID), summarise,Amp=mean(Amp)) %>% 
  ggplot(., aes(Group, Amp, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Peak amplitude (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_spon_amp_sta.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_spon_amp_sta
dev.off()

## plot trace to show opto inhibtiion of Pn neurons----

trace_pn_opto_before <- readABF("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/23518011.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 15000, length.out = 15000)) %>% 
  ggplot(., aes(Time, V3))+
  geom_line(color = "black")+
  theme_void()+
  scale_y_continuous(limits = c(-80, 40))+
  annotate(x=c(13000,13000,15000,15000), y=c(-50),"path")+ # 200ms
  annotate(x=c(15000), y=c(-50, -50, -30, -30),"path") # 20 mv


trace_pn_opto_light <- readABF("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/23518012.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 15000, length.out = 15000)) %>% 
  ggplot(., aes(Time, V3))+
  geom_line(color = "black")+
  theme_void()+
  scale_y_continuous(limits = c(-100, 40))+
  annotate(x=c(13000,13000,15000,15000), y=c(-50),"path")+ # 200ms
  annotate(x=c(15000), y=c(-50, -50, -30, -30),"path") # 20 mv

trace_pn_opto_after <- readABF("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/23518013.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, 15000, length.out = 15000)) %>% 
  ggplot(., aes(Time, V3))+
  geom_line(color = "black")+
  theme_void()+
  scale_y_continuous(limits = c(-80, 40))+
  annotate(x=c(13000,13000,15000,15000), y=c(-50),"path")+ # 200ms
  annotate(x=c(15000), y=c(-50, -50, -30, -30),"path") # 20 mv

p_opto_trace <- plot_grid(trace_pn_opto_before, trace_pn_opto_light, trace_pn_opto_after, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_opto_trace", width = 90/25.6, height = 60/25.6, family = "Arial")
p_opto_trace
dev.off()


## Plot light-induced voltage change-----

trace_opto_vol <- readABF("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/23518010.abf")[['data']] %>% 
  map(., function(x) x[,1]) %>% 
  do.call(cbind,. ) %>% 
  as_tibble() %>%
  slice(1000:4000) %>% 
  mutate(Time = seq(0, 300, length.out = 3001)) %>% 
  pivot_longer(-Time)

trace_opto_vol_mean <- trace_opto_vol %>% 
  ddply(., .(Time), summarise, mean= mean(value))

p_trace_vol <-
  ggplot()+
  geom_line(data = trace_opto_vol, aes(x = Time, y = value, group = name),color = "gray90")+
  geom_line(data= trace_opto_vol_mean, aes(x = Time, y = mean), color = "red")+
  theme_void()+
  scale_y_continuous(limits = c(-100, 300))+
  annotate(x=c(300,300,250,250), y=c(-50),"path")+ # 200ms
  annotate(x=c(300), y=c(-50, -50, 0, 0),"path") # 20 pA

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_opto_trace", width = 90/25.6, height = 60/25.6, family = "Arial")
p_opto_trace
dev.off()

## plot the current induced by light

p_light <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/Pn_ephys.xlsx", sheet = 1) %>% 
  ggplot(., aes(1, Value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, shape = 1)+
  labs(x="", y="Light-induced current (pA)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_light.pdf", width = 30/25.6, height = 60/25.6, family = "Arial")
p_light
dev.off()

## plot the Ap freq change 
p_light_ap <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/Pn_ephys.xlsx", sheet = 2) %>% 
  mutate(Treat= factor(Treat, levels= c("Before", "Light", "After"))) %>% 
  ggplot(., aes(Treat, Firing, group = Treat)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = Treat),width = 0.2,  size=2)+
  geom_line(aes(group = ID), colour="gray90")+
  labs(x="", y="AP frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_light_ap.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_light_ap
dev.off()


## statistical test

t_light_ap <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/ephys/Pn_opto_inhibition/Pn_ephys.xlsx", sheet = 2) %>% 
  mutate(Treat= factor(Treat, levels= c("Before", "Light", "After"))) %>% 
  aov(Firing~Treat,.)
summary(t_light_ap)
TukeyHSD(t_light_ap)

