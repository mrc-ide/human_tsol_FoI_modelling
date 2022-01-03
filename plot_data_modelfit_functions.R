#=============================================================================================================================#
#                           plot data & model fit functions                                                                   #
#=============================================================================================================================#


plot_ageprev_func <- function(data) {
  
  p <- ggplot() +    
  geom_point(data=data, aes(x=age, y=prev))+
  #geom_errorbar(data=predicted_simple,aes(x=age, y=prev, ymin=lower, ymax=upper, width=wd))+
  geom_errorbar(data=data,aes(x=age, y=prev, ymin=lower, ymax=upper))+
  #facet_wrap(~ref, scales = "free")+
  #facet_wrap(~ref)+
  ylim(0,1)+
  labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
  theme_bw() +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=11, face= "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title.x = element_text(size = 18, face= "bold"),
        axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
        legend.position = c(0.83, 0.15),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=16))
  
  return(p)
}

# for multiple datasets
plot_ageprev_func2 <- function(data) {
  
  p <- ggplot() +    
    geom_point(data=data, aes(x=age, y=prev))+
    #geom_errorbar(data=predicted_simple,aes(x=age, y=prev, ymin=lower, ymax=upper, width=wd))+
    geom_errorbar(data=data,aes(x=age, y=prev, ymin=lower, ymax=upper))+
    #facet_wrap(~ref, scales = "free")+
    facet_wrap(~ref)+
    ylim(0,1)+
    labs(x="Age (months) of human host", y="Seroprevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  return(p)
}