library(ggplot2)
library(survival)

cr.func <- function(data, model.fit, model.background.risk, data.title, preds.name, age.covariate=T){
  plot.data <- data.frame(expand.grid(age = c(0, 1),  hpv16=c(0, 1), cyt = c(1, 0)))
  plot.data$names <- levels(interaction(c("Age < 40, ", "Age \u2265 40, "), c("HPV16-, ", "HPV16+, "), 
                              c("abnormal cytology", "normal cytology"),sep=''))
  plot.data$names <- factor(plot.data$names, levels = plot.data$names[c(5, 1, 7, 3, 6, 2, 8, 4) ])
  
  if (!age.covariate) {
    plot.data <- plot.data[plot.data$age==0,]
    pi.covars <- c("hpv16", "cyt")
  }else{
    pi.covars <- c("age", "hpv16", "cyt")
  }
  
  plot.predict <- model.predict(c("hpv16"), c(), pi.covars , plot.data, fit=model.fit,
                                time.points=seq(0, 10, 0.2), include.h=model.background.risk, calc.CI = T)
  
  plot.predict <- cbind(do.call(rbind,plot.predict), as.data.frame(lapply(plot.data, rep, each=51)))
  
  # size: 1000 x 600
  cr.plot <- ggplot(plot.predict, aes(x=Time, y=CR, col=names, group=names, lty=names)) + 
    geom_line(lwd=0.6) + theme_classic() + ylim(c(0, 0.8)) + 
    ylab("Cumulative risk of CIN2+") + xlab("Time (years)") + ggtitle(data.title) + 
    scale_color_brewer("Baseline HPV & cytology result", palette = "Paired" ) +
    scale_linetype_manual("Baseline HPV & cytology result", values =  c(1,1,1,1, 2,2,2,2)) + #rep(c(2, 2, 1, 1), 2)) + 
    theme(legend.position = "bottom") + guides(colour = guide_legend(title.position = "top")) + 
    scale_x_continuous(breaks = c(0, 5, 10, 15,20)) 
  
  # table of predictons
  plot.predict.table <- model.predict(c("hpv16"), c(), pi.covars, plot.data, fit=model.fit, 
                                      time.points=c(0, 3, 5, 10), include.h=model.background.risk, calc.CI = T)

  rows <- list()
  for (i in 1:length(plot.predict.table)){
    rows[[i]] <- apply(plot.predict.table[[i]], 1, 
                       function(x) paste0(round(x*100, 2)[2], " (", round(x*100, 2)[4], " - ", round(x*100, 2)[5], ")"))
  }
  matrix(unlist(rows), ncol=4, byrow=T) %>% write.csv(., preds.name)
  
  return(list(matrix(unlist(rows), ncol=4, byrow=T), cr.plot))
}

cr.func(pob_int, pob.fit.h, model.background.risk = T, "POBASCAM (The Netherlands)", "pob_preds.csv", age.covariate=T)
cr.func(italy_int, italy.fit.h2, model.background.risk=T, "NTCC Study (Italy)", "italy_preds.csv", age.covariate=T)
cr.func(slov_int, slov.fit.h, model.background.risk=F, "Slovenian HPV Prevalence study (Slovenia)", "slov_preds.csv", age.covariate=T)
cr.func(swed_int, swed.fit.h, model.background.risk=F, "Swedescreen Intervention Arm (Sweden)", "swed_preds.csv", age.covariate=F)
