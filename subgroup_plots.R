library(ggplot2)
library(survival)
library(geomtextpath)

subgroup.plot.func <- function(data, age=T){
  
  predict.KM.plot <- function(data) {
    # function to fit KM non-parametric curve and our model for the same time points
    fit <-  model.fit(c(), c(), c(), data, silent = F, include.h = F, short.epsilon = 0.1, epsilon=1e-5)
    surv.fit <- survfit(Surv(data$left, data$right, type="interval2")~1)
    times <- surv.fit$time
    prob.surv <- 1- surv.fit$surv
    prob.model <- model.predict(c(), c(), c(), fit = fit, data=data[1,], time.points=times, calc.CI = F, include.h = F)[[1]]
    return(data.frame(times, prob.surv, prob.model))
  }
  
  dat_hpv1 <- predict.KM.plot(data[data$hpv16==1,])
  dat_hpv0 <- predict.KM.plot(data[data$hpv16==0,])
  dat_hpv <- vctrs::vec_c(dat_hpv1, dat_hpv0)
  dat_hpv$label <- c(rep("HPV16+", dim(dat_hpv1)[1]), rep("HPV16-", dim(dat_hpv0)[1]))
  dat_hpv$source <- "HPV16"

  dat_cyt1 <- predict.KM.plot(data[data$cyt_abnormal ==1,])
  dat_cyt0 <- predict.KM.plot(data[data$cyt_abnormal==0,])
  dat_cyt <- vctrs::vec_c(dat_cyt1, dat_cyt0)
  dat_cyt$label <- c(rep("Abnormal", dim(dat_cyt1)[1]), rep("Normal", dim(dat_cyt0)[1]))
  dat_cyt$source <- "Cytology"
  
  if (age==T){ 
    # age not estimated for sweden because all women are younger than 40
    dat_age1 <- predict.KM.plot(data[data$age==1,])
    dat_age0 <- predict.KM.plot(data[data$age==0,])
    dat_age <- vctrs::vec_c(dat_age1, dat_age0)
    dat_age$label <- c(rep(">40", dim(dat_age1)[1]), rep("<40", dim(dat_age0)[1]))
    dat_age$source <- "Age"
    
    dat_plot <- rbind(dat_hpv, dat_age, dat_cyt)
    dat_plot$source <- factor(dat_plot$source, levels=c("HPV16", "Age", "Cytology"))
  }else{
    dat_plot <- rbind(dat_hpv, dat_cyt)
    dat_plot$source <- factor(dat_plot$source, levels=c("HPV16", "Cytology"))
  }
  # size: 1000 x 350!!
  subgroup_plot <- ggplot(dat_plot, aes(x=times, y=CR, group=label)) + theme_bw() + facet_wrap(~source) +
    geom_line(aes(x=times, y=prob.surv), col='grey') +# geom_line(aes(x=times, y=CR), col='black', lwd=1.2) +
    xlab("Time (years)") + ylab(expression("Cumulative risk of CIN2+")) + ylim(c(0, 0.7)) +
    scale_hjust_manual(values = rep(0.9, 7)) +geom_textline(aes(label = label, hjust = label), lwd=1.2) + 
    scale_x_continuous(breaks=c(0, 5, 10, 15, 20)) + 
    theme(text = element_text(size = 20), panel.grid = element_blank()) 
  
  print(subgroup_plot)
  
  return(list(dat_plot, subgroup_plot))
}

pob_subgroup <- subgroup.plot.func(pob_int)  
italy_subgroup <- subgroup.plot.func(italy_int)  
slov_subgroup <- subgroup.plot.func(slov_int)  
swed_subgroup  <- subgroup.plot.func(swed_int, age=F)  
          