
if(getRversion() >= "2.15.1")
  utils::globalVariables(c("ci","cost","cost.se","event.prob","n.patients","n.patients.se",
                           "n.screened","n.screened.se","reduction","reduction.se"))

surv_plot_enrichment <- function (x, km.quantiles = c(0,0.25,0.5,0.75),
                                  km.range = NULL, alt.color = NULL){
  plot.error.bar <- as.numeric(!is.null(x$event.prob.se))
  reduc.error.bar <- as.numeric(!is.null(x$cost.reduction.se))

  if (is.null(x$event.prob.se))
    x$event.prob.se <- x$n.patients.se <- x$n.screened.se <- x$cost.se <-
      x$cost.reduction.se <- rep(0,length(x$selected.biomarker.quantiles))
  if (x$acc.fu)
    x$end.of.trial <- x$a+x$f
  if (x$acc.fu==FALSE){
    x$cost.reduction.se <- matrix(0,nrow=length(x$selected.biomarker.quantiles),ncol=length(x$end.of.trial))
    if (is.null(x$cost.se))
      x$cost.se <- matrix(0,nrow=length(x$selected.biomarker.quantiles),ncol=length(x$end.of.trial))
  }

  end.of.trial <- x$end.of.trial
  len.pos <- "bottom"
  if (length(end.of.trial)==1) len.pos <- "none"
  comp.cost <- !is.null(x$cost)

  #colors
  gg_color_hue <- function(n) {
    if (n==1) return("black")
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if (!is.null(alt.color) & length(alt.color)==length(end.of.trial))
    gg_color_hue <- function(n) alt.color

  # 1. K-M curve ################
  cols <- gray.colors(length(km.quantiles)+3)
  km.all <- survfit(x$response~1, error="greenwood")
  dat <- as.data.frame(seq(0,max(km.all$time),by=max(km.all$time)/500))
  if (!is.null(km.range)){
    km.range<=max(km.all$time)
    dat <- as.data.frame(seq(0,km.range,by=km.range/500))
  }
  colnames(dat) <- "time"
  for (j in 1:length(km.quantiles)){
    q <- quantile(x$biomarker,prob=km.quantiles[j])
    sobj <- x$response[x$biomarker>=q]
    km <- survfit(sobj~1,error="greenwood")
    survfun <- stepfun(km$time, c(1, km$surv))
    dat <- cbind(dat, survfun(dat[,1]))
    colnames(dat)[j+1] <- paste(j,"surv",sep=".")
  }
  dat <- reshape(dat, direction = 'long', timevar = 'level.enrichment',
                 varying=list(grep("surv", colnames(dat), value=TRUE)),
                 times = as.character(km.quantiles),
                 v.names = c("surv"),
                 idvar='time')
  g <- ggplot(dat,aes(x=time, y=surv, colour=level.enrichment)) +
    geom_line(size=1) + ylim(0,1) +
    labs(title ="Kaplan-Meier survival curves",
         x = "time", y = "survival estimate", color = "enrichment level") +
    scale_color_manual(labels = as.character(km.quantiles), values = cols[1:j]) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
  if (x$acc.fu==FALSE){
    for (k in 1:length(x$end.of.trial))
      g <- g + geom_vline(xintercept = x$end.of.trial[k], colour = gg_color_hue(length(x$end.of.trial))[k])
  }
  if (x$acc.fu==TRUE){
    g <- g + geom_vline(xintercept = x$f, colour = cols[j])
    #g <- g + geom_vline(xintercept = x$a/2+x$f, colour = cols[j])
    g <- g + geom_vline(xintercept = x$a+x$f, colour = cols[j])
  }

  # 2. event rate and sd ###############
  dat <- as.data.frame(cbind(x$selected.biomarker.quantiles, x$event.prob, x$event.prob.se))
  colnames(dat)[1] <- "level.enrichment"
  for (j in 1:length(x$end.of.trial)){
    colnames(dat)[1+j] <- paste(j,"prob", sep=".")
    colnames(dat)[(1+length(x$end.of.trial)+j)] <- paste(j, "sd", sep=".")
  }
  dat <- reshape(dat, direction = 'long', timevar = 'end.of.trial',
                 varying=list(grep("prob", colnames(dat), value=TRUE), grep("sd", colnames(dat), value=TRUE)),
                 times = as.character(seq(1,length(x$end.of.trial))),
                 v.names = c("event.prob", "event.prob.se"),
                 idvar='level.enrichment')
  dat$ci <- 1.96*dat$event.prob.se
  #pd <- position_dodge(0) # move them .05 to the left and right
  pd <- position_jitter()
  if (x$acc.fu == FALSE){
      tt <- "Event rate"
  }
  if (x$acc.fu == TRUE)
    tt <- "Average event rate"
  g2 <- ggplot(dat, aes(x=100*level.enrichment, y=event.prob, colour=end.of.trial)) +
    geom_errorbar(aes(ymin=event.prob-ci, ymax=event.prob+ci),
                  width=.05*length(x$end.of.trial)*sd(x$selected.biomarker.quantiles*100)*plot.error.bar) +
                  #position=pd) +
    geom_line() +#position=pd) +
    geom_point() +#position=pd) +
    labs(title = tt, x = "level of enrichment", y = "event rate") +
    labs(color = "end of trial") +
    scale_color_manual(labels = as.character(x$end.of.trial), values = gg_color_hue(length(x$end.of.trial))) +
    ylim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5), legend.position=len.pos)

  # 3. total sample size
  dat <- as.data.frame(cbind(x$selected.biomarker.quantiles, x$n.patients, x$n.patients.se))
  colnames(dat)[1] <- "level.enrichment"
  for (j in 1:length(x$end.of.trial)){
    colnames(dat)[1+j] <- paste(j,"num", sep=".")
    colnames(dat)[(1+length(x$end.of.trial)+j)] <- paste(j, "sd", sep=".")
  }
  dat <- reshape(dat, direction = 'long', timevar = 'end.of.trial',
                 #varying = colnames(dat)[-1],
                 varying=list(grep("num", colnames(dat), value=TRUE), grep("sd", colnames(dat), value=TRUE)),
                 times = as.character(seq(1,length(end.of.trial))),
                 v.names = c("n.patients", "n.patients.se"),
                 idvar='level.enrichment')
  g3 <- ggplot(dat, aes(x=level.enrichment*100, y=n.patients, colour=end.of.trial)) +
    geom_errorbar(aes(ymin=n.patients-1.96*n.patients.se, ymax=n.patients+1.96*n.patients.se),
                  width=.05*length(end.of.trial)*sd(x$selected.biomarker.quantiles*100)*plot.error.bar) +
                  #position=pd) +
    geom_line() +#position=pd) +
    geom_point() +#position=pd) +
    expand_limits(y=0) +
    labs(title ="Clinical trial sample size",
         x = "level of enrichment", y = "total sample size", color = "end of trial") +
    scale_color_manual(labels = as.character(x$end.of.trial), values = gg_color_hue(length(end.of.trial))) +
    theme(plot.title = element_text(hjust = 0.5), legend.position=len.pos)

  # 4. total patients screened
  dat <- as.data.frame(cbind(x$selected.biomarker.quantiles, x$n.screened, x$n.screened.se))
  colnames(dat)[1] <- "level.enrichment"
  for (j in 1:length(x$end.of.trial)){
    colnames(dat)[1+j] <- paste(j,"num", sep=".")
    colnames(dat)[(1+length(x$end.of.trial)+j)] <- paste(j, "sd", sep=".")
  }
  dat <- reshape(dat, direction = 'long', timevar = 'end.of.trial',
                 varying=list(grep("num", colnames(dat), value=TRUE), grep("sd", colnames(dat), value=TRUE)),
                 times = as.character(seq(1,length(end.of.trial))),
                 v.names = c("n.screened", "n.screened.se"),
                 idvar='level.enrichment')
  g4 <- ggplot(dat, aes(x=level.enrichment*100, y=n.screened, colour=end.of.trial)) +
    geom_errorbar(aes(ymin=n.screened-1.96*n.screened.se, ymax=n.screened+1.96*n.screened.se),
                  width=.05*length(end.of.trial)*sd(x$selected.biomarker.quantiles*100)*plot.error.bar) +
                  #position=pd) +
    geom_line() +#position=pd) +
    geom_point() +#position=pd) +
    labs(title ="Number of patients screened",
         x = "level of enrichment", y = "total # screened", color = "end of trial") +
    scale_color_manual(labels = as.character(x$end.of.trial), values = gg_color_hue(length(end.of.trial))) +
    expand_limits(y=0) +
    theme(plot.title = element_text(hjust = 0.5), legend.position=len.pos)

  # plots for cost
  g5 <- g6 <- NULL

  if (comp.cost){
    # 5. total costs for screening and trial
    dat <- as.data.frame(cbind(x$selected.biomarker.quantiles, x$cost, x$cost.se))
    colnames(dat)[1] <- "level.enrichment"
    for (j in 1:length(x$end.of.trial)){
      colnames(dat)[1+j] <- paste(j,"num", sep=".")
      colnames(dat)[(1+length(x$end.of.trial)+j)] <- paste(j, "sd", sep=".")
    }
    dat <- reshape(dat, direction = 'long', timevar = 'end.of.trial',
                   varying=list(grep("num", colnames(dat), value=TRUE), grep("sd", colnames(dat), value=TRUE)),
                   times = as.character(seq(1,length(end.of.trial))),
                   v.names = c("cost", "cost.se"),
                   idvar='level.enrichment')
    g5 <- ggplot(dat, aes(x=level.enrichment*100, y=cost, colour=end.of.trial)) +
      geom_errorbar(aes(ymin=cost-1.96*cost.se, ymax=cost+1.96*cost.se),
                    width=.05*length(end.of.trial)*sd(x$selected.biomarker.quantiles*100)*plot.error.bar) +
      #position=pd) +
      geom_line() +#position=pd) +
      geom_point() +#position=pd) +
      labs(title ="Total screening + trial cost",
           x = "level of enrichment", y = "total cost", color = "end of trial") +
      scale_color_manual(labels = as.character(x$end.of.trial), values = gg_color_hue(length(end.of.trial))) +
      expand_limits(y=0) +
      theme(plot.title = element_text(hjust = 0.5), legend.position=len.pos)

    # 6. % reduction in total cost
    dat <- as.data.frame(cbind(x$selected.biomarker.quantiles, x$cost.reduction, x$cost.reduction.se))
    colnames(dat)[1] <- "level.enrichment"
    for (j in 1:length(x$end.of.trial)){
      colnames(dat)[1+j] <- paste(j,"num", sep=".")
      colnames(dat)[(1+length(x$end.of.trial)+j)] <- paste(j, "sd", sep=".")
    }
    dat <- reshape(dat, direction = 'long', timevar = 'end.of.trial',
                   varying=list(grep("num", colnames(dat), value=TRUE), grep("sd", colnames(dat), value=TRUE)),
                   times = as.character(seq(1,length(end.of.trial))),
                   v.names = c("reduction", "reduction.se"),
                   idvar='level.enrichment')
    g6 <- ggplot(dat, aes(x=level.enrichment*100, y=reduction*100, colour=end.of.trial)) +
      geom_errorbar(aes(ymin=reduction*100-196*reduction.se, ymax=reduction*100+196*reduction.se),
                    width=.05*length(end.of.trial)*sd(x$selected.biomarker.quantiles*100)*reduc.error.bar) +
      #position=pd) +
      geom_line() +#position=pd) +
      geom_point() +#position=pd) +
      labs(title ="Reduction (%) in total cost",
           x = "level of enrichment", y = "% reduction in cost", color = "end of trial") +
      scale_color_manual(labels = as.character(x$end.of.trial), values = gg_color_hue(length(x$end.of.trial))) +
      expand_limits(y=0) +
      theme(plot.title = element_text(hjust = 0.5), legend.position=len.pos)
  }

  # plot and return
  if (comp.cost)
    summary <- arrangeGrob(grid.arrange(g,g2,g3,g4,g5,g6, nrow=3))
  if (!comp.cost)
    summary <- arrangeGrob(grid.arrange(g,g2,g3,g4, nrow=2))
  return(list(km.plot=g, prob.plot=g2, ss.plot=g3,
              screen.plot=g4, cost.plot=g5, reduction.cost.plot=g6,
              summary=summary))
}

