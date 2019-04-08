library(data.table)
db.path <- fread("db-path.tsv")
db.problems <- fread("db-problems.tsv")
db.loss <- fread("db-loss.tsv")
db.err <- db.loss[, list(
  min.err=min(remove.errors),
  max.err=max(remove.errors)
  ), by=list(prob.id, peaks)]
db.err[min.err != max.err]#TODO investigate why min.err does not equal max.err
db.err <- db.loss[order(prob.id, time.computed), .SD[1], by=list(prob.id, peaks)]
db.path.err <- db.err[db.path, on=list(prob.id, peaks)]
peak.type <- "remove"
for(col.name in c("errors", "fp", "fn")){
  old.name <- paste0(peak.type, ".", col.name)
  db.path.err[[col.name]] <- db.path.err[[old.name]]
}
err.not.na <- db.path.err[!is.na(errors)]
pid <- 1155
prob <- db.problems[prob.id==pid]
one <- err.not.na[prob.id==pid]

approx.target <- function
### Compute target interval for a segmentation problem. This function
### repeated calls problem.PeakSegFPOP with different penalty values,
### until it finds an interval of penalty values with minimal label
### error. The calls to PeakSegFPOP are parallelized using mclapply if
### you set options(mc.cores).
### A time limit in minutes may be specified in a file
### problem.dir/target.minutes;
### the search will stop at a sub-optimal target interval
### if this many minutes has elapsed. Useful for testing environments
### with build time limits (travis). 
(verbose=0
){
  ## Compute the label error for one penalty parameter.
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    stopifnot(length(penalty.str) == 1)
    pen.num <- as.numeric(penalty.str)
    data.table(
      iteration,
      one[min.lambda <= pen.num & pen.num <= max.lambda])
  }  
  error.list <- list()
  next.pen <- c(0, Inf)
  iteration <- 0
  last.target.vec <- c(-Inf, Inf)
  target.result.list <- list()
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    error.list[next.str] <- lapply(next.str, getError)
    error.dt <- do.call(rbind, error.list)[order(penalty)]
    if(!is.numeric(error.dt$penalty)){
      stop("penalty column is not numeric -- check loss in _loss.tsv files")
    }
    error.dt[, errors := fp+fn]
    if(verbose)print(
      error.dt[,.(penalty, peaks, equality.constraints, fp, fn, errors)])
    unique.peaks <- error.dt[, data.table(
      .SD[which.min(iteration)],
      penalties=.N
    ), by=list(peaks)]
    path.dt <- data.table(penaltyLearning::modelSelection(
      unique.peaks, "total.loss", "peaks"))
    path.dt[, next.pen := max.lambda]
    path.dt[, already.computed := next.pen %in% names(error.list)]
    path.dt[, no.next := c(diff(peaks) == -1, NA)]
    path.dt[, done := already.computed | no.next]
    path.dt[, is.min := errors==min(errors)]
    path.dt[, min.err.interval := cumsum(ifelse(
      c(is.min[1], diff(is.min))==1, 1, 0))]
    other.candidates <- path.dt[which(0<diff(fn) & diff(fp)<0)]
    interval.dt <- path.dt[is.min==TRUE, {
      i <- if(1 == .N || 0 == errors[1]){
        ## No middle candidate if there is only one model in the
        ## interval, or if there are no errors.
        NA
      }else{
        d <- data.table(
          i=1:(.N-1),
          ## do not attempt to explore other.candidates -- try
          ## different ones!
          is.other=next.pen[-.N] %in% other.candidates$next.pen,
          dist=diff(max.log.lambda)+diff(min.log.lambda),
          done=done[-.N])
        d[is.other==FALSE & done==FALSE, i[which.max(dist)]]
      }
      if(length(i)==0)i <- NA
      data.table(
        min.lambda=min.lambda[1],
        min.log.lambda=min.log.lambda[1],
        mid.lambda=max.lambda[i],
        max.lambda=max.lambda[.N],
        max.log.lambda=max.log.lambda[.N],
        log.size=max.log.lambda[.N]-min.log.lambda[1]
        )
    }, by=list(min.err.interval)]
    largest.interval <- interval.dt[which.max(log.size)]
    target.vec <- largest.interval[, c(min.log.lambda, max.log.lambda)]
    diff.target.vec <- target.vec-last.target.vec
    last.target.vec <- target.vec
    target.result.list[[paste(iteration)]] <- largest.interval[, data.table(
      iteration,
      min.log.lambda,
      max.log.lambda)]
    target.lambda <- largest.interval[, c(min.lambda, max.lambda)]
    error.candidates <- path.dt[next.pen %in% target.lambda]
    stopping.candidates <- rbind(error.candidates, other.candidates)[done==FALSE]
    next.pen <- if(nrow(stopping.candidates)){
      lambda.vec <- interval.dt[, c(min.lambda, mid.lambda, max.lambda)]
      interval.candidates <- path.dt[next.pen %in% lambda.vec][done==FALSE]
      unique(rbind(stopping.candidates, interval.candidates)$next.pen)
    }
  }#while(!is.null(pen))
  list(
    target=target.vec,
    target.iterations=do.call(rbind, target.result.list),
    models=error.dt)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), target.iterations is a
### data.table with target intervals as a function of iteration,
### models is a data.table with one row per model for which the label
### error was computed.
}

target.list <- approx.target()

interval.dt <- penaltyLearning::targetIntervals(one, "prob.id")

mean.seconds <- mean(one$seconds)
iterations.ord <- target.list$models[, list(
  seconds=mean.seconds*.N,
  models=.N
  ), by=list(iteration)][target.list$target.iterations, on=list(iteration)][order(iteration)]
iterations.ord[, cum.seconds := cumsum(seconds)]

seg.var.vec <- c("errors", "log10.peaks", "fp", "fn")
one[, log10.peaks := log10(peaks)]
seg.dt.list <- list()
for(seg.var in seg.var.vec){
  value.vec <- one[[seg.var]]
  diff.after <- diff(value.vec)
  change.vec <- which(diff.after!=0)
  end.vec <- c(change.vec, nrow(one))
  start.vec <- c(1, change.vec+1)
  seg.dt.list[[seg.var]] <- data.table(
    variable=seg.var,
    value=value.vec[start.vec],
    min.log.lambda=one$min.log.lambda[start.vec],
    max.log.lambda=one$max.log.lambda[end.vec])
}
seg.dt <- do.call(rbind, seg.dt.list)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_log10()+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value),
    data=seg.dt[variable=="errors"])

seg.dt[, facet := fac(ifelse(
  variable %in% c("fp", "fn"), "errors", variable))]
tlev <- "iterations"
fac <- function(x){
  factor(x, c("log10.peaks", "errors", tlev))
}
seg.size <- 2
small.size <- 1.5
efac <- function(x){
  factor(x, c("approx", "errors", "fp", "fn"))
}
seg.dt[, Error.type := efac(variable)]
interval.color <- "violet"
peak.range <- one[errors==0][peaks %in% range(peaks)]
peak.range[, vjust := c(-0.5, 1.5)]
peak.range[, pen.size := c("smallest", "largest")]
approx.path <- penaltyLearning::modelSelection(
  target.list$models[, .(errors, total.loss, peaks)],
  "total.loss", "peaks")
fullpath.dt <- one[, data.table(
  days=sum(seconds)/60/60/24,
  min.log.lambda=peak.range$min.log.lambda[1],
  max.log.lambda=peak.range$max.log.lambda[2],
  models=.N,
  facet=fac("errors"),
  data=prob$bedGraph.lines)]
approx.info <- iterations.ord[, data.table(
  minutes=max(cum.seconds)/60,
  max.iterations=max(iteration),
  min.log.lambda=target.list$target[1],
  max.log.lambda=target.list$target[2],
  models=sum(models),
  facet=fac(tlev),
  data=prob$bedGraph.lines)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(facet ~ ., scales="free")+
  scale_x_continuous("log(penalty)")+
  geom_text(aes(
    10, 60, 
    label=sprintf("%.0f days to compute exact target interval, (%.4f, %.4f)
by computing path of %d models for
one chrom subset with %d data", days, min.log.lambda, max.log.lambda, models, data)),
color=interval.color,
hjust=0,
data=fullpath.dt)+
  geom_text(aes(
    10, max.iterations, 
    label=sprintf("%.0f minutes to compute approximate target interval, (%.4f, %.4f)
by computing approximate error function using %d models/penalties ", minutes, min.log.lambda, max.log.lambda, models)),
hjust=0,
vjust=1,
data=approx.info)+
  geom_vline(aes(
    xintercept=log.penalty),
    color=interval.color,
    data=data.table(log.penalty=target.list$target))+
  coord_cartesian(xlim=c(3,13))+
  geom_hline(aes(
    yintercept=log10(peaks)),
    color=interval.color,
    data=data.table(peak.range, facet=fac("log10.peaks")))+
  geom_text(aes(
    -Inf, log10(peaks),
    vjust=vjust,
    label=paste0(
      " ", peaks, " peaks in minimum error model with ",
      pen.size, " penalty")),
    hjust=0,
    color=interval.color,
    data=data.table(peak.range, facet=fac("log10.peaks")))+
  geom_segment(aes(
    min.log.lambda, value, 
    xend=max.log.lambda, yend=value),
    size=seg.size,
    data=seg.dt[variable=="log10.peaks"])+
  scale_color_manual(values=c(
    errors="black",
    fp="red",
    approx="grey50",
    fn="deepskyblue"))+
  scale_size_manual(values=c(
    errors=seg.size,
    fp=small.size,
    approx=1,
    fn=small.size))+
  geom_segment(aes(
    min.log.lambda, value, color=Error.type, size=Error.type,
    xend=max.log.lambda, yend=value),
    data=seg.dt[variable!="log10.peaks"])+
  ylab("")+
  geom_segment(aes(
    min.log.lambda, errors, color=Error.type, size=Error.type,
    xend=max.log.lambda, yend=errors),
    data=data.table(approx.path, facet=fac("errors"), Error.type=efac("approx")))+
  geom_point(aes(
    log(penalty), iteration),
    shape=21,
    fill="white",
    color="black",
    size=3,
    data=data.table(target.list$models, facet=fac(tlev)))+
  geom_segment(aes(
    min.log.lambda, iteration,
    xend=max.log.lambda, yend=iteration),
    size=seg.size,
    data=data.table(iterations.ord, facet=fac(tlev)))+
  geom_text(aes(
    3.5, iteration, label=sprintf(
      "  iteration=%d",
      iteration)),
    hjust=1,
    data=data.table(iterations.ord, facet=fac(tlev)))+
  geom_text(aes(
    4.25, iteration, label=sprintf(
      "  penalties=%d",
      cumsum(models))),
    hjust=1,
    data=data.table(iterations.ord, facet=fac(tlev)))+
  geom_text(aes(
    5, iteration, label=sprintf(
      "  minutes=%.1f",
      cum.seconds/60)),
    hjust=1,
    data=data.table(iterations.ord, facet=fac(tlev)))
png("figure-approx-target.png", 17, 8, units="in", res=100)
print(gg)
dev.off()
