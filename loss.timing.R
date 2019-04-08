library(data.table)
lt <- fread("loss.timing.tsv")
(loss.inc <- lt[, list(
  total.loss=min(total.loss)
), by=list(peaks)][order(peaks)])
loss.inc[, cummin := cummin(total.loss)]
loss.min <- loss.inc[total.loss==cummin]
path.dt <- data.table(penaltyLearning::modelSelection(loss.min, "total.loss", "peaks"))
R.res <- loss.min[, modelSelectionR(total.loss, peaks, seq_along(peaks))]


plot(log10(seconds) ~ log10(penalty), lt)
plot(log10(megabytes) ~ log10(penalty), lt)
plot(log10(mean.intervals) ~ log10(penalty), lt)

lt.stats <- lt[, list(
  median.intervals=median(mean.intervals),
  q25=quantile(mean.intervals, 0.25),
  q75=quantile(mean.intervals, 0.75),
  count=.N
  ), by=list(penalty=10^floor(log10(penalty)))]
library(ggplot2)
ggplot()+
  geom_ribbon(aes(
    penalty, ymin=q25, ymax=q75),
    data=lt.stats,
    alpha=0.5)+
  geom_line(aes(
    penalty, median.intervals),
    alpha=0.5,
    data=lt.stats)+
  geom_text(aes(
    penalty, median.intervals, label=count),
    data=lt.stats)+
  scale_x_log10()+
  scale_y_log10("mean intervals")
