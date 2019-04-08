library(data.table)
db.loss <- fread("db-loss.tsv")
db.path <- fread("db-path.tsv")
db.err <- db.loss[, list(
  min.err=min(remove.errors),
  max.err=max(remove.errors)
  ), by=list(prob.id, peaks)]
db.err[min.err != max.err]#TODO investigate why min.err does not equal max.err
db.err <- db.loss[order(prob.id, time.computed), .SD[1], by=list(prob.id, peaks)]
db.path.err <- db.err[db.path, on=list(prob.id, peaks)]
db.path.err[, list(models=.N), by=list(prob.id)][order(-models)]
db.path.err[, table(remove.errors-join.errors)]
db.path.err[, table(use.errors-join.errors)]
db.path.err[, table(use.errors-remove.errors)]


one <- db.path.err[prob.id==1155]
##one <- db.path.err[prob.id==1]
library(ggplot2)
ggplot()+
  scale_x_log10()+
  geom_point(aes(
    peaks, remove.errors),
    shape=1,
    data=one)

seg.var.vec <- c(
  "equality.constraints", "remove.errors", "use.errors", "join.errors")
for(seg.var in seg.var.vec){
  diff.vec <- diff(one[[seg.var]])
  data.table(one, diff.with.next=c(diff.vec, NA))
}

one.tall <- melt(
  one, measure.vars=c("remove.errors", "use.errors", "join.errors"))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_log10()+
  geom_point(aes(
    peaks, value, color=variable),
    shape=1,
    data=one.tall)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(facet ~ ., scales="free")+
  scale_x_log10()+
  geom_point(aes(
    peaks, value, color=variable),
    shape=1,
    data=data.table(facet="errors", one.tall))+
  geom_point(aes(
    peaks, log10(equality.constraints)),
    shape=1,
    data=data.table(facet="log10(equality.constraints)", one))


