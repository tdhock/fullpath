library(data.table)
db.loss <- fread("db-loss.tsv")
path <- db.loss[, {
  print(prob.id)
  loss.dt <- data.table(.SD)
  (loss.inc <- loss.dt[, list(
    total.loss=total.loss[which.min(time.computed)]
  ), by=list(peaks)][order(peaks)])
  loss.inc[, cummin := cummin(total.loss)]
  loss.min <- loss.inc[total.loss==cummin]
  path.dt <- data.table(penaltyLearning::modelSelection(
    loss.min, "total.loss", "peaks"))
  path.dt[, max.computed := paste(max.lambda) %in% loss.dt$pen.str]
  path.dt[, no.next := c(diff(peaks) == -1, NA)]
  path.dt[, done := max.computed | no.next]
  print(path.dt[, table(no.next, max.computed)])
  path.dt
}, by=list(prob.id)]
