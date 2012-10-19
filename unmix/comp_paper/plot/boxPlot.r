
library("clinfun")

source("R/unmix/eval.r")

# source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")

method.name = c("ep", "tp", "tpc", "pseudo")  # , "mf")
method.label = c("expectation propagation",
  "constrained pseudoinverse",
  "constr. pseudoinverse w/ correlation",
  "pseudoinverse")
#  ,"average fold-change")
method.color = c("darkgrey", "red", "purple", "blue")   # , "green")
# method.color = c("#80808080", "#80000080", "#00800080", "#00008080")

# the "actual" answers
expr = list(expr.cell = expr.cell,
  sim.cell.cor = sim.cell.cor,
  sim.expr.cell = sim.expr.cell,
  sim.expr.sym = sim.expr.sym)

# definition of when genes are "on" or "off" (somewhat hokey)
expr.on.off = list(expr.cell = expr.cell >= 2000,
  sim.cell.cor = sim.cell.cor >= 0.2,
  sim.expr.cell = sim.expr.cell >= 5,
  sim.expr.sym = sim.expr.sym >= 5)


xlim = c(1, length(num.reporters)+1)

method.offset = c(0:(length(method.name)-1)) / (length(method.name)+1)
names(method.offset) = method.name


# Plots "correlation" accuracy for several data sets w.r.t. number of reporters.
plot.correlation.accuracy = function() {
  accuracies = NULL

  # way to compute accuracy
  accuracy = function(x, x.d) {
    r = diag(cor(t(x), t(x.d)))
#    r = r[!is.na(r)]
    r
  }

  plot.it = function(data.set, data.set.name) {
    plot(0,0, xlim=xlim, ylim=c(0,1), type="n",
      main=data.set.name, xlab="Number of reporters", ylab="Correlation",
      cex.main=1.1, cex.axis=1, cex.lab=1, xaxt="n")
    axis(1, at=0.4+(1:length(num.reporters)), labels=num.reporters)
    if (data.set == "expr.cell") {
      legend("bottomright", legend=method.label, col=method.color,
        lwd=5, seg.len=1, cex=0.75)
#      mtext("a)", adj=0, cex=1)
    }

    for(i in 1:length(method.name)) {
      for(j in 1:length(num.reporters)) {
cat("data.set =", data.set, "", class(data.set), "\n")
        nr = num.reporters[[j]]
        par(new=TRUE)
        a = accuracy(expr[[data.set]],
          unmix.r[[method.name[i] ]][[data.set]][[nr]])
        boxplot(a, at = j + method.offset[[method.name[[i]] ]],
          boxwex = 1 / 3, lwd = 0.7,
          xlab="", ylab="", xaxt="n", yaxt="n",
          add=TRUE, col=method.color[[i]])
        f = data.frame(method.name = method.name[[i]],
          data.set = data.set,
          num.reporters = nr, gene = rownames(expr[[data.set]]),
          a = as.vector(a))
print(dim(accuracies))
        accuracies <<- rbind(accuracies, f)
      }
    }
  }

  plot.it("expr.cell", "Measured expression")
  plot.it("sim.cell.cor", "Simulated patterns based on correlation")
  plot.it("sim.expr.cell", "Simulated one-lineage patterns")
  plot.it("sim.expr.sym", "Simulated two-lineage symmetric patterns")

  accuracies
}

# Plots "area under curve" accuracy for several data sets
# w.r.t. number of reporters.
plot.roc.accuracy = function() {
  accuracies = NULL

  # Accuracy by AUC.
  accuracy = function(x, x.d) {
    r = auc(x, x.d)
#    r = r[!is.na(r)]
    r
  }
  ylim = c(0.5, 1)

  plot.it = function(data.set, data.set.name) {
    plot(0,0, xlim=xlim, ylim=ylim, type="n",
      main=data.set.name, xlab="Number of reporters", ylab="Area under curve",
      cex.main=1.1, cex.axis=1, cex.lab=1, xaxt="n")
    axis(1, at=0.4+(1:length(num.reporters)), labels=num.reporters)
    if (data.set == "expr.cell") {
      legend("bottomright", legend=method.label, col=method.color,
        lwd=5, seg.len=1, cex=0.75)
#      mtext("a)", adj=0, cex=1)
    }

    for(i in 1:length(method.name)) {
      for(j in 1:length(num.reporters)) {
cat("data.set =", data.set, "", class(data.set), "\n")
        nr = num.reporters[[j]]
        par(new=TRUE)
        a = accuracy(expr.on.off[[data.set]],
          unmix.r[[method.name[i] ]][[data.set]][[nr]])
        boxplot(a, at = j + method.offset[[method.name[[i]] ]],
          boxwex = 1 / 3, lwd = 0.7,
          xlab="", ylab="", xaxt="n", yaxt="n",
          add=TRUE, col=method.color[[i]])
        f = data.frame(method.name = method.name[[i]],
          data.set = data.set, num.reporters = nr,
          gene = rownames(expr.on.off[[data.set]]), a = as.vector(a))
print(dim(accuracies))
        accuracies <<- rbind(accuracies, f)
      }
    }

  }

  plot.it("expr.cell", "Measured expression")
  plot.it("sim.cell.cor", "Simulated patterns based on correlation")
  plot.it("sim.expr.cell", "Simulated one-lineage patterns")
  plot.it("sim.expr.sym", "Simulated two-lineage symmetric patterns")

  accuracies
}

plot.it = function() {
  pdf("git/unmix/comp_paper/plot/accuracyBoxplot.pdf",
    title="Boxplot of accuracy vs. number of reporters", width=11.5, height=7)
  par(mfrow=c(2,4))

  cor.accuracy = plot.correlation.accuracy()
  save(cor.accuracy, file="git/unmix/comp_paper/cor.accuracy.Rdata")

  auc.accuracy = plot.roc.accuracy()
  save(auc.accuracy, file="git/unmix/comp_paper/auc.accuracy.Rdata")

  dev.off()
}

plot.it()

