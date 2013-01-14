
library("clinfun")

source("R/unmix/eval.r")

source("git/plot/label_panel.r")

source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")

method.name = c("ep", "tp", "tpc", "pseudo")  # , "mf")
method.label = c("expectation propagation",
  "constrained pseudoinverse",
  "constrained pseudoinverse with correlation",
  "pseudoinverse")

# method.label = c("EP", "CPseudo", "CPseudoCorr", "Pseudo")


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
  ylim=c(0,1.01)
  # way to compute accuracy
  accuracy = function(x, x.d) {
    r = diag(cor(t(x), t(x.d)))
#    r = r[!is.na(r)]
    r
  }

  plot.it = function(data.set, data.set.name) {
    plot(0,0, xlim=xlim, ylim=ylim, type="n",
      main=data.set.name, xlab="Number of reporters", ylab="Correlation",
      cex.main=1.6, cex.axis=1.3, cex.lab=1.4, xaxt="n")
    axis(1, at=0.4+(1:length(num.reporters)), cex.axis=1.3, labels=num.reporters)
    if (data.set == "expr.cell") {
#      legend("bottomright", legend=method.label, col=method.color,
#        lwd=5, seg.len=1, cex=0.75)
#      mtext("a)", adj=0, cex=1)
      label.panel("b)")
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

    # significance stars
    if (data.set %in% c("sim.expr.cell", "sim.expr.sym")) {
      par(new=TRUE)
      plot(3, 1.01,
        xlim=xlim, ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n",
        pch=8, cex=1.25)
    }


    }
  }

  plot.it("expr.cell", "Measured expression")
  plot.it("sim.cell.cor", "Patterns based on correlation")
  plot.it("sim.expr.cell", "One-lineage patterns")
  plot.it("sim.expr.sym", "Two-lineage patterns")

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
  ylim = c(0.5, 1.025)

  plot.it = function(data.set, data.set.name) {
    plot(0,0, xlim=xlim, ylim=ylim, type="n",
      main=data.set.name, xlab="Number of reporters", ylab="Area under curve",
      cex.main=1.6, cex.axis=1.3, cex.lab=1.4, xaxt="n")
    axis(1, at=0.4+(1:length(num.reporters)), cex.axis=1.3, labels=num.reporters)
    if (data.set == "expr.cell") {
#      legend("bottomright", legend=method.label, col=method.color,
#        lwd=5, seg.len=1, cex=0.89)
      label.panel("a)")
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

      # significance stars
      if (data.set %in% c("expr.cell", "sim.expr.cell", "sim.expr.sym")) {
        par(new=TRUE)
        plot(3, 1.025,
          xlim=xlim, ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n",
          pch=8, cex=1.25)
      }
    }

  }

  plot.it("expr.cell", "Measured expression")
  plot.it("sim.cell.cor", "Patterns based on correlation")
  plot.it("sim.expr.cell", "One-lineage patterns")
  plot.it("sim.expr.sym", "Two-lineage patterns")

  accuracies
}

plot.it = function() {
  pdf("git/unmix/comp_paper/plot/accuracyBoxplot.pdf",
    title="Boxplot of accuracy vs. number of reporters", width=11.5, height=8)
#  par(mfrow=c(2,4))

  mat = matrix(1:12, nrow=3, byrow=TRUE)
  mat[3,] = c(9,9,9,10)
  layout(mat, heights=c(1,1,0.52))

  auc.accuracy = plot.roc.accuracy()
  save(auc.accuracy, file="git/unmix/comp_paper/auc.accuracy.Rdata")

  cor.accuracy = plot.correlation.accuracy()
  save(cor.accuracy, file="git/unmix/comp_paper/cor.accuracy.Rdata")

  plot.new()
  par(mar=c(0,0,0,0)+0.2)

  legend("topleft", legend=method.label, col=method.color,
    lwd=10, seg.len=0.3, cex=1.2, y.intersp=2.5)

  dev.off()

#  embedFonts("git/unmix/comp_paper/plot/accuracyBoxplot1.pdf",
#    outfile="git/unmix/comp_paper/plot/accuracyBoxplot.pdf")
}

plot.it()

