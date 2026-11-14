# Figure 1 and Figure S1. These drawings describe group membership, not enrollment sizes.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, args[1] == "--output-dir")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ink <- "#293c48"
fills <- c(trial="#f1f3f3", treated="#e8f2ea", internal="#e7eff9", external="#f3edf6")
borders <- c(trial="#677780", treated="#286642", internal="#315f90", external="#704789")
panel <- function() {
  plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs="i", yaxs="i")
}
node <- function(x, y, label, type="trial", w=.38, dotted=FALSE) {
  rect(x-w/2, y-.058, x+w/2, y+.058, col=if (dotted) "white" else fills[type],
       border=borders[type], lty=if (dotted) 3 else 1)
  text(x,y,label,cex=.85,col=ink)
}
arrow <- function(x,y,x1,y1) arrows(x,y,x1,y1,length=.07,col=ink)
grDevices::cairo_pdf(file.path(out,"trial-designs.pdf"), width=10,height=3.65)
par(mfrow=c(1,3), mar=c(.15,.4,.4,.4))
for (design in c("RCT","HCT","ECT")) {
  panel()
  text(.02,.96,switch(design,RCT="(a) RCT",HCT="(b) HCT",ECT="(c) ECT"),adj=0,font=2)
  node(.5,.82,if (design=="ECT") "Single-arm trial cohort" else "Trial cohort",w=.77)
  if (design!="ECT") text(.5,.685,"Randomize",cex=.82)
  lines(c(.5,.5,.24),c(if (design=="ECT") .762 else .65,.625,.625),col=ink)
  arrow(.24,.625,.24,.57)
  lines(c(.5,.77),c(.625,.625),col=ink,lty=if (design=="ECT") 3 else 1)
  if (design=="ECT") segments(.77,.625,.77,.57,col=ink,lty=3) else arrow(.77,.625,.77,.57)
  node(.24,.505,"Treated","treated")
  node(.77,.505,if (design=="ECT") "Trial control\nunobserved" else "Internal controls",
       if (design=="ECT") "trial" else "internal",dotted=design=="ECT")
  if (design=="HCT") {
    node(.24,.31,"External controls","external")
    node(.77,.31,"Hybrid control\nmean")
    arrow(.43,.31,.575,.31)
    arrow(.77,.447,.77,.37)
    lines(c(.05,.025,.025,.5),c(.505,.505,.19,.19),col=ink)
  } else lines(c(.24,.24,.5),c(.447,.19,.19),col=ink)
  if (design=="ECT") node(.77,.31,"External controls","external")
  lines(c(.77,.77,.5),c(if (design=="RCT") .447 else .252,.19,.19),col=ink)
  arrow(.5,.19,.5,.148)
  node(.5,.085,switch(design,RCT="Randomized comparison",HCT="Hybrid comparison",
                       ECT="External comparison"),w=.85)
}
invisible(dev.off())

# Hypothetical records. The symbols on the right denote observed outcomes in each group.
dat <- data.frame(ID=sprintf("P%02d",1:6), S=c(1,1,1,1,0,0), A=c(1,1,0,0,0,0),
                  X1=c(".20",".70",".30",".80",".25",".65"), X2=c(0,1,1,0,0,1),
                  Y=c("8.4","10.1","7.2","8.5","6.9","8.0"))
grDevices::cairo_pdf(file.path(out,"hct_observed_data.pdf"),width=7.3,height=3.05)
par(mar=c(.1,.1,.2,.1)); plot.new(); plot.window(c(0,7.7),c(-.7,7.15),xaxs="i",yaxs="i")
for (i in 0:6) {
  y <- 6.55-i
  type <- if (i==0) "trial" else rep(c("treated","internal","external"),each=2)[i]
  rect(.05,y-.46,6.25,y+.47,col=fills[type],border="white")
  for (j in 1:6) text(j-.35,y,if (i==0) c("ID","S","A","X1","X2","Y")[j] else as.character(dat[i,j]),
                        font=if (i==0 || j==1) 2 else 1, cex=.95,col=ink)
}
for (k in 1:3) {
  y <- 5.05-2*(k-1)
  lines(c(6.36,6.48,6.48,6.36), y+c(.85,.85,-.85,-.85),col=borders[k+1],lwd=2)
  text(6.65,y,c(expression(Y[1]),expression(Y["1,0"]),expression(Y[0]))[k],
       adj=0,col=borders[k+1],cex=1.1)
}
legend("bottom", c("Trial treated","Internal control","External control"),horiz=TRUE,
       fill=fills[-1],border=borders[-1],bty="n",cex=.9)
invisible(dev.off())
