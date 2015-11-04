source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/ToolBox_PERCH.R")
source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/Sparse_Models.R")
library(doMC)
library(foreach)
library(MASS)


num_core = min(detectCores(), 12)
registerDoMC(num_core)

load("SimDat_K5n300.RData")

# ------------------------------------------------------------------------------------------ #
# call LattE function integrate from R
#GenMatrices(K=5)$Lmatrix

Lmat = GenMatrices.Sparse(K=5, Smax=3)
MuMat = Lmat$MuMat
PiMat = Lmat$PiMat
J1 = Lmat$J1
mu.tmp = c(0.258,0.157,0.165,0.288,0.541)
pi.tmp = c(0.103,0.513,0.256,0.128)

d = ncol(MuMat)
Resolution = 1
ScaleUp = 1000

Equality = round(cbind(c(mu.tmp, pi.tmp[-(1:2)])*ScaleUp, -rbind(MuMat,PiMat[-1,]))*Resolution,0)
Inequality = round(cbind(rep(1,d)*ScaleUp, -diag(1,d))*Resolution,0)
# -------------------------------------------------------------- #

b_A = rbind(Equality, Inequality)
m = nrow(b_A)


n_eq = 7
b_A = rbind(Equality[1:n_eq,], Inequality)
m = nrow(b_A)

b_A = Inequality/Resolution
m = nrow(b_A)
# -------------------------------------------------------------- #
# write constraint file
sink("Constraints.txt")
cat(paste(m,d+1,"\n"))
for (i in 1:m)
{
    cat(c(paste(b_A[i,]),"\n"))     
}
cat(c("linearity",paste(c(n_eq,1:n_eq)),"\n"))
cat(c("nonnegative",paste(c(d,1:d)),"\n"))
sink()

# -------------------------------------------------------------- #
# write integrand file
Integrand = diag(rep(1,d))
sink("Integrand.txt")
cat(c("[[1,[1,[", paste(Integrand[i,],collapse=","), "]]]]","\n"))
sink()


PathToLattE = "/home/bst/student/ddeng/latte/bin/integrate"

int_command = "~/latte/bin/integrate --valuation=integrate --triangulate --linear-forms=Integrand.txt Constraints.txt "

vol_command = "~/latte/bin/integrate --valuation=volume --triangulate Constraints.txt"
vol_command = "~/latte/bin/integrate --valuation=volume --cone-decompose Constraints.txt"

system("cat Constraints.txt")
system(vol_command, intern = T)
system(int_command, intern = T)


