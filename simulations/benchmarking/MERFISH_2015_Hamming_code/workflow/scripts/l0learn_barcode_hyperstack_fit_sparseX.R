if (!requireNamespace("L0Learn", quietly = TRUE)) {
  # No conda package provides L0Learn for win-64 (only dnachun's linux-64 build exists),
  # so install the CRAN release, which ships precompiled Windows binaries, straight into
  # this conda env's R library on first use.
  install.packages("L0Learn", repos = "https://cloud.r-project.org")
}
library(L0Learn)
library(tiff)
library(png)
library(Matrix)
library(glmnet)


cpaths = read.csv(snakemake@input[[1]])
cand_dots = read.csv(snakemake@input[[2]])
codebook = read.csv(snakemake@input[[3]])
n = dim(codebook)[2]
q = length(unique(c(as.matrix(codebook[,2:n]))))

cat("n: ", n,'; q', q, '\n')

if (nrow(cpaths) == 0){
  ranked_cpaths = data.frame(cpath=c(NA),
                             cost=c(NA),
                             gene_number=c(NA),
                             x=c(NA),
                             y=c(NA),
                             z=c(NA),
                             stringency_rank=c(NA),
                             lambda=c(NA),
                             gene_name=c(NA)
  )
  lambda_stats = data.frame(lambda=c(NA),
                            gamma=c(NA),
                            suppSize=c(NA)
  )
  write.csv(ranked_cpaths, snakemake@output[[1]])
  write.csv(lambda_stats, snakemake@output[[2]])
  write.csv(ranked_cpaths, snakemake@output[[3]])
  write.csv(data.frame(lambda=numeric(0), step=integer(0), gene=integer(0),
                       gene_number=integer(0), x=numeric(0), y=numeric(0), z=numeric(0)),
            snakemake@output[[4]])
} else{

psf = function(x, y, z){
  exp(-(x^2 + y^2 + z^2)/2.0)
}

fov_width = snakemake@config$roi_width + 2 * snakemake@config$roi_pad
i = 1

Y = array(0, dim=n*q*fov_width*fov_width)
for (img_i in 4:length(snakemake@input)){
  fname = snakemake@input[[img_i]]
  fname_split = strsplit(fname, "_")[[1]]
  r = as.numeric(fname_split[4])
  pc = as.numeric(fname_split[6])
  img = readPNG(fname)
  for (y in 1:fov_width){
    for (x in 1:fov_width){
      Y[i] = round(2^16*(img[y,x]))
      i = i + 1
    }
  }
}

n_nz_pixels_ub = 4*nrow(cpaths)*sum(psf(-10:10, 0,0) >= 0.001)^2

is = array(0, dim = n_nz_pixels_ub)
js = array(0, dim = n_nz_pixels_ub)
vals = array(0, dim = n_nz_pixels_ub)

nz_val = 1
for (j in 1:nrow(cpaths)){
  temp = cpaths[j,"cpath"]
  temp2<-substr(temp,2,nchar(temp)-1)
  temp3<-strsplit(temp2,", ")[[1]]
  dot_array <-as.numeric(temp3)
  for (dot in dot_array){
    xdot = cand_dots[dot,"x"]
    ydot = cand_dots[dot,"y"]
    zdot = cand_dots[dot,"z"]
    block = cand_dots[dot,"block"]
    pseudocolor = cand_dots[dot,"pseudocolor"]
    if (pseudocolor == 0){
      s = (q-1)*block
    } else{
      s = (q-1)*(block-1) + pseudocolor
    }

    i = 1
    for (yarr in 1:fov_width){
      for (xarr in 1:fov_width){
        psf_val = psf(xarr-xdot, yarr-ydot, 0)
        if (psf_val > 0.001){
          is[nz_val] = i + (s-1)*fov_width*fov_width
          js[nz_val] = j
          vals[nz_val] = psf_val
          nz_val = nz_val + 1
        }
        i = i + 1
      }
    }
  }
}


is = as.integer(is[vals > 0])
js = as.integer(js[vals > 0])
vals = as.vector(vals[vals > 0])
cat("The max i is ",max(i),"\n ",sep="")
cat("The max j is ",max(j),"\n",sep="")
cat("The min i is ",min(i),"\n ",sep="")
cat("The min j is ",min(j),"\n",sep="")
dimX = c(length(Y), nrow(cpaths))
print("dimX")
print(dimX)

X = sparseMatrix(i=is, j = js, x=vals, dims=dimX)
glmnet_res = glmnet(X,Y,alpha=1, lower.limits=0)
lasso_selected_coefs = coef(glmnet_res)[2:(nrow(cpaths)+1),dim(coef(glmnet_res))[2]] > 0
sum(lasso_selected_coefs)
cpaths2 = cpaths[lasso_selected_coefs,]

X2 = X[,lasso_selected_coefs]
print("dim(X2)")
print(dim(X))

reged = L0Learn.fit(X2, Y, penalty="L0", lows=0, maxSuppSize=1000, atol= 1000)


cpath_list = array(1:sum(lasso_selected_coefs))

stringency_rank = array(dim=length(cpath_list))
lambda = array(dim=length(cpath_list))
coef_at_most_stringent = array(dim=length(cpath_list))
ndecodings = dim(coef(reged))[2]


for (i in 2:ndecodings){
  
  nonzeros = which(coef(reged)[-1,i] != 0)
  if (any(coef(reged)[,i] < 0)){
    # A negative coefficient (in practice the unconstrained intercept) marks an over-fit
    # solution at this lambda. Discard it but keep ranking subsequent steps, rather than
    # aborting: breaking here drops every true barcode that first enters on this step.
    next
  }
  for (nz in nonzeros) {
    if (!is.na(nz)){
      if (nz > 0){
        if(is.na(stringency_rank[nz])) {
          stringency_rank[nz] = i-1
          lambda[nz] = print(reged)[i, "lambda"]
          coef_at_most_stringent[nz] = coef(reged)[nz+1,i]
        }
      }
    }
  }
}


cpath_selection_stringency = data.frame(cpath = cpath_list, stringency_rank = stringency_rank, lambda = lambda,coef_at_most_stringent=coef_at_most_stringent)


cpaths2$stringency_rank = stringency_rank
cpaths2$lambda = lambda

ranked_cpaths = cpaths2[complete.cases(cpaths2),]

cpath_selection_stringency = cpath_selection_stringency[complete.cases(cpath_selection_stringency ),]

n_cpath_at_stringency = array(dim=nrow(cpath_selection_stringency))
for (i in 1:nrow(cpath_selection_stringency)){
  n_cpath_at_stringency[i] = sum(cpath_selection_stringency[,"stringency_rank"] <= cpath_selection_stringency[i,"stringency_rank"])
}

ranked_cpaths$gene_name = codebook$gene_name[ranked_cpaths$gene_number]

write.csv(ranked_cpaths, snakemake@output[[1]])

ranked_cpaths[ranked_cpaths$lambda > 0.001,]

lambda_stats = print(reged)
lambda_stats$cvMeans = reged$cvMeans[[1]]
lambda_stats$cvSDs = reged$cvSDs[[1]]
write.csv(lambda_stats, snakemake@output[[2]])
write.csv(ranked_cpaths[ranked_cpaths$lambda >= 0.01,], snakemake@output[[3]])

# --- actual per-lambda support over the whole L0 path ---
# The stringency ranking above assigns each codepath its FIRST-ENTRY lambda, which assumes the L0
# support grows monotonically as lambda decreases. It does not: a codepath can leave and re-enter
# the support. To let the scorer report the best zero-FDR sensitivity found at ANY lambda, emit the
# genes actually selected (non-zero coefficient) at each lambda step.
cm = as.matrix(coef(reged))            # (1 intercept + kept codepaths) x lambda steps
reged_tab = as.data.frame(print(reged))
nsteps = min(ncol(cm), nrow(reged_tab))
perlam_parts = list()
for (L in 1:nsteps){
  supp = which(cm[-1, L] != 0)         # codepaths in the actual support at this lambda
  if (length(supp) == 0) next
  perlam_parts[[length(perlam_parts) + 1]] = data.frame(
    lambda = reged_tab[L, "lambda"], step = L,
    gene = cpaths2$gene[supp], gene_number = cpaths2$gene_number[supp],
    x = cpaths2$x[supp], y = cpaths2$y[supp], z = cpaths2$z[supp]
  )
}
perlambda = if (length(perlam_parts) > 0) do.call(rbind, perlam_parts) else
  data.frame(lambda=numeric(0), step=integer(0), gene=integer(0),
             gene_number=integer(0), x=numeric(0), y=numeric(0), z=numeric(0))
write.csv(perlambda, snakemake@output[[4]])

}
