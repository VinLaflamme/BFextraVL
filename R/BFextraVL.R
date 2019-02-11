#' Model Matrix for a posterior sample chain from a BF
#' 
#' Returns the model matrix of the posterior sample chain of a BFmcmc object. Can be used to compute the predicted means.
#' 
#' @param factorLevels Numerical vector which stores the number of levels for each factor in the design. The order of the factors should follow that of the elements composing the name of the highest level interaction effect in the design. For example if the highest level interaction is Etendue:levD and the Etendue and levD factors have 3 and 4 levels respectively, then the factorLevels argument should be c(3,4). 
#' @return The desired model matrix.
#' @examples
#' BFchainsModelMatrix(c(3,4))
#' @export

BFchainsModelMatrix<-function(factorLevels){
  if(class(factorLevels)!="numeric") stop("Bad factorLevels parameter.")
  if(!all(factorLevels>1)) stop("Bad factorLevels parameter.")
  dimensions<-length(factorLevels)
  blab<-rep(0,times=dimensions)
  bella<-NULL
  for(iii in 1:(dimensions)){
    blab[iii]<-1
    bella<-rbind(bella,unique(gtools::permutations(n = dimensions,r = dimensions,v = blab,set=F,repeats.allowed = F)))
  }
  gabrin<-matrix(rep(1,times=prod(factorLevels)))
  
  for(manol in 1:dim(bella)[1]){
    drab<-bella[manol,]
    labex<-vector(mode = "list",length = dim(bella)[2])
    for(jj in 1:dim(bella)[2]){
      if(drab[jj]==0){
        labex[[jj]]<-rep(1,times=factorLevels[jj])
      }else{
        labex[[jj]]<-diag(x = 1,nrow = factorLevels[jj])
      }
    }
    gabby<-diag(1,nrow = 1)
    for(jj in 1:dim(bella)[2]){
      gabby<- kronecker(gabby,labex[[jj]])
    }
    gabrin<-cbind(gabrin,gabby)
  }
  return(gabrin)
}


#' Constrained BF
#' 
#' Returns a Bayes Factor comparing a constrained model with an unconstrained one. Useful for paired comparisons.
#' 
#' @param postDT Posterior samples. Object of BFmcmc class.
#' @param const A character vector describing the constraints. For example, if you want to test the hypothesis according to which level A > level B and level B > level C, then const = "A>B>C". For more complicated constraints, use the semi-colon. For example if you want to test the hypothesis that level A > level B and that level A > level C, then const = "A>B;A>C". . If you want to test the two former hypotheses at once, const = c("A>B>C","A>B;A>C").
#' @param newColumnContrasts If supplied, this list of matrices is used to create new column from linear contrasts of the columns of the posterior samples. Each element of this list must be a matrix. The names of the columns of these matrices must correspond to column names in postDT. Rows must also be named. The row names of these matrices must be distinct from each other and from the column names of postDT. 
#' @param returnData Boolean. Indicates whether the data table of the posterior samples (possibly with the new columns) should be included in the output or not.
#' @section References:
#' Morey, R. D., & Wagenmakers, E. J. (2014). Simple relation between Bayesian order-restricted and point-null hypothesis tests. Statistics and Probability Letters, 92, 121–124. http://doi.org/10.1016/j.spl.2014.05.010
#' Klugkist, I., Laudy, O., & Hoijtink, H. (2005). Inequality constrained analysis of variance: A Bayesian approach. Psychological Methods, 10(4), 477–493. http://doi.org/10.1037/1082-989X.10.4.477
#' @return A list of lists with the following items:
#' \itemize{
#' \item prior Prior probability
#' \item post Posterior probability
#' \item BF_c_e Bayes factor comparing the constrained model vs the unconstrained model. 
#' }
#' @examples
#' data("FPpostDT")
#' out<-constrBF(FPpostDT,c("Etendue-150>Etendue-900;Etendue-300>Etendue-900","levD-1<levD-2<levD-3<levD-4"))
#' # A more complex example.
#' 
#' ooo<-data.table(as.data.frame(FPpostDT))
#' names1<-c(paste0("levD-",1:4),grep(x=names(ooo),pattern = "^Etendue:levD",value = T,fixed = F))
#' mat1<-BFchainsModelMatrix(c(3,4))[,-c(1:4)]
#' rownames(mat1)<-paste0("c_",rep(paste0("l",1:4),times=3),"_",rep(paste0("e",c(150,300,900)),each=4))
#' colnames(mat1)<-names1
#' 
#' names2<-grep(x=names(ooo),pattern = "^Etendue:levD",value = T,fixed = F)
#' mat2<-kronecker(diag(1,nrow = 3),matrix(c(1,-1, 0, 0,
#'                                           0, 1,-1, 0,
#'                                           0, 0, 1, -1),nrow = 3,byrow = T))
#' colnames(mat2)<-names2
#' temp<-data.table(expand.grid(c("1-2","2-3","3-4"),c("e150","e300","e900")))
#' temp[,tempi:=paste0("c",Var1,"_",Var2)]
#' rownames(mat2)<-temp$tempi
#' matlist<-list(mat1,mat2)
#' 
#' out2<-BFextraVL::constrBF(postDT = FPpostDT,const = c("c1-2_e150<c1-2_e300<c1-2_e900;c2-3_e150<c2-3_e300<c2-3_e900","c_l1_e150<c_l2_e150<c_l3_e150<c_l4_e150;c_l1_e300<c_l2_e300<c_l3_e300<c_l4_e300"),newColumnContrasts = matlist)
#' @export

constrBF<-function(postDT,const,newColumnContrasts=NULL,returnData=T){
  ooo<-data.table(as.data.frame(postDT))
  dabi<-newColumnContrasts
  rnames<-c()
  if(!is.null(dabi)){
    if(class(x = dabi)!="list"){
      stop("Bad newColumnContrasts parameter.")
    }
    for(iii in seq_len(length(dabi))){
      if(class(dabi[[iii]])!="matrix"){
        stop("Bad newColumnContrasts parameter.")
      }
      if(any(is.null(colnames(dabi[[iii]])))){
        stop("Bad newColumnContrasts parameter.")
      }
      if(!all(colnames(dabi[[iii]]) %in% names(ooo))){
        stop("Bad newColumnContrasts parameter.")
      }
      if(length(colnames(dabi[[iii]]))!=length(unique(colnames(dabi[[iii]])))){
        stop("Bad newColumnContrasts parameter.")
      }
      if(any(is.null(rownames(dabi[[iii]])))){
        stop("Bad newColumnContrasts parameter.")
      }
      if(any(rownames(dabi[[iii]] %in% names(ooo)))){
        stop("Bad newColumnContrasts parameter.")
      }
      rnames<-c(rnames,rownames(dabi[[iii]]))
    }
    if(length(unique(rnames))!=length(rnames)){
      stop("Bad newColumnContrasts parameter.")
    }
  }
  
  debe<-vector(mode = "list",length = length(dabi))
  for(iii in seq_along(dabi)){
    debe[[iii]]<-as.matrix(ooo[,colnames(dabi[[iii]]),with=F]) %*% t(dabi[[iii]])
  }
  for(jjj in seq_along(debe)){
    ooo<-cbind(ooo,data.table(debe[[jjj]]))
  }
  if(!(class(postDT)=="BFmcmc")) stop("Bad postDT parameter.")
  splitted<-const
  tempp<-list()
  temppp<-list()
  tempppp<-c()
  returnL<-list()
  for(i in 1:length(splitted)){
    temp<-unlist(strsplit(splitted[i],split = ";"))
    tempp[[i]]<-list()
    temppp[[i]]<-list()
    for(j in 1:length(temp)){
      numSm<-length(grep(x=temp[j],fixed = T,pattern = "<",value = T))
      numGr<-length(grep(x=temp[j],fixed = T,pattern = ">",value = T))
      if(!(xor(x = numSm>0,y = numGr>0))) stop("Bad constraint param")
      if(numSm>0){
        toko<-unlist(strsplit(temp[j],split = "<"))
        pp<-length(toko)
        for(kk in 1:pp){
          toko[kk]<-gsub(pattern = " ",replacement = "",x = toko[kk])
        }
        if(!(length(unique(toko))==length(toko))) stop("Error!")
        tempp[[i]][[j]]<-toko
        temppp[[i]][[j]]<-paste0("(`",toko[-pp],"`<`",toko[-1],"`)",collapse = "&")
      }else{
        toko<-unlist(strsplit(temp[j],split = ">"))
        pp<-length(toko)
        for(kk in 1:pp){
          toko[kk]<-gsub(pattern = " ",replacement = "",x = toko[kk])
        }
        if(!(length(unique(toko))==length(toko))) stop("Error!")
        tempp[[i]][[j]]<-toko
        temppp[[i]][[j]]<-paste0("(`",toko[-pp],"`>`",toko[-1],"`)",collapse = "&")
      }
    }
    tempppp[i]<-paste0(unlist(temppp[[i]]),collapse = "&")
    oo<-unique(unlist(tempp[[i]]))
    perms<-permutations(n = length(oo),r = length(oo),v = oo)
    matVal<-matrix(data = rep(0,times = dim(perms)[1]*dim(perms)[2]),nrow = dim(perms)[1])
    for(ww in 1:(dim(perms)[1])){
      matVal[ww,]<-chmatch(x =perms[ww,],table = oo)*-1
    }
    DTval<-data.table(matVal)
    setnames(x = DTval,old = oo)
    priorp<-mean(DTval[,eval(parse(text = tempppp[i]))])
    postp<-mean(ooo[,eval(parse(text = tempppp[i]))])
    returnL[[i]]<-list(prior=priorp,post=postp,BF_c_e=postp/priorp)
  }
  names(returnL)<-splitted
  if(returnData) returnL$Data<-ooo
  return(returnL)
}


#' Example data for BFextraVL
#' 
#' A data set used for examples.
#' 
#'
#' @docType data
#' @name FPpostDT
#' @usage data("FPpostDT")
#' @format BFmcmc
#' 

NULL

#' @import data.table
#' @import BayesFactor
#' @import gtools
#' 

NULL