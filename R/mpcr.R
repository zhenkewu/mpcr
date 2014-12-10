
      ###----------------------------------------------------------------------###
      ### BUGS fitting functions: set prior on second-level variance
      ###----------------------------------------------------------------------###

#       if (BAYES_STATUS){
#         if (!require("R2OpenBUGS")){
#           install.packages('R2OpenBUGS',type='source')
#         }
#       }
#

      call.bugs <- function(data, inits, parameters, m.file,
                            niters=20000, nburn=12000, nthin=2,
                            dic=FALSE, is.debug=FALSE, workd=getwd(),
                            bugsmodel.dir="" ) {

        m.file <- paste(bugsmodel.dir, m.file, sep="");
        f.tmp <- function() {
          if ("windows" == .Platform$OS.type | FALSE) {
            ##winbugs
            rst.bugs <- bugs(data, inits, parameters, model.file = m.file,
                             n.iter = niters, n.burnin=nburn, n.thin=nthin,
                             n.chains=1, working.directory=workd,
                             DIC=dic, useWINE=FALSE, clearWD=FALSE, debug=is.debug);
          } else {
            ##openbugs
            rst.bugs <- bugs(data, inits, parameters, floor(niters/nthin),
                             model.file = m.file,
                             n.burnin=floor(nburn/nthin),
                             n.thin=nthin,
                             n.chains=1, working.directory=workd,
                             DIC=dic, useWINE=FALSE, clearWD=FALSE, debug=is.debug);
          }
          rst.bugs;
        }

        bugs.try  <- try(rst.bugs <- f.tmp(), silent=FALSE);
        if (class(bugs.try) == "try-error") {
          rst.bugs <- NULL;
        }
        rst.bugs
      }

      bugs.1 <- function(...) {
        inits      <- list(list(mu=0,
                                sigma2delta=0.5));
        data       <- list("YDATA", "WEIGHT", "NTOTAL","HARMN");
        parameters <- c('mu', 'sigma2delta');
        rst.bugs   <- call.bugs(data, inits, parameters, "gc_bugs_2.txt",...);
        rst.bugs
      }


      ##-----------------------------------------------------------------------------##
      ## main function
      ##-----------------------------------------------------------------------------##
      #' Estimation of treatment effects in matched-pair cluster randomized
      #' trials by calibrating covarite imbalance between clusters.
      #'
      #' @param datinput The data set should have columns corresponding to the
      #' primary outcome, treatment assignment, pair IDs, cluster IDs, covariates
      #' used for covariate-calibration;
      #'
      #' @param arm The name of treatment assignment indicator. For two-arm trials,
      #' this variable takes value in {0,1}: 0 for control, 1 for treatment;
      #' @param cluster The variable name for cluster IDs;
      #' @param pair The variable name for pair IDs;
      #' @param outcome The variable name for the primary outcome;
      #' @param X_nm_all The vector of covariate names that enter the
      #' covariate-calibrated analysis
      #' @param X_nm_binary The vector of binary covariate names;
      #' @param X_nm_cat The vector of categorical (>2 categories) covariate names;
      #' @param X_nm_cont The vector of continuous covariate names.
      #'
      #' @return
      #' \itemize{
      #' \item Tables:
      #' \itemize{
      #'         \item  Table 1 - cluster sample sizes;
      #'           calibrated and uncalibrated outcome comparisons;
      #'          \item Table 2 - check covariate imbalances within each pair;
      #'          \item Table 3 -
      #'              \itemize{
      #'                   \item 1st level analysis: maximum likelihood estimate (MLE);
      #'                    permutation tests;
      #'                   \item 1st and 2nd level analysis:\cr
      #'                       MLE; profile MLE;\cr
      #'                       Bayes estimate with uniform shrinkage prior [link to paper];\cr
      #'                       permutation tests.
      #'              }
      #' }
      #'
      #' \item Figures:\cr
      #'
      #'   Check second level dependence for crude analysis \cr
      #'                      \deqn{\sqrt{v^{crude}_p} vs \delta^{crude}_p,}
      #'                      and covariate-calibrated analysis
      #'                      \deqn{\sqrt{v^{calibr}_p} vs \delta^{calibr}_p.}
      #'
      #' }
      #'
      #'
      #'
      #' @references Wu, Z., Frangakis, C. E., Louis, T. A. and Scharfstein, D. O. (2014), Estimation of treatment effects in matched-pair cluster randomized trials by calibrating covariate imbalance between clusters. Biometrics. doi: 10.1111/biom.12214
      #' @import sandwich
      #'
      #' @export


      mpcr = function(datinput,
                      arm = "tx",
                      cluster = "team",
                      pair    = "sitenew",
                      outcome = "sf36pcs32",
                      X_nm_all = c("race","ageatint","hcc","livesalone","education","s10","sf36mcs","sf36pcs","gender","h1"),
                      X_nm_binary = c("livesalone","education","gender"),
                      X_nm_cat    = c("race","s10","h1"),
                      X_nm_cont   = c("ageatint","hcc","sf36mcs","sf36pcs"),
                      BAYES_STATUS = FALSE,
                      digit.round = 2){


        cat("=====primary outcome:",outcome,"=====","\n")
        ###--------------------------------------------------------------------------###
        ###                    Covariate-Calibration
        ###--------------------------------------------------------------------------###

        calibration = function(dat,arm,pair,cluster,outcome,cov_bin,cov_cat,cov_cont,
                               siteind = 1:length(UniqSites),max_nsite = length(UniqSites),
                               use.covariate=TRUE){
          #require(sandwich)# for robust standard error calculation
          nsites = length(siteind)
          n = matrix(0,length(siteind),2)
          colnames(n) = c("cont","trt")
          ## important for jackknife, so as not confuse pair index:
          k=1
          tmp = dat[,pair]
          for (i in 1:length(siteind)) {

            n[i,1] = sum(dat[,pair]==siteind[i] & dat[,arm]==0)
            n[i,2] = sum(dat[,pair]==siteind[i] & dat[,arm]==1)
            if (n[i,1]<2 | n[i,2]<2) {
              dat$inc[dat$site==siteind[i]]=0
            }
            else {
              dat[tmp==siteind[i],pair] = k
              k = k+1
            }
          }

          #cat("Site indices are: ", siteind,"\n")
          #cat("Cluster sizes are:","\n")
          #print(n)

          ## the data enters regression should be only those included sites
          ## (useful for Jackknife, where we iteratively leave one pair out,
          ## siteind = 1:6, or 2:7, etc.)
          dat = dat[tmp %in% siteind,]
          dat = dat[order(dat[,cluster]),]

          form_bin_cont = paste(c(cov_bin,cov_cont),collapse = "+")
          form_cat = paste(paste("as.factor(",c(cluster,cov_cat),")",sep=""),collapse="+")
          model.form = as.formula(paste0(outcome,"~-1+",form_cat,"+",form_bin_cont))
          cov.form = as.formula(paste0("~-1+",form_cat,"+",form_bin_cont))

          # total 28 coefficients: 14 cluster-specific intercepts+14 covariates
          pacic_res = glm(model.form,data=dat,subset=(is.na(outcome)==FALSE))
          X = model.matrix(cov.form,data=dat)

          ##############################################################
          ###obtain robust variance-covariance matrix for coefficients
          ##############################################################
          Vrobust = sandwich(pacic_res) #Heteroskedasticity-consistent

          ###########################################################
          ## check continuous covariates
          ###########################################################
          balance_chk_cont = function(dat,nm){
            X_nm = nm
            cat("continous covariate:",X_nm,"\n")
            datX = dat[,X_nm]
            std_bias = matrix(NA,nrow=length(X_nm),ncol=nsites)
            rownames(std_bias)=X_nm
            colnames(std_bias)=1:nsites
            t_stat = std_bias

            for (j in 1:length(X_nm)){
              for (i in 1:nsites){
                x1=datX[dat[,cluster]==i,X_nm[j]]
                x2=datX[dat[,cluster]==i+nsites,X_nm[j]]
                out.t = t.test(x1,x2)
                t_stat[j,i]=out.t$statistic
                #stddev assume common variances
                stddev  = sqrt(((n[i,1]-1)*var(x1)+(n[i,2]-1)*var(x2))/(n[i,1]+n[i,2]-2))
                std_bias[j,i] = (mean(x1)-mean(x2))/stddev
              }
            }
            return(list(effectsize=std_bias,t = t_stat))
          }
          ###########################################################
          ## check categorical covariates
          ###########################################################
          balance_chk_cat = function(dat,nm){
            X_nm = nm
            cat("categorical covariate:",X_nm,"\n")
            datX = as.factor(dat[,X_nm])
            std_bias = matrix(NA,nrow=nlevels(datX),ncol=nsites)
            rownames(std_bias)=paste(X_nm,attributes(datX)$levels,sep="_")
            colnames(std_bias)=1:nsites
            pval = rep(NA,nsites)
            correction = 0.5
            for(i in 1:nsites){
              x1=datX[dat[,cluster]==i]
              x2=datX[dat[,cluster]==i+nsites]
              p1=(table(x1)+correction)/sum(table(x1)+correction)
              p2=(table(x2)+correction)/sum(table(x2)+correction)
              M <- as.table(rbind(table(x1)+correction, table(x2)+correction))
              dimnames(M) <- list(cluster = c("cont","trt"),
                                  levels = paste(X_nm,1:nlevels(datX),sep="_"))
              outchi = chisq.test(M,simulate.p.value=TRUE,B=10000)
              for (j in 1:nlevels(datX)){
                std_bias[j,i] = exp(log(p1[j]/(1-p1[j]))-log(p2[j]/(1-p2[j])))
              }
              pval[i]=outchi$p.value
            }
            return(rbind(std_bias,pval))
          }



          ##---------------------------------end of covariate balance check------####

          if (length(siteind)==max_nsite & use.covariate==FALSE){
            res_cov_cont= balance_chk_cont(dat,cov_cont)
            res_cov_cat = list()
            for (i in 1:length(c(cov_bin,cov_cat))){
              nm_cat =  c(cov_bin,cov_cat)[i]
              res_cov_cat[[i]]=balance_chk_cat(dat,nm_cat)
            }
            names(res_cov_cat) = c(c(cov_bin,cov_cat))

            naive_mn = matrix(0,nsites,2)
            naive_var = matrix(0,nsites,2)
            naive_pool = rep(0,nsites)
            for (i in 1:nsites){
              naive_mn[i,1]=mean(dat[,outcome][dat[,cluster]==i])
              naive_var[i,1]=var(dat[,outcome][dat[,cluster]==i])/length(dat[,outcome][dat[,cluster]==i])
              naive_mn[i,2]=mean(dat[,outcome][dat[,cluster]==i+nsites])
              naive_var[i,2]=var(dat[,outcome][dat[,cluster]==i+nsites])/length(dat[,outcome][dat[,cluster]==i+nsites])
              naive_pool[i] =((n[i,1]-1)*var(dat[,outcome][dat[,cluster]==i])+(n[i,2]-1)*var(dat[,outcome][dat[,cluster]==i+nsites]))/(n[i,1]+n[i,2]-2)
            }
            naive_diff = naive_mn[,1]-naive_mn[,2]
            naive_se2  = naive_var[,1]+naive_var[,2]
            effectsize = naive_diff/sqrt(naive_pool)
            jn_naive=rep(0,7)
            for (i in 1:7){
              jn_naive[i]=mean(naive_diff[-i])
            }
            jn_naive_se=sqrt(6/7*sum((jn_naive-mean(jn_naive))^2)) # jackknife s.e.

            return(list(YDATA=naive_diff,WEIGHT=naive_se2,HARMN=nsites/sum(1/naive_se2),
                        NTOTAL=length(naive_diff),effectsize = effectsize,
                        naive_pool_sqrt = sqrt(naive_pool),
                        jn_naive_se=jn_naive_se,
                        meanmat=naive_mn,n=n,
                        res_cov_cont=res_cov_cont,
                        res_cov_cat = res_cov_cat))

          } else{
            mn = matrix(0,nsites,2)
            diff = rep(0,nsites)

            dmu = matrix(0,2*nsites,dim(X)[2])
            D = matrix(0,2*nsites,2*nsites)
            S = matrix(0,nsites,2*nsites)
            for (i in 1:nsites){
              Xi = X[dat$sitenew==i,(2*nsites+1):dim(X)[2]]
              ni = dim(Xi)[1]
              u1 = matrix(0,ni,2*nsites)
              u1[,i] = 1
              X1 = cbind(u1,Xi)
              pred1 = X1 %*% pacic_res$coef
              mn[i,1] = mean(pred1)
              dmu[i,] = apply(X1,2,sum)/ni
              u2 = matrix(0,ni,2*nsites)
              u2[,i+nsites] = 1
              X2 = cbind(u2,Xi)
              pred2 = X2 %*% pacic_res$coef
              mn[i,2] = mean(pred2)
              dmu[i+nsites,] = apply(X2,2,sum)/ni
              diff[i] = mn[i,1]-mn[i,2]
              S[i,] = c(rep(0,i-1),1,rep(0,nsites-i),rep(0,i-1),-1,rep(0,nsites-i))
            }

            ########### adjusting for cluster/robust standard errors ##########
            #   require(clerror)
            #   tt2=clustvc(pacic_res, cluster = dat[,cluster])
            #   varmu = dmu %*% tt2 %*% t(dmu)
            ##################################################################

            #varmu = dmu %*% vcovHC(pacic_res) %*% t(dmu)
            #varmu = dmu %*%sandwich(pacic_res)%*% t(dmu)
            #varmu = dmu %*% vcov(pacic_res) %*% t(dmu) # model-based (const var/indep)
            varmu = dmu %*% Vrobust %*% t(dmu)
            sigma = S %*% varmu %*% t(S)

            YDATA = diff
            WEIGHT = diag(sigma)
            HARMN = nsites/sum(1/WEIGHT);

            #cat("diff=",YDATA,"\n")
            #cat("sqrtWEIGHT=",sqrt(WEIGHT),"\n")

            NTOTAL <- length(YDATA);

            #rst.bugs <- bugs.1();

            #tmp1 = c(tmp1,rst.bugs$sims.list$mu)
            #tmp2 = c(tmp2,rst.bugs$sims.list$sigma2delta)

            return(list(YDATA=diff,WEIGHT=WEIGHT,HARMN=HARMN,NTOTAL=NTOTAL,
                        meanmat=mn,n=n))
          }
        }


        ###--------------------------------------------------------------------------##
        ### main results
        ###-------------------------------------------------------------------------###
        estimation.no.cov = function(outnaive,BAYES=FALSE){
          if (BAYES){
            require("R2OpenBUGS");
            ###### 1.without covariate
            #     YDATA  = outnaive$YDATA
            #     WEIGHT = outnaive$WEIGHT
            #     HARMN  = outnaive$HARMN
            #     NTOTAL = outnaive$NTOTAL
            rst.bugs.naive <- bugs.1();
          }

          nsites <- length(outnaive$WEIGHT)
          #permutation using only 1st level
          naive_perm = as.matrix(expand.grid(rep(list(c(-1,1)),nsites)))
          naive_perm_seq = rep(NA,nrow(naive_perm))
          for (iter in 1:nrow(naive_perm)){
            naive_perm_seq[iter]=mean(as.numeric(naive_perm[iter,]*outnaive$YDATA))
          }
          pv_naive=round(sum(naive_perm_seq>=mean(outnaive$YDATA))/nrow(naive_perm),3)



          naive_diff=outnaive$YDATA
          naive_se2 =outnaive$WEIGHT

          ## maximum likelihood
          lkd = function(para,obsdiff,Sigma){
            require(mvtnorm)
            ##para[1] is delta effect
            ##para[2] is log(tau^2)
            res = - dmvnorm(obsdiff, rep(para[1],length(obsdiff)),
                            Sigma+diag(exp(para[2]),length(obsdiff)),log=TRUE)
            return(res)
          }


          start0 = c(mean(naive_diff),log(var(naive_diff)))
          out0 <- try(optim(par = start0, fn = lkd,control=list(maxit=1000),
                            obsdiff=naive_diff,Sigma=diag(naive_se2),
                            method="BFGS",hessian=TRUE))
          delta0hat = out0$par[1]
          tau0hat  = exp(out0$par[2])
          thetase0 = sqrt(solve(out0$hessian)[1,1])
          CI0 = out0$par[1]+c(-1.96,1.96)*thetase0

          ## profile likelihood

          #------------- fix theta
          maxlkd_fixtheta = function(theta,obsdiff,Sigma){
            require(mvtnorm)
            lkd_fixtheta = function(logtau2){
              require(mvtnorm)
              res = - dmvnorm(obsdiff, rep(theta,length(obsdiff)),
                              Sigma+diag(exp(logtau2),length(obsdiff)),log=TRUE)
              return(res)
            }
            start = log(var(obsdiff))
            out <- try(optim(par = start, fn = lkd_fixtheta,control=list(maxit=1000),
                             method="BFGS",hessian=TRUE))
            return(list(maxlkd = out$value,tau2=exp(out$par)))
          }

          #######naive
          theta_seq = seq(-3,3,by=0.01)
          res = as.matrix(sapply(theta_seq,maxlkd_fixtheta,obsdiff=naive_diff,Sigma=diag(naive_se2)))

          lkd_seq =-unlist(res[1,])
          theta_range=theta_seq[which(lkd_seq>=max(lkd_seq)-1.92)]
          profileCI=c(min(theta_range),max(theta_range))

          ## permutation (different weights)
          naive_mn = outnaive$meanmat

          perm_ind=as.matrix(expand.grid(rep(list(c(1,2)),nsites)))
          perm_diff = list()
          for (permiter in 1:nrow(perm_ind)){
            perm_diff[[permiter]] = naive_mn[cbind(1:nrow(naive_mn),perm_ind[permiter,])]-
              naive_mn[cbind(1:nrow(naive_mn),3-perm_ind[permiter,])]
          }

          for(iter in 1){
            start = c(mean(perm_diff[[iter]]),log(var(perm_diff[[iter]])))
            out1 <- try(optim(par = start, fn = lkd,control=list(maxit=1000),
                              obsdiff=perm_diff[[iter]],Sigma=diag(naive_se2),
                              method="BFGS",hessian=TRUE))
          }

          perm_res = rep(NA,nrow(perm_ind))
          for(iter in 1:nrow(perm_ind)){
            start = c(mean(perm_diff[[iter]]),log(var(perm_diff[[iter]])))
            outtemp <- try(optim(par = start, fn = lkd,control=list(maxit=1000),
                                 obsdiff=perm_diff[[iter]],Sigma=diag(naive_se2),
                                 method="BFGS",hessian=TRUE))
            perm_res[iter]=outtemp$par[1]
          }

          pv_naive_2nd=round(sum(perm_res>=out1$par[1])/nrow(perm_ind),2)

          ###### collect results
          nc_1st_MLE = c(mean(outnaive$YDATA),
                         mean(outnaive$YDATA)-1.96*outnaive$jn_naive_se,
                         mean(outnaive$YDATA)+1.96*outnaive$jn_naive_se,
                         outnaive$jn_naive_se,
                         NA,
                         2-2*pnorm(mean(outnaive$YDATA)/outnaive$jn_naive_se))
          nc_1st_permutation =c(NA,NA,NA,NA,NA,pv_naive*2)
          nc_1st2nd_MLE = c(delta0hat,CI0[1],CI0[2],thetase0,
                            tau0hat,2-2*pnorm(delta0hat/thetase0))
          nc_1st2nd_pMLE = c(delta0hat,profileCI[1],profileCI[2],NA,tau0hat,NA)
          nc_1st2nd_permutation =c(NA,NA,NA,NA,NA,pv_naive_2nd*2)
          if (BAYES){
            nc_1st2nd_Bayes = c(rst.bugs.naive$summary[1,1],
                                rst.bugs.naive$summary[1,3],
                                rst.bugs.naive$summary[1,7],
                                rst.bugs.naive$summary[1,2],
                                rst.bugs.naive$summary[2,1],
                                2-2*pnorm(rst.bugs.naive$summary[1,1]/rst.bugs.naive$summary[1,2]))

            mat_nc = rbind(nc_1st_MLE,
                           nc_1st_permutation,
                           nc_1st2nd_MLE,
                           nc_1st2nd_pMLE,
                           nc_1st2nd_Bayes,
                           nc_1st2nd_permutation)
          }else{
            mat_nc = rbind(nc_1st_MLE,
                           nc_1st_permutation,
                           nc_1st2nd_MLE,
                           nc_1st2nd_pMLE,
                           nc_1st2nd_permutation)
          }
          colnames(mat_nc)=c("delta_hat","2.5%","97.5%","se(delta_hat)","tau^2_hat","pval")
          return(mat_nc)
        }

        estimation.cov = function(outcalibr,max_nsite=length(UniqSites),BAYES=FALSE){
          if (BAYES){

            require(R2OpenBUGS)
            ######## 2.with covariate
            #     YDATA  = outcalibr$YDATA
            #     WEIGHT = outcalibr$WEIGHT
            #     HARMN  = outcalibr$HARMN
            #     NTOTAL = outcalibr$NTOTAL

            rst.bugs.calibr <- bugs.1();
          }

          nsites <- length(outcalibr$WEIGHT)
          #permutation using only 1st level
          calibr_perm = as.matrix(expand.grid(rep(list(c(-1,1)),nsites)))
          calibr_perm_seq = rep(NA,nrow(calibr_perm))
          for (iter in 1:nrow(calibr_perm)){
            calibr_perm_seq[iter]=mean(as.numeric(calibr_perm[iter,]*outcalibr$YDATA))
          }
          pv_calibr=round(sum(calibr_perm_seq>=mean(outcalibr$YDATA))/nrow(calibr_perm),3)


          #############################
          ## 1st level by Jackknife
          #############################
          jn_mean = rep(0,length=max_nsite)
          for (deletesite in 1:max_nsite){
            remainsite = (1:max_nsite)[-deletesite]
            calibr = calibration(datinput,arm,pair,cluster,outcome,X_nm_binary,X_nm_cat,X_nm_cont,siteind=remainsite,use.covariate=TRUE)
            jn_mean[deletesite]=mean(calibr$YDATA)
          }
          #full=calibration(dat,siteind=1:7)$y
          njn=max_nsite
          jn_se = sqrt((njn-1)/njn*sum((jn_mean-mean(jn_mean))^2))
          #pseudo = njn*mean(full)-(njn-1)*jn_mean
          #theta_tilda = mean(pseudo)


          #####################################
          ## input previous calibration results
          #####################################

          diff=outcalibr$YDATA
          sigma=diag(outcalibr$WEIGHT)

          ## maximum likelihood
          lkd = function(para,obsdiff,Sigma){
            require(mvtnorm)
            ##para[1] is delta effect
            ##para[2] is log(tau^2)
            res = - dmvnorm(obsdiff, rep(para[1],length(obsdiff)),
                            Sigma+diag(exp(para[2]),length(obsdiff)),log=TRUE)
            return(res)
          }

          start = c(mean(diff),log(var(diff)))
          out <- try(optim(par = start, fn = lkd,control=list(maxit=1000),
                           obsdiff=diff,Sigma=sigma,
                           method="BFGS",hessian=TRUE))
          thetase = sqrt(solve(out$hessian)[1,1])
          CI=out$par[1]+c(-1.96,1.96)*thetase
          deltahat= out$par[1]
          tauhat = exp(out$par[2])

          ## profile likelihood

          #------------- fix theta
          maxlkd_fixtheta = function(theta,obsdiff,Sigma){
            require(mvtnorm)
            lkd_fixtheta = function(logtau2){
              require(mvtnorm)
              res = - dmvnorm(obsdiff, rep(theta,length(obsdiff)),
                              Sigma+diag(exp(logtau2),length(obsdiff)),log=TRUE)
              return(res)
            }
            start = log(var(obsdiff))
            out <- try(optim(par = start, fn = lkd_fixtheta,control=list(maxit=1000),
                             method="BFGS",hessian=TRUE))
            return(list(maxlkd = out$value,tau2=exp(out$par)))
          }

          theta_seq = seq(-0.4,3,by=0.01)
          res = as.matrix(sapply(theta_seq,maxlkd_fixtheta,obsdiff=diff,Sigma=sigma))

          lkd_seq =-unlist(res[1,])
          theta_range=theta_seq[which(lkd_seq>=max(lkd_seq)-1.92)]
          profileCI=c(min(theta_range),max(theta_range))

          ######################## permutation (change label of intervention/control
          ###### within each pair)
          ######################################################################
          mn=outcalibr$meanmat

          perm_ind=as.matrix(expand.grid(rep(list(c(1,2)),nsites)))
          perm_diff = list()
          for (permiter in 1:nrow(perm_ind)){
            perm_diff[[permiter]] = mn[cbind(1:nrow(mn),perm_ind[permiter,])]-
              mn[cbind(1:nrow(mn),3-perm_ind[permiter,])]
          }

          for(iter in 1){
            start = c(mean(perm_diff[[iter]]),log(var(perm_diff[[iter]])))
            out1 <- try(optim(par = start, fn = lkd,control=list(maxit=1000),
                              obsdiff=perm_diff[[iter]],Sigma=sigma,
                              method="BFGS",hessian=TRUE))
          }

          perm_res = rep(NA,nrow(perm_ind))
          for(iter in 1:nrow(perm_ind)){
            start = c(mean(perm_diff[[iter]]),log(var(perm_diff[[iter]])))
            outtemp <- try(optim(par = start, fn = lkd,control=list(maxit=1000),
                                 obsdiff=perm_diff[[iter]],Sigma=sigma,
                                 method="BFGS",hessian=TRUE))
            perm_res[iter]=outtemp$par[1]
          }

          pv_calibr_2nd=round(sum(perm_res>=out1$par[1])/nrow(perm_ind),3)

          ###### collect results (with covariates)
          c_1st_MLE = c(mean(outcalibr$YDATA),
                        mean(outcalibr$YDATA)-1.96*jn_se,
                        mean(outcalibr$YDATA)+1.96*jn_se,
                        jn_se,
                        NA,
                        2-2*pnorm(mean(outcalibr$YDATA)/jn_se))
          c_1st_permutation =c(NA,NA,NA,NA,NA,pv_calibr*2)
          c_1st2nd_MLE = c(deltahat,CI[1],CI[2],thetase,
                           tauhat,2-2*pnorm(deltahat/thetase))
          c_1st2nd_pMLE = c(deltahat,profileCI[1],profileCI[2],NA,tauhat,NA)
          c_1st2nd_permutation =c(NA,NA,NA,NA,NA,pv_calibr_2nd*2)
          if (BAYES){
            c_1st2nd_Bayes = c(rst.bugs.calibr$summary[1,1],
                               rst.bugs.calibr$summary[1,3],
                               rst.bugs.calibr$summary[1,7],
                               rst.bugs.calibr$summary[1,2],
                               rst.bugs.calibr$summary[2,1],
                               2-2*pnorm(rst.bugs.calibr$summary[1,1]/rst.bugs.calibr$summary[1,2]))

            mat_c = rbind(c_1st_MLE,
                          c_1st_permutation,
                          c_1st2nd_MLE,
                          c_1st2nd_pMLE,
                          c_1st2nd_Bayes,
                          c_1st2nd_permutation)
          } else{
            mat_c = rbind(c_1st_MLE,
                          c_1st_permutation,
                          c_1st2nd_MLE,
                          c_1st2nd_pMLE,
                          c_1st2nd_permutation)
          }
          colnames(mat_c)=c("delta_hat","2.5%","97.5%","se(delta_hat)","tau^2_hat","pval")

          return(mat_c)
        }




        ## arrays to store results in table 1 and table 3
        table1 = array(NA,c(10,7))
        table3 = array(NA,c(12,6))

        ## specify variable names for treatment arm, cluster, pair,
        ## primary outcome, and covariates:

        covariates = c(X_nm_binary,X_nm_cat,X_nm_cont)

        UniqSites   = sort(unique(datinput[,pair]))

        cat("Number of pairs:",length(UniqSites),"\n")
        res_full_cov = calibration(datinput,arm,pair,cluster,outcome,X_nm_binary,X_nm_cat,X_nm_cont,siteind=1:length(UniqSites),use.covariate=TRUE)
        res_full_no_cov=calibration(datinput,arm,pair,cluster,outcome,X_nm_binary,X_nm_cat,X_nm_cont,siteind=1:length(UniqSites),use.covariate=FALSE)


        ###------------------------------------------------###
        ### making Table 1: summary of primary outcome for uncalibrated
        ### vs calibrated approaches
        ### Table 1: Mean potential outcomes, d_p, v_p
        ###------------------------------------------------###

        #dir.create("tables")

        outnaive = res_full_no_cov
        outcalibr = res_full_cov
        table1df1 = round(data.frame(mu1=outnaive$meanmat[,1],mu2=outnaive$meanmat[,2],
                                     dp=outnaive$YDATA,#effectsize = outnaive$effectsize,
                                     sqrtvp=sqrt(outnaive$WEIGHT)),2)
        table1df2 = round(data.frame(mu1=outcalibr$meanmat[,1],mu2=outcalibr$meanmat[,2],
                                     dp=outcalibr$YDATA,#effectsize = outcalibr$YDATA/naive_pool_sqrt,
                                     sqrtwp=sqrt(outcalibr$WEIGHT)),3)

        TABLE1 = round(rbind(t(outcalibr$n),t(table1df1),t(table1df2)),digit.round)
        table1_nm = paste0("tables//table1",".csv")
        #write.csv(TABLE1,table1_nm)
        table1 = TABLE1
        ####------------------------------------------------------------------------####
        ### Table 2: covariate imbalance check
        ####------------------------------------------------------------------------####

        TABLE2_cont = outnaive$res_cov_cont
        TABLE2_cat  = outnaive$res_cov_cat

        ###-------------------------------------------------------------------------###
        ### Table 3 main results
        ###-------------------------------------------------------------------------###
        YDATA  <<- outnaive$YDATA
        WEIGHT <<- outnaive$WEIGHT
        HARMN  <<- outnaive$HARMN
        NTOTAL <<- outnaive$NTOTAL
        nc_table = estimation.no.cov(outnaive,BAYES=BAYES_STATUS)

        YDATA  <<- outcalibr$YDATA
        WEIGHT <<- outcalibr$WEIGHT
        HARMN  <<- outcalibr$HARMN
        NTOTAL <<- outcalibr$NTOTAL
        c_table  = estimation.cov(outcalibr,BAYES=BAYES_STATUS)

        TABLE3 = round(rbind(nc_table,c_table),digit.round)
        table3_nm = paste0("tables//table3",".csv")
        #write.csv(TABLE3,table3_nm)
        table3=TABLE3

        ###------------------------------------------------###
        ### Figure 3: 2nd level dependence_check
        ###  v^crude_p vs delta^crude_p, and v^calibr_p vs delta^calibr_p
        ###------------------------------------------------###
        #dir.create("figures")
        fig_nm=paste("figures//dependence_check","_.pdf")
        #pdf(fig_nm,width=12)
        par(mfrow=c(1,2))
        par(mar=c(6,6,4,2));
        plot(table1df1$dp,table1df1$sqrtvp,
             xlab=expression(delta[p]^crude),ylab=expression(sqrt(v[p]^crude)),
             ylim=c(0,3),xlim=c(-3,5),las=1,cex.lab=2,cex=2,xaxt="n",yaxt="n")
        axis(1,at=c(-2,0,2,4),labels=c(-2,0,2,4),cex.axis=2)
        axis(2,at=seq(0,3,by=0.5),labels=seq(0,3,by=0.5),cex.axis=2)
        fit1 = lm(table1df1$sqrtvp~table1df1$dp)
        mtext("no calibration",3,line=2,cex=2)
        Rsq1 = round(summary(fit1)$r.squared,2)
        text(0,1.25,paste0("R squared=",Rsq1),cex=2)

        plot(table1df2$dp,table1df2$sqrtwp,
             xlab=expression(delta[p]^calibr),ylab=expression(sqrt(v[p]^calibr)),
             ylim=c(0,3),xlim=c(-3,5),pch=20,cex=2,las=1,cex.lab=2,xaxt="n",yaxt="n")
        axis(1,at=c(-2,0,2,4),labels=c(-2,0,2,4),cex.axis=2)
        axis(2,at=seq(0,3,by=0.5),labels=seq(0,3,by=0.5),cex.axis=2)
        fit2 = lm(table1df2$sqrtwp~table1df2$dp)
        Rsq2 = round(summary(fit2)$r.squared,2)
        text(0,1.25,paste0("R squared=",Rsq2),cex=2)
        mtext("with calibration",3,line=2,cex=2)
        #dev.off()

        res = list(table1,TABLE2_cont,TABLE2_cat,table3)
        names(res) = c("primary outcome summaries",
                       "covariate balance check (continuous)",
                       "covariate balance check (categorical)",
                       "estimates")
        return(res)
      }












