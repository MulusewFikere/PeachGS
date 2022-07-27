#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# Complete script for linear mixed model analysis of peach GxE paper
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Preliminary
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.2. Load libraries
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    library(gtools)
    library(asreml)  
    
    options(scipen = 999)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.3. Source functions
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 1.3.2. Estimate vcov from fa model
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        fatocov.f <- function(fatocov_in)
        {
            fc.Mpsi <- diag(fatocov_in$gammas[grepl(fatocov_in$fac,names(fatocov_in$gammas),fixed=T) &
                                                  grepl(fatocov_in$E,names(fatocov_in$gammas),fixed=T) &
                                                  grepl("var",names(fatocov_in$gammas))])
            rownames(fc.Mpsi) <- colnames(fc.Mpsi) <- sub(paste(fatocov_in$E,".",sep=""),"",
                                                          sub(".var","",
                                                              sapply(
                                                                  strsplit(
                                                                      names(fatocov_in$gammas[
                                                                          grepl(fatocov_in$fac,names(fatocov_in$gammas),fixed=T) &
                                                                              grepl(fatocov_in$E,names(fatocov_in$gammas),fixed=T) &
                                                                              grepl("var",names(fatocov_in$gammas))]),
                                                                      split = "!",fixed = T),
                                                                  "[[",2)
                                                          )
            )
            
            fc.Mlam <- matrix(fatocov_in$gammas[grepl(fatocov_in$fac,names(fatocov_in$gammas),fixed=T) &
                                                    grepl(fatocov_in$E,names(fatocov_in$gammas),fixed=T) &
                                                    !grepl("var",names(fatocov_in$gammas))],
                              nrow = nrow(fc.Mpsi))
            
            rownames(fc.Mlam) <- rownames(fc.Mpsi)
            colnames(fc.Mlam) <- paste("fa",1:ncol(fc.Mlam),sep="")
            
            fc.Mvcov <- as.matrix(fc.Mlam%*%t(fc.Mlam) + fc.Mpsi)
            
            fc_out <- NULL
            fc_out$Mpsi <- as.matrix(fc.Mpsi)
            fc_out$Mlam <- as.matrix(fc.Mlam)
            fc_out$Mvcov <- as.matrix(fc.Mvcov)
            
            fc_out$Mr <- cov2cor(fc_out$Mvcov)
            
            return(fc_out)
            
        }
        
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 1.3.3. Test model differences
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        tmd.f <- function(tmd_in)
        {
            #<<<<<<<<<<<<<<<<<<<
            # 1.3.3.1. Input tmd format
            #
            # tmd_in <- NULL
            # tmd_in$analy_name
            # tmd_in$fullmod_name
            # tmd_in$redmod_name
            # tmd_in$fullmod.asr
            # tmd_in$redmod.asr
            # tmd_in$test ("boundary"/"other")
            #
            #<<<<<<<<<<<<<<<<<<<
            
            analysum_temp <- NULL
            analysum_temp$analy_name <- tmd_in$analy_name
            analysum_temp$full_name <- tmd_in$fullmod_name
            analysum_temp$red_name <- tmd_in$redmod_name
            
            analysum_temp$fullmod.logl <- tmd_in$fullmod.asr$loglik
            analysum_temp$redmod.logl <- tmd_in$redmod.asr$loglik
            
            analysum_temp$fullmod.df <- sum(summary(tmd_in$fullmod.asr)$var$bound%in%c("P","U"))
            analysum_temp$redmod.df <- sum(summary(tmd_in$redmod.asr)$var$bound%in%c("P","U"))
            
            analysum_temp$Dlogl <- round(2*(analysum_temp$fullmod.logl - analysum_temp$redmod.logl),3)
            analysum_temp$Ddf <- analysum_temp$fullmod.df - analysum_temp$redmod.df
            
            analysum_temp$test <- tmd_in$test
            
            if(tmd_in$test == 'boundary')
            {
                if(analysum_temp$Ddf <1)
                {
                    analysum_temp$pDlogl = (1 - pchisq(analysum_temp$Dlogl,1))/2
                } else
                {
                    analysum_temp$pDlogl <- (1 - pchisq(analysum_temp$Dlogl,analysum_temp$Ddf))/2
                }
            } else
            {
                if(analysum_temp$Ddf <1)
                {
                    analysum_temp$pDlogl = (1 - pchisq(analysum_temp$Dlogl,1))
                } else
                {
                    analysum_temp$pDlogl <- (1 - pchisq(analysum_temp$Dlogl,analysum_temp$Ddf))
                }
                
            }
            
            return(as.data.frame(analysum_temp))
        }
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
# 2. prepare data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.1. Load curated data
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    load("../../4. Fitting linear models/4.3. XGWX/XGWX_210918.RData")
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Buidl GRMS
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
AQ.giv <- AQ.giv_QID_tune
DQ.giv <- DQ.giv_QID_tune
AB.giv <- AB.giv_tune
DB.giv <- DB.giv_tune
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 5. Multi-location QRM + BRM
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.1. Prep apdata
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    unique(apdata$AE4)
    
    apdata$fAQ1 <- factor(apdata$AE4)
    apdata$fAB1 <- factor(apdata$AE4)
    
    apdata$AQID <- factor(apdata$QID,rownames(AQ.giv))
    apdata$ABID <- factor(apdata$GID,rownames(AB.giv))

    apdata$DQID <- factor(apdata$QID,rownames(DQ.giv))
    apdata$DBID <- factor(apdata$GID,rownames(DB.giv))
    
    apdata$fU <- factor(apdata$U)  
    
    unique(apdata$fAQ1)
    unique(apdata$fAB1)
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.2. Univariate models
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.1. XQBU01
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.1.1. Fit model
            #<<<<<<<<<<<<<<<<<<<

            asreml.options(dense = ~ vm(AQID, AQ.giv))
            asreml.options(dense = ~ vm(ABID, AB.giv))
            asreml.options(dense = ~ vm(DQID, DQ.giv))
            asreml.options(dense = ~ vm(DBID, DB.giv))
            
            XQBU01.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) + vm(AQID, AQ.giv):idv(fAQ1) +
                                                vm(ABID, AB.giv) + vm(ABID, AB.giv):idv(fAB1) +
                                                vm(DQID, DQ.giv) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU),
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
                
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.1.2. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU01.asr$loglik
            
            summary(XQBU01.asr)$var
            
            sum(summary(XQBU01.asr)$var[,'bound']%in%c('P','U'))
    
            XQBU01_tmd_temp <- NULL
            XQBU01_tmd_temp$analy_name <- 'XQBU01' 
            XQBU01_tmd_temp$full_name <- 'XQBU01' 
            XQBU01_tmd_temp$fullmod.logl <- XQBU01.asr$loglik
            XQBU01_tmd_temp$fullmod.df <- sum(summary(XQBU01.asr)$var$bound%in%c("P","U"))
            XQBU01_tmd_temp$ranmod <- as.character(XQBU01.asr$call)[3]
            XQBU01_tmd_temp$AQE <- toString(unique(apdata$fAQ1))
            XQBU01_tmd_temp$ABE <- toString(unique(apdata$fAB1))
            
            XQBU01_tmd_temp$Bestmod <- 'XQBU01'
            XQBU01_tmd_temp
            
            tmd_sum <- as.data.frame(XQBU01_tmd_temp)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.2. XQBU02 remove  giv(AQID):idv(fAQ1)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.1. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU02.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) + vm(AQID, AQ.giv):idv(fAQ1) +
                                                vm(ABID, AB.giv) + vm(ABID, AB.giv):idv(fAB1) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU),
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.2. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBU02.asr)
            
            summary(XQBU02.asr)$var
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.3. Test difference bn models
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU02_tmd_in <- NULL
            XQBU02_tmd_in$analy_name <- 'XQBU02 = XWBU01' 
            XQBU02_tmd_in$fullmod_name <- 'XQBU01' 
            XQBU02_tmd_in$redmod_name <- 'XQBU02'
            XQBU02_tmd_in$fullmod.asr <- XQBU01.asr
            XQBU02_tmd_in$redmod.asr <- XQBU02.asr
            XQBU02_tmd_in$test <- "boundary" #("boundary"/"conventional")
            
            XQBU02_tmd_out <- tmd.f(tmd_in = XQBU02_tmd_in)
            XQBU02_tmd_out
            
#CH: XQBU02 < XQBU01_pd               
#CH: Best Model = XQBU01 
            
            XQBU02_tmd_out$Bestmod <- 'XQBU02'
            XQBU02_tmd_out$ranmod <- as.character(XQBU02.asr$call)[3]
            XQBU02_tmd_out$AQE <- toString(unique(apdata$fAQ1))
            XQBU02_tmd_out$ABE <- toString(unique(apdata$fAB1))
            
            tmd_sum <- smartbind(tmd_sum,XQBU02_tmd_out)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.2. XQBU03 remove  DB
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.1. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU03.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) + vm(AQID, AQ.giv):idv(fAQ1) +
                                                vm(ABID, AB.giv) + vm(ABID, AB.giv):idv(fAB1) +
                                                at(L):idv(fU),
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.2. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBU03.asr)
            
            summary(XQBU03.asr)$var
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.3. Test difference bn models
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU03_tmd_in <- NULL
            XQBU03_tmd_in$analy_name <- 'XQBU03=XQBU02' 
            XQBU03_tmd_in$fullmod_name <- 'XQBU02' 
            XQBU03_tmd_in$redmod_name <- 'XQBU03'
            XQBU03_tmd_in$fullmod.asr <- XQBU02.asr
            XQBU03_tmd_in$redmod.asr <- XQBU03.asr
            XQBU03_tmd_in$test <- "boundary" #("boundary"/"conventional")
            
            XQBU03_tmd_out <- tmd.f(tmd_in = XQBU03_tmd_in)
            XQBU03_tmd_out
            
            XQBU03_tmd_out$Bestmod <- 'XQBU02'
            XQBU03_tmd_out$ranmod <- as.character(XQBU03.asr$call)[3]
            XQBU03_tmd_out$AQE <- toString(unique(apdata$fAQ1))
            XQBU03_tmd_out$ABE <- toString(unique(apdata$fAB1))
            
            tmd_sum <- smartbind(tmd_sum,XQBU03_tmd_out)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.2. XQBU04 remove  giv(AQID):idv(fAQ1)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.1. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU04.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) + 
                                                vm(ABID, AB.giv) + vm(ABID, AB.giv):idv(fAB1) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU),
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.2. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBU04.asr)
            
            summary(XQBU04.asr)$var
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.3. Test difference bn models
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU04_tmd_in <- NULL
            XQBU04_tmd_in$analy_name <- 'XQBU04=XQBU02' 
            XQBU04_tmd_in$fullmod_name <- 'XQBU02' 
            XQBU04_tmd_in$redmod_name <- 'XQBU04'
            XQBU04_tmd_in$fullmod.asr <- XQBU02.asr
            XQBU04_tmd_in$redmod.asr <- XQBU04.asr
            XQBU04_tmd_in$test <- "boundary" #("boundary"/"conventional")
            
            XQBU04_tmd_out <- tmd.f(tmd_in = XQBU04_tmd_in)
            XQBU04_tmd_out
            
            #CH: XQBU04 < XQBU01_pd               
            #CH: Best Model = XQBU01 
            
            XQBU04_tmd_out$Bestmod <- 'XQBU02'
            XQBU04_tmd_out$ranmod <- as.character(XQBU04.asr$call)[3]
            XQBU04_tmd_out$AQE <- toString(unique(apdata$fAQ1))
            XQBU04_tmd_out$ABE <- toString(unique(apdata$fAB1))
            
            tmd_sum <- smartbind(tmd_sum,XQBU04_tmd_out)
    
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.2. XQBU05 remove  giv(AQID):idv(fAQ1)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.1. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU05.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) + vm(AQID, AQ.giv):idv(fAQ1) +
                                                vm(ABID, AB.giv) + 
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU),
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.2. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBU05.asr)
            
            summary(XQBU05.asr)$var
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.2.2.3. Test difference bn models
            #<<<<<<<<<<<<<<<<<<<
            
            XQBU05_tmd_in <- NULL
            XQBU05_tmd_in$analy_name <- 'XQBU05=XQBU02' 
            XQBU05_tmd_in$fullmod_name <- 'XQBU02' 
            XQBU05_tmd_in$redmod_name <- 'XQBU05'
            XQBU05_tmd_in$fullmod.asr <- XQBU02.asr
            XQBU05_tmd_in$redmod.asr <- XQBU05.asr
            XQBU05_tmd_in$test <- "boundary" #("boundary"/"conventional")
            
            XQBU05_tmd_out <- tmd.f(tmd_in = XQBU05_tmd_in)
            XQBU05_tmd_out
            
            #CH: XQBU05 < XQBU01_pd               
            #CH: Best Model = XQBU01 
            
            XQBU05_tmd_out$Bestmod <- 'XQBU02'
            XQBU05_tmd_out$ranmod <- as.character(XQBU05.asr$call)[3]
            XQBU05_tmd_out$AQE <- toString(unique(apdata$fAQ1))
            XQBU05_tmd_out$ABE <- toString(unique(apdata$fAB1))
            
            tmd_sum <- smartbind(tmd_sum,XQBU05_tmd_out)
            
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.4. Multivariate - Background (B)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.3.1. XQBMB1
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<
            # 5.3.1.2. Fit model
            #<<<<<<<<<<<<<<<<<<<
        
            asreml.options(dense = ~ vm(AQID, AQ.giv))
            asreml.options(dense = ~ vm(ABID, AB.giv))
            asreml.options(dense = ~ vm(DBID, DB.giv))

            XQBMB1.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) +
                                                vm(ABID, AB.giv):fa(fAB1,1) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU), 
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.3.1.3. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBMB1.asr)
            
            summary(XQBMB1.asr)$var
                
                ##########
                # 5.3.2.3.2. AB
                ##########
                
                XQBMB1.Mvcov_AB_F2C_in <- NULL
                XQBMB1.Mvcov_AB_F2C_in$gammas <- XQBMB1.asr$vparameters
                XQBMB1.Mvcov_AB_F2C_in$fac <- "ABID"
                XQBMB1.Mvcov_AB_F2C_in$E <- "fAB1"
                
                XQBMB1.Mvcov_AB <- fatocov.f(fatocov_in = XQBMB1.Mvcov_AB_F2C_in) 
                XQBMB1.Mvcov_AB
                
                #<<<<<<<<<
                # 5.3.2.3.3. Cluster of rG
                #<<<<<<<<<
                
                XQBMB1.vAQ <- XQBMB1.asr$vparameters[grepl("AQID",names(XQBMB1.asr$vparameters))]
                XQBMB1.vDB <- XQBMB1.asr$vparameters[grepl("DBID",names(XQBMB1.asr$vparameters))]
                
                XQBMB1.MvcovG <- XQBMB1.Mvcov_AB$Mvcov + XQBMB1.vAQ + XQBMB1.vAQ
                XQBMB1.MrG <- cov2cor(XQBMB1.MvcovG)
                
                plot(as.dendrogram(hclust(dist(XQBMB1.MrG), method = 'ward.D2'),hang = -1), horiz = T, xlab = "AB")

                            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.4.2. XQBMB2
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             
            ######2#############
            # 5.4.1.1. set up data
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XQBMB1.Mvcov_AB$Mvcov)),hang = -1)
            
            apdata$AB2 <- as.character(apdata$fAB1)
            apdata$AB2[apdata$fAB1 == 'F' | apdata$fAB1 == 'K2'] <- 'FK2'
            apdata$fAB2 <- factor(apdata$AB2) 
            
            unique(apdata$AB2)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.4.2.2. Fit model
            #<<<<<<<<<<<<<<<<<<<
             
            XQBMB2.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) +
                                                vm(ABID, AB.giv):fa(fAB2,1) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU), 
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
             #<<<<<<<<<<<<<<<<<<<
             # 5.4.2.3. Review results
             #<<<<<<<<<<<<<<<<<<<
             
             wald(XQBMB2.asr)
             
             summary(XQBMB2.asr)$var
             
                 ##########
                 # 5.4.2.2.2. AB
                 ##########
                 
                 XQBMB2.Mvcov_AB_F2C_in <- NULL
                 XQBMB2.Mvcov_AB_F2C_in$gammas <- XQBMB2.asr$vparameters
                 XQBMB2.Mvcov_AB_F2C_in$fac <- "ABID"
                 XQBMB2.Mvcov_AB_F2C_in$E <- "fAB2"
                 
                 XQBMB2.Mvcov_AB <- fatocov.f(fatocov_in = XQBMB2.Mvcov_AB_F2C_in) 
                 XQBMB2.Mvcov_AB
             
                 ##########
                 # 5.2.1.2.3. Test difference bn models
                 ##########
                 
                 XQBMB2_tmd_in <- NULL
                 XQBMB2_tmd_in$analy_name <- 'XQBMB2 = XQBMB1' 
                 XQBMB2_tmd_in$fullmod_name <- 'XQBMB1' 
                 XQBMB2_tmd_in$redmod_name <- 'XQBMB2'
                 XQBMB2_tmd_in$fullmod.asr <- XQBMB1.asr
                 XQBMB2_tmd_in$redmod.asr <- XQBMB2.asr
                 XQBMB2_tmd_in$test <- "conventional" #("boundary"/"conventional")
                 
                 XQBMB2_tmd_out <- tmd.f(tmd_in = XQBMB2_tmd_in)
                 XQBMB2_tmd_out
                 
#CH: XQBMB2 = XQBMB2MB1  
#CH: Best model = XQBMB2
                 
                 XQBMB2_tmd_out$Bestmod <- 'XQBMB2'
                 XQBMB2_tmd_out$ranmod <- as.character(XQBMB2.asr$call)[3]
                 XQBMB2_tmd_out$AQE <- toString(unique(apdata$fAQ3))
                 XQBMB2_tmd_out$ABE <- toString(unique(apdata$fAB2))
                 
                 tmd_sum <- smartbind(tmd_sum,XQBMB2_tmd_out) 
                 
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.4.3. XQBMB3
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             
            #<<<<<<<<<<<<<<<<<<<
            # 5.4.3.1. set up data
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XQBMB1.Mvcov_AB$Mvcov)),hang = -1)
            
            apdata$AB3 <- as.character(apdata$fAB1)
            unique(apdata$AB3)
            apdata$AB3[apdata$fAB1 == 'F' | apdata$fAB1 == 'G'] <- 'FG'
            apdata$fAB3 <- factor(apdata$AB3)
            
            unique(apdata$fAB3)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.4.3.2. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XQBMB3.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) +
                                                vm(ABID, AB.giv):fa(fAB3,1) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU), 
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.4.3.3. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBMB3.asr)
            
            summary(XQBMB3.asr)$var
            
                # ##########
                # # 5.4.3.3.1. AQ
                # ##########
                # 
                # XQBMB3.Mvcov_AQ_F2C_in <- NULL
                # XQBMB3.Mvcov_AQ_F2C_in$gammas <- XQBMB3.asr$vparameters
                # XQBMB3.Mvcov_AQ_F2C_in$E <- "fAQ3"
                # 
                # XQBMB3.Mvcov_AQ <- fatocov.f(fatocov_in = XQBMB3.Mvcov_AQ_F2C_in) 
                # XQBMB3.Mvcov_AQ
                # 
                ##########
                # 5.4.3.3.2. AB
                ##########

                XQBMB3.Mvcov_AB_F2C_in <- NULL
                XQBMB3.Mvcov_AB_F2C_in$gammas <- XQBMB3.asr$vparameters
                XQBMB3.Mvcov_AB_F2C_in$fac <- "ABID"
                XQBMB3.Mvcov_AB_F2C_in$E <- "fAB3"

                XQBMB3.Mvcov_AB <- fatocov.f(fatocov_in = XQBMB3.Mvcov_AB_F2C_in)
                XQBMB3.Mvcov_AB
                
                ##########
                # 5.4.3.3.3. Test difference bn models
                ##########
                
                XQBMB3_tmd_in <- NULL
                XQBMB3_tmd_in$analy_name <- 'XQBMB3' 
                XQBMB3_tmd_in$fullmod_name <- 'XQBMB1' 
                XQBMB3_tmd_in$redmod_name <- 'XQBMB3'
                XQBMB3_tmd_in$fullmod.asr <- XQBMB1.asr
                XQBMB3_tmd_in$redmod.asr <- XQBMB3.asr
                XQBMB3_tmd_in$test <- "conventional" #("boundary"/"conventional")
                
                XQBMB3_tmd_out <- tmd.f(tmd_in = XQBMB3_tmd_in)
                XQBMB3_tmd_out
                
#CH: XQBMB3 = XQBMB2  
#CH: Best model = XQBMB3
                
                XQBMB3_tmd_out$Bestmod <- 'XQBMB2'
                XQBMB3_tmd_out$ranmod <- as.character(XQBMB3.asr$call)[3]
                XQBMB3_tmd_out$AQE <- toString(unique(apdata$fAQ3))
                XQBMB3_tmd_out$ABE <- toString(unique(apdata$fAB3))
                XQBMB3_tmd_out
                
                tmd_sum <- smartbind(tmd_sum,XQBMB3_tmd_out)           
    
                
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.4.3. XQBMB9
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             
            #<<<<<<<<<<<<<<<<<<<
            # 5.4.3.2. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XQBMB9.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ 
                                                vm(AQID, AQ.giv) +
                                                vm(ABID, AB.giv):fa(LY,1) +
                                                vm(DBID, DB.giv) +
                                                at(L):idv(fU), 
                                    residual = ~ dsum(~idv(units)| L),
                                    data = apdata, 
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 5.4.3.3. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XQBMB9.asr)
            
            summary(XQBMB9.asr)$var
            
                # ##########
                # # 5.4.3.3.1. AQ
                # ##########
                # 
                # XQBMB9.Mvcov_AQ_F2C_in <- NULL
                # XQBMB9.Mvcov_AQ_F2C_in$gammas <- XQBMB9.asr$vparameters
                # XQBMB9.Mvcov_AQ_F2C_in$E <- "fAQ3"
                # 
                # XQBMB9.Mvcov_AQ <- fatocov.f(fatocov_in = XQBMB9.Mvcov_AQ_F2C_in) 
                # XQBMB9.Mvcov_AQ
                # 
                ##########
                # 5.4.3.3.2. AB
                ##########

                XQBMB9.Mvcov_AB_F2C_in <- NULL
                XQBMB9.Mvcov_AB_F2C_in$gammas <- XQBMB9.asr$vparameters
                XQBMB9.Mvcov_AB_F2C_in$fac <- "ABID"
                XQBMB9.Mvcov_AB_F2C_in$E <- "LY"

                XQBMB9.Mvcov_AB <- fatocov.f(fatocov_in = XQBMB9.Mvcov_AB_F2C_in)
                XQBMB9.Mvcov_AB
             
                plot(as.dendrogram(hclust(dist(XQBMB9.Mvcov_AB$Mr), method = 'ward.D2'),hang = -1), horiz = T, xlab = "AB")

                #<<<<<<<<<
                # 5.3.2.3.3. Cluster of rG
                #<<<<<<<<<
                
                XQBMB9.vAQ <- XQBMB9.asr$vparameters[grepl("AQID",names(XQBMB9.asr$vparameters))]
                XQBMB9.vDB <- XQBMB9.asr$vparameters[grepl("DBID",names(XQBMB9.asr$vparameters))]
                
                XQBMB9.MvcovG <- XQBMB9.Mvcov_AB$Mvcov + XQBMB9.vAQ + XQBMB9.vAQ
                XQBMB9.MrG <- cov2cor(XQBMB9.MvcovG)
                
                plot(as.dendrogram(hclust(dist(XQBMB9.MrG), method = 'ward.D2'),hang = -1), horiz = T, xlab = "AB")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 6. Save workspace
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
save.image("XQBX_210918.RData")
