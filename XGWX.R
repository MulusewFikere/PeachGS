#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# Complete script for linear mixed model analysis of peach GxE paper
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Preliminary ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.2. Load libraries ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    library(asreml)
    library(gtools)
    
    options(scipen = 999)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.3. Source functions ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 1.3.2. Estimate vcov from fa model ####
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
        # 1.3.3. Test model differences ####
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
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. prepare data ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.1. Load curated data ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    load("../4.1. Make_GIVs/Make_GIVs.RData")
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.2. pdata ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                
    pdata <- pdata_uG
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.3. sdata ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 2.3.1. wsdata ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        wsdata <- wsdata_uG

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Build GRM ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
            
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.1. Additive ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    AW.giv <- AW.giv_tune
    attr(AW.giv, "INVERSE") = TRUE
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.3.1.2.2. Dominance ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    DW.giv <- DW.giv_tune
    
    attr(DW.giv, "INVERSE") = TRUE

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Fit models ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
tmd_sum <- NULL            
          
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 4.1. Univariate ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.1. XGWU01 ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<
            # 4.1.1.1. Set data ####
            #<<<<<<<<<<<<<<<<<<
         
            apdata <- pdata
    
            sort(unique(apdata$LY))
        
            apdata$AE1 <- apdata$LY
            apdata$AE1[apdata$L == 'F'] <- 'F'
            apdata$AE1[apdata$L == 'G'] <- 'G'
            apdata$AE1[apdata$LY == 'K0' | apdata$LY == 'K1'] <- 'K01'
            
            sort(unique(apdata$AE1))
            
            apdata$L <- factor(apdata$L)
            apdata$Y <- factor(apdata$Y)
            apdata$LY <- factor(apdata$LY)
            apdata$U <- factor(apdata$U)
            
            apdata$AWID <- factor(apdata$GID,levels=rownames(AW.giv))
            apdata$DWID <- factor(apdata$GID,levels=rownames(DW.giv))
            apdata$fAE1 <- factor(apdata$AE1)
            
            unique(apdata$AE1)
            
            #<<<<<<<<<<<<<<<<<<
            # 4.1.1.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<

            asreml.options(dense = ~ vm(DWID, DW.giv))
            asreml.options(dense = ~ vm(AWID, AW.giv))
            
            XGWU01.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ vm(AWID, AW.giv) + vm(AWID, AW.giv):idv(fAE1) + 
                                                vm(DWID, DW.giv) + vm(DWID, DW.giv):idv(L) + 
                                                at(L):idv(U),
                                    residual = ~ dsum(~idv(units)| L),
                                    data= apdata,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
                        
            #<<<<<<<<<<<<<<<<<<
            # 4.1.1.3. Review results ####
            #<<<<<<<<<<<<<<<<<<
            
            wald(XGWU01.asr)
            
            summary(XGWU01.asr)$var
            
            XGWU01_analysum_temp <- data.frame("Model"= "XGWU01")
            
            XGWU01_varsum_temp <- summary(XGWU01.asr)$var
            XGWU01_varsum_select <- XGWU01_varsum_temp[!grepl("!R",rownames(XGWU01_varsum_temp),fixed=T),]
            XGWU01_varsum <- data.frame(XGWU01_varsum_select[,1])
            colnames(XGWU01_varsum) <- 'XWU01_vc'
            rownames(XGWU01_varsum) <- rownames(XGWU01_varsum_select)
            XGWU01_varsum['logl',1] <- XGWU01.asr$loglik
            XGWU01_varsum['df',1] <- length(summary(XGWU01.asr)$var$bound%in%c("P","U"))
            XGWU01_varsum['AIC',1] <- 2*(XGWU01_varsum['df',1] - XGWU01_varsum['logl',1])
            XGWU01_varsum$Component <- rownames(XGWU01_varsum)
            
            XGWU01_tmd_temp <- NULL
            XGWU01_tmd_temp$analy_name <- 'XGWU01' 
            XGWU01_tmd_temp$full_name <- 'XGWU01' 
            XGWU01_tmd_temp$fullmod.logl <- XGWU01.asr$loglik
            XGWU01_tmd_temp$fullmod.df <- sum(summary(XGWU01.asr)$var$bound%in%c("P","U"))
            XGWU01_tmd_temp$ranmod <- as.character(XGWU01.asr$call)[3]
            XGWU01_tmd_temp$EA <- toString(unique(apdata$fAE1))
            
            XGWU01_tmd_temp$Bestmod <- 'XGWU01'
            
            tmd_sum <- as.data.frame(XGWU01_tmd_temp)
            tmd_sum
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.2. XGWU02 ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.2.1. Set data ####
            #<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.2.2. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            XGWU02.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ vm(AWID, AW.giv) + vm(AWID, AW.giv):idv(fAE1) + 
                                                vm(DWID, DW.giv) + 
                                                at(L):idv(U),
                                    residual = ~ dsum(~idv(units)| L),
                                    data= apdata,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.2.3. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWU02.asr)
            
            summary(XGWU02.asr)$var
            
            XGWU02_analysum_temp <- data.frame("Model"= "XGWU02")
            
            XGWU02_varsum_temp <- summary(XGWU02.asr)$var
            XGWU02_varsum_select <- XGWU02_varsum_temp[!grepl("!R",rownames(XGWU02_varsum_temp),fixed=T),]
            XGWU02_varsum <- data.frame(XGWU02_varsum_select[,1])
            colnames(XGWU02_varsum) <- 'XWU02_vc'
            rownames(XGWU02_varsum) <- rownames(XGWU02_varsum_select)
            XGWU02_varsum['logl',1] <- XGWU02.asr$loglik
            XGWU02_varsum['df',1] <- length(summary(XGWU02.asr)$var$bound%in%c("P","U"))
            XGWU02_varsum['AIC',1] <- 2*(XGWU02_varsum['df',1] - XGWU02_varsum['logl',1])
            XGWU02_varsum$Component <- rownames(XGWU02_varsum)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.2.4. Test difference bn models
            #<<<<<<<<<<<<<<<<<<<
            
            XGWU02_tmd_in <- NULL
            XGWU02_tmd_in$analy_name <- 'XGWU02=XGWU01' 
            XGWU02_tmd_in$fullmod_name <- 'XGWU01' 
            XGWU02_tmd_in$redmod_name <- 'XGWU02'
            XGWU02_tmd_in$fullmod.asr <- XGWU01.asr
            XGWU02_tmd_in$redmod.asr <- XGWU02.asr
            XGWU02_tmd_in$test <- "boundary" #("boundary"/"other")
            
            XGWU02_tmd_out <- tmd.f(tmd_in = XGWU02_tmd_in)
            XGWU02_tmd_out
            
            XGWU02_tmd_out$Bestmod <- 'XGWU02'
            XGWU02_tmd_out$ranmod <- as.character(XGWU02.asr$call)[3]
            XGWU02_tmd_out$EA <- toString(unique(apdata$fAE1))
            
            tmd_sum <- smartbind(tmd_sum,XGWU02_tmd_out)
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.3. XGWU03
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.3.1. Set data
            #<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.3.2. Fit model
            #<<<<<<<<<<<<<<<<<<<
            
            asreml.options(dense = ~ vm(AWID, AW.giv))
            
            XGWU03.asr <- asreml(snewy ~ 1 + LY,
                                    random = ~ vm(AWID, AW.giv) + vm(AWID, AW.giv):idv(fAE1) + 
                                                at(L):idv(U),
                                    residual = ~ dsum(~idv(units)| L),
                                    data= apdata,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)
                
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.3.3. Review results
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWU03.asr)
            
            summary(XGWU03.asr)$var
            
            XGWU03_analysum_temp <- data.frame("Model"= "XGWU03")
            
            XGWU03_varsum_temp <- summary(XGWU03.asr)$var
            XGWU03_varsum_select <- XGWU03_varsum_temp[!grepl("!R",rownames(XGWU03_varsum_temp),fixed=T),]
            XGWU03_varsum <- data.frame(XGWU03_varsum_select[,1])
            colnames(XGWU03_varsum) <- 'XWU02_vc'
            rownames(XGWU03_varsum) <- rownames(XGWU03_varsum_select)
            XGWU03_varsum['logl',1] <- XGWU03.asr$loglik
            XGWU03_varsum['df',1] <- length(summary(XGWU03.asr)$var$bound%in%c("P","U"))
            XGWU03_varsum['AIC',1] <- 2*(XGWU03_varsum['df',1] - XGWU03_varsum['logl',1])
            XGWU03_varsum$Component <- rownames(XGWU03_varsum)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.1.3.4. Test difference bn models
            #<<<<<<<<<<<<<<<<<<<
            
            XGWU03_tmd_in <- NULL
            XGWU03_tmd_in$analy_name <- 'XGWU03=XGWU02' 
            XGWU03_tmd_in$fullmod_name <- 'XGWU02' 
            XGWU03_tmd_in$redmod_name <- 'XGWU03'
            XGWU03_tmd_in$fullmod.asr <- XGWU02.asr
            XGWU03_tmd_in$redmod.asr <- XGWU03.asr
            XGWU03_tmd_in$test <- "boundary" #("boundary"/"other")
            
            XGWU03_tmd_out <- tmd.f(tmd_in = XGWU03_tmd_in)
            XGWU03_tmd_out
            
            XGWU03_tmd_out$Bestmod <- 'XGWU02'
            XGWU03_tmd_out$ranmod <- as.character(XGWU03.asr$call)[3]
            XGWU03_tmd_out$EA <- toString(unique(apdata$fAE1))
            
            tmd_sum <- smartbind(tmd_sum,XGWU03_tmd_out)  
            
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 4.2. Multi-variate ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.2.XGWM02. Multivariate XGWU02 ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM02.1. Review data ####
            #<<<<<<<<<<<<<<<<<<<
            
            tapply(apdata$snewy,apdata$LY, FUN = 'sd', na.rm = T)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM02.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM02.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE1,1):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM02.3. Review fit ####
            #<<<<<<<<<<<<<<<<<<<
            
            summary(XGWM02.asr)$var
            
            XGWM02.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM02.asr)$var$bound%in%c("P","U"))
            
                #<<<<<<<<<
                # 4.2.XGWM02.3.1. Estimate genetic parameters ####
                #<<<<<<<<<
            
                    #<<<<
                    # 4.2.XGWM02.3.1.1. Additive ####
                    #<<<<
            
                    XGWM02.Mvcov_AW3C_in <- NULL
                    XGWM02.Mvcov_AW3C_in$gammas <- XGWM02.asr$vparameters
                    XGWM02.Mvcov_AW3C_in$fac <- "AWID"
                    XGWM02.Mvcov_AW3C_in$E <- "fAE1"
                    
                    XGWM02.Mvcov_AW <- fatocov.f(fatocov_in = XGWM02.Mvcov_AW3C_in) 
                    XGWM02.Mvcov_AW
                    
                    #<<<<
                    # 4.2.XGWM02.3.1.3. Total genetic ####
                    #<<<<
                    
                    XGWM02.Mvcov_GG <- NULL
                    XGWM02.Mvcov_GG$Mvcov <- XGWM02.Mvcov_AW$Mvcov
                    XGWM02.Mvcov_GG$Mr <- cov2cor(XGWM02.Mvcov_GG$Mvcov)
                    
                    #<<<< Cluster of XGWM02.Mvcov_GG$Mr
                    XGWM02_hclust <- hclust(dist(XGWM02.Mvcov_GG$Mr))
                    
                    plot(as.dendrogram(XGWM02_hclust),horiz = T)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.2.XGWM03. fa(AWID,2) ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM03.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM03.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE1,2):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)

            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM03.3. Review fit ####
            #<<<<<<<<<<<<<<<<<<<
            
            summary(XGWM03.asr)$var
            
            XGWM03.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM03.asr)$var$bound%in%c("P","U"))
            
                #<<<<<<<<<
                # 4.2.XGWM03.3.1. Estimate genetic parameters ####
                #<<<<<<<<<
            
                    #<<<<
                    # 4.2.XGWM03.3.1.1. Additive ####
                    #<<<<
            
                    XGWM03.Mvcov_AW3C_in <- NULL
                    XGWM03.Mvcov_AW3C_in$gammas <- XGWM03.asr$vparameters
                    XGWM03.Mvcov_AW3C_in$fac <- "AWID"
                    XGWM03.Mvcov_AW3C_in$E <- "fAE1"
                    
                    XGWM03.Mvcov_AW <- fatocov.f(fatocov_in = XGWM03.Mvcov_AW3C_in) 
                    XGWM03.Mvcov_AW
                    
                    #<<<<
                    # 4.2.XGWM03.3.1.2. Total genetic ####
                    #<<<<
                    
                    XGWM03.Mvcov_GG <- NULL
                    XGWM03.Mvcov_GG$Mvcov <- XGWM03.Mvcov_AW$Mvcov
                    XGWM03.Mvcov_GG$Mr <- cov2cor(XGWM03.Mvcov_GG$Mvcov)
                    
                    #<<<< Cluster of XGWM03.Mvcov_GG$Mr
                    XGWM03_hclust <- hclust(dist(XGWM03.Mvcov_GG$Mr))
                    
                    plot(as.dendrogram(XGWM03_hclust),horiz = T)
            
                #<<<<<<<<<<<<<<<<<<<
                # 4.2.XGWM03.3.2. Test difference bn models
                #<<<<<<<<<<<<<<<<<<<
                
                XGWM03_tmd_in <- NULL
                XGWM03_tmd_in$analy_name <- 'XGWM02=XGM03' 
                XGWM03_tmd_in$fullmod_name <- 'XGWM03' 
                XGWM03_tmd_in$redmod_name <- 'XGWM02'
                XGWM03_tmd_in$fullmod.asr <- XGWM03.asr
                XGWM03_tmd_in$redmod.asr <- XGWM02.asr
                XGWM03_tmd_in$test <- "conventional" #("boundary"/"other")
                
                XGWM03_tmd_out <- tmd.f(tmd_in = XGWM03_tmd_in)
                XGWM03_tmd_out
                
    #CHNOTE: XGWM03 = XGWM02            
            
                XGWM03_tmd_out$Bestmod <- 'XGWM02'
                XGWM03_tmd_out$ranmod <- as.character(XGWM03.asr$call)[3]
                XGWM03_tmd_out$EA <- toString(unique(apdata$fAE1))
                
                tmd_sum <- smartbind(tmd_sum,XGWM03_tmd_out)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.2.XGWM04. fa(AWID,1) - DWID ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM04.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM04.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE1,1):vm(AWID, AW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)

            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM04.3. Review fit ####
            #<<<<<<<<<<<<<<<<<<<
            
            summary(XGWM04.asr)$var
            
            XGWM04.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM04.asr)$var$bound%in%c("P","U"))
            
                #<<<<<<<<<
                # 4.2.XGWM04.3.1. Estimate genetic parameters ####
                #<<<<<<<<<
            
                    #<<<<
                    # 4.2.XGWM04.3.1.1. Additive ####
                    #<<<<
            
                    XGWM04.Mvcov_AW3C_in <- NULL
                    XGWM04.Mvcov_AW3C_in$gammas <- XGWM04.asr$vparameters
                    XGWM04.Mvcov_AW3C_in$fac <- "AWID"
                    XGWM04.Mvcov_AW3C_in$E <- "fAE1"
                    
                    XGWM04.Mvcov_AW <- fatocov.f(fatocov_in = XGWM04.Mvcov_AW3C_in) 
                    XGWM04.Mvcov_AW
                    
                    #<<<<
                    # 4.2.XGWM04.3.1.2. Total genetic ####
                    #<<<<
                    
                    XGWM04.Mvcov_GG <- NULL
                    XGWM04.Mvcov_GG$Mvcov <- XGWM04.Mvcov_AW$Mvcov
                    XGWM04.Mvcov_GG$Mr <- cov2cor(XGWM04.Mvcov_GG$Mvcov)
                    
                    #<<<< Cluster of XGWM04.Mvcov_GG$Mr
                    XGWM04_hclust <- hclust(dist(XGWM04.Mvcov_GG$Mr))
                    
                    plot(as.dendrogram(XGWM04_hclust),horiz = T)
            
                #<<<<<<<<<<<<<<<<<<<
                # 4.2.XGWM04.3.2. Test difference bn models
                #<<<<<<<<<<<<<<<<<<<
                
                XGWM04_tmd_in <- NULL
                XGWM04_tmd_in$analy_name <- 'XGWM04=XGWM02' 
                XGWM04_tmd_in$fullmod_name <- 'XGWM02' 
                XGWM04_tmd_in$redmod_name <- 'XGWM04'
                XGWM04_tmd_in$fullmod.asr <- XGWM02.asr
                XGWM04_tmd_in$redmod.asr <- XGWM04.asr
                XGWM04_tmd_in$test <- "conventional" #("boundary"/"other")
                
                XGWM04_tmd_out <- tmd.f(tmd_in = XGWM04_tmd_in)
                
                XGWM04_tmd_out
                
#CHNOTE: XGWM04 = XGWM02            
            
                2*(sum(summary(XGWM04.asr)$var$bound%in%c("P","U")) - XGWM04.asr$loglik)
                2*(sum(summary(XGWM03.asr)$var$bound%in%c("P","U")) - XGWM03.asr$loglik)
                XGWM04_tmd_out$Bestmod <- 'XGWM04'
                XGWM04_tmd_out$ranmod <- as.character(XGWM04.asr$call)[3]
                XGWM04_tmd_out$EA <- toString(unique(apdata$fAE1))
                
                tmd_sum <- smartbind(tmd_sum,XGWM04_tmd_out)                
 
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
        # 4.2.XGWM05. ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM05.1. Set data ####
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XGWM02.Mvcov_AW$Mvcov)),hang = -1)
            
            apdata$AE2 <- as.character(apdata$AE1)
            apdata$AE2[apdata$AE1 == 'F' | apdata$AE1 == 'K2'] <- 'FK2'
            apdata$fAE2 <- factor(apdata$AE2) 
            
            unique(apdata$fAE2)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM05.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            asreml.options(dense = ~ vm(AWID, AW.giv))
            asreml.options(dense = ~ vm(DWID, DW.giv))
            
            XGWM05.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE2,1):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)

            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM05.3. Review results ####
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWM05.asr)
            
            summary(XGWM05.asr)$var
            
            XGWM05.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM05.asr)$var$bound%in%c("P","U"))
            
            XGWM05.Mvcov_AW2C_in <- NULL
            XGWM05.Mvcov_AW2C_in$gammas <- XGWM05.asr$vparameters
            XGWM05.Mvcov_AW2C_in$fac <- "AWID"
            XGWM05.Mvcov_AW2C_in$E <- "fAE2"
            
            XGWM05.Mvcov_AW <- fatocov.f(fatocov_in = XGWM05.Mvcov_AW2C_in) 
            XGWM05.Mvcov_AW
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM05.4. Test difference bn models ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM05_tmd_in <- NULL
            XGWM05_tmd_in$analy_name <- 'XGWM05=XGWM02' 
            XGWM05_tmd_in$fullmod_name <- 'XGWM02' 
            XGWM05_tmd_in$redmod_name <- 'XGWM05'
            XGWM05_tmd_in$fullmod.asr <- XGWM02.asr
            XGWM05_tmd_in$redmod.asr <- XGWM05.asr
            XGWM05_tmd_in$test <- "conventional" #("boundary"/"other")
            
            XGWM05_tmd_out <- tmd.f(tmd_in = XGWM05_tmd_in)
            XGWM05_tmd_out

#CHNOTE: XGWM05 < XGWM02
#CHNOTE: Best model = XGWM02

            XGWM05_tmd_out$Bestmod <- 'XGWM05'
            XGWM05_tmd_out$ranmod <- as.character(XGWM05.asr$call)[3]
            XGWM05_tmd_out$EA <- toString(unique(apdata$fAE2))
            
            tmd_sum <- smartbind(tmd_sum,XGWM05_tmd_out)
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.2.XGWM06. ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM06.1. Set data ####
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XGWM02.Mvcov_AW$Mvcov)),hang = -1)
            
            apdata$AE3 <- apdata$AE1
            apdata$AE3[apdata$AE1 == 'S0' | apdata$AE1 == 'S1'] <- 'S01'
            apdata$fAE3 <- factor(apdata$AE3) 
            
            unique(apdata$fAE3)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM06.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            asreml.options(dense = ~ vm(AWID, AW.giv))
            asreml.options(dense = ~ vm(DWID, DW.giv))
            
            
            XGWM06.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE3,1):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM06.3. Review results ####
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWM06.asr)
            
            summary(XGWM06.asr)$var
            
            XGWM06.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM06.asr)$var$bound%in%c("P","U"))
            
            XGWM06.Mvcov_AW2C_in <- NULL
            XGWM06.Mvcov_AW2C_in$gammas <- XGWM06.asr$vparameters
            XGWM06.Mvcov_AW2C_in$fac <- "AWID"
            XGWM06.Mvcov_AW2C_in$E <- "fAE3"
            
            XGWM06.Mvcov_AW <- fatocov.f(fatocov_in = XGWM06.Mvcov_AW2C_in) 
            XGWM06.Mvcov_AW
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM06.4. Test difference bn models ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM06_tmd_in <- NULL
            XGWM06_tmd_in$analy_name <- 'XGWM06=XGWM05' 
            XGWM06_tmd_in$fullmod_name <- 'XGWM05' 
            XGWM06_tmd_in$redmod_name <- 'XGWM06'
            XGWM06_tmd_in$fullmod.asr <- XGWM05.asr
            XGWM06_tmd_in$redmod.asr <- XGWM06.asr
            XGWM06_tmd_in$test <- "conventional" #("boundary"/"other")
            
            XGWM06_tmd_out <- tmd.f(tmd_in = XGWM06_tmd_in)
            XGWM06_tmd_out

#CHNOTE: XGWM06 < XGWM02
#CHNOTE: Best model = XGWM02

            XGWM06_tmd_out$Bestmod <- 'XGWM06'
            XGWM06_tmd_out$ranmod <- as.character(XGWM06.asr$call)[3]
            XGWM06_tmd_out$EA <- toString(unique(apdata$fAE3))
            
            tmd_sum <- smartbind(tmd_sum,XGWM06_tmd_out)            
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
        # 4.2.XGWM07. ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM07.1. Set data ####
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XGWM02.Mvcov_AW$Mvcov)),hang = -1)
            
            apdata$AE4 <- apdata$AE1
            apdata$AE4[apdata$AE3 == 'K01' | apdata$AE3 == 'S2'] <- 'K01S2'
            apdata$fAE4 <- factor(apdata$AE4) 
            
            unique(apdata$fAE4)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM07.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            asreml.options(dense = ~ vm(AWID, AW.giv))
            
            XGWM07.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE4,1):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM07.3. Review results ####
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWM07.asr)
            
            summary(XGWM07.asr)$var
            
            XGWM07.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM07.asr)$var$bound%in%c("P","U"))
            
            XGWM07.Mvcov_AW2C_in <- NULL
            XGWM07.Mvcov_AW2C_in$gammas <- XGWM07.asr$vparameters
            XGWM07.Mvcov_AW2C_in$fac <- "AWID"
            XGWM07.Mvcov_AW2C_in$E <- "fAE4"
            
            XGWM07.Mvcov_AW <- fatocov.f(fatocov_in = XGWM07.Mvcov_AW2C_in) 
            XGWM07.Mvcov_AW
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM07.4. Test difference bn models ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM07_tmd_in <- NULL
            XGWM07_tmd_in$analy_name <- 'XGWM02 = XGWM07' 
            XGWM07_tmd_in$fullmod_name <- 'XGWM02' 
            XGWM07_tmd_in$redmod_name <- 'XGWM07'
            XGWM07_tmd_in$fullmod.asr <- XGWM02.asr
            XGWM07_tmd_in$redmod.asr <- XGWM07.asr
            XGWM07_tmd_in$test <- "conventional" #("boundary"/"other")
            
            XGWM07_tmd_out <- tmd.f(tmd_in = XGWM07_tmd_in)
            XGWM07_tmd_out

#CHNOTE: XGWM07 < XGWM02
#CHNOTE: Best model = XGWM02

            XGWM07_tmd_out$Bestmod <- 'XGWM07'
            XGWM07_tmd_out$ranmod <- as.character(XGWM07.asr$call)[3]
            XGWM07_tmd_out$EA <- toString(unique(apdata$fAE4))
            
            tmd_sum <- smartbind(tmd_sum,XGWM07_tmd_out)             
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
        # 4.2.XGWM08. ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM08.1. Set data ####
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XGWM07.Mvcov_AW$Mvcov)),hang = -1)
            
            apdata$AE5 <- apdata$AE4
            apdata$AE5[apdata$AE4 == 'F' | apdata$AE4 == 'K2'] <- 'FK2'
            apdata$fAE5 <- factor(apdata$AE5) 
            
            unique(apdata$fAE5)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM08.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            asreml.options(dense = ~ vm(AWID, AW.giv))
            
            XGWM08.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE5,1):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM08.3. Review results ####
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWM08.asr)
            
            summary(XGWM08.asr)$var
            
            XGWM08.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM08.asr)$var$bound%in%c("P","U"))
            
            XGWM08.Mvcov_AW2C_in <- NULL
            XGWM08.Mvcov_AW2C_in$gammas <- XGWM08.asr$vparameters
            XGWM08.Mvcov_AW2C_in$fac <- "AWID"
            XGWM08.Mvcov_AW2C_in$E <- "fAE5"
            
            XGWM08.Mvcov_AW <- fatocov.f(fatocov_in = XGWM08.Mvcov_AW2C_in) 
            XGWM08.Mvcov_AW
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM08.4. Test difference bn models ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM08_tmd_in <- NULL
            XGWM08_tmd_in$analy_name <- 'XGWM08 = XGWM07' 
            XGWM08_tmd_in$fullmod_name <- 'XGWM07' 
            XGWM08_tmd_in$redmod_name <- 'XGWM08'
            XGWM08_tmd_in$fullmod.asr <- XGWM07.asr
            XGWM08_tmd_in$redmod.asr <- XGWM08.asr
            XGWM08_tmd_in$test <- "conventional" #("boundary"/"other")
            
            XGWM08_tmd_out <- tmd.f(tmd_in = XGWM08_tmd_in)
            XGWM08_tmd_out

#CHNOTE: XGWM08 < XGWM02
#CHNOTE: Best model = XGWM02

            XGWM08_tmd_out$Bestmod <- 'XGWM07'
            XGWM08_tmd_out$ranmod <- as.character(XGWM08.asr$call)[3]
            XGWM08_tmd_out$EA <- toString(unique(apdata$fAE5))
            
            tmd_sum <- smartbind(tmd_sum,XGWM08_tmd_out)             

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
        # 4.2.XGWM09. ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM09.1. Set data ####
            #<<<<<<<<<<<<<<<<<<<
            
            plot(hclust(dist(XGWM08.Mvcov_AW$Mvcov)),hang = -1)
            
            apdata$AE6 <- apdata$AE4
            apdata$AE6[apdata$AE4 == 'F' | apdata$AE4 == 'G'] <- 'FG'
            apdata$fAE6 <- factor(apdata$AE6) 
            
            unique(apdata$fAE6)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM09.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<
            
            asreml.options(dense = ~ vm(AWID, AW.giv))
            
            XGWM09.asr <- asreml(snewy ~ 1 + LY,
                                 random = ~ fa(fAE6,1):vm(AWID, AW.giv) +
                                            vm(DWID, DW.giv) +
                                            at(L):idv(U),
                                 residual = ~ dsum(~idv(units)| L),
                                 data= apdata,
                                 na.action = na.method(y = "include", x = "include"),
                                 workspace = 6e+08,
                                 maxiter=50)
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM09.3. Review results ####
            #<<<<<<<<<<<<<<<<<<<
            
            wald(XGWM09.asr)
            
            summary(XGWM09.asr)$var
            
            XGWM09.asr$loglik
            
            #<<<< Degrees of freedom
            sum(summary(XGWM09.asr)$var$bound%in%c("P","U"))
            
            XGWM09.Mvcov_AW2C_in <- NULL
            XGWM09.Mvcov_AW2C_in$gammas <- XGWM09.asr$vparameters
            XGWM09.Mvcov_AW2C_in$fac <- "AWID"
            XGWM09.Mvcov_AW2C_in$E <- "fAE6"
            
            XGWM09.Mvcov_AW <- fatocov.f(fatocov_in = XGWM09.Mvcov_AW2C_in) 
            XGWM09.Mvcov_AW
            
            #<<<<<<<<<<<<<<<<<<<
            # 4.2.XGWM09.4. Test difference bn models ####
            #<<<<<<<<<<<<<<<<<<<
            
            XGWM09_tmd_in <- NULL
            XGWM09_tmd_in$analy_name <- 'XGWM09 = XGWM07' 
            XGWM09_tmd_in$fullmod_name <- 'XGWM07' 
            XGWM09_tmd_in$redmod_name <- 'XGWM09'
            XGWM09_tmd_in$fullmod.asr <- XGWM07.asr
            XGWM09_tmd_in$redmod.asr <- XGWM09.asr
            XGWM09_tmd_in$test <- "conventional" #("boundary"/"other")
            
            XGWM09_tmd_out <- tmd.f(tmd_in = XGWM09_tmd_in)
            XGWM09_tmd_out

#CHNOTE: XGWM09 < XGWM02
#CHNOTE: Best model = XGWM02

            XGWM09_tmd_out$Bestmod <- 'XGWM07'
            XGWM09_tmd_out$ranmod <- as.character(XGWM09.asr$call)[3]
            XGWM09_tmd_out$EA <- toString(unique(apdata$fAE6))
            
            tmd_sum <- smartbind(tmd_sum,XGWM09_tmd_out)             
            
            write.csv(tmd_sum,"tmd_sum_211128.csv",row.names = F)
                    
        
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 5. Save workspace
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
save.image("XGWX_210918.RData")
            
        
            