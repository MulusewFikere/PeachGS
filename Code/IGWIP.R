#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# Complete script for linear mixed model analysis of peach GxE paper
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Preliminary ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.1. Load libraries ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    library(asreml) 
    library(gtools)
    options(scipen = 999)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.2. Source functions ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 1.2.1. Estimate vcov from fa model ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        fatocov <- function(fatocov_in)
        {
            fc.Mpsi <- diag(fatocov_in$vparameters[grepl(fatocov_in$fac,names(fatocov_in$vparameters),fixed=T) &
                                                  grepl(fatocov_in$E,names(fatocov_in$vparameters),fixed=T) &
                                                  grepl("var",names(fatocov_in$vparameters))])
            rownames(fc.Mpsi) <- colnames(fc.Mpsi) <- sub(paste(fatocov_in$E,".",sep=""),"",
                                                          sub(".var","",
                                                              sapply(
                                                                  strsplit(
                                                                      names(fatocov_in$vparameters[
                                                                          grepl(fatocov_in$fac,names(fatocov_in$vparameters),fixed=T) &
                                                                              grepl(fatocov_in$E,names(fatocov_in$vparameters),fixed=T) &
                                                                              grepl("var",names(fatocov_in$vparameters))]),
                                                                      split = "!",fixed = T),
                                                                  "[[",2)
                                                          )
            )
            
            fc.Mlam <- matrix(fatocov_in$vparameters[grepl(fatocov_in$fac,names(fatocov_in$vparameters),fixed=T) &
                                                    grepl(fatocov_in$E,names(fatocov_in$vparameters),fixed=T) &
                                                    !grepl("var",names(fatocov_in$vparameters))],
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
        # 1.2.2. Test model differences ####
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
        
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. Prepare data ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.1. Load RData ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    load("../4.1. Make_GIVs/Make_GIVs.RData")
                
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.2. Prepare pdata ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 2.2.1. pdata ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        apdata_prep <- pdata_uG
                    
        tapply(apdata_prep$value,apdata_prep$LY,FUN = var,na.rm = T)
        tapply(apdata_prep$snewy,apdata_prep$LY,FUN = var,na.rm = T)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 2.2.2. GID list ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         
        GID_in_apdata <- unique(apdata_prep$GID)
        length(GID_in_apdata)
        
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.3. Select and modify sdata ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 2.3.1.  ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
            #<<<<<<<<<<<<<<<<<<<
            # 2.3.1.1. Prepare wsdata ####
            #<<<<<<<<<<<<<<<<<<<
         
            wsdata_prep <- wsdata_uG
            nrow(wsdata_prep)
        
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Individual location models with location specific GRM ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.1. Fresno (F) ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.1.1. Prepare data ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.1.1.1 Select data ####
            #<<<<<<<<<<<<<<<<<<<
            
            apdata_F <- apdata_prep[apdata_prep$L%in%c("F"),]
            
            GID_apdata_F <- unique(apdata_F$GID)
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.1.1.2. Build AGRM and DGRM ####
            #<<<<<<<<<<<<<<<<<<<
            
                #<<<<<<<<<
                # 3.1.1.2.1. Additive
                #<<<<<<<<<
            
                    #<<<<#
                    # 3.1.1.2.1.1 AG_F
                    #<<<<#
            
                    AG_F.giv <- AW.giv_F_tune
                    
                #<<<<<<<<<
                # 3.1.1.2.2. Dominance ####
                #<<<<<<<<<
                    
                    #<<<<
                    # 3.1.1.2.2.1. DG_F
                    #<<<<
                
                    DG_F.giv <- DW.giv_F_tune
                    
            #<<<<<<<<<<<<<<<<<<<
            # 3.1.1.3. Set factors ####
            #<<<<<<<<<<<<<<<<<<<
            
            ### Factorial
            apdata_F$Y <- factor(apdata_F$Y)
            apdata_F$AGID <- factor(apdata_F$GID,levels = colnames(AG_F.giv))
            apdata_F$DGID <- factor(apdata_F$GID,levels = colnames(DG_F.giv))
            apdata_F$U <- factor(apdata_F$U)
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.1.2. Univariate models ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.1.2.1. FGWFU01 ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.1.2.1.1. Fit model ####
                #<<<<<<<<<
            
                options(scipen = 999)
                FGWFU01.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_F.giv) + vm(AGID,AG_F.giv):idv(Y) + 
                                                    vm(DGID,DG_F.giv) + vm(DGID,DG_F.giv):idv(Y) + 
                                                    U, 
                                        data= apdata_F,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.1.2.1.2. Review model fit ####
                #<<<<<<<<<
                
                    #<<<<
                    # 3.1.2.1.2.1. Residuals
                    #<<<<
                
                    plot(FGWFU01.asr)
                
                    apdata_res_F <- cbind(apdata_F,FGWFU01.asr$residuals)
                    colnames(apdata_res_F)[ncol(apdata_res_F)] <-  "res"
                    
                    apdata_res_F$stdres <- apdata_res_F$res/sqrt(var(apdata_res_F$res,na.rm = T))
                    hist(apdata_res_F$stdres, main = "apdata_stdres_F", xlab = "std residuals distribution", col= "sky blue")
                    
                    apdata_res_F[abs(apdata_res_F$stdres) > 3.5 & !is.na(apdata_res_F$stdres),]
                    
                    #<<<<
                    # 3.1.2.1.2.2. Model parameters
                    #<<<<
                
                    wald(FGWFU01.asr)
                    
                    FGWFU01_vc <- summary(FGWFU01.asr)$var
                    FGWFU01_vc
                    
                #<<<<<<<<<
                # 3.1.2.1.3. Collate results ####
                #<<<<<<<<<
                    
                FGWFU01_analysum <- data.frame("Model_ID"="FGWFU01")
                    
                FGWFU01_analysum$logl <-   FGWFU01.asr$logl
                FGWFU01_analysum$df <-  sum(FGWFU01_vc[,'bound']%in%c("U","P"))
                
                FGWFU01_analysum$pY <- round(wald(FGWFU01.asr)[grepl("Y",rownames(wald(FGWFU01.asr))),4],3)
                
                FGWFU01_vc <- summary(FGWFU01.asr)$var
                
                FGWFU01_analysum$vAGID <- round(FGWFU01_vc[grepl("AGID",rownames(FGWFU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU01_vc)),1],3)
                    
                FGWFU01_analysum$vAGIDxY <- round(FGWFU01_vc[grepl("AGID",rownames(FGWFU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU01_vc)),1],3)
                
                FGWFU01_analysum$vDGID <- round(FGWFU01_vc[grepl("DGID",rownames(FGWFU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU01_vc)),1],3)
                    
                FGWFU01_analysum$vDGIDxY <- round(FGWFU01_vc[grepl("DGID",rownames(FGWFU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU01_vc)),1],3)
                
                FGWFU01_analysum$vU <- round(FGWFU01_vc[grepl("U",rownames(FGWFU01_vc),fixed = T),1],3)
                
                FGWFU01_analysum$vR <- round(FGWFU01_vc[grepl("R",rownames(FGWFU01_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- FGWFU01_analysum
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.1.2.2. FGWFU02 ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.1.2.2.1. Fit model
                #<<<<<<<<<
            
                FGWFU02.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_F.giv) + vm(AGID,AG_F.giv):idv(Y) + 
                                                    vm(DGID,DG_F.giv) + 
                                                    U, 
                                        data= apdata_F,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.1.2.2.2. Review model fit ####
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.1.2.2.2.2. Model parameters
                    #<<<<
                
                    wald(FGWFU02.asr)
                    
                    FGWFU02_vc <- summary(FGWFU02.asr)$var
                    FGWFU02_vc
                    
                #<<<<<<<<<
                # 3.1.2.2.3. Collate results ####
                #<<<<<<<<<
                
                FGWFU02_analysum <- data.frame("Model_ID"="FGWFU02")
                    
                FGWFU02_analysum$logl <-   FGWFU02.asr$logl
                FGWFU02_analysum$df <-  sum(FGWFU01_vc[,'bound']%in%c("U","P"))
                
                FGWFU02_analysum$pY <- round(wald(FGWFU02.asr)[grepl("Y",rownames(wald(FGWFU02.asr))),4],3)
                
                FGWFU02_analysum$vAGID <- round(FGWFU02_vc[grepl("AGID",rownames(FGWFU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU02_vc)),1],3)
                    
                FGWFU02_analysum$vAGIDxY <- round(FGWFU02_vc[grepl("AGID",rownames(FGWFU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU02_vc)),1],3)
                
                FGWFU02_analysum$vDGID <- round(FGWFU02_vc[grepl("DGID",rownames(FGWFU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU02_vc)),1],3)
                    
                FGWFU02_analysum$vDGIDxY <- round(FGWFU02_vc[grepl("DGID",rownames(FGWFU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU02_vc)),1],3)
                
                FGWFU02_analysum$vU <- round(FGWFU02_vc[grepl("U",rownames(FGWFU02_vc),fixed = T),1],3)
                
                FGWFU02_analysum$vR <- round(FGWFU02_vc[grepl("R",rownames(FGWFU02_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,FGWFU02_analysum)
                
                #<<<<<<<<<
                # 3.1.2.2.4. test models if required  ####
                #<<<<<<<<<
                
                FGWFU02_tmd_in <- NULL
                FGWFU02_tmd_in$analy_name <- "FU01=FU02"
                FGWFU02_tmd_in$fullmod_name <- "FGWFU01"
                FGWFU02_tmd_in$redmod_name <- "FGWFU02"
                FGWFU02_tmd_in$fullmod.asr <- FGWFU01.asr
                FGWFU02_tmd_in$redmod.asr <- FGWFU02.asr
                FGWFU02_tmd_in$test <- "boundary" #("boundary"/"other")
                
                FGWFU02_tmd_out <- tmd.f(tmd_in =  FGWFU02_tmd_in)
                FGWFU02_tmd_out
            
                FGWFU02_tmd_out$Bestmod <- 'FGWFU02'
                
                IGWIPU_tmd_sum <- FGWFU02_tmd_out
                
                #<<<<<<<<<
                # 3.1.2.2.5. 
                #<<<<<<<<<
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,FGWFU02_analysum)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.1.3.2. FGWFU03 ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.1.3.2.1. Fit model ####
                #<<<<<<<<<
            
                FGWFU03.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_F.giv) +  
                                                    vm(DGID,DG_F.giv) + 
                                                    U, 
                                        data= apdata_F,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.1.3.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.1.3.2.2.2. Model parameters
                    #<<<<
                
                    wald(FGWFU03.asr)
                    
                    FGWFU03_vc <- summary(FGWFU03.asr)$var
                    FGWFU03_vc
                    
                #<<<<<<<<<
                # 3.1.3.2.3. Collate results
                #<<<<<<<<<
                
                FGWFU03_analysum <- data.frame("Model_ID"="FGWFU03")
                    
                FGWFU03_analysum$logl <-   FGWFU03.asr$logl
                FGWFU03_analysum$df <-  sum(FGWFU01_vc[,'bound']%in%c("U","P"))
                
                FGWFU03_analysum$pY <- round(wald(FGWFU03.asr)[grepl("Y",rownames(wald(FGWFU03.asr))),4],3)
                
                FGWFU03_analysum$vAGID <- round(FGWFU03_vc[grepl("AGID",rownames(FGWFU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU03_vc)),1],3)
                    
                FGWFU03_analysum$vAGIDxY <- round(FGWFU03_vc[grepl("AGID",rownames(FGWFU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU03_vc)),1],3)
                
                FGWFU03_analysum$vDGID <- round(FGWFU03_vc[grepl("DGID",rownames(FGWFU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU03_vc)),1],3)
                    
                FGWFU03_analysum$vDGIDxY <- round(FGWFU03_vc[grepl("DGID",rownames(FGWFU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU03_vc)),1],3)
                
                FGWFU03_analysum$vU <- round(FGWFU03_vc[grepl("U",rownames(FGWFU03_vc),fixed = T),1],3)
                
                FGWFU03_analysum$vR <- round(FGWFU03_vc[grepl("R",rownames(FGWFU03_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,FGWFU03_analysum)
                
                #<<<<<<<<<
                # 3.1.3.2.4. Test models #### 
                #<<<<<<<<<
                
                FGWFU03_tmd_in <- NULL
                FGWFU03_tmd_in$analy_name <- "FU02=FU03"
                FGWFU03_tmd_in$fullmod_name <- "FGWFU02"
                FGWFU03_tmd_in$redmod_name <- "FGWFU03"
                FGWFU03_tmd_in$fullmod.asr <- FGWFU02.asr
                FGWFU03_tmd_in$redmod.asr <- FGWFU03.asr
                FGWFU03_tmd_in$test <- "boundary" #("boundary"/"other")
                
                FGWFU03_tmd_out <- tmd.f(tmd_in =  FGWFU03_tmd_in)
                FGWFU03_tmd_out
            
                FGWFU03_tmd_out$Bestmod <- 'FGWFU03'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,FGWFU03_tmd_out)

            #<<<<<<<<<<<<<<<<<<<
            # 3.1.4.2. FGWFU04
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.1.4.2.1. Fit model
                #<<<<<<<<<
            
                FGWFU04.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_F.giv) +
                                                    U, 
                                        data= apdata_F,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.1.4.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.1.4.2.2.2. Model parameters
                    #<<<<
                
                    wald(FGWFU04.asr)
                    
                    FGWFU04_vc <- summary(FGWFU04.asr)$var
                    FGWFU04_vc
                    
                #<<<<<<<<<
                # 3.1.3.2.3. Collate results
                #<<<<<<<<<
                
                FGWFU04_analysum <- data.frame("Model_ID"="FGWFU04")
                    
                FGWFU04_analysum$logl <-   FGWFU04.asr$logl
                FGWFU04_analysum$df <-  sum(FGWFU01_vc[,'bound']%in%c("U","P"))
                
                FGWFU04_analysum$pY <- round(wald(FGWFU04.asr)[grepl("Y",rownames(wald(FGWFU04.asr))),4],3)
                
                FGWFU04_analysum$vAGID <- round(FGWFU04_vc[grepl("AGID",rownames(FGWFU04_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU04_vc)),1],3)
                    
                FGWFU04_analysum$vAGIDxY <- round(FGWFU04_vc[grepl("AGID",rownames(FGWFU04_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU04_vc)),1],3)
                
                FGWFU04_analysum$vDGID <- round(FGWFU04_vc[grepl("DGID",rownames(FGWFU04_vc),fixed = T)
                                                                    & !grepl('Y',rownames(FGWFU04_vc)),1],3)
                    
                FGWFU04_analysum$vDGIDxY <- round(FGWFU04_vc[grepl("DGID",rownames(FGWFU04_vc),fixed = T)
                                                                    & grepl('Y',rownames(FGWFU04_vc)),1],3)
                
                FGWFU04_analysum$vU <- round(FGWFU04_vc[grepl("U",rownames(FGWFU04_vc),fixed = T),1],3)
                
                FGWFU04_analysum$vR <- round(FGWFU04_vc[grepl("R",rownames(FGWFU04_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,FGWFU04_analysum)
                
                #<<<<<<<<<
                # 3.1.4.2.4. test models if required 
                #<<<<<<<<<
                
                FGWFU04_tmd_in <- NULL
                FGWFU04_tmd_in$analy_name <- "FU03=FU04"
                FGWFU04_tmd_in$fullmod_name <- "FGWFU03"
                FGWFU04_tmd_in$redmod_name <- "FGWFU04"
                FGWFU04_tmd_in$fullmod.asr <- FGWFU03.asr
                FGWFU04_tmd_in$redmod.asr <- FGWFU04.asr
                FGWFU04_tmd_in$test <- "boundary" #("boundary"/"other")
                
                FGWFU04_tmd_out <- tmd.f(tmd_in =  FGWFU04_tmd_in)
                FGWFU04_tmd_out
            
                FGWFU04_tmd_out$Bestmod <- 'FGWFU04'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,FGWFU04_tmd_out)
                
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.2. colleGe station (G)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.2.1. Prepare data
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.2.1.1 Select data
            #<<<<<<<<<<<<<<<<<<<
            
            apdata_G <- apdata_prep[apdata_prep$L%in%c("G"),]
            
            GID_apdata_G <- unique(apdata_G$GID)
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.2.1.2. Build AGRM and DGRM
            #<<<<<<<<<<<<<<<<<<<
            
            #<<<< Additive
            AG_G.giv <- AW.giv_G_tune
                    
            #<<<< Dominance
            DG_G.giv <- DW.giv_G_tune
                    
            #<<<<<<<<<<<<<<<<<<<
            # 3.2.1.3. Set factors
            #<<<<<<<<<<<<<<<<<<<
            
            ### Factorial
            apdata_G$Y <- factor(apdata_G$Y)
            apdata_G$AGID <- factor(apdata_G$GID,levels = colnames(AG_G.giv))
            apdata_G$DGID <- factor(apdata_G$GID,levels = colnames(DG_G.giv))
            apdata_G$U <- factor(apdata_G$U)
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.2.2. Univariate models
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.2.2.1. GGWGU01
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.2.2.1.1. Fit model
                #<<<<<<<<<
            
                options(scipen = 999)
                GGWGU01.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_G.giv) + vm(AGID,AG_G.giv):idv(Y) + 
                                                    vm(DGID,DG_G.giv) + vm(DGID,DG_G.giv):idv(Y) + 
                                                    U, 
                                        data= apdata_G,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.2.2.1.2. Review model fit
                #<<<<<<<<<
                
                    #<<<<
                    # 3.2.2.1.2.1. Residuals
                    #<<<<
                
                    plot(GGWGU01.asr)
                
                    apdata_res_G <- cbind(apdata_G,GGWGU01.asr$residuals)
                    colnames(apdata_res_G)[ncol(apdata_res_G)] <-  "res"
                    
                    apdata_res_G$stdres <- apdata_res_G$res/sqrt(var(apdata_res_G$res,na.rm = T))
                    hist(apdata_res_G$stdres, main = "apdata_stdres_G", xlab = "std residuals distribution", col= "sky blue")
                    
                    apdata_res_G[abs(apdata_res_G$stdres) > 3.5 & !is.na(apdata_res_G$stdres),]
                    
                    #<<<<
                    # 3.2.2.1.2.2. Model parameters
                    #<<<<
                
                    wald(GGWGU01.asr)
                    
                    GGWGU01_vc <- summary(GGWGU01.asr)$var
                    GGWGU01_vc
                    
                #<<<<<<<<<
                # 3.2.2.1.3. Collate results
                #<<<<<<<<<
                    
                GGWGU01_analysum <- data.frame("Model_ID"="GGWGU01")
                    
                GGWGU01_analysum$logl <-   GGWGU01.asr$logl
                GGWGU01_analysum$df <-  sum(GGWGU01_vc[,'bound']%in%c("U","P"))
                
                GGWGU01_analysum$pY <- round(wald(GGWGU01.asr)[grepl("Y",rownames(wald(GGWGU01.asr))),4],3)
                
                GGWGU01_vc <- summary(GGWGU01.asr)$var
                
                GGWGU01_analysum$vAGID <- round(GGWGU01_vc[grepl("AGID",rownames(GGWGU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU01_vc)),1],3)
                    
                GGWGU01_analysum$vAGIDxY <- round(GGWGU01_vc[grepl("AGID",rownames(GGWGU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU01_vc)),1],3)
                
                GGWGU01_analysum$vDGID <- round(GGWGU01_vc[grepl("DGID",rownames(GGWGU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU01_vc)),1],3)
                    
                GGWGU01_analysum$vDGIDxY <- round(GGWGU01_vc[grepl("DGID",rownames(GGWGU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU01_vc)),1],3)
                
                GGWGU01_analysum$vU <- round(GGWGU01_vc[grepl("U",rownames(GGWGU01_vc),fixed = T),1],3)
                
                GGWGU01_analysum$vR <- round(GGWGU01_vc[grepl("R",rownames(GGWGU01_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,GGWGU01_analysum)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.2.2.2. GGWGU02
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.2.2.2.1. Fit model
                #<<<<<<<<<
            
                GGWGU02.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_G.giv) + vm(AGID,AG_G.giv):idv(Y) + 
                                                    vm(DGID,DG_G.giv) + 
                                                    U, 
                                        data= apdata_G,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.2.2.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.2.2.2.2.2. Model parameters
                    #<<<<
                
                    wald(GGWGU02.asr)
                    
                    GGWGU02_vc <- summary(GGWGU02.asr)$var
                    GGWGU02_vc
                    
                #<<<<<<<<<
                # 3.2.2.2.3. Collate results
                #<<<<<<<<<
                
                GGWGU02_analysum <- data.frame("Model_ID"="GGWGU02")
                    
                GGWGU02_analysum$logl <-   GGWGU02.asr$logl
                GGWGU02_analysum$df <-  sum(GGWGU01_vc[,'bound']%in%c("U","P"))
                
                GGWGU02_analysum$pY <- round(wald(GGWGU02.asr)[grepl("Y",rownames(wald(GGWGU02.asr))),4],3)
                
                GGWGU02_analysum$vAGID <- round(GGWGU02_vc[grepl("AGID",rownames(GGWGU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU02_vc)),1],3)
                    
                GGWGU02_analysum$vAGIDxY <- round(GGWGU02_vc[grepl("AGID",rownames(GGWGU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU02_vc)),1],3)
                
                GGWGU02_analysum$vDGID <- round(GGWGU02_vc[grepl("DGID",rownames(GGWGU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU02_vc)),1],3)
                    
                GGWGU02_analysum$vDGIDxY <- round(GGWGU02_vc[grepl("DGID",rownames(GGWGU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU02_vc)),1],3)
                
                GGWGU02_analysum$vU <- round(GGWGU02_vc[grepl("U",rownames(GGWGU02_vc),fixed = T),1],3)
                
                GGWGU02_analysum$vR <- round(GGWGU02_vc[grepl("R",rownames(GGWGU02_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,GGWGU02_analysum)
                
                #<<<<<<<<<
                # 3.2.2.2.4. test models if required 
                #<<<<<<<<<
                
                GGWGU02_tmd_in <- NULL
                GGWGU02_tmd_in$analy_name <- "GU01=GU02"
                GGWGU02_tmd_in$fullmod_name <- "GGWGU01"
                GGWGU02_tmd_in$redmod_name <- "GGWGU02"
                GGWGU02_tmd_in$fullmod.asr <- GGWGU01.asr
                GGWGU02_tmd_in$redmod.asr <- GGWGU02.asr
                GGWGU02_tmd_in$test <- "boundary" #("boundary"/"other")
                
                GGWGU02_tmd_out <- tmd.f(tmd_in =  GGWGU02_tmd_in)
                GGWGU02_tmd_out
            
                GGWGU02_tmd_out$Bestmod <- 'GGWGU02'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,GGWGU02_tmd_out)
                
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.2.3.2. GGWGU03
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.2.3.2.1. Fit model
                #<<<<<<<<<
            
                GGWGU03.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_G.giv) +  
                                                    vm(DGID,DG_G.giv) + 
                                                    U, 
                                        data= apdata_G,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.2.3.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.2.3.2.2.2. Model parameters
                    #<<<<
                
                    wald(GGWGU03.asr)
                    
                    GGWGU03_vc <- summary(GGWGU03.asr)$var
                    GGWGU03_vc
                    
                #<<<<<<<<<
                # 3.2.3.2.3. Collate results
                #<<<<<<<<<
                
                GGWGU03_analysum <- data.frame("Model_ID"="GGWGU03")
                    
                GGWGU03_analysum$logl <-   GGWGU03.asr$logl
                GGWGU03_analysum$df <-  sum(GGWGU01_vc[,'bound']%in%c("U","P"))
                
                GGWGU03_analysum$pY <- round(wald(GGWGU03.asr)[grepl("Y",rownames(wald(GGWGU03.asr))),4],3)
                
                GGWGU03_analysum$vAGID <- round(GGWGU03_vc[grepl("AGID",rownames(GGWGU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU03_vc)),1],3)
                    
                GGWGU03_analysum$vAGIDxY <- round(GGWGU03_vc[grepl("AGID",rownames(GGWGU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU03_vc)),1],3)
                
                GGWGU03_analysum$vDGID <- round(GGWGU03_vc[grepl("DGID",rownames(GGWGU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU03_vc)),1],3)
                    
                GGWGU03_analysum$vDGIDxY <- round(GGWGU03_vc[grepl("DGID",rownames(GGWGU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU03_vc)),1],3)
                
                GGWGU03_analysum$vU <- round(GGWGU03_vc[grepl("U",rownames(GGWGU03_vc),fixed = T),1],3)
                
                GGWGU03_analysum$vR <- round(GGWGU03_vc[grepl("R",rownames(GGWGU03_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,GGWGU03_analysum)
                
                #<<<<<<<<<
                # 3.2.3.2.4. Test models #### 
                #<<<<<<<<<
                
                GGWGU03_tmd_in <- NULL
                GGWGU03_tmd_in$analy_name <- "GU02=GU03"
                GGWGU03_tmd_in$fullmod_name <- "GGWGU02"
                GGWGU03_tmd_in$redmod_name <- "GGWGU03"
                GGWGU03_tmd_in$fullmod.asr <- GGWGU02.asr
                GGWGU03_tmd_in$redmod.asr <- GGWGU03.asr
                GGWGU03_tmd_in$test <- "boundary" #("boundary"/"other")
                
                GGWGU03_tmd_out <- tmd.f(tmd_in =  GGWGU03_tmd_in)
                GGWGU03_tmd_out
            
                GGWGU03_tmd_out$Bestmod <- 'GGWGU03'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,GGWGU03_tmd_out)

            #<<<<<<<<<<<<<<<<<<<
            # 3.2.4.2. GGWGU04
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.2.4.2.1. Fit model
                #<<<<<<<<<
            
                GGWGU04.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_G.giv) +
                                                    U, 
                                        data= apdata_G,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.2.4.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.2.4.2.2.2. Model parameters
                    #<<<<
                
                    wald(GGWGU04.asr)
                    
                    GGWGU04_vc <- summary(GGWGU04.asr)$var
                    GGWGU04_vc
                    
                #<<<<<<<<<
                # 3.2.3.2.3. Collate results
                #<<<<<<<<<
                
                GGWGU04_analysum <- data.frame("Model_ID"="GGWGU04")
                    
                GGWGU04_analysum$logl <-   GGWGU04.asr$logl
                GGWGU04_analysum$df <-  sum(GGWGU01_vc[,'bound']%in%c("U","P"))
                
                GGWGU04_analysum$pY <- round(wald(GGWGU04.asr)[grepl("Y",rownames(wald(GGWGU04.asr))),4],3)
                
                GGWGU04_analysum$vAGID <- round(GGWGU04_vc[grepl("AGID",rownames(GGWGU04_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU04_vc)),1],3)
                    
                GGWGU04_analysum$vAGIDxY <- round(GGWGU04_vc[grepl("AGID",rownames(GGWGU04_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU04_vc)),1],3)
                
                GGWGU04_analysum$vDGID <- round(GGWGU04_vc[grepl("DGID",rownames(GGWGU04_vc),fixed = T)
                                                                    & !grepl('Y',rownames(GGWGU04_vc)),1],3)
                    
                GGWGU04_analysum$vDGIDxY <- round(GGWGU04_vc[grepl("DGID",rownames(GGWGU04_vc),fixed = T)
                                                                    & grepl('Y',rownames(GGWGU04_vc)),1],3)
                
                GGWGU04_analysum$vU <- round(GGWGU04_vc[grepl("U",rownames(GGWGU04_vc),fixed = T),1],3)
                
                GGWGU04_analysum$vR <- round(GGWGU04_vc[grepl("R",rownames(GGWGU04_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,GGWGU04_analysum)
                
                #<<<<<<<<<
                # 3.2.4.2.4. Test model differences
                #<<<<<<<<<
                
                GGWGU04_tmd_in <- NULL
                GGWGU04_tmd_in$analy_name <- "GU03=GU04"
                GGWGU04_tmd_in$fullmod_name <- "GGWGU03"
                GGWGU04_tmd_in$redmod_name <- "GGWGU04"
                GGWGU04_tmd_in$fullmod.asr <- GGWGU03.asr
                GGWGU04_tmd_in$redmod.asr <- GGWGU04.asr
                GGWGU04_tmd_in$test <- "boundary" #("boundary"/"other")
                
                GGWGU04_tmd_out <- tmd.f(tmd_in =  GGWGU04_tmd_in)
                GGWGU04_tmd_out
            
                GGWGU04_tmd_out$Bestmod <- 'GGWGU03'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,GGWGU04_tmd_out)                               
                
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.3. clarKsville (K)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.3.1. Prepare data
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.1.1 Select data
            #<<<<<<<<<<<<<<<<<<<
            
            apdata_K <- apdata_prep[apdata_prep$L%in%c("K"),]
            
            GID_apdata_K <- unique(apdata_K$GID)
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.1.2. Build AGRM and DGRM
            #<<<<<<<<<<<<<<<<<<<
            
            #<<<< Additive
            AG_K.giv <- AW.giv_K_tune
                    
            #<<<< Dominance
            DG_K.giv <- DW.giv_K_tune
                    
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.1.3. Set factors
            #<<<<<<<<<<<<<<<<<<<
            
            ### Factorial
            apdata_K$Y <- factor(apdata_K$Y)
            apdata_K$AGID <- factor(apdata_K$GID,levels = colnames(AG_K.giv))
            apdata_K$DGID <- factor(apdata_K$GID,levels = colnames(DG_K.giv))
            apdata_K$U <- factor(apdata_K$U)
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.3.2. Univariate models
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.2.1. KGWKU01
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.2.1.1. Fit model
                #<<<<<<<<<
            
                options(scipen = 999)
                KGWKU01.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv) + vm(AGID,AG_K.giv):idv(Y) + 
                                                    vm(DGID,DG_K.giv) + vm(DGID,DG_K.giv):idv(Y) + 
                                                    U, 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.2.1.2. Review model fit
                #<<<<<<<<<
                
                    #<<<<
                    # 3.3.2.1.2.1. Residuals
                    #<<<<
                
                    plot(KGWKU01.asr)
                
                    apdata_res_K <- cbind(apdata_K,KGWKU01.asr$residuals)
                    colnames(apdata_res_K)[ncol(apdata_res_K)] <-  "res"
                    
                    apdata_res_K$stdres <- apdata_res_K$res/sqrt(var(apdata_res_K$res,na.rm = T))
                    hist(apdata_res_K$stdres, main = "apdata_stdres_K", xlab = "std residuals distribution", col= "sky blue")
                    
                    apdata_res_K[abs(apdata_res_K$stdres) > 3.5 & !is.na(apdata_res_K$stdres),]
                    
                    #<<<<
                    # 3.3.2.1.2.2. Model parameters
                    #<<<<
                
                    wald(KGWKU01.asr)
                    
                    KGWKU01_vc <- summary(KGWKU01.asr)$var
                    KGWKU01_vc
                    
                #<<<<<<<<<
                # 3.3.2.1.3. Collate results
                #<<<<<<<<<
                    
                KGWKU01_analysum <- data.frame("Model_ID"="KGWKU01")
                    
                KGWKU01_analysum$logl <-   KGWKU01.asr$logl
                KGWKU01_analysum$df <-  sum(KGWKU01_vc[,'bound']%in%c("U","P"))
                
                KGWKU01_analysum$pY <- round(wald(KGWKU01.asr)[grepl("Y",rownames(wald(KGWKU01.asr))),4],3)
                
                KGWKU01_vc <- summary(KGWKU01.asr)$var
                
                KGWKU01_analysum$vAGID <- round(KGWKU01_vc[grepl("AGID",rownames(KGWKU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU01_vc)),1],3)
                    
                KGWKU01_analysum$vAGIDxY <- round(KGWKU01_vc[grepl("AGID",rownames(KGWKU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU01_vc)),1],3)
                
                KGWKU01_analysum$vDGID <- round(KGWKU01_vc[grepl("DGID",rownames(KGWKU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU01_vc)),1],3)
                    
                KGWKU01_analysum$vDGIDxY <- round(KGWKU01_vc[grepl("DGID",rownames(KGWKU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU01_vc)),1],3)
                
                KGWKU01_analysum$vU <- round(KGWKU01_vc[grepl("U",rownames(KGWKU01_vc),fixed = T),1],3)
                
                KGWKU01_analysum$vR <- round(KGWKU01_vc[grepl("R",rownames(KGWKU01_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,KGWKU01_analysum)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.2.2. KGWKU02
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.2.2.1. Fit model
                #<<<<<<<<<
            
                KGWKU02.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv) + vm(AGID,AG_K.giv):idv(Y) + 
                                                    vm(DGID,DG_K.giv) + 
                                                    U, 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.2.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.3.2.2.2.2. Model parameters
                    #<<<<
                
                    wald(KGWKU02.asr)
                    
                    KGWKU02_vc <- summary(KGWKU02.asr)$var
                    KGWKU02_vc
                    
                #<<<<<<<<<
                # 3.3.2.2.3. Collate results
                #<<<<<<<<<
                
                KGWKU02_analysum <- data.frame("Model_ID"="KGWKU02")
                    
                KGWKU02_analysum$logl <-   KGWKU02.asr$logl
                KGWKU02_analysum$df <-  sum(KGWKU01_vc[,'bound']%in%c("U","P"))
                
                KGWKU02_analysum$pY <- round(wald(KGWKU02.asr)[grepl("Y",rownames(wald(KGWKU02.asr))),4],3)
                
                KGWKU02_analysum$vAGID <- round(KGWKU02_vc[grepl("AGID",rownames(KGWKU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU02_vc)),1],3)
                    
                KGWKU02_analysum$vAGIDxY <- round(KGWKU02_vc[grepl("AGID",rownames(KGWKU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU02_vc)),1],3)
                
                KGWKU02_analysum$vDGID <- round(KGWKU02_vc[grepl("DGID",rownames(KGWKU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU02_vc)),1],3)
                    
                KGWKU02_analysum$vDGIDxY <- round(KGWKU02_vc[grepl("DGID",rownames(KGWKU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU02_vc)),1],3)
                
                KGWKU02_analysum$vU <- round(KGWKU02_vc[grepl("U",rownames(KGWKU02_vc),fixed = T),1],3)
                
                KGWKU02_analysum$vR <- round(KGWKU02_vc[grepl("R",rownames(KGWKU02_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,KGWKU02_analysum)
                
                #<<<<<<<<<
                # 3.3.2.2.4.  Test model differences 
                #<<<<<<<<<
                
                KGWKU02_tmd_in <- NULL
                KGWKU02_tmd_in$analy_name <- "KU01=KU02"
                KGWKU02_tmd_in$fullmod_name <- "KGWKU01"
                KGWKU02_tmd_in$redmod_name <- "KGWKU02"
                KGWKU02_tmd_in$fullmod.asr <- KGWKU01.asr
                KGWKU02_tmd_in$redmod.asr <- KGWKU02.asr
                KGWKU02_tmd_in$test <- "boundary" #("boundary"/"other")
                
                KGWKU02_tmd_out <- tmd.f(tmd_in =  KGWKU02_tmd_in)
                KGWKU02_tmd_out
            
                KGWKU02_tmd_out$Bestmod <- 'KGWKU02'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,KGWKU02_tmd_out)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.3.3. KGWKU03
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.3.3.1. Fit model
                #<<<<<<<<<
            
                KGWKU03.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv) +  
                                                    vm(DGID,DG_K.giv) + 
                                                    U, 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.3.3.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.3.3.3.2.2. Model parameters
                    #<<<<
                
                    wald(KGWKU03.asr)
                    
                    KGWKU03_vc <- summary(KGWKU03.asr)$var
                    KGWKU03_vc
                    
                #<<<<<<<<<
                # 3.3.3.3.3. Collate results
                #<<<<<<<<<
                
                KGWKU03_analysum <- data.frame("Model_ID"="KGWKU03")
                    
                KGWKU03_analysum$logl <-   KGWKU03.asr$logl
                KGWKU03_analysum$df <-  sum(KGWKU01_vc[,'bound']%in%c("U","P"))
                
                KGWKU03_analysum$pY <- round(wald(KGWKU03.asr)[grepl("Y",rownames(wald(KGWKU03.asr))),4],3)
                
                KGWKU03_analysum$vAGID <- round(KGWKU03_vc[grepl("AGID",rownames(KGWKU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU03_vc)),1],3)
                    
                KGWKU03_analysum$vAGIDxY <- round(KGWKU03_vc[grepl("AGID",rownames(KGWKU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU03_vc)),1],3)
                
                KGWKU03_analysum$vDGID <- round(KGWKU03_vc[grepl("DGID",rownames(KGWKU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU03_vc)),1],3)
                    
                KGWKU03_analysum$vDGIDxY <- round(KGWKU03_vc[grepl("DGID",rownames(KGWKU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU03_vc)),1],3)
                
                KGWKU03_analysum$vU <- round(KGWKU03_vc[grepl("U",rownames(KGWKU03_vc),fixed = T),1],3)
                
                KGWKU03_analysum$vR <- round(KGWKU03_vc[grepl("R",rownames(KGWKU03_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,KGWKU03_analysum)
                
                #<<<<<<<<<
                # 3.3.3.2.4.  Test model differences #### 
                #<<<<<<<<<
                
                KGWKU03_tmd_in <- NULL
                KGWKU03_tmd_in$analy_name <- "KU02=KU03"
                KGWKU03_tmd_in$fullmod_name <- "KGWKU02"
                KGWKU03_tmd_in$redmod_name <- "KGWKU03"
                KGWKU03_tmd_in$fullmod.asr <- KGWKU02.asr
                KGWKU03_tmd_in$redmod.asr <- KGWKU03.asr
                KGWKU03_tmd_in$test <- "boundary" #("boundary"/"other")
                
                KGWKU03_tmd_out <- tmd.f(tmd_in =  KGWKU03_tmd_in)
                KGWKU03_tmd_out
            
                KGWKU03_tmd_out$Bestmod <- 'KGWKU02'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,KGWKU03_tmd_out)

            #<<<<<<<<<<<<<<<<<<<
            # 3.3.4.2. KGWKU05
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.4.2.1. Fit model
                #<<<<<<<<<
            
                KGWKU05.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv) + vm(AGID,AG_K.giv):idv(Y) +
                                                    U, 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.4.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.3.4.2.2.2. Model parameters
                    #<<<<
                
                    wald(KGWKU05.asr)
                    
                    KGWKU05_vc <- summary(KGWKU05.asr)$var
                    KGWKU05_vc
                    
                #<<<<<<<<<
                # 3.3.3.3.3. Collate results
                #<<<<<<<<<
                
                KGWKU05_analysum <- data.frame("Model_ID"="KGWKU05")
                    
                KGWKU05_analysum$logl <-   KGWKU05.asr$logl
                KGWKU05_analysum$df <-  sum(KGWKU01_vc[,'bound']%in%c("U","P"))
                
                KGWKU05_analysum$pY <- round(wald(KGWKU05.asr)[grepl("Y",rownames(wald(KGWKU05.asr))),4],3)
                
                KGWKU05_analysum$vAGID <- round(KGWKU05_vc[grepl("AGID",rownames(KGWKU05_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU05_vc)),1],3)
                    
                KGWKU05_analysum$vAGIDxY <- round(KGWKU05_vc[grepl("AGID",rownames(KGWKU05_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU05_vc)),1],3)
                
                KGWKU05_analysum$vDGID <- round(KGWKU05_vc[grepl("DGID",rownames(KGWKU05_vc),fixed = T)
                                                                    & !grepl('Y',rownames(KGWKU05_vc)),1],3)
                    
                KGWKU05_analysum$vDGIDxY <- round(KGWKU05_vc[grepl("DGID",rownames(KGWKU05_vc),fixed = T)
                                                                    & grepl('Y',rownames(KGWKU05_vc)),1],3)
                
                KGWKU05_analysum$vU <- round(KGWKU05_vc[grepl("U",rownames(KGWKU05_vc),fixed = T),1],3)
                
                KGWKU05_analysum$vR <- round(KGWKU05_vc[grepl("R",rownames(KGWKU05_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,KGWKU05_analysum)
                
                #<<<<<<<<<
                # 3.3.4.2.4. Test model differences ####
                #<<<<<<<<<
                
                KGWKU05_tmd_in <- NULL
                KGWKU05_tmd_in$analy_name <- "KU02=KU05"
                KGWKU05_tmd_in$fullmod_name <- "KGWKU02"
                KGWKU05_tmd_in$redmod_name <- "KGWKU05"
                KGWKU05_tmd_in$fullmod.asr <- KGWKU02.asr
                KGWKU05_tmd_in$redmod.asr <- KGWKU05.asr
                KGWKU05_tmd_in$test <- "boundary" #("boundary"/"other")
                
                KGWKU05_tmd_out <- tmd.f(tmd_in =  KGWKU05_tmd_in)
                KGWKU05_tmd_out
            
                KGWKU05_tmd_out$Bestmod <- 'KGWKU05'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,KGWKU05_tmd_out)                
                
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.3.3. Multivariate models ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        unique(apdata_K$Y)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.3.1. KGWKM02 ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.3.1.1. Fit model
                #<<<<<<<<<
                
                KGWKM02.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv):fa(Y,1) + vm(DGID,DG_K.giv) +
                                                    idv(U), 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.3.1.2. Review model
                #<<<<<<<<<
                    
                KGWKM02.asr$loglik
                
                summary(KGWKM02.asr)$var
                
                sum(summary(KGWKM02.asr)$var$bound%in%c("P","U"))
                
                KGWKM02.Mvcov_AG_F2C_in <- NULL
                KGWKM02.Mvcov_AG_F2C_in$vparameters <- KGWKM02.asr$vparameters
                KGWKM02.Mvcov_AG_F2C_in$fac <- "AGID"
                KGWKM02.Mvcov_AG_F2C_in$E <- "Y"
                
                KGWKM02.Mvcov_AG <- fatocov(fatocov_in = KGWKM02.Mvcov_AG_F2C_in) 
                KGWKM02.Mvcov_AG
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.3.1a. KGWKM02aa ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.3.1a.1. Fit model
                #<<<<<<<<<
                
                KGWKM02a.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv):corgh(Y) + vm(DGID,DG_K.giv) +
                                                    idv(U), 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.3.1a.2. Review model
                #<<<<<<<<<
                    
                KGWKM02a.asr$loglik
                
                summary(KGWKM02a.asr)$var
                
                sum(summary(KGWKM02a.asr)$var$bound%in%c("P","U"))
                
                KGWKM02a.Mvcov_AG_F2C_in <- NULL
                KGWKM02a.Mvcov_AG_F2C_in$vparameters <- KGWKM02a.asr$vparameters
                KGWKM02a.Mvcov_AG_F2C_in$fac <- "AGID"
                KGWKM02a.Mvcov_AG_F2C_in$E <- "Y"
                
                KGWKM02a.Mvcov_AG <- fatocov(fatocov_in = KGWKM02a.Mvcov_AG_F2C_in) 
                KGWKM02a.Mvcov_AG
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.3.3.2. KGWKM06
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.3.3.2.1. Prepare data
                #<<<<<<<<<
                
                apdata_K$KAE1 <- as.character(apdata_K$Y)

                apdata_K$KAE1[apdata_K$Y == 0] <- 'K01'
                apdata_K$KAE1[apdata_K$Y == 1] <- 'K01'

                unique(apdata_K$KAE1)

                apdata_K$KAE1 <- factor(apdata_K$KAE1)
                
                #<<<<<<<<<
                # 3.3.3.2.2. Fit model
                #<<<<<<<<<
                
                KGWKM06.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_K.giv):corgh(KAE1) + vm(DGID,DG_K.giv) + 
                                                    idv(U), 
                                        data= apdata_K,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.3.3.2.3. Review model
                #<<<<<<<<<
                    
                summary(KGWKM06.asr)$var
                
                #<<<<<<<<<
                # 3.3.3.2.4. Test model differences ####
                #<<<<<<<<<
                
                KGWKM06_tmd_in <- NULL
                KGWKM06_tmd_in$analy_name <- "KM05=KM06"
                KGWKM06_tmd_in$fullmod_name <- "KGWKM05"
                KGWKM06_tmd_in$redmod_name <- "KGWKM06"
                KGWKM06_tmd_in$fullmod.asr <- KGWKM05.asr
                KGWKM06_tmd_in$redmod.asr <- KGWKM06.asr
                KGWKM06_tmd_in$test <- "other" #("boundary"/"other")
                
                KGWKM06_tmd_out <- tmd.f(tmd_in =  KGWKM06_tmd_in)
                KGWKM06_tmd_out
            
                KGWKM06_tmd_out$Bestmod <- 'KGWKM06'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,KGWKM06_tmd_out)  
                
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.4. Seneca (S)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.4.1. Prepare data
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.1.1 Select data
            #<<<<<<<<<<<<<<<<<<<
            
            apdata_S <- apdata_prep[apdata_prep$L%in%c("S"),]
            
            GID_apdata_S <- unique(apdata_S$GID)
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.1.2. Build AGRM and DGRM
            #<<<<<<<<<<<<<<<<<<<
            
            #<<<< Additive
            AG_S.giv <- AW.giv_S_tune
                    
            #<<<< Dominance
            DG_S.giv <- DW.giv_S_tune
                    
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.1.3. Set factors
            #<<<<<<<<<<<<<<<<<<<
            
            ### Factorial
            apdata_S$Y <- factor(apdata_S$Y)
            apdata_S$AGID <- factor(apdata_S$GID,levels = colnames(AG_S.giv))
            apdata_S$DGID <- factor(apdata_S$GID,levels = colnames(DG_S.giv))
            apdata_S$U <- factor(apdata_S$U)
            
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.4.2. Univariate models
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.2.1. SGWSU01
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.2.1.1. Fit model
                #<<<<<<<<<
            
                SGWSU01.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv) + vm(AGID,AG_S.giv):idv(Y) + 
                                                    vm(DGID,DG_S.giv) + vm(DGID,DG_S.giv):idv(Y) + 
                                                    U, 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.2.1.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.4.2.1.2.1. Residuals
                    #<<<<
                
                    plot(SGWSU01.asr)
                
                    apdata_res_S <- cbind(apdata_S,SGWSU01.asr$residuals)
                    colnames(apdata_res_S)[ncol(apdata_res_S)] <-  "res"
                    
                    apdata_res_S$stdres <- apdata_res_S$res/sqrt(var(apdata_res_S$res,na.rm = T))
                    hist(apdata_res_S$stdres, main = "apdata_stdres_S", xlab = "std residuals distribution", col= "sky blue")
                    
                    apdata_res_S[abs(apdata_res_S$stdres) > 3.5 & !is.na(apdata_res_S$stdres),]
                    
                    #<<<<
                    # 3.4.2.1.2.2. Model parameters
                    #<<<<
                
                    wald(SGWSU01.asr)
                    
                    SGWSU01_vc <- summary(SGWSU01.asr)$var
                    SGWSU01_vc
                    
                #<<<<<<<<<
                # 3.4.2.1.3. Collate results
                #<<<<<<<<<
                    
                SGWSU01_analysum <- data.frame("Model_ID"="SGWSU01")
                    
                SGWSU01_analysum$logl <-   SGWSU01.asr$logl
                SGWSU01_analysum$df <-  sum(SGWSU01_vc[,'bound']%in%c("U","P"))
                
                SGWSU01_analysum$pY <- round(wald(SGWSU01.asr)[grepl("Y",rownames(wald(SGWSU01.asr))),4],3)
                
                SGWSU01_vc <- summary(SGWSU01.asr)$var
                
                SGWSU01_analysum$vAGID <- round(SGWSU01_vc[grepl("AGID",rownames(SGWSU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU01_vc)),1],3)
                    
                SGWSU01_analysum$vAGIDxY <- round(SGWSU01_vc[grepl("AGID",rownames(SGWSU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU01_vc)),1],3)
                
                SGWSU01_analysum$vDGID <- round(SGWSU01_vc[grepl("DGID",rownames(SGWSU01_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU01_vc)),1],3)
                    
                SGWSU01_analysum$vDGIDxY <- round(SGWSU01_vc[grepl("DGID",rownames(SGWSU01_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU01_vc)),1],3)
                
                SGWSU01_analysum$vU <- round(SGWSU01_vc[grepl("U",rownames(SGWSU01_vc),fixed = T),1],3)
                
                SGWSU01_analysum$vR <- round(SGWSU01_vc[grepl("R",rownames(SGWSU01_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,SGWSU01_analysum)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.2.2. SGWSU02
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.2.2.1. Fit model
                #<<<<<<<<<
            
                SGWSU02.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv) + vm(AGID,AG_S.giv):idv(Y) + 
                                                    vm(DGID,DG_S.giv) + 
                                                    U, 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.2.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.4.2.2.2.2. Model parameters
                    #<<<<
                
                    wald(SGWSU02.asr)
                    
                    SGWSU02_vc <- summary(SGWSU02.asr)$var
                    SGWSU02_vc
                    
                #<<<<<<<<<
                # 3.4.2.2.3. Collate results
                #<<<<<<<<<
                
                SGWSU02_analysum <- data.frame("Model_ID"="SGWSU02")
                    
                SGWSU02_analysum$logl <-   SGWSU02.asr$logl
                SGWSU02_analysum$df <-  sum(SGWSU02_vc[,'bound']%in%c("U","P"))
                
                SGWSU02_analysum$pY <- round(wald(SGWSU02.asr)[grepl("Y",rownames(wald(SGWSU02.asr))),4],3)
                
                SGWSU02_analysum$vAGID <- round(SGWSU02_vc[grepl("AGID",rownames(SGWSU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU02_vc)),1],3)
                    
                SGWSU02_analysum$vAGIDxY <- round(SGWSU02_vc[grepl("AGID",rownames(SGWSU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU02_vc)),1],3)
                
                SGWSU02_analysum$vDGID <- round(SGWSU02_vc[grepl("DGID",rownames(SGWSU02_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU02_vc)),1],3)
                    
                SGWSU02_analysum$vDGIDxY <- round(SGWSU02_vc[grepl("DGID",rownames(SGWSU02_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU02_vc)),1],3)
                
                SGWSU02_analysum$vU <- round(SGWSU02_vc[grepl("U",rownames(SGWSU02_vc),fixed = T),1],3)
                
                SGWSU02_analysum$vR <- round(SGWSU02_vc[grepl("R",rownames(SGWSU02_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,SGWSU02_analysum)
                
                #<<<<<<<<<
                # 3.4.2.2.4.  Test model differences 
                #<<<<<<<<<
                
                SGWSU02_tmd_in <- NULL
                SGWSU02_tmd_in$analy_name <- "SU01=SU02"
                SGWSU02_tmd_in$fullmod_name <- "SGWSU01"
                SGWSU02_tmd_in$redmod_name <- "SGWSU02"
                SGWSU02_tmd_in$fullmod.asr <- SGWSU01.asr
                SGWSU02_tmd_in$redmod.asr <- SGWSU02.asr
                SGWSU02_tmd_in$test <- "boundary" #("boundary"/"other")
                
                SGWSU02_tmd_out <- tmd.f(tmd_in =  SGWSU02_tmd_in)
                SGWSU02_tmd_out
            
                SGWSU02_tmd_out$Bestmod <- 'SGWSU01'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,SGWSU02_tmd_out)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.2.4. SGWSU03 ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.2.4.1. Fit model
                #<<<<<<<<<
            
                SGWSU03.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv) +  
                                                    vm(DGID,DG_S.giv) + 
                                                    U, 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.2.4.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.4.2.4.2.2. Model parameters
                    #<<<<
                
                    wald(SGWSU03.asr)
                    
                    SGWSU03_vc <- summary(SGWSU03.asr)$var
                    SGWSU03_vc
                    
                #<<<<<<<<<
                # 3.4.2.4.3. Collate results
                #<<<<<<<<<
                
                SGWSU03_analysum <- data.frame("Model_ID"="SGWSU03")
                    
                SGWSU03_analysum$logl <-   SGWSU03.asr$logl
                SGWSU03_analysum$df <-  sum(SGWSU03_vc[,'bound']%in%c("U","P"))
                
                SGWSU03_analysum$pY <- round(wald(SGWSU03.asr)[grepl("Y",rownames(wald(SGWSU03.asr))),4],3)
                
                SGWSU03_analysum$vAGID <- round(SGWSU03_vc[grepl("AGID",rownames(SGWSU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU03_vc)),1],3)
                    
                SGWSU03_analysum$vAGIDxY <- round(SGWSU03_vc[grepl("AGID",rownames(SGWSU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU03_vc)),1],3)
                
                SGWSU03_analysum$vDGID <- round(SGWSU03_vc[grepl("DGID",rownames(SGWSU03_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU03_vc)),1],3)
                    
                SGWSU03_analysum$vDGIDxY <- round(SGWSU03_vc[grepl("DGID",rownames(SGWSU03_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU03_vc)),1],3)
                
                SGWSU03_analysum$vU <- round(SGWSU03_vc[grepl("U",rownames(SGWSU03_vc),fixed = T),1],3)
                
                SGWSU03_analysum$vR <- round(SGWSU03_vc[grepl("R",rownames(SGWSU03_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,SGWSU03_analysum)
                
                #<<<<<<<<<
                # 3.4.2.2.4.  Test model differences #### 
                #<<<<<<<<<
                
                SGWSU03_tmd_in <- NULL
                SGWSU03_tmd_in$analy_name <- "SU02=SU03"
                SGWSU03_tmd_in$fullmod_name <- "SGWSU02"
                SGWSU03_tmd_in$redmod_name <- "SGWSU03"
                SGWSU03_tmd_in$fullmod.asr <- SGWSU02.asr
                SGWSU03_tmd_in$redmod.asr <- SGWSU03.asr
                SGWSU03_tmd_in$test <- "boundary" #("boundary"/"other")
                
                SGWSU03_tmd_out <- tmd.f(tmd_in =  SGWSU03_tmd_in)
                SGWSU03_tmd_out
            
                SGWSU03_tmd_out$Bestmod <- 'SGWSU02'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,SGWSU03_tmd_out)

            #<<<<<<<<<<<<<<<<<<<
            # 3.4.2.2. SGWSU05
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.2.2.1. Fit model
                #<<<<<<<<<
            
                SGWSU05.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv) + vm(AGID,AG_S.giv):idv(Y) +
                                                    U, 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.4.2.2. Review model fit
                #<<<<<<<<<
                
                options(scipen = 999)
                
                    #<<<<
                    # 3.4.4.2.2.2. Model parameters
                    #<<<<
                
                    wald(SGWSU05.asr)
                    
                    SGWSU05_vc <- summary(SGWSU05.asr)$var
                    SGWSU05_vc
                    
                #<<<<<<<<<
                # 3.4.2.4.3. Collate results
                #<<<<<<<<<
                
                SGWSU05_analysum <- data.frame("Model_ID"="SGWSU05")
                    
                SGWSU05_analysum$logl <-   SGWSU05.asr$logl
                SGWSU05_analysum$df <-  sum(SGWSU05_vc[,'bound']%in%c("U","P"))
                
                SGWSU05_analysum$pY <- round(wald(SGWSU05.asr)[grepl("Y",rownames(wald(SGWSU05.asr))),4],3)
                
                SGWSU05_analysum$vAGID <- round(SGWSU05_vc[grepl("AGID",rownames(SGWSU05_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU05_vc)),1],3)
                    
                SGWSU05_analysum$vAGIDxY <- round(SGWSU05_vc[grepl("AGID",rownames(SGWSU05_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU05_vc)),1],3)
                
                SGWSU05_analysum$vDGID <- round(SGWSU05_vc[grepl("DGID",rownames(SGWSU05_vc),fixed = T)
                                                                    & !grepl('Y',rownames(SGWSU05_vc)),1],3)
                    
                SGWSU05_analysum$vDGIDxY <- round(SGWSU05_vc[grepl("DGID",rownames(SGWSU05_vc),fixed = T)
                                                                    & grepl('Y',rownames(SGWSU05_vc)),1],3)
                
                SGWSU05_analysum$vU <- round(SGWSU05_vc[grepl("U",rownames(SGWSU05_vc),fixed = T),1],3)
                
                SGWSU05_analysum$vR <- round(SGWSU05_vc[grepl("R",rownames(SGWSU05_vc),fixed = T),1],3)
                
                IGWIPU_analysum <- smartbind(IGWIPU_analysum,SGWSU05_analysum)
                
                #<<<<<<<<<
                # 3.4.4.2.4. Test model differences ####
                #<<<<<<<<<
                
                SGWSU05_tmd_in <- NULL
                SGWSU05_tmd_in$analy_name <- "SU02=SU05"
                SGWSU05_tmd_in$fullmod_name <- "SGWSU02"
                SGWSU05_tmd_in$redmod_name <- "SGWSU05"
                SGWSU05_tmd_in$fullmod.asr <- SGWSU02.asr
                SGWSU05_tmd_in$redmod.asr <- SGWSU05.asr
                SGWSU05_tmd_in$test <- "boundary" #("boundary"/"other")
                
                SGWSU05_tmd_out <- tmd.f(tmd_in =  SGWSU05_tmd_in)
                SGWSU05_tmd_out
            
                SGWSU05_tmd_out$Bestmod <- 'SGWSU02'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,SGWSU05_tmd_out)                
                
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 3.4.3. Multivariate models ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
        unique(apdata_S$Y)
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.3.1. SGWSM02 ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.3.1.1. Fit model
                #<<<<<<<<<
                
                SGWSM02.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv):fa(Y,1) + vm(DGID,DG_S.giv) + 
                                                    vm(DGID,DG_S.giv) +
                                                    idv(U), 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.3.1.2. Review model
                #<<<<<<<<<
                    
                SGWSM02.asr$loglik
                
                SGWSM02_vcomp <- summary(SGWSM02.asr)$var
                
                sum(summary(SGWSM02.asr)$var[,'bound']%in%c("U","P"))
                
                SGWSM02.Mvcov_AG_F2C_in <- NULL
                SGWSM02.Mvcov_AG_F2C_in$vparameters <- SGWSM02_vcomp[,1]
                names(SGWSM02.Mvcov_AG_F2C_in$vparameters) <- rownames(SGWSM02_vcomp)
                SGWSM02.Mvcov_AG_F2C_in$fac <- "AGID"
                SGWSM02.Mvcov_AG_F2C_in$E <- "Y"
                
                SGWSM02.Mvcov_AG <- fatocov(fatocov_in = SGWSM02.Mvcov_AG_F2C_in) 
                SGWSM02.Mvcov_AG
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.2.1a. SGWSM02a ####
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.2.1a.1. Fit model
                #<<<<<<<<<
                
                SGWSM02a.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv):corgh(Y) + vm(DGID,DG_S.giv) + 
                                                    vm(DGID,DG_S.giv) +
                                                    idv(U), 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.2.1a.2. Review model
                #<<<<<<<<<
                    
                SGWSM02a.asr$loglik
                
                summary(SGWSM02a.asr)$var
                
                sum(summary(SGWSM02a.asr)$var[,'bound']%in%c("U","P"))
                
            #<<<<<<<<<<<<<<<<<<<
            # 3.4.3.2. SGWSM06
            #<<<<<<<<<<<<<<<<<<<
                
                #<<<<<<<<<
                # 3.4.3.2.1. Prepare data
                #<<<<<<<<<
                
                apdata_S$SAE1 <- as.character(apdata_S$Y)

                apdata_S$SAE1[apdata_S$Y == 0] <- 'S01'
                apdata_S$SAE1[apdata_S$Y == 1] <- 'S01'

                unique(apdata_S$SAE1)

                apdata_S$SAE1 <- factor(apdata_S$SAE1)
                
                #<<<<<<<<<
                # 3.4.3.2.2. Fit model
                #<<<<<<<<<
                
                SGWSM06.asr <- asreml(snewy ~ 1 + Y,
                                        random = ~ vm(AGID,AG_S.giv):corgh(SAE1) + vm(DGID,DG_S.giv) +
                                                    vm(DGID,DG_S.giv) +
                                                    idv(U), 
                                        data= apdata_S,
                                        na.action = na.method(y = "include", x = "include"),
                                        workspace = 6e+08,
                                        maxiter=50)
                
                #<<<<<<<<<
                # 3.4.3.2.3. Review model
                #<<<<<<<<<
                    
                summary(SGWSM06.asr)$var
                
                #<<<<<<<<<
                # 3.4.3.2.4. Test model differences ####
                #<<<<<<<<<
                
                SGWSM06_tmd_in <- NULL
                SGWSM06_tmd_in$analy_name <- "SM02=SM06"
                SGWSM06_tmd_in$fullmod_name <- "SGWSM02"
                SGWSM06_tmd_in$redmod_name <- "SGWSM06"
                SGWSM06_tmd_in$fullmod.asr <- SGWSM02.asr
                SGWSM06_tmd_in$redmod.asr <- SGWSM06.asr
                SGWSM06_tmd_in$test <- "other" #("boundary"/"other")
                
                SGWSM06_tmd_out <- tmd.f(tmd_in =  SGWSM06_tmd_in)
                SGWSM06_tmd_out
            
                SGWSM06_tmd_out$Bestmod <- 'SGWSM02'
                
                IGWIPU_tmd_sum <- smartbind(IGWIPU_tmd_sum,SGWSM06_tmd_out)        
                
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Write out summaries
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                
write.csv(IGWIPU_tmd_sum,"IGWIPU_tmd_sum.csv")
write.csv(IGWIPU_analysum,"IGWIPU_analysum.csv")
