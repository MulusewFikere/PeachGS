#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
##### Cross validation of models ####
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
    library(ASRgenomics)
    
    options(scipen = 999)
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 1.2. functions
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. prepare data ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # # 2.1. Load curated data ####
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 
    # load("../../../4. Fitting linear models/4.4. XQBX/XQBX_210918.RData")
    # 
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # # 2.2. Set up validation populations in data ####
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 
    # nCV_vpops = 5
    # 
    #     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #     # 2.2.1. across trial validation, i.e. exclude accession from all data  ####
    #     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 
    #     XTCV_vpop_list <- as.data.frame(unique(apdata$GID),stringsAsFactors = F)
    #     colnames(XTCV_vpop_list) <- 'GID'
    # 
    #     XTCV_vpop_list$XTCV_vpop <- ceiling(sample(1:nrow(XTCV_vpop_list))/(nrow(XTCV_vpop_list)/nCV_vpops))
    # 
    #     apdata_XTCV_temp <- merge(apdata,XTCV_vpop_list,all.x = T)
    # 
    #     apdata_XTCV_nna <- apdata_XTCV_temp[!is.na(apdata_XTCV_temp$snewy),]
    # 
    #     apdata_XTCV_final <- apdata_XTCV_nna
    # 
    #     table(apdata_XTCV_nna$XTCV_vpop,apdata_XTCV_nna$LY)
    # 
    #     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #     # 2.2.2. Within trial validation, i.e. exclude accession from individual trial data  ####
    #     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 
    #         #<<<<<<<<<<<<<<<<<<<
    #         # 2.2.2.1. F  ####
    #         #<<<<<<<<<<<<<<<<<<<
    # 
    #         apdata_XTCV_F <- apdata_XTCV_nna[apdata_XTCV_nna$L == 'F',]
    # 
    #         WTCV_tvpop_list_F <- as.data.frame(unique(apdata_XTCV_F$GID),stringsAsFactors = F)
    #         colnames(WTCV_tvpop_list_F) <- 'GID'
    # 
    #         WTCV_tvpop_list_F$WTCV_tvpop <- ceiling(sample(1:nrow(WTCV_tvpop_list_F))/(nrow(WTCV_tvpop_list_F)/nCV_vpops))
    # 
    #         apdata_WTCV_F <- merge(apdata_XTCV_F,WTCV_tvpop_list_F,all.x = T)
    # 
    #         #<<<<<<<<<<<<<<<<<<<
    #         # 2.2.2.2. G ####
    #         #<<<<<<<<<<<<<<<<<<<
    # 
    #         apdata_XTCV_G <- apdata_XTCV_nna[apdata_XTCV_nna$L == 'G',]
    # 
    #         WTCV_tvpop_list_G <- as.data.frame(unique(apdata_XTCV_G$GID),stringsAsFactors = F)
    #         colnames(WTCV_tvpop_list_G) <- 'GID'
    # 
    #         WTCV_tvpop_list_G$WTCV_tvpop <- ceiling(sample(1:nrow(WTCV_tvpop_list_G))/(nrow(WTCV_tvpop_list_G)/nCV_vpops))
    # 
    #         apdata_WTCV_G <- merge(apdata_XTCV_G,WTCV_tvpop_list_G,all.x = T)
    # 
    #         #<<<<<<<<<<<<<<<<<<<
    #         # 2.2.2.3. K ####
    #         #<<<<<<<<<<<<<<<<<<<
    # 
    #         apdata_XTCV_K <- apdata_XTCV_nna[apdata_XTCV_nna$L == 'K',]
    # 
    #         WTCV_tvpop_list_K <- as.data.frame(unique(apdata_XTCV_K$GID),stringsAsFactors = F)
    #         colnames(WTCV_tvpop_list_K) <- 'GID'
    # 
    #         WTCV_tvpop_list_K$WTCV_tvpop <- ceiling(sample(1:nrow(WTCV_tvpop_list_K))/(nrow(WTCV_tvpop_list_K)/nCV_vpops))
    # 
    #         apdata_WTCV_K <- merge(apdata_XTCV_K,WTCV_tvpop_list_K,all.x = T)
    # 
    #         #<<<<<<<<<<<<<<<<<<<
    #         # 2.2.2.4. S ####
    #         #<<<<<<<<<<<<<<<<<<<
    # 
    #         apdata_XTCV_S <- apdata_XTCV_nna[apdata_XTCV_nna$L == 'S',]
    # 
    #         WTCV_tvpop_list_S <- as.data.frame(unique(apdata_XTCV_S$GID),stringsAsFactors = F)
    #         colnames(WTCV_tvpop_list_S) <- 'GID'
    # 
    #         WTCV_tvpop_list_S$WTCV_tvpop <- ceiling(sample(1:nrow(WTCV_tvpop_list_S))/(nrow(WTCV_tvpop_list_S)/nCV_vpops))
    # 
    #         apdata_WTCV_S <- merge(apdata_XTCV_S,WTCV_tvpop_list_S,all.x = T)
    # 
    #         #<<<<<<<<<<<<<<<<<<<
    #         # 2.2.2.5. apdata_CV_temp ####
    #         #<<<<<<<<<<<<<<<<<<<
    # 
    #         apdata_CV_temp <- rbind(apdata_WTCV_F,apdata_WTCV_G,apdata_WTCV_K,apdata_WTCV_S)
    #         nrow(apdata_CV_temp)
    # 
    #         apdata_CV_temp_nna <- apdata_CV_temp[!is.na(apdata_CV_temp$snewy),]
    # 
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # # 2.3. Output workspace ####
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 
    # uGLX <- unique(apdata_CV_temp_nna[,c("GID","L","XTCV_vpop")])
    # table(uGLX$XTCV_vpop,uGLX$L)
    #
    # uGLW <- unique(apdata_CV_temp_nna[,c("GID","L","WTCV_tvpop")])
    # table(uGLW$WTCV_tvpop,uGLW$L)
    #         
    # save.image("apdata_CV_temp.RData")
            
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 2.4. Input workspace ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    load("apdata_CV_temp.RData")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Predict non-genomic validation phenotypes ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.1. Prepare  data ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    apdata_NG <- apdata_CV_temp_nna

    apdata_NG$fL <- factor(apdata_NG$L)
    apdata_NG$fLY <- factor(apdata_NG$LY)
    apdata_NG$fU <- factor(apdata_NG$U)
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.2. Fit non-genomic model ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    NG.asr <- asreml(snewy ~ 1 + fLY,
                            random = ~ at(fL):idv(fU),
                            residual = ~ dsum(~idv(units)| fL),
                            data = apdata_NG,
                            na.action = na.method(y = "include", x = "include"),
                            workspace = 6e+08,
                            maxiter=50)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.3. Prepare validation phenotypes ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<
        # 3.3.1. bU_NG ####
        #<<<<<<<<<<<<<<<<<<<

        bU_L_tmp <- as.data.frame(NG.asr$coefficients$random)
        colnames(bU_L_tmp) <- 'bU'

        bU_L_tmp$L <- sub("at(fL, ","",sapply(strsplit(rownames(bU_L_tmp),split = "):"),'[',1),fixed = T)
        bU_L_tmp$U  <- sub("fU_","",sapply(strsplit(rownames(bU_L_tmp),split = "):"),'[',2),fixed = T)

        bU_L <- bU_L_tmp[bU_L_tmp$L == substr(bU_L_tmp$U,1,1),]
        table(bU_L$L)

        rownames(bU_L) <- NULL

        #<<<<<<<<<<<<<<<<<<<
        # 3.3.2. bUe_NG ####
        #<<<<<<<<<<<<<<<<<<<

        LUYe_tmp <- cbind(apdata_NG[,c("L","Y","LY","U","GID","QID","fAE4","fAB1","snewy","XTCV_vpop","WTCV_vpop")],NG.asr$residuals)

        ##### Merge U and e
        apdata_CV <- merge(LUYe_tmp,bU_L)

        #<<<< Create Ue
        apdata_CV$bUe <- apdata_CV$bU + apdata_CV$e

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.3. Create new fields ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    apdata_CV$AE <- as.character(apdata_CV$fAE4)
    apdata_CV$LAE <- paste(apdata_CV$L, apdata_CV$AE,sep = "_")
    apdata_CV$AB1 <- as.character(apdata_CV$fAB1)
    apdata_CV$msnewy <- apdata_CV$snewy
    apdata_CV$WTCV_tvpop <- paste(apdata_CV$L,apdata_CV$WTCV_vpop,sep="_")

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.5. Create lists and keys
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    QID_list <- unique(apdata_CV$QID)

    QID_GID_key <- unique(apdata_CV[,c("QID","GID")])

    LY_L_key <- unique(apdata_CV[,c("LY","L")])

    LY_AE_key <- unique(apdata_CV[,c("L","Y","LY","AE")])

    L_AE_AB1_LAE_key <- unique(apdata_CV[,c("L","AE","AB1","LAE")])

    L_AE_LAE_key <- unique(apdata_CV[,c("L","AE","LAE")])
    
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # # 3.6. Alternative adjusted phenotypes
    # #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 
    # mLY <- predict(XQBMB1.asr,classify = "LY")$pvals
    # colnames(mLY)[2] <-"mLY"    
    # 
    # apdata_CV_merge <- merge(apdata_CV,mLY[,c("LY","mLY")])
    #     
    # apdata_CV_merge$snewy_mLY <- apdata_CV_merge$snewy - apdata_CV_merge$mLY   
    # 
    # plot(apdata_CV_merge$snewy_mLY,apdata_CV_merge$bUe)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Full models
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 4.1. Individual models ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.1. Fresno ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.1.1. Prepare data
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            apdata_CV_F <- apdata_CV[apdata_CV$L == "F",]

            AW.giv_F_tune <- AW.giv_F_tune
            DW.giv_F_tune <- DW.giv_F_tune

            apdata_CV_F$fL <- factor(apdata_CV_F$L)
            apdata_CV_F$fU <- factor(apdata_CV_F$U)
            apdata_CV_F$AWID <- factor(apdata_CV_F$GID, levels = rownames(AW.giv_F_tune))
            apdata_CV_F$DWID <- factor(apdata_CV_F$GID, levels = rownames(DW.giv_F_tune))
            apdata_CV_F$fY <- factor(apdata_CV_F$Y)
            apdata_CV_F$fLY <- factor(apdata_CV_F$LY)

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.1.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            FGW.asr <- asreml(msnewy ~ 1 + fLY,
                                    random = ~ vm(AWID,AW.giv_F_tune) +
                                                vm(DWID,DW.giv_F_tune) +
                                                fU,
                                    residual = ~ idv(units),
                                    data= apdata_CV_F,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.2. colleGe Station ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.2.1. Prepare data
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            apdata_CV_G <- apdata_CV[apdata_CV$L == "G",]

            AW.giv_G_tune <- AW.giv_G_tune
            DW.giv_G_tune <- DW.giv_G_tune

            apdata_CV_G$fL <- factor(apdata_CV_G$L)
            apdata_CV_G$fU <- factor(apdata_CV_G$U)
            apdata_CV_G$AWID <- factor(apdata_CV_G$GID, levels = rownames(AW.giv_G_tune))
            apdata_CV_G$DWID <- factor(apdata_CV_G$GID, levels = rownames(DW.giv_G_tune))
            apdata_CV_G$fU <- factor(apdata_CV_G$U)
            apdata_CV_G$fY <- factor(apdata_CV_G$Y)
            apdata_CV_G$fLY <- factor(apdata_CV_G$LY)

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.2.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            GGW.asr <- asreml(msnewy ~ 1 + fLY,
                                    random = ~ vm(AWID,AW.giv_G_tune) +
                                                vm(DWID,DW.giv_G_tune) +
                                                fU,
                                    residual = ~ idv(units),
                                    data= apdata_CV_G,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.3. clarKsville ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.3.1. Prepare data
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            apdata_CV_K <- apdata_CV[apdata_CV$L == "K",]

            AW.giv_K_tune <- AW.giv_K_tune
            DW.giv_K_tune <- DW.giv_K_tune

            apdata_CV_K$fL <- factor(apdata_CV_K$L)
            apdata_CV_K$fU <- factor(apdata_CV_K$U)
            apdata_CV_K$AWID <- factor(apdata_CV_K$GID, levels = rownames(AW.giv_K_tune))
            apdata_CV_K$DWID <- factor(apdata_CV_K$GID, levels = rownames(DW.giv_K_tune))
            apdata_CV_K$fY <- factor(apdata_CV_K$Y)
            apdata_CV_K$fLY <- factor(apdata_CV_K$LY)
            apdata_CV_K$fAE <- factor(apdata_CV_K$AE)

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.3.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            KGW.asr <- asreml(msnewy ~ 1 + fLY,
                                    random = ~ vm(AWID,AW.giv_K_tune):corgh(fAE) +
                                                vm(DWID,DW.giv_K_tune) +
                                                fU,
                                    residual = ~ idv(units),
                                    data= apdata_CV_K,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.1.4. Seneca ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.4.1. Prepare data
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            apdata_CV_S <- apdata_CV[apdata_CV$L == "S",]

            AW.giv_S_tune <- AW.giv_S_tune
            DW.giv_S_tune <- DW.giv_S_tune

            apdata_CV_S$fL <- factor(apdata_CV_S$L)
            apdata_CV_S$fU <- factor(apdata_CV_S$U)
            apdata_CV_S$AWID <- factor(apdata_CV_S$GID, levels = rownames(AW.giv_S_tune))
            apdata_CV_S$DWID <- factor(apdata_CV_S$GID, levels = rownames(DW.giv_S_tune))
            apdata_CV_S$fY <- factor(apdata_CV_S$Y)
            apdata_CV_S$fLY <- factor(apdata_CV_S$LY)
            apdata_CV_S$fAE <- factor(apdata_CV_S$AE)

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.1.4.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            SGW.asr <- asreml(msnewy ~ 1 + fLY,
                                    random = ~ vm(AWID,AW.giv_S_tune):fa(fAE) +
                                                vm(DWID,DW.giv_S_tune) +
                                                fU,
                                    residual = ~ idv(units),
                                    data= apdata_CV_S,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 4.2. Multi-trial models
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.2.1. XGW ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.2.1.1. Prepare data
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            apdata_CV_X <- apdata_CV

            AW.giv_X_tune <- AW.giv_tune
            DW.giv_X_tune <- DW.giv_tune

            apdata_CV_X$fL <- factor(apdata_CV_X$L)
            apdata_CV_X$fU <- factor(apdata_CV_X$U)
            apdata_CV_X$AWID <- factor(apdata_CV_X$GID, levels = rownames(AW.giv_X_tune))
            apdata_CV_X$DWID <- factor(apdata_CV_X$GID, levels = rownames(DW.giv_X_tune))
            apdata_CV_X$fY <- factor(apdata_CV_X$Y)
            apdata_CV_X$fLY <- factor(apdata_CV_X$LY)
            apdata_CV_X$fAE <- factor(apdata_CV_X$AE)

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.2.1.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            XGW.asr <- asreml(msnewy ~ 1 + fLY,
                                    random = ~ vm(AWID,AW.giv_X_tune):fa(fAE) +
                                                vm(DWID,DW.giv_X_tune) +
                                                at(fL):idv(fU),
                                    residual = ~ dsum(~idv(units)| fL),
                                    data= apdata_CV_X,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 4.2.2. XQB ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.2.2.1. Prepare data
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            apdata_CV_X <- apdata_CV_X

            AQ.giv_QID_X_tune <- AQ.giv_QID_tune
            AB.giv_X_tune <- AB.giv_tune
            DW.giv_X_tune <- DW.giv_tune

            apdata_CV_X$fL <- factor(apdata_CV_X$L)
            apdata_CV_X$fU <- factor(apdata_CV_X$U)
            apdata_CV_X$AQID <- factor(apdata_CV_X$QID, levels = rownames(AQ.giv_QID_X_tune))
            apdata_CV_X$ABID <- factor(apdata_CV_X$GID, levels = rownames(AB.giv_X_tune))
            apdata_CV_X$DWID <- factor(apdata_CV_X$GID, levels = rownames(DW.giv_X_tune))
            apdata_CV_X$fY <- factor(apdata_CV_X$Y)
            apdata_CV_X$fLY <- factor(apdata_CV_X$LY)
            apdata_CV_X$fAE <- factor(apdata_CV_X$AE)

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 4.2.2.2. Fit model ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            XQB.asr <- asreml(msnewy ~ 1 + fLY,
                                    random = ~  vm(AQID,AQ.giv_QID_X_tune) +
                                                vm(ABID,AB.giv_X_tune):fa(fAE) +
                                                vm(DWID,DW.giv_X_tune) +
                                                at(fL):idv(fU),
                                    residual = ~ dsum(~idv(units)| fL),
                                    data= apdata_CV_X,
                                    na.action = na.method(y = "include", x = "include"),
                                    workspace = 6e+08,
                                    maxiter=50)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 5. Loop for LOO cross trial cross validation ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.1. Set up pre-loop objects ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    colnames_b_LOOCV_g_template <- c("v_sampmeth","vpop_predmeth","g","model","mod_id","Env","bAW","bAQ","bAB","bDW")
    
    b_LOOCV_g_template <- as.data.frame(matrix(NA,ncol=length(colnames_b_LOOCV_g_template),nrow=1))
    colnames(b_LOOCV_g_template) <- colnames_b_LOOCV_g_template
            
    b_LOOCV_g_template$v_sampmeth <- 'LOO'
    b_LOOCV_g_template$vpop_predmeth <- 'mask'

    b_LOOCV_out <- NULL

    g_list <- GID_list
    #g_list <- c("galaxy","elberta","a_657")
    
    n = 0
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.2. Open loop ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    for(g in g_list)
    {
        #g = g_list[1]

        print(g)
        print(n+1)
        
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.1. Set up summary object ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        b_LOOCV_g <-b_LOOCV_g_template
        b_LOOCV_g$g <- g

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.2. Set up data set ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        apdata_LOOCV_g <- apdata_CV_X
        apdata_LOOCV_g$msnewy[apdata_LOOCV_g$GID == g] <- NA

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.3. Individual trial models ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 5.2.3_F. Fresno ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            b_LOOCV_mask_FGW_g <- b_LOOCV_g
            b_LOOCV_mask_FGW_g$model <- "FGWFU03"
            b_LOOCV_mask_FGW_g$mod_id <- "STGW"

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_F.1. Data set up ####
                #<<<<<<<<<<<<<<<<<<<

                #<<<<< tpop
                apdata_LOOCV_mask_F_tpop_g <- apdata_LOOCV_g[apdata_LOOCV_g$L == "F",]
            
            if(g %in% unique(apdata_LOOCV_mask_F_tpop_g$GID))
            {
                apdata_LOOCV_mask_F_tpop_g$fU <- factor(apdata_LOOCV_mask_F_tpop_g$U)
                apdata_LOOCV_mask_F_tpop_g$AWID <- factor(apdata_LOOCV_mask_F_tpop_g$GID,levels = rownames(AW.giv_F_tune))
                apdata_LOOCV_mask_F_tpop_g$DWID <- factor(apdata_LOOCV_mask_F_tpop_g$GID,levels = rownames(DW.giv_F_tune))
                apdata_LOOCV_mask_F_tpop_g$fLY <- factor(apdata_LOOCV_mask_F_tpop_g$LY)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_F.mask.2. Fit model ####
                #<<<<<<<<<<<<<<<<<<<

                LOOCV_mask_FGW_g.asr <- update(FGW.asr,
                                                data = apdata_LOOCV_mask_F_tpop_g)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_F.mask.3. Extract AW and DW for g ####
                #<<<<<<<<<<<<<<<<<<<

                b_LOOCV_mask_FGW_g$Env <- unique(apdata_LOOCV_mask_F_tpop_g$L)
                b_LOOCV_mask_FGW_g$bAW <- LOOCV_mask_FGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_FGW_g.asr$coefficients$random)) &
                                                                    grepl("AWID",rownames(LOOCV_mask_FGW_g.asr$coefficients$random)),]
                b_LOOCV_mask_FGW_g$bDW <- LOOCV_mask_FGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_FGW_g.asr$coefficients$random)) &
                                                                    grepl("DWID",rownames(LOOCV_mask_FGW_g.asr$coefficients$random)),]
            } 
            
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 5.2.3_G. College Station ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            b_LOOCV_mask_GGW_g <- b_LOOCV_g
            b_LOOCV_mask_GGW_g$model <- "GGWGU03"
            b_LOOCV_mask_GGW_g$mod_id <- "STGW"

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_G.1. Data set up ####
                #<<<<<<<<<<<<<<<<<<<

                #<<<<< tpop
                apdata_LOOCV_mask_G_tpop_g <- apdata_LOOCV_g[apdata_LOOCV_g$L == "G",]
            
            if(g %in% unique(apdata_LOOCV_mask_G_tpop_g$GID))
            {
                apdata_LOOCV_mask_G_tpop_g$fU <- factor(apdata_LOOCV_mask_G_tpop_g$U)
                apdata_LOOCV_mask_G_tpop_g$AWID <- factor(apdata_LOOCV_mask_G_tpop_g$GID,levels = rownames(AW.giv_G_tune))
                apdata_LOOCV_mask_G_tpop_g$DWID <- factor(apdata_LOOCV_mask_G_tpop_g$GID,levels = rownames(DW.giv_G_tune))
                apdata_LOOCV_mask_G_tpop_g$fLY <- factor(apdata_LOOCV_mask_G_tpop_g$LY)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_G.mask.2. Fit model ####
                #<<<<<<<<<<<<<<<<<<<

                LOOCV_mask_GGW_g.asr <- update(GGW.asr,
                                                data = apdata_LOOCV_mask_G_tpop_g)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_G.mask.3. Extract AW and DW for g ####
                #<<<<<<<<<<<<<<<<<<<

                b_LOOCV_mask_GGW_g$Env <- unique(apdata_LOOCV_mask_G_tpop_g$L)
                b_LOOCV_mask_GGW_g$bAW <- LOOCV_mask_GGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_GGW_g.asr$coefficients$random)) &
                                                                    grepl("AWID",rownames(LOOCV_mask_GGW_g.asr$coefficients$random)),]
                b_LOOCV_mask_GGW_g$bDW <- LOOCV_mask_GGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_GGW_g.asr$coefficients$random)) &
                                                                    grepl("DWID",rownames(LOOCV_mask_GGW_g.asr$coefficients$random)),]
            }
            
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 5.2.3_K. Clarksville ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            b_LOOCV_mask_KGW_g <- b_LOOCV_g
            b_LOOCV_mask_KGW_g$model <- "KWKM06"
            b_LOOCV_mask_KGW_g$mod_id <- "STGW"

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_K.1. Data set up ####
                #<<<<<<<<<<<<<<<<<<<

                #<<<<< tpop
                apdata_LOOCV_mask_K_tpop_g <- apdata_LOOCV_g[apdata_LOOCV_g$L == "K",]
            
            if(g %in% unique(apdata_LOOCV_mask_K_tpop_g$GID))
            {
                apdata_LOOCV_mask_K_tpop_g$fU <- factor(apdata_LOOCV_mask_K_tpop_g$U)
                apdata_LOOCV_mask_K_tpop_g$AWID <- factor(apdata_LOOCV_mask_K_tpop_g$GID,levels = rownames(AW.giv_K_tune))
                apdata_LOOCV_mask_K_tpop_g$DWID <- factor(apdata_LOOCV_mask_K_tpop_g$GID,levels = rownames(DW.giv_K_tune))
                apdata_LOOCV_mask_K_tpop_g$fLY <- factor(apdata_LOOCV_mask_K_tpop_g$LY)
                apdata_LOOCV_mask_K_tpop_g$fAE <- factor(apdata_LOOCV_mask_K_tpop_g$fAE)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_K.mask.2. Fit model ####
                #<<<<<<<<<<<<<<<<<<<

                LOOCV_mask_KGW_g.asr <- update(KGW.asr,
                                                data = apdata_LOOCV_mask_K_tpop_g)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_K.mask.3. Extract AW and DW for g ####
                #<<<<<<<<<<<<<<<<<<<

                b_LOOCV_mask_KGW_g$bDW <- LOOCV_mask_KGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_KGW_g.asr$coefficients$random)) &
                                                                    grepl("DWID",rownames(LOOCV_mask_KGW_g.asr$coefficients$random)),]

                b_LOOCV_mask_KGW_g_bAW <- as.data.frame(LOOCV_mask_KGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_KGW_g.asr$coefficients$random)) &
                                                                    grepl("AWID",rownames(LOOCV_mask_KGW_g.asr$coefficients$random)),])
                colnames(b_LOOCV_mask_KGW_g_bAW) <- "bAW"
                b_LOOCV_mask_KGW_g_bAW$Env <-  sapply(strsplit(rownames(b_LOOCV_mask_KGW_g_bAW),split = "fAE_",fixed = T),"[[",2)
                
                b_LOOCV_mask_KGW_g <- merge(b_LOOCV_mask_KGW_g[,!colnames(b_LOOCV_mask_KGW_g)%in%c("Env","bAW")],b_LOOCV_mask_KGW_g_bAW, all = T) 
            } 

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 5.2.3_S. Seneca ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            b_LOOCV_mask_SGW_g <- b_LOOCV_g
            b_LOOCV_mask_SGW_g$model <- "SWSM06"
            b_LOOCV_mask_SGW_g$mod_id <- "STGW"

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_S.1. Data set up ####
                #<<<<<<<<<<<<<<<<<<<

                #<<<<< tpop
                apdata_LOOCV_mask_S_tpop_g <- apdata_LOOCV_g[apdata_LOOCV_g$L == "S",]
            
            if(g %in% unique(apdata_LOOCV_mask_S_tpop_g$GID))
            {
                apdata_LOOCV_mask_S_tpop_g$fU <- factor(apdata_LOOCV_mask_S_tpop_g$U)
                apdata_LOOCV_mask_S_tpop_g$AWID <- factor(apdata_LOOCV_mask_S_tpop_g$GID,levels = rownames(AW.giv_S_tune))
                apdata_LOOCV_mask_S_tpop_g$DWID <- factor(apdata_LOOCV_mask_S_tpop_g$GID,levels = rownames(DW.giv_S_tune))
                apdata_LOOCV_mask_S_tpop_g$fLY <- factor(apdata_LOOCV_mask_S_tpop_g$LY)
                apdata_LOOCV_mask_S_tpop_g$fAE <- factor(apdata_LOOCV_mask_S_tpop_g$fAE)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_S.mask.2. Fit model ####
                #<<<<<<<<<<<<<<<<<<<

                LOOCV_mask_SGW_g.asr <- update(SGW.asr,
                                                data = apdata_LOOCV_mask_S_tpop_g)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.3_S.mask.3. Extract AW and DW for g ####
                #<<<<<<<<<<<<<<<<<<<

                b_LOOCV_mask_SGW_g$bDW <- LOOCV_mask_SGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_SGW_g.asr$coefficients$random)) &
                                                                    grepl("DWID",rownames(LOOCV_mask_SGW_g.asr$coefficients$random)),]

                b_LOOCV_mask_SGW_g_bAW <- as.data.frame(LOOCV_mask_SGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_SGW_g.asr$coefficients$random)) &
                                                                    grepl("AWID",rownames(LOOCV_mask_SGW_g.asr$coefficients$random)) &
                                                                    !grepl("Comp1",rownames(LOOCV_mask_SGW_g.asr$coefficients$random)),])
                colnames(b_LOOCV_mask_SGW_g_bAW) <- "bAW"
                b_LOOCV_mask_SGW_g_bAW$Env <-  sapply(strsplit(rownames(b_LOOCV_mask_SGW_g_bAW),split = "fa(fAE)_",fixed = T),"[[",2)
                
                b_LOOCV_mask_SGW_g <- merge(b_LOOCV_mask_SGW_g[,!colnames(b_LOOCV_mask_SGW_g)%in%c("Env","bAW")],b_LOOCV_mask_SGW_g_bAW, all = T) 
            } 
       
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.4. Mulit-trial models ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 5.2.4._XGW ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            b_LOOCV_mask_XGW_g <- b_LOOCV_g
            b_LOOCV_mask_XGW_g$model <- "XGWM07"
            b_LOOCV_mask_XGW_g$mod_id <- "XGW"

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.4_XGW.1. Data set up ####
                #<<<<<<<<<<<<<<<<<<<

                apdata_LOOCV_mask_X_tpop_g <- apdata_LOOCV_g
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.4_XGW.mask.2. Fit model ####
                #<<<<<<<<<<<<<<<<<<<

                LOOCV_mask_XGW_g.asr <- update(XGW.asr,
                                                data = apdata_LOOCV_mask_X_tpop_g)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.4_XGW.mask.3. Extract predictions
                #<<<<<<<<<<<<<<<<<<<
            
                b_LOOCV_mask_XGW_g$bDW <- LOOCV_mask_XGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_XGW_g.asr$coefficients$random)) &
                                                                    grepl("DWID",rownames(LOOCV_mask_XGW_g.asr$coefficients$random)),]

                b_LOOCV_mask_XGW_g_bAW <- as.data.frame(LOOCV_mask_XGW_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_XGW_g.asr$coefficients$random)) &
                                                                    grepl("AWID",rownames(LOOCV_mask_XGW_g.asr$coefficients$random)) &
                                                                    !grepl("Comp1",rownames(LOOCV_mask_XGW_g.asr$coefficients$random)),])
                colnames(b_LOOCV_mask_XGW_g_bAW) <- "bAW"
                b_LOOCV_mask_XGW_g_bAW$Env <-  sapply(strsplit(rownames(b_LOOCV_mask_XGW_g_bAW),split = "fa(fAE)_",fixed = T),"[[",2)
                
                b_LOOCV_mask_XGW_g <- merge(b_LOOCV_mask_XGW_g[,!colnames(b_LOOCV_mask_XGW_g)%in%c("Env","bAW")],b_LOOCV_mask_XGW_g_bAW, all = T)
                
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # 5.2.4.2. XQB ####
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            b_LOOCV_mask_XQB_g <- b_LOOCV_g
            b_LOOCV_mask_XQB_g$model <- "XQBMB1"
            b_LOOCV_mask_XQB_g$mod_id <- "XQB"

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.4_XQB.1. Data set up ####
                #<<<<<<<<<<<<<<<<<<<

                #<<<<< tpop
                apdata_LOOCV_mask_X_tpop_g <- apdata_LOOCV_mask_X_tpop_g

                #<<<<<<<<<<<<<<<<<<<
                # 5.2.4_XQB.mask.2. Fit model ####
                #<<<<<<<<<<<<<<<<<<<

                LOOCV_mask_XQB_g.asr <- update(XQB.asr,
                                                data = apdata_LOOCV_mask_X_tpop_g)
            
                #<<<<<<<<<<<<<<<<<<<
                # 5.2.4_XGW.mask.3. Extract predictions
                #<<<<<<<<<<<<<<<<<<<
            
                b_LOOCV_mask_XQB_g$bDW <- LOOCV_mask_XQB_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_XQB_g.asr$coefficients$random)) &
                                                                    grepl("DWID",rownames(LOOCV_mask_XQB_g.asr$coefficients$random)),]

                GID_QID <- unique(apdata_LOOCV_mask_X_tpop_g[apdata_LOOCV_mask_X_tpop_g$GID == g,c("GID","QID")])
                
                b_LOOCV_mask_XQB_g_bAQ <- LOOCV_mask_XQB_g.asr$coefficients$random[grepl(GID_QID$QID,rownames(LOOCV_mask_XQB_g.asr$coefficients$random)) &
                                                                    grepl("AQID",rownames(LOOCV_mask_XQB_g.asr$coefficients$random)),]
                if(length(b_LOOCV_mask_XQB_g_bAQ) >1)
                {
                    b_LOOCV_mask_XQB_g$bAQ <- b_LOOCV_mask_XQB_g_bAQ[order(names(b_LOOCV_mask_XQB_g_bAQ))][1]
                } else
                {
                    b_LOOCV_mask_XQB_g$bAQ <- b_LOOCV_mask_XQB_g_bAQ
                }
                b_LOOCV_mask_XQB_g_bAB <- as.data.frame(LOOCV_mask_XQB_g.asr$coefficients$random[grepl(g,rownames(LOOCV_mask_XQB_g.asr$coefficients$random)) &
                                                                    grepl("ABID",rownames(LOOCV_mask_XQB_g.asr$coefficients$random)) &
                                                                    !grepl("Comp1",rownames(LOOCV_mask_XQB_g.asr$coefficients$random)),])
                colnames(b_LOOCV_mask_XQB_g_bAB) <- "bAB"
                b_LOOCV_mask_XQB_g_bAB$Env <-  sapply(strsplit(rownames(b_LOOCV_mask_XQB_g_bAB),split = "fa(fAE)_",fixed = T),"[[",2)
                
                b_LOOCV_mask_XQB_g <- merge(b_LOOCV_mask_XQB_g[,!colnames(b_LOOCV_mask_XQB_g)%in%c("Env","bAB")],b_LOOCV_mask_XQB_g_bAB, all = T)

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 5.2.5. Collate summaries ####
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        b_LOOCV_out <- rbind(   b_LOOCV_out,
                                b_LOOCV_mask_FGW_g,
                                b_LOOCV_mask_GGW_g,
                                b_LOOCV_mask_KGW_g,
                                b_LOOCV_mask_SGW_g,
                                b_LOOCV_mask_XGW_g,
                                b_LOOCV_mask_XQB_g)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.3. Close loop ####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    }
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 6. estimate pacc
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.1. Set up summary objects
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    apdata_CV$TE <- paste(apdata_CV$L,apdata_CV$AB1,sep = "__")
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.2 Loop for trial-by-environment
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    for(te in unique(apdata_CV$TE))
    {
        te <- unique(apdata_CV$TE)[1]
        
        
        
    }


    