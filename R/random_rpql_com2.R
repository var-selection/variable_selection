# step 2 of random rpql

random_rpql_com2  <- function(q2_fix = 5,
                              q2_random=3,
                              bt=10,
                              data=data,
                              fixed_effect,
                              random_effect,
                              ycol,
                              id_vector,
                              no_scale,
                              cluster_id,
                              time,
                              family = "gaussian",
                              #cov.groups=NULL,
                              #categorical_levels = c(0),
                              pen.type="adl",
                              result=res1,
                              lam=exp(seq(from=log(0.55),to=log(0.001),length.out=70)),
                              verbose = TRUE,
                              random_draw_sample = 0.5,
                              ci_criteria,
                              use_lmer=T,# whether the weights of ALasso from lmer
                              use_fre_fun1=T # use the frequency list to calculate the important index from function1
){#start function2

  if(is.null(q2_fix)){
    q2_fix=floor(length(fixed_effect)/2)
  }

  if(is.null(q2_random)){
    q2_random=length(t(unique(data%>%dplyr::select(time))))-1
  }
  num_w= length(t(unique(data%>%dplyr::select(time))))
  if(num_w<=q2_random){
    stop("number of random effect used in each bootscrapping step need to be less than the within group size")
  } else{

    result_matrix<- vector(mode = "list", length = bt)

    var_list_x            <- c(fixed_effect, random_effect, id_vector)
    var_list              <- c(var_list_x, ycol)
    var_model_selection_x <- setdiff(var_list_x, id_vector)
    fl                    <- length(fixed_effect)
    rl                    <- length(random_effect)
    rlcov                 <- (rl*(rl+1))/2-rl

    rep=floor(q2_fix/q2_random) # proportion between fix and random


    M    <- result$est%>%
      purrr::map_dfc(., as.character)%>%
      purrr::map_dfc(., as.numeric)#from function1
    m    <-dim(M)[2]


    p    <-dim(data[,var_list_x])[2]
    n    <-dim(data[,var_list_x])[1]

    #x_name<-result$x_name #from function1
    x_name<-c(fixed_effect, random_effect) #this name is for prob
    x_name_rand <- outer(random_effect, random_effect, FUN = "paste0") #add on Feb 2022
    x_name_full <- c(fixed_effect, diag(x_name_rand), x_name_rand[lower.tri(x_name_rand, diag = F)]) #add on Feb 2022, the order of first 13 should be matached to N

    v            <-dim(data[,var_model_selection_x])[2]
    # N            <-array(0,dim=c(bt,v,1))
    # NN           <-matrix(0,bt,v+rlcov)


    #prepare names
    check_name <- x_name %>% as.data.frame() %>%
      dplyr::rename(name = ".") %>%
      tibble::rownames_to_column(., "number") %>%
      dplyr::mutate(
        V = paste0("V", number),
        no_number = gsub("\\d{1}$","\\" ,name),
      )

    PO_list <- list() #for probability result table
    CC_list <- list() #for coefficient result table


    check_name_full <-  x_name_full %>% as.data.frame() %>%
      dplyr::rename(name = ".") %>%
      tibble::rownames_to_column(., "number") %>%
      dplyr::mutate(
        V = paste0("V", number),
        no_number = sub("(.).*", "\\1",  name),
      )

    check_namename<-as.character(check_name$name)

    #prepare data
    M <- as.data.frame(M) %>% purrr::map_dfc(., as.numeric)
    M <- data.table::setnames(M, old=colnames(M), new=check_namename) %>% as.matrix()

    #IMPORTANT!!! #Calculate the importance index
    converge_bt  <- as.numeric(dim(na.omit(M))[1]) #any missing indicated non-converge from function1

    #update calculate the importance index by using actual frequency
    if(use_fre_fun1){

      Impms        <-as.data.frame(colSums(abs(M), na.rm = T))%>%
        mutate(rownames(.))
      colnames(Impms)=c("sum","name")
      Impms=
        Impms%>%
        left_join(result$fre,by="name")%>%
        dplyr::select(name,sum,total)%>%
        rowwise()%>%
        mutate(prob=sum/total)%>%
        mutate(prob=if_else(total==0,0,prob))%>%
        ungroup()%>%
        dplyr::select(name,prob)



      nms=Impms$name
      Impms=as.matrix(Impms[,2])
      rownames(Impms)=nms

    } else{ # using total number of bootstrapping

      Impms        <-as.matrix(colSums(abs(M), na.rm = T)/converge_bt )
    }

    y            <- as.matrix(data[,ycol])
    scale_x_name <- setdiff(var_list_x, no_scale)
    xs <- data %>%
      dplyr::select(all_of(var_list_x) )%>%
      mutate_at(., vars(var_list_x), as.numeric) %>% #update July 2022 --> assuming no categorical predictor
      mutate_at(., vars(scale_x_name), scale) #standardize all predictors
    fc2=0

    # full variable lists
    name_res=as.data.frame(c(fixed_effect,random_effect))
    names(name_res)="name"
    fre_list=name_res # list of number of variables selected

    j=1
    while(j<=bt){
      if (verbose) message(sprintf("Starting bootstrap % 3i of %i", j, bt))
      unique_id2 <- as.data.frame(unique(data[, cluster_id]))
      names(unique_id2)="group_id"
      random_unique <- as.numeric(sample(unique_id2$group_id, round(dim(unique_id2)[1])*random_draw_sample,  replace = T)) ##!!changed

      random_unique1<- random_unique %>%
        as.data.frame() %>%
        dplyr::rename(cluster_id = ".") %>%
        dplyr::arrange(cluster_id) %>%
        dplyr::mutate(
          cluster_id = as.numeric(cluster_id),
        ) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::mutate(
          time = as.numeric(rep(1:dplyr::n()))
        ) %>%
        dplyr::summarise(
          weight = max(time) #weight is to summarize how many times of being sampled
        ) %>%
        dplyr::ungroup()


      IVs2_x <- as.data.frame(xs) %>%
        dplyr::rename(time = time,
                      cluster_id = cluster_id)

      smp2 <- IVs2_x %>% #update the id infor
        dplyr::select(time, cluster_id) %>%
        dplyr::right_join(random_unique1, by="cluster_id") %>%  #depends on the random_unique1
        dplyr::group_split(cluster_id) %>%
        purrr::map(function(data){
          rep(data$time, times=max(data$weight)) %>% #replicate by weight
            as.data.frame() %>%
            dplyr::rename(
              time = "."
            ) %>%
            dplyr::mutate(cluster_id = rep(unique(data$cluster_id))) %>%
            dplyr::group_by(time) %>%
            dplyr::mutate(
              seq = rep(1:dplyr::n()),
              b_cluster_id = as.integer(paste0(cluster_id, seq))
            ) %>%
            dplyr::select(
              time, b_cluster_id, cluster_id
            )}) %>% do.call(rbind.data.frame, .)%>%
        group_by(b_cluster_id)%>%
        mutate(cluster_id_rpql=cur_group_id())


      #if (class(data[,ycol] ) == "numeric" | family == "gaussian")
      if (family == "gaussian"){
        outcome2 <- data %>%
          dplyr::select(cluster_id, time, ycol) %>%
          purrr::map_dfc(., as.numeric) %>% #for continuous outcome, can do this
          dplyr::mutate_at(., vars(ycol), scale) %>% #DV needs to be scale
          dplyr::arrange(cluster_id) %>%
          dplyr::rename(time = time, cluster_id = cluster_id)

      }else{
        outcome2 <- data %>%
          dplyr::select(cluster_id, time, ycol) %>%
          dplyr::mutate_at(., vars(c(cluster_id, time)), as.numeric) %>%
          arrange(cluster_id) %>%
          dplyr::rename(time = time, cluster_id = cluster_id) #DV no scale
      }#end of dv scale issue


      IVs2 <- smp2 %>% dplyr::left_join(IVs2_x, by=c("cluster_id", "time"))
      DV2  <- smp2 %>% dplyr::left_join(outcome2, by=c("cluster_id", "time"))

      #if(dim(IVs2)[1] != dim(DV2)[1]){warning("The number of rows of x is not as same as the number of the outcome")}

      x1   <- IVs2 %>% dplyr::arrange(cluster_id, time) %>% as.matrix()
      y1   <- DV2 %>% dplyr::arrange(cluster_id, time) %>%
        dplyr::select(b_cluster_id, time, ycol) %>%
        as.matrix()

      x2_data <- IVs2 %>%
        dplyr::ungroup() %>%
        dplyr::left_join(DV2, by=c("cluster_id", "time", "b_cluster_id","cluster_id_rpql")) %>%
        dplyr::arrange(cluster_id, b_cluster_id,cluster_id_rpql, time) %>%
        dplyr::select(-time)

      #the following condition is for all continuous predictors random draw
      #if (is.null(categorical_var)){ #start continuous predictors random draw;
      x2_data_selective_x_only <- x2_data %>%
        dplyr::select(all_of(var_model_selection_x))%>%
        dplyr::mutate_all(., as.numeric)%>%
        # change 9/5 choose w1 and z1 also
        dplyr::select(starts_with("x")) #, -x1) #CHANGE for package, specify the intercept

      x2_data_selective_z_only <- x2_data %>%
        dplyr::select(all_of(var_model_selection_x))%>%
        dplyr::mutate_all(., as.numeric)%>%
        # change 9/5 choose w1 and z1 also
        dplyr::select(starts_with("z"))#,-z1 ) #CHANGE for package, specify the intercept

      x2_lmer_use <- x2_data %>%
        dplyr::ungroup() %>%
        dplyr::select(all_of(var_model_selection_x), ycol, b_cluster_id,cluster_id_rpql) %>%
        dplyr::mutate_all(., as.numeric) #assume all continuous predictors

      ## method 2 half sample from fixed effect half from random effect
      # fix

      # change 9/5 choose w1 and z1 also
      name_fix    <- gtools::mixedsort(grep("x",check_name$name,value=T)) #update 06/30, replace rownames function
      #loc_x1      <- which(rownames(Impms)=="x1") #F1 and F2 different part
      #name_fix    <- name_fix[name_fix!="x1"]
      fix_ind     <- grep("x",rownames(Impms))
      #fix_ind     <- fix_ind[fix_ind!=loc_x1] #F1 and F2 different part
      impms_fix   <- Impms[fix_ind]
      # random
      name_random <- gtools::mixedsort(grep("z",rownames(Impms),value=T))
      #loc_z1      <- which(rownames(Impms)=="z1")
      #name_random <- name_random[name_random!="z1"]
      random_ind  <- grep("z",rownames(Impms))
      #random_ind  <- random_ind[random_ind!=loc_z1]
      impms_random<- Impms[random_ind]

      # number of non-zero parameters
      # size_fixed  <- ceiling(length(impms_fix)/(length(impms_fix)+length(impms_random))*q2_fix) #q2 also needs to set up a default value for package
      # size_random <-ceiling(length(impms_random)/(length(impms_fix)+length(impms_random))*q2_random)
      size_fixed  <- q2_fix #q2 also needs to set up a default value for package
      size_random <- q2_random

      w1          <- sample(name_fix,size=size_fixed,prob=impms_fix)#update by name rather than order
      w1=gtools::mixedsort(w1)

      # for (i in 1:rep){
      w2          <- sample(name_random,size=size_random,prob=impms_random)
      w2=gtools::mixedsort(w2)

      fre=cbind(as.data.frame(c(w1,w2)),rep(1,q2_fix+q2_random))
      names(fre)=c("name","fre")

      fre <- name_res %>%
        dplyr::left_join(fre, by="name") %>%
        dplyr::mutate(
          unselect = rep(0),
          fre = dplyr::coalesce(as.numeric(fre), unselect)
        ) %>%
        dplyr::select(-unselect)






      x2_data2_x          <- cbind(x2_data_selective_x_only[, w1])
      x2_data2_z         <- cbind(x2_data_selective_z_only[, w2])


      # change 9/5 choose w1 and z1 also
      x2_data3 <- names(cbind(x2_data_selective_x_only,(x2_data_selective_z_only)))%>%


        as.data.frame() %>%
        dplyr::rename(name = ".") %>%
        tibble::rownames_to_column(., "number")


      x2_data4 <- names(cbind(x2_data2_x,x2_data2_z)) %>%
        as.data.frame() %>%
        dplyr::rename(name = ".") %>%
        left_join(x2_data3,by="name")



      b_fix_effect    <- x2_data4 %>% dplyr::filter(name %in% fixed_effect) %>% dplyr::select(name) %>% dplyr::pull(name)
      b_random_effect <- x2_data4 %>% dplyr::filter(name %in% random_effect) %>% dplyr::select(name) %>% dplyr::pull(name)

      x <- as.matrix(x2_lmer_use[, b_fix_effect])
      y <- x2_lmer_use[, ycol] %>% dplyr::pull(y)
      z <- as.matrix(x2_lmer_use[, b_random_effect])


      list2 <- list(id = list(cluster_id_rpql = x2_lmer_use$cluster_id_rpql), X=x, Z=list(cluster_id_rpql=z), Y=y)
      if(pen.type=="adl" & use_lmer==T){ # use lmer to calculate initial value
        if (family == "gaussian"){#for lmer DV continue
          fit2  <- lme4::lmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use)
        } # # gaussian
        if (family %in% "binomial"){#for glmer DV binary
          fit2  <- lme4::glmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use,family = "binomial")} #binomial
        xx    <- x %>% as.data.frame()

        fit_sat2 <- build.start.fit(fit2, gamma = 2, cov.groups = NULL) #no categorical predictor --> CHANGE for package
        fit_final2 <- tryCatch(rpqlseq(y = list2$Y, X = list2$X, Z = list2$Z, id = list2$id,
                                       lambda = lam, pen.type = "adl", hybrid.est = TRUE,
                                       pen.weights = fit_sat2$pen.weights, start = fit_sat2,restarts = 50),
                               error=function(e){return(NA)})
      } else if(pen.type=="adl" & use_lmer==F){ # do not use lmer
        weight_f=cbind(names=rownames(Impms),as.data.frame(Impms))
        colnames(weight_f)[2]="prob"
        weight_f=weight_f%>%
          filter(names %in% colnames(x))%>%
          dplyr::select(prob)%>%
          t()%>%as.data.frame()%>%
          dplyr::select(colnames(x))%>%
          as.vector()
        weight_f=weight_f^-1
        weight_f=as.vector(t(weight_f))
        weight_r=cbind(names=rownames(Impms),as.data.frame(Impms))
        colnames(weight_r)[2]="prob"
        weight_r=weight_r%>%
          # rename(prob=V1)%>%
          filter(names %in% colnames(z))%>%
          dplyr::select(prob)%>%
          t()%>%as.data.frame()%>%
          dplyr::select(colnames(z))
        weight_r=weight_r^-1
        weight_r=as.vector(t(weight_r))

        wt=list(fixed=weight_f,random=list(cluster_id_rpql=weight_r))


        fit_final2 <- tryCatch(rpqlseq(y = list2$Y, X = list2$X, Z = list2$Z, id = list2$id,
                                       lambda = lam, pen.type = "adl", hybrid.est = TRUE,pen.weights=wt,
                                       restarts = 50),
                               error=function(e){return(NA)})

      } else if(pen.type=="lasso"){
        fit_final1 <- tryCatch(rpqlseq(y = list1$Y, X = list1$X, Z = list1$Z, id = list1$id,
                                       lambda = lam, pen.type = "lasso", hybrid.est = TRUE,
                                       restarts = 50),
                               error=function(e){return(NA)})
      }



      if(is.na(fit_final2)){j=j} else{
        j=j+1
      }


      if(!is.na(fit_final2)){
        final_result <- fit_final2$best.fit[[ci_criteria]] #save all results
        result2      <- summary(fit_final2$best.fit[[ci_criteria]]) #save summary result
        fixcoef2     <- as.data.frame(result2$fixef) %>%  tibble::rownames_to_column(., "name") %>% dplyr::rename(value = `result2$fixef`)
        randcov2     <- as.matrix(result2$ran.cov$cluster_id_rpql)
        randcoef2    <- diag(randcov2)%>%
          as.data.frame() %>%
          cbind(name=as.character(b_random_effect),.)%>%
          dplyr::rename(value = ".")
        randcoef3    <-  randcov2[lower.tri(randcov2, diag = T)]
        togethercoef <- rbind(fixcoef2, randcoef2)
        fre_list=left_join(fre_list,fre,by="name")

        name_res=as.data.frame(c(fixed_effect,random_effect))
        names(name_res)="name"
        x2_data5_cc <- name_res %>% #this coefficient table does not include the covariance values
          dplyr::left_join(togethercoef, by="name") %>%
          dplyr::mutate(
            unselect = rep(NA),
            value = dplyr::coalesce(as.numeric(value), unselect)
          ) %>%
          dplyr::select(name, value) %>%
          # dplyr::arrange(name) %>%
          t() %>%
          as.data.frame() %>%
          header.true()


        x2_data5_po <- name_res %>%
          dplyr::left_join(togethercoef, by="name") %>%
          #purrr::map_dfc(., as.numeric) %>%
          dplyr::mutate(
            unselect = rep(0),
            value1 = dplyr::coalesce(as.numeric(value), unselect),
            value = dplyr::if_else(value1 != 0, 1, 0)
          ) %>%
          dplyr::select(name, value) %>%
          # dplyr::arrange(name) %>%
          t() %>%
          as.data.frame()%>%
          header.true()




      }else{
        final_result <- list() #null list
        fc2=fc2+1
        x2_data5 <- x2_data3 %>%
          dplyr::mutate(
            value = rep(NA_real_)
          ) %>%
          dplyr::select(name, value) %>%
          dplyr::arrange(name) %>%
          t() %>%
          as.data.frame() %>%
          header.true()
        x2_data5_cc <- x2_data5
        x2_data5_po <- x2_data5


      }#end of managing result2

      #rm(random_unique, random_unique1, w1,w2) #delete the random draw id

      PO_list[[j]] <- x2_data5_po
      CC_list[[j]] <- x2_data5_cc


      # }
    } #end of boostrapping


    PO <- do.call(rbind, PO_list) %>% as.data.frame() %>% purrr::map_dfc(., as.character)%>%purrr::map_dfc(., as.numeric)
    CC <- do.call(rbind, CC_list) %>% as.data.frame() %>% purrr::map_dfc(., as.character)%>%purrr::map_dfc(., as.numeric)
    converge_bt  <- as.numeric(dim(na.omit(PO))[1]) #any missing indicated non-converge from function1





    fre_res=fre_list%>%
      rowwise()%>%
      mutate(total=sum(c_across(where(is.numeric)),na.rm=T))%>%
      ungroup()%>%
      dplyr::select(name,total)



    coefficient_matrix   <- as.matrix(abs(colSums(CC, na.rm = T)/converge_bt) )
    probability_matrix   <- as.matrix(colSums(PO, na.rm = T)/converge_bt)



    #make a little bit improvement:
    ds_prob <- as.data.frame(probability_matrix) %>%
      data.table::setDT(keep.rownames = T) %>%
      dplyr::rename(prob = `V1`, name = rn) %>%
      dplyr::arrange(desc(prob))

    ds_coef <- as.data.frame(coefficient_matrix) %>%
      data.table::setDT(keep.rownames = T) %>%
      dplyr::rename(estimate = `V1`, name = rn) %>%
      dplyr::arrange(desc(abs(estimate)))






    round2_res<-
      list(probability = list(matrix = PO, summary = ds_prob),
           coefficient = list(matrix = CC, summary = ds_coef),
           raw_rpql_result  = final_result,
           frequency=fre_res,
           fre_matrix=fre_list)





    return(round2_res)

  }

}#end of function2
