# step 1 of random rpql

random_rpql_com1 <- function(q1_fix = 10, # number of fixed effect used in each bootscrap step
                             q1_random=10,# number of random effect used in each bootscrap step
                             bt=10, # time of bootscrapping
                             data=df,
                             fixed_effect,
                             random_effect,
                             ycol,
                             id_vector,
                             no_scale,
                             cluster_id,
                             time,
                             family ,
                             lam=exp(seq(from=log(0.55),to=log(0.001),length.out=70)),
                             verbose = TRUE,
                             ci_criteria,# select whether ci used to select the optimal model
                             pen.type, # penalty methods
                             random_draw_sample = 0.5 #new: manipulate how much sample size of id to draw, if it is 1 means use the orginal ss
){
  if(is.null(q1_fix)){
    q1_fix=floor(length(fixed_effect)/2)
  }

  if(is.null(q1_random)){
    q1_random=length(t(unique(data%>%dplyr::select(time))))-1
  }
  num_w= length(t(unique(data%>%dplyr::select(time))))
  if(num_w<=q1_random){
    stop("number of random effect used in each bootscrapping step need to be less than the within group size")
  } else{


    var_list_x            <- c(fixed_effect, random_effect, id_vector)
    var_list              <- c(var_list_x, ycol)
    var_model_selection_x <- setdiff(var_list_x, id_vector)
    fl                    <- length(fixed_effect) #original fixed effect length
    rl                    <- length(random_effect) #original random effect length
    rlcov                 <- (rl*(rl+1))/2-rl

    #no categorical predictors in this simulation

    p    <-dim(data[,var_list_x])[2] #all variables
    n    <-dim(data[,var_list_x])[1]
    v    <-dim(data[,var_model_selection_x])[2] #only the selection candidates


    x_name<-c(fixed_effect, random_effect)

    check_name <- x_name %>% as.data.frame() %>%
      dplyr::rename(name = ".") %>%
      tibble::rownames_to_column(., "number") %>%
      dplyr::mutate(
        V = paste0("V", number),
        no_number = gsub("\\d{1}$","\\" ,name),)


    M_list <- list() #build a null list for function1 result

    # full variable lists
    name_res=as.data.frame(c(fixed_effect,random_effect))
    names(name_res)="name"
    fre_list=name_res # list of number of variables selected


    rep=floor(q1_fix/q1_random)
    bt_rep=ceiling(bt/rep)
    for (j in 1:bt_rep){#start bootstrapping
      y            <- as.matrix(data[,ycol])
      scale_x_name <- setdiff(var_list_x, no_scale)

      xs <- data %>%
        dplyr::select(all_of(var_list_x) )%>%
        mutate_at(., vars(var_list_x), as.numeric)%>%
        mutate_at(., vars(scale_x_name), scale)


      fc2=0

      if (verbose) message(sprintf("Starting bootstrap % 3i of %i", j, bt))
      unique_id2 <- as.data.frame(unique(data[, cluster_id]))##!!changed
      names(unique_id2)="group_id"
      #IMPORTANT: randomly sample half of the people rather than half of the data.
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


      smp2 <- IVs2_x %>% #update the id info
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
        mutate(cluster_id_rpql=cur_group_id()) # update the cluster_id with a sequential number due to the rpql function




      #if (class(data[,ycol] ) == "numeric" | family == "gaussian")
      if (family == "gaussian"){
        outcome2 <- data %>%
          dplyr::select(cluster_id, time, ycol) %>%
          purrr::map_dfc(., as.numeric) %>% #for continuous outcome, can do this
          dplyr::mutate_at(., vars(ycol), scale) %>% #DV needs to be scale
          dplyr::arrange(cluster_id) %>%
          dplyr::rename(time = time, cluster_id = cluster_id)


        # outcome2_lmer=data %>%
        #   dplyr::select(cluster_id, time, ycol) %>%
        #   purrr::map_dfc(., as.numeric) %>% #for continuous outcome, can do this
        #   dplyr::arrange(cluster_id) %>%
        #   dplyr::rename(time = time, cluster_id = cluster_id)

      }else{ # if the outcome variable is discrete, DV do not scale
        outcome2 <- data %>%
          dplyr::select(cluster_id, time, ycol) %>%
          dplyr::mutate_at(., vars(c(cluster_id, time)), as.numeric) %>%
          arrange(cluster_id) %>%
          dplyr::rename(time = time, cluster_id = cluster_id) #DV no scale
        outcome2_lmer=outcome2

      }

      DV2  <- smp2 %>% dplyr::left_join(outcome2, by=c("cluster_id", "time"))
      IVs2 <- smp2 %>% dplyr::left_join(IVs2_x, by=c("cluster_id", "time"))



      x1   <- IVs2 %>% dplyr::arrange(cluster_id, time) %>% as.matrix()
      y1   <- DV2 %>% dplyr::arrange(cluster_id, time) %>%
        dplyr::select(b_cluster_id, time,cluster_id_rpql, ycol) %>%
        as.matrix()

      x2_data <- IVs2 %>%
        dplyr::ungroup() %>%
        dplyr::left_join(DV2, by=c("cluster_id", "time", "b_cluster_id","cluster_id_rpql")) %>%
        dplyr::arrange(cluster_id, b_cluster_id,cluster_id_rpql,time) %>%
        dplyr::select(-time)




      #TODO start continuous predictors random draw; No if condition here


      x2_data_selective_x_only <- x2_data %>%
        dplyr::select(all_of(var_model_selection_x))%>%
        purrr::map_dfc(., as.numeric) %>% #for no categorical var
        dplyr::select(starts_with("x"))

      x2_data_selective_z_only <- x2_data %>%
        dplyr::select(all_of(var_model_selection_x))%>%
        dplyr::mutate_all(., as.numeric)%>%
        dplyr::select(starts_with("z"))

      x2_lmer_use <- x2_data %>%
        dplyr::ungroup() %>%
        dplyr::select(all_of(var_model_selection_x), ycol, b_cluster_id,cluster_id_rpql) %>%
        purrr::map_dfc(., as.numeric) #assume all continuous predictors



      # change 9/5 choose w1 and z1 also
      # fix
      name_fix  <- gtools::mixedsort(grep("x",check_name$name,value=T)) #replace rownames function
      #loc_x1    <- which(check_name$name=="x1") # used for simulation
      #name_fix  <- name_fix[name_fix!="x1"]
      fix_ind   <- grep("x", check_name$name)
      #fix_ind   <- fix_ind[fix_ind!=loc_x1]

      # random
      name_random  <- gtools::mixedsort(grep("z",check_name$name,value=T))
      #loc_z1       <- which(check_name$name=="z1") #only for simulation, need CHANGE for package
      #name_random  <- name_random[name_random!="z1"]
      random_ind   <- grep("z",check_name$name)
      #random_ind   <- random_ind[random_ind!=loc_z1]



      # proportional of fixed and random which is non-zero
      # number of non-zero parameters
      fix_len <- length(fix_ind)
      ran_len <- length(random_ind)
      total_len <- fix_len+ran_len


      size_fixed  <- q1_fix #q1 needs to set up a default value for package
      size_random <- q1_random


      w1          <- gtools::mixedsort(sample(name_fix,size=size_fixed)) #update by name rather than order
      # the proportion of fix and random with different number of selection times
      for (i in 1:rep){

        w2          <- gtools::mixedsort(sample(name_random,size=size_random))

        # update
        # number of time vairables selected

        fre=cbind(as.data.frame(c(w1,w2)),rep(1,q1_fix+q1_random))

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




        x2_data3 <- names(cbind(x2_data_selective_x_only,(x2_data_selective_z_only)))%>%

          as.data.frame() %>%
          dplyr::rename(name = ".") %>%
          tibble::rownames_to_column(., "number")


        x2_data4 <- names(cbind(x2_data2_x,x2_data2_z)) %>%
          as.data.frame() %>%
          dplyr::rename(name = ".") %>%
          left_join(x2_data3,by="name")


        b_fix_effect <- x2_data4 %>% dplyr::filter(name %in% fixed_effect) %>% dplyr::select(name) %>% dplyr::pull(name)
        b_random_effect <- x2_data4 %>% dplyr::filter(name %in% random_effect) %>% dplyr::select(name) %>% dplyr::pull(name)

        x <- as.matrix(x2_lmer_use[, b_fix_effect])
        y <- x2_lmer_use[, ycol] %>% dplyr::pull(y)
        z <- as.matrix(x2_lmer_use[, b_random_effect])



        list1 <- list(id = list(cluster_id_rpql = x2_lmer_use$cluster_id_rpql), X=x, Z=list(cluster_id_rpql=z), Y=y)

        if(pen.type=="adl"){# use lmer to calculate the initial value
          if (family %in% "gaussian"){#for lmer DV continue
            fit1  <- lme4::lmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use)
          } #cate needs to be factor
          if (family %in% "binomial"){#for glmer DV binary
            fit1  <- lme4::glmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use,family = "binomial")} #cate needs to be factor}

          xx    <- x %>% as.data.frame()
          fit_sat1 <- build.start.fit(fit1, gamma = 2, cov.groups = NULL) #no categorical predictor --> CHANGE for package
          fit_final1 <- tryCatch(rpqlseq(y = list1$Y, X = list1$X, Z = list1$Z, id = list1$id,
                                         lambda = lam, pen.type = "adl", hybrid.est = TRUE,
                                         pen.weights = fit_sat1$pen.weights, start = fit_sat1,restarts = 50),
                                 error=function(e){return(NA)})
        } else if(pen.type=="lasso"){ # do not use lmer as initial value
          fit_final1 <- tryCatch(rpqlseq(y = list1$Y, X = list1$X, Z = list1$Z, id = list1$id,
                                         lambda = lam, pen.type = "lasso", hybrid.est = TRUE,
                                         restarts = 50),
                                 error=function(e){return(NA)})

        }






        if(is.na(fit_final1)){fc2=fc2+1 }
        if(fc2>bt*0.05){break} # if 5% bootscraping cannot converage, break the loop and give NA
        if(!is.na(fit_final1)){
          result      <- summary(fit_final1$best.fit[[ci_criteria]])
          fixcoef     <- as.data.frame(result$fixef) %>%  tibble::rownames_to_column(., "name") %>% dplyr::rename(value = `result$fixef`)
          randcov     <- as.matrix(result$ran.cov$cluster_id_rpql)
          randcoef    <- diag(randcov)%>% as.data.frame() %>% cbind(b_random_effect,.)  %>% dplyr::rename(value = ".", name=b_random_effect)
          togethercoef<- rbind(fixcoef, randcoef)
          # frequency list of variables selected  if the function is converage combined the fre into fre_list otherwise ignore the fre

          fre_list=left_join(fre_list,fre,by="name")

          name_res=as.data.frame(c(fixed_effect,random_effect))
          names(name_res)="name"
          x2_data5 <- name_res %>%
            dplyr::left_join(togethercoef, by="name") %>%
            dplyr::mutate(
              unselect = rep(0),
              value = dplyr::coalesce(as.numeric(value), unselect)
            ) %>%
            dplyr::select(name, value) %>%
            #dplyr::arrange(name) %>%
            t() %>%
            as.data.frame() %>%
            header.true()


        } #end of non-null fit_final1
        else{
          fc2=fc2+1
          x2_data5 <- x2_data3 %>%
            dplyr::mutate(
              value = rep(NA)
            ) %>%
            dplyr::select(name, value) %>%
            dplyr::arrange(name) %>%
            t() %>%
            as.data.frame() %>%
            header.true()
        }#end of null fit_final1



        M_list[[(j-1)*rep+i]] <- x2_data5

      }}#end of bootstrapping


    M <- do.call(rbind, M_list)

    fre_res=fre_list%>%
      rowwise()%>%
      mutate(total=sum(c_across(where(is.numeric)),na.rm=T))%>%
      ungroup()%>%
      dplyr::select(name,total)



    round1_res <- list(est=as.data.frame(M), fre=fre_res,family=family)

    return(round1_res)
  }

}#end of mirl_rpql_selection_com1
