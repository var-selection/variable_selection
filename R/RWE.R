change_rank=function(res2=res2,fixed_effect,random_effect){

  if(!is.na(res2)){

    if(length(res2)==5){
      co_matrix=res2[[2]]
      fre_list=t(res2[[5]][,-1])
      rownames(fre_list)=NULL

      rank_x=as.data.frame(fixed_effect)
      rank_z=as.data.frame(random_effect)
      names(rank_x)=names(rank_z)="name"

      tmp_cof=co_matrix$matrix%>%
        na.omit()
      if(nrow(tmp_cof)==0){
        co_matrix=co_matrix$matrix
      } else{
        co_matrix=tmp_cof
      }
      colnames(fre_list)=colnames(co_matrix)
      bt_time=dim(co_matrix)[1]


      rk=purrr::reduce(lapply(1:bt_time,change_rank_fun,co_matrix=co_matrix,fre_list=fre_list,fixed_effect,random_effect), dplyr::left_join, by = 'name')
      colnames(rk)=c("name",paste0("rank_",seq(1,ncol(rk)-1,1)))

      # rlist::list.save(rk,paste0(path_rank,n_group,"_",n_time,"_",type,id,".rds"))


    }
  }
  return(rk)
}


change_na=function(co,fre){
  if(fre==1){
    co=co
  } else{
    co=NA
  }}


change_rank_fun=function(i,co_matrix,fre_list,fixed_effect,random_effect){

  rk_name_x=as.data.frame(fixed_effect)
  rk_name_z=as.data.frame(random_effect)
  num_var_fix=length(fixed_effect)
  num_var_random=length(random_effect)

  names(rk_name_x)=names(rk_name_z)="name"

  if(i<=nrow(fre_list)){
    if(sum(is.na(co_matrix))==0){
      tmp_co=purrr::map2_dfr(co_matrix[i,],fre_list[i,],change_na)
    } else{
      tmp_co=co_matrix[i,]
    }
    tmp_co_t_x=tmp_co%>%

      dplyr::select({{fixed_effect}})%>%
      data.table::transpose()%>%
      mutate(name=colnames(tmp_co)[1:num_var_fix])%>%
      arrange(abs(V1))%>%
      mutate(rank=seq_len(num_var_fix))%>%
      mutate(rank=ifelse(is.na(V1),V1,(ifelse(V1==0,0,rank))))%>%
      mutate(rank1=max(rank,na.rm=T),
             rank=rank1-rank+1)%>%
      dplyr::select(name,rank)
    tmp_co_t_x=dplyr::left_join(rk_name_x,tmp_co_t_x,by="name")


    tmp_co_t_z=tmp_co%>%

      dplyr::select({{random_effect}})%>%
      data.table::transpose()%>%
      mutate(name=colnames(tmp_co)[(num_var_fix+1):(ncol(fre_list))])%>%
      arrange(abs(V1))%>%
      mutate(rank=seq_len(num_var_random))%>%
      mutate(rank=ifelse(is.na(V1),V1,(ifelse(V1==0,0,rank))))%>%
      mutate(rank1=max(rank,na.rm=T),
             rank=rank1-rank+1)%>%
      dplyr::select(name,rank)
    tmp_co_t_z=dplyr::left_join(rk_name_z,tmp_co_t_z,by="name")

    tmp_co_t=rbind(tmp_co_t_x,tmp_co_t_z)
    print(i)


  }
  return(tmp_co_t)

}

#------------change_prob--------------------------------

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))


rp_change_fun=
  function(rank,plot_fixed_effect=T,plot_random_effect=T,fixed_effect,random_effect,fixed_effect_name,random_effect_name){

    rank_x=as.data.frame(fixed_effect)
    rank_z=as.data.frame(random_effect)
    names(rank_x)=names(rank_z)="name"
    nm_x=rank_x
    nm_z=rank_z
    #id=6
    num_var_fix=length(fixed_effect)
    num_var_random=length(random_effect)


    if(!all(is.na(rank))){


      n_missing_t=sum(colSums(is.na(rank))==num_var_fix+num_var_random)
      uncon_index=which(colSums(is.na(rank))==num_var_fix+num_var_random)
      colnames(rank)=c("name",paste0("rank_",seq_len(ncol(rank)-1)))
      ini=colnames(rank)[2]
      end=colnames(rank)[ncol(rank)]
      #rank[,-1] <- mutate_all(rank[,-1], ~coalesce(.,max(rank[,-1],na.rm=T)+1))
      n_miss=rank%>%
        mutate(across({{ini}}:{{end}},is.na))%>%
        rowwise()%>%
        dplyr::summarize(miss=sum(c_across({{ini}}:{{end}})))%>%
        mutate(name=c(fixed_effect,random_effect))%>%
        arrange(desc(miss))%>%
        filter(miss>=(ncol(rank)-1))
      n_miss_x=n_miss%>%
        filter(name %in% fixed_effect)

      n_miss_z=n_miss%>%
        filter(name %in% random_effect)
      #------------x----------------------------
      ini=colnames(rank)[2]
      end=colnames(rank)[ncol(rank)]
      # remove candidate fixed-effects which were never selected
      # remove the bootscrapping unconverage
      rank_r <- rank %>%
        # dplyr::select(
        #   where(
        #     ~sum(!is.na(.x)) > 0
        #   )
        # )%>%
        filter(name %in% fixed_effect)%>%

        dplyr::select(-name)%>%
        data.table::transpose()%>%
        data.table::setnames(fixed_effect)%>%
        dplyr::select(!n_miss_x$name)


      # calculate mean of rank
      mean_rank=rank %>%
        filter(name %in% fixed_effect)%>%
        dplyr::select(-name)%>%
        mutate_at(vars(one_of(names(rank)[-1])),funs(if_else(is.na(.),0,1)))%>%
        as.data.frame()%>%
        dplyr::rowwise()%>%
        dplyr::summarise(fre=sum(c_across({{ini}}:{{end}})))%>%
        dplyr::mutate(name=fixed_effect)%>%
        dplyr::left_join(rank%>%
                           filter(name %in% fixed_effect)%>%
                           dplyr::select(-name)%>%
                           as.data.frame()%>%
                           dplyr::rowwise()%>%
                           dplyr::summarise(rank=sum(c_across({{ini}}:{{end}}),na.rm=T))%>%
                           mutate(name=fixed_effect),by="name")%>%
        mutate(mean_rank=rank/fre)%>%
        #dplyr::filter(!name %in%n_miss_x$name)%>%
        dplyr::select(name,mean_rank,fre)%>%
        arrange(mean_rank)
      bt_time=dim(rank)[2]-1
      # empirical probability for fixed effect
      prob=rank %>%
        filter(name %in% fixed_effect)%>%
        dplyr::select(-name)%>%
        mutate_at(vars(one_of(names(rank)[-1])),funs(if_else(is.na(.),0,if_else(.==max(rank[,-1],na.rm=T),0,1))))%>%
        as.data.frame()%>%
        dplyr::rowwise()%>%
        dplyr::summarise(fre=sum(c_across({{ini}}:{{end}})))%>%
        dplyr::mutate(name=fixed_effect)%>%
        dplyr::mutate(prob=fre/bt_time)%>%
        dplyr::select(name,prob)%>%
        arrange(desc(prob))
      #%>%
      #  dplyr::filter(!name %in%n_miss_x$name)


      nms_zero=prob%>%
        dplyr::filter(!name %in%n_miss_x$name)
      max_rank=max(rank_r,na.rm=T)


      nm_x=as.data.frame(colnames(rank_r))
      names(nm_x)="name"

      # change rank as ranking format
      rank_r=rank_r%>%
        #dplyr::select((!(mean_rank$name[mean_rank$fre<=4]))&(!nms_zero$name[nms_zero$prob==0]))%>%
        PlackettLuce::as.rankings(.)

      # use PL model to calculate the worth paramters

      res_pl <- tryCatch(PlackettLuce::PlackettLuce(rank_r),error=function(e){return(NA)})
      if(!all(is.na(res_pl))){
        cf_t=tryCatch(PlackettLuce::qvcalc(res_pl,ref=NULL),error=function(e){return(NA)})
        if(!all(is.na(cf_t))){

          plot_x_name=
            fixed_effect_name%>%
            as.data.frame()%>%
            dplyr::filter(!new %in%n_miss_x$name)


          cf_t$qvframe=cf_t$qvframe%>%
            mutate(x=plot_x_name[,1])%>%
            arrange(desc(estimate))%>%
            `rownames<-`(.[,5])%>%
            dplyr::select(-x)

          if(plot_fixed_effect==T){
            plot(cf_t,ylab = "Worth (log)", main = "importance for fixed effects")
          }
        } else{
          cat("Do not have enough data to plot, please decrease the value of lambda")
        }

        rp=tryCatch(itempar(res_pl,ref=NULL)%>% as.vector(),error=function(e){return(NA)})
        rp1=coef(res_pl,ref=NULL)
        tie_ind=grep("tie", names(rp1))
        if(length(tie_ind)!=0){
          rp1=rp1[-tie_ind]
        }
        rp2=exp(rp1)/sum(exp(rp1))
        #--------plots-----------------------------
        rp=
          rp%>%
          as.data.frame()%>%
          setNames("rp")%>%
          arrange(desc(rp))



        rp1=
          rp1%>%
          as.data.frame()%>%
          setNames("rp")%>%
          arrange(desc(rp))

        rp=data.frame(name=rownames(rp1),rp=rp)




        rp=left_join(prob,rp,by="name")%>%
          dplyr::select(name,rp)
        rp[is.na(rp)]=0

        if(!all(is.na(rp))){


          # rp=rp[1:num_var]
          rank_x_rep=rp%>%
            arrange(desc(rp))%>%
            dplyr::left_join(prob,by="name")%>%
            dplyr::rename("new"=name)%>%
            dplyr::left_join(fixed_effect_name%>%as.data.frame(),by="new")%>%
            dplyr::rename("variable"=old)%>%
            dplyr::select(-new)%>%
            dplyr::rename(!!("RWE"):=rp)%>%
            dplyr::rename(!!("EP"):=prob)%>%
            dplyr::select(variable,RWE,EP)
        } else{
          rank_x_rep=cbind(nm_x,rp2)%>%
            arrange(desc(rp2))%>%
            dplyr::left_join(prob,by="name")%>%
            dplyr::rename("new"=name)%>%
            dplyr::left_join(fixed_effect_name%>%as.data.frame(),by="new")%>%
            dplyr::rename("variable"=old)%>%
            dplyr::select(-new)%>%
            dplyr::rename(!!("RWE"):=rp2)%>%
            dplyr::rename(!!("EP"):=prob)%>%
            dplyr::select(variable,RWE,EP)
        }
        #rank_x=dplyr::left_join(rank_x,rank_x_rep,by="name")
      } else{ rank_x_rep=NA}
      #----------z----------------------------
      if(!all(is.na(rank_x_rep))){

        rank_r <- rank %>%
          filter(name %in% random_effect)%>%
          dplyr::select(-name)%>%
          data.table::transpose()%>%
          data.table::setnames(random_effect)%>%
          dplyr::select(!n_miss_z$name)

        ini=colnames(rank)[2]
        end=colnames(rank)[ncol(rank)]

        mean_rank=rank %>%
          filter(name %in% random_effect)%>%
          dplyr::select(-name)%>%
          mutate_at(vars(one_of(names(rank)[-1])),funs(if_else(is.na(.),0,1)))%>%
          as.data.frame()%>%
          dplyr::rowwise()%>%
          dplyr::summarise(fre=sum(c_across({{ini}}:{{end}})))%>%
          dplyr::mutate(name=random_effect)%>%
          dplyr::left_join(rank%>%
                             filter(name %in% random_effect)%>%
                             dplyr::select(-name)%>%
                             as.data.frame()%>%
                             dplyr::rowwise()%>%
                             dplyr::summarise(rank=sum(c_across({{ini}}:{{end}}),na.rm=T))%>%
                             mutate(name=random_effect),by="name")%>%
          mutate(mean_rank=rank/fre)%>%

          #dplyr::filter(!name %in%n_miss_z$name)%>%
          dplyr::select(name,mean_rank,fre)

        prob=rank %>%
          filter(name %in% random_effect)%>%
          dplyr::select(-name)%>%
          mutate_at(vars(one_of(names(rank)[-1])),funs(if_else(is.na(.),0,if_else(.==max(rank[,-1],na.rm=T),0,1))))%>%
          as.data.frame()%>%
          dplyr::rowwise()%>%
          dplyr::summarise(fre=sum(c_across({{ini}}:{{end}})))%>%
          dplyr::mutate(name=random_effect)%>%
          dplyr::mutate(prob=fre/400)%>%
          dplyr::select(name,prob)%>%
          arrange(desc(prob))
        nms_zero=prob%>%
          dplyr::filter(!name %in%n_miss_z$name)

        max_rank=max(rank_r,na.rm=T)
        nm_z=as.data.frame(colnames(rank_r))
        names(nm_z)="name"



        rank_r=rank_r%>%
          PlackettLuce::as.rankings(.)

        res_pl <- tryCatch(PlackettLuce::PlackettLuce(rank_r),error=function(e){return(NA)})


        if(!all(is.na(res_pl))){
          cf_t=tryCatch(PlackettLuce::qvcalc(res_pl,ref=NULL),error=function(e){return(NA)})
          if(!all(is.na(cf_t))){

            plot_z_name=
              random_effect_name%>%
              as.data.frame()%>%
              dplyr::filter(!new %in%n_miss_z$name)


            cf_t$qvframe=cf_t$qvframe%>%
              mutate(x=plot_z_name[,1])%>%
              arrange(desc(estimate))%>%
              `rownames<-`(.[,5])%>%
              dplyr::select(-x)

            if(plot_random_effect==T){
              p_random= plot(cf_t,ylab = "Worth (log)", main = "importance for random effects")
            }

          } else{
            cat("Do not have enough data to plot, please decrease the value of lambda")
          }

          rp=tryCatch(itempar(res_pl,ref=NULL)%>% as.vector(),error=function(e){return(NA)})
          rp1=coef(res_pl,ref=NULL)
          tie_ind=grep("tie", names(rp1))

          if(length(tie_ind!=0)){
            rp1=rp1[-tie_ind]
          }
          rp2=exp(rp1)/sum(exp(rp1))
          #--------plots-----------------------------
          rp=
            rp%>%
            as.data.frame()%>%
            setNames("rp")%>%
            arrange(desc(rp))



          rp1=
            rp1%>%
            as.data.frame()%>%
            setNames("rp")%>%
            arrange(desc(rp))

          rp=data.frame(name=rownames(rp1),rp=rp)




          rp=left_join(prob,rp,by="name")%>%
            dplyr::select(name,rp)
          rp[is.na(rp)]=0

          if(!all(is.na(rp))){


            # rp=rp[1:num_var]
            rank_z_rep=rp%>%
              arrange(desc(rp))%>%
              dplyr::left_join(prob,by="name")%>%
              dplyr::rename("new"=name)%>%
              dplyr::left_join(random_effect_name%>%as.data.frame(),by="new")%>%
              dplyr::rename("variable"=old)%>%
              dplyr::select(-new)%>%
              dplyr::rename(!!("RWE"):=rp)%>%
              dplyr::rename(!!("EP"):=prob)%>%
              dplyr::select(variable,RWE,EP)

          }else{
            rank_z_rep=cbind(nm_z,rp2)%>%
              arrange(desc(rp2))%>%
              dplyr::left_join(prob,by="name")%>%
              dplyr::rename("new"=name)%>%
              dplyr::left_join(random_effect_name%>%as.data.frame(),by="new")%>%
              dplyr::rename("variable"=old)%>%
              dplyr::select(-new)%>%
              dplyr::rename(!!("RWE"):=rp2)%>%
              dplyr::rename(!!("EP"):=prob)%>%
              dplyr::select(variable,RWE,EP)
          }
        } else{rank_z_rep=NA}



        #rank_z=dplyr::left_join(rank_z,rank_z_rep,by="name")

        return(list(importance_fixed_effect=rank_x_rep,importance_random_effect=rank_z_rep))
      }

    }
  }
