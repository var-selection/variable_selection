header.true <- function(data) {
  names(data) <- as.character(unlist(data[1,]))
  data[-1,]
}

data_change_format=function(data,fixed_effect,random_effect,ycol,cluster_id,time){

  dt=data%>%
    dplyr::select(fixed_effect,random_effect,ycol,cluster_id,time)
  x1=rep(1,nrow(dt))
  z1=rep(1,nrow(dt))
  # len of fixed and random effect
  num_fixed=length(fixed_effect)
  num_random=length(random_effect)
  # get new name of fixed and random effect
  random_effect_new <-paste0("z",seq(1,num_random+1,1))
  fixed_effect_new <- paste0("x",seq(1,num_fixed+1,1))

  dt=cbind(x1,z1,dt)


  random_effect_name <-cbind(old=c("intercept",random_effect),new=random_effect_new)

  fixed_effect_name <-cbind(old=c("intercept",fixed_effect),new=fixed_effect_new)

  #-------change name of fix and random effect--------------------

  dt=dt%>%
    dplyr::select(fixed_effect_new,random_effect_new,ycol,cluster_id,time)
  names(dt)=c(fixed_effect_new,random_effect_new,"y","participant_id","wave")

  return(list(dt,fixed_effect_name,random_effect_name))

}

# d=data_change_format(data)
# df=d[[1]]
# fixed_effect_name=d[[2]]
# random_effect_name=d[[3]]

