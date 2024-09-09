#' random rPQL
#'@description
#'
#'random_rpql() allows the user to do variable selection for fixed effect and random effect in GLMMs model by using rpql(Hui,et al,2017) and random lasso(Wang, et al., 2011)
#' @usage random_rpql(data,q1_fix,q1_random,bt1,fixed_effect,random_effect,ycol,id_vector,no_scale,cluster_id,within_cluster_id,family="gaussian",lam1,verbose=T,ci_criteria1=5,pen.type1="adl",random_draw_sample1=0.5,q2_fix,q2_random,bt2,lam2,pen.type2="adl",ci_criteria2=5,random_draw_sample2=0.5,use_lmer=F,use_fre_fun1=T)
#' @return variable importance plots and tables for fixed effects and random effects
#' @param data dataset without missing data
#' @param q1_fix number of fixed effects used in each bootstrapping step in first step of random rpql
#' @param q1_random number of random effects used in each bootstrapping step in first step of random rpql
#' @param bt1  number of bootstrapping in first step of random rpql
#' @param fixed_effect candidate variables used to fixed effects
#' @param random_effect candidate variables used to random effects
#' @param ycol  outcome variables
#' @param id_vector column names for groups (cross-sectional design)/participants id(for longitudinal design) and within-units id
#' @param no_scale  variable names of discrete variables
#' @param cluster_id variable name for participants id/group id
#' @param within_cluster_id variable name for time/id for participant within each unit(participant/group)
#' @param family The distribution for the outcome in GLMM.  Supported arguments include: gaussian(), binomial().
#' @param lam1 a vector for tuning parameters used in rPQL for each bootstrapping step in first step
#' @param ci_criteria1 criteria to choose the optimal tuning parameters for step 1, details refers to Hui,et al,2017
#' @param pen.type1 the type of penalty used in random rpql first step. supported argument include: "lasso" for standard lasso(Tibshirani,1996),"adl"for adaptive lasso(Zou,2006)
#' @param random_draw_sample1 the proportion sample used in each bootstrapping in first step
#' @param q2_fix  number of fixed effects used in each bootstrapping step in second step of random rpql
#' @param q2_random number of random effects used in each bootstrapping step in second step of random rpql
#' @param bt2 number of bootstrapping in second step of random rpql
#' @param lam2 a vector for tuning parameters used in rPQL for each bootstrapping step in first step
#' @param pen.type2 the type of penalty used in random rpql second step
#' @param ci_criteria2 criteria to choose the optimal tuning parameters for step 2, details refers to Hui,et al,2017
#' @param random_draw_sample2 the proportion sample used in each bootstrapping in second step
#' @param use_lmer if pen.type2 is true, whether use lmer to derived the weights of penalty. if use_lmer=T, using lmer to calculate weights for "adl". otherwise, the weights for "adl" are calculated by the importance index obtained by first step
#' @param use_fre_fun1 whether used frequency list to calculate importance index
#' @param verbose show the running details
#'
#' @return
#' @export
#' @import rpql
#' @import dplyr
#' @import lme4
#' @import purrr
#' @import gtools
#' @import PlackettLuce
#' @import data.table
#' @import tibble
#' @importFrom methods new
#' @importFrom stats coef df na.omit setNames



#'
#'
#' @examples
#' \dontrun{
#' id_vector<- c("participant_id", "wave")
#'
#' no_scale <- c("participant_id", "wave")
#' cluster_id="participant_id"
#' within_cluster_id="wave"
#' ycol <- "y"
#' family<-"gaussian"
#' lam1=exp(seq(from=log(1e-4),to=log(1e-8),length.out=100))
#' lam2=exp(seq(from=log(1e-1),to=log(1e-8),length.out=100))
#'  random_rpql(data,q1_fix=4,q1_random=3,bt1=200,fixed_effect,random_effect,ycol,id_vector,no_scale,cluster_id,within_cluster_id,family="gaussian",lam1,verbose=T,ci_criteria1=5,pen.type1="adl",random_draw_sample1=0.5,q2_fix=4,q2_random=3,bt2=200,lam2,pen.type2="adl",ci_criteria2=4,random_draw_sample2=0.5,use_lmer=F,use_fre_fun1=T)
#'  }
#' @references
#' Hui, F. K. C., Mueller, S., and Welsh, A.H. (2017). Hierarchical Selection of Fixed and Random Effects in Generalized Linear Mixed Models. Statistica Sinica, 27, 501-518.
#' Wang, S., Nan, B., Rosset, S., & Zhu, J. (2011). Random lasso. The annals of applied statistics, 5(1), 468.
#' @examples
random_rpql=
  function(data,q1_fix,q1_random,bt1,
           fixed_effect,random_effect,
           ycol,id_vector,no_scale,cluster_id,
           within_cluster_id,family="gaussian",lam1,
           verbose=T,ci_criteria1=5,pen.type1="adl",random_draw_sample1=0.5,
           q2_fix,q2_random,bt2,lam2,pen.type2="adl",ci_criteria2=5,
           random_draw_sample2=0.5,use_lmer=F,use_fre_fun1=T
  ){
    options(warn=-1)
    d=data_change_format(data,fixed_effect,random_effect,ycol,cluster_id,time=within_cluster_id)
    df=d[[1]]
    fixed_effect_name=d[[2]]
    random_effect_name=d[[3]]
    fixed_effect_n=d[[2]][,2]
    random_effect_n=d[[3]][,2]


    res1=random_rpql_com1 (q1_fix=q1_fix, # number of fixed effect used in each bootscrap step
                           q1_random=q1_random,# number of random effect used in each bootscrap step
                           bt=bt1, # time of bootscrapping
                           data=df,
                           fixed_effect_n,
                           random_effect_n,
                           ycol,
                           id_vector,
                           no_scale,
                           cluster_id,
                           time=within_cluster_id,
                           family =family,
                           lam=lam1,
                           verbose,
                           ci_criteria=ci_criteria1,# select whether ci used to select the optimal model
                           pen.type=pen.type1,
                           random_draw_sample =random_draw_sample1 #new: manipulate how much sample size of id to draw, if it is 1 means use the orginal ss
    )

    res2=random_rpql_com2(q2_fix = q2_fix,
                          q2_random=q2_random,
                          bt=bt2,
                          data=df,
                          fixed_effect_n,
                          random_effect_n,
                          ycol,
                          id_vector,
                          no_scale,
                          cluster_id,
                          time=within_cluster_id,
                          family = family,

                          pen.type=pen.type2,
                          result=res1,
                          lam=lam2,
                          verbose = TRUE,
                          random_draw_sample = random_draw_sample2,
                          ci_criteria=ci_criteria2,
                          use_lmer,# whether the weights of ALasso from lmer
                          use_fre_fun1) # use the frequency list to calculate the important index from function1}


    rk=change_rank(res2=res2,fixed_effect=fixed_effect_n,random_effect=random_effect_n)
    p=rp_change_fun(rank=rk,fixed_effect=fixed_effect_n,random_effect=random_effect_n,
                    fixed_effect_name=fixed_effect_name,random_effect_name=random_effect_name)

    return(p)}
