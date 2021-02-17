library(purrr)
library(gtools)
separate_pvals<-function(pvals) #returns #1 lnamed list of sig pvalues #2 named list of non-sig pvalues
{
  sig=list()
  non_sig=list()
  
  for( i in 1:length(pvals))
  {
    if(pvals[i]<0.05)
    {
      value<-pvals[i]
      names(value)<-names(pvals)[i]
      sig<-append(sig, value)
    }
    else
    {
      value_2<-pvals[i]
      names(value_2)<-names(pvals)[i]
      non_sig<-append(non_sig, value_2)
    }
  }
  #print(sig)
  #this needs to be optional for first case or take it out of the first 
  
  sig$comp<-names(sig)
  non_sig$comp<-names(non_sig)
  
  out<-list(sig, non_sig)
  names(out)<-c("sig", "non_sig")
  return(out)
}

get_comp<-function(pvals)
{
  pvals$sig$comp<-lapply(pvals$sig$comp, strsplit , "_vs_")
  pvals$sig$comp<-lapply(rapply( pvals$sig$comp, enquote, how="unlist"), eval)
  #print(sig)
  pvals$non_sig$comp<-lapply(pvals$non_sig$comp, strsplit , "_vs_")
  pvals$non_sig$comp<-lapply(rapply(pvals$non_sig$comp, enquote, how="unlist"), eval)
  return(pvals)
}

mean_tests<-function(df, continuous, categorial1, categorial2)
{
  df$merged_cat<-paste(unlist(df[categorial1]), unlist(df[categorial2]), sep="_") #df[categorial2], sep="_"
  df$merged_cat<-as.factor(df$merged_cat)
  full_comb <- combn(levels(df$merged_cat), 2, FUN=list)
  
  
  #
  #Shapiro Wilk Test for Normality Distribution
  #
  shapiro_vec <- map(.x =levels(df$merged_cat),
                     
                     .f = ~ subset(df,merged_cat==.x))
  shapiro_vec_list=list()
  for(i in 1:length(shapiro_vec))
  {
    cp<-shapiro_vec[i]
    shapiro_vec_list[[i]]<-cp[[1]][[continuous]]
  }
  
  shapiro_test<-map(.x=shapiro_vec_list, .f=~shapiro.test(.x))
  
  shapiro_pvals <- map_dbl(.x = shapiro_test, .f  = "p.value")
  names(shapiro_pvals)<-levels(df$merged_cat)
  
  shapiro_pvals<-p.adjust(shapiro_pvals, method="BH")
  norm_eval<-separate_pvals(shapiro_pvals)
  
  norm_comb <- combn(unlist(norm_eval$non_sig$comp), 2, FUN=list) 
  
  #
  #Var Test for Homoscedasticity, when data normally distributed
  #
  var_try<-map(.x= norm_comb,
               .f=~subset(df, merged_cat %in% .x))
  
  var_test <- map(.x = norm_comb, 
                  .f = ~ var.test(formula= as.formula(paste0(continuous,  "~ merged_cat")), 
                                  data=df[which(df$merged_cat %in% unlist(.x)),]))
  
  var_pvals <- map_dbl(.x = var_test, .f  = "p.value")
  
  names(var_pvals) <- map_chr(.x = norm_comb, .f = ~ paste0(.x, collapse = "_vs_"))
  var_pvals<-p.adjust(var_pvals, method="BH")
  
  var_eval<-separate_pvals(var_pvals)
  
  var_eval_2<-get_comp(var_eval) 
  
  
  #
  # t-test comparing means when variance is equal
  #
  t_test <- map(.x = var_eval_2$non_sig$comp, 
                .f = ~ t.test(formula= as.formula(paste0(continuous,  "~ merged_cat")), 
                              var.equal=T,
                              data= subset(df, merged_cat %in% .x)))
  
  ttest_pvals <- map_dbl(.x = t_test, .f  = "p.value")
  names(ttest_pvals) <- map_chr(.x = var_eval_2$non_sig$comp, .f = ~ paste0(.x, collapse = "_vs_"))
  ttest_pvals<-p.adjust(ttest_pvals, method="BH")
  
  ttest_eval<-separate_pvals(ttest_pvals)
  
  var_neq_comb<-setdiff(norm_comb, var_eval_2$non_sig$comp)
  
  #
  # welch test comparing means when variance is equal
  #
  welch_test <- map(.x = var_neq_comb, 
                    .f = ~ t.test(formula= as.formula(paste0(continuous,  "~ merged_cat")), 
                                  var.equal=F,
                                  data= subset(df, merged_cat %in% .x)))
  welch_pvals <- map_dbl(.x = welch_test, .f  = "p.value")
  names(welch_pvals) <- map_chr(.x = var_neq_comb, .f = ~ paste0(.x, collapse = "_vs_"))
  welch_pvals<-p.adjust(welch_pvals, method="BH")
  
  welch_eval<-separate_pvals(welch_pvals)
  welch_eval<-get_comp(welch_eval)
  not_norm_comb<-setdiff(full_comb, norm_comb)
  if(length(not_norm_comb)>0)
  {
    #
    # welch test comparing means when data is not normally distributed
    #
    wilcox_test <- map(.x = not_norm_comb, 
                       .f = ~ wilcox.test(as.formula(paste0(continuous, "~ merged_cat")), exact=F,
                                          data    = subset(subsetted_df, comb %in% .x)))
    wilcox_pvals<-lapply(wilcox_test, function (x) x["p.value"][[1]])#
    print("wilcox done")
    names(wilcox_pvals) <- map_chr(.x = not_norm_comb, .f = ~ paste0(.x, collapse = "_vs_"))
    wilcox_pvals<-p.adjust(wilcox_pvals, method="BH")
    wilcox_eval<-separate_pvals(wilcox_pvals)
    wilcox_eval<-get_comp(wilcox_eval)
  }
  else(wilcox_eval<-list())
  cat("Groups with Normal Distribution and Homoscedasticity, subjected to T-Test and rejecting the Null-Hypothesis of equal means: \n")
  ttest_eval$sig<-ttest_eval$sig[-length(ttest_eval$sig)] #ugly design error, last element is comparison list
  for(i in names(ttest_eval$sig))
  {
    cat(names(ttest_eval$sig[i]), " ", ttest_eval$sig[[i]], " ", stars.pval(ttest_eval$sig[[i]]), "\n") #" ", stars.pval(ttest_eval$sig[[i]]),
  }
  cat("Groups with Normal Distribution without equal Variances, subjected to Welch-Test and rejecting the Null-Hypothesis of equal means: \n")
  welch_eval$sig<-welch_eval$sig[-length(welch_eval$sig)]
  for(i in names(welch_eval$sig))
  {
    cat(names(welch_eval$sig[i]), " ", welch_eval$sig[[i]], " ", stars.pval(welch_eval$sig[[i]]), "\n") #" ", stars.pval(welch_eval$sig[[i]]),
  }
  cat("Groups without Normal Distribution, subjected to Wilcox Rang Sum Test and rejecting the Null-Hypothesis of equal means: \n")
  wilcox_eval$sig<-wilcox_eval$sig[-length(wilcox_eval$sig)]
  for(i in names(wilcox_eval$sig))
  {
    cat(names(wilcox_eval$sig[i]), " ", wilcox_eval$sig[[i]],  " ", stars.pval(wilcox_eval$sig[[i]]), "\n") #" ", stars.pval(wilcox_eval$sig[[i]]),
  }
}
setwd("C:/Users/...")
ex18<-read.csv(file="EX30_7d_ns_hypocotyl.csv")

colnames(ex18)<-c("genotype", "number","area", "perim", "angle", "size", "treatment", "image")
summary(ex18)

mean_tests(ex18, "size", "genotype", "treatment")




