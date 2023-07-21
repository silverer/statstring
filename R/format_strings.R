#' Helper function to format p-value
#'
#' This is the base function that format_pval_apa calls.
#'
#' @param p_val P-value to be formatted
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @return A formatted p-value string
#' @export

pval_fmt <- function(p_val, markdown=T){
  if(is.na(p_val)){
    return("")
  }else{
    tmp = as.numeric(p_val)
    if(tmp < 0.001){
      return_str = "_p_ < .001"
    }
    else if(tmp < 0.01){
      tmp = scales::number(tmp, accuracy = 0.001)
      tmp = stringr::str_replace(tmp, "0.", ".")
      return_str = paste0("_p_ = ", tmp)
    }
    else if(scales::number(tmp, accuracy = .01)=="1.00"){
      return_str = paste0("_p_ > .99")
    }
    else{
      tmp = scales::number(tmp, accuracy = 0.01)
      tmp = stringr::str_replace(tmp, "0.", ".")
      return_str = paste0("_p_ = ", tmp)
    }
  }
  if(markdown==F){
    return(stringr::str_remove_all(return_str, "[_]"))
  }else{
    return(return_str)
  }

}


#' Format a p-value
#'
#' This function takes a p-value (numeric or character) or list of p-values and formats it so that
#' values less than 0.001 are treated as <0.001. Values less than 0.01 are treated as
#' < 0.01
#'
#' @param p_val P-value to be formatted
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @return A formatted p-value string
#' @export

format_pval_apa <- function(p_val, as.markdown=T){
  if(length(p_val)>1){
    return(mapply(pval_fmt, p_val=p_val, markdown=as.markdown, USE.NAMES = F))
  }else{
    return(pval_fmt(p_val))
  }
}

#' Get p-value significance stars
#'
#' This function takes a p-value and returns a string containing significance asterisks
#' Returns nothing if p > 0.05
#'
#' @param p_val P-value to get stars for
#' @return A string of significance asterisks
#' @export

format_sig_stars <- function(p_val){
  tmp = as.numeric(p_val)
  if(tmp < 0.001){
    return("***")
  }else if(tmp < 0.01){
    return("**")
  }else if(tmp < 0.05){
    return("*")
  }else{
    return("")
  }
}


#' Format results of a t-test
#'
#' This function takes the results of a t-test and formats it for
#' easy markdown outputs
#'
#' @param t_result t-test result to be formatted
#' @param cohen_d optional, Cohen's d standardized effect size
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @return A formatted t-statistic
#' @export

format_tstat_apa <- function(t_result, cohen_d="none", as.markdown=T){
  t_stat = scales::number(t_result$statistic,
                          accuracy = 0.01)
  dof = scales::number(t_result$parameter,
                       big.mark = ",")
  p_val = format_pval_apa(t_result$p.value)
  if(length(t_result$estimate) == 2){
    mdiff = t_result$estimate[1] - t_result$estimate[2]

  }
  else{
    mdiff = t_result$estimate
  }
  if(mdiff < 0){
    comparison_dir = "negative"
  }
  else{
    comparison_dir = "positive"
  }
  if(scales::number(mdiff, accuracy = .01)=="0.00"){
    lci = scales::number(t_result$conf.int[1], accuracy = 0.001)
    uci = scales::number(t_result$conf.int[2], accuracy = 0.001)
    mdiff = scales::number(mdiff,
                           accuracy = 0.001)
  }else{
    lci = scales::number(t_result$conf.int[1], accuracy = 0.01)
    uci = scales::number(t_result$conf.int[2], accuracy = 0.01)
    mdiff = scales::number(mdiff,
                           accuracy = 0.01)
  }

  return_str = paste0("_t_(", dof, ") = ", t_stat,
                      ", ", p_val, ", M~diff~ (95% CI) = ",
                      mdiff, " (", lci, ", ", uci, ")")
  if(as.markdown==F){
    return_str = stringr::str_remove_all(return_str, "[_]")
    return_str = stringr::str_remove_all(return_str, "[~]")
  }
  if(is.character(cohen_d)){
    return(return_str)
  }
  else{
    #make sure cohen's d is the same direction as the t-statistic
    if(((comparison_dir == "positive" & cohen_d < 0)|
       (comparison_dir == "negative" & cohen_d > 0))){
      cohen_d = -1*cohen_d
    }
    cohen_d = scales::number(cohen_d, accuracy = 0.01)
    return_str = paste0(return_str, ", _d_ = ", cohen_d)
    if(as.markdown==F){
      return_str = stringr::str_remove_all(return_str, "[_]")
      return_str = stringr::str_remove_all(return_str, "[~]")
    }
    return(return_str)
  }

}

#' Format results of a pairwise comparison
#'
#' This function takes the results of a pairwise comparison and formats it for
#' easy markdown outputs
#'
#' @param t_stat t-statistic
#' @param df degrees of freedom
#' @param p_val p-value of t-statistic
#' @param mdiff optional, the mean difference of the comparison
#' @param lci optional, the lower end of the confidence interval on the comparison
#' @param uci optional, the upper end of the confidence interval on the comparison
#' @param cohen_d optional, Cohen's d standardized effect size
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @return A formatted t-statistic
#' @export

format_pairwise_comparison <- function(t_stat, df, p_val,
                                       mdiff = "none", lci = "none", uci = "none",
                                       cohen_d="none",
                                       as.markdown=T){
  t_stat = scales::number(t_stat,
                          accuracy = 0.01)
  dof = scales::number(df,
                       big.mark = ",")
  p_val = format_pval_apa(p_val)
  return_str = paste0("_t_(", dof, ") = ", t_stat,
                      ", ", p_val)
  if(is.numeric(mdiff) & is.numeric(lci) & is.numeric(uci)){
    if(mdiff < 0){
      comparison_dir = "negative"
    }
    else{
      comparison_dir = "positive"
    }
    if(scales::number(mdiff, accuracy = .01)=="0.00"){
      lci = scales::number(lci, accuracy = 0.001)
      uci = scales::number(uci, accuracy = 0.001)
      mdiff = scales::number(mdiff,
                             accuracy = 0.001)
    }else{
      lci = scales::number(lci, accuracy = 0.01)
      uci = scales::number(uci, accuracy = 0.01)
      mdiff = scales::number(mdiff,
                             accuracy = 0.01)
    }
  }
  return_str = paste0(return_str, ", M~diff~ (95% CI) = ",
                      mdiff, " (", lci, ", ", uci, ")")
  if(as.markdown==F){
    return_str = stringr::str_remove_all(return_str, "[_]")
    return_str = stringr::str_remove_all(return_str, "[~]")
  }
  if(is.character(cohen_d)){
    return(return_str)
  }
  else{
    #make sure cohen's d is the same direction as the t-statistic
    if(((comparison_dir == "positive" & cohen_d < 0)|
        (comparison_dir == "negative" & cohen_d > 0))){
      cohen_d = -1*cohen_d
    }
    cohen_d = scales::number(cohen_d, accuracy = 0.01)
    return_str = paste0(return_str, ", _d_ = ", cohen_d)
    if(as.markdown==F){
      return_str = stringr::str_remove_all(return_str, "[_]")
      return_str = stringr::str_remove_all(return_str, "[~]")
    }
    return(return_str)
  }

}

#' Format results of a bivariate correlation
#'
#' This function takes the results of a bivariate correlation and formats it for
#' easy markdown outputs
#'
#' @param r value of correlation to be formatted
#' @param dof degrees of freedom for the correlation
#' @param p p-value of correlation
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @return A formatted bivariate correlation statistic
#' @export

format_corr_apa <- function(r, dof, p, as.markdown=T){
  r = scales::number(r,accuracy = 0.01)
  r = stringr::str_replace(r, "0.", ".")
  dof = scales::number(dof,
                       big.mark = ",")
  p_val = format_pval_apa(p)

  return_str = paste0("_r_(", dof, ") = ", r,
                      ", ", p_val)
  if(as.markdown){
    return(return_str)
  }else{
    return_str = stringr::str_remove_all(return_str, "[_]")
    return(return_str)
  }

}

#' Format ANOVA results for APA reporting
#'
#' This function takes the outputs from an ANOVA object (the F-ratio, the degrees of freedom,
#' the p-value, and, optionally, the partial eta squared) and formats them as a string for
#' easy copying and pasting into a document or piping into a markdown script.
#' To return a string that is not formatted for markdown, set as.markdown to FALSE.
#' All arguments can be passed either as strings or as numbers.
#'
#' @param f_stat The F ratio from an ANOVA
#' @param dfn The degrees of freedom (numerator) from an ANOVA
#' @param dfd The degrees of freedom (denominator) from an ANOVA
#' @param p_val The p-value from an ANOVA
#' @param partial_eta Optional, the effect size from an ANOVA. Default is NA
#' @param as.markdown Return a string formatted for markdown output? Default is TRUE
#' @return An APA-formatted string of the ANOVA results
#' @export

aov_wrapper <- function(f_stat, dfn, dfd, p_val,
                        partial_eta=NA,
                        as.markdown = TRUE){
  stat_string <- ""
  if(is.na(partial_eta)|partial_eta=="NA"|partial_eta==""){
    stat_string = paste0('_F_(', scales::number(as.numeric(dfn), big.mark = ","),
                         ', ',
                         scales::number(as.numeric(dfd), big.mark = ","),
                         ')',
                         ' = ',
                         scales::number(as.numeric(f_stat), accuracy = 0.01),
                         ', ',
                         format_pval_apa(as.numeric(p_val)))
  }else{
    if(scales::number(as.numeric(partial_eta), accuracy = 0.01)=="0.00"){
      partial_eta_str = scales::number(as.numeric(partial_eta), accuracy = 0.001)
    }else{
      partial_eta_str = scales::number(as.numeric(partial_eta), accuracy = 0.01)
    }
    stat_string = paste0('_F_(',
                         scales::number(as.numeric(dfn), big.mark = ","),
                         ', ',
                         scales::number(as.numeric(dfd), big.mark = ","), ')',
                         ' = ',
                         scales::number(as.numeric(f_stat), accuracy = 0.01),
                         ', ', format_pval_apa(as.numeric(p_val)),
                         ', ~partial~ $\\eta^2$ = ',
                         partial_eta_str)

  }
  if(as.markdown == FALSE){
    stat_string = stringr::str_remove_all(stat_string, "[_]")
    stat_string = stringr::str_replace(stat_string, "~partial~ $\\eta^2$", "partial eta sq.")
  }
  return(stat_string)
}

#' Extract statistics and format raw ANOVA results for APA reporting in markdown
#'
#' This function takes the summary object from an ANOVA object and feeds it
#' to the format_anova_string function.
#' Accepted ANOVA objects include:
#' stats::aov()
#' summary(stats::aov())
#' rstatix::get_anova_table()
#' apaTables::apa.aov.table()
#'
#'
#' @param aov_summary_obj The ANOVA summary object
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @param get.all deprecated, no longer an option but retaining in the function prevents errors from prior versions
#' @return An APA-formatted string (or strings) of ANOVA results
#' @export

format_anova_string <- function(aov_summary_obj, as.markdown=T, get.all=T){
  if (length(aov_summary_obj)==1) {
    #summary(stats::aov()) output
    res = aov_summary_obj[[1]]
    sstring = mapply(aov_wrapper, f_stat = res$`F value`,
                     dfn = res$Df,
                     dfd = res$Df[nrow(res)],
                     p_val = res$`Pr(>F)`,
                     as.markdown=as.markdown)
    sstring[length(sstring)] = NA #last row is just the error term

  } else if (length(aov_summary_obj)==7) {
    #rstatix output
    sstring = mapply(aov_wrapper, f_stat = aov_summary_obj$`F`,
                     dfn = aov_summary_obj$`DFn`,
                     dfd = aov_summary_obj$`DFd`,
                     p_val = aov_summary_obj$`p`,
                     partial_eta = aov_summary_obj[,ncol(aov_summary_obj)],
                     as.markdown=as.markdown)
  } else if  (length(aov_summary_obj)==4) {
    #APA Tables output
    aov_summary_obj = aov_summary_obj$table_body
    sstring = mapply(aov_wrapper, f_stat = aov_summary_obj$`F`,
                     dfn = aov_summary_obj$`df`,
                     dfd = aov_summary_obj$`df`[nrow(aov_summary_obj)],
                     p_val = aov_summary_obj$`p`,
                     partial_eta = aov_summary_obj$partial_eta2,
                     as.markdown=as.markdown)

  } else if (length(aov_summary_obj)==13) {
    #stats::aov() output
    aov_summary_obj = summary(aov_summary_obj)
    res = aov_summary_obj[[1]]
    sstring = mapply(aov_wrapper, f_stat = res$`F value`,
                     dfn = res$Df,
                     dfd = res$Df[nrow(res)],
                     p_val = res$`Pr(>F)`,
                     as.markdown=as.markdown)
    sstring[length(sstring)] = NA #last row is just the error term
  } else {
    sstring = ""
  }
  return(sstring)
}

#' Wrapper to format raw point estimate and lower and upper CI
#'
#' This function takes a point estimate, lower CI, and upper CI and returns
#' a formatted string with the results
#'
#' @param num A point estimate
#' @param lower.ci Lower bound of a confidence interval
#' @param upper.ci Upper bound of a confidence interval
#' @return An APA-formatted string of CI results
#' @export
format_confint <- function(num, lower.ci, upper.ci){
  num = scales::number(num, accuracy = 0.01)
  lower.ci = scales::number(lower.ci, accuracy = 0.01)
  upper.ci = scales::number(upper.ci, accuracy = 0.01)

  return(paste0(num, " [", lower.ci, ", ", upper.ci, "]"))
}

#' Wrapper to extract R-squared from a linear regression model
#'
#' This function takes an object produced by stats::lm() and returns
#' the model R-squared
#'
#' @param mod A linear regression model or a summary(lm()) object
#' @return Model R-squared
#' @export
extract_r2 <- function(mod){
  if(length(mod) == 13){
    #Check if summary has already been called
    mod.sum = stats::summary.lm(mod)
  }else{
    mod.sum = mod
  }
  return(mod.sum$r.squared)
}


#' Wrapper to extract R-squared from a linear regression model
#'
#' This function takes an object produced by stats::lm() and returns
#' the model R-squared as a string formaatted for markdown output
#'
#' @param mod A linear regression model or the R-squared from a regression model
#' @param as.markdown optional, indicates whether the resulting string should be formatted for markdown. Default is TRUE
#' @return Formatted string containing model R-squared
#' @export
format_r2 <- function(mod,as.markdown=T){
  if(is.numeric(mod)){
    r2 = mod
  }else{
    r2 = extract_r2(mod)
  }
  if(as.markdown){
    r2.txt = paste0("_R_^2^ = ",scales::number(r2, accuracy = 0.01))
    return(r2.txt)
  }else{
    return(paste0("R2 = ",scales::number(r2, accuracy = 0.01)))
  }

}

