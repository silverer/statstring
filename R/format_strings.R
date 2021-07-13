
#' Format a p-value
#'
#' This function takes a p-value (numeric or character) and formats it so that
#' values less than 0.001 are treated as <0.001. Values less than 0.01 are treated as
#' < 0.01
#'
#' @param p_val P-value to be formatted
#' @return A formatted p-value string
#' @export

format_pval_apa <- function(p_val){
  if(is.na(p_val)){
    return("")
  }else{
    tmp = as.numeric(p_val)
    if(tmp < 0.001){
      return("_p_ < .001")
    }else if(tmp < 0.01){
      return("_p_ < .01")
    }else{
      return(paste0("_p_ = ", scales::number(tmp, accuracy = 0.01)))
    }
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
#' @return A formatted t-statistic
#' @export

format_tstat_apa <- function(t_result){
  t_stat = scales::number(t_result$statistic,
                          accuracy = 0.01)
  dof = scales::number(t_result$parameter,
                       big.mark = ",")
  p_val = format_pval_apa(t_result$p.value)

  if(length(t_result$estimate) == 2){
    mdiff = t_result$estimate[1] - t_result$estimate[2]
    if(mdiff < 0.001){
      lci = scales::number(t_result$conf.int[1], accuracy = 0.001)
      uci = scales::number(t_result$conf.int[2], accuracy = 0.001)
      mdiff = scales::number(t_result$estimate[1] - t_result$estimate[2],
                             accuracy = 0.001)
    }else{
      lci = scales::number(t_result$conf.int[1], accuracy = 0.01)
      uci = scales::number(t_result$conf.int[2], accuracy = 0.01)
      mdiff = scales::number(t_result$estimate[1] - t_result$estimate[2],
                             accuracy = 0.01)
    }


  }else{
    mdiff = scales::number(t_result$estimate[1], accuracy = 0.01)
  }
  return_str = paste0("_t_(", dof, ") = ", t_stat,
                      ", ", p_val, ", M~diff~ (95% CI) = ",
                      mdiff, " (", lci, ", ", uci, ")")
  return(return_str)
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
#' @param as_markdown Return a string formatted for markdown output? Default it TRUE
#' @return An APA-formatted string of the ANOVA results
#' @export

aov_wrapper <- function(f_stat, dfn, dfd, p_val,
                        partial_eta=NA,
                        as_markdown = TRUE){
  stat_string <- ""
  if(is.na(partial_eta)){
    stat_string = paste0('_F_(', scales::number(as.numeric(dfn)),
                         ', ',
                         scales::number(as.numeric(dfd)), ')',
                         ' = ', scales::number(as.numeric(f_stat), accuracy = 0.01),
                         ', ', format_pval_apa(as.numeric(p_val)))
  }else{
    stat_string = paste0('_F_(', scales::number(as.numeric(dfn)),
                         ', ',
                         scales::number(as.numeric(dfd)), ')',
                         ' = ',
                         scales::number(as.numeric(f_stat), accuracy = 0.01),
                         ', ', format_pval_apa(as.numeric(p_val)),
                         ', ~partial~ $\\eta^2$ = ',
                         scales::number(as.numeric(partial_eta), accuracy = 0.01))

  }
  if(as_markdown == FALSE){
    stat_string = stringr::str_replace_all(stat_string, "_", "")
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
#' @param predictor The row or variable name to format output for (defaults to last row)
#' @param get.all Return a list of all F-values (TRUE) or just one row (FALSE, default, overrides predictor)
#' @return An APA-formatted string (or strings) of ANOVA results
#' @export

format_anova_string <- function(aov_summary_obj, predictor = NA,
                                get.all = FALSE){
  if (length(aov_summary_obj)==1) {
    #summary(stats::aov()) output
    res = aov_summary_obj[[1]]
    if(get.all==TRUE){
      sstring = mapply(aov_wrapper, f_stat = res$`F value`,
                       dfn = res$Df,
                       dfd = res$Df[nrow(res)],
                       p_val = res$`Pr(>F)`)
      sstring[length(sstring)] = NA #last row is just the error term
    } else if (is.numeric(predictor)){
      sstring = aov_wrapper(f_stat = res$`F value`[predictor],
                            dfn = res$Df[predictor],
                            dfd = res$Df[nrow(res)],
                            p_val = res$`Pr(>F)`[predictor])
    } else if (is.character(predictor)){
      #remove whitespace from rownames
      rownames(res) = unlist(lapply(rownames(res), function(x){x=gsub(" ", "", x); x}))
      tmp = res[predictor,]
      print(tmp)
      sstring = aov_wrapper(f_stat = tmp$`F`[1],
                            dfn = tmp$`Df`[1],
                            dfd = res$`Df`[nrow(res)],
                            p_val = tmp$`Pr(>F)`[1])
    } else{
      sstring = aov_wrapper(f_stat = res$`F value`[nrow(res)-1],
                            dfn = res$Df[nrow(res)-1],
                            dfd = res$Df[nrow(res)],
                            p_val = res$`Pr(>F)`[nrow(res)-1])
    }
  } else if (length(aov_summary_obj)==7) {
    #rstatix output
    if(get.all==TRUE){
      sstring = mapply(aov_wrapper, f_stat = aov_summary_obj$`F`,
                       dfn = aov_summary_obj$`DFn`,
                       dfd = aov_summary_obj$`DFd`,
                       p_val = aov_summary_obj$`p`,
                       partial_eta = aov_summary_obj[,ncol(aov_summary_obj)])
    } else if (is.numeric(predictor)){
      sstring = aov_wrapper(f_stat = aov_summary_obj$`F`[predictor],
                            dfn = aov_summary_obj$`DFn`[predictor],
                            dfd = aov_summary_obj$`DFd`[predictor],
                            p_val = aov_summary_obj$`p`[predictor],
                            partial_eta = aov_summary_obj[predictor,
                                                          ncol(aov_summary_obj)])
    } else if (is.character(predictor)){
      tmp = aov_summary_obj[aov_summary_obj$Effect == predictor,]
      sstring = aov_wrapper(f_stat = tmp$`F`[1],
                            dfn = tmp$`DFn`[1],
                            dfd = tmp$`DFd`[1],
                            p_val = tmp$`p`[1],
                            partial_eta = tmp[1,ncol(tmp)])
    } else{
      sstring = aov_wrapper(f_stat = aov_summary_obj$`F`[nrow(aov_summary_obj)],
                            dfn = aov_summary_obj$`DFn`[nrow(aov_summary_obj)],
                            dfd = aov_summary_obj$`DFd`[nrow(aov_summary_obj)],
                            p_val = aov_summary_obj$`p`[nrow(aov_summary_obj)],
                            partial_eta = aov_summary_obj[nrow(aov_summary_obj),
                                                          ncol(aov_summary_obj)])
    }
  } else if  (length(aov_summary_obj)==4) {
    #APA Tables
    aov_summary_obj = aov_summary_obj$table_body
    if (get.all==TRUE){
      sstring = mapply(aov_wrapper, f_stat = aov_summary_obj$`F`,
                       dfn = aov_summary_obj$`df`,
                       dfd = aov_summary_obj$`df`[nrow(aov_summary_obj)],
                       p_val = aov_summary_obj$`p`,
                       partial_eta = aov_summary_obj$partial_eta2)
    } else if (is.numeric(predictor)){
      sstring = aov_wrapper(f_stat = aov_summary_obj$`F`[predictor],
                            dfn = aov_summary_obj$`df`[predictor],
                            dfd = aov_summary_obj$`df`[nrow(aov_summary_obj)],
                            p_val = aov_summary_obj$`p`[predictor],
                            partial_eta = aov_summary_obj$partial_eta2[predictor])

    } else if (is.character(predictor)){
      tmp = aov_summary_obj[aov_summary_obj$Predictor == predictor,]
      sstring = aov_wrapper(f_stat = tmp$`F`[1],
                            dfn = tmp$`df`[1],
                            dfd = aov_summary_obj$`df`[nrow(aov_summary_obj)],
                            p_val = tmp$`p`[1],
                            partial_eta = tmp$partial_eta2[1])
    } else{
      sstring = aov_wrapper(f_stat = aov_summary_obj$`F`[nrow(aov_summary_obj)-1],
                            dfn = aov_summary_obj$`df`[nrow(aov_summary_obj)-1],
                            dfd = aov_summary_obj$`df`[nrow(aov_summary_obj)],
                            p_val = aov_summary_obj$`p`[nrow(aov_summary_obj)-1],
                            partial_eta = aov_summary_obj$partial_eta2[nrow(aov_summary_obj)-1])
    }
  } else if (length(aov_summary_obj)==13) {
    #stats::aov() output
    aov_summary_obj = summary(aov_summary_obj)
    res = aov_summary_obj[[1]]
    if(get.all==TRUE){
      sstring = mapply(aov_wrapper, f_stat = res$`F value`,
                       dfn = res$Df,
                       dfd = res$Df[nrow(res)],
                       p_val = res$`Pr(>F)`)
      sstring[length(sstring)] = NA #last row is just the error term
    }else if (is.numeric(predictor)){
      sstring = aov_wrapper(f_stat = res$`F value`[predictor],
                            dfn = res$Df[predictor],
                            dfd = res$Df[nrow(res)],
                            p_val = res$`Pr(>F)`[predictor])
    } else if (is.character(predictor)){
      rownames(res) = unlist(lapply(rownames(res), function(x){x=gsub(" ", "", x); x}))
      tmp = res[predictor,]
      sstring = aov_wrapper(f_stat = tmp$`F`[1],
                            dfn = tmp$`Df`[1],
                            dfd = res$`Df`[nrow(res)],
                            p_val = tmp$`Pr(>F)`[1])
    } else{
      sstring = aov_wrapper(f_stat = res$`F value`[nrow(res)-1],
                            dfn = res$Df[nrow(res)-1],
                            dfd = res$Df[nrow(res)],
                            p_val = res$`Pr(>F)`[nrow(res)-1])
    }
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
#' @return Formatted string containing model R-squared
#' @export
format_r2 <- function(mod){
  if(is.numeric(mod)){
    r2 = mod
  }else{
    r2 = extract_r2(mod)
  }
  r2.txt = paste0("_R_^2^ = ",scales::number(r2, accuracy = 0.01))
  return(r2.txt)
}


