# Group fitting adapted from MODGEE_Singapore.R; k remains disabled.
.compute_neg_LL_detect_param_groups <- function(parms,l_M,l_t_hats)
  {
  l_n =  dim(l_M)[1]
  l_TT = dim(l_M)[2]
  l_ds = c(1.0,exp(parms[1:(l_n-1)]))
  l_hs = exp(parms[(1:l_TT)+l_n-1])
  Ls = lapply(1:l_n,function(i)
    {
    t_max = l_t_hats[i]
    xs_i = l_ds[i]*l_hs[1:t_max]
    L_i = compute_L_i_no_extinction_Cpp(l_M[i,1:t_max],xs_i)
    return(L_i)
    })
  
  LL = sum(log(unlist(Ls)))
  stopifnot(!is.nan(LL))

  return(-LL)
  }

.compute_neg_LL_extinct_param_groups <- function(parms,l_M,l_H,l_ds,l_hs,do_k=T,do_l=T)
  {
  l_n =  dim(l_M)[1]
  l_TT = dim(l_M)[2]

  if ( do_l) l_n_H = dim(l_H)[2]
    
  l_h_T = exp(unlist(parms[1]))
  l_ms = c(0,0,unlist(exp(parms[1+(1:(l_TT-3))])),0)
  if ( do_k ) l_k=parms[1+(l_TT-3)+1]/1 else l_k=0
  if ( do_l ) l_ls=c(0,unlist(parms[1+(l_TT-3)+do_k+(1:(l_n_H-1))])) else l_ls=NULL

  zeros = rep(0,l_TT)
  ones = lapply(1:l_TT,function(i) {x=zeros;x[i]=1;return(x)})
  Ls = lapply(1:l_n,function(i)
    {
    l_ms_prime = l_ms
    if ( do_k ) l_ms_prime = l_ds[i]^l_k*l_ms_prime
    if ( do_l ) l_ms_prime = l_ms_prime*exp(sum(l_H[i,]*l_ls))
    P_i_unconditional = compute_L_i_Cpp(l_M[i,],l_ds[i]*c(l_hs,l_h_T),l_ms_prime)
    P_i_M_zero = compute_L_i_Cpp(zeros,l_ds[i]*c(l_hs,l_h_T),l_ms_prime)
    P_i_M_one = unlist(lapply(ones,function(ones_i) compute_L_i_Cpp(ones_i,l_ds[i]*c(l_hs,l_h_T),l_ms_prime)))
    L_i_conditional = P_i_unconditional/(1-P_i_M_zero-sum(P_i_M_one))
    return(L_i_conditional)
    })
  LL = sum(log(unlist(Ls)))
  stopifnot(!is.nan(LL))

  return(-LL)
  }

.compute_neg_LL_detect_param_1_obs_groups <- function(parm,l_M,l_H,l_ms,l_hs,l_k,l_ls,do_k=T,do_l=T)
  {
  l_TT = length(l_M)

  l_d = exp(parm)

  zeros = rep(0,l_TT)
  ones = lapply(1:l_TT,function(i) {x=zeros;x[i]=1;return(x)})
  l_ms_prime = l_ms
  if ( do_k ) l_ms_prime = l_d^l_k*l_ms_prime
  if ( do_l ) l_ms_prime = l_ms_prime*exp(sum(l_H*l_ls))
  P_i_unconditional = compute_L_i_Cpp(l_M,l_d*l_hs,l_ms_prime)
  P_i_M_one = unlist(lapply(ones,function(ones_i) compute_L_i_Cpp(ones_i,l_d*l_hs,l_ms_prime)))
  LL = .modgee_singleton_loglik(P_i_unconditional, sum(P_i_M_one))

  return(-LL)
  }


.est_modgee_grouped <- function(M,H,ds_start,hs_start,ms_start,k_start=NULL,ls_start=NULL)
  {
  optimizer_runs <- list()
  optimx <- function(...) {
    ans <- tryCatch(optimx::optimx(...), error = function(e) {
      optimizer_runs[[length(optimizer_runs) + 1L]] <<-
        data.frame(convcode = 999L, value = NA_real_, row.names = "BFGS")
      stop(e)
    })
    optimizer_runs[[length(optimizer_runs) + 1L]] <<- ans
    if (any(is.na(ans$convcode) | ans$convcode != 0))
      warning("An optimization did not report normal convergence.", call. = FALSE)
    ans
  }
  stopifnot(is.null(k_start), ncol(H) >= 2L)
  n_spp = dim(M)[1]
  TT = dim(M)[2]

  do_k = !is.null(k_start)
  do_l = !is.null(ls_start)

  if ( do_l) n_H = dim(H)[2]
      
  
  M_prime = t(apply(M,1,function(temp)
    {
    temp2 = 0*temp
    iii = which(temp==1)
    if ( length(iii) ) if ( max(iii)>1 ) temp2[1:(max(iii)-1)] = 1
    return(temp2)
    }))
  
  n_obs = rowSums(M)
  iii = which(n_obs>=2)
  not_iii = which(n_obs<2)
  n_fit = length(iii)
  t_hats = apply(M_prime,1,function(M_prime_i) {x=utils::tail(which(M_prime_i==1),1);ifelse(length(x),x,NA)})

  parms = c(as.numeric(log(ds_start[iii][-1])),as.numeric(log(hs_start[-TT])))
  
  cat('Fitting detection parameters...\n')
  best_parms = try(optimx(parms,
                      .compute_neg_LL_detect_param_groups,l_M=M[iii,-TT,drop=FALSE],l_t_hats=t_hats[iii],
                      method="BFGS"))
  
  if ( 'try-error'%in%class(best_parms) | length(best_parms)==1 | any(is.na(best_parms[1:length(parms)])) )
    {
    cat('Could not fit detection parameters!\n')
    return(NULL)
    }
   
  ds_est = unlist(c(1.0,exp(best_parms[1:(n_fit-1)])))
  hs_est = unlist(c(exp(best_parms[n_fit-1+(1:(TT-1))]),NA))
  
  parms = c(log(hs_start[TT]),log(ms_start[-c(1,2,TT)]))
  if ( do_k ) parms=c(parms,k_start)
  if ( do_l ) parms=c(parms,ls_start[-1])

  cat('Fitting extinction parameters...\n')
  best_parms = try(optimx(parms,
                      .compute_neg_LL_extinct_param_groups,l_M=M[iii,,drop=FALSE],l_H=H[iii,,drop=FALSE],l_ds=ds_est,l_hs=hs_est[-TT],
                      do_k=do_k,do_l=do_l,
                      method="BFGS"))

  if ( 'try-error'%in%class(best_parms) | length(best_parms)==1 | any(is.na(best_parms[1:length(parms)])) )
    {
    cat('Could not fit extinction parameters!\n')
    return(NULL)
    }
  
  hs_est[TT] = unlist(exp(best_parms[1]))
  ms_est = c(0,0,exp(unlist(best_parms[1+(1:(TT-3))])),0)
  if ( do_k ) k_est = best_parms[1+(TT-3)+1] else k_est = NULL
  if ( do_l ) ls_est = c(0,unlist(best_parms[1+(TT-3)+do_k+(1:(n_H-1))])) else ls_est = NULL
  LL = -unlist(best_parms["value"])

  ds_est_final = rep(NA,n_spp)
  ds_est_final[iii] = ds_est

  small_d = 1e-5
  big_d = 100
  
  cat('Fitting detection parameters for species with only one detection...\n')
  if ( length(not_iii) )
    {
    for ( j in 1:length(not_iii) )
      {
      i = not_iii[j]
      n_guess=1000
      d_guesses = exp(seq(log(small_d),log(big_d),length.out=n_guess))
      neg_LL_guesses = rep(NA,n_guess)
      for ( k in 1:n_guess )
        {
        neg_LL_guesses[k] = .modgee_singleton_initial_value(
          .compute_neg_LL_detect_param_1_obs_groups(log(d_guesses[k]),
          l_M=M[i,],l_H=H[i,],l_ms=ms_est,l_hs=hs_est,l_k=k_est,l_ls=ls_est,do_k=do_k,do_l=do_l), species = rownames(M)[i])
        }

      parm = log(d_guesses[which.min(neg_LL_guesses)])
      best_parm = try(optimx(parm,
                          .compute_neg_LL_detect_param_1_obs_groups,l_M=M[i,],l_H=H[i,],
                          l_ms=ms_est,l_hs=hs_est,l_k=k_est,l_ls=ls_est,
                          do_k=do_k,do_l=do_l,
                          method="BFGS"))
    
      if ( 'try-error'%in%class(best_parm) | is.na(best_parm[[1]]) )
        {
        ds_est_final[i] = d_guesses[which.max(-neg_LL_guesses)]
        } else
        {
        ds_est_final[i] = exp(best_parm[[1]])
        }
      }
    
    big_d = 100
    iii_too_big = which(ds_est_final[not_iii]>big_d)
    ds_est_final[not_iii][iii_too_big] = big_d
    }
  
  r_val = list(M=M,H=H,ds_est=ds_est_final,hs_est=hs_est,ms_est=ms_est,k_est=k_est,ls_est=ls_est)
  
  r_val$optimizer_runs <- optimizer_runs
  r_val$LL = LL
  return(r_val);
  }

