#Wrapper function around rvinecopulib (spline-based)
library(kde1d)
library(rvinecopulib)

#object for estimation and methods contour(), plot() and simulate()

dat = cop_data_MIMIC 
variables_of_interest = c("anchor_age", "weight", "Glucose finger stick (range 70-100)",
                          "Hematocrit (serum)", "Hemoglobin") #names(cop_data_MIMIC)
keep_data = FALSE
family_set = "parametric"
cores = 40

variable <- "anchor_age"

estimate_vinecopula_from_data <- function(dat, variables_of_interest = NULL, 
                                          polynomial = FALSE, ID_name = NULL, 
                                          time_name = NULL, keep_data = FALSE, ...) {
  
  #与えられた共変量（covariate）のためのマージナルスプラインを推定
  #関数は、共変量のデータに対してカーネル密度推定（KDE）を実行し、
  #その結果を使用してさまざまな統計関数を提供するリストを生成
  estimate_spline_marginal <- function(covariate, xmin = NaN) { # covariate = dat_out[, variable]
    covariate <- covariate[!is.na(covariate)] # NAの除外
    param <- kde1d(covariate, xmin = xmin) # kde1d 関数を使用して共変量の1次元カーネル密度推定
    param$x <- NULL # 後続の計算には必要ないため，param から $x コンポーネントを削除。
    marg <- list(pdf = function(u) qkde1d(u, param), # 確率密度関数
                 pit = function(x) pkde1d(x, param), # 確率積分変換
                 rdist = function(n) rkde1d(n, param), # ランダム分布, KDEからのランダムサンプルを生成
                 density = function(x) dkde1d(x, param), # 特定の点における密度を計算
                 dist = param) # KDEのパラメータ param 
    return(marg)
  }
  
  
  if(is.null(variables_of_interest)) {
    #変数名を抽出するためのコード
    variables_of_interest <- setdiff(colnames(dat), # 列名を取得
                                     c(ID_name, time_name)) # ただし，IDと時間の列はsetdiffで除外する
  }
  
  if (polynomial) { # 多項式のランダム効果モデルを使うかどうかを指定する
    #多項式のランダム効果モデルを使用する場合に必要な変数が指定されていない時にエラーを発生
    if (is.null(ID_name) | is.null(time_name)) {
      stop("ID_name and time_name should be specified when polynomial = TRUE")
    }
    
    time_range <- range(dat[, time_name])
    
    #estimate individual coefficients from polynomial random effect model
    dat_out <- as.data.frame(matrix(NA, nrow = length(unique(dat[, ID_name])), 
                                    ncol = 3*length(variables_of_interest)))
    names(dat_out) <- paste0("b", 0:2, "_", rep(variables_of_interest, each = 3))
    for (variable in variables_of_interest) {
      formula_v <- as.formula(paste(variable, "~ poly(", time_name, ", 2, raw = TRUE) +", 
                                    "(1 + poly(", time_name,",2, raw = TRUE)|", ID_name,")"))
      re_lm <- lmer(data = dat, formula_v, REML = TRUE)
      #extract individual coefficients
      dat_out[grep(paste0("_", variable), names(dat_out))] <-  coef(re_lm)[[ID_name]]
    }
  } else {# 多項式のランダム効果モデルを使うかどうかを指定していなければ
    dat_out <- dat[, variables_of_interest] # 指定した列名のあるデータを抽出
    time_range <- NULL
  }
  
  #データセットの各変数に対してマージナルスプラインを推定
  marginals <- list()
  dat_unif <- dat_out
  for (variable in names(dat_unif)) { # names(dat_unif)は関心のある変数，確率積分変換を行なった結果を返す
    marginals[[variable]] <- estimate_spline_marginal(dat_out[, variable]) # margとして推定した周辺分布の情報を返す
    ind_na <- is.na(dat_unif[, variable]) # variable 列におけるNA（非数値）の位置を特定
    dat_unif[!ind_na, variable] <- marginals[[variable]]$pit(dat_unif[!ind_na, variable])
    # NAでない要素に対して、marginals[[variable]]$pit 関数を適用し、その結果を dat_unif の対応する列に格納します。pit は、確率積分変換を行う
    # 
  }
  
  #離散変数（discrete variables）の確率積分変換を行う
  #おそらく，var_types = ['c', 'c', 'd', 'c', 'd']みたいに渡してあげることで以下のような処理が可能になる
  if (hasArg(var_types)) {
    arg_list <- list(...) # 最初のfunctionoの「...」の部分
    var_types <- arg_list[["var_types"]] # 「...」の部分にvar_typesの引数があればそれを渡す
    rm(arg_list)
    if (any("d" %in% var_types)) { # var_types が "d" である変数の名前を選択
      discrete_vars <- names(dat_out[var_types == "d"])
      discrete_data <- dat_out[, discrete_vars, drop = F]
      for (variable in discrete_vars) {
        ind_na <- is.na(dat_out[, variable])
        discrete_data[!ind_na, variable] <- marginals[[variable]]$pit(dat_out[!ind_na, variable] - 1)
      }
      names(discrete_data) <- paste0(names(discrete_data), "_d")
      dat_unif <- cbind.data.frame(dat_unif, discrete_data)
    }
  }
  
  
  #estimate copula
  #　https://www.rdocumentation.org/packages/rvinecopulib/versions/0.6.3.1.1/topics/vinecop
  #　5次元のデータだと5分くらいかかる
  vine_coefs <- vinecop(dat_unif) # 元コード：vinecop(dat_unif, ...)
  
  #create output object
  vine_output <- list(vine_copula = vine_coefs, marginals = marginals, 
                      polynomial = polynomial, 
                      names = c(ID_name = ID_name, time_name = time_name), 
                      time_range = time_range, 
                      variables_of_interest = variables_of_interest)
  
  if (keep_data) { # keep_data = TRUEで，データをkeepできる
    vine_output <- append(vine_output, list(original_data = dat_out, uniform_data = dat_unif), after = 3)
  }
  
  class(vine_output) <- "estVineCopula"
  # オブジェクトにクラスを割り当てることで、そのオブジェクトに特定の振る舞いやメソッドを関連付けることができる
  
  return(vine_output)
}

plot.estVineCopula <- function(vine_output, ...) {
  plot(vine_output$vine_copula, ...) # vinecopulaのtree構造を返してくれる
}

contour.estVineCopula <- function(vine_output, ...) {
  contour(vine_output$vine_copula, ...) # 等高線を返してくれる
}

#simulation from estimated vine copula
#value_only = TRUE for longitudinal predictions, or FALSE for list with 
#   parameters and longitudinal predictions
#   nは生成したいsimulationの個数
simulate.estVineCopula <- function(vine_output, n, value_only = TRUE) {
  
  #Simulation
  if (any(vine_output$vine_copula$var_types == "d")) { # 離散変数の場合
    vine_distribution <- vinecop_dist(vine_output$vine_copula$pair_copulas,
                                      vine_output$vine_copula$structure, var_types = vine_output$vine_copula$var_types)
    dat_sim <- as.data.frame(rvinecop(n, vine_distribution))
    names(dat_sim) <- vine_output$variables_of_interest[1:ncol(dat_sim)] # 列名の作成
  } else {
    # 連続変数だけなら
    dat_sim <- as.data.frame(rvinecop(n, vine_output$vine_copula))
  }
  
  for (variable in names(dat_sim)) {
    ind_na <- is.na(dat_sim[, variable])
    dat_sim[!ind_na, variable] <- vine_output$marginals[[variable]]$pdf(dat_sim[!ind_na, variable])
  }
  if (any(vine_output$vine_copula$var_types == "d")) {
    dat_sim[, vine_output$vine_copula$var_types == "d"] <- 
      round(dat_sim[, vine_output$vine_copula$var_types == "d"])
  }
  
  #polynomials 多項式かどうか
  if (vine_output$polynomial) {
    gest_times <- seq(vine_output$time_range[1], vine_output$time_range[2], length.out = 100)
    time_data <- as.data.frame(expand_grid(rownames(dat_sim), gest_times))
    names(time_data) <- vine_output$names
    suppressMessages(df_sim <- dat_sim %>% 
                       rownames_to_column(vine_output$names["ID_name"]) %>%
                       right_join(time_data))
    
    for (variable in vine_output$variables_of_interest) {
      col_ind <- grep(paste0("_", variable), names(df_sim))
      df_sim[, variable] <- df_sim[, col_ind[1]] + df_sim[, col_ind[2]]*df_sim[, vine_output$names["time_name"]] + 
        df_sim[, col_ind[3]]*df_sim[, vine_output$names["time_name"]]^2
    }
    output <- df_sim %>% dplyr::select(all_of(c(as.character(vine_output$names), vine_output$variables_of_interest)))
    if (value_only) {
      return(output)
    } else {
      return(list(values = output, parameters = dat_sim))
    }
  }
  
  return(dat_sim)
}

