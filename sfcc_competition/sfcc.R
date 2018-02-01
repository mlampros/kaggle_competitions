
# code script for the San Francisco Crime Classification competition 

# [ it runs in 2.89 hours (using 6 cores in parallel), test-error : 2.2324, leaderboard-score : 2.238 ]

# refer to the blog-post on how to download the data


train <- read.csv("train.csv", stringsAsFactors = F)

sort(table(train$Category), decreasing = T)

head(train)

test <- read.csv("test.csv", stringsAsFactors = F)
head(test)

train = train[, -c(3,6)]

train$Id = sort(seq(-nrow(train), -1, 1), decreasing = T)

test$Category = rep('none', nrow(test))

train = train[, c(8, 1:7)]

test = test[, c(1:2,8,3:7)]

train = rbind(train, test)

addr = train$Address

library(parallel)

rem_sl = unlist(mclapply(addr, function(x) stringr::str_replace(x, "/", ""), mc.cores = 4))

rem_sl1 = unlist(mclapply(rem_sl, function(x) stringr::str_replace_all(x, pattern=" ", repl=""), mc.cores = 4))

rem_sl2 = as.vector(sapply(rem_sl1, tolower))

train$Address = rem_sl2

date = train$Dates           


library(lubridate)

date1 = ymd_hms(date)

Year = year(date1)

Month = month(date1)

YDay = yday(date1)

WDay = wday(date1)

char_wdays = weekdays(date1)

Day = day(date1)

Hour = hour(date1)

Minute = minute(date1)

remov = data.frame(date_tmp = date1, order_dat = 1:length(date1))

remov1 = remov[order(remov$date_tmp, decreasing = F), ]

remov2 = cbind(remov1, order_out = 1:nrow(remov1))

remov3 = remov2[order(remov2$order_dat, decreasing = F), ]

ORD_rows = remov3$order_out

library(zoo)

yq <- as.yearqtr(as.yearmon(as.Date(train$Dates), "%m/%d/%Y") + 1/12)

Season <- factor(format(yq, "%q"), levels = 1:4, labels = c("winter", "spring", "summer", "fall"))

Season = as.numeric(Season)

newy = which(Month == 1 & Day == 1)

newy1 = which(Month == 12 & Day == 31)

newal = rep(0, length(Month))

newal[newy] = 1

newal[newy1] = 1

train1 = data.frame(year = Year, month = Month, yday = YDay, weekday = WDay, day = Day, hour = Hour, minutes = Minute, season = Season, newy = newal, ORD_date = ORD_rows)

train1 = cbind(train, train1)

library(dplyr)

DISTRICTS = lapply(unique(train1$PdDistrict), function(x) filter(train1, PdDistrict == x))

median_outliers = function(sublist) {
  
  if (max(sublist$X) == -120.5 || max(sublist$Y) == 90.00) {
    
    sublist$X[which(sublist$X == -120.5)] = median(sublist$X)
    sublist$Y[which(sublist$Y == 90.00)] = median(sublist$Y)
  }
  
  sublist
}

distr = lapply(DISTRICTS, function(x) median_outliers(x))

distr1 = do.call(rbind, distr)

address_frequency = function(sublist) {
  
  tmp_df = data.frame(table(sublist$Address))
  
  tmp_df = tmp_df[order(tmp_df$Freq, decreasing = T), ]
  
  tmp_df[1, ]$Var1
}


gcd.hf <- function(long1, lat1, long2, lat2) {                                                                             # http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
  
  R <- 6371                                               # Earth mean radius [km]
  
  delta.long <- (long2 - long1)
  
  delta.lat <- (lat2 - lat1)
  
  a <- sin(delta.lat/2) ^ 2 + cos(lat1) * cos(lat2) * sin(delta.long/2) ^ 2
  
  c <- 2 * asin(min(1,sqrt(a)))
  
  d = R * c
  
  return(d)                                                # Distance in km
}


get_reference_address = function(initial_data, split_column) {                                                               # function to calculate km-distances
  
  s_col = lapply(unique(initial_data[, split_column]), function(x) initial_data[initial_data[, split_column] == x, ])
  
  reference_address = lapply(s_col, function(x) as.character(address_frequency(x)))
  
  reference_lon_lat = lapply(1:length(s_col), function(x) filter(s_col[[x]], Address == reference_address[[x]])[1, c('X','Y')])
  
  Distance = lapply(1:length(s_col), function(f) sapply(1:nrow(s_col[[f]]), function(x) gcd.hf(s_col[[f]][x, 7], s_col[[f]][x, 8], reference_lon_lat[[f]]$X, reference_lon_lat[[f]]$Y)))
  
  tmp_id = do.call(rbind, s_col)$Id
  
  tmp_df = data.frame(id = tmp_id, unlist(Distance))
  
  colnames(tmp_df) = c('Id', paste('dist_', stringr::str_trim(split_column, side = 'both' )))
  
  return(tmp_df)
}


lst_out = list()

for (i in c('PdDistrict', 'weekday', 'day', 'hour', 'season')) {
  
  cat(i, '\n')
  
  lst_out[[i]] = get_reference_address(distr1, i)
}


merg = merge(lst_out[[1]], lst_out[[2]], by.x = 'Id', by.y = 'Id')

merg = merge(merg, lst_out[[3]], by.x = 'Id', by.y = 'Id')

merg = merge(merg, lst_out[[4]], by.x = 'Id', by.y = 'Id')

merg = merge(merg, lst_out[[5]], by.x = 'Id', by.y = 'Id')

ndf = merge(distr1, merg, by.x = 'Id', by.y = 'Id')

ndf$`dist_ weekday` = log(ndf$`dist_ weekday` + 1)

ndf$`dist_ day` = sqrt(ndf$`dist_ day` + 1)

ndf$`dist_ hour` = 2 * sqrt(ndf$`dist_ hour` + 3/8)

ndf = ndf[, -c(2, 4)]

table(ndf$PdDistrict)

pdD = as.factor(ndf$PdDistrict)

mdM = model.matrix(~.-1, data.frame(pdD))

ndf$PdDistrict = NULL

ndf = cbind(ndf, mdM)

ndf$Address = as.numeric(as.factor(ndf$Address))

ntrain = filter(ndf, Id < 0)

ntrain$Id = NULL

response = ntrain$Category

y = c(0:38)[ match(response, sort(unique(response))) ]

ntrain$Category = NULL

ntest = filter(ndf, Id >= 0)

ID_TEST = as.integer(ntest$Id)

ntest$Id = NULL

ntest$Category = NULL

library(Matrix)

ntrain = Matrix(as.matrix(ntrain), sparse = T)

ntest = Matrix(as.matrix(ntest), sparse = T)

VALID_FUNC = function(EVAL_METRIC, arg_actual, arg_predicted, inverse_order = FALSE) {
  
  if (inverse_order == TRUE) {
    
    args_list = list(arg_predicted, arg_actual)
  }
  
  else {
    
    args_list = list(arg_actual, arg_predicted)
  }
  
  result = do.call(EVAL_METRIC, args_list)
  
  result
}




MultiLogLoss = function (y_true, y_pred) {
  
  if (is.factor(y_true)) {
    
    y_true_mat <- matrix(0, nrow = length(y_true), ncol = length(levels(y_true)))
    
    sample_levels <- as.integer(y_true)
    
    for (i in 1:length(y_true)) y_true_mat[i, sample_levels[i]] <- 1
    
    y_true <- y_true_mat
  }
  
  eps <- 1e-15
  
  N <- dim(y_pred)[1]
  
  y_pred <- pmax(pmin(y_pred, 1 - eps), eps)
  
  MultiLogLoss <- (-1/N) * sum(y_true * log(y_pred))
  
  return(MultiLogLoss)
}



xgboost_cv = function(RESP, data, TEST, repeats, Folds, idx_train = NULL, param, num_rounds, print_every_n  = 10, 
                      
                      early_stop = 10, maximize = FALSE, verbose = 1, EVAL_METRIC, set_seed = 2) {
  
  start = Sys.time()
  
  library(caret)
  library(xgboost)
  library(Metrics)
  
  out_ALL = list()
  
  for (j in 1:repeats) {
    
    cat('REPEAT', j, '\n')
    
    TEST_lst = PARAMS = PREDS_tr = PREDS_te = list()
    PARAMS = PREDS_tr = PREDS_te = list()
    
    if (is.numeric(Folds)) {
      
      if (is.null(set_seed)) {
        
        sample_seed = sample(seq(1, 1000000, 1), 1)}
      
      else {
        
        sample_seed = set_seed
      }
      
      set.seed(sample_seed)
      folds = createFolds(RESP, k = Folds, list = TRUE)}
    
    else {
      
      if (is.null(idx_train)) stop(simpleError('give index of train data in form of a vector'))
      
      out_idx = 1:dim(data)[1]
      folds = lapply(1:length(Folds), function(x) out_idx[which(idx_train %in% Folds[[x]])])
    }
    
    tr_er <- tes_er <- rep(NA, length(folds))
    
    for (i in 1:length(folds)) {
      
      cat('fold', i, '\n')
      
      dtrain <- xgb.DMatrix(data = data[unlist(folds[-i]), ], label = RESP[unlist(folds[-i])])
      
      dtest <- xgb.DMatrix(data = data[unlist(folds[i]), ], label = RESP[unlist(folds[i])])
      
      watchlist <- list(train = dtrain, test = dtest)
      
      fit = xgb.train(param, dtrain, nround = num_rounds, print.every.n  = print_every_n, watchlist = watchlist, 
                      
                      early.stop.round = early_stop, maximize = maximize, verbose = verbose)
      
      PARAMS[[i]] = list(param = param, bst_round = fit$bestInd)
      
      pred_tr = predict(fit, data[unlist(folds[-i]), ], ntreelimit = fit$bestInd)
      pred_tr = matrix(pred_tr, nrow =  dim(data[unlist(folds[-i]), ])[1], ncol = length(unique(y)), byrow = TRUE)
      
      pred_te = predict(fit, data[unlist(folds[i]), ], ntreelimit = fit$bestInd)
      pred_te = matrix(pred_te, nrow =  dim(data[unlist(folds[i]), ])[1], ncol = length(unique(y)), byrow = TRUE)
      
      tr_er[i] = VALID_FUNC(EVAL_METRIC, as.factor(RESP[unlist(folds[-i])]), pred_tr)
      tes_er[i] = VALID_FUNC(EVAL_METRIC, as.factor(RESP[unlist(folds[i])]), pred_te)
      
      tmp_TEST = matrix(predict(fit, TEST, ntreelimit = fit$bestInd), nrow =  dim(TEST)[1], ncol = length(unique(y)), byrow = TRUE)
      
      TEST_lst[[paste0('preds_', i)]] = tmp_TEST
      
      cat('---------------------------------------------------------------------------', '\n')
      
      save(tmp_TEST, file = paste('sfcc_', paste(sample(1:1000000000, 1), '_REPEAT_save.RDATA', sep = ""), sep = ""))
      
      gc()
    }
    
    out_ALL[[j]] = list(TEST_lst = TEST_lst, PARAMS = PARAMS, sample_seed = sample_seed, tr_er = tr_er,PREDS_tr = PREDS_tr, PREDS_te = PREDS_te, tes_er = tes_er)
    
    cat('================================================================================================================', '\n')
    
    gc()
  }
  
  end = Sys.time()
  
  return(list(res = out_ALL, time = end - start))
}



params = list("objective" = "multi:softprob", "eval_metric" = "mlogloss", "num_class" = 39, "booster" = "gbtree", "bst:eta" = 0.245, 
              
              "subsample" = 0.7, "max_depth" = 7, "colsample_bytree" = 0.7, "nthread" = 6, "scale_pos_weight" = 0.0, 
              
              "min_child_weight" = 0.0, "max_delta_step" = 1.0) 


fit = xgboost_cv(y, ntrain, ntest, repeats = 1, Folds = 4, idx_train = NULL, params, num_rounds = 145, print_every_n = 5, 
                 
                 early_stop = 10, maximize = FALSE, verbose = 1, MultiLogLoss)




tr_er = unlist(lapply(fit$res, function(x) x$tr_er))

tes_er = unlist(lapply(fit$res, function(x) x$tes_er))

cat('log loss of train is :', mean(tr_er), '\n')

cat('log loss of train is :',  mean(tes_er), '\n')

lap = unlist(lapply(fit$res, function(x) x$TEST_lst), recursive = FALSE)

avg_dfs = (lap[[1]] + lap[[2]] + lap[[3]] + lap[[4]])/4

subms = data.frame(ID_TEST, avg_dfs)

sampleSubmission <- read.csv("~/Desktop/kaggle_gpu/SFCC/sampleSubmission.csv", check.names = F)

colnames(subms) = colnames(sampleSubmission)

subms = subms[order(subms$Id, decreasing = F), ]

write.csv(subms, "xgb_post_submission_train_error_1_95_test_error_2_2324.csv", row.names=FALSE, quote = FALSE)


# plot important variables

ntr = as.matrix(ntrain)

colnames(ntr) = make.names(colnames(ntr))


library(FeatureSelection)

params_xgboost = list(params = list("objective" = "multi:softprob", "eval_metric" = "mlogloss", "num_class" = 39, "booster" = "gbtree", "bst:eta" = 0.5, 
                      
                      "subsample" = 0.65, "max_depth" = 5, "colsample_bytree" = 0.65, "nthread" = 6, "scale_pos_weight" = 0.0, 
                      
                      "min_child_weight" = 0.0, "max_delta_step" = 1.0), 
                      
                      nrounds = 50, print.every.n = 5, verbose = 1, maximize = F)

params_ranger = list(write.forest = TRUE, probability = TRUE, num.threads = 6, num.trees = 50, verbose = FALSE, classification = TRUE, 
                     
                     mtry = 4, min.node.size = 5, importance = 'impurity')

params_features = list(keep_number_feat = NULL, union = F)

feat = wrapper_feat_select(ntr, y, params_glmnet = NULL, params_xgboost = params_xgboost, params_ranger = params_ranger, xgb_sort = NULL,
                           
                           CV_folds = 4, stratified_regr = FALSE, cores_glmnet = 2, params_features = params_features)


params_barplot = list(keep_features = 28, horiz = TRUE, cex.names = 0.8)

barplot_feat_select(feat, params_barplot, xgb_sort = 'Cover')
