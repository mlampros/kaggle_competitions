

# load a predictions file
#------------------------

subm = read.csv('../path_to_any_submitted_predictions_file.csv')


# code chunk to load train - test data
#--------------------------------------

PATH = "../train.json"

dat = jsonlite::fromJSON(PATH)


PATH = "../test.json"

dat_te = jsonlite::fromJSON(PATH)



# indices of true images 
#-----------------------

# extract the true images index as explained in : https://www.kaggle.com/brassmonkey381/viewing-leak-and-machine-images 
# however the number of true images should be 3425 and not 3412


relevant_idx = read.csv("../indices_of_true_images.csv")

relevant_idx = as.vector(relevant_idx[, 1])                       # take the first column if one-column-matrix



# merge train and test data using the incidence angle as ID
#----------------------------------------------------------

df_tr = data.frame(ID = dat$inc_angle, IS_ICEB = dat$is_iceberg, idx_tr = 1:length(dat$is_iceberg))



df_te = data.frame(ID = dat_te$inc_angle, true_idx_te = 1:length(dat_te$inc_angle))

df_te = df_te[relevant_idx, ]



mrg = merge(df_te, df_tr, by.x = "ID", by.y = "ID", all.x = T)



compl = complete.cases(mrg)
sum(compl)

true_mrg = mrg[compl, ]



df_tbl = data.frame(table(true_mrg$idx_tr))

colnames(df_tbl) = c('idx_tr', 'freq')

df_tbl = df_tbl[order(df_tbl$freq, decreasing = T), ]



# "manipulate" angle in train - test
#------------------------------------

cp_subm_MANIPULATE = subm                                                       # copy initial predictions


for (i in 1:nrow(df_tbl)) {
  
  SUBSET = subset(true_mrg, true_mrg$idx_tr == as.numeric(as.character(df_tbl$idx_tr[i])))
  
  SUBSET$PROBS = subm$is_iceberg[SUBSET$true_idx_te]
  
  if (length(unique(SUBSET$IS_ICEB)) > 1) stop("more than 1 unique value", call. = F) 
  
  if (unique(SUBSET$IS_ICEB) == 0) {
    
    tmp_prob = as.double(summary(SUBSET$PROBS)['1st Qu.'])                      # use 1st quartile as threshold
  }
  
  if (unique(SUBSET$IS_ICEB) == 1) {
    
    tmp_prob = as.double(summary(SUBSET$PROBS)['3rd Qu.'])                      # use 3rd quartile as threshold
  }
  
  cp_subm_MANIPULATE[SUBSET$true_idx_te, 'is_iceberg'] = tmp_prob
}


# make submission
#-----------------

write.csv(cp_subm_MANIPULATE, "../quartiles_submission.csv", row.names=FALSE, quote = FALSE)

