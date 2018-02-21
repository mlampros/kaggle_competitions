
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and build features
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# volume and degrees to radians   [ both appeared in one of the kaggle-kernels ]
#------------------------------

volume_structure = function(lattice_vector_1_ang, lattice_vector_2_ang, lattice_vector_3_ang,
                            radians_lattice_angle_alpha_degree, radians_lattice_angle_beta_degree, radians_lattice_angle_gamma_degree) {

  tmp_calc = lattice_vector_1_ang * lattice_vector_2_ang * lattice_vector_3_ang * sqrt(1.0 + 2.0 * cos(radians_lattice_angle_alpha_degree) * cos(radians_lattice_angle_beta_degree) * cos(radians_lattice_angle_gamma_degree)
                                                                                       - (cos(radians_lattice_angle_alpha_degree)^2) - (cos(radians_lattice_angle_beta_degree)^2) - (cos(radians_lattice_angle_gamma_degree)^2))

  return (tmp_calc)
}


degrees_2radians = function(lattice_vector_alpha_beta_gamma_degree) {

  tmp_out = pi * (lattice_vector_alpha_beta_gamma_degree / 180.0)

  return(tmp_out)
}


# normalize x, y, z
#------------------

normalized = function(x) {

  out = (x - min(x))/(max(x) - min(x))

  out = out / (sum(out))

  out
}



# build features for train-test data
#-----------------------------------

DAT = c('train', 'test')

for (DAT_ITER in DAT) {

  tr_te = DAT_ITER

  FOLDER_PATH = paste(c(".../nomad2018/", tr_te, "/"), collapse = "")                             # adjust to the folder where train and test are saved

  train_lst = list.files(FOLDER_PATH)

  #--------------------------------


  NUM_CORES = 5         # number of cores to use in parallel


  # function to use both for train and test
  #-----------------------------------------


  tbl_last_col = do.call(rbind, parallel::mclapply(1:length(train_lst), function(ITEM) {

    dat = read.delim( paste(c(FOLDER_PATH, train_lst[ITEM], "/geometry.xyz"), collapse = ""), header = F, skip = 3, na.strings = c(""), sep = " ")


    #---------------------------------------------------------------------------------------

    # unlist 'lattice-vector' and 'pca-components' for each category in atoms

    subs_lattice = t(as.matrix(as.vector(as.matrix(subset(dat, dat$V1 == 'lattice_vector')[, -c(1, ncol(dat))]))))

    colnames(subs_lattice) = paste("subs_lattice_vec_", 1:ncol(subs_lattice), sep = "")

    split_nams = c("Al", "Ga", "In", "O")

    summary_lst = list()

    for (x in 1:length(split_nams)) {

      tmp1 = subset(dat, dat$V1 == 'atom')[, 2:5]

      tmp_mt1 = subset(tmp1, tmp1$V5 == split_nams[x])

      if (nrow(tmp_mt1) == 0) {

        summary_lst[[x]] = rep(0, 18)}

      else {

        summary_lst[[x]] = as.vector(apply(tmp_mt1[, 1:3], 2, summary))
      }
    }

    df_summary = matrix(unlist(summary_lst), nrow = 1, ncol = 72)

    colnames(df_summary) = unlist(lapply(c("Al_", "Ga_", "In_", "O_"), function(x) paste(x, paste(c('MIN_', '1stQUART_', 'MEDIAN_', 'MEAN_', '3rdQUART_', 'MAX_'), c(rep('x', 6), rep('y', 6), rep('z', 6)), sep = ""), sep = "")))


    subs_pca_atom = unlist(lapply(1:length(split_nams), function(x) {

      tmp = subset(dat, dat$V1 == 'atom')[, 2:5]

      tmp_mt = subset(tmp, tmp$V5 == split_nams[x])

      if (nrow(tmp_mt) == 0) {

        out = as.matrix(rep(0.0, 6))}

      else {

        tmp_mt = t(tmp_mt[, -ncol(tmp_mt)])

        out = tryCatch({

          data_pca = prcomp(tmp_mt, scale = TRUE)

          as.matrix(as.vector(data_pca$x[, 1:2]))

        },

        error=function(cond) {

          as.matrix(rep(0.0, 6))
        }
        )
      }

      out
    }))

    subs_pca_atom = t(as.matrix(subs_pca_atom))

    colnames(subs_pca_atom) = unlist(lapply(c("Al_", "Ga_", "In_", "O_"), function(x) paste(x, 1:6, sep = "")))

    #---------------------------------------------------------------------------------------


    f0 = data.frame(matrix(as.numeric(train_lst[ITEM]), nrow = 1, ncol = 1))

    colnames(f0) = "idx"

    f1 = data.frame(matrix(nrow(dat), nrow = 1, ncol = 1))

    colnames(f1) = "nrows"

    lattice = dat[1:3, 2:4]

    f2 = data.frame(matrix(colSums(lattice) / sum(colSums(lattice)), nrow = 1, ncol = 3))

    colnames(f2) = paste("lattice", 1:3, sep = "_")

    x_y_z = dat[4:nrow(dat), 2:4]

    x_y_z_norm = do.call(rbind, lapply(1:nrow(x_y_z), function(x) normalized(x_y_z[x, ])))

    f3 = data.frame(matrix(colSums(x_y_z) / sum(colSums(x_y_z)), nrow = 1, ncol = 3))

    colnames(f3) = paste("atom", 1:3, sep = "_")

    f4 = data.frame(matrix(0, nrow = 1, ncol = 4))

    colnames(f4) = c("Al", "Ga", "In", "O")

    tbl = table(dat$V5)

    nam_tbl = names(tbl)

    val_tbl = as.vector(tbl)

    val_tbl = val_tbl / sum(val_tbl)

    for (i in 1:length(nam_tbl)) {

      f4[, nam_tbl[i]] = val_tbl[i]
    }

    al_ga_in = as.vector(dat[4:nrow(dat), 5])

    accum_value = 0                                    # accumulate the values for Al, Ga, In

    for (j in 1:length(al_ga_in)) {

      if (al_ga_in[j] == "Al") {

        accum_value = accum_value + x_y_z[j, 1]
      }

      if (al_ga_in[j] == "Ga") {

        accum_value = accum_value + x_y_z[j, 2]
      }

      if (al_ga_in[j] == "In") {

        accum_value = accum_value + x_y_z[j, 3]
      }
    }

    f5 = data.frame(matrix(accum_value, nrow = 1, ncol = 1))

    colnames(f5) = "simple_feature"

    cbind(subs_lattice, subs_pca_atom, f0, f1, f2, f3, f4, f5, df_summary)
  }, mc.cores = NUM_CORES))


  # exclude the 'O' column     [ it consists of 0.6 or 3 * N whereas Al, Ga, In consist of 0.4 or 2 * N ]

  tbl_last_col = tbl_last_col[, -which(colnames(tbl_last_col) == 'O')]

  tbl_last_col = tbl_last_col[order(tbl_last_col$idx, decreasing = F), ]

  tbl_last_col = tbl_last_col[, -which(colnames(tbl_last_col) %in% c('idx', 'nrows', "Al", "Ga", "In"))]


  # cbind initial data
  #--------------------

  all_dat <- read.csv(paste(c(".../nomad2018/", tr_te, ".csv"), collapse = ""), header = T)                     # adjust to the folder where train and test are saved


  # calculate global features
  #--------------------------

  newd_glob = all_dat[, c('percent_atom_al', 'percent_atom_ga')] * all_dat[, 'number_of_total_atoms']

  sm = colSums(newd_glob)

  newd_glob$percent_atom_al = newd_glob$percent_atom_al / as.vector(sm)[1]
  newd_glob$percent_atom_ga = newd_glob$percent_atom_ga / as.vector(sm)[2]

  colnames(newd_glob) = paste("global_atom_al_ga_", 1:2, sep = "")


  # cbind data
  #-----------

  all_dat = cbind(all_dat, newd_glob, tbl_last_col)


  # calculate volume of data
  #-------------------------

  res_volume = volume_structure(all_dat$lattice_vector_1_ang, all_dat$lattice_vector_2_ang, all_dat$lattice_vector_3_ang,
                                degrees_2radians(all_dat$lattice_angle_alpha_degree), degrees_2radians(all_dat$lattice_angle_beta_degree),
                                degrees_2radians(all_dat$lattice_angle_gamma_degree))

  all_dat$volume_structure = res_volume


  # atomic density
  #---------------

  all_dat$atomic_density = all_dat$number_of_total_atoms / all_dat$volume_structure


  if (tr_te == "train") {

    y1 = all_dat$formation_energy_ev_natom

    y2 = all_dat$bandgap_energy_ev


    all_dat$formation_energy_ev_natom = NULL

    all_dat$bandgap_energy_ev = NULL

    all_dat$id = NULL

    ntrain = all_dat
  }

  if (tr_te == "test") {

    IDX_TEST = all_dat$id

    all_dat$id = NULL

    ntest = all_dat
  }
}

print(dim(ntrain))
print(dim(ntest))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-validation and shuffling function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# shuffling function for cross-validation folds
#-----------------------------------------------


func_shuffle = function(vec, times = 10) {

  for (i in 1:times) {

    out = sample(vec, length(vec))
  }
  out
}


# regression cross-validation folds
#-----------------------------------


regr_folds = function(folds, RESP, stratified = FALSE) {

  if (is.factor(RESP)) {

    stop(simpleError("this function is meant for regression for classification use the 'class_folds' function"))
  }

  samp_vec = rep(1/folds, folds)

  sort_names = paste0('fold_', 1:folds)


  if (stratified == TRUE) {

    stratif = cut(RESP, breaks = folds)

    clas = lapply(unique(stratif), function(x) which(stratif == x))

    len = lapply(clas, function(x) length(x))

    prop = lapply(len, function(y) sapply(1:length(samp_vec), function(x) round(y * samp_vec[x])))

    repl = unlist(lapply(prop, function(x) sapply(1:length(x), function(y) rep(paste0('fold_', y), x[y]))))

    spl = suppressWarnings(split(1:length(RESP), repl))}

  else {

    prop = lapply(length(RESP), function(y) sapply(1:length(samp_vec), function(x) round(y * samp_vec[x])))

    repl = func_shuffle(unlist(lapply(prop, function(x) sapply(1:length(x), function(y) rep(paste0('fold_', y), x[y])))))

    spl = suppressWarnings(split(1:length(RESP), repl))
  }

  spl = spl[sort_names]

  if (length(table(unlist(lapply(spl, function(x) length(x))))) > 1) {

    warning('the folds are not equally split')            # the warning appears when I divide the unique labels to the number of folds and instead of an ingeger I get a float
  }

  if (length(unlist(spl)) != length(RESP)) {

    stop(simpleError("the length of the splits are not equal with the length of the response"))
  }

  spl
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# regularized greedy forest (5-fold cross-validation repeated 6 times each time with different seed)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

START = Sys.time()

# sample different seeds
#-----------------------

SEEDS = c(473L, 691L, 783L, 799L, 984L, 999L)


# parameters for regularized greedy forest
#-----------------------------------------

# parameters RGF
#---------------

params_rgf_y1 = list(LEAF = c(1500, 1500, 1500, 1500, 1000, 1000), ALGORITHM = rep("RGF_Sib", 6), RATE = c(0.3, 0.3, 0.3, 0.3, 0.1, 0.1), INTERVALL = c(1000, 1500, 1000, 1000, 1000, 1500),
                     NORMALIZE = rep(F, 6), NITER = c(30, 30, 10, 30, 10, 10), SAMPLES = c(50, 50, 50, 10, 10, 10), OPT_INTERVALL = c(100, 100, 50, 100, 100, 100))

params_rgf_y2 = list(LEAF = c(1500, 1500, 1000, 1000, 1500, 1500), ALGORITHM = c("RGF_Opt", "RGF_Opt", "RGF_Opt", "RGF_Sib", "RGF_Opt", "RGF_Opt"),
                     RATE = rep(0.1, 6), INTERVALL = c(1000, 1500, 1000, 1000, 1000, 1000), NORMALIZE = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
                     NITER = c(10, 10, 10, 30, 10, 30), SAMPLES = c(10, 10, 10, 10, 50, 20), OPT_INTERVALL = c(100, 100, 100, 100, 50, 100))


# set the number of folds
#-------------------------

NUM_FOLDS = 5


library(RGF)


# regularized greedy forest       [ the function requires also the 'MLmetrics' library ]
#--------------------------

FUNC_RGF = function(x, y, MAX_LEAF, ALGORITHM, RATE, INTERVALL, NORMALIZE, NITER, SAMPLES, OPT_INTERVALL, TEST_DATA, folds, THREADS = 5, SEED = 1) {

  lst_preds = parallel::mclapply(1:length(folds), function(i) {

    rgf = RGF_Regressor$new(max_leaf = MAX_LEAF,
                            algorithm = ALGORITHM,
                            learning_rate = RATE,
                            test_interval = INTERVALL,
                            normalize = NORMALIZE,
                            loss = "LS",
                            n_iter = NITER,
                            min_samples_leaf = SAMPLES,
                            opt_interval = OPT_INTERVALL,
                            verbose = 0)

    rgf$fit(as.matrix(x[unlist(folds[-i]), ]), log(y[unlist(folds[-i])] + 1))

    pr_tr = exp(rgf$predict(as.matrix(x[unlist(folds[-i]), ]))) - 1

    pr_te = exp(rgf$predict(as.matrix(x[unlist(folds[i]), ]))) - 1

    preds_TEST = exp(rgf$predict(as.matrix(TEST_DATA))) - 1

    er_tr = MLmetrics::RMSLE(y[unlist(folds[-i])], pr_tr)

    er_te = MLmetrics::RMSLE(y[unlist(folds[i])], pr_te)


    list(er_tr = er_tr, er_te = er_te, MAX_LEAF = MAX_LEAF, ALGORITHM = ALGORITHM,

         RATE = RATE, INTERVALL = INTERVALL, NORMALIZE = NORMALIZE, NITER = NITER, SAMPLES = SAMPLES, OPT_INTERVALL = OPT_INTERVALL, pr_te, preds_TEST)

  }, mc.cores = THREADS)

  return(lst_preds)
}


# loop over parameter-settings and save predictions  [ here is necessary that each sublist of the params is equal in length for both y1 and y2 ]
#--------------------------------------------------

res_RGF_y1 = res_RGF_y2 = list()


for (ITER in 1:length(params_rgf_y1[['LEAF']])) {

  cat('iteration : ', ITER, '\n')

  set.seed(SEEDS[ITER])
  FOLDS = regr_folds(folds = NUM_FOLDS, y1, stratified = F)           # stratified FALSE because I split the data using y1 and I use the folds for both y1 and y2

  res_RGF_y1[[ITER]] = FUNC_RGF(ntrain, y1, params_rgf_y1[['LEAF']][ITER], params_rgf_y1[['ALGORITHM']][ITER], params_rgf_y1[['RATE']][ITER],
                                params_rgf_y1[['INTERVALL']][ITER], params_rgf_y1[['NORMALIZE']][ITER], params_rgf_y1[['NITER']][ITER], params_rgf_y1[['SAMPLES']][ITER],
                                params_rgf_y1[['OPT_INTERVALL']][ITER], TEST = ntest, folds = FOLDS, THREADS = 5)

  res_RGF_y2[[ITER]] = FUNC_RGF(ntrain, y2, params_rgf_y2[['LEAF']][ITER], params_rgf_y2[['ALGORITHM']][ITER], params_rgf_y2[['RATE']][ITER],
                                params_rgf_y2[['INTERVALL']][ITER], params_rgf_y2[['NORMALIZE']][ITER], params_rgf_y2[['NITER']][ITER], params_rgf_y2[['SAMPLES']][ITER],
                                params_rgf_y2[['OPT_INTERVALL']][ITER], TEST = ntest, folds = FOLDS, THREADS = 5)
}


END = Sys.time()

print(END - START)              # Time difference : 3.745848 mins    [ for rgf loop ]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make submission
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mean RMSLE
#-----------

avg_rmsle_y1 = mean(unlist(lapply(res_RGF_y1, function(x) lapply(x, function(y) y$er_te))))

avg_rmsle_y2 = mean(unlist(lapply(res_RGF_y2, function(x) lapply(x, function(y) y$er_te))))

print(avg_rmsle_y1)
print(avg_rmsle_y2)


# average all predictions from all repeats
#-----------------------------------------

avg_PREDS_y1 = rowMeans(do.call(cbind, lapply(res_RGF_y1, function(x) do.call(cbind, lapply(x, function(y) y[[12]])))))

avg_PREDS_y2 = rowMeans(do.call(cbind, lapply(res_RGF_y2, function(x) do.call(cbind, lapply(x, function(y) y[[12]])))))



# submission
#-----------

sample_submission = read.csv(".../nomad2018/sample_submission.csv", header = T)                      # adjust to the folder where the sample-submission is saved
head(sample_submission)

sample_submission$formation_energy_ev_natom = avg_PREDS_y1
sample_submission$bandgap_energy_ev = avg_PREDS_y2

head(sample_submission)
summary(sample_submission)

# neg_idx = which(sample_submission$formation_energy_ev_natom < 0)                 # it might be the case that some predictions have negative values
# sample_submission[neg_idx, ]
# sample_submission$formation_energy_ev_natom[neg_idx] = 0.0                       # if negative set to 0.0

write.csv(sample_submission, ".../late_submission.csv", row.names=FALSE, quote = FALSE)

