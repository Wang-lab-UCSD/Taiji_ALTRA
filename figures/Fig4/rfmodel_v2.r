# v2 also changes split of train_df and test_df
library(caret)
library(dplyr)
library(parallel)
setwd('/home/jupyter/output_20230310/post-analysis/')

# data2 <- read.csv('data_for_rf_model.csv', row.names = 1)
# data2 <- read.csv('data_for_rf_model_70_cytokines.csv', row.names = 1)
# data2 <- read.csv('data_for_rf_model_no_receptor.csv', row.names = 1)
# data2 <- read.csv('data_for_rf_model_270_genes.csv', row.names = 1)
data2 <- read.csv('data_for_rf_model_converter.csv', row.names = 1)
# data2 <- read.csv('data_for_rf_model_converter_both_C4.csv', row.names = 1)
# data2 <- read.csv('data_antibody_for_rf_model.csv', row.names = 1)
# data2 <- read.csv('data_antibody_for_rf_model_converter.csv', row.names = 1)

# max_feat = 67
# max_feat = 34
max_feat = 26
# max_feat = 4
numCores = detectCores()
print(numCores)

# from Barton's code
randomForestModelFitter = function(outer_seed, max_featureCount=max_feat, dependent_var='class'){
    
    # set random seed
    set.seed(outer_seed)
    
    # split data in train and test data
    smp_size <- floor(0.7 * nrow(data2))
    train_ind <- sample(seq_len(nrow(data2)), size = smp_size)
    train_df <- data2[train_ind, ]
    train_df$class <- factor(train_df$class)
    test_df <- data2[-train_ind, ]
    test_df$class <- factor(test_df$class)

    modelData_list = list()
    runName = paste0("seed", outer_seed)

    predictors = paste(names(train_df)[-c(1,2)], collapse = " + ")
    model_formula = as.formula(paste0(dependent_var, " ~ ", predictors))

    ### Rank most important predictors
    # Set RFE control
    ctrl = rfeControl(functions = rfFuncs,
                      method = "repeatedcv", repeats = 10,
                      saveDetails = TRUE)
    # Set a sequence of feature-space sizes to search over:
    sizes = seq(1, max_featureCount, by = 1)
    # Use caret's rfe function to fit RF models to these different feature spaces
    rfeResults = rfe(model_formula, data = train_df, 
                     sizes = sizes, 
                     rfeControl = ctrl, 
                     verbose = TRUE)
    # get features ranked by importance
    varimp_df = data.frame(feature = row.names(varImp(rfeResults)), importance = varImp(rfeResults)[, 1])
    varimp_df_name = paste0(runName, "_varimpDF")
    modelData_list[[varimp_df_name]] = varimp_df  ### STORE ALL PREDICTORS BY IMPORTANCE
    # remove negative importance predictors
    varimp_df = varimp_df[which(varimp_df$importance > 0), ]
    max_featureCount = min(max_featureCount, dim(varimp_df)[1])
    ### Train Random Forest Model
    for (featureCount in seq(1, max_featureCount, 1)){        
        featureCount = min(featureCount, dim(varimp_df)[1])
        # select features
        selectedFeatures_vector = varimp_df[1:featureCount, "feature"]
        # define model formula
        predictors = paste(selectedFeatures_vector, collapse = " + ")
        model_formula = as.formula(paste0(dependent_var, " ~ ", predictors))
        # Define parameters
        control = trainControl(method='repeatedcv', number=10, repeats=5)
        metric = "Accuracy"
        mtry_upper = ceiling(sqrt(length(selectedFeatures_vector)))
        grid = expand.grid(mtry=seq(1, mtry_upper, 1))
        # train model
        set.seed(outer_seed)
        rfClassModel = train(model_formula, data=train_df, method='rf', trControl=control, metric=metric, tuneGrid=grid)
        rfClassModel_name = paste0(runName, "_predictors", featureCount, "_rfClass_model")
        modelData_list[[rfClassModel_name]] = rfClassModel  ### STORE MODEL
        ### Test performance
        # confusion matrix
        predictions = predict(object=rfClassModel, newdata=test_df)
        rfClassModel_testPredictions_name = paste0(runName, "_predictors", featureCount, "_rfClass_testPredictions")
        modelData_list[[rfClassModel_testPredictions_name]] = predictions  ### STORE PREDICTIONS DF
        confusion_mtrx = confusionMatrix(data=predictions, test_df[[dependent_var]])
        rfClassModel_testConfusionMtrx_name = paste0(runName, "_predictors", featureCount, "_rfClass_testConfusionMtrx")
        modelData_list[[rfClassModel_testConfusionMtrx_name]] = confusion_mtrx  ### STORE CONFUSION MATRIX
        # AUC
        roc_score = pROC::roc(test_df[[dependent_var]], predict(rfClassModel, newdata=test_df, type = "prob")[['1']]) #AUC score
        rfClassModel_testAUC_name = paste0(runName, "_predictors", featureCount, "_rfClass_testAUC")
        modelData_list[[rfClassModel_testAUC_name]] = roc_score  ### STORE AUC
        }
    saveRDS(modelData_list, file = paste0("rfmodels/rfmodel_converter_C2/seed",outer_seed,"_rfClass_list.RDS"))
    }

tmp = mclapply(1:100,randomForestModelFitter,mc.cores=numCores)
# tmp = lapply(1:100,randomForestModelFitter)
