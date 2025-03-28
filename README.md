# Description  
This code implements elastic net regression (EN) on metabolomic data from two nested case-control studies using the R package glmnet. We applied 10-fold cross-validation, ensuring that samples from each matched case-control pair were assigned to the same fold, while samples from different studies were evenly distributed across the folds. The trained model was then validated on the other cohort, and the area under the curve (AUC) was computed using the R package pROC for (1) EN-selected metabolites, (2) known cancer risk factors (pack-years, diabetes, and BMI), and (3) all features combined. We used a bootstrap method implemented in pROC to compare AUC values across models.    
# Code
elasticnet_2cohorts.R: Elastic Net (EN) model for metabolomic data, using one cohort as the discovery dataset and the other as the validation dataset.  
elasticnet_2cohorts_exclude1year.R: EN model for metabolomic data, with one cohort as the discovery dataset and the other as the validation dataset, excluding samples with follow-up within the first year.    
elasticnet_2cohorts_5years.R: EN model for metabolomic data, using one cohort for discovery and the other for validation, restricting analysis to samples with follow-up of less than 5 years.  
compute_AUC_includingPRS.R: Computes AUC and p-values for different models incorporating known cancer risk factors, metabolomic markers, and polygenic risk scores (PRS).    

# Reference  
Ting Zhang, Steven C. Moore, Sheng Fu, et al. "Association between prediagnostic serum metabolic profiles and pancreatic ductal adenocarcinoma risk in two prospective cohorts." Submitted to International Journal of Cancer.  
