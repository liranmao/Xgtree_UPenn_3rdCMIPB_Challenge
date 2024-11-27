# Xgtree_UPenn_3rdCMIPB_Challenge

Machine learning model implementation for predicting vaccine responses as part of the CMI-PB Challenge. The model predicts three primary metrics and their fold changes:
- Monocyte levels (Day 1)
- CCL3 expression (Day 3) 
- IgG-PT antibody titers (Day 14)

## Data Processing

Input features:
- Demographics: age at boost, biological sex, infancy vaccination status
- Baseline measurements (Day 0): Monocytes, CCL3, IgG-PT

Preprocessing steps:
- Gene expression data: log2(TPM + 1) transformation
- Categorical variables: One-hot encoding
- Missing values: Median imputation
- Response variables: Log2 fold change from baseline

## Model Architecture

We evaluated many models and found that XGBoost tree performs best so in the end we submitted XGBoost tree results. Here is a detailed description about the whole evaluation step:

Model evaluation framework includes:

Linear Models:
- Linear Regression
- Lasso (α=1)
- Ridge (α=0) 
- Elastic Net (α=0.1-0.9)

Non-linear Models:
- Random Forest
- SVM (Linear/Radial kernels)
- Gradient Boosting
- XGBoost (Linear/Tree)

Training approach:
- Leave-One-Out Cross-Validation (LOOCV)
- Optimization metric: RMSE
- Performance evaluation: Correlation between predictions and observations


## Files

- `lm_age_baseline_model.R`: Model training implementation
- `lm_age_baseline_model_submission.R`: Prediction generation


## Contact
Liran Mao, Liran.Mao@Pennmedicine.upenn.edu, University of Pennsylvania

Lexi Li, lexi.li@pennmedicine.upenn.edu, University of Pennsylvania

Yuhao Tan: Yuhao.Tan@pennmedicine.upenn.edu, University of Pennsylvania

Tram Anh Nguyen: tram.nguyen1@pennmedicine.upenn.edu, University of Pennsylvania

