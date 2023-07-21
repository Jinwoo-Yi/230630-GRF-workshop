### ----------------------------------------------------
### GRF Workshop: Hands-on Tutorials 
### date: 2023/06/30
### ----------------------------------------------------
### written by
### Jinwoo, Yi 
### adem1997@snu.ac.kr
### SNU Connectome Lab
### ----------------------------------------------------

library(caret)
library(lmtest)
library(sandwich)
library(grf)
library(glmnet)
library(splines)
library(ggplot2)
library(dplyr)
library(rpart)
library(MASS)
library(tibble)
library(lsr)   # ADD!
library(stats) # ADD!

### @@@ PART I. Data Loading and Preprocessing @@@ ####
## > I-1. to import dataset ####
data <- read.csv("/Users/Jinwoo_1/Desktop/grf_tutorial_data.csv")
summary(data)

## > I-2. to define variable types ####
factor.idx <- c(1, 2, 4) 
for (i in factor.idx) {data[[i]] <- as.factor(data[[i]])}

## > I-3. to check near-zero-variance ####
data.nzv <- nearZeroVar(data, saveMetrics = TRUE) # there is not any nzv variable!

## > I-4. standardization for continuous variables ####
scale.model <- preProcess(data, method = c("center", "scale")) 
data.std <- predict(scale.model, data)
summary(data.std)

## > I-5. to make all variables numeric ####
for (i in c(1:ncol(data.std))) {data.std[[i]] <- as.numeric(data.std[[i]])}
data.std$dep_diag <- data.std$dep_diag - 1
data.std$drug <- data.std$drug - 1
data.std$sex_M <- data.std$sex_M - 1

### @@@ PART II. Fitting Generalized Random Forests @@@ ####
## > II-1. to define X, Y, and W ####
X <- data.std[, c(3:8)]
Y <- data.std$dep_diag
W <- data.std$drug

## > II-2. to fit causal forest and assess model fit ####
set.seed(20230630, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
cf <- causal_forest(X, Y, W, Y.hat = NULL, W.hat = NULL, tune.parameters = "all", num.trees = 2000, tune.num.trees = 100)

cf.fit <- test_calibration(cf)
cf.fit

### @@@ PART III. ATE and CATE estimation @@@ ####
## > III-1. ATE estimation ####
ate.est <- average_treatment_effect(cf, method = "AIPW", num.trees.for.weights = 1000)[1] 
ate.se <- average_treatment_effect(cf, method = "AIPW", num.trees.for.weights = 1000)[2]
ate.tstat <- ate.est / ate.se
ate.pvalue <- 1.96 * (1 - pnorm(abs(ate.tstat)))
ate.summary <- c(estimate = ate.est, std.error = ate.se, t.stat = ate.tstat, p.value = ate.pvalue)
ate.summary

## > III-2. CATE estimation ####
cate <- get_scores(cf)
data[, "cate"] <- cate
hist(data$cate)


### @@@ PART IV. Assessment of Heterogeneity: Median Comparison @@@ ####
data$group <- ifelse(data$cate > median(cate), "high", "low")
data$group <- as.factor(data$group)
summary(data)

t.test(data[data$group == 'high', 'cate'], data[data$group == 'low', 'cate'])
cohensD(data$cate[data$group == "high"], data$cate[data$group == "low"])


### @@@ PART V. Risk/protective factor analysis ####
## > V-1. variable importance ####
varimp <- c(variable_importance(cf)) 
names(varimp) <- colnames(X)
varimp <- sort(varimp, decreasing = TRUE)
varimp <- data.frame(varimp)
varimp[, "rank"] <- c(1:length(colnames(X)))
varimp <- rownames_to_column(varimp, var = "feature")

## > V-2. group comparison between 'high' and 'low' groups ####
# to define variable types again for precise test
for (i in factor.idx) {data[[i]] <- as.factor(data[[i]])}
summary(data)

sum(varimp$varimp)

# to create empty archive dataframe
compare.x <- data.frame(colnames(X))
names(compare.x) <- "feature"
compare.x[, "high_mean"] <- 0
compare.x[, "high_sd"] <- 0
compare.x[, "low_mean"] <- 0
compare.x[, "low_sd"] <- 0
compare.x[, "test_t"] <- 0
compare.x[, "test_p"] <- 0

# to perform Welch's T-test for group comparsion only for continuous variables first
for (i in 3:8) {
  compare.x[i-2, "high_mean"] <- ifelse(is.factor(data[[i]]) == TRUE, NA, 
                                      mean(data[data$group == 'high', colnames(data)[i]]))
  compare.x[i-2, "high_sd"] <- ifelse(is.factor(data[[i]]) == TRUE, NA, 
                                    sd(data[data$group == 'high', colnames(data)[i]]))
  compare.x[i-2, "low_mean"] <- ifelse(is.factor(data[[i]]) == TRUE, NA, 
                                     mean(data[data$group == 'low', colnames(data)[i]]))
  compare.x[i-2, "low_sd"] <- ifelse(is.factor(data[[i]]) == TRUE, NA, 
                                   sd(data[data$group == 'low', colnames(data)[i]]))
  compare.x[i-2, "test_t"] <- ifelse(is.factor(data[[i]]) == TRUE, NA, 
                                   t.test(data[data$group == 'high', i], data[data$group == 'low', i])$statistic[[1]])
  compare.x[i-2, "test_p"] <- ifelse(is.factor(data[[i]]) == TRUE, NA, 
                                   t.test(data[data$group == 'high', i], data[data$group == 'low', i])$p.value)
}

# to perform Fisher's Exact Test for discrete variable (i.e., sex_M)
nrow(data[data$group == 'high' & data$sex_M == '1', ])
nrow(data[data$group == 'high' & data$sex_M == '0', ])
nrow(data[data$group == 'low' & data$sex_M == '1', ])
nrow(data[data$group == 'low' & data$sex_M == '0', ])

sex_table <- data.frame("high" = c(nrow(data[data$group == 'high' & data$sex_M == '1', ]), nrow(data[data$group == 'high' & data$sex_M == '0', ])),
                        "low"  = c(nrow(data[data$group == 'low' & data$sex_M == '1', ]), nrow(data[data$group == 'low' & data$sex_M == '0', ])),
                        row.names = c("1", "2"))
sex_table

compare.x[2, 'test_p'] <- fisher.test(sex_table)$p.value
compare.x[2, 'high_mean'] <- nrow(data[data$group == 'high' & data$sex_M == '1', ]) / nrow(data)
compare.x[2, 'low_mean'] <- nrow(data[data$group == 'low' & data$sex_M == '1', ]) / nrow(data)

# FDR correction for p-values
compare.x[, "fdr_p"] <- p.adjust(compare.x$test_p, method = 'fdr')

## > V-3. best linear projection ####
blp <- best_linear_projection(cf, X)
blp <- as.data.frame(blp[, ])
colnames(blp) <- c("estimate", "std.error", "t.stat", "p.value")
blp <- rownames_to_column(blp, var = "feature")

## > V-4. summary of feature analysis ####
feature.summary <- cbind(compare.x, blp[-1, -1])

# tiering features
feature.summary <- feature.summary %>%
  mutate(tier = case_when(fdr_p < 0.05 & p.value < 0.05 ~ '1',
                          fdr_p < 0.05 | p.value < 0.05 ~ '2',
                          TRUE ~ '3'))

# typing features
feature.summary <- feature.summary %>%
  mutate(type = case_when(high_mean > low_mean & estimate > 0 ~ 'risk',
                          high_mean < low_mean & estimate < 0 ~ 'protective',
                          TRUE ~ 'mixed'))
