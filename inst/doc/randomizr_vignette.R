## ------------------------------------------------------------------------
# Load built-in dataset
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)

# Transform so each row is a subject
# Columns describe subject's hair color, eye color, and gender
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                        times = HairEyeColor$Freq), 1:3]

N <- nrow(hec)

# Fix the rownames
rownames(hec) <- 1:N

## ------------------------------------------------------------------------
# Set a seed for reproducability
set.seed(343)

# Create untreated and treated outcomes for all subjects
hec <- within(hec,{
  Y0 <- rnorm(n = N,mean = (2*as.numeric(Hair) + -4*as.numeric(Eye) + -6*as.numeric(Sex)), sd = 5)
  Y1 <- Y0 + 6*as.numeric(Hair) + 4*as.numeric(Eye) + 2*as.numeric(Sex)
})

# Calculate true ATE
with(hec, mean(Y1 - Y0))

## ------------------------------------------------------------------------
library(randomizr)
Z <- simple_ra(N = N)
table(Z)

## ------------------------------------------------------------------------
Z <- simple_ra(N = N, prob = 0.30)
table(Z)

## ------------------------------------------------------------------------
Z <- simple_ra(N = N, num_arms = 3)
table(Z)

## ------------------------------------------------------------------------
Z <- simple_ra(N = N, prob_each = c(.2, .2, .6))
table(Z)

## ------------------------------------------------------------------------
Z <- simple_ra(N = N, prob_each = c(.2, .2, .6),
               condition_names=c("control", "placebo", "treatment"))
table(Z)

## ------------------------------------------------------------------------
Z <- complete_ra(N = N)
table(Z)

## ------------------------------------------------------------------------
Z <- complete_ra(N = N, m=200)
table(Z)

## ------------------------------------------------------------------------
Z <- complete_ra(N = N, num_arms = 3)
table(Z)

## ------------------------------------------------------------------------
Z <- complete_ra(N = N, m_each = c(100, 200, 292))
table(Z)

## ------------------------------------------------------------------------
Z <- complete_ra(N = N, m_each = c(100, 200, 292),
               condition_names=c("control", "placebo", "treatment"))
table(Z)

## ------------------------------------------------------------------------
sims <- 1000

# Set up empty vectors to collect results
simple_ests <- rep(NA, sims)
complete_ests <- rep(NA, sims)

# Loop through simulation 2000 times
for(i in 1:sims){
  hec <- within(hec,{
    
    # Conduct both kinds of random assignment
    Z_simple <- simple_ra(N = N)
    Z_complete <- complete_ra(N = N)
    
    # Reveal observed potential outcomes
    Y_simple <- Y1*Z_simple + Y0*(1-Z_simple)
    Y_complete <- Y1*Z_complete + Y0*(1-Z_complete)
    })
  
  # Estimate ATE under both models
  fit_simple <- lm(Y_simple ~ Z_simple, data=hec)
  fit_complete <- lm(Y_complete ~ Z_complete, data=hec)
  
  # Save the estimates
  simple_ests[i] <- coef(fit_simple)[2]
  complete_ests[i] <- coef(fit_complete)[2]
}

## ------------------------------------------------------------------------
sd(simple_ests)
sd(complete_ests)

## ------------------------------------------------------------------------
Z <- block_ra(block_var = hec$Hair)
table(Z, hec$Hair)

## ------------------------------------------------------------------------
Z <- block_ra(block_var = hec$Hair, num_arms=3)
table(Z, hec$Hair)
Z <- block_ra(block_var = hec$Hair, condition_names=c("Control", "Placebo", "Treatment"))
table(Z, hec$Hair)

## ------------------------------------------------------------------------
Z <- block_ra(block_var = hec$Hair,block_prob = c(.3, .7))
table(Z, hec$Hair)

## ------------------------------------------------------------------------
sort(unique(hec$Hair))
block_m <- cbind(c(78, 186, 51, 87),c(30, 100, 20, 40))
block_m
Z <- block_ra(block_var = hec$Hair,block_m=block_m, condition_names = c(0, 1))
table(Z, hec$Hair)

## ------------------------------------------------------------------------
hec <- within(hec,{
  Z_blocked <- block_ra(block_var = hec$Hair,block_m=block_m, condition_names = c(0, 1))
  Y_blocked <- Y1*(Z_blocked) + Y0*(1-Z_blocked)
  prob_Z <- rep(NA)
  prob_Z[Hair %in% c("Black")] <- 30/108
  prob_Z[Hair %in% c("Brown")] <- 100/286
  prob_Z[Hair %in% c("Red")] <- 20/71
  prob_Z[Hair %in% c("Blond")] <- 40/127
  IPW_weights <- 1/(prob_Z*Z_blocked + (1 - prob_Z)*(1-Z_blocked))
})

fit_control_for_blocks <- lm(Y_blocked ~ Z_blocked + Hair, data=hec)
fit_IPW <- lm(Y_blocked ~ Z_blocked, weights=IPW_weights, data=hec)

## ------------------------------------------------------------------------
suppressMessages(library(dplyr))
block_id <- id(hec[,c("Hair", "Eye", "Sex")])
block_var <- paste0("block_", sprintf("%02d", block_id))
table(block_var)

## ------------------------------------------------------------------------
Z_blocked <- block_ra(block_var = block_var)
head(table(block_var, Z_blocked))

## ------------------------------------------------------------------------
library(blockTools)

# BlockTools requires that all variables be numeric
numeric_mat <- model.matrix(~Hair+Eye+Sex, data=hec)[,-1]

# BlockTools also requres an id variable
df_forBT <- data.frame(id_var = 1:nrow(numeric_mat), numeric_mat)

# Conducting the actual blocking: let's make trios
out <- block(df_forBT, n.tr = 3, id.vars = "id_var", 
             block.vars = colnames(df_forBT)[-1])

# Extact the block_ids
hec$block_id <- createBlockIDs(out, df_forBT, id.var = "id_var")

# Conduct actual random assignment with randomizr
Z_blocked <- block_ra(block_var = hec$block_id, num_arms = 3)
head(table(hec$block_id, Z_blocked))

## ------------------------------------------------------------------------
clust_id <- id(hec[,c("Hair", "Eye", "Sex")])
clust_var <- paste0("clust_", sprintf("%02d", clust_id))

Z_clust <- cluster_ra(clust_var = clust_var)
head(table(clust_var, Z_clust))

## ------------------------------------------------------------------------
Z_clust <- cluster_ra(clust_var=clust_var, num_arms=3)
head(table(clust_var, Z_clust))

## ------------------------------------------------------------------------
Z_clust <- cluster_ra(clust_var=clust_var, 
                      condition_names=c("Control", "Placebo", "Treatment"))
head(table(clust_var, Z_clust))

## ------------------------------------------------------------------------
Z_clust <- cluster_ra(clust_var=clust_var, m_each=c(5, 15, 12))
head(table(clust_var, Z_clust))

## ------------------------------------------------------------------------
# Create cluster variable
clust_id <- id(hec[,c("Hair", "Eye", "Sex")])
clust_var <- paste0("clust_", sprintf("%02d", clust_id))

# Sort cluster names by increasing size
cluster_names_by_size <- names(sort(table(clust_var)))

# Create block variable
block_var <- rep(NA, nrow(hec))
for(i in 1:16){
  block_var[clust_var %in% cluster_names_by_size[((2*i)-1):(2*i)]] <- 
    paste0("block_", sprintf("%02d", i))
}

special_ra <- function(){
  assign <- rep(NA, length(clust_var))
  unique_blocks <- unique(block_var)
  for(i in 1:length(unique_blocks)){
    assign[block_var==unique_blocks[i]] <-
      cluster_ra(clust_var[block_var==unique_blocks[i]])
  }
  return(assign)
}
 
Z_special <- special_ra()
head(table(clust_var, Z_special))
head(table(block_var, Z_special))

## ------------------------------------------------------------------------
# 400 families have 1 child in the lottery, 100 families have 2
family_id <- c(sprintf("%03d", 1:500), sprintf("%03d", 1:100))

school_ra <- function(m){
  N <- length(family_id)
  random_number <- sample(1:N, replace=FALSE)
  Z <- rep(0, N)
  i <- 1
  while(sum(Z) <m){
    Z[family_id==family_id[random_number[i]]] <- 1
    i <- i + 1
  }
  return(Z)
}

Z <- school_ra(200)
table(Z)

## ------------------------------------------------------------------------
Z_matrix <- replicate(1000, school_ra(200))
plot(rowMeans(Z_matrix))

## ----,eval=FALSE---------------------------------------------------------
#  hec <- within(hec,{
#    Z_blocked <- complete_ra(N = N, m_each = c(100, 200, 292),
#                 condition_names=c("control", "placebo", "treatment"))
#    id_var <- 1:nrow(hec)
#  })
#  write.csv(hec[,c("id_var", "Z_blocked")], file="MyRandomAssignment.csv")

