

run_sim = function(seed = NA,
                   nsim = 500,
                   n = 600,
                   p = 20,
                   relative_weights = c(1,1,1)) {
  
  require(MASS);
  require(tidyverse);
  require(broom);
  require(glue);
  require(logistf);
  require(glmnet);
  require(glmnetUtils);
  
  if(!is.na(seed)) {
    set.seed(seed);
  }
  #Will split into equal parts
  proportions = cumsum(relative_weights / sum(relative_weights));
  training_subset = 1:round(proportions[1] * n);
  validation_subset = (round(proportions[1] * n) + 1):round(proportions[2]*n);
  testing_subset = (round(sum(proportions[2])*n)+1):n;
  intercept = log(0.6/0.4);
  beta <- c(rep(0.25, floor(p/2)), numeric(ceiling(p/2)));
  
  model_order = 
    c("(truth)","null","full","forward","forward_pruned", "backward", "firth", "lasso");
  results = 
    expand_grid(sim = 1:nsim, 
                model_name = model_order, 
                step = c("training", "validation", "testing"), 
                mspe = NA_real_, 
                absolute = NA_real_,
                zero_one = NA_real_,
                deviance = NA_real_) %>%
    mutate(row_number = 1:n()) %>%
    mutate(model_name = factor(model_name, levels = model_order), 
           step = factor(step, levels = c("training", "validation", "testing")))
  
  #simulating all data one time, outside of the for-loop, is faster than 
  #simulating at each step
  all_x <- matrix(rnorm(nsim*n*p), 
                  nrow = n * nsim);
  all_x <- all_x %*% chol(0.1 + diag(0.9, p));
  colnames(all_x) <- glue("x{1:p}");
  all_true_probs <- 1/(1 + exp(-intercept - all_x%*%beta));
  all_y <- rbinom(n*nsim, 1, all_true_probs);
  all_data <- bind_cols(y = all_y, data.frame(all_x));
  full_fmla <- 
    glue("y~{glue_collapse(glue('x{1:p}'),sep='+')}") %>%
    as.formula();
  
  for(i in 1:nsim) {
    curr_y = all_y[(i-1)*n + (1:n)]
    
    null_model = glm(y ~ 1,
                     data = all_data,
                     subset = (i-1)*n + training_subset, 
                     family = "binomial");
    
    forward_model <- 
      stepAIC(null_model,
              scope = list(upper = full_fmla),
              direction = "forward",
              trace = F);
    
    forward_pruned_model <- 
      stepAIC(null_model,
              scope = list(upper = full_fmla),
              direction = "forward",
              steps = max(1, floor(proportions[1] * n / 25)),
              trace = F);
    
    full_model = glm(full_fmla,
                     data = all_data,
                     subset = (i-1)*n + training_subset, 
                     family = "binomial");
    
    backward_model <- 
      stepAIC(full_model,
              direction = "backward",
              trace = F);
    
    firth_model <-
      logistf(full_fmla, 
              data = slice(all_data, training_subset))
    
    lasso_model <- 
      cv.glmnet(full_fmla, 
                data = all_data,
                subset = training_subset, 
                family = "binomial", alpha = 1)
    
    
    predict_models <- 
      list(null = null_model, 
           full = full_model, 
           forward = forward_model, 
           forward_pruned = forward_pruned_model, 
           backward = backward_model, 
           firth = firth_model, 
           lasso = lasso_model) %>%
      map(predict, newdata = all_data[(i-1)*n + (1:n),], type = 'resp') %>%
      map(as.numeric) %>%
      bind_cols() %>%
      mutate(`(truth)` = all_true_probs[(i-1)*n + (1:n)],
             y = curr_y, 
             training = row_number() %in% training_subset,
             validation = row_number() %in% validation_subset, 
             testing = row_number() %in% testing_subset, 
             sim = i) %>%
      pivot_longer(c(`(truth)`, null:lasso),
                   names_to = "model_name") %>%
      mutate(model_name = factor(model_name, levels = model_order)) %>%
      group_by(sim, model_name, training, validation, testing) %>% 
      summarize(mspe = mean((y - value)^2),
                absolute = mean(abs(y - value)), 
                zero_one = mean((1-y) * (value >= 0.5) + y * (value < 0.5)), 
                deviance = mean(-2 * (1-y) * log(1-value) - 2 * y * log(value)),
                .groups = "drop")
    
    #Training
    curr_training_index = 
      results %>% 
      filter(sim == i, step == "training") %>%
      pull(row_number)
    
    results[curr_training_index, c("mspe", "absolute", "zero_one", "deviance")] =
      predict_models %>% 
      filter(training) %>%
      arrange(model_name) %>% 
      dplyr::select(mspe, absolute, zero_one, deviance)
    
    
    #Validation
    curr_validation_index = 
      results %>% 
      filter(sim == i, step == "validation") %>%
      pull(row_number)
    
    results[curr_validation_index, c("mspe", "absolute", "zero_one", "deviance")] =
      predict_models %>% 
      filter(validation) %>%
      arrange(model_name) %>% 
      dplyr::select(mspe, absolute, zero_one, deviance)
    
    #Testing
    curr_testing_index = 
      results %>% 
      filter(sim == i, step == "testing") %>%
      pull(row_number)
    
    results[curr_testing_index, c("mspe", "absolute", "zero_one", "deviance")] =
      predict_models %>% 
      filter(testing) %>%
      arrange(model_name) %>% 
      dplyr::select(mspe, absolute, zero_one, deviance)
    
    #Selection / ranking based on validation step
    #results[curr_row_index[-c(1,2,7,8,13,14)],"ranking"] = 
    #  rep(rank(pull(results[curr_row_index[9:12],"value"]), ties.method = "random"), times = 3)
    
  }
  
  results <- 
    results %>% 
    group_by(sim) %>% 
    mutate(mspe_ranking = rank(mspe / (!model_name %in% c("(truth)", "null")) / (step == "validation"), ties.method = "average")) %>%
    group_by(sim, model_name) %>%
    mutate(mspe_ranking = 
             case_when(
               model_name %in% c("(truth)", "null") ~ Inf,
               TRUE ~ as.numeric(min(mspe_ranking)))) %>%
    ungroup()
  
  results;
  
}