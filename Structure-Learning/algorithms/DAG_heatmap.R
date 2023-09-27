# plot DAGs
DAG_heatmap <- function(W, text_bound = 0.05, low_col = "white", high_col = "red", true_graph = NULL) {
  
  p <- nrow(W)
  x <- 1:p
  data <- expand.grid(X=x, Y=x)
  data$Z <- as.vector(t(W))
  heatmap_W <- data %>% mutate(X = factor(X),
                               Y = factor(Y),
                               Z_value = round(ifelse(Z > text_bound, Z, NA), 2)) %>% 
    ggplot(aes(X, Y, fill= Z)) + 
    geom_tile() +
    scale_y_discrete(limits=rev, breaks = as.character(round(seq(1, p, length.out = 10), digits = 0))) +
    ylab("Parents nodes") +
    scale_x_discrete(position = "top", breaks = as.character(round(seq(1, p, length.out = 10), digits = 0))) +
    xlab("Children nodes") +
    scale_fill_gradient(low=low_col, high=high_col, name = "PEP") + 
    geom_text(aes(label = Z_value), na.rm = TRUE, size=2) 
  
  
  if (!is.null(true_graph)) {
    
    data_ture <- data %>% 
      mutate(W_true = as.vector(t(true_graph)))
    
    heatmap_W <- heatmap_W +
      geom_tile(data = data_ture, fill=NA,
                aes(factor(X), factor(Y),
                    fill = NA, colour = as.factor(W_true))) +
      scale_colour_manual(values=c("white", "black"), 
                          guide="none")
    
  }
  
  return(heatmap_W)
  
}











DAG_heatmap_cor <- function(W, text_bound = 0.05, low_col = "white", high_col = "red", true_graph = NULL) {
  
  p <- nrow(W)
  x <- 1:p
  data <- expand.grid(X=x, Y=x)
  data$Z <- as.vector(t(W))
  heatmap_W <- data %>% mutate(X = factor(X),
                               Y = factor(Y),
                               Z_value = round(ifelse(abs(Z) > text_bound, Z, NA), 2)) %>% 
    ggplot(aes(X, Y, fill= Z)) + 
    geom_tile() +
    scale_y_discrete(limits=rev, breaks = as.character(round(seq(1, p, length.out = 10), digits = 0))) +
    ylab("x") +
    scale_x_discrete(position = "top", breaks = as.character(round(seq(1, p, length.out = 10), digits = 0))) +
    xlab("X") +
    scale_fill_gradient(low=low_col, high=high_col, name = "correlation") + 
    geom_text(aes(label = Z_value), na.rm = TRUE) 
  
  
  if (!is.null(true_graph)) {
    
    data_ture <- data %>% 
      mutate(W_true = as.vector(t(true_graph)))
    
    heatmap_W <- heatmap_W +
      geom_tile(data = data_ture, fill=NA,
                aes(factor(X), factor(Y),
                    fill = NA, colour = as.factor(W_true))) +
      scale_colour_manual(values=c("white", "black"), 
                          guide="none")
    
  }
  
  return(heatmap_W)
  
}








DAG_heatmap_true <- function(W, low_col = "white", high_col = "red") {
  
  p <- nrow(W)
  x <- 1:p
  data <- expand.grid(X=x, Y=x)
  data$Z <- as.vector(t(W))
  heatmap_W <- data %>% mutate(X = factor(X),
                               Y = factor(Y),
                               Z = factor(Z, levels = c(0,1), labels = c("disconnected", "connected"))) %>% 
    ggplot(aes(X, Y, fill= Z)) + 
    geom_tile() +
    scale_y_discrete(limits=rev, breaks = as.character(round(seq(1, p, length.out = 10), digits = 0))) +
    ylab("Parents nodes") +
    scale_x_discrete(position = "top", breaks = as.character(round(seq(1, p, length.out = 10), digits = 0))) +
    xlab("Children nodes") +
    scale_fill_manual(values=c(low_col,high_col), name = "Direct Edges") +
    # theme(legend.position="none") +
    geom_tile(fill=NA, aes(X, Y,
                  fill = NA, colour = Z)) +
    scale_colour_manual(values=c("white", "black"), 
                        guide="none")
  
  
  return(heatmap_W)
  
}





