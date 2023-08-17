



CPU_time <- matrix(0, 3,3)
colnames(CPU_time) <- c("PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")
rownames(CPU_time) <- c("Logistic", "Cox's PH", "Weibull")

CPU_time[1,1] <- as.numeric(list.load("study6/logistic_aALA.rds")$CPU_time[1])
CPU_time[1,2] <- as.numeric(list.load("study6/logistic_LA.rds")$CPU_time[1])
CPU_time[1,3] <- as.numeric(list.load("study6/logistic_ALA.rds")$CPU_time[1])

CPU_time[2,1] <- as.numeric(list.load("study6/cox_aALA.rds")$CPU_time[1])
CPU_time[2,2] <- as.numeric(list.load("study6/cox_LA.rds")$CPU_time[1])
CPU_time[2,3] <- as.numeric(list.load("study6/cox_ALA.rds")$CPU_time[1])

CPU_time[3,1] <- as.numeric(list.load("study6/weibull_aALA.rds")$CPU_time[1])
CPU_time[3,2] <- as.numeric(list.load("study6/weibull_LA.rds")$CPU_time[1])
CPU_time[3,3] <- as.numeric(list.load("study6/weibull_ALA.rds")$CPU_time[1])



CPU_time <- data.frame(CPU_time) %>%
  mutate(Model = rownames(CPU_time))

CPU_time_cleaned <- CPU_time %>%
  pivot_longer(cols=starts_with("PARNI"), names_to="Algorithm") %>%
  mutate(Algorithm = gsub("\\.", "-", Algorithm)) %>%
  mutate(Algorithm = factor(Algorithm, levels = c("PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")),
         Model = factor(Model, levels = c("Logistic", "Cox's PH", "Weibull"))) 


p2 <- CPU_time %>%
  pivot_longer(cols=starts_with("PARNI"), names_to="Algorithm") %>%
  mutate(Algorithm = gsub("\\.", "-", Algorithm)) %>%
  mutate(Algorithm = factor(Algorithm, levels = c("PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")),
         Model = factor(Model, levels = c("Logistic", "Cox's PH", "Weibull")),
         value10 = value * 1.2) %>%
  ggplot(aes(x = Algorithm, y = value)) +
  facet_grid(rows = vars(Model), scales = "free_y") +
  geom_blank(aes(x = Algorithm, y = value10)) +
  geom_col() +
  geom_text(aes(label = round(value,2)), vjust = -0.5) + 
  scale_y_continuous(position = "right") +
  ylab("CPU time (in min)") +
  xlab("")



iterations <- 2000:10000

traces_logistic <- matrix(0, nrow = length(iterations), ncol = 5)
colnames(traces_logistic) <- c("Model", "Iterations", "PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")
traces_logistic <- data.frame(traces_logistic)

traces_logistic$Model <- "Logistic"
traces_logistic$Iterations <- iterations
traces_logistic$PARNI.adaptiveALA <- list.load("study6/logistic_aALA.rds")$log_post_trace[-(1:2000),1]
traces_logistic$PARNI.LA <- list.load("study6/logistic_LA.rds")$log_post_trace[-(1:2000),1]
traces_logistic$PARNI.ALA <- list.load("study6/logistic_ALA.rds")$log_post_trace[-(1:2000),1]


traces_cox <- matrix(0, nrow = length(iterations), ncol = 5)
colnames(traces_cox) <- c("Model", "Iterations", "PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")
traces_cox <- data.frame(traces_cox)

traces_cox$Model <- "Cox's PH"
traces_cox$Iterations <- iterations
traces_cox$PARNI.adaptiveALA <- list.load("study6/cox_aALA.rds")$log_post_trace[-(1:2000),1]
traces_cox$PARNI.LA <- list.load("study6/cox_LA.rds")$log_post_trace[-(1:2000),1]
traces_cox$PARNI.ALA <- list.load("study6/cox_ALA.rds")$log_post_trace[-(1:2000),1]



traces_weibull <- matrix(0, nrow = length(iterations), ncol = 5)
colnames(traces_weibull) <- c("Model", "Iterations", "PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")
traces_weibull <- data.frame(traces_weibull)

traces_weibull$Model <- "Weibull"
traces_weibull$Iterations <- iterations
traces_weibull$PARNI.adaptiveALA <- list.load("study6/weibull_aALA.rds")$log_post_trace[-(1:2000),1]
traces_weibull$PARNI.LA <- list.load("study6/weibull_LA.rds")$log_post_trace[-(1:2000),1]
traces_weibull$PARNI.ALA <- list.load("study6/weibull_ALA.rds")$log_post_trace[-(1:2000),1]



traces <- traces_logistic %>%
  add_case(traces_cox) %>%
  add_case(traces_weibull)


p1 <- traces %>% pivot_longer(cols=starts_with("PARNI"), names_to="Algorithm") %>%
  mutate(Algorithm = gsub("\\.", "-", Algorithm)) %>%
  mutate(Algorithm = factor(Algorithm, levels = c("PARNI-adaptiveALA", "PARNI-LA", "PARNI-ALA")),
         Model = factor(Model, levels = c("Logistic", "Cox's PH", "Weibull"))) %>%
  ggplot(aes(x = Iterations, y = value)) +
  geom_line() +
  facet_grid(rows = vars(Model), cols = vars(Algorithm), scales = "free_y", switch = "y") +
  ylab("Log posterior model probaiblity") +
  xlab("Iteration")
  # +
  #  scale_x_continuous(labels = scales::scientific)


library(gridExtra)

grid.arrange(p1, p2, 
             widths = c(3, 1),
             layout_matrix = matrix(c(1, 2), nrow = 1))


grid.arrange(p1+theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")),
             p2 +
               theme(plot.margin = margin(22, 5.5, 5.5, 5.5, "pt"),
                     strip.text = element_blank()),
             widths = c(2.7, 1),
             layout_matrix = matrix(c(1, 2), nrow = 1))







