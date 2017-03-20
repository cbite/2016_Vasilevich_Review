library(tidyverse)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(magrittr)
library(purrr)
library(MASS)
set.seed(1)

generate_clusters <- function(ndim, nclusters, cluster_sizes) {
  list(n = cluster_sizes, 
       mu = rerun(nclusters, rnorm(ndim) * 3), 
       Sigma = rerun(nclusters, cov(matrix(rnorm(ndim * 50), 50, ndim)))) %>%
    pmap(mvrnorm) %>%
    map(as.data.frame) %>%
    bind_rows()
}

x <- 
  generate_clusters(4, 3, list(5, 15, 10)) %>%
  sample_n(nrow(.)) %>%
  as.matrix()

x <- x - min(x)

f <- 5

y <- x %*% matrix(c(1, rnorm(1)/(f*1), rnorm(1)/(f*.2), rnorm(1)/(f*.2), 
                    1, rnorm(1)/(f*1), 1, rnorm(1)/(f*.01), 
                    rnorm(1)/(f*1), 1, rnorm(1)/(f*1), rnorm(1)/(f*1)), 4, 3)

err <- matrix(rnorm(dim(y)[1] * dim(y)[2]) * f, dim(y)[1], dim(y)[2])

#y <- y + err

df <- data_frame(x = x[,1], y = y[,1])

df$small_sample <- FALSE

df$small_sample[sample(30, 4)] <- TRUE

df<-df[!df$small_sample=="TRUE",]

m <- lm(y ~ x + 1, df)

p <- 
  ggplot(df, aes(x, y)) + 
  geom_point(aes(shape = small_sample)) + 
  geom_smooth(method=lm) + 
  coord_fixed(ratio = .45, 
              xlim = c(min(df$x), max(df$x)),
              ylim = c(min(df$y), max(df$y))) +
  xlab('grain size (a.u.)') + 
  ylab('cell proliferation (a.u.)') + 
  theme_classic() + 
  scale_color_manual(values = c("#000000", "#00FF00")) + 
  guides(color=FALSE, shape=FALSE)

xs <- 6

pred <- predict(m, data.frame(x = xs), se.fit = TRUE)

ys <- pred$fit

se <- 2 * pred$se.fit

# m_small <- lm(y ~ x + 1, df %>% filter(small_sample == TRUE))
# 
# y_small <- predict(m_small, data.frame(x = c(min(df$x), max(df$x))))

p <- p + 
  geom_segment(aes(x = xs, y = min(y), xend = xs, yend = ys), linetype = 2, color = "red") + 
  geom_segment(aes(x = min(x),  y = ys, xend = xs, yend = ys), linetype = 2, color = "red") +
  geom_segment(aes(x = min(x),  y = ys+se, xend = xs, yend = ys+se), linetype = 2, color = "darkgray") +
  geom_segment(aes(x = min(x),  y = ys-se, xend = xs, yend = ys-se), linetype = 2, color = "darkgray")
#  +
#  geom_segment(aes(x = min(x),  y = y_small[[1]], xend = max(x), yend = #y_small[[2]]), linetype = 1, color = "purple", size = 1)


p

t <- textGrob(expression(Y == a + b * X))

ggsave("grain_size_v_cell_proliferation.png", 
       grid.arrange(p, t, ncol=2), 
       width = 6, height = 4)

ggsave("grain_size_v_cell_proliferation.pdf", 
       grid.arrange(p, t, ncol=2), 
       width = 6, height = 4)

format_matrix <- function(mat) {
  df <- as.data.frame(mat)
  names(df) <- seq(ncol(mat))
  df  %<>% 
    rownames_to_column() %>% 
    rename(row = rowname) %>%
    mutate(row = as.integer(row))
  df %<>% gather(col, val, -row)
  df
}

plot_matrix <- function(df, palette = "BuGn", 
                        xlab = "features",
                        ylab = "samples") {
  p <- 
    ggplot(df, aes(row, col, fill=val)) + 
    geom_tile(color="white", size = 1) + 
    coord_equal() + 
    theme(line = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_text(angle = 0)) + 
    guides(fill=FALSE) +
    xlab(xlab) + ylab(ylab)
  
  if (length(palette) == 1) {
    p <- p + scale_fill_brewer(palette = palette)
  } else {
    p <- p + scale_fill_manual(values = palette)
  }
  
  p + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
}

discrete <- function(s, breaks = 8) {
  s %>% 
    cut(breaks, 
        labels = seq(breaks), 
        ordered_result = TRUE) %>%
    matrix(dim(s)[1], 
           dim(s)[2])
}
p1 <- 
  plot_matrix(format_matrix(discrete(as.matrix(y[,1]))), 
              palette = "BuGn",
              ylab = "cell\nproliferation\n(1 feature)",
              xlab = "samples (n=30)")

p2 <- 
  plot_matrix(format_matrix(as.matrix(discrete(x))), 
              palette = "BuPu",
              ylab = "biomaterial\ndesign\n(4 features)",
              xlab = "samples (n=30)")

t <- textGrob(expression(Y == a + b[1] * X[1] + b[2] * X[2] + b[3] * X[3] + b[4] * X[4]))

p <- grid.arrange(p2, p1, t, nrow=3, heights = c(.01, .01, .01))

p 

ggsave("all_params_v_cell_proliferation.png", p, width = 8, height = 5)

ggsave("all_params_v_cell_proliferation.pdf", p, width = 8, height = 5)

p1 <- 
  plot_matrix(format_matrix(discrete(as.matrix(y[,1]), 3)), 
              palette =  c("#FDC086", "#AAAAAA", "#7FC97F"),
              ylab = "cell\nproliferation\n(1 feature)",
              xlab = "samples (n=30)")

t <- textGrob(expression(Y == f(X[1], X[2], X[3], X[4])))

p <- grid.arrange(p2, p1, t, nrow=3, heights = c(.01, .01, .01))

p 

ggsave("all_params_v_cell_proliferation_classification.png", p, width = 8, height = 5)

ggsave("all_params_v_cell_proliferation_classification.pdf", p, width = 8, height = 5)

x_df <- as.data.frame(x)
x_df$class <- discrete(as.matrix(y[,1]), 3)
x_df$class <- factor(x_df$class, labels = c("high", "medium", "low"))

x_new <- data.frame(V1 = 6, V2 = 3.5 )

library(FNN)
k <- knn(x_df %>% dplyr::select(V1, V2), 
         x_new,
         x_df$class, 
         k = 5, algorithm="cover_tree")
indices <- attr(k, "nn.index")

x_df_nn <- x_df[as.vector(indices),]

max_neighbor_dist <- tail(as.vector(attr(k, "nn.dist")), 1)

p <- 
  ggplot(data=x_df, aes(x=V1, y=V2)) +
  geom_point(
    aes(fill=class),
    pch=21,
    size = 4,
    color="#333333") + 
  geom_point(data = x_new, color = "red", size = 4) +
  coord_equal() + 
  scale_fill_manual(values = c("#FDC086", "#AAAAAA", "#7FC97F")) + 
  xlab("Proliferation") + ylab("Protein level")

p 
ggsave("knn.png", p, width = 4, height = 4)

ggsave("knn.pdf", p, width = 4, height = 4)

library(plotly)
library(RColorBrewer)

m <- kmeans(y, 3)

y_df <- as.data.frame(y)

y_df %<>% mutate(cluster = as.factor(m$cluster))

p <- plot_ly(y_df, 
             x = ~V1, 
             y = ~V2, 
             z = ~V3, 
             color = ~cluster, 
             colors = brewer.pal(3, "Set1")) %>%
  add_markers() %>%
  layout(scene = list(zaxis = list(title = 'ECM'),
                      yaxis = list(title = 'Protein\n level'),
                      xaxis = list(title = 'Proliferation')))
p
library(ggdendro)

model <- cluster::agnes(y, metric = "manhattan", stand = TRUE)
dg <- as.dendrogram(model)
p <- 
  ggdendrogram(dg) + 
  geom_hline(yintercept = 5, color = "red", linetype = 2)+
  geom_hline(yintercept = 3.5, color = "blue", linetype = 2)
p
ggsave("hier_clustering.png", p, width = 4, height = 4)
ggsave("hier_clustering.pdf", p, width = 4, height = 4)

p1 <- 
  plot_matrix(format_matrix(discrete(as.matrix(scale(y)))), 
              palette = "BuGn",
              ylab = "cell\nproliferation\n(3 features)",
              xlab = "samples (n=30)")

ggsave("all_params_v_cell_proliferation_discrete.pdf", p1, width = 8, height = 5)
ggsave("all_params_v_cell_proliferation_discrete.png", p1, width = 8, height = 5)


p2 <- qplot() + 
  annotation_custom(rasterGrob(png::readPNG("../results/manual/3dcrop.png"))) + theme_void()

p3 <- qplot() + 
  annotation_custom(rasterGrob(png::readPNG("hier_clustering.png"))) + theme_void()

p <- grid.arrange(p1, 
                  grid.arrange(p2, p3, 
                               ncol=2, 
                               widths = c(1, 1)), 
                  nrow = 2, 
                  heights = c(.5, 1))

ggsave("clustering.png", p, width = 8, height = 5)

p