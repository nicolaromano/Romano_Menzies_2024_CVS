library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(emmeans)
library(ggpmisc)
library(tidyr)
library(tibble)
library(gridExtra)
library(purrr)
library(corrplot)
library(pheatmap)

my_theme <- theme_minimal() +
    theme(
        text = element_text(size = 16),
        legend.position = "none"
    )

effect_size <- read.csv("effect_size.csv")

get_cohen_d <- function(n1, n2, mu1, mu2, s1, s2) {
    #' Calculate Cohen's d
    #' Parameters:
    #' n1: sample size of group 1
    #' n2: sample size of group 2
    #' mu1: mean of group 1
    #' mu2: mean of group 2
    #' s1: standard deviation of group 1
    #' s2: standard deviation of group 2
    #' Returns:
    #' d: Cohen's d

    # Convert all parameters to numeric
    n1 <- as.numeric(n1)
    n2 <- as.numeric(n2)
    mu1 <- as.numeric(mu1)
    mu2 <- as.numeric(mu2)
    s1 <- as.numeric(s1)
    s2 <- as.numeric(s2)

    s_pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    d <- (mu1 - mu2) / s_pooled

    return(d)
}

tests <- c("FST", "SPT", "TST", "EPM", "OFT")

for (test in tests) {
    # Calculate Cohen's d for the current test
    effect_size[[paste0(test, "_Cohen_d")]] <- get_cohen_d(
        effect_size[[paste0(test, "_Test_n_used")]],
        effect_size[[paste0(test, "_Control_n_used")]],
        effect_size[[paste0(test, "_Test_mean")]],
        effect_size[[paste0(test, "_Control_mean")]],
        effect_size[[paste0(test, "_Test_back_calculated_SD")]],
        effect_size[[paste0(test, "_Control_back_calculated_SD")]]
    )
}

effect_size %>%
    select(
        PMID, First_Author, Year, Species, Duration, Burden, Diversity, contains("Cohen_d")
    ) %>%
    mutate(Reference = paste(First_Author, Year, sep = ", "), .after = PMID) %>%
    select(-First_Author, -Year) %>%
    # All columns containing "effect size"
    pivot_longer(cols = c(contains("Cohen_d")), names_to = "Test", values_to = "Effect_size") %>%
    mutate(Test = gsub(paste0("_", "Cohen_d"), "", Test)) %>%
    mutate(Species = factor(Species)) -> effect_size_long

write.csv(effect_size_long, "effect_size_long.csv", row.names = FALSE)

correlation_data <- effect_size_long %>%
    group_by(Test) %>%
    summarize(
        correlation_dur = cor(Duration, Effect_size, use = "complete.obs", method = "pearson"),
        correlation_test_dur = cor.test(Duration, Effect_size, method = "pearson")$p.value,
        correlation_bur = cor(Burden, Effect_size, use = "complete.obs", method = "pearson"),
        correlation_test_bur = cor.test(Burden, Effect_size, method = "pearson")$p.value,
        correlation_div = cor(Diversity, Effect_size, use = "complete.obs", method = "pearson"),
        correlation_test_div = cor.test(Diversity, Effect_size, method = "pearson")$p.value
    )

# Figure 2 - Bubble chart of size, duration and burden
jpeg("bubble_chart.jpg", width = 1024, height = 768, res = 300, quality = 100)
effect_size %>%
    # Remove duplicated entries based on First_Author, Journal, Year
    distinct(First_Author, Year, .keep_all = TRUE) %>%
    ggplot() +
    geom_point(aes(
        x = Duration, y = Burden, size = Diversity,
        col = Species
    ), alpha = 0.3) +
    scale_color_manual(
        values = c("purple", "green"), name = "Species",
        labels = c("Mouse", "Rat")
    ) +
    scale_size_continuous(range = c(1, 30)) +
    xlab("Duration (days)") +
    xlim(0, 60) +
    ylim(0, 160) +
    my_theme +
    theme(legend.position = "right") +
    guides(col = guide_legend(override.aes = list(size = 15)))
dev.off()

# Figure 3 - Correlation between effect size and duration, burden, and diversity for each test, both species combined; gray band indicates 95% confidence intervals.


plot_correlation <- function(split_by_species = FALSE, show_R2 = FALSE) {
    if (split_by_species) {
        regression_aes <- aes(group = Species, lty = Species)
        points_aes <- aes(shape = Species)
    } else {
        regression_aes <- aes()
        points_aes <- aes()
    }

    p1 <- ggplot(effect_size_long, aes(x = Duration, y = Effect_size, color = Test)) +
        geom_point(points_aes, size = 2) +
        scale_shape_manual(values = c(21, 20)) +
        stat_poly_line(regression_aes, se = !split_by_species) +
        xlim(0, 60) +
        ylab("Effect size") +
        facet_wrap(~Test, ncol = 1) +
        theme_minimal() +
        my_theme

    p2 <- ggplot(effect_size_long, aes(x = Burden, y = Effect_size, color = Test)) +
        geom_point(points_aes, size = 2) +
        scale_shape_manual(values = c(21, 20)) +
        stat_poly_line(regression_aes, se = !split_by_species) +
        xlim(0, 150) +
        ylab("Effect size") +
        facet_wrap(~Test, ncol = 1) +
        my_theme

    p3 <- ggplot(effect_size_long, aes(x = Diversity, y = Effect_size, color = Test)) +
        geom_point(points_aes, size = 2) +
        scale_shape_manual(values = c(21, 20)) +
        stat_poly_line(regression_aes, se = !split_by_species) +
        xlim(0, 15) +
        ylab("Effect size") +
        facet_wrap(~Test, ncol = 1) +
        my_theme

    if (show_R2) {
        p1 <- p1 + stat_poly_eq(use_label(c("R2", "p")), col = "black", size = 5)
        p2 <- p2 + stat_poly_eq(use_label(c("R2", "p")), col = "black", size = 5)
        p3 <- p3 + stat_poly_eq(use_label(c("R2", "p")), col = "black", size = 5)
    }

    do.call("grid.arrange", c(list(p1, p2, p3), ncol = 3))
}

jpeg("Figure 3.jpg", width = 2048, height = 2048, res = 300, quality = 100)

plot_correlation(split_by_species = FALSE, show_R2 = FALSE)
dev.off()

jpeg("Figure 3_by_species.jpg", width = 2048, height = 2048, res = 300, quality = 100)
plot_correlation(split_by_species = TRUE, show_R2 = FALSE)
dev.off()

# Fit the models

models <- lapply(unique(effect_size_long$Test), function(x) {
    #  We only have mice in the TST test
    if (x != "TST") {
        effect_size_long %>%
            subset(Test == x & !is.na(Effect_size)) %>%
            subset(complete.cases(.)) %>%
            lm(Effect_size ~ Duration * Burden * Diversity + Species, data = .)
    } else {
        effect_size_long %>%
            subset(Test == x & !is.na(Effect_size)) %>%
            subset(complete.cases(.)) %>%
            lm(Effect_size ~ Duration * Burden * Diversity, data = .)
    }
})

names(models) <- unique(effect_size_long$Test)

for (i in seq_along(models)) {
    m <- models[[i]]
    print(summary(m))

    # Diagnostic plots
    # Residuals vs fitted
    # p1 <- ggplot(m, aes(.fitted, .resid)) +
    #     geom_point() +
    #     geom_hline(yintercept = 0, linetype = "dashed") +
    #     xlab("Fitted values") +
    #     ylab("Residuals") +
    #     ggtitle(paste("Test: ", names(models)[i])) +
    #     my_theme +
    #     theme(plot.title = element_text(size = 20))
    # # Residuals QQ plot
    # p2 <- ggplot(m, aes(sample = .stdresid)) +
    #     geom_qq() +
    #     geom_abline() +
    #     ggtitle(paste("Test: ", names(models)[i])) +
    #     my_theme +
    #     theme(plot.title = element_text(size = 20))

    # grid.arrange(p1, p2, ncol = 2)
}

# Figure 4 - Significant predictors of effect size
significant_interactions <- NULL

for (i in seq_along(models)) {
    m <- models[[i]]
    test <- names(models)[i]
    pvals <- coef(summary(m))[, "Pr(>|t|)"]
    pvals <- pvals[pvals < 0.05]
    pvals <- pvals[!grepl("Intercept", names(pvals))]
    factors <- unique(unlist(strsplit(names(pvals), ":")))
    if (length(pvals)) {
        significant_interactions <- rbind(
            significant_interactions,
            data.frame(
                Test = test,
                Factors = factors
            )
        )
    }
}

significant_interactions$Test <- factor(significant_interactions$Test, levels = names(models))

jpeg("Figure 4.jpg", width = 2048, height = 2048, res = 300, quality = 100)
ggplot(significant_interactions, aes(x = Factors, y = Test)) +
    geom_point(size = 10) +
    # Make sure all tests are shown, even if they have no significant interactions
    scale_y_discrete(limits = levels(significant_interactions$Test)) +
    xlab("") +
    ylab("") +
    my_theme +
    theme(axis.text = element_text(angle = 45, hjust = 1, size = 20))
dev.off()

### Figure 5 - Distance in the duration/burden/diversity space vs. difference in effect size

test <- "FST"
effect_size_long_test <- effect_size_long %>%
    filter(Test == test) %>%
    filter(complete.cases(.)) %>%
    mutate(Reference = make.unique(Reference))

# Get the distance matrix in the duration/burden/diversity space
distance_matrix <- effect_size_long_test %>%
    select(Duration, Burden, Diversity) %>%
    dist() %>%
    as.matrix()
rownames(distance_matrix) <- effect_size_long_test$Reference
colnames(distance_matrix) <- effect_size_long_test$Reference

pheatmap(distance_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE
)

distance_matrix[lower.tri(distance_matrix)] <- NA

distance_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "Paper_1") %>%
    pivot_longer(-Paper_1, names_to = "Paper_2", values_to = "Distance") %>%
    filter(!is.na(Distance)) %>%
    mutate(
        Duration_1 = effect_size_long_test$Duration[match(Paper_1, effect_size_long_test$Reference)],
        Duration_2 = effect_size_long_test$Duration[match(Paper_2, effect_size_long_test$Reference)],
        Burden_1 = effect_size_long_test$Burden[match(Paper_1, effect_size_long_test$Reference)],
        Burden_2 = effect_size_long_test$Burden[match(Paper_2, effect_size_long_test$Reference)],
        Diversity_1 = effect_size_long_test$Diversity[match(Paper_1, effect_size_long_test$Reference)],
        Diversity_2 = effect_size_long_test$Diversity[match(Paper_2, effect_size_long_test$Reference)],
        Species_1 = effect_size_long_test$Species[match(Paper_1, effect_size_long_test$Reference)],
        Species_2 = effect_size_long_test$Species[match(Paper_2, effect_size_long_test$Reference)]
    ) %>%
    mutate(
        Effect_size_P1 = effect_size_long_test$Effect_size[match(Paper_1, effect_size_long_test$Reference)],
        Effect_size_P2 = effect_size_long_test$Effect_size[match(Paper_2, effect_size_long_test$Reference)]
    ) %>%
    mutate(Effect_size_diff = abs(Effect_size_P1 - Effect_size_P2)) -> distance_data


jpeg("Figure 5.jpg", width = 3000, height = 3000, res = 300, quality = 100)

distance_data %>%
    ggplot(aes(x = Distance, y = Effect_size_diff)) +
    # geom_point(aes(col = Species_1 == Species_2)) +
    geom_point() +
    # Annotate a few points
    geom_text_repel(
        data = distance_data %>%
            filter(str_detect(Paper_1, "Eid RS") & str_detect(Paper_2, "Yu X") |
                str_detect(Paper_1, "Walton NL") & str_detect(Paper_2, "Qin L") |
                str_detect(Paper_1, "Wang F") & str_detect(Paper_2, "Verma H")),
        aes(label = paste(
            "Du:", Duration_1, "/", Duration_2, "\n",
            "B:", Burden_1, "/", Burden_2, "\n",
            "Di:", Diversity_1, "/", Diversity_2, "\n",
            "ES:", format(Effect_size_P1, digits = 2),
            "/", format(Effect_size_P2, digits = 2)
        )),
        nudge_x = 1, nudge_y = 1.5, size = 4, color = rgb(.4, .4, .4)
    ) +
    annotate("text",
        x = c(3, 89, 128),
        y = c(4.5, 5.9, 0.6),
        label = c("A", "B", "C"),
        size = 8, color = rgb(.7, .2, .2)
    ) +
    xlab("Distance in the Duration/Burden/Diversity space") +
    ylab("Difference in effect size") +
    my_theme
dev.off()

jpeg("Figure 5 not annotated.jpg", width = 3000, height = 3000, res = 300, quality = 100)

distance_data %>%
    ggplot(aes(x = Distance, y = Effect_size_diff)) +
    # geom_point(aes(col = Species_1 == Species_2)) +
    geom_point() +
    xlab("Distance in the Duration/Burden/Diversity space") +
    ylab("Difference in effect size") +
    my_theme
dev.off()

cor.test(distance_data$Distance, distance_data$Effect_size_diff,
    method = "pearson"
)

# Figure 6 - Correlations between tests
jpeg("Figure 6.jpg", width = 2048, height = 2048, res = 300, quality = 100)
effect_size %>%
    select(contains("Effect_size")) %>%
    # Rename columns to remove the "Effect_size" part
    rename_all(~ gsub("_Effect_size", "", .)) %>%
    cor(use = "pairwise.complete.obs", method = "pearson") %>%
    corrplot.mixed(
        upper = "ellipse", lower = "number",
        order = "hclust", upper.col = COL2("PRGn"),
        lower.col = "black", tl.col = "black"
    )
dev.off()


distance_data %>%
    filter(Distance > 80 &
        Effect_size_diff > 5 & Effect_size_diff < 6) %>%
    pull(Paper_1, Paper_2)
