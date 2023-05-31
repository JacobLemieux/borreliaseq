heatMapPAMplas <- function(p, data, col_colours = "black", null_colour = "grey90",
                       presence_label = "Present", absence_label = "Absent",
                       border_colour = "white", cluster_cols = FALSE,
                       cluster_method = "single", cluster_distance = "binary",
                       rev_cols = FALSE, colnames = TRUE,
                       colnames_position = "bottom", colnames_angle = 0,
                       colnames_level = NULL, set_label_colours = FALSE,
                       colnames_offset_x = 0, colnames_offset_y = 0,
                       font.size = 4, hjust = 0.5, offset = 0, width = 1,
                       show_legend = FALSE) {
  # The first two packages are dependencies of the package ggtree.
  require(ggplot2)
  require(tidyr)  # for the function gather
  require(tibble)
  require(magrittr)  # for operators "%<>%" and "%>%"(github.com/GuangchuangYu/ggtree/blob/master/R/operator.R)
  require(ggtree)
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert relative width of the PAM to width of each cell
  width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x, na.rm = TRUE) + offset
  
  # Column-wise clustering
  if (cluster_cols && ncol(data) > 2) {
    hc <- hclust(d = dist(t(data), method = cluster_distance),
                 method = cluster_method)
    data <- data[, hc$order]  # reorder columns
    clustered <- TRUE
  } else {
    clustered <- FALSE
  }
  
  # Inverse columns for improving visualisation
  if (rev_cols) {
    data <- data[, ncol(data) : 1]
  }
  
  # weight cells before converting the PAM (dd) into a data frame
  n_colours <- length(col_colours)
  is_multi_colour <- n_colours > 1 && !is.null(names(col_colours))
  if (is_multi_colour) {
    colours_uniq <- sort(unique(as.character(col_colours)), decreasing = FALSE)  # colours for positive values in PAM. Column names are discarded here.
    colour_codes <- 1 : length(colours_uniq)  # codes for colours
    names(colour_codes) <- colours_uniq  # e.g., colour_codes becomes c("red" = 1, "blue" = 2, ...)
    column_names <- colnames(data)  # Since matrix product loses column names (allele names), we need to save them beforehand.
    if (clustered) {
      col_colours <- col_colours[column_names]  # rearrange colour codes when columns of the PAM have been clustered
    } else if (n_colours > length(column_names))  {
      col_colours <- col_colours[column_names]
    }
    data <- data %*% diag(as.integer(colour_codes[as.character(col_colours)]))  # convert colour characters into integer codes
    colnames(data) <- column_names
    names(colours_uniq) <- as.character(colour_codes)
    colours_uniq <- append(c("0" = null_colour), colours_uniq)
    
    # Convert values in the matrix to factors so as to map colours to the levels
    dd <- as.data.frame(data)
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, as.integer(colour_codes)),
                        labels = c("0", as.character(colour_codes)))
    }
  } else {
    dd <- as.data.frame(data)
    
    # Here, I convert binary values to factors in order to show the legend as binary rather than continuous.
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, 1), labels = c("0", "1"))
    }
  }
  
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  
  # tibble::set_tidy_names(dd) solves the problem of the error: "Can't bind data because some arguments have the same name"
  # see https://github.com/tidyverse/tidyr/issues/472
  dd <- gather(data = tibble::set_tidy_names(dd), key = variable, value = value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  
  # Create a data frame for label attributes
  mapping <- data.frame(from = as.character(dd$variable), to = V2, stringsAsFactors = FALSE)  # from: label texts
  mapping <- unique(mapping)
  paint_labels <- set_label_colours && is_multi_colour
  if (paint_labels) {
    mapping$col <- as.character(col_colours[mapping$from])  # variable label colour
    mapping$col <- factor(mapping$col, levels = as.character(colours_uniq),
                          labels = names(colours_uniq))
  }
  mapping$from = as.factor(mapping$from)  # change back to factors
  
  # Create the coloured tile matrix
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  if (is.null(border_colour)) {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        inherit.aes = FALSE)
  } else {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        colour = border_colour, inherit.aes = FALSE)
  }
  
  # Print column names
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    if (paint_labels) {
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from, colour = col),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    } else {  # Labels are printed in black.
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    }
  }
  
  # Scale colours for both tiles and column labels
  if (is_multi_colour) {
    # This command does not scale label colours when paint_labels = FALSE.
    # Class: title of the legend
    # scale_fill_manual sets tile colours; scale_colour_manual sets label colours
    p2 <- p2 + scale_fill_manual(name = "Class",
                                 breaks = c("0", as.character(colour_codes)),
                                 values = colours_uniq, na.value = NA) +
      scale_colour_manual(name = "Class",
                          breaks = c("0", as.character(colour_codes)),
                          values = colours_uniq, na.value = NA)
  } else {
    p2 <- p2 + scale_fill_manual(name = "Plasmid", breaks = c("1", "0"),
                                 values = c("1" = col_colours, "0" = null_colour),
                                 labels = c(presence_label, absence_label),
                                 na.value = NA)
    #p2 <- p2 + scale_fill_gradient(low = null_colour, high = col_colours, na.value = NA)  # reserve this command for continuous variable
  }
  
  # Print legend
  if (show_legend) {
    p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  }
  
  attr(p2, "mapping") <- mapping
  
  return(list(p = p2, mapping = mapping, data = dd))
}

heatMapPAMpf32 <- function(p, data, col_colours = "black", null_colour = "grey90",
                           presence_label = "Present", absence_label = "Absent",
                           border_colour = "white", cluster_cols = FALSE,
                           cluster_method = "single", cluster_distance = "binary",
                           rev_cols = FALSE, colnames = TRUE,
                           colnames_position = "bottom", colnames_angle = 0,
                           colnames_level = NULL, set_label_colours = FALSE,
                           colnames_offset_x = 0, colnames_offset_y = 0,
                           font.size = 4, hjust = 0.5, offset = 0, width = 1,
                           show_legend = FALSE) {
  # The first two packages are dependencies of the package ggtree.
  require(ggplot2)
  require(tidyr)  # for the function gather
  require(tibble)
  require(magrittr)  # for operators "%<>%" and "%>%"(github.com/GuangchuangYu/ggtree/blob/master/R/operator.R)
  require(ggtree)
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert relative width of the PAM to width of each cell
  width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x, na.rm = TRUE) + offset
  
  # Column-wise clustering
  if (cluster_cols && ncol(data) > 2) {
    hc <- hclust(d = dist(t(data), method = cluster_distance),
                 method = cluster_method)
    data <- data[, hc$order]  # reorder columns
    clustered <- TRUE
  } else {
    clustered <- FALSE
  }
  
  # Inverse columns for improving visualisation
  if (rev_cols) {
    data <- data[, ncol(data) : 1]
  }
  
  # weight cells before converting the PAM (dd) into a data frame
  n_colours <- length(col_colours)
  is_multi_colour <- n_colours > 1 && !is.null(names(col_colours))
  if (is_multi_colour) {
    colours_uniq <- sort(unique(as.character(col_colours)), decreasing = FALSE)  # colours for positive values in PAM. Column names are discarded here.
    colour_codes <- 1 : length(colours_uniq)  # codes for colours
    names(colour_codes) <- colours_uniq  # e.g., colour_codes becomes c("red" = 1, "blue" = 2, ...)
    column_names <- colnames(data)  # Since matrix product loses column names (allele names), we need to save them beforehand.
    if (clustered) {
      col_colours <- col_colours[column_names]  # rearrange colour codes when columns of the PAM have been clustered
    } else if (n_colours > length(column_names))  {
      col_colours <- col_colours[column_names]
    }
    data <- data %*% diag(as.integer(colour_codes[as.character(col_colours)]))  # convert colour characters into integer codes
    colnames(data) <- column_names
    names(colours_uniq) <- as.character(colour_codes)
    colours_uniq <- append(c("0" = null_colour), colours_uniq)
    
    # Convert values in the matrix to factors so as to map colours to the levels
    dd <- as.data.frame(data)
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, as.integer(colour_codes)),
                        labels = c("0", as.character(colour_codes)))
    }
  } else {
    dd <- as.data.frame(data)
    
    # Here, I convert binary values to factors in order to show the legend as binary rather than continuous.
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, 1), labels = c("0", "1"))
    }
  }
  
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  
  # tibble::set_tidy_names(dd) solves the problem of the error: "Can't bind data because some arguments have the same name"
  # see https://github.com/tidyverse/tidyr/issues/472
  dd <- gather(data = tibble::set_tidy_names(dd), key = variable, value = value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  
  # Create a data frame for label attributes
  mapping <- data.frame(from = as.character(dd$variable), to = V2, stringsAsFactors = FALSE)  # from: label texts
  mapping <- unique(mapping)
  paint_labels <- set_label_colours && is_multi_colour
  if (paint_labels) {
    mapping$col <- as.character(col_colours[mapping$from])  # variable label colour
    mapping$col <- factor(mapping$col, levels = as.character(colours_uniq),
                          labels = names(colours_uniq))
  }
  mapping$from = as.factor(mapping$from)  # change back to factors
  
  # Create the coloured tile matrix
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  if (is.null(border_colour)) {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        inherit.aes = FALSE)
  } else {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        colour = border_colour, inherit.aes = FALSE)
  }
  
  # Print column names
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    if (paint_labels) {
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from, colour = col),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    } else {  # Labels are printed in black.
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    }
  }
  
  # Scale colours for both tiles and column labels
  if (is_multi_colour) {
    # This command does not scale label colours when paint_labels = FALSE.
    # Class: title of the legend
    # scale_fill_manual sets tile colours; scale_colour_manual sets label colours
    p2 <- p2 + scale_fill_manual(name = "Class",
                                 breaks = c("0", as.character(colour_codes)),
                                 values = colours_uniq, na.value = NA) +
      scale_colour_manual(name = "Class",
                          breaks = c("0", as.character(colour_codes)),
                          values = colours_uniq, na.value = NA)
  } else {
    p2 <- p2 + scale_fill_manual(name = "PF-32", breaks = c("1", "0"),
                                 values = c("1" = col_colours, "0" = null_colour),
                                 labels = c(presence_label, absence_label),
                                 na.value = NA)
    #p2 <- p2 + scale_fill_gradient(low = null_colour, high = col_colours, na.value = NA)  # reserve this command for continuous variable
  }
  
  # Print legend
  if (show_legend) {
    p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  }
  
  attr(p2, "mapping") <- mapping
  
  return(list(p = p2, mapping = mapping, data = dd))
}


heatMapPAM <- function(p, data, col_colours = "black", null_colour = "grey90",
                           presence_label = "Present", absence_label = "Absent",
                           border_colour = "white", cluster_cols = FALSE,
                           cluster_method = "single", cluster_distance = "binary",
                           rev_cols = FALSE, colnames = TRUE,
                           colnames_position = "bottom", colnames_angle = 0,
                           colnames_level = NULL, set_label_colours = FALSE,
                           colnames_offset_x = 0, colnames_offset_y = 0,
                           font.size = 4, hjust = 0.5, offset = 0, width = 1,
                           show_legend = FALSE) {
  # The first two packages are dependencies of the package ggtree.
  require(ggplot2)
  require(tidyr)  # for the function gather
  require(tibble)
  require(magrittr)  # for operators "%<>%" and "%>%"(github.com/GuangchuangYu/ggtree/blob/master/R/operator.R)
  require(ggtree)
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert relative width of the PAM to width of each cell
  width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x, na.rm = TRUE) + offset
  
  # Column-wise clustering
  if (cluster_cols && ncol(data) > 2) {
    hc <- hclust(d = dist(t(data), method = cluster_distance),
                 method = cluster_method)
    data <- data[, hc$order]  # reorder columns
    clustered <- TRUE
  } else {
    clustered <- FALSE
  }
  
  # Inverse columns for improving visualisation
  if (rev_cols) {
    data <- data[, ncol(data) : 1]
  }
  
  # weight cells before converting the PAM (dd) into a data frame
  n_colours <- length(col_colours)
  is_multi_colour <- n_colours > 1 && !is.null(names(col_colours))
  if (is_multi_colour) {
    colours_uniq <- sort(unique(as.character(col_colours)), decreasing = FALSE)  # colours for positive values in PAM. Column names are discarded here.
    colour_codes <- 1 : length(colours_uniq)  # codes for colours
    names(colour_codes) <- colours_uniq  # e.g., colour_codes becomes c("red" = 1, "blue" = 2, ...)
    column_names <- colnames(data)  # Since matrix product loses column names (allele names), we need to save them beforehand.
    if (clustered) {
      col_colours <- col_colours[column_names]  # rearrange colour codes when columns of the PAM have been clustered
    } else if (n_colours > length(column_names))  {
      col_colours <- col_colours[column_names]
    }
    data <- data %*% diag(as.integer(colour_codes[as.character(col_colours)]))  # convert colour characters into integer codes
    colnames(data) <- column_names
    names(colours_uniq) <- as.character(colour_codes)
    colours_uniq <- append(c("0" = null_colour), colours_uniq)
    
    # Convert values in the matrix to factors so as to map colours to the levels
    dd <- as.data.frame(data)
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, as.integer(colour_codes)),
                        labels = c("0", as.character(colour_codes)))
    }
  } else {
    dd <- as.data.frame(data)
    
    # Here, I convert binary values to factors in order to show the legend as binary rather than continuous.
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, 1), labels = c("0", "1"))
    }
  }
  
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  
  # tibble::set_tidy_names(dd) solves the problem of the error: "Can't bind data because some arguments have the same name"
  # see https://github.com/tidyverse/tidyr/issues/472
  dd <- gather(data = tibble::set_tidy_names(dd), key = variable, value = value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  
  # Create a data frame for label attributes
  mapping <- data.frame(from = as.character(dd$variable), to = V2, stringsAsFactors = FALSE)  # from: label texts
  mapping <- unique(mapping)
  paint_labels <- set_label_colours && is_multi_colour
  if (paint_labels) {
    mapping$col <- as.character(col_colours[mapping$from])  # variable label colour
    mapping$col <- factor(mapping$col, levels = as.character(colours_uniq),
                          labels = names(colours_uniq))
  }
  mapping$from = as.factor(mapping$from)  # change back to factors
  
  # Create the coloured tile matrix
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  if (is.null(border_colour)) {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        inherit.aes = FALSE)
  } else {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        colour = border_colour, inherit.aes = FALSE)
  }
  
  # Print column names
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    if (paint_labels) {
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from, colour = col),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    } else {  # Labels are printed in black.
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    }
  }
  
  # Scale colours for both tiles and column labels
  if (is_multi_colour) {
    # This command does not scale label colours when paint_labels = FALSE.
    # Class: title of the legend
    # scale_fill_manual sets tile colours; scale_colour_manual sets label colours
    p2 <- p2 + scale_fill_manual(name = "Class",
                                 breaks = c("0", as.character(colour_codes)),
                                 values = colours_uniq, na.value = NA) +
      scale_colour_manual(name = "Class",
                          breaks = c("0", as.character(colour_codes)),
                          values = colours_uniq, na.value = NA)
  } else {
    p2 <- p2 + scale_fill_manual(name = "ORF", breaks = c("1", "0"),
                                 values = c("1" = col_colours, "0" = null_colour),
                                 labels = c(presence_label, absence_label),
                                 na.value = NA)
    #p2 <- p2 + scale_fill_gradient(low = null_colour, high = col_colours, na.value = NA)  # reserve this command for continuous variable
  }
  
  # Print legend
  if (show_legend) {
    p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  }
  
  attr(p2, "mapping") <- mapping
  
  return(list(p = p2, mapping = mapping, data = dd))
}

heatMapPAMlipo <- function(p, data, col_colours = "black", null_colour = "grey90",
                       presence_label = "Present", absence_label = "Absent",
                       border_colour = "white", cluster_cols = FALSE,
                       cluster_method = "single", cluster_distance = "binary",
                       rev_cols = FALSE, colnames = TRUE,
                       colnames_position = "bottom", colnames_angle = 0,
                       colnames_level = NULL, set_label_colours = FALSE,
                       colnames_offset_x = 0, colnames_offset_y = 0,
                       font.size = 4, hjust = 0.5, offset = 0, width = 1,
                       show_legend = FALSE) {
  # The first two packages are dependencies of the package ggtree.
  require(ggplot2)
  require(tidyr)  # for the function gather
  require(tibble)
  require(magrittr)  # for operators "%<>%" and "%>%"(github.com/GuangchuangYu/ggtree/blob/master/R/operator.R)
  require(ggtree)
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert relative width of the PAM to width of each cell
  width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x, na.rm = TRUE) + offset
  
  # Column-wise clustering
  if (cluster_cols && ncol(data) > 2) {
    hc <- hclust(d = dist(t(data), method = cluster_distance),
                 method = cluster_method)
    data <- data[, hc$order]  # reorder columns
    clustered <- TRUE
  } else {
    clustered <- FALSE
  }
  
  # Inverse columns for improving visualisation
  if (rev_cols) {
    data <- data[, ncol(data) : 1]
  }
  
  # weight cells before converting the PAM (dd) into a data frame
  n_colours <- length(col_colours)
  is_multi_colour <- n_colours > 1 && !is.null(names(col_colours))
  if (is_multi_colour) {
    colours_uniq <- sort(unique(as.character(col_colours)), decreasing = FALSE)  # colours for positive values in PAM. Column names are discarded here.
    colour_codes <- 1 : length(colours_uniq)  # codes for colours
    names(colour_codes) <- colours_uniq  # e.g., colour_codes becomes c("red" = 1, "blue" = 2, ...)
    column_names <- colnames(data)  # Since matrix product loses column names (allele names), we need to save them beforehand.
    if (clustered) {
      col_colours <- col_colours[column_names]  # rearrange colour codes when columns of the PAM have been clustered
    } else if (n_colours > length(column_names))  {
      col_colours <- col_colours[column_names]
    }
    data <- data %*% diag(as.integer(colour_codes[as.character(col_colours)]))  # convert colour characters into integer codes
    colnames(data) <- column_names
    names(colours_uniq) <- as.character(colour_codes)
    colours_uniq <- append(c("0" = null_colour), colours_uniq)
    
    # Convert values in the matrix to factors so as to map colours to the levels
    dd <- as.data.frame(data)
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, as.integer(colour_codes)),
                        labels = c("0", as.character(colour_codes)))
    }
  } else {
    dd <- as.data.frame(data)
    
    # Here, I convert binary values to factors in order to show the legend as binary rather than continuous.
    for (i in names(dd)) {
      dd[, i] <- factor(dd[, i], levels = c(0, 1), labels = c("0", "1"))
    }
  }
  
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  
  # tibble::set_tidy_names(dd) solves the problem of the error: "Can't bind data because some arguments have the same name"
  # see https://github.com/tidyverse/tidyr/issues/472
  dd <- gather(data = tibble::set_tidy_names(dd), key = variable, value = value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  
  # Create a data frame for label attributes
  mapping <- data.frame(from = as.character(dd$variable), to = V2, stringsAsFactors = FALSE)  # from: label texts
  mapping <- unique(mapping)
  paint_labels <- set_label_colours && is_multi_colour
  if (paint_labels) {
    mapping$col <- as.character(col_colours[mapping$from])  # variable label colour
    mapping$col <- factor(mapping$col, levels = as.character(colours_uniq),
                          labels = names(colours_uniq))
  }
  mapping$from = as.factor(mapping$from)  # change back to factors
  
  # Create the coloured tile matrix
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  if (is.null(border_colour)) {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        inherit.aes = FALSE)
  } else {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                        colour = border_colour, inherit.aes = FALSE)
  }
  
  # Print column names
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    if (paint_labels) {
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from, colour = col),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    } else {  # Labels are printed in black.
      p2 <- p2 + geom_text(data = mapping,
                           aes(x = to, y = y, label = from),
                           size = font.size, inherit.aes = FALSE,
                           angle = colnames_angle, nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y, hjust = hjust)
    }
  }
  
  # Scale colours for both tiles and column labels
  if (is_multi_colour) {
    # This command does not scale label colours when paint_labels = FALSE.
    # Class: title of the legend
    # scale_fill_manual sets tile colours; scale_colour_manual sets label colours
    p2 <- p2 + scale_fill_manual(name = "Class",
                                 breaks = c("0", as.character(colour_codes)),
                                 values = colours_uniq, na.value = NA,
                                 guide = guide_legend(direction = "horizontal", title.position = "left", title.theme = element_text(angle = 90, size = 20),
                                                      label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, label.theme = element_text(angle = 90, size = 20))) +
      scale_colour_manual(name = "Class",
                          breaks = c("0", as.character(colour_codes)),
                          values = colours_uniq, na.value = NA,
                          guide = guide_legend(title.theme = element_text(angle = 90)))
  } else {
    p2 <- p2 + scale_fill_manual(name = "Surface Lipoprotein", breaks = c("1", "0"),
                                 values = c("1" = col_colours, "0" = null_colour),
                                 labels = c(presence_label, absence_label),
                                 na.value = NA, 
                                 guide = guide_legend(direction = "horizontal", title.position = "left", title.theme = element_text(angle = 90, size = 20),
                                                      label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, label.theme = element_text(angle = 90, size = 20)))
    #p2 <- p2 + scale_fill_gradient(low = null_colour, high = col_colours, na.value = NA)  # reserve this command for continuous variable
  }
  
  # Print legend
  if (show_legend) {
    p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  }
  
  attr(p2, "mapping") <- mapping
  
  return(list(p = p2, mapping = mapping, data = dd))
}