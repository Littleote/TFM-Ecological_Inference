wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(glue)

R <- 2
table <- read.csv(file = glue("out/louisiana_stats.csv"), check.names = FALSE)

if (!dir.exists("../latex/taules")) {
  dir.create("../latex/taules", recursive = TRUE)
}

name.map <- list()
name.map[["logit covariate"]] <- "Our model"
name.map[["logit covariate without random effects"]] <- "Our model without random effects"
name.map[["cluster"]] <- "Cluster model"
name.map[["spatial"]] <- "Spatial model"
name.map[["covariate"]] <- "Covariate model"
name.map[["lphom"]] <- "Linear programming model"
name.map[["ecolRxC"]] <- "Latent structure model"

# =================================================================================

file.name <- "../latex/taules/louisiana_bias_time.tex"
Y.index <- ""
for(r in 1:R) {
  Y.index <- sprintf("%s & $Y_{%d, 1}$", Y.index, r)
}
cat(
  "\\begin{table}\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{75pt}>{\\raggedright\\arraybackslash}m{55pt}",
  rep("r", R + 1),
  "}\n",
  "\t\t\\hline\n",
  "\t\tMethod & Model",
  Y.index,
  "& Execution time \\\\ \\hline\n",
  file = file.name
)
last.item <- ""
for(i in 1:nrow(table)) {
  index <- strsplit(table[i, "Model"], ",")[[1]]
  name <- "~"
  if(last.item != index[1]) {
    last.item <- index[1]
    count <- sum(mapply(function(x) { startsWith(x, glue("{last.item},")) } , table[,1]))
    name <- glue("\\SetCell[r={count}]{{}}{name.map[[last.item]]}")
    cat(
      "\t\t\\hline\n",
      file = file.name, append = TRUE
    )
  }
  b.err <- ""
  for(r in 1:R) {
    bias <- table[i, glue("Bias [{r}, 1]")]
    err <- 1.96 * table[i, glue("Deviation [{r}, 1]")]
    b.err <- sprintf("%s & $%.3f\\pm%.3f$", b.err, bias, err)
  }
  time <- sprintf("$%.2f$s", table[i, "Execution Time"])
  cat(
    "\t\t",
    glue("{name} & {index[2]}{b.err} & {time}"),
    "\\\\\n",
    file = file.name, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\t\\caption{Louisiana bias and execution time for each model.}\n",
  "\t\\label{tbl:louisiana-bias-time}\n",
  "\\end{table}",
  file = file.name, append = TRUE
)

# =================================================================================

file.name <- "../latex/taules/louisiana_bounds_error.tex"
Y.index <- ""
Y.bound <- ""
for(r in 1:R) {
  Y.index <- sprintf("%s & \\SetCell[c=2]{}{$Y_{%d, 1}$ bounds} & ", Y.index, r)
  Y.bound <- sprintf("%s & Size & Acc.", Y.bound)
}
cat(
  "\\begin{table}\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{75pt}>{\\raggedright\\arraybackslash}m{55pt}",
  rep("r", 2 * R + 2 * 2),
  "}\n",
  "\t\t\\hline\n",
  "\t\t\\SetCell[r=2]{}{Method} & \\SetCell[r=2]{}{Model}",
  Y.index,
  "& \\SetCell[c=2]{}{Local Error} & & \\SetCell[c=2]{}{Global Error} \\\\ \\hline\n",
  "\t\t~&~",
  Y.bound,
  "& MAE & MSE & MAE & MSE \\\\ \\hline\n",
  file = file.name
)
last.item <- ""
for(i in 1:nrow(table)) {
  index <- strsplit(table[i, 1], ",")[[1]]
  name <- "~"
  if(last.item != index[1]) {
    last.item <- index[1]
    count <- sum(mapply(function(x) { startsWith(x, glue("{last.item},")) } , table[,1]))
    name <- glue("\\SetCell[r={count}]{{}}{name.map[[last.item]]}")
    cat(
      "\t\t\\hline\n",
      file = file.name, append = TRUE
    )
  }
  bounds <- ""
  for(r in 1:R) {
    size <- table[i, glue("Bounds Size [{r}, 1]")]
    percent <- 100 * table[i, glue("Accuracy [{r}, 1]")]
    bounds <- sprintf("%s & $%.3f$ & $%.0f\\%%$", bounds, size, percent)
  }
  local.MAE <- sprintf("$%.3f$", table[i, "Local MAE"])
  local.MSE <- sprintf("$%.3f$", table[i, "Local MSE"])
  global.MAE <- sprintf("$%.3f$", table[i, "Global MAE"])
  global.MSE <- sprintf("$%.3f$", table[i, "Global MSE"])
  cat(
    "\t\t",
    glue("{name} & {index[2]}{bounds} & {local.MAE} & {local.MSE} & {global.MAE} & {global.MSE}"),
    "\\\\\n",
    file = file.name, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\t\\caption{Louisiana bounds and error for each model.}\n",
  "\t\\label{tbl:louisiana-bounds-error}\n",
  "\\end{table}",
  file = file.name, append = TRUE
)

# =================================================================================

df <- read.csv("../data/Louisiana-2020_presidential_election.csv", skip = 1)
true.count <- unlist(df[6, paste0(c("WHITE", "BLACK"), ".", c(1, 1, 2, 2))])
idx <- paste0("Example Counts [", 1:2, ", ", c(1, 1, 2, 2), "]")
format <- function(counts, bf = FALSE) {
  if (bf) {
    return(
      do.call(sprintf, c(
        list("$\\mathbf{%d}$ & $\\mathbf{%d}$ & $\\mathbf{%d}$ & $\\mathbf{%d}$"),
        as.list(counts)
      ))
    )
  }
  do.call(sprintf, c(list("$%d$ & $%d$ & $%d$ & $%d$"), as.list(counts)))
}

file.name <- "../latex/taules/louisiana_counts.tex"
cat(
  "\\begin{table}\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{lrrrr",
  "}\n",
  "\t\t\\hline\n",
  "\t\tModel & $Y_{1,\\, 1, 1}$ & $Y_{1,\\, 2, 1}$ & $Y_{1,\\, 1, 2}$ & $Y_{1,\\, 2, 2}$ \\\\ \\hline\n",
  "\t\t\\textbf{True counts} & ",
  format(true.count, bf = TRUE),
  " \\\\ \n",
  file = file.name
)

for(i in c(4, 5, 10, 11, 12, 13, 20)) {
  name <- strsplit(table[i, "Model"], ", ")[[1]]
  cat(
    "\t\t",
    name.map[[name[1]]],
    ": ",
    name[2],
    " & ",
    format(round(table[i, idx])),
    " \\\\ \n",
    file = file.name, append = TRUE
  )
}

cat(
  "\t\t\\hline\n",
  "\t\\end{tblr}\n",
  "\t\\caption{Counts predicted for Acadia for the best models of each method}\n",
  "\t\\label{tbl:louisiana-example}\n",
  "\\end{table}",
  file = file.name, append = TRUE
)