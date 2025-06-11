wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(glue)

R <- 2
table <- read.csv(file = glue("out/louisiana_stats.csv"))

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
name.map[["ecolRxC"]] <- "Latent space model"

# =================================================================================

file.name <- "../latex/taules/louisiana_bias_time.tex"
Y.index <- ""
for(r in 1:R) {
  Y.index <- sprintf("%s & $Y_{%d, 1}$", Y.index, r)
}
cat(
  "\\begin{table}[!ht]\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{80pt}>{\\raggedright\\arraybackslash}m{60pt}",
  rep("r", R + 1),
  "}\n",
  "\t\t\\hline\n",
  "\t\tMethod & Parameters",
  Y.index,
  "& Execution time \\\\ \\hline\n",
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
  b.err <- ""
  for(r in 1:R) {
    bias <- table[i, glue("Individual.bias..{r}..1.")]
    err <- 1.96 * table[i, glue("Individual.deviation..{r}..1.")]
    b.err <- sprintf("%s & $%.3f\\pm%.3f$", b.err, bias, err)
  }
  time <- sprintf("%.2fs", table[i, "Execution.time"])
  cat(
    "\t\t",
    glue("{name} & {index[2]}{b.err} & {time}"),
    "\\\\\n",
    file = file.name, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\t\\caption{Louisiana bias and execution time.}\n",
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
  "\\begin{table}[!ht]\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{60pt}>{\\raggedright\\arraybackslash}m{60pt}",
  rep("r", 2 * R + 2 * 2),
  "}\n",
  "\t\t\\hline\n",
  "\t\t\\SetCell[r=2]{}{Method} & \\SetCell[r=2]{}{Parameters}",
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
    size <- table[i, glue("Bounds.size..{r}..1.")]
    percent <- 100 * table[i, glue("True.in.bounds..{r}..1.")]
    bounds <- sprintf("%s & $%.3f$ & $%.0f\\%%$", bounds, size, percent)
  }
  local.MAE <- sprintf("%.3f", table[i, "local.MAE"])
  local.MSE <- sprintf("%.3f", table[i, "local.MSE"])
  global.MAE <- sprintf("%.3f", table[i, "global.MAE"])
  global.MSE <- sprintf("%.3f", table[i, "global.MSE"])
  cat(
    "\t\t",
    glue("{name} & {index[2]}{bounds} & {local.MAE} & {local.MSE} & {global.MAE} & {global.MSE}"),
    "\\\\\n",
    file = file.name, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\t\\caption{Louisiana bounds and error.}\n",
  "\t\\label{tbl:louisiana-bounds-error}\n",
  "\\end{table}",
  file = file.name, append = TRUE
)