wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(glue)

R <- 2
table <- read.csv(file = glue("out/louisiana_stats.csv"))

if (!dir.exists("../latex/taules")) {
  dir.create("../latex/taules", recursive = TRUE)
}

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
    name <- glue("\\SetCell[r={count}]{{}}{last.item}")
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
  "\\end{table}",
  file = file.name, append = TRUE
)

# =================================================================================

file.name <- "../latex/taules/louisiana_bounds_error.tex"
Y.index <- ""
Y.bound <- ""
for(r in 1:R) {
  Y.index <- sprintf("%s & \\SetCell[c=2]{}{$Y_{%d, 1}$ bounds} & ", Y.index, r)
  Y.bound <- sprintf("%s & Size & Accuracy", Y.bound)
}
cat(
  "\\begin{table}[!ht]\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{80pt}>{\\raggedright\\arraybackslash}m{60pt}",
  rep("r", 2 * R + 2),
  "}\n",
  "\t\t\\hline\n",
  "\t\t\\SetCell[r=2]{}{Method} & \\SetCell[r=2]{}{Parameters}",
  Y.index,
  "& \\SetCell[c=2]{}{Local Error} & \\SetCell[c=2]{}{Global Error} \\\\ \\hline\n",
  "\t\t~&~",
  Y.bound,
  "& EI & EQ & EI & EQ \\\\ \\hline\n",
  file = file.name
)
last.item <- ""
for(i in 1:nrow(table)) {
  index <- strsplit(table[i, 1], ",")[[1]]
  name <- "~"
  if(last.item != index[1]) {
    last.item <- index[1]
    count <- sum(mapply(function(x) { startsWith(x, glue("{last.item},")) } , table[,1]))
    name <- glue("\\SetCell[r={count}]{{}}{last.item}")
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
  local.EI <- sprintf("%.3f", table[i, "local.EI"])
  local.EQ <- sprintf("%.3f", table[i, "local.EQ"])
  global.EI <- sprintf("%.3f", table[i, "global.EI"])
  global.EQ <- sprintf("%.3f", table[i, "global.EQ"])
  cat(
    "\t\t",
    glue("{name} & {index[2]}{bounds} & {local.EI} & {local.EQ} & {global.EI} & {global.EQ}"),
    "\\\\\n",
    file = file.name, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\\end{table}",
  file = file.name, append = TRUE
)