wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(glue)

file <- "louisiana"
R <- 2
table <- read.csv(file = glue("{file}_stats.csv"))

f <- glue("{file}_bias.tex")
Y <- ""
for(r in 1:R) {
  Y <- sprintf("%s & $Y_{%d, 1}$", Y, r)
}
cat(
  "\\begin{table}[!ht]\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{80pt}>{\\raggedright\\arraybackslash}m{60pt}",
  rep("r", R + 1),
  "}\n",
  "\t\t\\hline\n",
  "\t\tMethod & Parameters",
  Y,
  "& Execution time \\\\ \\hline\n",
  file = f
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
      file = f, append = TRUE
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
    file = f, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\\end{table}",
  file = f, append = TRUE
)

# =================================================================================
f <- glue("{file}_bounds.tex")
Y <- ""
Y.2 <- ""
for(r in 1:R) {
  Y <- sprintf("%s & \\SetCell[c=2]{}{$Y_{%d, 1}$ bounds} & ", Y, r)
  Y.2 <- sprintf("%s & Size & Accuracy", Y.2)
}
cat(
  "\\begin{table}[!ht]\n",
  "\t\\centering\n",
  "\t\\begin{tblr}{>{\\raggedright\\arraybackslash}m{80pt}>{\\raggedright\\arraybackslash}m{60pt}",
  rep("r", 2 * R + 2),
  "}\n",
  "\t\t\\hline\n",
  "\t\t\\SetCell[r=2]{}{Method} & \\SetCell[r=2]{}{Parameters}",
  Y,
  "& \\SetCell[c=2]{}{Error} \\\\ \\hline\n",
  "\t\t~&~",
  Y.2,
  "& EI & EQ \\\\ \\hline\n",
  file = f
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
      file = f, append = TRUE
    )
  }
  bounds <- ""
  for(r in 1:R) {
    size <- table[i, glue("Bounds.size..{r}..1.")]
    percent <- 100 * table[i, glue("True.in.bounds..{r}..1.")]
    bounds <- sprintf("%s & $%.3f$ & $%.0f\\%%$", bounds, size, percent)
  }
  EI <- sprintf("%.3f", table[i, "EI"])
  EQ <- sprintf("%.3f", table[i, "EQ"])
  cat(
    "\t\t",
    glue("{name} & {index[2]}{bounds} & {EI} & {EQ}"),
    "\\\\\n",
    file = f, append = TRUE
  )
}

cat(
  "\t\\end{tblr}\n",
  "\\end{table}",
  file = f, append = TRUE
)