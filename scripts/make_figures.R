# Regenerates every figure and data file the portfolio site needs, from the
# underlying models and datasets. Nothing on the site is hand-drawn.
#
#   Rscript scripts/make_figures.R
#
# Raw surveillance data is read from a local path and never committed; only the
# derived aggregates written into assets/ ship with the repository.

suppressPackageStartupMessages({
  library(deSolve); library(ggplot2); library(dplyr); library(tidyr)
})

ROOT <- if (dir.exists("assets")) normalizePath(".") else "C:/Users/Zander/Claude Code/portfolio"
OUT  <- file.path(ROOT, "assets", "img")
DAT  <- file.path(ROOT, "assets", "data")
DATA <- "C:/Users/Zander/Downloads/Files/06_Data-Spreadsheets"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DAT, showWarnings = FALSE, recursive = TRUE)

ink <- "#12233a"; accent <- "#c2410c"
base <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 13, colour = ink),
        plot.subtitle = element_text(size = 10, colour = "#5b6b7f"),
        axis.title    = element_text(size = 10, colour = "#5b6b7f"),
        legend.position = "top", legend.title = element_blank())

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, name), p, width = w, height = h, dpi = 150, bg = "white")
  cat("wrote", name, "\n")
}

# =========================================================================
# 1. SIR, parameters exactly as in basicSIRsimulation.R
# =========================================================================
sir <- function(t, y, p) with(as.list(c(y, p)),
  list(c(-beta*S*I, beta*S*I - gamma*I, gamma*I)))

out <- ode(c(S = 999, I = 1, R = 0), seq(0, 100, 1), sir,
           c(beta = 0.0003, gamma = 0.1)) |> as.data.frame()

p1 <- out |> pivot_longer(-time) |>
  mutate(name = factor(name, c("S","I","R"), c("Susceptible","Infected","Recovered"))) |>
  ggplot(aes(time, value, colour = name)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("#2563eb", accent, "#059669")) +
  labs(title = "SIR model with mass-action incidence",
       subtitle = "beta = 0.0003, gamma = 0.1, N = 1000",
       x = "Day", y = "Individuals") + base
save_fig(p1, "fig-sir.png", 7, 4.2)

# =========================================================================
# 2. SEIR, structure as in ShinyApp_SEIRmodel.R
# =========================================================================
seir <- function(t, y, p) with(as.list(c(y, p)),
  list(c(-beta*S*E, beta*S*E - sigma*E, sigma*E - gamma*I, gamma*I)))

out2 <- ode(c(S = 999, E = 1, I = 0, R = 0), seq(0, 150, 1), seir,
            c(beta = 0.0004, sigma = 0.2, gamma = 0.1)) |> as.data.frame()

p2 <- out2 |> pivot_longer(-time) |>
  mutate(name = factor(name, c("S","E","I","R"),
                       c("Susceptible","Exposed","Infected","Recovered"))) |>
  ggplot(aes(time, value, colour = name)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("#2563eb", "#d97706", accent, "#059669")) +
  labs(title = "SEIR model with a latent compartment",
       subtitle = "The exposed class delays and flattens the infection peak",
       x = "Day", y = "Individuals") + base
save_fig(p2, "fig-seir.png", 7, 4.2)

# =========================================================================
# 3. Cholera SEIRB, equations and defaults taken from Final_Cholera_Shiny.R
#    Environmental transmission through a bacterial reservoir B.
# =========================================================================
cholera <- function(t, y, p) with(as.list(c(y, p)), {
  inf_env    <- omega * alpha1 * S * B / K
  inf_direct <- phi * alpha2 * S * I
  list(c(
    Lambda - inf_env - inf_direct - (nu + mu) * S,
    inf_env + inf_direct - (kappa + mu) * E,
    kappa * E - (gamma + delta + mu) * I,
    nu * S + gamma * I - mu * R,
    epsilon * I - (tau + eta) * B
  ))
})

cp <- c(Lambda = 0.00913 * 1007, mu = 5.48e-5, nu = 0.02, eta = 4,
        kappa = 0.7143, gamma = 0.2, delta = 0.02, omega = 0.00162,
        alpha1 = 0.00015, phi = 0.0073, alpha2 = 0.00096, K = 500,
        epsilon = 10, tau = 0.033)

out3 <- ode(c(S = 1000, E = 5, I = 2, R = 0, B = 20),
            seq(0, 400, 1), cholera, cp) |> as.data.frame()

p3 <- out3 |> select(time, S, E, I, R) |> pivot_longer(-time) |>
  mutate(name = factor(name, c("S","E","I","R"),
                       c("Susceptible","Exposed","Infectious","Recovered"))) |>
  ggplot(aes(time, value, colour = name)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("#2563eb", "#d97706", accent, "#059669")) +
  labs(title = "Cholera SEIRB model with an environmental reservoir",
       subtitle = "Transmission by contaminated water as well as direct contact",
       x = "Day", y = "Individuals") + base
save_fig(p3, "fig-cholera.png", 7, 4.2)

# =========================================================================
# 4. Kuramoto synchronisation transition
#    Order parameter r against coupling strength K.
# =========================================================================
set.seed(42)
kuramoto_r <- function(Kc, n = 400, dt = 0.05, steps = 3000, burn = 2000) {
  omega <- rcauchy(n, 0, 0.5)
  theta <- runif(n, 0, 2 * pi)
  acc <- numeric(0)
  for (s in seq_len(steps)) {
    z <- mean(exp(1i * theta))
    theta <- theta + dt * (omega + Kc * Mod(z) * sin(Arg(z) - theta))
    if (s > burn) acc <- c(acc, Mod(z))
  }
  mean(acc)
}

Ks <- seq(0, 4, by = 0.2)
kur <- data.frame(K = Ks, r = vapply(Ks, kuramoto_r, numeric(1)))

p4 <- ggplot(kur, aes(K, r)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#94a3b8") +
  annotate("text", x = 1.08, y = .07, label = "critical coupling",
           hjust = 0, size = 3, colour = "#64748b") +
  geom_line(colour = accent, linewidth = 1.1) +
  geom_point(colour = accent, size = 2) +
  ylim(0, 1) +
  labs(title = "Kuramoto transition to collective synchronisation",
       subtitle = "400 oscillators with Cauchy-distributed natural frequencies",
       x = "Coupling strength K", y = "Order parameter r") + base
save_fig(p4, "fig-kuramoto.png", 7, 4.2)

# =========================================================================
# 5. Kaprekar triples, counted directly
#    k is an n-Kaprekar triple when k^3 = p*N^2 + q*N + r with k = p+q+r,
#    N = 10^n, 0 <= q,r < N and p > 0. Cubes stay exact below 2^53, so this
#    brute force is trustworthy up to n = 5.
# =========================================================================
count_triples <- function(n) {
  N <- 10^n
  # k must lie strictly below N; k = N satisfies the split trivially with
  # q = r = 0 and is not counted as a triple.
  k <- 1:(N - 1)
  cube <- as.numeric(k)^3
  p <- floor(cube / N^2)
  rem <- cube - p * N^2
  q <- floor(rem / N)
  r <- rem - q * N
  hit <- (p > 0) & (p + q + r == k)
  list(count = sum(hit), members = k[hit])
}

kap <- lapply(1:5, count_triples)
kdf <- data.frame(n = 1:5, count = vapply(kap, \(x) x$count, numeric(1)))
for (i in 1:5) cat(sprintf("  K(10^%d) = %d  %s\n", i, kap[[i]]$count,
                           paste(head(kap[[i]]$members, 6), collapse = ", ")))

# ---- the structure worth showing: everything depends on n mod 81 ---------
# Class counts come from kt_period_table() in the thesis PARI/GP program. That
# program is re-run here when the interpreter is present, and the committed CSV
# is used otherwise, so the figure never silently drifts from the source.
GPDIR <- "C:/Users/Zander/Downloads/Files/01_School/Thesis-Kaprekar/05_computations/pari-gp"
GPEXE <- file.path(GPDIR, "gp64-2-17-4.exe")
CLSCSV <- file.path(DAT, "kaprekar-classes.csv")

if (file.exists(GPEXE)) {
  ok <- tryCatch({
    owd <- setwd(GPDIR)
    raw <- system2(GPEXE, c("-q", "-f", "kaprekar_triples.gp"),
                   stdout = TRUE, stderr = TRUE, input = "quit;")
    setwd(owd)
    rowsx <- grep("^ +[0-9]+ : [0-9]+", raw, value = TRUE)
    # Empty classes print a trailing "<-- K(10^n) EMPTY" annotation, so the two
    # integers are captured explicitly rather than split off the line. Splitting
    # on the colon swallows the annotation into the count and silently yields NA.
    m <- regmatches(rowsx, regexec("^ +([0-9]+) +: +([0-9]+)", rowsx))
    got <- vapply(m, length, integer(1)) == 3
    parsed <- if (all(got)) {
      data.frame(class = as.integer(vapply(m, `[`, character(1), 2)),
                 count = as.integer(vapply(m, `[`, character(1), 3)))
    } else NULL
    if (length(rowsx) == 54 && !is.null(parsed) && !anyNA(parsed)) {
      write.csv(parsed, CLSCSV, row.names = FALSE)
      cat(sprintf("  kaprekar classes: re-derived from PARI/GP (%d classes, %d empty)\n",
                  nrow(parsed), sum(parsed$count == 0)))
      TRUE
    } else FALSE
  }, error = function(e) FALSE)
  if (!ok) cat("  kaprekar classes: PARI/GP run did not return 54 rows, using committed CSV\n")
}

cls <- read.csv(CLSCSV, stringsAsFactors = FALSE)
grid <- data.frame(class = 1:81) |>
  left_join(cls, by = "class") |>
  mutate(col = (class - 1) %% 9 + 1,
         row = (class - 1) %/% 9 + 1,
         excluded = class %% 3 == 0,
         label = ifelse(excluded, "", as.character(count)))

n_empty <- sum(grid$count == 0, na.rm = TRUE)
peaks   <- grid$class[which(grid$count == max(grid$count, na.rm = TRUE))]
cat(sprintf("  kaprekar classes: %d covered, %d empty, max %d at n = %s (mod 81)\n",
            sum(!grid$excluded), n_empty, max(grid$count, na.rm = TRUE),
            paste(peaks, collapse = " and ")))

p5 <- ggplot(grid, aes(col, -row)) +
  geom_tile(aes(fill = count), colour = "white", linewidth = 1.4) +
  # width and height must be stated: with only two highlighted cells ggplot
  # infers the tile size from the data resolution and draws them eight rows tall.
  geom_tile(data = subset(grid, class %in% peaks),
            fill = NA, colour = ink, linewidth = 1.1, width = 1, height = 1) +
  geom_text(aes(label = label, colour = count > 2), size = 3.5, fontface = "bold") +
  geom_text(data = subset(grid, excluded), aes(label = "3|n"),
            size = 2.5, colour = "#9aa7b5") +
  scale_fill_gradient(low = "#fdf1ea", high = "#7c2d12", na.value = "#f1efec",
                      breaks = 0:4, name = "Triples") +
  scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white"), guide = "none") +
  coord_equal() +
  labs(title = "Every power of ten is decided by n mod 81",
       subtitle = sprintf("One cell per residue class. %d classes carry the theorem, %d of them admit no triple at all,\nand the count never exceeds %d, reached only at n = %s. Grey cells fall outside it.",
                          sum(!grid$excluded), n_empty, max(grid$count, na.rm = TRUE),
                          paste(peaks, collapse = " and ")),
       x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, colour = ink),
        plot.subtitle = element_text(size = 9.2, colour = "#5b6b7f", lineheight = 1.25),
        legend.position = "right")
save_fig(p5, "fig-kaprekar.png", 7.6, 5.2)

# =========================================================================
# 6. Dengue surveillance, Baguio City
# =========================================================================
wk <- list.files(DATA, pattern = "^INC PER 1000 - [A-Z]", full.names = TRUE)
mon <- c(Jan=1,Feb=2,Mar=3,Apr=4,May=5,Jun=6,Jul=7,Aug=8,Sep=9,Oct=10,Nov=11,Dec=12)

# The record is not continuous: only 23 of 52 weeks are present. Plotting the
# labels as evenly spaced categories would draw a straight line across a
# four-month hole, so the series is placed on a real date axis and the gaps are
# left as gaps. REFYEAR only fixes weekday spacing; the true year is unconfirmed
# and is never printed.
REFYEAR <- 2023
rows <- lapply(wk, function(f) {
  lab <- sub("[.]csv$", "", sub("^.*INC PER 1000 - ", "", f))
  tok <- strsplit(lab, "[ -]+")[[1]]
  m <- mon[substr(tok[1], 1, 3)]; d <- suppressWarnings(as.integer(tok[2]))
  if (is.na(m) || is.na(d)) return(NULL)
  v <- read.csv(f, stringsAsFactors = FALSE)
  i <- suppressWarnings(as.numeric(v$INC))
  data.frame(label = lab,
             date = as.Date(sprintf("%d-%02d-%02d", REFYEAR, m, d)),
             total = sum(i, na.rm = TRUE),
             affected = sum(i > 0, na.rm = TRUE))
}) |> bind_rows() |> arrange(date)

# split into runs of consecutive weeks so no line is drawn across a gap
rows$run <- cumsum(c(0, as.numeric(diff(rows$date)) > 7))
gaps <- which(as.numeric(diff(rows$date)) > 7)
cat(sprintf("  coverage: %d of 52 weeks (%.0f%%), %d gap(s), longest %.0f weeks\n",
            nrow(rows), nrow(rows) / 52 * 100, length(gaps),
            max(c(0, (as.numeric(diff(rows$date))[gaps] - 7) / 7))))

p6 <- ggplot(rows, aes(date, total)) +
  geom_area(aes(group = run), fill = accent, alpha = .13) +
  geom_line(aes(group = run), colour = accent, linewidth = 1.1) +
  geom_point(colour = accent, size = 1.8) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(title = "Baguio City dengue cases by surveillance week",
       subtitle = sprintf("%d barangays. Only %d of 52 weeks were supplied; the line breaks where the data does.",
                          nrow(read.csv(wk[1])), nrow(rows)),
       x = NULL, y = "Reported cases") +
  base
save_fig(p6, "fig-dengue-weekly.png", 8, 4.4)

brgy <- lapply(wk, function(f) {
  v <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(BRGY = v$BRGY, INC = suppressWarnings(as.numeric(v$INC)))
}) |> bind_rows() |> group_by(BRGY) |>
  summarise(mean_inc = mean(INC, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_inc))

# ---- how concentrated is the burden? -------------------------------------
# A ranked bar chart of anonymised barangays says almost nothing. The question
# worth answering is whether dengue is spread evenly across the city or carried
# by a small number of areas, which is what decides where control effort goes.
per_brgy <- lapply(wk, function(f) {
  v <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(code = v$BRGY.CODE, inc = suppressWarnings(as.numeric(v$INC)))
}) |> bind_rows() |> group_by(code) |>
  summarise(total = sum(inc, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(total))

gini <- function(x) { x <- sort(x); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }

# Route 1: rank by raw case count, x is share of barangays.
cnt <- per_brgy |> arrange(desc(total)) |>
  mutate(x = row_number() / n(), y = cumsum(total) / sum(total))
g_cnt   <- gini(per_brgy$total)
c10     <- cnt$y[which.min(abs(cnt$x - 0.10))]

# Route 2: normalise by residential population and rank by incidence, x is
# share of population. Residential denominators are wrong for the commercial
# core, so those barangays are held out rather than silently distorting it.
POPFILE <- file.path("C:/Users/Zander/Downloads/Files/10_Uncategorized/QGIS",
                     "QGIS Files-20240102T115653Z-001", "QGIS Files", "CSV",
                     "DENGUE_CASES_IN_BAGUIO_CITY - INCIDENCE_RATE.csv")
have_pop <- file.exists(POPFILE)

if (have_pop) {
  pop <- read.csv(POPFILE, stringsAsFactors = FALSE) |>
    transmute(code = CODE, popln = suppressWarnings(as.numeric(POPLN))) |>
    filter(!is.na(popln), popln > 0, nzchar(code))
  jn <- inner_join(per_brgy, pop, by = "code") |> mutate(rate = total / popln * 1000)
  r_all <- cor(jn$total, jn$popln)
  small <- jn |> filter(popln < 1000)
  res <- jn |> filter(popln >= 1000) |> arrange(desc(rate)) |>
    mutate(x = cumsum(popln) / sum(popln), y = cumsum(total) / sum(total))
  g_rate <- gini(jn$rate[jn$popln >= 1000])
  r10    <- res$y[which.min(abs(res$x - 0.10))]
  cat(sprintf("  cases vs population correlation: %.2f (counts are not a population artefact)\n", r_all))
  cat(sprintf("  by count : worst 10%% of barangays  = %.1f%% of cases, Gini %.2f\n", c10 * 100, g_cnt))
  cat(sprintf("  by rate  : worst 10%% of population = %.1f%% of cases, Gini %.2f (%d barangays, popln >= 1000)\n",
              r10 * 100, g_rate, nrow(res)))
  cat(sprintf("  held out : %d barangays under 1,000 residents carrying %.0f%% of cases; max rate there %.0f per 1,000\n",
              nrow(small), sum(small$total) / sum(jn$total) * 100, max(small$rate)))

  lines <- bind_rows(
    cnt |> transmute(x, y, basis = "By case count, share of barangays"),
    res |> transmute(x, y, basis = "By incidence, share of population")
  )
} else {
  cat("  population file not found, plotting the count route only\n")
  lines <- cnt |> transmute(x, y, basis = "By case count, share of barangays")
  r10 <- c10
}

p7 <- ggplot(lines, aes(x, y, colour = basis)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#94a3b8") +
  geom_line(linewidth = 1.2) +
  annotate("segment", x = .10, xend = .10, y = 0, yend = max(c10, r10),
           colour = "#64748b", linetype = "dotted") +
  annotate("text", x = .125, y = max(c10, r10) - .10,
           label = sprintf("worst tenth of the city:\n%.0f%% and %.0f%% of cases", c10 * 100, r10 * 100),
           hjust = 0, size = 3.1, colour = "#475569") +
  scale_colour_manual(values = c(accent, "#0f766e")) +
  scale_x_continuous(labels = function(v) paste0(v * 100, "%")) +
  scale_y_continuous(labels = function(v) paste0(v * 100, "%")) +
  labs(title = "The same concentration, reached two different ways",
       subtitle = "Dashed line is an even spread. Both routes put roughly a third of cases in the worst tenth of the city.",
       x = "Share of the city, ranked worst first", y = "Share of all reported cases") +
  base
save_fig(p7, "fig-dengue-concentration.png", 7.8, 4.8)

# ---- choropleth, deliberately unlabelled ---------------------------------
# Boundaries arrive as CAD linestrings, so they are polygonised and dissolved
# by PSGC code before joining. No barangay is named on the published map.
suppressPackageStartupMessages(library(sf))
BND <- file.path("C:/Users/Zander/Downloads/Files/10_Uncategorized/QGIS",
                 "BAGUIO BARANGAY BOUNDARY-20240101T132904Z-001",
                 "BAGUIO BARANGAY BOUNDARY", "BRGY BOUNDARY.shp")

if (file.exists(BND)) {
  ln <- suppressWarnings(st_read(BND, quiet = TRUE))
  ln <- ln[!is.na(ln$BRGY.code), c("BRGY.code")]
  pg <- suppressWarnings(st_polygonize(st_geometry(ln)))
  keep <- !st_is_empty(pg)
  shp <- st_sf(code = ln$BRGY.code[keep],
               geometry = st_collection_extract(pg[keep], "POLYGON"))
  shp <- shp |> group_by(code) |> summarise(.groups = "drop")

  wks <- nrow(rows)
  mp <- left_join(shp, per_brgy, by = "code") |>
    mutate(rate = total / wks)

  p8 <- ggplot(mp) +
    geom_sf(aes(fill = rate), colour = "white", linewidth = .18) +
    scale_fill_gradient(low = "#fdf1ea", high = "#7c2d12", na.value = "#eceae6",
                        name = "Mean cases\nper week") +
    labs(title = "Dengue burden across Baguio City",
         subtitle = sprintf("%d barangay polygons built from survey boundaries and joined on PSGC code. Names withheld.",
                            nrow(mp))) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13, colour = ink),
          plot.subtitle = element_text(size = 9.5, colour = "#5b6b7f"),
          legend.position = "right")
  save_fig(p8, "fig-dengue-map.png", 7.4, 5.6)
} else {
  cat("  boundary shapefile not found, skipping the map\n")
}

# ---- aggregates for the interactive chart on the site --------------------
# Only these derived numbers leave the machine; the raw files stay put.
jstr <- function(x) paste0('"', gsub('"', '', x), '"')
jnum <- function(x) formatC(x, format = "f", digits = 2, drop0trailing = TRUE)

# Ranked labels only. No barangay name reaches the published data file.
top <- head(brgy, 12) |> mutate(label = sprintf("Barangay %02d", row_number()))
json <- paste0(
  '{\n  "weeks": [', paste(jstr(as.character(rows$label)), collapse = ", "), '],\n',
  '  "cases": [', paste(jnum(rows$total), collapse = ", "), '],\n',
  '  "days": [', paste(as.numeric(rows$date - min(rows$date)), collapse = ", "), '],\n',
  '  "barangays": [', paste(jstr(top$label), collapse = ", "), '],\n',
  '  "means": [', paste(jnum(top$mean_inc), collapse = ", "), '],\n',
  '  "n_barangays": ', nrow(read.csv(wk[1])), ',\n',
  '  "n_weeks": ', nrow(rows), '\n}'
)
writeLines(json, file.path(DAT, "dengue.json"))
cat("wrote dengue.json\n")

cat("\nall figures and data written\n")
