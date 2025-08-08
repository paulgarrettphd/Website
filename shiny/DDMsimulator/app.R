# DIFFUSION DECISION MODEL
library(shiny)
library(rtdists)    # exact Wiener RT simulator

# --------------------------------------------------------------------
#  DDM: draw RT + accuracy analytically via rwiener -------------------
simulate_ddm_rt <- function(n, v, s, a, zprop, t0 = 0) {
  rdiffusion(
    n     = n,
    a     = a,
    v     = v,        # drift rate  (positive → upper bound more likely)
    t0    = t0,
    z     = a * zprop,
    s     = s         # diffusion coefficient (noise)
  )
}

# --------------------------------------------------------------------
#  Generate a *sample path* that matches a given RT & choice ----------
#  Uses Brownian bridge conditioned to hit the correct bound.
ddm_path <- function(rt, a, z, correct, v, s, dt = 0.002) {
  if (is.na(rt)) return(NULL)
  n  <- ceiling(rt / dt)
  t  <- seq(0, rt, length.out = n + 1)
  
  # Unconstrained Brownian motion with drift v and noise s
  bm <- c(0, cumsum(rnorm(n, mean = v * dt, sd = s * sqrt(dt))))
  path <- z + bm
  
  # Brownian bridge correction to force hitting the bound at rt
  target <- if (correct) a else 0
  correction <- (target - path[length(path)]) / rt
  path <- path + correction * t
  
  data.frame(t = t, x = pmax(pmin(path, a), 0))  # trim minor jitter
}

# --------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Drift Diffusion Model Simulator"),
  sidebarLayout(
    sidebarPanel(
      numericInput("ntrials", "N-Trials:", 1000, 1, 10000, 1),
      sliderInput("v", "Drift rate  v:",  -3, 3,  0.5, 0.05),
      sliderInput("s", "Noise σ  (s):",   0.1, 2, 1,   0.05),
      sliderInput("a", "Threshold  a:",   0.1, 3, 1,   0.05),
      sliderInput("z", "Starting-point proportion  z/a:", 0, 1, 0.5, 0.05),
      actionButton("go", "Run Simulation")
    ),
    mainPanel(
      plotOutput("histPlot", height = 350),
      plotOutput("pathPlot", height = 350),
      plotOutput("densPlot", height = 350)
    )
  )
)

# --------------------------------------------------------------------
server <- function(input, output, session) {
  
  sim_data <- eventReactive(input$go, {
    raw <- simulate_ddm_rt(input$ntrials,
                           v = input$v, s = input$s,
                           a = input$a, zprop = input$z)
    data.frame(
      RT      = raw$rt,
      correct = raw$response == "upper",   # 'upper' vs 'lower'
      start   = input$a * input$z
    )
  })
  
  # ---- Histograms ---------------------------------------------------
  output$histPlot <- renderPlot({
    df <- sim_data(); par(mfrow = c(1, 2), mar = c(4,4,3,1))
    if (any(df$correct))
      hist(df$RT[df$correct], breaks = 30, col = "steelblue",
           border = "white", main = "Correct RTs", xlab = "RT (s)")
    else { plot.new(); title("Correct RTs\n(no trials)", col.main="grey50") }
    
    if (any(!df$correct))
      hist(df$RT[!df$correct], breaks = 30, col = "salmon",
           border = "white", main = "Error RTs", xlab = "RT (s)")
    else { plot.new(); title("Error RTs\n(no trials)", col.main="grey50") }
  })
  
  # ---- Trajectory panels (≤ 100 paths) ------------------------------
  output$pathPlot <- renderPlot({
    df <- sim_data()
    par(mfrow = c(1,2), mar = c(4,4,3,1))
    nshow <- min(100, nrow(df))
    df    <- df[seq_len(nshow), ]
    
    draw_set <- function(sub, col_line, title_txt) {
      xmax <- max(sub$RT)
      plot(NA, xlim = c(0, xmax), ylim = c(0, input$a),
           xlab = "Time", ylab = "Evidence", main = title_txt)
      abline(h = input$a, lty = 3)
      for (i in seq_len(nrow(sub))) {
        p <- ddm_path(sub$RT[i], input$a, sub$start[i], sub$correct[i],
                      v = input$v, s = input$s)
        if (!is.null(p))
          lines(p$t, p$x, col = adjustcolor(col_line, 0.15))
      }
    }
    
    draw_set(df[df$correct, ], "steelblue",
             paste("Correct paths (", sum(df$correct), ")", sep=""))
    draw_set(df[!df$correct, ], "salmon",
             paste("Error paths (", sum(!df$correct), ")", sep=""))
  })
  
  # ---- Density curves ----------------------------------------------
  output$densPlot <- renderPlot({
    df <- sim_data()
    par(mfrow = c(1, 2), mar = c(4,4,3,1))
    dens_panel <- function(idx, col_fill, label) {
      n <- sum(idx)
      if (n >= 2) {
        d <- density(df$RT[idx])
        plot(d, col = col_fill, lwd = 2,
             main = paste(label, "(n=", n, ")", sep=""),
             xlab = "RT (s)")
        polygon(d, col = adjustcolor(col_fill, 0.3), border = NA)
      } else {
        plot.new(); title(paste(label, "\n(< 2 trials)"), col.main="grey50")
        if (n == 1) rug(df$RT[idx], col = col_fill, lwd = 2)
      }
    }
    dens_panel(df$correct,     "steelblue", "Correct RT Density")
    dens_panel(!df$correct,    "salmon",    "Error RT Density")
  })
}

shinyApp(ui, server)
