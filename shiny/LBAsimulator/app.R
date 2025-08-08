# LINEAR BALLISTIC ACCUMULATOR MODEL
# app.R ───────────────────────────────────────────────────────────────
library(shiny)

rtruncnorm <- function(n, mean, sd, a=0) {
  x <- rnorm(n, mean, sd)
  while (any(bad <- x <= a))          # resample negatives
    x[bad] <- rnorm(sum(bad), mean, sd)
  x
}

# ── LBA simulation ───────────────────────────────────────────────────
simulate_lba <- function(n, v, s, b, a) {
  start_correct <- runif(n, 0, a)
  start_error   <- runif(n, 0, a)
  drift_correct <- rtruncnorm(n, v,     s)
  drift_error   <- rtruncnorm(n, 1-v,   s)
  
  time_correct  <- (b - start_correct) / drift_correct
  time_error    <- (b - start_error)   / drift_error
  correct_win   <- time_correct < time_error
  
  ## winning accumulator’s start & drift (for path plot)
  win_start <- ifelse(correct_win, start_correct, start_error)
  win_drift <- ifelse(correct_win, drift_correct, drift_error)
  win_rt    <- ifelse(correct_win, time_correct,  time_error)
  
  data.frame(
    RT      = win_rt,
    correct = correct_win,
    start   = win_start,
    drift   = win_drift
  )
}

# ── UI ────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Linear Ballistic Accumulator Simulator"),
  sidebarLayout(
    sidebarPanel(
      numericInput("ntrials", "N-Trials:", 1000, 1, 10000, 1),
      sliderInput("v", "Drift Rate (v):", 0.01, 1, 0.5, 0.01),
      sliderInput("s", "Drift SD (s):",   0.01, 1, 0.1, 0.01),
      sliderInput("b", "Threshold (b):",   1,   10,  3,  0.1),
      sliderInput("a", "Starting Point (a):", 0, 10, 1, 0.1),
      actionButton("simulate", "Run Simulation")
    ),
    mainPanel(
      plotOutput("histPlot", height = 350),
      plotOutput("pathPlot", height = 350), 
      plotOutput("densPlot", height = 350) 
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  sim_data <- eventReactive(input$simulate, {
    validate(need(input$a <= input$b,
                  "Starting point (a) must be ≤ threshold (b)."))
    simulate_lba(input$ntrials, input$v, input$s, input$b, input$a)
  })
  
  # Histograms ---------------------------------------------------------
  output$histPlot <- renderPlot({
    df <- sim_data()
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    # -------- Correct RTs --------
    if (any(df$correct)) {
      hist(df$RT[df$correct],
           breaks = 30, col = "steelblue", border = "white",
           main = "Correct RTs", xlab = "RT (s)")
    } else {
      plot.new(); title("Correct RTs\n(no trials)", col.main = "grey50")
    }
    
    # -------- Error RTs ----------
    if (any(!df$correct)) {
      hist(df$RT[!df$correct],
           breaks = 30, col = "salmon", border = "white",
           main = "Error RTs",   xlab = "RT (s)")
    } else {
      plot.new(); title("Error RTs\n(no trials)", col.main = "grey50")
    }
  })
  
  
  # ── Trajectory plots ---------------------------------------
  output$pathPlot <- renderPlot({
    df    <- sim_data()
    nshow <- min(250, nrow(df))
    df    <- df[seq_len(nshow), ]
    
    xmax <- max(df$RT)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    # ---------- Correct ----------
    plot(NA, xlim = c(0, xmax), ylim = c(0, input$b),
         xlab = "Time", ylab = "Accumulated Evidence",
         main = paste("Correct (n =", sum(df$correct), ")"))
    abline(h = input$b, lty = 3)
    
    if (any(df$correct)) {
      apply(df[df$correct, ], 1, function(tr) {
        lines(c(0, tr["RT"]),
              c(tr["start"], input$b),
              col = adjustcolor("steelblue", alpha.f = 0.15))
      })
    } else {
      text(xmax/2, input$b/2, "No correct\ntrials", col = "grey50")
    }
    
    # ---------- Error ----------
    plot(NA, xlim = c(0, xmax), ylim = c(0, input$b),
         xlab = "Time", ylab = "Accumulated Evidence",
         main = paste("Error (n =", sum(!df$correct), ")"))
    abline(h = input$b, lty = 3)
    
    if (any(!df$correct)) {
      apply(df[!df$correct, ], 1, function(tr) {
        lines(c(0, tr["RT"]),
              c(tr["start"], input$b),
              col = adjustcolor("salmon", alpha.f = 0.15))
      })
    } else {
      text(xmax/2, input$b/2, "No error\ntrials", col = "grey50")
    }
  })
  
  # Density curves ------------------------------------------------------
output$densPlot <- renderPlot({
  df <- sim_data()
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # ---- Correct density ----
  n_cor <- sum(df$correct)
  if (n_cor >= 2) {
    d <- density(df$RT[df$correct])
    plot(d, col = "steelblue", lwd = 2,
         main = paste("Correct RT Density (n =", n_cor, ")"),
         xlab = "RT (s)")
    polygon(d, col = adjustcolor("steelblue", .30), border = NA)
  } else {
    plot.new(); title("Correct RT Density\n(< 2 trials)", col.main = "grey50")
    if (n_cor == 1) rug(df$RT[df$correct], col = "steelblue", lwd = 2)
  }

  # ---- Error density ----
  n_err <- sum(!df$correct)
  if (n_err >= 2) {
    d <- density(df$RT[!df$correct])
    plot(d, col = "salmon", lwd = 2,
         main = paste("Error RT Density (n =", n_err, ")"),
         xlab = "RT (s)")
    polygon(d, col = adjustcolor("salmon", .30), border = NA)
  } else {
    plot.new(); title("Error RT Density\n(< 2 trials)", col.main = "grey50")
    if (n_err == 1) rug(df$RT[!df$correct], col = "salmon", lwd = 2)
  }
})

  


}

shinyApp(ui, server)
