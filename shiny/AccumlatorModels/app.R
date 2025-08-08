################################################################################
#  ONE-STOP ACCUMULATOR SIMULATOR (LBA | DDM | CIRCULAR DDM)
################################################################################
library(shiny)
library(rtdists)     # for DDM
library(circular)    # for circular DDM
library(Rcpp)        # for fast CDM RTs

# ────────────────────────────────────────────────────────────────────────────
#  CDM core in Rcpp
cppFunction('
  Rcpp::NumericMatrix sim_cdm_cpp(int    n,
                                  double v,
                                  double s,
                                  double R,
                                  double dt,
                                  double cosT,
                                  double sinT) {
    Rcpp::NumericMatrix out(n, 2);              // rt, angle
    double vx = v * cosT,  vy = v * sinT;
    for (int i = 0; i < n; ++i) {
      double x = 0.0, y = 0.0, t = 0.0;
      while (true) {
        t += dt;
        x += vx * dt + R::rnorm(0.0, s * sqrt(dt));
        y += vy * dt + R::rnorm(0.0, s * sqrt(dt));
        if (x*x + y*y >= R*R) {
          out(i,0) = t;
          out(i,1) = atan2(y, x);
          break;
        }
        if (t > 20.0) { out(i,0)=NA_REAL; out(i,1)=NA_REAL; break; }
      }
    }
    return out;
  }
')

# ────────────────────────────────────────────────────────────────────────────
#  Helper functions for each model
rtruncnorm <- function(n, mean, sd, a = 0) {
  x <- rnorm(n, mean, sd)
  while (any(bad <- x <= a))
    x[bad] <- rnorm(sum(bad), mean, sd)
  x
}

simulate_lba <- function(n, v, s, b, a) {
  start_c <- runif(n, 0, a)
  start_e <- runif(n, 0, a)
  drift_c <- rtruncnorm(n, v,   s)
  drift_e <- rtruncnorm(n, 1-v, s)
  t_c <- (b - start_c) / drift_c
  t_e <- (b - start_e) / drift_e
  winC <- t_c < t_e
  data.frame(RT = ifelse(winC, t_c, t_e), correct = winC,
             start = ifelse(winC, start_c, start_e))
}

simulate_ddm <- function(n, v, s, a, zprop) {
  raw <- rdiffusion(n, a = a, v = v, s = s, z = a*zprop, t0 = 0)
  data.frame(RT = raw$rt, correct = raw$response == "upper",
             start = a * zprop)
}

simulate_cdm <- function(n, v, s, R, theta_deg) {
  th <- theta_deg*pi/180
  as.data.frame(sim_cdm_cpp(n, v, s, R, 0.002, cos(th), sin(th)))
}

# ────────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Accumulator-Model Simulator"),
  sidebarLayout(
    sidebarPanel(
      selectInput("model", "Choose model:",
                  c("Linear Ballistic Accumulator" = "lba",
                    "Drift Diffusion Model"        = "ddm",
                    "Circular Diffusion Model"     = "cdm")),
      hr(),
      # common
      numericInput("n", "N-trials", 1000, 1, 10000, 1),
      uiOutput("paramUI"),
      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      plotOutput("histPlot",   height = 300),
      plotOutput("extraPlot1", height = 300),
      plotOutput("extraPlot2", height = 300),
      plotOutput("extraPlot3", height = 400)
    )
  )
)

# ────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # dynamic parameter panel -------------------------------------------
  output$paramUI <- renderUI({
    switch(input$model,
           lba = tagList(
             sliderInput("v", "Drift rate  v", 0.01, 1, 0.5, 0.01),
             sliderInput("s", "Drift SD      s", 0.01, 1, 0.1, 0.01),
             sliderInput("b", "Threshold     b",   1, 10, 3, 0.1),
             sliderInput("a", "Start max     a",   0, 10, 1, 0.1)
           ),
           ddm = tagList(
             sliderInput("v", "Drift rate  v", -3, 3, 0.5, 0.05),
             sliderInput("s", "Noise σ      s", 0.1, 2, 1, 0.05),
             sliderInput("a", "Threshold   a", 0.1, 3, 1, 0.05),
             sliderInput("z", "Start prop  z/a", 0, 1, 0.5, 0.05)
           ),
           cdm = tagList(
             sliderInput("v", "Drift speed v", 0.01, 1, 0.25, 0.01),
             sliderInput("s", "Noise SD    s", 0.01, 1, 0.15, 0.01),
             sliderInput("R", "Radius      R", 0.5, 2, 1, 0.05),
             sliderInput("theta", "Target angle (°)", 0, 360, 0, 5)
           )
    )
  })
  
  # ── central data store: NULL until user runs ------------------------
  simData <- reactiveVal(NULL)
  
  # Clear plots by clearing data whenever model changes
  observeEvent(input$model, {
    simData(NULL)
  }, ignoreInit = TRUE)
  
  # Run button: compute and store data in simData ----------------------
  observeEvent(input$run, {
    df <- switch(input$model,
                 lba = simulate_lba(input$n, input$v, input$s, input$b, input$a),
                 ddm = simulate_ddm(input$n, input$v, input$s, input$a, input$z),
                 cdm = {
                   dat <- simulate_cdm(input$n, input$v, input$s, input$R, input$theta)
                   names(dat) <- c("RT","angle")
                   target <- input$theta*pi/180
                   err <- ((dat$angle - target + pi) %% (2*pi)) - pi
                   dat$correct <- abs(err) <= pi/6   # used only to reuse some helpers; not "errors"
                   dat
                 }
    )
    simData(df)
  })
  
  # ── Histogram (per model) -------------------------------------------
  output$histPlot <- renderPlot({
    df <- simData()
    if (is.null(df)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    
    if (input$model == "cdm") {
      par(mfrow = c(1,1), mar = c(4,4,3,1))
      hist(df$RT, breaks = 40, col = "grey85", border = "white",
           main = "RT distribution", xlab = "RT (s)")
      abline(v = mean(df$RT, na.rm = TRUE), col = "steelblue", lwd = 3)
      legend("topright", legend = "mean RT", lwd = 3, col = "steelblue", bty = "n")
    } else {
      par(mfrow = c(1,2), mar = c(4,4,3,1))
      # Correct
      if (any(df$correct, na.rm = TRUE)) {
        hist(df$RT[df$correct], breaks = 30, col = "steelblue", border = "white",
             main = "Correct RTs", xlab = "RT (s)")
      } else {
        plot.new(); title("Correct RTs\n(no trials)", col.main = "grey50")
      }
      # Error
      if (any(!df$correct, na.rm = TRUE)) {
        hist(df$RT[!df$correct], breaks = 30, col = "salmon", border = "white",
             main = "Error RTs", xlab = "RT (s)")
      } else {
        plot.new(); title("Error RTs\n(no trials)", col.main = "grey50")
      }
    }
  })
  
  # ── Route the three secondary panels based on model -----------------
  output$extraPlot1 <- renderPlot({
    df <- simData()
    if (is.null(df)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (input$model %in% c("lba","ddm")) {
      pathPlot(df, input$model)
    } else {
      polarErrPlot(df, input)
    }
  })
  
  output$extraPlot2 <- renderPlot({
    df <- simData()
    if (is.null(df)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (input$model %in% c("lba","ddm")) {
      densPlot(df)
    } else {
      rtPolarPlot(df, input)
    }
  })
  
  output$extraPlot3 <- renderPlot({
    df <- simData()
    if (is.null(df)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (input$model %in% c("lba","ddm")) {
      plot.new(); title("")  # intentionally blank for LBA/DDM
    } else {
      cdmPathPlot(df, input) # CDM Brownian paths
    }
  })
  
  # ---------- helper plotting functions ------------------------------
  # LBA: ballistic paths; DDM: Brownian paths conditioned on RT & bound
  pathPlot <- function(df, model) {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    xmax <- max(df$RT, na.rm = TRUE)
    
    if (model == "ddm") {
      # ----- DDM: Brownian trajectories -----
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, input$a),
             xlab = "Time", ylab = "Evidence",
             main = paste(ifelse(flag, "Correct", "Error"), "paths"))
        abline(h = input$a, lty = 3)
        
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow == 0L) { text(xmax/2, input$a/2, "No trials", col = "grey50"); next }
        cols <- adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.18)
        
        for (k in seq_len(nshow)) {
          j <- idx[k]
          p <- ddm_path(
            rt       = df$RT[j],
            a        = input$a,
            z        = df$start[j],
            to_upper = df$correct[j],
            v        = input$v,
            s        = input$s,
            dt       = 0.002
          )
          if (!is.null(p)) lines(p$t, p$x, col = cols)
        }
      }
      
    } else { # LBA
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, input$b),
             xlab = "Time", ylab = "Accumulated Evidence",
             main = paste(ifelse(flag, "Correct", "Error"), "paths"))
        abline(h = input$b, lty = 3)
        if (any(df$correct == flag)) {
          apply(df[df$correct == flag, ], 1, function(r)
            lines(c(0, r["RT"]), c(r["start"], input$b),
                  col = adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.15)))
        } else {
          text(xmax/2, input$b/2, "No trials", col = "grey50")
        }
      }
    }
  }
  
  # --- DDM Brownian path conditioned on RT & bound --------------------
  ddm_path <- function(rt, a, z, to_upper, v, s, dt = 0.002) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n  <- max(1L, ceiling(rt / dt))
    t  <- seq(0, rt, length.out = n + 1)
    
    # Drift + noise (Brownian motion)
    drift <- v * t
    W     <- c(0, cumsum(rnorm(n, mean = 0, sd = s * sqrt(dt))))
    x     <- z + drift + W
    
    # Brownian-bridge correction to force end at threshold
    target <- if (to_upper) a else 0
    corr   <- (target - x[length(x)]) / rt
    x      <- x + corr * t
    
    # Clamp tiny numeric overshoots
    x <- pmin(pmax(x, 0), a)
    
    data.frame(t = t, x = x)
  }
  
  densPlot <- function(df){
    par(mfrow = c(1,2), mar = c(4,4,3,1))
    for (flag in c(TRUE,FALSE)){
      n <- sum(df$correct==flag, na.rm = TRUE)
      if (n>=2){
        d <- density(df$RT[df$correct==flag])
        plot(d, col=ifelse(flag,"steelblue","salmon"), lwd=2,
             main=paste(ifelse(flag,"Correct","Error"),"RT Density"))
      } else {
        plot.new(); title("Not enough", col.main="grey50")
      }
    }
  }
  
  polarErrPlot <- function(df, inp){
    target <- inp$theta*pi/180
    ang  <- circular(df$angle, units="radians")
    dens <- density.circular(ang, bw = 20)
    maxR <- max(1,dens$y)*1.2
    plot(dens, shrink=1, xlim=c(-maxR,maxR), ylim=c(-maxR,maxR),
         main="Angular-error density", col="darkorange", lwd=2)
    arrows(0,0, cos(target), sin(target), length=.08, lwd=2)
  }
  
  rtPolarPlot <- function(df, inp){
    if (nrow(df) < 2) { plot.new(); title("Not enough data"); return(invisible()) }
    
    # Angles in [0, 2π)
    ang <- (df$angle + 2*pi) %% (2*pi)
    
    # 12 × 30° bins (keep empty bins so lengths stay consistent)
    sector <- factor(
      cut(ang, breaks = seq(0, 2*pi, length.out = 13),
          include.lowest = TRUE, labels = FALSE),
      levels = 1:12
    )
    mu <- tapply(df$RT, sector, mean, na.rm = TRUE)
    
    # Plot at BIN CENTERS (nicer than edges)
    edges   <- seq(0, 2*pi, length.out = 13)
    centers <- edges[-length(edges)] + diff(edges)[1]/2
    theta   <- c(centers, centers[1])   # wrap
    r       <- c(mu, mu[1])             # wrap
    
    if (all(is.na(r))) { plot.new(); title("No RTs in any sector", col.main="grey50"); return(invisible()) }
    
    maxR <- max(r, na.rm = TRUE) * 1.15
    plot(0, 0, type = "n",
         xlim = c(-maxR, maxR), ylim = c(-maxR, maxR),
         asp = 1, axes = FALSE, main = "Mean RT by angle (12 bins)")
    
    # radial grid
    rgrid <- pretty(c(0, maxR), 4); rgrid <- rgrid[rgrid > 0]
    if (length(rgrid)) {
      symbols(rep(0, length(rgrid)), rep(0, length(rgrid)),
              circles = rgrid, inches = FALSE, add = TRUE, fg = "grey90")
      text(rgrid, 0, labels = round(rgrid, 2), col = "grey60", pos = 4, cex = 0.7)
    }
    
    # polar line
    xs <- r * cos(theta);  ys <- r * sin(theta)
    lines(xs, ys, col = "purple", lwd = 2)
    points(xs, ys, pch = 16, col = "purple")
    
    # target arrow
    tar <- inp$theta * pi/180
    arrows(0, 0, maxR*cos(tar), maxR*sin(tar), length = 0.06, lwd = 2)
  }
  
  cdmPathPlot <- function(df, inp){
    nshow  <- min(100, nrow(df))
    target <- inp$theta * pi/180
    par(mar = c(1,1,3,1))
    plot(NA, xlim=c(-inp$R*1.1, inp$R*1.1), ylim=c(-inp$R*1.1, inp$R*1.1),
         asp=1, axes=FALSE, main=paste("First", nshow, "diffusion paths"))
    symbols(0,0, circles=inp$R, inches=FALSE, add=TRUE, fg="grey60")
    ramp <- colorRampPalette(c("steelblue","grey80"))
    dt <- 0.02
    for (i in seq_len(nshow)) {
      x <- y <- 0; xs <- ys <- numeric(0)
      while (x^2 + y^2 < inp$R^2) {
        xs <- c(xs, x); ys <- c(ys, y)
        x  <- x + inp$v*dt*cos(target) + rnorm(1, 0, inp$s*sqrt(dt))
        y  <- y + inp$v*dt*sin(target) + rnorm(1, 0, inp$s*sqrt(dt))
      }
      ang  <- atan2(y, x)
      delta<- ((ang - target + pi) %% (2*pi)) - pi
      errA <- abs(delta)                             # [0, π]
      col  <- adjustcolor(ramp(101)[floor(errA/pi*100)+1], 0.35)
      lines(xs, ys, col = col)
    }
    points(0,0,pch=3)
  }
  
}

shinyApp(ui, server)
