################################################################################
#  ONE-STOP ACCUMULATOR SIMULATOR (snapshot-based reactivity)
#  - Plots update ONLY when you click "Run simulation"
#  - Parameter values are snapshotted on click and stored with the results
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
#  Helpers: RNGs & simulators
rtruncnorm <- function(n, mean, sd, a = 0) {
  x <- rnorm(n, mean, sd)
  while (any(bad <- x <= a))
    x[bad] <- rnorm(sum(bad), mean, sd)
  x
}

simulate_rw <- function(n, p, h, a, zprop, dt){
  # n - ntrials, p = prob step, h = step size, a = boundary, zprop = bias, dt = drift time step
  RT <- numeric(n); correct <- logical(n); start <- numeric(n); steps <- integer(n)
  z <- a * zprop
  max_steps <- max(1L, ceiling(20 / dt))  # cap ~20s like others
  for (i in seq_len(n)) {
    x <- z; t <- 0L
    for (k in seq_len(max_steps)) {
      t <- t + 1L
      x <- x + if (runif(1) < p) h else -h
      if (x >= a) { RT[i] <- t*dt; correct[i] <- TRUE;  start[i] <- z; steps[i] <- t; break }
      if (x <= 0) { RT[i] <- t*dt; correct[i] <- FALSE; start[i] <- z; steps[i] <- t; break }
    }
    if (steps[i] == 0L) { RT[i] <- NA_real_; correct[i] <- NA; start[i] <- z; steps[i] <- NA_integer_ }
  }
  data.frame(RT=RT, correct=correct, start=start, steps=steps)
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
#  UI
ui <- fluidPage(
  titlePanel("Accumulator-Model Simulator"),
  sidebarLayout(
    sidebarPanel(
      selectInput("model", "Choose model:",
                  c("Random Walk" = "rw",
                    "Linear Ballistic Accumulator" = "lba",
                    "Drift Diffusion Model"        = "ddm",
                    "Circular Diffusion Model"     = "cdm")),
      hr(),
      # common
      numericInput("n", "N-trials", 1000, 1, 10000, 1),
      uiOutput("paramUI"),
      actionButton("run", "Run simulation", class = "btn-primary"),
      downloadButton("savecsv", "Save CSV", class = "btn-success")
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
  
  # dynamic parameter panel (reacts immediately, but this does NOT trigger plots)
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
           ),
           rw = tagList(
             sliderInput("p",  "Step up prob  p", 0, 1, 0.55, 0.01),
             sliderInput("h",  "Step size     h",  0.01, 0.2, 0.05, 0.01),
             sliderInput("a",  "Threshold     a",  0.1, 3, 1, 0.05),
             sliderInput("z",  "Start prop   z/a", 0, 1, 0.5, 0.05),
             sliderInput("dt", "Step time     dt", 0.001, 0.02, 0.005, 0.001)
           )
    )
  })
  
  # ── central state store: NULL until user runs ----------------------
  simState <- reactiveVal(NULL)
  
  # Clear plots by clearing state whenever model changes
  observeEvent(input$model, {
    simState(NULL)
  }, ignoreInit = TRUE)
  
  # Run button: snapshot params and compute once ----------------------
  observeEvent(input$run, {
    mdl <- isolate(input$model)  # snapshot which model
    n   <- isolate(input$n)
    
    if (mdl == "lba") {
      p <- list(v = isolate(input$v),
                s = isolate(input$s),
                b = isolate(input$b),
                a = isolate(input$a),
                n = n)
      df <- simulate_lba(p$n, p$v, p$s, p$b, p$a)
      
    } else if (mdl == "ddm") {
      p <- list(v = isolate(input$v),
                s = isolate(input$s),
                a = isolate(input$a),
                z = isolate(input$z),
                n = n)
      df <- simulate_ddm(p$n, p$v, p$s, p$a, p$z)
      
    } else if (mdl == "rw") {
      p <- list(p = isolate(input$p),
                h = isolate(input$h),
                a = isolate(input$a),
                z = isolate(input$z),
                dt = isolate(input$dt),
                n = n)
      df <- simulate_rw(p$n, p$p, p$h, p$a, p$z, p$dt)
      
    } else { # cdm
      p <- list(v     = isolate(input$v),
                s     = isolate(input$s),
                R     = isolate(input$R),
                theta = isolate(input$theta),
                n     = n)
      dat <- simulate_cdm(p$n, p$v, p$s, p$R, p$theta)
      names(dat) <- c("RT","angle")
      target <- p$theta*pi/180
      err <- ((dat$angle - target + pi) %% (2*pi)) - pi
      dat$correct <- abs(err) <= pi/6   # pseudo "correct" to reuse helpers
      df <- dat
    }
    
    simState(list(model = mdl, params = p, df = df))
  })
  
  # ── Histogram (per model) -----------------------------------------
  output$histPlot <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    
    df <- st$df
    if (st$model == "cdm") {
      par(mfrow = c(1,1), mar = c(4,4,3,1))
      hist(df$RT, breaks = 40, col = "grey85", border = "white",
           main = "RT distribution", xlab = "RT (s)")
      abline(v = mean(df$RT, na.rm = TRUE), lwd = 3)
      legend("topright", legend = "mean RT", lwd = 3, bty = "n")
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
  
  # ── Route the three secondary panels based on model ----------------
  output$extraPlot1 <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (st$model %in% c("lba","ddm","rw")) {
      pathPlot(st$df, st$model, st$params)
    } else {
      polarErrPlot(st$df, st$params)
    }
  })
  
  output$extraPlot2 <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (st$model %in% c("lba","ddm","rw")) {
      densPlot(st$df)
    } else {
      rtPolarPlot(st$df, st$params)
    }
  })
  
  output$extraPlot3 <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (st$model %in% c("lba","ddm","rw")) {
      plot.new(); title("")  # intentionally blank for LBA/DDM
    } else {
      cdmPathPlot(st$df, st$params) # CDM Brownian paths
    }
  })
  
  # ── Save CSV of current simulation ----------------------------------------
  output$savecsv <- downloadHandler(
    filename = function(){
      st <- simState()
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      if (is.null(st)) return(paste0("accumulator_no-data_", ts, ".csv"))
      mdl <- st$model
      p   <- st$params
      num <- function(x) format(round(x, 3), trim = TRUE, scientific = FALSE)
      key <- switch(mdl,
                    lba = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_b", num(p$b), "_a", num(p$a)),
                    ddm = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_a", num(p$a), "_z", num(p$z)),
                    cdm = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_R", num(p$R), "_theta", num(p$theta)),
                    paste0("n", p$n)
      )
      paste0("accumulator_", mdl, "_", key, "_", ts, ".csv")
    },
    content = function(file){
      st <- simState()
      if (is.null(st)) {
        write.csv(data.frame(note = "No data yet. Click 'Run simulation' first."), file, row.names = FALSE)
      } else {
        write.csv(st$df, file, row.names = FALSE)
      }
    }
  )
  
  # ---------- helper plotting functions ------------------------------
  # LBA: ballistic paths; DDM: Brownian paths conditioned on RT & bound
  pathPlot <- function(df, model, params) {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    xmax <- max(df$RT, na.rm = TRUE)
    
    if (model == "ddm") {
      # ----- DDM: Brownian trajectories -----
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, params$a),
             xlab = "Time", ylab = "Evidence",
             main = paste(ifelse(flag, "Correct", "Error"), "paths"))
        abline(h = params$a, lty = 3)
        
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow == 0L) { text(xmax/2, params$a/2, "No trials", col = "grey50"); next }
        cols <- adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.18)
        
        for (j in idx[seq_len(nshow)]) {
          p <- ddm_path(
            rt       = df$RT[j],
            a        = params$a,
            z        = df$start[j],
            to_upper = df$correct[j],
            v        = params$v,
            s        = params$s,
            dt       = 0.002
          )
          if (!is.null(p)) lines(p$t, p$x, col = cols)
        }
      }
      
    } else if (model == "rw") {
      # ----- Random Walk: discrete steps paths reconstructed from RT & bound -----
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, params$a),
             xlab = "Time", ylab = "Evidence",
             main = paste(ifelse(flag, "Correct", "Error"), "paths (RW)"))
        abline(h = params$a, lty = 3)
        
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow == 0L) { text(xmax/2, params$a/2, "No trials", col = "grey50"); next }
        cols <- grDevices::adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.18)
        
        for (j in idx[seq_len(nshow)]) {
          # fallback if df$steps is missing/NA
          stps <- df$steps[j]
          if (is.null(stps) || is.na(stps)) {
            stps <- max(1L, round(df$RT[j] / params$dt))
          }
          
          p <- rw_path(
            rt       = df$RT[j],
            a        = params$a,
            z        = df$start[j],
            to_upper = isTRUE(df$correct[j]),
            h        = params$h,
            dt       = params$dt,
            steps    = stps
          )
          if (!is.null(p)) lines(p$t, p$x, col = cols)
        }
      }
    } else { # LBA
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, params$b),
             xlab = "Time", ylab = "Accumulated Evidence",
             main = paste(ifelse(flag, "Correct", "Error"), "paths"))
        abline(h = params$b, lty = 3)
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow > 0L) {
          for (j in idx[seq_len(nshow)]) {
            lines(c(0, df$RT[j]), c(df$start[j], params$b),
                  col = adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.15))
          }
        } else {
          text(xmax/2, params$b/2, "No trials", col = "grey50")
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
  
  polarErrPlot <- function(df, params){
    target <- params$theta*pi/180
    ang  <- circular(df$angle, units="radians")
    dens <- density.circular(ang, bw = 20)
    maxR <- max(1,dens$y)*1.2
    plot(dens, shrink=1, xlim=c(-maxR,maxR), ylim=c(-maxR,maxR),
         main="Angular-error density", col="darkorange", lwd=2)
    arrows(0,0, cos(target), sin(target), length=.08, lwd=2)
  }
  
  rtPolarPlot <- function(df, params){
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
    tar <- params$theta * pi/180
    arrows(0, 0, maxR*cos(tar), maxR*sin(tar), length = 0.06, lwd = 2)
  }
  
  cdmPathPlot <- function(df, params){
    get_num <- function(x, default) {
      y <- suppressWarnings(as.numeric(x))
      if (length(y) != 1 || is.na(y) || !is.finite(y)) default else y
    }
    R     <- get_num(params[["R"]], 1)
    theta <- get_num(params[["theta"]], 0) * pi/180
    v     <- get_num(params[["v"]], 0.25)
    s     <- max(0, get_num(params[["s"]], 0.15))
    
    if (!is.finite(R) || R <= 0 || !is.finite(v) || !is.finite(s)) {
      plot.new(); title("Invalid CDM params. Click Run again.", col.main = "grey50"); return(invisible())
    }
    
    nshow  <- min(100, nrow(df))
    par(mar = c(1,1,3,1))
    plot(NA, xlim=c(-R*1.1, R*1.1), ylim=c(-R*1.1, R*1.1),
         asp=1, axes=FALSE, main=paste("First", nshow, "diffusion paths"))
    symbols(0,0, circles=R, inches=FALSE, add=TRUE, fg="grey60")
    ramp <- grDevices::colorRampPalette(c("steelblue","grey80"))
    dt <- 0.02
    for (i in seq_len(nshow)) {
      x <- y <- 0; xs <- ys <- numeric(0)
      iter <- 0L; iter_max <- 1e6L
      while ((x^2 + y^2) < (R^2)) {
        if (iter >= iter_max) break
        iter <- iter + 1L
        xs <- c(xs, x); ys <- c(ys, y)
        x  <- x + v*dt*cos(theta) + rnorm(1, 0, s*sqrt(dt))
        y  <- y + v*dt*sin(theta) + rnorm(1, 0, s*sqrt(dt))
      }
      ang  <- atan2(y, x)
      delta<- ((ang - theta + pi) %% (2*pi)) - pi
      errA <- abs(delta)
      idx  <- min(101, max(1, floor(errA/pi*100)+1))
      col  <- grDevices::adjustcolor(ramp(101)[idx], 0.35)
      if (length(xs) > 1) lines(xs, ys, col = col)
    }
    points(0,0,pch=3)
  }
  
  
  rw_path <- function(rt, a, z, to_upper, h, dt, steps = NULL) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n  <- if (is.null(steps)) max(1L, round(rt / dt)) else as.integer(steps)
    t  <- (0:n) * dt
    
    target <- if (isTRUE(to_upper)) a else 0
    net    <- target - z
    k      <- round(net / h)
    U      <- (n + k) / 2
    D      <- n - U
    
    if (U < 0 || D < 0 || abs(U - round(U)) > 1e-8) {
      x <- z + (net / (n*dt)) * t
      x <- pmin(pmax(x, 0), a)
      return(data.frame(t = t, x = x))
    }
    
    U <- as.integer(round(U)); D <- as.integer(round(D))
    stepsigns <- sample(c(rep(1L, U), rep(-1L, D)))
    x <- numeric(n + 1); x[1] <- z
    for (i in 1:n) x[i + 1] <- x[i] + h * stepsigns[i]
    x <- pmin(pmax(x, 0), a)
    data.frame(t = t, x = x)
  }
}

shinyApp(ui, server)
