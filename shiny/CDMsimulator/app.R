# CIRCULAR DIFFUSION MODEL
library(shiny)
library(circular)
library(Rcpp)

# ── Rcpp: fast RT sampler toward arbitrary angle ──────────────────────
cppFunction('
  Rcpp::NumericMatrix sim_cdm_cpp(int    n,
                                  double v,
                                  double s,
                                  double R,
                                  double dt,
                                  double cosT,
                                  double sinT) {
    Rcpp::NumericMatrix out(n, 2);              // rt, angle
    double vx = v * cosT,  vy = v * sinT;       // drift components
    for (int i = 0; i < n; ++i) {
      double x = 0.0, y = 0.0, t = 0.0;
      while (true) {
        t += dt;
        x += vx * dt + R::rnorm(0.0, s * sqrt(dt));
        y += vy * dt + R::rnorm(0.0, s * sqrt(dt));
        if (x*x + y*y >= R*R) {                 // crossed circular bound
          out(i,0) = t;
          out(i,1) = atan2(y, x);               // global angle
          break;
        }
        if (t > 20.0) {                         // failsafe
          out(i,0) = NA_REAL; out(i,1) = NA_REAL; break;
        }
      }
    }
    return out;
  }
')

simulate_cdm <- function(n, v, s, R, theta_deg, dt = 0.002) {
  th <- theta_deg * pi/180
  as.data.frame(sim_cdm_cpp(n, v, s, R, dt, cos(th), sin(th)))
}

# ── UI ────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Circular Diffusion Model Simulator"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n", "N-trials", 100, 1, 10000, 1),
      sliderInput("v", "Drift  v",  0.01, 1, 0.5, 0.01),
      sliderInput("s", "Noise  s",  0.01, 1, 0.1, 0.01),
      sliderInput("R", "Threshold radius R", 0.5, 3, 1.5, 0.05),
      sliderInput("theta", "Target angle (°)", 0, 360, 0, 5),
      actionButton("go", "Simulate")
    ),
    mainPanel(
      plotOutput("rtHist",   height = 300),
      plotOutput("polarErr", height = 300),
      plotOutput("rtPolar",  height = 300),
      plotOutput("pathPlot", height = 400)
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  cdm <- eventReactive(input$go, {
    df <- simulate_cdm(input$n, input$v, input$s, input$R, input$theta)
    names(df) <- c("RT","angle")
    df <- subset(df, !is.na(RT))
    target <- input$theta * pi/180
    err    <- ((df$angle - target + pi) %% (2*pi)) - pi  # wrap to (-π,π]
    df$err <- err
    df$correct <- abs(err) <= pi/6
    df
  })
  
  # 1 ─ RT histogram ---------------------------------------------------
  output$rtHist <- renderPlot({
    df <- cdm()
    hist(df$RT, breaks = 40, col = "grey85", border = "white",
         main = "RT distribution", xlab = "RT (s)")
    if (any(df$correct))
      abline(v = mean(df$RT[df$correct]), col = "steelblue", lwd = 3)
    legend("topright", legend = "mean correct", lwd = 3,
           col = "steelblue", bty = "n")
  })
  
  # 2 ─ Polar density (centred & not clipped) ----------------------------
  output$polarErr <- renderPlot({
    df <- cdm()
    target <- input$theta * pi/180
    
    ang  <- circular(df$angle, units = "radians")
    dens <- density.circular(ang, bw = 20)
    
    maxR <- max(1, dens$y) * 1.2      # enlarge frame if density > 1
    
    plot(dens, shrink = 1,
         xlim = c(-maxR, maxR), ylim = c(-maxR, maxR),
         main = "Angular-error density",
         col = "darkorange", lwd = 2)
    
    arrows(0, 0, cos(target), sin(target), length = 0.08, lwd = 2)
    text(1.05*cos(target), 1.05*sin(target), "Target",
         adj = c(sign(cos(target))*0.4, sign(sin(target))*0.4))
  })
  
  # 3 ─ Random-walk paths (≤ 100) ---------------------------------------
  output$pathPlot <- renderPlot({
    df <- cdm(); nshow <- min(100, nrow(df))
    target <- input$theta * pi/180
    par(mar = c(1,1,3,1))
    
    plot(NA,
         xlim = c(-input$R*1.1, input$R*1.1),
         ylim = c(-input$R*1.1, input$R*1.1),
         asp  = 1, axes = FALSE,
         main = paste("First", nshow, "diffusion paths"))
    symbols(0, 0, circles = input$R, inches = FALSE, add = TRUE, fg = "grey60")
    
    ramp <- colorRampPalette(c("steelblue", "grey80"))
    dt   <- 0.02                            # coarser step for speed
    
    for (i in seq_len(nshow)) {
      x <- y <- 0
      xs <- ys <- numeric(0)
      
      # simulate one path
      while (x^2 + y^2 < input$R^2) {
        xs <- c(xs, x);  ys <- c(ys, y)
        x  <- x + input$v*dt*cos(target) + rnorm(1, 0, input$s*sqrt(dt))
        y  <- y + input$v*dt*sin(target) + rnorm(1, 0, input$s*sqrt(dt))
      }
      
      ang  <- atan2(y, x)
      delta<- ((ang - target + pi) %% (2*pi)) - pi   # wrap to (-π, π]
      errA <- abs(delta)                             # absolute error ∈ [0, π]
      col  <- adjustcolor(ramp(101)[floor(errA/pi*100)+1], 0.35)
      
      lines(xs, ys, col = col)
    }
    points(0, 0, pch = 3)
  })
  
  # 4 ─ Polar RT profile -------------------------------------------------
  output$rtPolar <- renderPlot({
    df <- cdm()
    if (nrow(df) < 5) { plot.new(); title("Not enough data"); return() }
    
    # Bin into 12 × 30° sectors
    ang    <- (df$angle + 2*pi) %% (2*pi)          # wrap to [0,2π)
    sector <- cut(ang, breaks = seq(0, 2*pi, length.out = 13),
                  include.lowest = TRUE, labels = FALSE)
    medRT  <- tapply(df$RT, sector, mean, na.rm = TRUE)
    
    theta  <- seq(0, 2*pi, length.out = 13)        # sector centres + wrap
    rvals  <- c(medRT, medRT[1])                   # close the loop
    
    maxR   <- max(rvals, na.rm = TRUE) * 1.15
    plot(0, 0, type = "n",
         xlim = c(-maxR, maxR), ylim = c(-maxR, maxR),
         asp  = 1, axes = FALSE, xlab = "", ylab = "",
         main = "Mean RT by angle (12 bins)")
    
    # radial grid --------------------------------------------------------
    rgrid <- pretty(c(0, maxR), 4)[-1]
    symbols(rep(0, length(rgrid)), rep(0, length(rgrid)),
            circles = rgrid, inches = FALSE, add = TRUE,
            fg = "grey90")
    text(rgrid, 0,  labels = round(rgrid, 2), col = "grey60", pos = 4, cex = 0.7)
    
    # polar line ---------------------------------------------------------
    xs <- rvals * cos(theta)
    ys <- rvals * sin(theta)
    lines(xs, ys, col = "purple", lwd = 2)
    points(xs, ys, pch = 16, col = "purple")
    
    # target arrow -------------------------------------------------------
    tar <- input$theta * pi/180
    arrows(0, 0, maxR*cos(tar), maxR*sin(tar),
           length = 0.06, lwd = 2)
    text(1.05*maxR*cos(tar), 1.05*maxR*sin(tar), "Target",
         adj = c(sign(cos(tar))*0.4, sign(sin(tar))*0.4))
  })
}

shinyApp(ui, server)
