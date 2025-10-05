library(shiny)
library(ggplot2)

# Import functions
source("functions.R")

ui <- fluidPage(
  titlePanel("AMISE Bandwidth Selection — Simulation Study"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Sample size n", min = 50, max = 5000, value = 3000, step = 50),
      sliderInput("alpha", "Beta parameter α", min = 0.1, max = 8, value = 2, step = 0.5),
      sliderInput("beta", "Beta parameter β", min = 0.1, max = 8, value = 2, step = 0.5),
      sliderInput("sigma2", "Error variance σ²", min = 0.1, max = 3, value = 1, step = 0.1),
      checkboxInput("useOptN", "Use optimal number of blocks (via Cp criterion)", TRUE),
      uiOutput("Nchoice"),
      hr(),
      sliderInput("n_min", "Minimum sample size (for h vs n curve)", min = 50, max = 2000, value = 200, step = 100),
      sliderInput("n_max", "Maximum sample size (for h vs n curve)", min = 200, max = 5000, value = 2000, step = 100),
      numericInput("n_step", "Sample size step", value = 200, min = 50, step = 50)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", tableOutput("summary")),
        tabPanel("h vs N", plotOutput("plot_hN")),
        tabPanel("h vs n", plotOutput("plot_hn")),
        tabPanel("Beta & h_AMISE", 
                 fluidRow(
                   column(12, plotOutput("plot_beta", height = "400px")),
                 ),
                 fluidRow(
                   column(12, plotOutput("plot_hBeta", height = "400px"))  
                 )
        ),
        
        tabPanel("Data & Fit", plotOutput("plot_datafit"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive data
  dat <- reactive({
    simulate_data(input$n, input$alpha, input$beta, input$sigma2)
  })
  
  # N choice
  output$Nchoice <- renderUI({
    Nmax <- max(min(floor(input$n/20), 5), 1)
    if (input$useOptN) {
      helpText(paste("N_max =", Nmax))
    } else {
      sliderInput("N", "Number of blocks (N)", min = 1, max = Nmax, value = Nmax, step = 1)
    }
  })
  
  # Summary: N_opt or chosen N
  output$summary <- renderTable({
    X <- dat()$X; Y <- dat()$Y
    if (input$useOptN) {
      opt <- optimal_N(X, Y)
      est <- blockwise_fit(X, Y, opt$N_opt)
      data.frame(N_used = opt$N_opt,
                 Cp_value = round(opt$Cp_opt,3),
                 sigma2_hat = round(est$sigma2,3),
                 theta22_hat = round(est$theta22,3),h_AMISE = round(h_amise(input$n, est$sigma2, est$theta22),3))
    } else {
      est <- blockwise_fit(X, Y, input$N)
      data.frame(N_used = input$N,
                 sigma2_hat = round(est$sigma2,3),
                 theta22_hat = round(est$theta22,3),
                 h_AMISE = round(h_amise(input$n, est$sigma2, est$theta22),3))
    }
  })
  
  # h vs N
  output$plot_hN <- renderPlot({
    X <- dat()$X; Y <- dat()$Y
    n <- length(Y)
    Nmax <- max(min(floor(n/20), 5), 1)
    h_vals <- sapply(1:Nmax, function(N) {
      est <- blockwise_fit(X, Y, N)
      h_amise(n, est$sigma2, est$theta22)
    })
    opt <- optimal_N(X, Y)
    
    ggplot(data.frame(N=1:Nmax, h=h_vals), aes(N, h)) +
      geom_line() + geom_point() +
      geom_vline(xintercept = opt$N_opt, linetype = 2, color="green") +
      annotate("text", x = opt$N_opt, y = max(h_vals, na.rm=TRUE),
               label = "N optimal", vjust = -0.5, hjust = 0, color="green") +
      labs(title="h_AMISE vs N", y="h_AMISE", x="Number of blocks (N)") +
      theme_minimal()
    
  })
  
  
  # h vs n
  output$plot_hn <- renderPlot({
    n_seq <- seq(input$n_min, input$n_max, by=input$n_step)
    h_vals <- sapply(n_seq, function(n){
      d <- simulate_data(n, input$alpha, input$beta, input$sigma2)
      if (input$useOptN) {
        opt <- optimal_N(d$X, d$Y)
        est <- blockwise_fit(d$X, d$Y, opt$N_opt)
      } else {
        N <-input$N
        est <- blockwise_fit(d$X, d$Y, N)
      }
      h_amise(n, est$sigma2, est$theta22)
    })
    ggplot(data.frame(n=n_seq, h=h_vals), aes(n,h)) +
      geom_line() + geom_point() +
      labs(title="h_AMISE vs sample size n", y="h_AMISE") +
      theme_minimal()
  })
  
  #Beta distribution
  output$plot_beta <- renderPlot({
    X <- rbeta(5000, input$alpha, input$beta)
    
    ggplot(data.frame(X), aes(X)) +
      geom_histogram(bins=10, fill="skyblue", color="black") +
      labs(title=paste0("Beta(", input$alpha, ", ", input$beta, ") distribution"),
           x="x", y="Density") +
      theme_minimal()
  })
  
  
  # h vs (alpha, beta)
  output$plot_hBeta <- renderPlot({
    alpha_seq <- seq(0.1,8,by=1)
    beta_seq  <- seq(0.1,8,by=1)
    grid <- expand.grid(alpha=alpha_seq, beta=beta_seq)
    
    h_vals <- apply(grid, 1, function(p){
      d <- simulate_data(input$n, p[1], p[2], input$sigma2)
      if (input$useOptN) {
        opt <- optimal_N(d$X,d$Y)
        est <- blockwise_fit(d$X,d$Y,opt$N_opt)
      } else {
        est <- blockwise_fit(d$X,d$Y,input$N)
      }
      h_amise(input$n, est$sigma2, est$theta22)
    })
    
    grid$h <- h_vals
    
    ggplot(grid, aes(alpha, beta, fill=h)) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(title="h_AMISE across (alpha, beta)", 
           x=expression(alpha), y=expression(beta), fill="h") +
      theme_minimal()
  })
  
  
  # Data & Fit
  output$plot_datafit <- renderPlot({
    d <- dat()
    if (input$useOptN) {
      opt <- optimal_N(d$X,d$Y)
      N <- opt$N_opt
    } else {
      N <- input$N
    }
    fit <- blockwise_fit(d$X,d$Y,N)
    
    # Grid for polynomial fits
    xs <- seq(0,1,length.out=200)
    block_id <- cut(xs, breaks=seq(0,1,length.out=N+1), include.lowest=TRUE, labels=FALSE)
    yhat <- numeric(length(xs))
    for (j in 1:N) {
      idx <- which(block_id==j)
      df <- data.frame(x=d$X[d$X>=min(xs[idx]) & d$X<=max(xs[idx])], y=d$Y[d$X>=min(xs[idx]) & d$X<=max(xs[idx])])
      if(nrow(df)>=6){
        f <- lm(y ~ x + I(x^2)+I(x^3)+I(x^4), data=df)
        yhat[idx] <- predict(f, newdata=data.frame(x=xs[idx]))
      } else {
        yhat[idx] <- NA
      }
    }
    ggplot(d, aes(X,Y)) +
      geom_point(alpha=0.4) +
      stat_function(fun=m_true, colour="blue", size=1) +
      geom_line(data=data.frame(x=xs,y=yhat), aes(x,y), colour="red") +
      labs(title=paste("Data with true m(x) and blockwise fit (N=",N,")"), y="Y") +
      theme_minimal()
  })
}

shinyApp(ui, server)
