# Save as app.R and run in RStudio
library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Sensitivity Analysis of Compound Poisson Process S(t)"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Parameters"),
      sliderInput("lambda", "Inter-arrival Rate (lambda):", min=0.1, max=10, value=2, step=0.1),
      sliderInput("theta", "Jump Size Rate (theta):", min=0.1, max=5, value=1, step=0.1),
      
      h4("Time Settings"),
      selectInput("time_t", "Select Time (t) to View:", 
                  choices = c(10, 100, 1000, 10000), selected = 100),
      
      sliderInput("num_sims", "Number of Simulations:", min=1000, max=50000, value=10000),
      
      hr(),
      helpText("Note: As 't' increases, observe the convergence to Normal Distribution (CLT).")
    ),
    
    mainPanel(
      plotOutput("distPlot", height = "500px"),
      verbatimTextOutput("statsOutput")
    )
  )
)

server <- function(input, output) {
  
  sim_data <- reactive({
    t <- as.numeric(input$time_t)
    lambda <- input$lambda
    theta <- input$theta
    n_sims <- input$num_sims
    
    # 1. Simulate N(t) ~ Poisson(lambda * t)
    n_jumps <- rpois(n_sims, lambda * t)
    
    # 2. Simulate S(t)
    # If N=0, S=0. If N>0, S ~ Gamma(n, rate=theta)
    # We use vectorization where possible
    s_values <- numeric(n_sims)
    
    # Indices where jumps occurred
    nonzero_idx <- which(n_jumps > 0)
    
    if(length(nonzero_idx) > 0) {
      # Draw Gamma variables for all non-zero cases
      # rgamma(n, shape, rate)
      s_values[nonzero_idx] <- rgamma(length(nonzero_idx), 
                                      shape = n_jumps[nonzero_idx], 
                                      rate = theta)
    }
    
    return(s_values)
  })
  
  output$distPlot <- renderPlot({
    data <- data.frame(S_t = sim_data())
    t <- as.numeric(input$time_t)
    
    # Theoretical Mean and Variance for Normal Approx
    mu <- (input$lambda * t) / input$theta
    sigma2 <- (2 * input$lambda * t) / (input$theta^2)
    
    ggplot(data, aes(x=S_t)) +
      geom_histogram(aes(y=..density..), bins=50, fill="skyblue", color="black", alpha=0.7) +
      stat_function(fun = dnorm, args = list(mean = mu, sd = sqrt(sigma2)), 
                    color = "red", size = 1.5) +
      labs(title = paste("Distribution of S(t) at t =", t),
           subtitle = "Blue: Simulated | Red: Normal Approximation",
           x = "S(t)", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(size=18, face="bold"))
  })
  
  output$statsOutput <- renderText({
    s_val <- sim_data()
    t <- as.numeric(input$time_t)
    theo_mean <- (input$lambda * t) / input$theta
    theo_var <- (2 * input$lambda * t) / (input$theta^2)
    
    paste0(
      "--- Statistics ---\n",
      "Simulated Mean: ", round(mean(s_val), 4), "\n",
      "Theoretical Mean: ", round(theo_mean, 4), "\n",
      "Simulated Variance: ", round(var(s_val), 4), "\n",
      "Theoretical Variance: ", round(theo_var, 4)
    )
  })
}

shinyApp(ui = ui, server = server)
