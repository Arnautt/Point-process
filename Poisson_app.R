library(shiny)
library(ggplot2)


# Fonctions utiles à l'application


Sliding_window <- function(x) {                        # Prend en arguments x, liste de réels
  t <- rep(0, length(x))                               # Renvoie un vecteur de 0.5 si -1<=x=<+1 sinon 0
  test = x <= 1 & x>=-1
  t[test == 1] <- (1/2)
  t
}

Gaussien <- function(x) {
  dnorm(x)
}

Epanechnikov <- function(x) {
  t <- rep(0, length(x))
  test = x <= 1 & x>=-1
  xt <- x[test == 1]
  t[test == 1] <- (3/4)*(1-xt**2)
  t
}

Triangulaire <- function(x) {
  t <- rep(0, length(x))                               
  test = x <= 1 & x>=-1
  t[test == 1] <- 1 - abs(t[test == 1])
  t
}


gen_one_sample <- function(Tmax, lambd) {        # Génère un processus de poisson homogène
  Tn = vector(mode = 'numeric')                  # Renvoie une LISTE contenant les temps, le param, T max, estimation de lambda  
  Tn[1] <- 0                                     
  i = 1
  while (Tn[i] < Tmax) {                   
    Tn[i+1] <- Tn[i] + rexp(1, lambd)   
    i <- i+1
  }
  N_T <- length(Tn)
  if (Tn[N_T] > Tmax) {                    
    Tn <- Tn[1:length(Tn)-1]
    N_T <- N_T - 1
  }
  N_T <- N_T - 1                                                     # car T[1] = 0 ne doit pas compter dans nos ?venements
  list(Tn = Tn, lambd = lambd, Tmax = Tmax, lamb_hat = N_T/Tmax)     # Tn contient le temps 0 comme premier ?lement
}


gen_by_thin <- function(f, I) {                                  # Génère une trajectoire d'intensité f sur I=[0, Tmax]
  Tmax <- I[2]                                                   # Renvoie un vecteur contenant les temps Tn (sans le temps 0)
  majorant <- optimize(f, interval = I, maximum = T)$objective
  Tn_barre <- gen_one_sample(Tmax, majorant)$Tn
  Tn_barre <- Tn_barre[2:length(Tn_barre)]
  N <- length(Tn_barre)
  Un <- runif(N)
  test <- Un <= f(Tn_barre)/majorant
  To_n <- (1:N)[test == 1]
  Tn <- Tn_barre[To_n]
  Tn
}  


nsamp_thin <- function(n, f, Tmax) {                            # Génère n processus de Poisson inhomogène d'intensité f
  N_agreg <- vector(mode = 'numeric')                           # Renvoie l'ensemble des temps de chaque simulation N_agreg
  for (i in 1:n) {
    N_agreg <- c(N_agreg, gen_by_thin(f, I = c(0, Tmax)))
  }
  N_agreg
}


KDE_fun <- function(Data_aggreg, n, h, K, x_list, Tmax, f) {     # Permet de représenter l'estimation par noyau 
  Kh_x <- vector(mode = 'numeric', length = length(x_list))
  Kh <- function(t) {(1/h)*K(t/h)}                                          
  for (i in 1:length(x_list)) {                         
    Kh_x[i] <- (1/n) * sum(Kh(x_list[i] - Data_aggreg))
  }
  df <- data.frame(xs = x_list, ys = Kh_x)
  j <- ggplot(df, aes(xs, ys)) +
    geom_line() +
    stat_function(fun = f, color = "blue") +
    xlim(0, Tmax) +
    xlab("") +
    ylab("") +
    ggtitle("Fonction d'intensité (bleu) et estimation par noyau") +
    theme(plot.title = element_text(hjust = 0.5))
    
  print(j)
}


ui <- fluidPage(
  
  titlePanel("Processus de Poisson non homogène"),

  fluidRow(column(width = 4, textInput("fun", "Choisir la fonction d'intensité : ", "-x^2 + 10*x")),
           column(width = 4, sliderInput("Tmax", "Jusqu'au temps Tmax :", min = 5, max = 50, value = 10)),
           column(width = 4, sliderInput("n", "Nombre d'observations :", min = 5, max = 100, value = 50))
           ),
  
      
  fluidRow(column(width = 4, selectInput("kernel", "Choisir le noyau :", c("Sliding_window", "Gaussien", "Epanechnikov","Triangulaire"), "Epanechnikov")),
           column(width = 4, numericInput("h", "Paramètre de lissage h :", min = -2, max = 2, value = 1, step = 0.05)),
           column(width = 4, actionButton("go", "Lancer le processus")),
           column(width = 3, helpText("Remarque : il faut que Tmax soit divisible par h pour éviter des erreurs dans la représentation"))
           ),
      
  
  fluidRow(column(width = 12, 
                  tabsetPanel(
                    tabPanel("Intensité", plotOutput("intensity")), 
                    tabPanel("Histogramme", plotOutput("histo")),
                    tabPanel("Noyau", plotOutput("kernel"))
                    )
                )
            )
                      
  
    
)




server <- function(input, output) {
  
  
  observeEvent(input$go, {
    
    input$go
    eval(parse(text = paste('f <<- function(x)', input$fun, sep="")))
    eval(parse(text = paste('kern <<- function(x)', input$kernel, "(x)", sep="")))
    N_agreg <- nsamp_thin(input$n, f, input$Tmax)
    N_1 <- gen_by_thin(f, c(0, input$Tmax))
    
    mydf = data.frame(xaxis = N_1, yaxis = rep(0, length(N_1)))
    
    p <- ggplot(mydf, aes(xaxis, yaxis)) +
      geom_point(shape=3) +
      stat_function(fun = f, color = "blue") +
      xlim(0, input$Tmax) +
      xlab("") + 
      ylab("") +
      ggtitle("Fonction d'intensité et temps de sauts pour un processus (n=1)") +
      theme(plot.title = element_text(hjust = 0.5))
    

    
    df2 = data.frame(tps = seq(0, input$Tmax-input$h, input$h), cnt = (hist(N_agreg, breaks = seq(0,input$Tmax,input$h), plot = FALSE)$counts)/input$n)
    
    j <- ggplot(df2, aes(tps, cnt)) +
      geom_step() +
      stat_function(fun = f, color = "blue") +
      ggtitle(paste("Estimation par histogramme avec h = ", input$h, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(0, input$Tmax) +
      xlab("") +
      ylab("")
    
    
    


    
    
    output$intensity <- renderPlot({
      print(p)
    })
    
    
    output$histo <- renderPlot ({
      print(j)
    })
    
    
    output$kernel <- renderPlot ({
      KDE_fun(N_agreg, input$n, input$h, K=kern, x_list=seq(0, input$Tmax, by=0.2), input$Tmax, f)
    })
  })
  
  
  
}




shinyApp(ui = ui, server = server)

