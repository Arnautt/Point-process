library(shiny)
library(ggplot2)


hawkes_simul <- function(lambda, Tmax, noyau) {
  
  
  Ti <- 0
  T_list <- vector(mode = 'numeric')
  T_list[1] <- 0
  eps <- 10^(-10)
  
  
  
  L <- function(x, list) {           
    if (length(list)==0) {lambda}    
    else {
      to_sum <- list[list < x]
      lambda + sum(noyau(x - to_sum))
    }
  }
  
  
  i = 2                        # G?n?re le processus de Hawkes jusqu'? Tmax
  while (Ti < Tmax) {          # Par thinning
    M <- L(Ti + eps, T_list)
    E <- rexp(1,M)
    Ti <- Ti + E
    U <- runif(1,0,M)
    if ((Ti < Tmax) & (U <= L(Ti, T_list))) {
      T_list[i] <- Ti
      i = i + 1
    }
  }
  
  T_list <- T_list[-1]
  T_list
  
}


find_inter <- function(x, vec) {
  N <- length(vec)
  idx = ((1:N)[vec > x])[1] - 1
  if (is.na(idx)) {idx = 0}
  idx
}



intgr_nlk <- function(Mat1, Mat2) {
  newMat1 <- Mat1
  newMat1[newMat1[,1] < 0,1] <- 0
  newMat1[newMat1[,1] > Tmax,1] <- Tmax
  
  newMat2 <- Mat2
  newMat2[newMat2[,1] < 0,1] <- 0
  newMat2[newMat2[,1] > Tmax,1] <- Tmax
  
  res = 0
  
  tps1 <- newMat1[,1]
  tps2 <- newMat2[,1]
  if (length(tps1) == length(tps2) && sum(tps1 - tps2) == 0) {T_tilde <- tps1}   # lorsque Mat1=Mat2, cas particulier.
  else {T_tilde <- sort(unique(c(tps1, tps2)))}                                          # sinon, on trie l'ensemble des temps o? il y a des sauts
  
  N <- length(T_tilde)                                                           # pour chaque intervalle, on regarde la valeur des 2 fct cpm
  trans_1 <- 0                                                                   # on fait le produit, et on intègre
  trans_2 <- 0
  
  
  for (i in 1:(N-1)) {
    t1 <- T_tilde[i]     
    t2 <- T_tilde[i+1]
    idx1 <- find_inter(t1, tps1)
    idx2 <- find_inter(t1, tps2)
    if (isTRUE(idx1 == 0) || isTRUE(idx2 == 0)) {}
    else {
      trans_1 <- newMat1[idx1, 2]
      trans_2 <- newMat2[idx2, 2]
      res = res + (t2 - t1)*trans_1*trans_2
    }
  }
  
  res
}





nlk <- function(l, k) {       #delta
  Tps <- Ts_list[[paste("T",l, sep = "")]]   
  N_tps <- length(Tps)
  col1 <- c()
  col2 <- c()
  col1 <- unique(c(Tps + (k-1)*delta, Tps + k*delta))
  col1 <- sort(col1)          # col1 = ensemble de tous les temps où "?a saute" de +/-1
  N_res <- length(col1)
  Res <- matrix(nrow = N_res, ncol = 2)                        # résultat stocké dans une matrice à deux colonnes et 2*N lignes
  Res[,1] <- col1
  col2 <- rep(0, N_res)
  
  for (i in 1:N_tps) { 
    for (j in 1:(N_res-1)) {
      if (col1[j] >= Tps[i]+(k-1)*delta && col1[j] < Tps[i]+k*delta) {col2[j] = col2[j] + 1}
    }
  }
  
  col2[N_res] <- 0
  Res[,2] <- col2
  Res
}




intgr_cpm <- function(Mat) {
  newMat <- Mat
  newMat[newMat[,1] < 0,1] <- 0
  newMat[newMat[,1] > Tmax,1] <- Tmax
  res = 0
  tps <- newMat[,1]
  pts <- newMat[,2]
  
  
  for (i in 1:(length(pts)-1)) {
    res = res + (tps[i+1]-tps[i])*pts[i]
  }
  
  
  res
}




nmlk <- function(m, l, k) {
  res = 0
  tps_l <- Ts_list[[paste("T",l, sep = "")]]
  tps_m <- Ts_list[[paste("T",m, sep = "")]]
  tps_m <- tps_m[tps_m > 0 && tps_m <= Tmax]
  N <- length(tps_m)
  
  for (t in 1:N) {
    to_sum <- tps_l >= (tps_m[t] - k*delta) & tps_l < (tps_m[t] - (k-1)*delta)
    res = res + sum(to_sum)
  }
  
  res
}


b_m <- function(m) {
  res <- rep(0, 1+K*M)
  tps <- Ts_list[[paste("T",m, sep = "")]]
  res[1] <- length(tps)
  
  for (l in 1:M) {
    for (k in 1:K) {
      res[1+ (l-1)*K + k] <- 1/sqrt(delta) * nmlk(m, l, k)
    }
  }
  
  res
}




nlk_t <- function(l, k, x) {
  tps <- Ts_list[[paste("T",l, sep = "")]]
  to_sum <- tps >= (x - (k*delta)) & tps < (x - (k-1)*delta)
  sum(to_sum)
}


reconstruct_fct <- function(a_hat_m) {
  
  f <<- function(x) {
    res = a_hat_m[1]
    for (l in 1:M) {
      for (k in 1:K) {
        res = res + a_hat_m[1 + (l-1)*K + k]*(1/sqrt(delta))*nlk_t(l, k, x)
      }
    }
    res
  }
  
  return(f)
  
}


np_estimation <- function(K, delta, Tmax) {
  
  A <- delta*K
  G <- matrix(ncol = 1+K, nrow = 1+K)
  G[1,1] <- Tmax
  
  
  for (k in 1:K) {
    g <- (1/sqrt(delta))*intgr_cpm(nlk(1,k))
    if (length(g)==0) {G[1, 1+k] = G[1+k, 1] = 0}
    else {
      G[1, 1+k] = g
      G[1+k, 1] = g
    }
  }
  
  for (k in 1:K) {
    for (k_prime in 1:K) {
      trans1 <- nlk(1,k)
      trans2 <- nlk(1, k_prime)
      g <- intgr_nlk(trans1, trans2)
      G[1+k, 1+k_prime] <- (1/sqrt(delta))*g
    }
  }
  
  
  G_inv <- solve(G)
  b1 <- b_m(1)
  G_inv %*% b1
  
  
}






ui <- fluidPage(
  
  
  titlePanel("Processus de Hawkes univarié"),
  
  fluidRow(column(width = 4, textInput("phi", "Choisir la fonction d'intensité phi : ", "1.3*exp(-1.4*x)"), numericInput("lambda", "Paramètre lambda :", min = 1, max = 5, value = 1, step = 0.1)),
           column(width = 4, sliderInput("Tmax", "Jusqu'au temps :", min = 5, max = 100, value = 20), helpText("Remarque : l'estimation non paramétrique est coûteuse en terme de temps de calcul, cela peut prendre un petit moment.")),
           column(width = 4, sliderInput("K", "Paramètre K :", min = 5, max = 20, value = 10), numericInput("delta", "Paramètre delta :", min = 0.2, max = 2, value = 1, step = 0.1), actionButton("go", "Lancer le processus"))
  ),
  
  
  fluidRow(column(width = 12, plotOutput("comparaison")))
  
  # comparaison : vraie, vraisemblance, np, lasso ?
  
)




server <- function(input, output) {
  
  
  
  observeEvent(input$go, {
    
    input$go
    
    ## Génération des temps ##
    
    eval(parse(text = paste('phi <<- function(x)', input$phi, sep="")))
    #phi <<- Vectorize(phi)
    Tps <<- hawkes_simul(input$lambda, input$Tmax, phi)
    true_intensity <<- function(x, list) {
      if (length(list)==0) {input$lambda}
      else {
        to_sum <- list[list < x]
        return(input$lambda + sum(phi(x - to_sum)))
      }
    }
    
    delta <<- input$delta
    Tmax <<- input$Tmax
    M <<- 1
    K <<- input$K
    lambda <<- input$lambda
    
    ## Estimation non paramétrique ##
    
    Lambda_list <<- list()
    Ts_list <<- list()
    Lambda_list[[paste("L",1, sep = "")]] <<- true_intensity
    Ts_list[[paste("T",1, sep = "")]] <<- Tps
    
    #a_hat <- list()
    #a_hat[[paste("a",1, sep = "")]] <- np_estimation(input$K, input$delta, input$Tmax)
    a_hat_m <- np_estimation(input$K, input$delta, input$Tmax)
    np_func_est <- reconstruct_fct(a_hat_m)
    
    
    ## Estimation par vraisemblance ##
    
    N <- length(Tps)
    
    log_likelihood <- function(theta) {
      A_list <- rep(0,N)
      A_list[1] <- 0
      for (i in 2:N) {
        A_list[i] <- (1 + A_list[i-1]) * exp(-theta[3]*(Tps[i]-Tps[i-1]))
      }
      sum(log(theta[1] + (theta[2] * A_list))) - theta[1]*Tps[N] + (theta[2]/theta[3])*sum(exp(-theta[3]*(Tps[N] - Tps)) - 1)
    }
    
    mlog_likelihood <- function(theta) {
      -log_likelihood(theta)
    }
    
    est_par <- optim(c(2,2,2), mlog_likelihood)$par
    vrais_func_est <- function(x, list) {          
      if (length(list)==0) {est_par[1]}    
      else {
        to_sum <- list[list < x]
        est_par[1] + sum(est_par[2]*exp(-est_par[3]*(x - to_sum)))
      }
    }
    
    
    ## Représentation ##
    
    
    inter = seq(0, Tmax, by=0.1)
    to_plot1 = rep(0, length(inter))
    to_plot2 = rep(0, length(inter))
    to_plot3 = rep(0, length(inter))
    for (i in 1:length(inter)) {
      x <- inter[i]
      to_plot1[i] <- true_intensity(x, Tps)
      to_plot2[i] <- vrais_func_est(x, Tps)
      to_plot3[i] <- np_func_est(x)
    }
    
    
    
    mydf = data.frame(interv = inter, true_int = to_plot1, vraisemblance = to_plot2, np_estimation = to_plot3)
    
    p <- ggplot(mydf, aes(interv)) +
      geom_line(aes(interv, true_int, color = "red")) +
      geom_line(aes(interv, vraisemblance, color = "green")) +
      geom_line(aes(interv, np_estimation, color = "blue")) +
      xlab("") + 
      ylab("") +
      ggtitle("Comparaison des estimations") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_colour_manual(values = c("blue", "green", "red"), name = "Intensités",
                          labels = c("Non paramétrique", "Vraisemblance", "Réelle"),guide = guide_legend(reverse = TRUE))
    
    
    
    
    output$comparaison <- renderPlot({
      print(p)
    })
    
    
  })
  
  
  
  
}






shinyApp(ui = ui, server = server)

