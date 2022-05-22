#156Final
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(stats4)
library(ggfortify)
stylesheet <- tags$head(tags$style(HTML('
    .main-header .logo {
      font-family: "Georgia", Times, "Times New Roman", serif;
      font-weight: bold;
      font-size: 24px;
    }
    
    .large-emphasis {
      font-weight: bold;
      text-decoration: underline;
      font-family: "Georgia", Times, "Times New Roman", serif;
    }
    
    .small-emphasis {
      font-weight: bold;
    }
    
    .bact {
      color: red;
      font-style: italic;
    }
    
    .dis {
      color: blue;
    }
    
    .ADHD {
      color: red;
    }
    
    .Bul {
      color: green;
    }
    
    .Epi {
      color: blue;
    }
    
    .Dep {
      color: purple;
    }
  ')
))
#Import Data
agdata <- read.csv("modAGdata.csv")[-1]
#Fix Age
bad <- which(agdata$Age > 100)
agdata <- agdata[-bad, ]
#Rearrange Columns
agdata <- agdata[,c("Sex", "Age", "Weight", "BMI", "ADHD", "Bulimia", "Epilepsy", "Depression", "Bacteroidaceae", "Ruminococcaceae", "Intrasporangiaceae", "Enterobacteriaceae", "Erysipelotrichaceae", "Bacillaceae", "Kineosporiaceae")]
#Identify Factors
facs <- c(1,5,6,7,8)
#Modify Sex Data
fs <- which(agdata$Sex)
agdata$Sex[fs] <- "F"
agdata$Sex[-fs] <- "M"
#Identify Bacteria Counts and Create Log Counts Data
bact <- c(9,10,11,12,13,14,15)
lagdata <- agdata
lagdata[bact] <- log(lagdata[bact]+1)
bactnames <- c("Bacter.", "Rumin.", "Intra.", "Enter.", "Erysi.", "Bacil.", "Kineo.")
#Bayes Variables
abBay <- data.frame(replicate(4, rep(1, 2)))
dis <- c(5, 6, 7, 8)
colnames(abBay) <- c("ADHD", "Bulimia", "Epilepsy", "Depression")
ntot <- length(agdata[,1])
pDis <- c(sum(agdata$ADHD), sum(agdata$Bulimia), sum(agdata$Epilepsy), sum(agdata$Depression))
pDis <- pDis/ntot
bmin <- 0
bmax <- 1
#t variables
toohigh <- 0
toolow <- 0
justright <- 0

#The user interface
header <- dashboardHeader(title = "American Gut Project Data",
                          titleWidth = 500)
sidebar <- dashboardSidebar(width = 150,
                            radioButtons("dset","Bacterial Count Type", choices = c(
                                "Standard Counts" = "std",
                                "Log Counts" = "log"
                            ))
)
body <- dashboardBody(
    fluidRow(stylesheet,
             column(width=2,
                    h3("American Gut Data", class = "large-emphasis"),
                    h4("Summary", class = "small-emphasis"),
                    h5("15 Columns, 13620 Rows"),
                    h5("5 Factors, 10 Numeric"),
                    h5("7 Bacteria, 4 Disorders"),
                    h4("Factor Columns", class = "small-emphasis"),
                    tags$ul(tags$li("Sex"),
                            tags$li("ADHD", class = "dis"),
                            tags$li("Bulimia", class = "dis"),
                            tags$li("Epilepsy", class = "dis"),
                            tags$li("Depression", class = "dis")),
                    h4("Numeric Columns", class = "small-emphasis"),
                    tags$ul(tags$li("Age"),
                            tags$li("Weight"),
                            tags$li("BMI"),
                            tags$li("Bacter.", class = "bact"), 
                            tags$li("Rumin.", class = "bact"), 
                            tags$li("Intra.", class = "bact"), 
                            tags$li("Enter.", class = "bact"), 
                            tags$li("Erysi.", class = "bact"), 
                            tags$li("Bacil.", class = "bact"), 
                            tags$li("Kineo.", class = "bact"))
             ),
             column(width = 5,
                    h3("Principal Component Analysis (PCA)", class = "large-emphasis"),
                    plotOutput("PCA"),
                    br(),
                    selectInput("PCAcol", "Coloring", colnames(agdata)),
                    br(),
                    actionBttn("PCAgo", "Update Plot", color = "success"),
                    br(),
                    h3("Student t Confidence Intervals", class = "large-emphasis"),
                    plotOutput("Confidence"),
                    br(),
                    selectInput("Connum", "Numeric", colnames(agdata)[-facs]),
                    br(),
                    sliderInput("Consamp", "Sample Size", 5, 100, 10),
                    br(),
                    sliderInput("Conconf", "Confidence", 75, 95, 90),
                    br(),
                    actionBttn("Tgo", "Update Plot", color = "success"),
                    br(),
                    h4("Interval Testing:", class = "small-emphasis"),
                    uiOutput("TooHigh"),
                    uiOutput("TooLow"),
                    uiOutput("JustRight")
             ),
             column(width = 4,
                    h3("Regression", class = "large-emphasis"),
                    plotOutput("Regress"),
                    br(),
                    selectInput("Regresstype", "Regression Type", c("Linear", "Logistic")),
                    br(),
                    selectInput("Regressfac", "Factor", colnames(agdata)[dis]),
                    br(),
                    selectInput("Regressnum", "Numeric", colnames(agdata)[-facs]),
                    br(),
                    actionBttn("Reggo", "Update Plot", color = "success"),
                    br(),
                    h3("Bayesian Methods", class = "large-emphasis"),
                    plotOutput("Bayes", brush = "plot_brush"),
                    br(),
                    actionBttn("BayinitR", "Initialize Plot (Random)"),
                    br(),
                    actionBttn("BayinitU", "Initialize Plot (Uniform)"),
                    br(),
                    actionBttn("Baygo", "Update Plot", color = "success"),
                    br(),
                    h4("Full Data Probabilities:", class = "small-emphasis"),
                    h5("ADHD Probability: ", pDis[1], class = "ADHD"),
                    h5("Bulimia Probability: ", pDis[2], class = "Bul"),
                    h5("Epilepsy Probability: ", pDis[3], class = "Epi"),
                    h5("Depression Probability: ", pDis[4], class = "Dep")
            )
    )
)
ui <- dashboardPage(header, sidebar, body, skin = "black") #other colors available

#Functions that implement the mathematics
#This file must go into the same directory as app.R
#source(".R")

#Additional functions are OK here, but no variables

updatePCA <- function (f1) {
    renderPlot({
        X <- dta[bact]
        X.pca <- prcomp(X, center = TRUE, scale. = FALSE)
        autoplot(X.pca, data = dta, colour = f1)
    })
}

updateReg <- function (f, n, l) {
    renderPlot({
        fac <- dta[,f]
        num <- dta[,n] 
        tbl <- table(num, fac)
        x <- as.numeric(rownames(tbl))
        prob <- numeric(length(x))
        for (i in 1:length(x)) {
            prob[i] = sum(fac[which(num == x[i])])/sum(num == x[i])
        }
        plot(x, prob)
        if(l == "Linear") {
            lm <- lm(prob~x)
            abline(lm, col = "red")
        } else {
            try({
                    MLL <- function(alpha, beta) -sum(log(exp(alpha+beta*num)/(1+exp(alpha+beta*num)))*fac + log(1/(1+exp(alpha+beta*num)))*(1-fac))
                    results <- mle(MLL, start = list(alpha = 1, beta = 1))
                    curve(exp(results@coef[1]+results@coef[2]*x)/ (1+exp(results@coef[1]+results@coef[2]*x)),col = "blue", add=TRUE)
                }, silent = TRUE)
            
        }
    })
}

updateBayes <- function () {
    renderPlot({
        curve(dbeta(x, abBay[1,1], abBay[2,1]), col = "Red", xlab = "Disorder Probability", ylab = "Likelihood", xlim = c(bmin, bmax))
        curve(dbeta(x, abBay[1,2], abBay[2,2]), col = "Green", add = TRUE)
        curve(dbeta(x, abBay[1,3], abBay[2,3]), col = "Blue", add = TRUE)
        curve(dbeta(x, abBay[1,4], abBay[2,4]), col = "Purple", add = TRUE)
        abline(v = pDis[1], col = "Red", lty = 2)
        abline(v = pDis[2], col = "Green", lty = 2)
        abline(v = pDis[3], col = "Blue", lty = 2)
        abline(v = pDis[4], col = "Purple", lty = 2)
    })
}

updateConfidencePlot <- function (nmr, smps, conf) {
    renderPlot({
        nmr <<- dta[, nmr]
        toohigh <<- 0
        toolow <<- 0
        justright <<- 0
        mu <<- mean(nmr)
        alph <<- 1 - conf/100
        plot(x = c(0, max(nmr)), y = c(1,100), type = "n", xlab = "", ylab = "")
        for (i in 1:1000) {
            x <<- sample(nmr, smps)
            L <<- mean(x) + qt(alph/2, smps - 1) * sd(x)/sqrt(smps)
            U <<- mean(x) + qt(1 - alph/2, smps - 1) * sd(x)/sqrt(smps)
            if (L > mu) toohigh <<- toohigh + 1
            else if (U < mu) toolow <<- toolow + 1
            else justright <<- justright + 1
            if(i <= 100) segments(L, i, U, i)
        }
        abline (v = mu, col = "red")
    })
}

resetValues <- function () {
    bmin <<- 0
    bmax <<- 1
    abBay[,] <<- 1
    toohigh <<- 0
    toolow <<- 0
    justright <<- 0
}

server <- function(session, input, output) {
    observeEvent(input$dset,{
        dta <<- switch(input$dset,
                       std = agdata,
                       log = lagdata
        )
        resetValues()
        output$Bayes <<- updateBayes()
        output$Confidence <<- updateConfidencePlot(isolate(input$Connum), isolate(input$Consamp), isolate(input$Conconf))
        output$TooHigh <<- renderUI(h5("Confidence Interval Too High: ", toohigh)) 
        output$TooLow <<- renderUI(h5("Confidence Interval Too Low: ", toolow)) 
        output$JustRight <<- renderUI(h5("Confidence Interval Correct: ", justright)) 
        output$PCA <<- updatePCA(isolate(input$PCAcol))
        output$Regress <<- updateReg(isolate(input$Regressfac), isolate(input$Regressnum), isolate(input$Regresstype))
    })
    observeEvent(input$Reggo, {
        output$Regress <<- updateReg(isolate(input$Regressfac), isolate(input$Regressnum), isolate(input$Regresstype))
    })
    observeEvent(input$PCAgo, {
        output$PCA <<- updatePCA(isolate(input$PCAcol))
    })
    observeEvent(input$BayinitR, {
        resetValues()
        b <<- runif(1, 0, 100)
        abBay[2,] <<- b
        abBay[1,] <<- runif(4, 0, b)
        output$Bayes <<- updateBayes()
    })
    observeEvent(input$BayinitU, {
        resetValues()
        output$Bayes <<- updateBayes()
    })
    observeEvent(input$Baygo, {
        abBay[2,] <<- abBay[2,] + 100
        succ <<- numeric(4)
        for (i in 1:4){
            succ[i] <<- sum(sample(dta[,dis[i]], 100))
        }
        abBay[1,] <<- abBay[1,] + succ
        output$Bayes <<- updateBayes()
    })
    observeEvent(input$plot_brush,{
        bmin <<- input$plot_brush$xmin
        bmax <<- input$plot_brush$xmax
        output$Bayes <<- updateBayes()
        session$resetBrush("plot_brush")
    })
    
    observeEvent(input$Tgo,{
        output$Confidence <<- updateConfidencePlot(isolate(input$Connum), isolate(input$Consamp), isolate(input$Conconf))
        output$TooHigh <<- renderUI(h5("Confidence Interval Too High: ", toohigh)) 
        output$TooLow <<- renderUI(h5("Confidence Interval Too High: ", toolow)) 
        output$JustRight <<- renderUI(h5("Confidence Interval Correct: ", justright)) 
    })
}

#Run the app
shinyApp(ui = ui, server = server)

#https://abimm117.shinyapps.io/156FinalProject/
