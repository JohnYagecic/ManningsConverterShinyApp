library(shiny)
library(ggplot2)

ui<-fluidPage(
  sidebarLayout(
    sidebarPanel(
      
  tags$h3("Manning's Input Variables"),
  sliderInput("Depth", "Expected Max. Depth (m)", min=0.1, max=20, value=10),
  sliderInput("n", "Manning's n", min=0.000, max=0.300, value=0.04),
  sliderInput("WidthB", "Bottom Width (m)", min=1, max=1000, value=50),
  sliderInput("SideS", "Side Slope (m/m)", min=0.01, max=1, value=0.5),
  sliderInput("SlopeB", "Bed Slope (m/m)", min=0.0001, max=0.5, value=0.003),
  
  img(src="Mannings.jpg"),
  img(src="MultExp.jpg")
    ),
  
  mainPanel(
    tags$h2("Convert Manning's Equation Input to Multiplier Exponent Expressions"),
  
  #tableOutput("table"),
  "An app for fitting muliplier and exponent values for V=aQ^b and Depth=cQ^d using Mannings equation as input.",
  
  plotOutput("scatter"),
  
  "Programmed by John Yagecic, P.E.  (JYagecic@gmail.com)",
  tags$br(),
  tags$br(),
  tags$a(href="https://en.wikipedia.org/wiki/Manning_formula", "More about Manning's Equation."),
  tags$br(),
  tags$a(href="https://www.epa.gov/exposure-assessment-models/water-quality-analysis-simulation-program-wasp", "More about the WASP water quality model."),
  tags$br(),
  tags$br(),
  "WASP and other flow and water quality models ask users to input multiplier and exponent values relating velocity and depth to discharge.",
  "These relationships take the form of V=aQ^b and depth=cQ^d, where the values a, b, c, and d describe the curve that approximates paired points from other sources.",
  "When paired field measurements are lacking, Manning's equation provides an estimate of open channel flow based on channel characteristics.",
  tags$br(),
  "This app takes Manning's equation input variables and fits nonlinear least squares parameters to estimate a, b, c, and d over 20 equal increments of depth ",
  "from 0.05 meters to the user specified maximum.",
  "The computations assume a symmetrical trapezoidal channel with continuous sloped sides.",
  tags$br(),
  tags$br(),
  "If you use this product or the underlying code in any professional or academic product, please consider ",
  "using a citation such as:",
  tags$br(),
  tags$br(),
  "Yagecic, John, August 2016.  A web app fitting multiplier and exponent parameters relating velocity and depth to discharge using Manning's equation as input.",
  "https://johnyagecic.shinyapps.io/ManningsConverterApp/",
  tags$br(),
  tags$br(),
  tags$a(href="http://www.nj.gov/drbc/", "Get the script"),
  tags$br(),
  tags$br(),
  "There is also a batch script for experienced R users for converting multiple reaches from a .csv file.  That script is available",
  tags$a(href="https://github.com/JohnYagecic/ManningsConverterBatch", "at this link.")
  )
  )
)

server<-function(input, output){
  
  
  Depth <- reactive({seq(0.05, input$Depth[1], by=(input$Depth/20))})
  SideRun<-reactive({Depth()/input$SideS})
  AreaSide <- reactive({0.5*SideRun()*Depth()})
  Hyp<-reactive({sqrt(Depth()^2 + SideRun()^2)})
  AreaTotal <- reactive({(Depth() * input$WidthB) + 2*AreaSide()})
  WettedPerim <- reactive({input$WidthB + 2*Hyp()})
  RadiusHyd <- reactive({AreaTotal() / WettedPerim()})
  Velm <- reactive({(1/input$n)*((RadiusHyd())^(2/3))*(input$SlopeB^0.5)})
  Q<-reactive({Velm()*AreaTotal()})
  
  ManningDF <- reactive({
    data.frame(Depth=Depth(),Area=AreaTotal(),WettedPerimeter=WettedPerim(),
               HydraulicRadius=RadiusHyd(),Velocity=Velm(),Discharge=Q())
  })
  
  output$scatter<-renderPlot({
    
    model.1<-lm(log(Velm()) ~ Q())
    start <- list(a=exp(coef(model.1)[1]), b=coef(model.1)[2]) # this is the part the allows nls function to work
    nlmod.VgivenQ<-nls(Velm() ~ a * Q() ^ b, start=start)
    
    mya<-coef(nlmod.VgivenQ)[1]
    myb<-coef(nlmod.VgivenQ)[2]
    
    model.2<-lm(log(Depth()) ~ Q())
    start2 <- list(c=exp(coef(model.2)[1]), d=coef(model.2)[2]) # this is the part the allows nls function to work
    nlmod.DgivenQ<-nls(Depth() ~ c * Q() ^ d, start=start2)
    
    myc<-coef(nlmod.DgivenQ)[1]
    myd<-coef(nlmod.DgivenQ)[2]
    
    par(mfrow=c(1,2))
    plot(Q(), Velm(), xlab="Q (cubic meters per second)", ylab="Velocity (m/s)")
    points(Q(), mya*Q()^myb, type="l", col="red")
    text(quantile(Q(),probs=0.75), quantile(Velm(), probs=0.35), paste("a =", round(mya,4)), cex=1.5, font=2)
    text(quantile(Q(),probs=0.75), quantile(Velm(), probs=0.25), paste("b =", round(myb,4)), cex=1.5, font=2)
    legend(x="topleft", legend=c("Mannings", "V=aQ^b"), col=c("black", "red"), pch=c(1,NA), lty=c(NA,1))
    
    
    plot(Q(), Depth(), xlab="Q (cubic meters per second)", ylab="Depth (m)")
    points(Q(), myc*Q()^myd, type="l", col="red")
    text(quantile(Q(),probs=0.75), quantile(Depth(), probs=0.35), paste("c =", round(myc,4)), cex=1.5, font=2)
    text(quantile(Q(),probs=0.75), quantile(Depth(), probs=0.25), paste("d =", round(myd,4)), cex=1.5, font=2)
    legend(x="topleft", legend=c("Mannings", "depth=cQ^d"), col=c("black", "red"), pch=c(1,NA), lty=c(NA,1))
    
    
  })
  
}

shinyApp(ui=ui, server=server)