library(shiny)
library(latex2exp)
library(DT)

#Aquí empieza el modelo en sí
#Leemos las proteínas modelo (Brahms, 1979)
bturn = read.csv("./data/beta-turn.csv",dec=",")
bsheet = read.csv("./data/beta-sheet.csv",dec=",")
ahelix = read.csv("./data/alpha-helix.csv",dec=",")
random = read.csv("./data/random-coil.csv",dec=",")

#Acotamos la lambda a la que trabajamos
lambda = seq(170,250,len=1000)

#Generamos un ajuste suave para cada dataset
bturnpol = loess(bturn$Curve1~bturn$x,span=0.2)
bsheetpol = loess(bsheet$Curve1~bsheet$x,span=0.2)
ahelixpol = loess(ahelix$Curve1~ahelix$x,span=0.2)
randompol = loess(random$Curve1~random$x,span=0.2)

#Now we define the function to minimise, a linear combination of all spectra
linear = function(wavelength,BS,BT,AH,RA) {
  #BS = % of beta sheet
  #BT = % of beta turn
  #AH = % of alpha helix
  #RA = % of random coil
  return(BS*predict(bsheetpol,wavelength) + BT*predict(bturnpol,wavelength) +
           AH*predict(ahelixpol,wavelength) + RA*predict(randompol,wavelength))
}

# ui object
ui <- fluidPage(
  titlePanel(p("Spatial app", style = "color:#3474A7")),
  sidebarLayout(
    sidebarPanel(
      #Aquí va algo de texto con instrucciones y tal.
      fileInput("upload", "Subir un archivo.",multiple=F,
                buttonLabel="buscar..."),accept = c(".tsv", ".txt"),
      DT::dataTableOutput("rawdata"),
      actionButton("compute", "Intentar ajuste")
    ),
    mainPanel(
      p("¡Hola! Bienvenid@ a PAULA, acrónimo de Protein Analysis Using a Least-squares Approach.
        Aquí puedes obtener una estimación razonablemente buena de la estructura secundaria de tu proteína
        problema a partir de su espectro de dicroísmo circular en el ultravioleta lejano. Para que el análisis
        funcione correctamente es necesario que tu archivo cumpla ciertos requisitos.
        En primer lugar, debe ser un archivo .txt con dos columnas separadas por una espacio o por un tabulador.
        Los decimales deben estar indicados con puntos, no comas. La primera columna debe contener longitudes de onda expresadas en nanómetros. La segunda columna debe
        contener la cantidad de giro de la luz registrada por el espectrómetro DC (las unidades son irrelevantes).
        También es necesario que la longitud de onda esté acotada en el rango [170-250] nm. Una vez subas un 
        archivo y hagas click en Intentar análisis, obtendrás una imagen con tu resultado. Si te parece que la estimación
        no es suficientemente buena puedes reintentar el análisis tantas veces como consideres. Si quieres
        saber cuál es el fundamento de PAULA o estás experimentando algún problema, puedes escribirme a
        malmriv [@] correo.ugr.es."),
      img(src = "http://shiny.rstudio-staging.com/tutorial/written-tutorial/lesson2/images/image-in-app.png"),
      plotOutput(
        "raw",
        width = "100%",
        height = "400px",
        click = NULL,
        dblclick = NULL,
        hover = NULL,
        brush = NULL,
        inline = FALSE
      ),
      plotOutput(
        "results",
        width = "100%",
        height = "400px",
        click = NULL,
        dblclick = NULL,
        hover = NULL,
        brush = NULL,
        inline = FALSE
      )
      
    )
  )
)

# server()
server <- function(input, output) {
  #Leer el archivo subido
  df_products_upload <- reactive({
    inFile <- input$upload
    if (is.null(inFile))
      return(NULL)
    df <- read.table(inFile$datapath, header = T)
    return(df)
  })
  #Mostrarlo como una tabla
  output$rawdata<- DT::renderDataTable({
    df <- df_products_upload()
    DT::datatable(df)
  })
  #Hacemos el ajuste
  observeEvent(input$compute, {
    df <- df_products_upload()
    wv = df[,1]
    cd = df[,2]
    seed = runif(4,min=0,max=10)
    results = nls(cd~linear(wv,BS,BT,AH,RA),data=list(df),
                  start=list(BS=seed[1],BT=seed[2],AH=seed[3],RA=seed[4]))
    #The coefficients need to be normalised:
    coefs = summary(results)$coefficients
    total = sum(abs(coefs[,1]))
    BS = abs(coefs[1,1]/total)
    BT = abs(coefs[2,1]/total)
    AH = abs(coefs[3,1]/total)
    RA = abs(coefs[4,1]/total)
    
    #Plot the result
    output$results = renderPlot({
      plot(df,type="l",col="red",lwd=2,xlab="wavelength (nm)",
           ylab=TeX("$\\theta_{MWE}$ (deg $cm^2$ $dmol^{-1}$)"),xlim=c(175,240),
           main="CD spectra: experimental data & best fit.")
      lines(wv,linear(wv,BS*total,BT*total,AH*total,RA*total),col="blue",lwd=2)
      mtext(paste("α-helix:",round(AH*100,0),"%, β-sheet:",round(BS*100,0),
                  "%, random coil:",round(RA*100,0),"%, β-turn: ",round(BT*100,0),"%"))
      legend("topright",c("Experimental","Best fit"),lty=c(1,1),lwd=c(2,2),
             col=c("red","blue"))
    })
    })
  
  #Graficamos los resultados fuera
  output$raw = renderPlot({
    df <- df_products_upload()
    #Hay que andarse con cuidado porque al tratar con objetos responsive
    #vamos a obtener un error cada vez que abramos la página sin un
    #archivo y acargado; usamos tryCatch para que, en caso de este "error",
    #el usuario no vea nada.
    tryCatch(
      expr = plot(df,main="Espectro suministrado (sin procesar).",
                  xlab="lambda (nm)",ylab=""),
      error = function(e) {}
    )
    
    
  })
}

# shinyApp()
shinyApp(ui = ui, server = server)