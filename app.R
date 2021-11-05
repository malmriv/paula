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
ui <- fluidPage(title = "PAULA",
  titlePanel(img(src = "https://github.com/malmriv/paula/blob/main/www/logo.png?raw=true",width="31%"),
             tags$head(tags$link(rel="shortcut icon", href="https://github.com/malmriv/paula/blob/main/www/favicon.png?raw=true"))),
  sidebarLayout(
    sidebarPanel(
      #Aquí va algo de texto con instrucciones y tal.
      fileInput("upload", "Subir un archivo.",multiple=F,
                buttonLabel="Buscar..."),accept = c(".tsv", ".txt"),
      actionButton("compute", "Analizar"),
      DT::dataTableOutput("rawdata"),
    ),
    mainPanel(
      column(10,p("¡Hola! Bienvenid@ a PAULA (Protein Analysis Using a Least-squares Approach).
        Aquí puedes obtener una estimación razonablemente buena de la estructura secundaria de tu proteína
        problema a partir de su espectro de dicroísmo circular en el ultravioleta lejano. Para que el análisis
        funcione correctamente es necesario que tu archivo cumpla ciertos requisitos.
        En primer lugar, debe ser un archivo .txt con dos columnas separadas por un tabulador o por un espacio.
        Los decimales deben estar indicados con puntos, no comas. La primera columna debe contener longitudes de
        onda expresadas en nanómetros. La segunda columna debe contener la cantidad de giro de la luz registrada por el espectrómetro DC (las unidades son irrelevantes).
        También es necesario que ",strong("la longitud de onda esté acotada en el rango [170-250] nm. "),
        a(href="https://raw.githubusercontent.com/malmriv/paula/main/sample_proteins/cytochrome_c.txt","Esto")," es un ejemplo
        de un archivo válido (se permite la notación científica). Una vez subas un 
        archivo y hagas click en Analizar, obtendrás una imagen con tu resultado. Si te parece que la estimación
        no es suficientemente buena puedes reintentar el análisis tantas veces como consideres. ",
        a(href="https://github.com/malmriv/paula/blob/main/README.md","En este enlace"),"pueden encontrarse
        algunos detalles sobre el funcionamiento de esta aplicación. Si el servicio deja de funcionar o hay cualquier
        duda acerca de su funcionamiento se me puede contactar mediante mi correo académico, ",a(href="mailto:malmriv@correo.ugr.es","malmriv [@] correo.ugr.es"),
        br(),div("Manuel Almagro Rivas, 2021. Licencia Creative Commons Zero v1.0 Universal.",style="font-size:80%; color: #808080;"))),
      plotOutput(
        "results",
        width = "90%",
        height = "500px",
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
}

# shinyApp()
shinyApp(ui = ui, server = server)