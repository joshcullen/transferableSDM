
## Shiny web application to interactively explore maps of green turtle tracks ##

## Plots Leaflet map w/ inputs and reactivity
## Created 2023-07-13 by Josh Cullen (joshcullen10@gmail.com)

library(shiny)
library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(leaflet)
library(leafem)
library(cmocean)
library(MetBrewer)
library(htmltools)

source("leaflet_utils.R")  #helper functions


# Load turtle tracks

gom.tracks <- read_csv("GoM_Cm_Tracks_behav.csv") %>%
  arrange(date)
br.tracks <- read_csv("Brazil_Cm_Tracks_behav.csv") %>%
  arrange(date)
qa.tracks <- read_csv("Qatar_Cm_Tracks_behav.csv") %>%
  arrange(date)

tracks <- list(GoM = gom.tracks,
               Brazil = br.tracks,
               Qatar = qa.tracks)

# Load metadata for life stage by ID
dat.age <- read_csv("Turtle tag metadata.csv") %>%
  dplyr::select(Age, Ptt) %>%
  distinct()





# Create indexing column "month.year" and only relevant cols
tracks <- tracks %>%
  map(., ~{.x %>%
      mutate(month.year = as_date(date), .after = 'date') %>%
      mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"),
             id = as.character(id)) %>%
      dplyr::select(Region, id, date, month.year, x, y, behav) %>%
      st_changeCoords(data = ., coords = c('x','y'), crs_in = 3395, crs_out = 4326) %>%
      left_join(y = dat.age, by = c("id" = "Ptt"))
  })



# # Load environmental covariates
#
# cov_list_gom <- readRDS(gzcon(url("https://github.com/joshcullen/transferableSDM/raw/main/green_turtle_tracks/GoM_covars.rds"))) %>%
#   map(rast)
# cov_list_br <- readRDS(gzcon(url("https://github.com/joshcullen/transferableSDM/raw/main/green_turtle_tracks/Brazil_covars.rds"))) %>%
#   map(rast)
# cov_list_qa <- readRDS(gzcon(url("https://github.com/joshcullen/transferableSDM/raw/main/green_turtle_tracks/Qatar_covars.rds"))) %>%
#   map(rast)
#
#
# cov_list <- list(GoM = cov_list_gom,
#                  Brazil = cov_list_br,
#                  Qatar = cov_list_qa)
#
# # Modify units for depth and NPP
# for (i in seq_along(cov_list)) {
#   cov_list[[i]]$bathym <- cov_list[[i]]$bathym * -1  #convert depth to positive values
#   cov_list[[i]]$npp <- cov_list[[i]]$npp / 1000  #convert from mg to g C
# }










### UI ###

ui <- fluidPage(title = "Green Turtle Observations Compared to Depth, NPP, and SST",
                header = tags$head(includeHTML("gtag.html")),

                leafletOutput("leaflet_map", width = "100%", height = "750px"),

                tags$style("#controls {
                             background-color: #f0f0f0;
                             opacity: 0.5;
                             }
                            #controls:hover{
                             opacity: 1;
                             }"),

                absolutePanel(id = "controls",
                              class = "panel panel-default",
                              top = 200,
                              left = 25,
                              width = 250,
                              fixed = TRUE,
                              draggable = TRUE,
                              height = "auto",

                              h3("Choose region and data to map"),
                              selectInput(inputId = "region",
                                          label = "Select geographic region",
                                          choices = names(tracks),
                                          selected = names(tracks)[1]),
                              selectInput(inputId = "covar",
                                          label = "Select environmental variable",
                                          choices = c("bathym", "npp", "sst"),
                                          selected = "sst"),
                              selectInput(inputId = "month_year",
                                          label = "Select month-year of observations",
                                          choices = unique(tracks[[1]]$month.year),
                                          selected = unique(tracks[[1]]$month.year)[1])

                )  #close absolutePanel

)  #close fluidPage





### Server ###

server <- function(input, output, session) {

  # Define dynamic UI for month-year selection by Region
  region <- reactive({
    tracks[[input$region]]
  })

  observeEvent(region(), {
    freezeReactiveValue(input, "month_year")

    choices <- unique(region()$month.year)
    updateSelectInput(inputId = "month_year", choices = choices)
  })



  # Load environmental covariates as reactive object by region

  cov_list <- reactive({

    covar_rast <- switch(input$region,
           "GoM" = readRDS(gzcon(url("https://github.com/joshcullen/transferableSDM/raw/main/green_turtle_tracks/GoM_covars.rds"))),
           "Brazil" = readRDS(gzcon(url("https://github.com/joshcullen/transferableSDM/raw/main/green_turtle_tracks/Brazil_covars.rds"))),
           "Qatar" = readRDS(gzcon(url("https://github.com/joshcullen/transferableSDM/raw/main/green_turtle_tracks/Qatar_covars.rds")))
    ) %>%
      map(rast)


    # Keep only the selected covariate and month-year
    covar_rast <- covar_rast[[input$covar]]

    if (input$covar %in% c('npp','sst')) {
      covar_rast <- covar_rast[[input$month_year]]
    }


    # Modify values to adjust units
    if (input$covar == "bathym") {
      covar_rast <- covar_rast * -1  #convert depth to positive values
    }

    if (input$covar == "npp") {
      covar_rast <- covar_rast / 1000  #convert from mg to g C
    }


    covar_rast
  })





  # Generate color palette for tracks
  set.seed(123)
  turt.id <- tracks %>%
    map(pull, id) %>%
    map(unique) %>%
    unlist()

  col.pal <- sapply(names(MetBrewer::MetPalettes),
                    function(x) MetBrewer::met.brewer(palette_name = x, n = 5,
                                                      type = "continuous")) %>%
    as.vector() %>%
    sample(., size = length(turt.id))
  tracks.pal <- colorFactor(col.pal, factor(turt.id))



  # Generate Leaflet basemap
    output$leaflet_map <- renderLeaflet({
      leaflet() %>%
        addProviderTiles(provider = providers$Esri.WorldImagery, group = "World Imagery",
                         options = tileOptions(zIndex = -10)) %>%
        addProviderTiles(provider = providers$Esri.OceanBasemap, group = "Ocean Basemap",
                         options = tileOptions(zIndex = -10)) %>%
        addProviderTiles(provider = providers$OpenStreetMap, group = "Open Street Map",
                         options = tileOptions(zIndex = -10)) %>%
        addLayersControl(baseGroups = c("World Imagery", "Ocean Basemap", "Open Street Map"),
                         overlayGroups = c("Raster", "Tracks"),
                         options = layersControlOptions(collapsed = TRUE,
                                                        autoZIndex = FALSE)) %>%
        setView(lng = -87, lat = 25.5, zoom = 6) %>%
        addMouseCoordinates() %>%
        addScaleBar(position = "bottomleft") %>%
        addMeasure(position = "topleft",
                   primaryLengthUnit = "kilometers",
                   primaryAreaUnit = "hectares",
                   activeColor = "#3D535D",
                   completedColor = "#7D4479")
    })  #close renderLeaflet


    region.my <- reactive({
      tracks[[input$region]] %>%
        filter(month.year == input$month_year)
    }) %>%
      debounce(500)  #to delay reactive invalidation signal causing downstream issues w/ subsetting rasters after changing Region


    # Update Leaflet map w/ reactive layers
    observeEvent(c(input$region, input$covar, input$month_year), {

      req(region.my)


      proxy <- leafletProxy("leaflet_map") %>%
        clearMarkers() %>%
        clearImages() %>%
        fitBounds(lng1 = min(region.my()$x) - 1,
                  lng2 = max(region.my()$x) + 1,
                  lat1 = min(region.my()$y) - 1,
                  lat2 = max(region.my()$y) + 1) %>%
        addCircleMarkers(data = region.my(),
          lng = ~x,
          lat = ~y,
          radius = 5,
          group = "Tracks",
          stroke = FALSE,
          fillColor = ~tracks.pal(id),
          fillOpacity = 0.75,
          label = ~paste0("<b>ID:</b> ", id,
                          "<br> <b>Date:</b> ", date,
                          "<br> <b>Behavior:</b> ", behav,
                          "<br> <b>Life Stage:</b> ", Age) %>%
            lapply(htmltools::HTML),
          labelOptions = labelOptions(style = list(
                                        "font-size" = "14px"))
          )


      # Remove any existing legend and replace w/ updated one
      proxy %>%
        clearControls()


      # Add raster and associated legend by covar
      if (input$covar == "sst") {

        sst.pal <- colorNumeric(cmocean("thermal")(256),
                                domain = as.vector(values(cov_list())),
                                na.color = "transparent")

        proxy %>%
          addRasterImage(x = cov_list(),
                         colors = sst.pal,
                         opacity = 0.9,
                         group = "Raster") %>%
          addImageQuery(cov_list(),
                        group = "Raster",
                        project = TRUE) %>%
          addLegend_decreasing(pal = sst.pal,
                               values = as.vector(
                                 values(
                                   cov_list()
                                   )
                                 ),
                               title = "SST (\u00B0C)",
                               decreasing = TRUE) %>%
          addLegend(pal = tracks.pal,
                    values = region.my()$id,
                    title = "ID")

      } else if (input$covar == "bathym") {

        depth.pal <- colorNumeric(cmocean("deep")(256),
                                  domain = values(cov_list()),
                                  na.color = "transparent")

        proxy %>%
          addRasterImage(x = cov_list(),
                         colors = depth.pal,
                         opacity = 0.9,
                         group = "Raster") %>%
          addImageQuery(cov_list(),
                        group = "Raster",
                        project = TRUE) %>%
          addLegend(pal = depth.pal,
                    values = as.vector(values(cov_list())),
                    title = "Depth (m)") %>%
          addLegend(pal = tracks.pal,
                    values = region.my()$id,
                    title = "ID")

      } else if (input$covar == "npp") {

        npp.pal <- colorNumeric(cmocean("algae")(256),
                                domain = as.vector(values(cov_list())),
                                na.color = "transparent")

        proxy %>%
          addRasterImage(x = cov_list(),
                         colors = npp.pal,
                         opacity = 0.9,
                         group = "Raster") %>%
          addImageQuery(cov_list(),
                        group = "Raster",
                        project = TRUE) %>%
          addLegend_decreasing(pal = npp.pal,
                               values = as.vector(
                                 values(
                                   cov_list()
                                 )
                               ),
                               title = "NPP (g C m^-2 d^-1)",
                               decreasing = TRUE) %>%
          addLegend(pal = tracks.pal,
                    values = region.my()$id,
                    title = "ID")
      }
    })  #close observeEvent for leafletProxy mapping
}  #close server

# Run the Shiny app
shinyApp(ui = ui, server = server)
