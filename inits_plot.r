

plot_init_Had_multiple_noAA <- function(){
    name_style<<-"noAA"

    paper<<-c(7.2,3)
    margins<<-c(0,0,0,0)

    outer_cut<<-10
    inner_cut<<-33

    xAusschnitt<<-c(-168.75,191.75)    
    yAusschnitt<<-c(-60,90)

    greyLand<<-TRUE

    asp<<-1
    pointsize<<-0.3
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.1

    region<<-NA

    season_auswahl<<-1:5
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    color_legend<<-"seperate"
    layout_mat<<-c(NA)

    cexIndex<<-2
    #indexTopRight<<-c("a","A","b","B","c","C","d","D","e","E")
    indexTopRight<<-c("")
    posTopRight<<-c(172,80)

    indexBottomLeft<<-c("MAM","MAM","JJA","JJA","SON","SON","DJF","DJF","Annual","Annual")
    posBottomLeft<<-c(-175,-45)

    nbcol<<-101

    closePlot<<-TRUE
}

plot_init_multi_SH <- function(){
    name_style<<-"SH"

    paper<<-c(8,1.55)
    margins<<-c(0,0,0,5)

    outer_cut<<-10
    inner_cut<<-33

    greyLand<<-TRUE

    yAusschnitt<<-c(-65,0)
    xAusschnitt<<-c(-120,240)
    asp<<-1
    pointsize<<-0.3
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.1

    region<<-NA

    season_auswahl<<-1:5
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    cexIndex<<-2
    indexTopRight<<-c("")
    posTopRight<<-c(172,80)

    indexBottomLeft<<-c("")
    posBottomLeft<<-c(-175,-45)

    col_row<<-c(1,1)
    color_legend<<-"right"
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")

    nbcol<<-101

    border_col<<-"black"
    land_col<<-NA #rgb(0.5,0.5,0.5,0.5)
    water_col<<-rgb(0,0.2,0.8,0.0)

    closePlot<<-TRUE
}

plot_init_EU <- function(){
    name_style<<-"EU"

    paper<<-c(6,3.5)
    margins<<-c(0,0,0,0)

    outer_cut<<-10
    inner_cut<<-1
    
    yAusschnitt<<-c(25,75)
    xAusschnitt<<-c(-13,80)
    asp<<-1
    pointsize<<-0.1
    pch_points<<-c(1,NA,0.25,0.25)

    greyLand<<-FALSE

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.03

    region<<-NA

    season_auswahl<<-c(1,2,3,4,5)
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    color_legend<<-"seperate"
    layout_mat<<-c(NA)

    cexIndex<<-1
    indexTopRight<<-c(NA)
    posTopRight<<-c(175,80)

    indexBottomLeft<<-c(NA)
    posBottomLeft<<-c(-155,-45)
    
    nbcol<<-101

    border_col<<-"black"
    land_col<<-NA #rgb(0.5,0.5,0.5,0.5)
    water_col<<-rgb(0,0.2,0.8,0.0)

    closePlot<<-TRUE
}