

plot_init_midlat <- function(){
    paper<<-c(8,2)
    margins<<-c(0,0,0,5)

    outer_cut<<-9
    inner_cut<<-2

    yAusschnitt<<-c(0,90)
    xAusschnitt<<-c(-180,180)
    asp<<-1
    pointsize<<-1
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.1

    region<<-NA

    season_auswahl<<-1:5
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    color_legend<<-"right"
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

plot_init_EU <- function(){
    paper<<-c(8,5)
    margins<<-c(2,3,2,5)

    outer_cut<<-0
    inner_cut<<-0
    
    yAusschnitt<<-c(20,80)
    xAusschnitt<<-c(-30,80)
    asp<<-1
    pointsize<<-0.44
    pch_points<<-c(1,NA,0.25,0.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.03

    region<<-NA

    season_auswahl<<-c(1,2,3,4,5)
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    color_legend<<-"right"
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

plot_init_Had <- function(){
    paper<<-c(8,3.5)
    margins<<-c(0,0,0,5)

    outer_cut<<-9
    inner_cut<<-2

    yAusschnitt<<-c(-80,80)
    xAusschnitt<<-c(-180,180)
    asp<<-1
    pointsize<<-1
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.1

    region<<-NA

    season_auswahl<<-1:5
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    color_legend<<-"right"
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

plot_init_Had_multiple <- function(){
    paper<<-c(7.2,3.5)
    margins<<-c(0,0,0,0)

    outer_cut<<-9
    inner_cut<<-0.3

    yAusschnitt<<-c(-90,90)
    xAusschnitt<<-c(-168,192)
    asp<<-1
    pointsize<<-1
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
    subIndex<<-c("a","b")
}

plot_init_multi_midlat <- function(){
    paper<<-c(7,2)
    margins<<-c(0,0,0,5)

    outer_cut<<-9
    inner_cut<<-2

    yAusschnitt<<-c(0,90)
    xAusschnitt<<-c(-168,192)
    asp<<-1
    pointsize<<-1
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
    subIndex<<-c("a","b")
}