

plot_init_midlat <- function(){
    paper<<-c(8,2)
    margins<<-c(0,0,0,5)

    outer_cut<<-5
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
    mat<<-NA
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

plot_init_midlat()