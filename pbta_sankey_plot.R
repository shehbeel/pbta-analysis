library(plotly)

fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  node = list(
    label = c("PBTA", "ATRT", "Cranio", "Ependy", "Gangli", "GNT", "HGAT", "LGAT", "Medullo", "Schwan", "Tert",
              "G1", "G2"),
    color = c("black", "blue", "red", "green", "purple", "orange", "lightblue", "pink", "lightgreen", "violet", "yellow"),
    pad = 15,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = c(0,0,0,0,0,0,0,0,0,0,4,4),
    target = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
    value =  c(17,21,36,19,3,34,102,35,2,1,10,9)
  )
)
fig <- fig %>% layout(
  title = "Pediatric Brain Tumor Subtyping by miRNA",
  font = list(
    size = 10
  )
)

fig
