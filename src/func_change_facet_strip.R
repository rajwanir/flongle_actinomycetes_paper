require(gtable)
require(grid)
require(gridExtra)


change_strip_color = function(p1,dummy) {

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(dummy)

gtable_select <- function (x, ...) 
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip_r", g2$layout$name)
g2$layout$r[panels] <- g2$layout$r[panels] + 1
g2$layout$l[panels] <- g2$layout$l[panels] + 1


stript <- grepl(pattern="strip-r", g2$layout$name)



new_strips <- gtable_select(g2, panels | strips | stript)
grid.newpage()
grid.draw(new_strips)

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}
## ideally you'd remove the old strips, for now they're just covered
new_plot <- gtable_stack(g1, new_strips)

return(new_plot)
}