library(eulerr)

# Create the Euler diagram
set.seed(42)
euler.plot <- plot(
  euler(c("AUTOTYP" = 99, "WALS" = 168, "Lexibank" = 20, "PHOIBLE" = 22, "WALS&AUTOTYP" = 15, "WALS&Lexibank" = 6, "WALS&PHOIBLE" = 19),
        shape = "ellipse"), quantities = T, fill = c("lavender","gainsboro","lavenderblush1","bisque1"), labels = list(font = 2))

# Save the plot as a png file
png("Fig1_typlink_logical_euler_diagram.png", 
    bg = "transparent", 
    res = 300, width = 5.5, height = 5.5, units = "in")
plot(euler.plot)
dev.off()  # Close the png file

# Create the Euler diagram
set.seed(41)
euler.plot <- plot(
  euler(c("AUTOTYP" = 91, "WALS" = 163, "Lexibank" = 18, "PHOIBLE" = 21, "WALS&AUTOTYP" = 15, "WALS&Lexibank" = 6, "WALS&PHOIBLE" = 19),
        shape = "ellipse"), quantities = T, fill = c("lavender","gainsboro","lavenderblush1","bisque1"), labels = list(font = 2))
euler.plot

# Save the plot as a png file
png("Fig1_typlink_statistical_euler_diagram.png", 
    bg = "transparent", 
    res = 300, width = 5.5, height = 5.5, units = "in")
plot(euler.plot)
dev.off()  # Close the png file
