# Load the required libraries
library(tidyverse)
library(grid)

#colors
colors <- c("#55AA56", "#438244", "#346b34", "#F04442")

#data to plot
my_species <- c('ERR2929116', 'ERR2929116', 'ERR2929116', 'ERR2929116', 'ERR2929117', 'ERR2929117', 'ERR2929117', 'ERR2929117', 'ERR2929118', 'ERR2929118', 'ERR2929118', 'ERR2929118', 'ERR2929122', 'ERR2929122', 'ERR2929122', 'ERR2929122', 'ERR2929123', 'ERR2929123', 'ERR2929123', 'ERR2929123', 'ERR2929124', 'ERR2929124', 'ERR2929124', 'ERR2929124')
my_species <- factor(my_species)
my_values <- c(23184254, 1899721, 239505, 286784, 27807848, 2562758, 268152, 394315, 24159449, 2086408, 252929, 322817, 28847983, 2406690, 292088, 380744, 20793725, 1646617, 1108062, 569792, 29799326, 2338747, 308449, 364107)
my_percentage <- c(90.53, 7.42, 0.94, 1.12, 89.61, 8.26, 0.86, 1.27, 90.07, 7.78, 0.94, 1.20, 90.35, 7.54, 0.91, 1.19, 86.22, 6.83, 4.59, 2.36, 90.82, 7.13, 0.94, 1.11)

#build the graph
category <- c(rep(c("S","D","F","M"),c(1)))
category <-factor(category)
category = factor(category,levels(category)[c(4,1,2,3)])
df = data.frame(my_species,my_percentage,my_values,category)

figure <- ggplot() + 
  
	geom_bar(aes(y = my_percentage, x = my_species, fill = category), position = position_stack(reverse = TRUE), data = df, stat="identity", width=0.75) + 
	coord_flip() + 
	theme_gray(base_size = 8) + 
	scale_y_continuous(labels = c("0%","20%","40%","60%","80%","100%"), breaks = c(0,20,40,60,80,100), expand = c(0, 0)) +
	scale_fill_manual(values = colors,labels =c("Both surviving",
                                                 "Forward only surviving",
                                                 "Reverse only surviving",
                                                 "Dropped")) +   
	ggtitle("Run") +
	xlab("") + 
	ylab("\n%reads") + 

	theme(legend.title = element_blank(),
		legend.text = element_text(size = rel(1.2)), 
		panel.background = element_rect(color="#FFFFFF", fill="white"),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		plot.title = element_text(hjust = -0.15, size = rel(1.8), face = "bold"),
		axis.title.x = element_text(face="bold", size=rel(1.8)),
		axis.text.y = element_text(colour = "black", size = rel(1.66)),
		axis.text.x = element_text(colour = "black", size = rel(1.66)))

  
  figure

  
#ggsave("trimmomatic_stats.pdf", plot=last_plot(), device="pdf", dpi="retina")
