#!/opt/R/4.3.1/lib/R/bin//R


#import data
selected_spots <- read.csv("./data/spot_coords/out1.csv")
total_spots <- read.csv("./data/spot_coords/Spatial-Projection.csv")

#Select selected spot coords
selected_coords <- subset(total_spots, Barcode %in% selected_spots$Barcode) # get coords of selected spots
selected_coords <- cbind(selected_spots$out1, selected_coords) # add in region names to coords
selected_coords <- subset(selected_coords, select = -Barcode) # remove temp barcodes
selected_coords <- selected_coords[order(selected_coords$`selected_spots$out1`),] # group regions together



#rename region names to spot names
index =1
name = 1
current <- selected_coords$`selected_spots$out1`[1]
for (spot in selected_coords$`selected_spots$out1`){
  if (spot != current ){
    current <- spot
    name = 1
    selected_coords$`selected_spots$out1`[index] <- paste0(spot, "_spot_", name)
  }
  else {
    selected_coords$`selected_spots$out1`[index] <- paste0(spot, "_spot_", name)
  }
  
  index = index + 1
  name = name + 1
  current <- spot
}

#export
write.csv(selected_coords, file = "./data/spot_coords/regions_coords.csv", row.names = FALSE)
