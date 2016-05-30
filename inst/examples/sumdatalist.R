  
  # Start with two data frames
  mydata1 <- data.frame(
    name = "A",
    time = 0:1,
    value = 1:2,
    sigma = .1,
    compound = c("DEM", "APAP"),
    dose = "0.1"
  )
  
  mydata2 <- data.frame(
    name = "A",
    time = 0:1,
    value = 3:4,
    sigma = .1,
    compound = c("APAP", "DCF"),
    dose = "0.1"
  )
 
  # Create datalists from dataframes
  data1 <- as.datalist(mydata1, split.by = c("compound", "dose")) 
  data2 <- as.datalist(mydata2, split.by = c("compound", "dose")) 
  
  # Direct sum of datalists
  data <- data1 + data2
   print(data)
  
  # Check the condition.grid (if available)
  condition.grid <- attr(data, "condition.grid")
   print(condition.grid)

