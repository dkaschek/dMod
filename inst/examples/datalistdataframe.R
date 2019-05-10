   
   
   Eample: Interconversion between datalist and data.frame
    
   
   library(dMod)
  
  # * Setting 1: no condition column available but covariates
  df1 <- expand.grid(name = letters[1:2], time = 1:2, value = 1, sigma = 1, lloq = 0, cov1 = LETTERS[3:4], cov2 = LETTERS[5:6], cov3 = 1:2, stringsAsFactors = FALSE)
  dl1 <- as.datalist(df1)
  
  # * Setting 2: convert back to data.frame -> now has condition column available
  df2.1 <- as.data.frame(dl1)
  dl2.1 <- as.datalist(df2.1)
  identical(df2.1, df1)
  #   * From the point on when there is a condition-column, we have identical(x, as.datalist(as.data.frame(x)))
  df2.2 <- as.data.frame(dl2.1)
  dl2.2 <- as.datalist(df2.2)
  identical(df2.1, df2.2)
  
  # * Setting 3: additional options with keep.covariates and split.by
  df3.1 <- as.datalist(df1, split.by = "cov1", keep.covariates = FALSE)
  df3.2 <- as.datalist(df1, split.by = "cov1", keep.covariates = "cov2")
  df3.3 <- as.datalist(df1, split.by = c("cov1", "cov2"), keep.covariates = TRUE)
  
  covariates(df3.1)
  covariates(df3.2)
  covariates(df3.3)
