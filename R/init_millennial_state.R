init_millennial_state <- function(Detritivore = F) {
  # Units: g C m^-2
  if(Detritivore){
    c(
      Litter  = 200,
      CWD     = 1000,
      Organic = 2000,
      DOM     = 10,
      MIC     = 0.1,
      P       = 400,
      L       = 10,
      A       = 1000,
      M       = 3000,
      B       = 40,
      Detritivore = 10
    )
  }else{
    c(
      Litter  = 200,
      CWD     = 1000,
      Organic = 2000,
      DOM     = 10,
      MIC     = 0.1,
      P       = 400,
      L       = 10,
      A       = 1000,
      M       = 3000,
      B       = 40
    )
  }
}