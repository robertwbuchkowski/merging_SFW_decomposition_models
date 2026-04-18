init_millennial_state <- function(Detritivore = F, Earthworm = F, HerbNem = F) {
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
    if(Earthworm){
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
        Earthworm = 10
      )
    }else{
      if(HerbNem){
        c(C_shoot = 100,
          C_root = 100,
          Litter  = 200,
          CWD     = 0,
          Organic = 100,
          DOM     = 10,
          MIC     = 0.1,
          P       = 400,
          L       = 10,
          A       = 1000,
          M       = 3000,
          B       = 40,
          RootHerb = 0.01
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
  }
}