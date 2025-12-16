protect_table <- function(df,
                          method = c("ckm",
                                     "sdc"),
                          freq_col = "freq",
                          threshold = 3,
                          noise_range = -2:2) {
  
  method <- match.arg(method)
  
  library(dplyr)
  library(purrr)
  
  # ============================================================
  # CKM-LIKNANDE BRUSNING (inbyggd funktion)
  # ============================================================
  if (method == "ckm") {
    
    library(digest)
    
    brusa_ckm_like <- function(df, cols = freq_col, brus = noise_range) {
      
      # Skapa deterministisk brusvektor baserat på dimensionerna
      noise_vec <- df %>%
        select(-all_of(cols)) %>%
        pmap_chr(~ paste(..., sep = "_")) %>%
        map_chr(~ {
          h <- digest(.x, algo = "xxhash32", serialize = FALSE)
          idx <- (strtoi(substr(h, 1, 4), base = 16) %% length(brus)) + 1
          as.character(brus[idx])
        }) %>%
        as.integer()
      
      df %>%
        mutate(across(all_of(cols), ~ pmax(.x + noise_vec, 0)))
    }
    
    return(brusa_ckm_like(df))
  }
  
  # ============================================================
  # SDC-TABLE BASERAT SKYDD
  # ============================================================
  library(sdcTable)
  library(sdcHierarchies)
  
  dims <- df %>%
    select(-all_of(freq_col)) %>%
    names()
  
  # Skapa riktiga hierarkier med purrr
  dimList <- dims %>%
    map(~ hier_create(
      root  = "Total",
      nodes = sort(unique(df[[.x]]))
    )) %>%
    set_names(dims)
  
  # Skapa SDC-problem
  sdc <- makeProblem(
    data = df,
    dimList = dimList,
    freqVar = freq_col
  )
  
  if (method == "sdc") {
    
    sdc <- primarySuppression(
      sdc,
      type = "freq",
      maxN = threshold - 1
    )
    
    sdc <- protectTable(sdc, method = "HITAS")
    
    df_sdc <- getInfo(sdc, type = "finalData") %>% 
      mutate(Freq = ifelse(sdcStatus != "s", "x", Freq))
    # %>% 
    #   filter(!if_any(all_of(dims), ~ .x == "Total"))
    
    summary(sdc)
    print("Maskerar celler med sdcStatus != s")
    
    return(df_sdc)
  }
  
  stop("Okänd metod.")
}
