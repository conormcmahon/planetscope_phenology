# This code sourced from John Baums at: https://gist.github.com/johnbaums/61c8062938e05a4c6b92
# Original discussion of this code stemmed from a StackOverflow question here: https://stackoverflow.com/questions/33700755/how-can-i-find-the-pixel-wise-standard-deviation

gdal_sd <- function(infile, outfile, quiet=TRUE) {
  require(rgdal)
  # infile:   The multiband raster file (or a vector of paths to multiple raster
  #           files) for which to calculate cell standard deviations.
  # outfile:  Path to raster output file.
  # quiet:    Logical. Should gdal_calc.py output be silenced?
  gdal_calc <- Sys.which('gdal_calc.py')
  if(gdal_calc=='') stop('gdal_calc.py not found on system.')
  if(file.exists(outfile)) stop('outfile already exists.')
  nbands <- sapply(infile, function(x) nrow(attr(GDALinfo(x), 'df')))
  if(length(infile) > 26 || nbands > 26) stop('Maximum number of inputs is 26.')
  if(length(nbands) > 1 & any(nbands > 1)) 
    warning('One or more rasters have multiple bands. First band used.')
  
  if(length(infile)==1) {
    inputs <- paste0('-', LETTERS[seq_len(nbands)], ' ', infile, ' --', 
                     LETTERS[seq_len(nbands)], '_band ', seq_len(nbands), collapse=' ')
    n <- nbands
  } else {
    inputs <- paste0('-', LETTERS[seq_along(nbands)], ' ', infile, ' --', 
                     LETTERS[seq_along(nbands)], '_band 1', collapse=' ')
    n <- length(infile)
  }
  
  message('Calculating standard deviation and writing to ', basename(outfile))
  cmd <- 'python %s %s --outfile=%s --calc="std([%s], 0, ddof=1)" --co="COMPRESS=LZW"'
  out <- system(
    sprintf(cmd, gdal_calc, inputs, outfile, 
            paste0(LETTERS[seq_len(n)], collapse=',')),
    show.output.on.console=!quiet, intern=TRUE
  )
  if(any(grepl('Error', out))) stop(out, call.=FALSE) else NULL
}