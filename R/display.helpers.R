#' @keywords internal
format_windows <- function(windows, .Object) {
  if (is(windows, 'character'))
    windows = BiocGenerics::unlist(parse.grl(windows, seqlengths(seqinfo(.Object))))

  if (is(windows, 'Seqinfo'))
    windows = gUtils::si2gr(windows)

  if (is(windows, 'GRangesList'))
    windows = BiocGenerics::unlist(windows)

  windows = windows[width(windows)>0]  ## otherwise will get non-matching below

  if (is.null(windows$col))
    windows$col = 'gray98'

  #    if (is.list(windows))
  #        windows = do.call('GRangesList', windows)

  ## collapse and match metadata back to original
  tmp = reduce(gr.stripstrand(windows))
  ix = gUtils::gr.match(tmp, windows)
  values(tmp) = values(windows)[ix, ]
  windows = tmp
  ##    if (!inherits(windows, 'GRangesList')) ## GRangesList windows deprecated
  ##        windows = GenomicRanges::GRangesList(windows)

  if (sum(as.numeric(width(gUtils::grl.unlist(windows))))==0)
  {

    if (length(seqinfo(.Object))) {
      windows = gUtils::si2gr(seqinfo(.Object))
    } else {
      warning("no windows provided and no seqinfo. Drawing blank plot")
      return(GRanges())
    }
  }

  return (gUtils::grl.unlist(windows))
}

#' @keywords internal
prep_defaults_for_plotting <- function(.Object) {

  # layout y coordinates
  sumh = sum(.Object$height + .Object$ygap)
  .Object$height  = .Object$height/sumh
  .Object$ygap  = .Object$ygap/sumh
  #            .Object$ywid  = .Object$ywid/sumh

  if (is.null(.Object$max.ranges))
    .Object$max.ranges = NA

  return(.Object)
}

#' @keywords internal
extract_data_from_tmp_dat <- function(.Object, j, this.windows) {

  tmp.dat = dat(.Object)[[j]]

  if (is.character(tmp.dat))
  {
    if (file.exists(tmp.dat))
    {
      bed.style = FALSE
      if (grepl('(\\.bw)|(\\.bigwig)', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::BigWigFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.wig', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::WIGFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.bed', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::BEDFile(normalizePath(tmp.dat))
        #                            si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
        bed.style = TRUE
      }
      else if (grepl('\\.gff', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::GFFFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.2bit', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::TwoBitFile(normalizePath(tmp.dat))
        si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
      }
      else if (grepl('\\.bedgraph', tmp.dat, ignore.case = TRUE))
      {
        f = rtracklayer::BEDGraphFile(normalizePath(tmp.dat))
        #                           si = tryCatch(seqinfo(f), error = function(tmp.dat) NULL)
        bed.style = TRUE
      }

      if (!bed.style) ## bed.style file objects do not have seqinfo option
      {
        si = tryCatch(seqinfo(f), error = function(x) '')
        pre.filtered = T

        if (!is.character(si)) {
          sel = gUtils::gr.fix(this.windows, si, drop = TRUE)

          if (length(this.windows)>0 & length(sel)==0)
            warning('zero ranges intersecting UCSC file sequences selection. Perhaps "chr" prefix mismatch?')

          tmp.dat = rtracklayer::import(f, selection = sel, asRangedData = FALSE)

          if (is.na(.Object$y.field[j]))
            .Object$y.field[j] = 'score'

        }
      }
      else
      {
        tmp.dat = rtracklayer::import(f, asRangedData = FALSE)

        if (is.na(.Object$y.field[j]))
          .Object$y.field[j] = 'score'
      }
    }
  }
  else if (is(tmp.dat, 'ffTrack'))
  {
    tmp.dat = tmp.dat[this.windows, gr = TRUE]
    .Object$y.field[j] = 'score'
    pre.filtered = T
  }

  if (is.character(tmp.dat)) ## file was not found
  {
    warnings('Track bigwig file not found')
    tmp.dat = GRanges()
  }

  return(list(o=.Object, t=tmp.dat))
}

#' @keywords internal
enforce_max_ranges <- function(.Object, pre.filtered, j, tmp.dat, this.windows) {

  ## enforce max.ranges (everything should be a GRanges at this point)
  if (!pre.filtered & nrow(edgs(.Object)[[j]])==0 & length(.Object@vars[[j]])==0) ## assume any track associated with edges is pre-filtered
  {
    if (length(tmp.dat)>.Object$max.ranges[j])
      if (inherits(tmp.dat, 'GRangesList'))
      {
        vals = values(tmp.dat)
        nm = names(tmp.dat)
        tmp2 = gUtils::grl.unlist(tmp.dat)
        tmp2 = tmp2[gr.in(tmp2, this.windows)]
        ##tmp2 = tmp2[tmp2 %over% this.windows]
        tmp.dat = GenomicRanges::split(tmp2, tmp2$grl.ix)
        values(tmp.dat) = vals[as.numeric(names(tmp.dat)), ]
        names(tmp.dat) = nm[as.numeric(names(tmp.dat))]
      }
    else if (length(tmp.dat)>.Object$max.ranges[j])
    {
      tmp.dat = tmp.dat[gr.in(tmp.dat, this.windows)]
      ##tmp.dat = tmp.dat[tmp.dat %over% this.windows]
      pre.filtered = TRUE
    }
  }

  if (length(tmp.dat)>.Object$max.ranges[j] &
      nrow(edgs(.Object)[[j]])==0 & !.Object@formatting$triangle[j]) ## don't sample if there are edges or triangle
  {
    tmp.dat = sample(tmp.dat, ceiling(.Object$max.ranges[j]))
  }

  return(list(p=pre.filtered, t=tmp.dat))

}

#' @keywords internal
smooth_yfield <- function(.Object, j, tmp.dat) {

  tmp = S4Vectors::runmean(GenomicRanges::coverage(tmp.dat, weight = values(tmp.dat)[, .Object$y.field[j]]),
                           k = floor(.Object$smooth[j]/2)*2+1, endrule = 'constant', na.rm = TRUE)

  if (!is.na(.Object$round[j]))
    tmp = round(tmp, .Object$round[j])

  tmp = as(tmp, 'GRanges')
  tmp = tmp[gr.in(tmp, tmp.dat)]
  ##tmp = tmp[tmp %over% tmp.dat]
  tmp.val = tmp$score
  values(tmp) = values(tmp.dat)[gUtils::gr.match(tmp, tmp.dat), , drop = F]
  values(tmp)[, .Object$y.field[j]] = tmp.val
  tmp.dat = tmp

}

#' @keywords internal
format_yfield_limits <- function(.Object, j, tmp.dat, pre.filtered, this.windows) {

  strand(this.windows) <- rep("*", length(this.windows))

  if (!(.Object$y.field[j] %in% names(values(tmp.dat))))
    stop('y.field missing from input granges')

  y0.global = min(values(tmp.dat)[, .Object$y.field[j]], na.rm = TRUE)
  y1.global = max(values(tmp.dat)[, .Object$y.field[j]], na.rm = TRUE)

  if (!pre.filtered)
    tmp.dat.r <- tmp.dat[gr.in(tmp.dat, this.windows)]
  ##tmp.dat.r <- tmp.dat[tmp.dat %over% this.windows]
  else
    tmp.dat.r = tmp.dat

  val = values(tmp.dat.r)[, .Object$y.field[j]]
  #                          r = range(val[!is.infinite(val)], na.rm = TRUE)

  p.quantile = 0.01
  quant = .Object@yaxis[[j]]$quantile
  if (!is.na(quant))
    p.quantile = pmin(pmax(pmin(quant, 1-quant), 0), 1)

  r = quantile(val[!is.infinite(val)], probs = c(p.quantile, 1-p.quantile), na.rm = TRUE)

  if (is.na(diff(r)))
    r = c(0, 0)

  y0 = .Object@yaxis[[j]]$min
  y1 = .Object@yaxis[[j]]$max
  y0.bar = .Object@yaxis[[j]]$min.bar

  ## adjust y0 if not specified
  if (is.na(y0)) {
    y0 =  pmax(y0.global, r[1] - diff(r)*0.05)
    if (diff(r) == 0 & !is.na(y1)[j])
      y0 = r[1] - diff(c(r[1], y1))*0.05

  }

  ## adjust y1 if not specified
  if (is.na(y1)) {
    y1 = pmin(y1.global, r[2] + diff(r)*0.05)
    if (diff(r) == 0 & !is.na(y0))
      y1 = r[2] + diff(c(y0, r[1]))*0.05
  }

  if (!is.null(y0.bar) && !is.na(y0.bar))
    y0 = min(y0.global, y0.bar)

  return(.Object)
}
