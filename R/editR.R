editrversion <- "1.0.10"

GetGuideMatch <- function(guide, input.seq)
{
  align.f  <- pairwiseAlignment(pattern=guide, subject=input.seq, type="overlap")

  guide.coord <- list(match = "forward", start = align.f@subject@range@start,
                      end = align.f@subject@range@start + align.f@subject@range@width-1)

  if(align.f@pattern@range@width <= 1){
    stop("Failed to find a forward match of the gRNA -- should it be a reverse complement?")
    return(NA)
  }

  is.match.guide.length <- (guide.coord$end - guide.coord$start + 1) == length(guide)

  if(is.match.guide.length != TRUE){

    startdiff <- 1 - align.f@pattern@range@start
    enddiff <- 20 - (align.f@pattern@range@start + align.f@pattern@range@width-1)
    guide.coord$start <- guide.coord$start + startdiff
    guide.coord$end <- guide.coord$end + enddiff

  }

  return(guide.coord)
}


CreateSangs <- function(peakAmp, basecalls){

  sangs <- as.data.frame(peakAmp)
  names(sangs) <- c("A.area","C.area","G.area","T.area")

  sangs %<>%
    mutate(Tot.area = A.area + C.area + G.area + T.area,
           A.perc = 100*A.area / Tot.area,
           C.perc = 100*C.area / Tot.area,
           G.perc = 100*G.area / Tot.area,
           T.perc = 100*T.area / Tot.area)

  sangs$base.call <- strsplit(x = toString(basecalls@primarySeq), split = "") %>%
    unlist

  sangs$index <- seq_along(sangs$base.call)
  return(sangs)
}


GetNullDistModel <- function(sangs.filt, guide.coord)
{
  sangs.filt %<>% filter(!(index %in% (guide.coord$start:guide.coord$end)) )

  nvals <- list()
  nvals$t <- sangs.filt %>% filter(base.call != "T") %>% select(T.perc) %>% unlist()
  nvals$c <- sangs.filt %>% filter(base.call != "C") %>% select(C.perc) %>% unlist()
  nvals$g <- sangs.filt %>% filter(base.call != "G") %>% select(G.perc) %>% unlist()
  nvals$a <- sangs.filt %>% filter(base.call != "A") %>% select(A.perc) %>% unlist()

  replacement_zaga = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
                       0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
                       0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)

  n.models <- lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = replacement_zaga; message("Replacement vector used for low noise.")}
    tryCatch(gamlss((x)~1, family = ZAGA), error=function(e)
      tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e)
        tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e)
          tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e)
            gamlss((x)~1, family = ZAGA, mu.start = mean(x))
          )
        )
      )
    )
  })

  null.m.params <- lapply(n.models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]

    return(data.frame(mu= mu, sigma = sigma, nu = nu, fillibens = fillibens))
  })

  return(null.m.params)
}


CreateEditingDF <- function(guide.coord, guide, sangs, null.m.params){

  guide.df <- sangs[guide.coord$start:guide.coord$end,]

  guide.df$guide.seq <- guide %>% toString() %>% strsplit(. , "") %>% unlist

  calcBaseProb <- function(params, perc){
    outprobs <- pZAGA(q = perc,
                      mu = params$mu,
                      sigma = params$sigma,
                      nu = params$nu,
                      lower.tail = FALSE)
  }

  guide.df$T.pval <- calcBaseProb(params = null.m.params$t, perc = guide.df$T.perc)
  guide.df$C.pval <- calcBaseProb(params = null.m.params$c, perc = guide.df$C.perc)
  guide.df$G.pval <- calcBaseProb(params = null.m.params$g, perc = guide.df$G.perc)
  guide.df$A.pval <- calcBaseProb(params = null.m.params$a, perc = guide.df$A.perc)

  guide.df$guide.position <- seq_along(guide.df$guide.seq)

  editing.df <- guide.df
  return(editing.df)
}


makeEditingBarPlot <- function(edit.long, null.m.params, base, pval, editing.df){

  cutoff.line <- qZAGA(p = pval, mu = null.m.params$mu,
                       sigma = null.m.params$sigma,
                       nu = null.m.params$nu, lower.tail = FALSE)

  edit.filt <- edit.long %>% filter(guide.seq != base, focal.base == base)

  edit.filt %>% filter(measure == "perc") %>%
    ggplot(aes(x = guide.position, y = value)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(ylim = c(0,100)) +
    scale_x_continuous(breaks = seq_along(editing.df$guide.seq)) +
    annotate(geom = "segment", x = 0, xend = length(editing.df$guide.seq),
             y = cutoff.line, yend = cutoff.line) +
    labs(title = paste0("Percent ", base, " editing"),
         x = "Base postion in guide RNA",
         y = "Percent Total peak area") +
    theme(axis.text.x = element_text(size = 8))

}


run_editr <- function(ab1_file = NULL,
                      guide = NULL,
                      pval_cutoff = 0.01,
                      trim5 = NA,
                      trim3 = NA,
                      guide_is_rc = FALSE,
                      example = FALSE,
                      output_dir = NULL) {

  if(example) {
    ab1_file <- system.file("example.ab1", package = "editR")
    guide <- "CACTGGAATGACACACGCCC"
  }

  if(is.null(ab1_file)) {
    stop("Please provide an .ab1 file path or set example = TRUE")
  }
  if(is.null(guide)) {
    stop("Please provide a guide RNA sequence")
  }

  input.seq <- readsangerseq(ab1_file)
  input.basecalls <- makeBaseCalls(input.seq)
  input.peakampmatrix <- peakAmpMatrix(input.basecalls)
  sangs <- CreateSangs(input.peakampmatrix, input.basecalls)

  guide_dna <- DNAString(guide)
  if(guide_is_rc) {
    guide_dna <- reverseComplement(guide_dna)
  }

  guide.coord <- GetGuideMatch(guide_dna, input.basecalls@primarySeq)

  if(is.na(trim5) & is.na(trim3)){
    sangs.filt <- sangs %>% filter(index > 20)
    peakTotAreaCutoff <- mean(sangs.filt$Tot.area)/10
    sangs.filt %<>% filter(Tot.area > peakTotAreaCutoff)
  } else {
    if(is.na(trim5)) {
      stop("A value to trim the 5' end is required when specifying trim ranges")
    }
    if(is.na(trim3)) {
      stop("A value to trim the 3' end is required when specifying trim ranges")
    }
    sangs.filt <- sangs[trim5:trim3, ]
  }

  null.m.params <- GetNullDistModel(sangs.filt, guide.coord)

  editing.df <- CreateEditingDF(guide.coord, guide_dna, sangs, null.m.params)

  base.info <- NULL
  if(!is.null(sangs.filt) && !is.null(null.m.params) && !is.null(pval_cutoff)) {
    avg.base <- sangs.filt %>% gather(key = focal.base, value = value,
                                      A.area:T.area, A.perc:T.perc) %>%
      separate(col = focal.base, into = c("focal.base", "measure")) %>%
      spread(key = measure, value = value) %>%
      filter(base.call == focal.base) %>%
      group_by(focal.base) %>%
      summarize(avg.percsignal = mean(perc),
                avg.areasignal = mean(area))

    crit.vals <- c(
      a = qZAGA(p = pval_cutoff, mu = null.m.params$a$mu,
                sigma = null.m.params$a$sigma,
                nu = null.m.params$a$nu,
                lower.tail = FALSE),
      c = qZAGA(p = pval_cutoff, mu = null.m.params$c$mu,
                sigma = null.m.params$c$sigma,
                nu = null.m.params$c$nu,
                lower.tail = FALSE),
      g = qZAGA(p = pval_cutoff, mu = null.m.params$g$mu,
                sigma = null.m.params$g$sigma,
                nu = null.m.params$g$nu,
                lower.tail = FALSE),
      t = qZAGA(p = pval_cutoff, mu = null.m.params$t$mu,
                sigma = null.m.params$t$sigma,
                nu = null.m.params$t$nu,
                lower.tail = FALSE)
    )

    fil <- lapply(null.m.params, FUN = function(x){x$fillibens})
    filvec <- c(a = fil$a, c = fil$c, g = fil$g, t = fil$t)

    mul <- lapply(null.m.params, FUN = function(x){x$mu})
    mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)

    base.info <- data.frame(avg.base, crit.perc.area = crit.vals, mu = mulvec, fillibens = filvec)
  }

  report_data <- list(
    null.m.params = null.m.params,
    guide = guide_dna,
    guide.coord = guide.coord,
    editing.df = editing.df,
    p.val.cutoff = pval_cutoff,
    input.seq = input.seq,
    input.basecalls = input.basecalls,
    sangs = sangs,
    sangs.filt = sangs.filt,
    base.info = base.info,
    editrversion = editrversion
  )

  if(!is.null(output_dir)) {
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    report_file <- file.path(output_dir, "editr_report.html")
    rmarkdown::render(
      system.file("report.Rmd", package = "editR"),
      output_file = report_file,
      params = report_data,
      envir = new.env(parent = asNamespace("editR")),
      intermediates_dir = tempdir()
    )
    message("Report generated at: ", report_file)
  }

  invisible(report_data)
}