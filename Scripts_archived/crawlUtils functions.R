### Modifications of {crawlUtils} functions to source in for now


# cu_crw_argos <- function(data_list, bm=FALSE, fixPar=NULL){
#   i <- datetime <- type <- const <- NULL #handle 'no visible binding...'
#   p <- progressr::progressor(length(data_list))
#   fits <- foreach(i=1:length(data_list), .packages="sf") %dorng% {
#     dat <- data_list[[i]] %>% dplyr::arrange(datetime)
#     alsg <- all(dat$type%in%c("Argos_ls","FastGPS","known"))
#     akfg <- all(dat$type%in%c("Argos_kf","FastGPS","known"))
#     if(!alsg & !akfg) {
#       pct_alsg <- (dat %>% filter(type=="Argos_ls") %>% nrow()) / nrow(dat)
#       pct_akfg <- (dat %>% filter(type=="Argos_kf") %>% nrow()) / nrow(dat)
#       argos_pick <- if_else(pct_alsg > pct_akfg, "LS", "KF")
#       if (argos_pick == "LS") {
#         dat <- dat %>% filter(type %in% c("Argos_ls","FastGPS","known"))
#       }
#       if (argos_pick == "KF") {
#         dat <- dat %>% filter(type %in% c("Argos_kf","FastGPS","known"))
#       }
#       warning("Animal ", i, " has both LS and KF Argos location types or other unknown types!\n",
#               "Keeping ", argos_pick, " becuase ", argos_pick, " represents ",
#               "the larger percentage of the observations.")
#     }
#     if(alsg & akfg) alsg <- FALSE
#     if(alsg){
#       err.model <- list(x =  ~0+ln.sd.x+aq0+aqA+aqB)
#       fixPar <- c(1,NA,NA,NA,NA,NA)
#       constr <- list(
#         lower=c(rep(0,3), -Inf, log(-log(1-1.0e-4))),
#         upper=c(rep(Inf,3), Inf, log(-log(1.0e-4)))
#       )
#       theta <- c(rep(log(1.2),3),9,0.5)
#       prior <- function(x) {
#         dnorm(x[5], -4, 2, log = TRUE)
#       }
#
#     } else{
#       err.model <- list(x =  ~0+ln.sd.x, y = ~0+ln.sd.y, rho= ~error.corr)
#       fixPar <- c(1,1,NA,NA)
#       constr <- list(
#         lower=c(-Inf, log(-log(1-1.0e-4))),
#         upper=c(Inf, log(-log(1.0e-4)))
#       )
#       theta <- c(9,.5)
#     }
#     if(bm){
#       fixPar[length(fixPar)] <- log(-log(1.0e-4))
#       constr$lower <- const$lower[-length(const$lower)]
#       constr$lower <- const$lower[-length(const$lower)]
#       theta <- theta[-length(theta)]
#     }
#     # Fit ctcrw model
#     suppressMessages(
#       out <- crawl::crwMLE(
#         mov.model = ~1, err.model = err.model, data = dat, Time.name="datetime",
#         fixPar = fixPar, constr = constr, theta = theta, prior = prior,
#         control = list(maxit=2000, trace = 2), initialSANN = list(maxit=1500, temp=10),
#         attempts=10, method = "L-BFGS-B", retrySD = 2)
#     )
#     p()
#     out
#   }
#   return(fits)
# }

#-----------------------------

# cu_batch_predict <- function(fit_list, predTime, barrier=NULL, vis_graph=NULL, crs){
#   i <- NULL #handle 'no visible binding...'
#   p <- progressr::progressor(length(fit_list))
#   plist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr")) %dorng% {
#     pred <- crawl::crwPredict(fit_list[[i]], predTime=predTime, return.type="flat")
#     if(!is.null(barrier) & !is.null(vis_graph)){
#       if (!requireNamespace("pathroutr", quietly = TRUE)) stop("Please install {pathroutr}: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
#       pred <- pred %>% crw_as_sf(ftype="POINT", locType="p", crs = crs)
#       pred <- pred %>% pathroutr::prt_trim(barrier)
#       fix <- pathroutr::prt_reroute(pred, barrier, vis_graph, blend=FALSE)
#       pred <- pathroutr::prt_update_points(fix, pred)
#     }
#     p()
#     pred
#   }
#   return(plist)
# }

#-----------------------------

# crw_as_sf = function(data, ftype, locType = c("p", "o", "f"), group = NULL, crs, ...) {
#   locType <- enquo(locType)
#   stopifnot(!missing(ftype), ftype %in% c("POINT", "LINESTRING"))
#   crw_crs <- crs
#   # if (is.null(crw_crs) || is.na(crw_crs))
#   #   crw_crs <- attr(data, "proj4")
#   if (ftype == "POINT" && is.null(group)) {
#     data <- crw_as_tibble(data) %>% dplyr::filter(.data$locType %in%
#                                                     !!locType) %>% dplyr::arrange(.data$TimeNum) %>%
#       sf::st_as_sf(coords = c("mu.x", "mu.y"))
#     data = data %>% sf::st_set_crs(crw_crs)
#   }
#   if (ftype == "POINT" && !is.null(group)) {
#     warning("group argument not applicable for 'POINT' type. ignoring")
#   }
#   if (ftype == "LINESTRING" && is.null(group)) {
#     data <- crw_as_tibble(data) %>% dplyr::filter(.data$locType %in%
#                                                     !!locType) %>% dplyr::arrange(.data$TimeNum) %>%
#       sf::st_as_sf(coords = c("mu.x", "mu.y"))
#     data = data %>% sf::st_set_crs(crw_crs)
#     data = data %>% summarise(id = 1, do_union = FALSE) %>%
#       sf::st_cast("LINESTRING")
#   }
#   if (ftype == "LINESTRING" && !is.null(group)) {
#     data <- crw_as_tibble(data) %>% dplyr::filter(.data$locType %in%
#                                                     !!locType) %>% dplyr::arrange(.data$TimeNum) %>%
#       sf::st_as_sf(coords = c("mu.x", "mu.y"))
#     data = data %>% sf::st_set_crs(crw_crs)
#     data = data %>% dplyr::group_by(group) %>% dplyr::summarise(do_union = FALSE) %>%
#       sf::st_cast("LINESTRING")
#   }
#   return(data)
# }

#-----------------------------

cu_crw_sample <- function(size=8, fit_list, predTime, barrier=NULL, vis_graph=NULL){
  i <- j <- NULL #handle 'no visible binding...'
  p <- progressr::progressor(length(fit_list))
  slist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"))%dorng%{
    simObj <- crawl::crwSimulator(fit_list[[i]], parIS = 0, predTime=predTime)
    out <- foreach(j=1:size)%do%{
      samp <- crawl::crwPostIS(simObj, fullPost = FALSE)

      samp$datetime <- lubridate::as_datetime(simObj$TimeNum*3600)  #modify based on helper function from Johnson et al 2021

      if(!is.null(barrier) & !is.null(vis_graph)){
        if (! requireNamespace("pathroutr", quietly = TRUE)) stop("Please install {pathroutr}: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
        samp <- samp %>% crawl::crw_as_sf(ftype="POINT", locType="p")
        samp <- samp %>% pathroutr::prt_trim(barrier)
        fix <- pathroutr::prt_reroute(samp, barrier, vis_graph, blend=FALSE)
        samp <- pathroutr::prt_update_points(fix, samp) %>% dplyr::mutate(rep=j)
      }
      samp
    }
    p()
    out
  }
  return(slist)
}

