predict <- function(object, ...)
UseMethod("predict")


#' Predict function for spm
#'@param object An object of class \code{spm}
#'@param data A data frame with current values of \code{Y} and \code{t}
#'@param dt A time interval. Default: 8.
#'@return A data frame consisting of predicted values of \code{Y}, 
#'@return hazard \code{H}, \code{t}, \code{t.next}.
#'@export
#'@examples \dontrun{
#'library(stpm)
#'data <- simdata_discr(N=100, format="long")
#'res <- spm_discrete(data)
#'splitted <- split(data, data$id)
#'df <- data.frame()
#'lapply(1:100, function(i) {df<<-rbind(df,splitted[[i]][dim(splitted[[i]])[1],c("id", "xi", "t1", "y1")])})
#'names(df) <- c("id", "xi", "t", "y")
#'predicted <- predict(object=res, data=df, dt=3)
#'head(predicted)
#'}
#'
predict.spm <- function(object, ...) {
    predict.default(object, data, dt)
}


predict.default <- function(object, data, dt=8) {
    if(!(class(object) %in% c("spm", "spm.discrete", "spm.continuous"))) {
        stop("Wrong class of object. It should be spm")
    }
    
    
    # Next value of Y:
    pred.y <- c()
    haz <- c()
    for(i in 1:length(data$y)) {
        pred.y <- c(pred.y, getNextY.cont(data$y[i], data$t[i], data$t[i]+dt, 
                            object$cmodel$a, object$cmodel$f1, 
                            object$cmodel$Q, object$cmodel$f, 
                            object$cmodel$b, object$cmodel$theta))
        haz <- c(haz, mu(data$y[i], object$cmodel$mu0, object$cmodel$b, 
                object$cmodel$Q, object$cmodel$theta, time))
    }
    
    
    
    res <- data.frame(cbind(pred.y = pred.y, 
                            H=haz,
                            t=data$t,
                            t.next=data$t+dt))
    
    return(res)
}
