\name{predict_fpqr}
\alias{predict_fpqr}
\title{Prediction for a function-on-function linear quantile regression model based on functional partial quantile regression}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted function-on-function linear quantile regression model in the output of \code{\link{fpqr}}.}
\usage{
predict_fpqr(object, xnew)
}
\arguments{
   \item{object}{An output object obtained from \code{\link{fpqr}}.}
  \item{xnew}{A matricx consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fpqr}}.}

}
\value{
A matrix of predicted values of the functional response variable for the given set of new functional predictor \code{xnew}.
}

\author{Muge Mutis, Ufuk Beyaztas, Filiz Karaman, and Han Lin Shang}
\examples{
gpx <- (1:50)/50
gpy <- (1:60)/60
data <- dgp(n = 250, gpy = gpy, gpx = gpx, sd.error = 0.01)
Y <- data$y
X <- data$x

data.test <- dgp(n = 100, gpy = gpy, gpx = gpx, sd.error = 0.01)
X.test <- data.test$x
Y.test <- data.test$y

fpqr.model.li <- fpqr(y=Y, x=X, h=5, tau=0.5, nby=20,
nbx=20, gpy=gpy, gpx=gpx, qc.type = "li")

predictions <- predict_fpqr(object = fpqr.model.li, xnew = X.test)
}
