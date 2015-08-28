testMatrixMultiply <- function(a,b) {
  c <- .Call("testMatrixMultiply", a, b)
  c
}

testMatrixAdd <- function(a,b) {
  c <- .Call("testMatrixAdd", a, b)
  c
}

testMatrixMultiplyByConst <- function(a,b) {
  c <- .Call("testMatrixMultiplyByConst", a, b)
  c
}

testPowerMatrix <- function(a,b) {
  c <- .Call("testPowerMatrix", a, b)
  c
}