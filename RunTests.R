library(RUnit)

test.suite <- defineTestSuite("Sumrep Test Suite",
                              "Tests",
                              testFileRegexp="Test.R")
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result, showDetails=TRUE)
