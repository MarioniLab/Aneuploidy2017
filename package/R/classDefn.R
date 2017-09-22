setClass("ploidytest",
         contains = "matrix",
         slots = c(counts = "matrix",
                   cpm = "matrix",
                   chrs = "character",
                   cellNames = "character",
                   cellGroups = "character",
                   scores = "data.frame",
                   aneuploidies = "data.frame",
                   params = "list",
                   knownAneu = "data.frame",
                   metrics = "list")
)
