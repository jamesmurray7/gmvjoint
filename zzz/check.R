.CHECK = function(){
cat("Checking gmvjoint   **************")
cat("\n")
cat("Run as CRAN (y/n)?")
as.cran = readline()
cat("\n")
cat("Run with timings (y/n)?")
w.timings = readline()
cat("\n")
cat("Run dontttest{...} examples (y/n)?")
dont.test = readline()
cat("\n")
cat("Run with valgrind enabled (y/n)?")
val.grind = readline()
cat("\n")
cat("Run with gctorture enabled (y/n)?")
gct = readline()
cat("\n\n")

if(tolower(w.timings) == 'y') w.timings = '--timings' else w.timings = NULL
if(tolower(as.cran) == 'y') as.cran = '--as-cran' else as.cran = NULL
if(tolower(dont.test) == 'y') dont.test = '--run-donttest' else dont.test = NULL
if(tolower(val.grind) == 'y') val.grind = '--use-valgrind' else val.grind = NULL
if(tolower(gct) == 'y') gct = '--use-gct' else gct = NULL

args = c(as.cran, dont.test, w.timings, val.grind, gct)
devtools::check(args = args)
}
