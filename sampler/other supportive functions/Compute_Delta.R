compute_Delta <- function(PIPs) {
  2 * (sum(PIPs[PIPs >= 0.5]) + sum(1 - PIPs[PIPs < 0.5]))
}