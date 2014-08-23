# Draws a sequence logo as an SVG file, hopefully
# more compactly than, e.g., the current version
# of webLogo.
# XXX writes raw SVG, which is a bit hokey.

source("git/utils.r")

# from MotIV package:
# get information content profile from PWM
pwm2ic<-function(pwm) {
    npos<-ncol(pwm)
    ic<-numeric(length=npos)
    for (i in 1:npos) {
        ic[i]<-2 + sum(sapply(pwm[, i], function(x) { 
            if (x > 0) { x*log2(x) } else { 0 }
        }))
    }    
    ic
}


# Pads a motif to have at least some number of characters.
pad.motif = function(m, len) {
  if (ncol(m) < len) {
    m = cbind(m,
      matrix(0, nrow=nrow(m), ncol=len - ncol(m)))
  }
  m
}

# Scales a motif by information content.
scale.by.ic = function(m) {
  a = pwm2ic(m) / 2
  t( t(m) * a )
}

default.dna.colors =
  c(A="#00ff00", C="#0000ff", G="#ffa500", T="#ff0000")

# Args:
#   m - motif, as a matrix with one column per letter,
#      and one row per letter.
#   letter.colors - what color to draw each letter
#   title - what to call the layer
#   font.family - the font family to use
# Returns: SVG for that motif, as a string.
motif.svg = function(m,
  letter.colors = default.dna.colors,
  title = "motif PWM logo",
  font.family = "Sans-serif") {

  # dimensions of a letter, in pixels
  letter.width = 17
  letter.height = 21
  letter.base = -2.5

  s = paste0("<svg width=\"300\" height=\"20\" ",
    "xmlns=\"http://www.w3.org/2000/svg\" ",
    "xmlns:svg=\"http://www.w3.org/2000/svg\"> ",
    "<g><title>", title, "</title>")

  # SVG for one letter
  svg.for.letter = function(letter, x, y, rel.height) {
    letter.color = letter.colors[letter]
    paste(
  "\n<text transform=\"matrix(1, 0, 0, ", rel.height, ",",
  letter.width * (x-1), ",",
  letter.base + letter.height * (1-y),
  ")\" fill=\"", letter.color,
  "\" stroke=\"#000000\" stroke-width=\"0\" ",
  "x=\"0\" y=\"0\" id=\"svg_2\" font-size=\"24\" ",
  "font-family=\"Sans-serif\" text-anchor=\"bottom\" ",
  "xml:space=\"preserve\">", letter, "</text>", sep="")
  }

  for(j in 1:ncol(m)) {
    # this is where the current letter starts
    y = 0

    # order the letters, smallest first
    p = rev(order(m[,j], decreasing=TRUE, na.last=TRUE))
    for(i in 1:nrow(m)) {

      if (m[p[i], j] > 0) {
        s = paste0(s,
          svg.for.letter(rownames(m)[ p[i] ],
            j, y,
            m[ p[i] , j ]))
      }

      y = y + m[ p[i] , j ]
    }
  }

  paste0(s, "</g></svg>")
}

# Writes an SVG motif as a file.
write.motif.svg = function(m, output.file,
  width=8, height=2, min.motif.width=15) {

  m1 = pad.motif(scale.by.ic(m), min.motif.width)
  s = motif.svg(m1)
  cat(s, file=output.file)
}

if (FALSE) {
m = matrix(c(
  1, 0.25, 0.8, 0,
  0, 0.25, 0.1, 0,
  0, 0.25, 0.05, 1,
  0, 0.25, 0.05, 0), nrow=4, byrow=TRUE)
rownames(m) = c("A", "C", "G", "T")

m = scale.by.ic(pad.motif(m, 7))
}

