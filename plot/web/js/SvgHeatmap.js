// Creates a heatmap using SVG.
// ??? this might be simpler using a Transform node

svgns = "http://www.w3.org/2000/svg";
xlns = "http://www.w3.org/1999/xlink";

// various color utilities.
function rgbToColor(r,g,b) {
  function f(x) {
    x = x < 0 ? 0 : x;
    x = x > 255 ? 255 : x;
    x = Math.floor(x);
    return x;
  }
  return "rgb(" + f(r) + "," + f(g) + "," + f(b) + ")";
}

function scaleColor(range, x) {
  r = x > 0 ? range[1] : range[0];
  x = Math.abs(x);
  return rgbToColor(r[0] * x, r[1] * x, r[2] * x);
}

/* Constructs an SVG heatmap.
  Args:
  svg - SVG object to add onto
  m, n - dimensions (rows & columns) of the grid (in cells)
  cellSize - size of each square (in pixels) */
function SvgHeatmap(svg, m, n, cellSize) {

  // dimensions
  this.m = m;
  this.n = n;

  // size of each grid square
  this.cellSize = cellSize;

  // save the grid, for future use
  this.grid = [];

  // color range for the heatmap (currently fixed)
  this.colorRange = [[0,0,255], [255,255,0]];

  // Adds the labels to the given SVG object.
  for(var i=0; i<m; i++) {
    this.grid[i] = [];

    for(var j=0; j<n; j++) {

      // add a grid square
      var r = document.createElementNS(svgns, "rect");
      r.setAttributeNS(null, "x", this.cellSize * j);
      r.setAttributeNS(null, "y", this.cellSize * i);
      r.setAttributeNS(null, "width", this.cellSize);
      r.setAttributeNS(null, "height", this.cellSize);
      r.setAttributeNS(null, "fill", "white");
      r.setAttributeNS(null, "stroke", "black");
      r.setAttributeNS(null, "stroke-width", "1");

      // save for updating in the future
      this.grid[i][j] = r;

      svg.appendChild(r);
    }
  }
 
  // Sets the i'th row.
  // i - the row
  // x - array of data
  this.setRow = function(i, x) {
    for(var j=0; j<this.n; j++) {
      var s = this.grid[i][j];
      var c1 = scaleColor(this.colorRange, x[j]);
      s.setAttributeNS(null, "fill", scaleColor(this.colorRange, x[j]));
      s.setAttributeNS(null, "alt", x[j]);
    }
  }
};

