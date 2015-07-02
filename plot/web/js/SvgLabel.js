// Creates a row of labels using SVG.
// ??? possibly using some framework such as D3 would simplify this
// Wishlist:
// - buttons changing colors on click

svgns = "http://www.w3.org/2000/svg";
xlns = "http://www.w3.org/1999/xlink";

function SvgLabel(n, rowHeight, width) {

  // number of labels
  this.n = n;

  // height of each row
  this.rowHeight = rowHeight;

  // the group node containing this
  this.g = document.createElementNS(svgns, "g");
  this.g.setAttributeNS(null, "x", 0);
  this.g.setAttributeNS(null, "y", 0);
  this.g.setAttributeNS(null, "width", width);
  this.g.setAttributeNS(null, "height", n * rowHeight);
  this.g.setAttributeNS(null, "clip-path", "rect(0,0,"
    + n * rowHeight + "," + width + ")");

  // the SVG Text objects (so as to be able to update them)
  this.svgText = [];

  // similarly for the SVG Link objects
  this.svgLink = [];

  // Adds the labels to the given SVG object.
  for(var i=0; i<n; i++) {
    // background rectangle
    var r = document.createElementNS(svgns, "rect");
    r.setAttributeNS(null, "x", 0);
    r.setAttributeNS(null, "y", this.rowHeight * i);
    r.setAttributeNS(null, "width", width);
    r.setAttributeNS(null, "height", this.rowHeight);
    r.setAttributeNS(null, "fill", i % 2 ? "#fff" : "#bbb");

    // the actual text
    var t = document.createElementNS(svgns, "text");
    t.setAttributeNS(null, "x", 0);
    // XXX there's probably a better way to set this location
    t.setAttributeNS(null, "y", this.rowHeight * (i + 0.75));
    t.setAttributeNS(null, "clip-path",
      "rect(0,0," + width + "," + this.rowHeight + ")");

    // the hyperlink
    var link = document.createElementNS(svgns, "a");
    link.appendChild(r);
    link.appendChild(t);

    // save these for future use
    this.svgText[i] = t;
    this.svgLink[i] = link;

    this.g.appendChild(link);
  }
 
  // Sets the i'th text object.
  this.setText = function(i, s) {
    this.svgText[i].textContent = s;
  }

  // Sets the i'th link href and title (hover text.)
  this.setLink = function(i, href, title) {
    this.svgLink[i].setAttributeNS(xlns, 'xlink:href', href);
    this.svgLink[i].setAttributeNS(xlns, 'xlink:title', title);
  }
};

