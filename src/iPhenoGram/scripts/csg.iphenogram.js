/**
* @module csg/iphenogram
* @license
* Copyright (c) Marylyn Ritchie 2015
*
* This file is part of iPhenoGram.
*
* iPhenoGram is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* iPhenoGram is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with iPhenoGram.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 *@fileOverview
 *@version 1.0
 *
 * @namespace csg.iphenogram
 */

/**
* Implements the jQuery widget factory interface.
*/
$.widget('csg.iphenogram',{


	/**
	* @property {object} options The jQuery widget factory options
	* @property {number} options.height Overall height in pixels
	* @property {number} options.width Overall width in pixels
	* @property {number} options.xmax Horizontal max coordinate
	* @property {number} options.ymax Vertical max coordinate
	* @property {number} options.chromheight Height of a chromosome in coordinate scale
	* @property {number} options.chromwidth Width of a chromosome in coordinate scale
	* @property {number} options.circlesize Diameter of phenotype circles in coordinate scale
	* @property {object} options.margin Margin sizes in coordinate scale for top, bottom, right and left
	* @property {number} options.circleradius Radius of phenotype circles in coordinate scale
	* @property {number} options.linelength Length of line from chromosome to phenotype circle
	* @property {number} options.chromboxheight Height including margin for chromosome
	* @property {number} options.chromsperrow Chromosomes per row
	* @property {number} options.chromrows Number of rows of chromosomes
	* @property {number} options.capheight Height of chromosome cap
	* @property {array} options.chroms Array of objects containing chromosome information
	* @property {string} options.chroms.name Chromosome name
	* @property {number} options.chroms.size Chromosome size in bp
	* @property {number} options.chroms.centromere_start Chromosome centromere start location in bp
	* @property {number} options.chroms.centromere_end Chromosome centromere end location in bp
	* @property {array} options.phenos Array of objects containing phenotype information
	* @property {string} options.phenos.chrom Chromosome name for phenotype position
	* @property {number} options.phenos.position Position of phenotype on chromosome in bp
	* @property {string} options.phenos.id Position identifier on chromosoome
	* @property {string} options.phenos.pheno Phenotype name
	* @property {number} options.phenos.endPosition Optional end position of phenotype region on chromosome in bp
	* @property {string} options.phenos.posColor Optional color for line (or region) across chromosome
	* @property {string} options.phenos.group Optional group for defining shape on the plot
	* @property {array} options.cytoBands Cytogenetic bands on chromosomes
	* @property {string} options.cytoBands.type Type of band
	* @property {number} options.cytoBands.startBand Start position in bp
	* @property {number} options.cytoBands.finishBand End position in bp
	* @property {string} options.cytoBands.chrom Chromosome ID
	* @property {object} options.cytoColor Specifies colors for each cytogenetic band type
	* @property {boolean} options.zoom_enabled Allows zooming when true
	* @property {number} options.lineOpacity Set opacity of lines and regions across chromosomes marking phenotype positions
	* @property {number} options.backgroundColor Plot background color
  * @property {number} options.centromereIndentation Horizontal size of centromere indentation
  * @property {boolean} options.chromsOnly If true, no phenotype circles are shown on the plot
  * @property {boolean} options.drawCyto Only display cytobands when true
  * @property {number} options.minChromsRow Minimum number of chromosomes allowed in row
  * @property {number} options.maxChromsRow Maximum number of chromosomes allowed im row
	*/
	options: {
		height:  680,
		width:  732,
		xmax: 2928,
		ymax: 2720,
		chromheight: 1200,
		chromwidth: 64,
		circlesize: 32,
		margin: {top: 80, right: 0, bottom: 80, left: 20},
		circleradius: 16,
		linelength: 64,
		chromboxwidth: 244,
		chromboxheight: 1360,
		chromsperrow: 12,
		chromrows: 2,
		capheight: 40,
		chroms: [],
		phenos: [],
		cytoBands: [],
		cytoColor: {gneg: 'white', gpos25: '#ECECEC', gpos50: '#C7C7C7',
			gpos75: '9F9F9F', gpos100: '#787878', stalk: '#63B8FF',
			gvar: '#4F94CD', acen: '#0000A0'},
		posColorMap: {0: 'black', 1: 'blue', 2: 'red', 3: 'green', 4: 'orange', 5: 'purple',
			6: 'pink', 7: 'gold'},
		zoom_enabled: true,
		tip_enabled: true,
		include_map: false,
		zoom_map: false,
		lineOpacity: 0.5,
		backgroundColor: "white",
		centromereIndentation: 4,
		chromsOnly: false,
		drawCyto: true,
		minChromsRow: 6,
		maxChromsRow: 12,
		includedChroms: [],
		pixelscale: 4, // scale for pixels to xmax, ymax conversion
	},


	/**
    * Standard jQuery widget factory _create function.  Utilizes the options object for
    * parameters and creating the widget.
    * @private
  */
	_create: function() {
		// change the ymax based on the height/width ratio
		this.options.ymax = this.options.xmax * this.options.height/this.options.width;

		this.options.pixelscale = this.options.xmax / this.options.width;

		this.doefunc = function(val){return true;}

		// set up drawing canvas
		this.widgetID = this.element.attr('id');
		this.canvasID = 'csg-iphenogram-'+this.widgetID;
// 		this.highlighted_pheno = "";
		this.highlightedPhenos = {};
		this.highlightedShapes = {};
		this.highlightedFills = {};

		this.canvasjquery = $('<svg id="'+this.canvasID+'"></svg>').appendTo(this.element);
		this.canvas = d3.select('#'+this.canvasID)
			.attr("width", this.options.width)
			.attr("height", this.options.height)
			.attr("viewBox", "0 0 " + this.options.xmax + " " + this.options.ymax)
			.append("g");

		var zrect = this.canvas.append("rect")
    	.attr("class", "overlay")
    	.attr("width", this.options.xmax+2000)
	    .attr("height", this.options.ymax+2000)
	    .attr("pointer-events", "all")
	    .attr("x", -1000)
	    .attr("y", -1000)
	    .style("fill", this.options.backgroundColor)
		this._fillChroms();
		this._setIncludedChroms();
		this._placeChroms();
		var widget=this;


		if(this.options.tip_enabled){
			// Define 'div' for tooltips
			widget.tooltip = d3.select("body")
				.append("div")  // declare the tooltip div
				.attr("class", "tooltip")              // apply the 'tooltip' class
				.attr("id", "iPhenoTip")
				.style("position", "absolute")
  			.style("text-align", "center")
  			.style("padding", "10px")
  			.style("font", "0.8em sans-serif")
  			.style("background", "lightsteelblue")
				.style("border-radius", "8px")
				.style("border", "0px")
				.style("opacity", 0);                  // set the opacity to nil
		}

		this._addCytoBands();
		this._addPhenotypes();

    var w = this;
    
    w.zoom = d3.zoom().scaleExtent([1, 10]).on("zoom", function(event){w._zoom(w, event)});
    if(this.options.zoom_enabled){
// 			w.zm = d3.zoom().scaleExtent([1, 10]).on("zoom", function(){w._zoom(w)});
// 			this.canvas.call(w.zm);

// 			this.canvas.call(d3.zoom().scaleExtent([1, 10]).on("zoom", function(){
// 				w.canvas.attr("transform", d3.event.transform)
// 			}));

// 			this.canvas.call(d3.zoom().scaleExtent([1, 10]).on("zoom", function(){w._zoom(w)}));
			w.zoomMain=this.canvas.call(w.zoom);
		}

		if(this.options.drawCyto){
			this.drawCytoBands();
		}
		this._drawPhenos();

		this._drawChroms();
		if(this.options.include_map)
			this.maprectangle=this.canvas.append("rect")
  	  	.attr("class", "mapping")
    		.attr("width", this.options.xmax)
	    	.attr("height", this.options.ymax)
		    .attr("stroke", "black")
		    .attr("stroke-width", 8)
		    .attr("fill", "lightgray")
		    .attr("pointer-events", "all")
		    .attr("fill-opacity", 0.2);

	},

	/**
    * Standard jQuery widget factory _setOption function
    * @private
  */
	_setOption: function( key, value ) {
		var self = this;
		fnMap = {
			'chroms': function() {
				self._setIncludedChroms();
			},
      'includedChroms': function(){
				self._removePhenos(self);
				self._addPhenotypes();
				self._setIncludedChroms();
      	self._removeChroms();
      	self._placeChroms();
      	self.drawCytoBands();
      	self._drawPhenos();
      	self._drawChroms();
      },
      'chromsOnly': function(){
      	self._changePhenoVisibility();
      },
      'phenos': function(){
      	self._sortPhenos();
      }
    };
    this.options[ key ] = value;

		if (key in fnMap) {
	    fnMap[key]();
	  }
  },

	/**
    * Set which chromosomes are included on plot
    * @private
  */
  _setIncludedChroms: function(){
  	this.includedChroms = this.options.includedChroms;
  },

	/**
    * Removes all phenotypes
    * @private
  */
  _removePhenos: function(self){
  	self.canvas.selectAll(".phenocircle")
  		.remove();
  	self.canvas.selectAll(".lineLinker")
  		.remove();
  	self.canvas.selectAll(".chromLinker")
  		.remove();
  },

	/**
    * Removes all chromosomes
    * @private
  */
  _removeChroms:function(){
		this.chromd3.remove();
	},


	/**
    * Destroys this object
    * @private
  */
	_destroy: function() {
		return this._superApply(arguments);
	},

	/**
    * Returns overlay rectangle for use in mapping zoom and location.
    * @return {object} rectangle D3 rectangle
  */
	getMapRectangle: function(){
		return this.maprectangle;
	},

	/**
    * Set rectangle to use for mapping and reflecting zoom of plot.
    * @param {object} rectangle D3 rectangle
  */
	setZoomRectangle: function(rectangle){
		this.zoomrect = rectangle;
	},

	/**
    * Draw chromosomes on canvas
    * @private
  */
	_drawChroms: function(){
		var options = this.options;
		var widget = this;

		var chrompath = function(chromosome){
			var startchrom = chromosome.offset;
			var returnstr;
			if(chromosome.centromere_start != undefined){
				returnstr= "M0," + Math.round(widget.yconversion(chromosome.offset)) +" C0,"+ (Math.round(widget.yconversion(chromosome.offset))-widget.options.capheight).toString() + " " + widget.chromWidth + ","+ (Math.round(widget.yconversion(chromosome.offset))-widget.options.capheight).toString() + " " + widget.chromWidth + "," + Math.round(widget.yconversion(chromosome.offset)) +
			" V" + Math.round(widget.yconversion(chromosome.centromere_start+chromosome.offset)) +
			" L" + Math.round((widget.chromWidth-widget.chromWidth/widget.options.centromereIndentation)).toString() + "," + Math.round(widget.yconversion((chromosome.centromere_end-chromosome.centromere_start)/2+chromosome.centromere_start+chromosome.offset)) +
			" L" + widget.chromWidth + "," + Math.round(widget.yconversion(chromosome.centromere_end+chromosome.offset)) +
			" V" + Math.round(widget.yconversion(chromosome.size+chromosome.offset)) +
			" C" + widget.chromWidth + "," + (Math.round(widget.yconversion(chromosome.size+chromosome.offset))+widget.options.capheight).toString() + " 0," + (Math.round(widget.yconversion(chromosome.size+chromosome.offset))+widget.options.capheight).toString()  +
			" 0," + Math.round(widget.yconversion(chromosome.size+chromosome.offset)) +
			" V" + Math.round(widget.yconversion(chromosome.centromere_end+chromosome.offset)) +
			" L" + Math.round(widget.chromWidth/widget.options.centromereIndentation) + "," +   Math.round(widget.yconversion((chromosome.centromere_end-chromosome.centromere_start)/2+chromosome.centromere_start+chromosome.offset)) +
			" L0," + Math.round(widget.yconversion(chromosome.centromere_start+chromosome.offset)) +
			" V" + Math.round(widget.yconversion(chromosome.offset));
			}
			else{
				returnstr= "M0," + Math.round(widget.yconversion(chromosome.offset)) +" C0,"+ (Math.round(widget.yconversion(chromosome.offset))-widget.options.capheight).toString() +
					" " + widget.chromWidth + ","+ (Math.round(widget.yconversion(chromosome.offset))-widget.options.capheight).toString() +
					" " + widget.chromWidth + "," + Math.round(widget.yconversion(chromosome.offset)) +
					" V" + Math.round(widget.yconversion(chromosome.size+chromosome.offset)) +
					" C" + widget.chromWidth + "," + (Math.round(widget.yconversion(chromosome.size+chromosome.offset))+widget.options.capheight).toString() + " 0," + (Math.round(widget.yconversion(chromosome.size+chromosome.offset))+widget.options.capheight).toString() +
						" 0," + Math.round(widget.yconversion(chromosome.size+chromosome.offset)) +
					" V" + Math.round(widget.yconversion(chromosome.offset));
			}
			return returnstr;
		}


		var leftCent = function(chromosome){
			var points;
			points = "0," + Math.round(widget.yconversion(chromosome.centromere_start + chromosome.offset))
				+ " " +  widget.chromWidth/widget.options.centromereIndentation + "," +
				Math.round(widget.yconversion((chromosome.centromere_end-chromosome.centromere_start)/2+chromosome.centromere_start+chromosome.offset)) +
				" 0," + Math.round(widget.yconversion(chromosome.centromere_end + chromosome.offset));
			return points;
		}

		var rightCent = function(chromosome){
			var points;
			points = widget.chromWidth + "," + Math.round(widget.yconversion(chromosome.centromere_start + chromosome.offset)) +
			" " + Math.round((widget.chromWidth - widget.chromWidth/widget.options.centromereIndentation)).toString() +
			"," + Math.round(widget.yconversion((chromosome.centromere_end-chromosome.centromere_start)/2+chromosome.centromere_start+chromosome.offset)) +
			" " + widget.chromWidth + "," +  Math.round(widget.yconversion(chromosome.centromere_end + chromosome.offset));
			return points;
		}

		// when there is a centromere draw a white box over that location to clear any colored bands
		this.chromd3.select(function(d,i){return d.centromere_start != undefined ? this : null;})
			.append("polygon")
			.attr("points", function(d){return leftCent(d);})
			.attr("fill", widget.options.backgroundColor)
			.attr("stroke", widget.options.backgroundColor);

		this.chromd3.select(function(d,i){return d.centromere_start != undefined ? this : null;})
			.append("polygon")
			.attr("points", function(d){return rightCent(d);})
			.attr("fill", widget.options.backgroundColor)
			.attr("stroke", widget.options.backgroundColor);


		this.chromd3.append("path")
			.attr("d",function(d){ return chrompath(d);})
			.attr("class", "chromOutline")
			.attr("stroke-width",8)
			.attr("stroke","gray")
			.attr("fill","none");

		this.chromd3.append("text")
			.text(function(d){return d.name;})
			.attr("x", Math.round(widget.chromWidth/2))
			.attr("y", function(d){return Math.round(widget.yconversion(d.size+d.offset) + widget.options.capheight*2.5);})
			.attr("font-family", "sans-serif")
			.attr("font-size", "56")
			.attr("fill", "black")
			.attr("text-anchor", "middle");
	},

	/**
    * Place phenotypes using constrained force directed layout and draw phenotypes on plot.
    * @private
  */
	_drawPhenos: function(){
		var widget = this;
		for(i=0; i<widget.options.chroms.length; i++){

		var links = [], nodes =[];
		// set x and y for each node based on positions on the chromosome
		for(j=0; j<widget.options.chroms[i].positions.length;j++){

			widget.options.chroms[i].positions[j].count=0;
			widget.options.chroms[i].positions[j].x = widget.options.chroms[i].positions[j].fx = widget.chromWidth+widget.options.linelength;
			widget.options.chroms[i].positions[j].y = Math.round(this.yconversion(widget.options.chroms[i].positions[j].pos + widget.options.chroms[i].offset));
			widget.options.chroms[i].positions[j].orig_y = widget.options.chroms[i].positions[j].y;
			// add a node only if phenotypes associated with this position
			if(widget.options.chroms[i].positions[j].hasPhenos){
				nodes.push({x: widget.chromWidth, y:widget.options.chroms[i].positions[j].y, fx: widget.chromWidth, fy:widget.options.chroms[i].positions[j].y});
				nodes.push(widget.options.chroms[i].positions[j]);
				links.push({source: nodes[nodes.length-2], target: nodes[nodes.length-1]});
			}
		}			
			
		
		widget.options.chroms[i].layout = d3.forceSimulation()
											.nodes(nodes)
											.force("charge_force", d3.forceManyBody())
// 											.force("center_force", d3.forceCenter(widget,widget.options.chromheight/2))
											.force("links", d3.forceLink(links).strength(0.3))
											.on("tick", function(){
												for(var ni=0; ni<nodes.length; ni++){
													if(nodes[ni].y <  widget.options.circleradius/2){
 															nodes[ni].y = widget.options.circleradius/2;
 													}
 													else if(nodes[ni].y > widget.options.chromheight){
															nodes[ni].y = widget.options.chromheight - widget.options.circleradius/2;
													}
												}
											})
											.stop()
		
				
		for(var k=0; k<50; k++) widget.options.chroms[i].layout.tick();
		widget.options.chroms[i].layout.stop();
		var nodepositions=[];			
		for(var z=1; z<nodes.length; z+=2){
			nodepositions.push(nodes[z].y);
		}
		nodepositions.sort(function(a, b){return a - b;});
		// loop through nodepositions and set the positions to match in order for those that
		// have phenotypes
		var posIndx=0;
		for(z=0; z<nodepositions.length; z++){
			while(!widget.options.chroms[i].positions[posIndx].hasPhenos){
				posIndx++;
			}
				widget.options.chroms[i].positions[posIndx].y = Math.round(nodepositions[z]);
			posIndx++;
		}
		// update the x and y locations for the phenotypes
		// when multiple phenotypes at the same position, increment the x position by value
		// of the count variable in the position
		for(ph=0; ph<widget.options.chroms[i].phenos.length; ph++){
			if(widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].count >= widget.options.circlesPerChrom){
				widget.options.chroms[i].phenos[ph].y = widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].y +
					widget.options.circleradius;
				var posIndex = 2*widget.options.circlesPerChrom - widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].count-1;
				widget.options.chroms[i].phenos[ph].x = widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].x +
					posIndex * widget.options.circleradius*2;
			}
			else{
				widget.options.chroms[i].phenos[ph].y = widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].y
				widget.options.chroms[i].phenos[ph].x = widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].x +
					widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].count * widget.options.circleradius*2;
			}
			widget.options.chroms[i].posMap[widget.options.chroms[i].phenos[ph].position.toString()].count++;
		}

		}

		// add lines across chromosomes
		this._drawTransverseLines();

		var circleVis = 'hidden';
		if(!this.options.chromsOnly){
			circleVis = 'visible';
		}

		var circles = widget.chromd3.selectAll("circle")
				.data(function(d){return d.phenos;})
				.enter()
				.append("path");
		elementsize = (widget.options.circleradius*2) * (widget.options.circleradius*2);
		circles.attr("class","phenocircle")
			.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
			.attr("d", d3.symbol()
				.size(function(d){return elementsize;})
				.type(function(d){return widget.groupMap[d.group];}))
			.style("fill", function(d){
				return widget.fillMap[d.shading][d.pheno];})
			.style("stroke",  function(d){return widget.phenoColors(d.pheno);})
			.attr("visibility", circleVis)
			.on("click", function(event,d){
				var ph={};
				ph[d.pheno]=1;
				widget.highlightPhenos(ph, false);
				widget._triggerPhenoSelection(d.pheno);
			});				
			
		var regex = /[a-z0-9]/i;
		if(widget.options.tip_enabled){
			circles.on("mouseover", function(event,d){
				var message = "<ul class='no-bullet'><li>" + d.pheno + "</li>";
				if(regex.test(d.group))
					message += "<li>" + d.group + "</li>";
				if(d.pvalue > 0.0)
					message += "<li>pval: " + d.pvalue.toString() + "</li>";
				if(d.doe)
					message += "<li>beta: " + d.doe.toString() + "</li>";
				message += "</ul>";
				widget.tooltip.html(message)
					.style("left", (event.pageX) + "px")
					.style("top", (event.pageY + Math.round(widget.options.circlesize/widget.options.pixelscale)) + "px");
				widget.tooltip.transition()
					.duration(200)
					.style("opacity", 0.9);
			})
			.on("mouseout", function(d){
				widget.tooltip.transition()
					.duration(750)
					.style("opacity", 0.0);
			});
			}

	},

	/**
    * Toggle phenotype visibility
    * @private
  */
	_changePhenoVisibility: function(){
		var circleVis;
		this.options.chromsOnly ? circleVis = "hidden" : circleVis = "visible";
		this.chromd3.selectAll(".phenocircle")
			.attr("visibility", circleVis);
		this.chromd3.selectAll(".lineLinker")
			.attr("visibility", circleVis);
	},


	/**
    * Draw lines and regions on chromosomes for positions in plot.
    * @private
  */
	_drawTransverseLines: function(){
		var widget = this;
		this.options.chromsOnly ? lineVis = "hidden" : lineVis = "visible";

		var lineLength = widget.chromWidth +  widget.options.linelength - widget.options.circleradius;
		var linkLines = this.chromd3.selectAll("line")
			.data(function(d){
				var posWithPhenos = [];
				for(var i=0; i<d.positions.length; i++){
					if(d.positions[i].count > 0){
						posWithPhenos.push(d.positions[i]);
					}
				}
				return posWithPhenos;})
			.enter()
			.append("line")
			.attr("class", "lineLinker")
			.attr("x1", widget.chromWidth)
			.attr("y1", function(d){return d.orig_y;})
			.attr("x2", lineLength)
			.attr("y2", function(d){return d.y;})
			.attr("stroke", "gray")
			.attr("stroke-width", 4)
			.attr("visibility", lineVis);

		var transversePath = function(position, widget){
			var chromOffset = widget.chromMap[position.chrom].offset;
			var dpath = "M0," + position.orig_y + " H" + widget.chromWidth;
			// add lower portion when needed
			if(position.endPos != undefined){
				var ending_y = widget.yconversion(chromOffset + position.endPos);
				if( ending_y - position.orig_y >= 1.0){
					dpath += " V" + Math.round(widget.yconversion(chromOffset + position.endPos)) +
						" H0 V" + position.orig_y + " Z";
				}
			}
			return dpath;
		}

		// for drawing lines to phenotypes from positions on chromosomes
		var transverseLines = this.chromd3.selectAll(".chromLinker")
			.data(function(d){return d.positions;})
			.enter()
			.append("path")
			.attr("d", function(d){return transversePath(d, widget);})
			.attr("class", "chromLinker")
			.attr("stroke", function(d){
				var strokecolor = 'gray';
				if(d.endPos != undefined){
					d.posColor != undefined ? strokecolor=d.posColor: strokecolor='gray';
				}
				return strokecolor;})
			.attr("stroke-width", 4)
			.attr("fill", function(d){return d.posColor != undefined? d.posColor:"red";})
			.attr("opacity", widget.options.lineOpacity);

		if(this.options.tip_enabled){
			if(!this.options.chromsOnly){
				linkLines.on("mouseover", function(event, d){
					widget.tooltip.html("<ul class='no-bullet' style='text-align:left;'><li><span id='positionSpan' style='color:red;'>ID: " +
						d.id + "</span></li><li>pos: " + d.pos + "</li></ul>")
						.style("left", (event.pageX) + "px")
						.style("top", (event.pageY + Math.round(widget.options.circlesize/widget.options.pixelscale)) + "px")
						.style('opacity',0.8)
						.style('pointer-events','none');
				})
				.on("mouseout", function(event, d){
					widget.tooltip.transition()
						.duration(250)
						.delay(3000)
						.style("opacity", 0.0)
						.style("pointer-events", "none");
				});
			}

			transverseLines.on("mouseover", function(event,d){
					widget.tooltip.html("<ul class='no-bullet' style='text-align:left;'><li><span id='positionSpan' style='color:red;'>ID: " +
						d.id + "</span></li><li>pos: " + d.pos + "</li></ul>")
					.style("left", (event.pageX) + "px")
					.style("top", (event.pageY + Math.round(widget.options.circlesize/widget.options.pixelscale)) + "px")
					.style("opacity", 0.8)
					.style("pointer-events", "none");
			})
			.on("mouseout", function(d){
				widget.tooltip.transition()
					.duration(250)
					.delay(3000)
					.style("opacity", 0.0)
					.style("pointer-events", "none");
			});

		}

	},

	/**
    * Calculate and return maximum x coordinate
    * @private
    * @return {integer} Maximum x coordinate
  */
	_calcXmax:function(){
		return this.options.pixelscale * this.options.width;
	},

	/**
    * Calculate and return maximum y coordinate
    * @private
    * @return {integer} Maximum y coordinate
  */
	_calcYmax:function(){
		return this.options.pixelscale * this.options.height;
	},

	/**
    * Clear old chromosome information and add new
    * @private
  */
	_fillChroms:function(){

		this.chromMap = {};
		for(var i=0; i<this.options.chroms.length; i++){
			// map the chromosome name to the object
			this.chromMap[this.options.chroms[i].name.toString()]=this.options.chroms[i];
			this.options.chroms[i].links = [];
			this.options.chroms[i].phenos = [];
			this.options.chroms[i].positions = [];
			this.options.chroms[i].posMap = {};
			this.options.chroms[i].cytoBands = [];
			this.options.visible = false;
		}

	},

	/**
    * Place chromosomes on plot
    * @private
  */
	_placeChroms: function(){
		if(this.options.chromsOnly){
			this.options.minChromsRow = 8;
			this.options.maxChromsRow = 16;
		}

		var nChromsIncluded = this.options.includedChroms.length;
		var includedChromMap = {};
		if(nChromsIncluded ==  0){
			nChromsIncluded = this.options.chroms.length;
			for(var i=0; i<this.options.chroms.length; i++){
				this.options.chroms[i].visible = true;
				this.options.includedChroms.push(this.options.chroms[i].name);
				includedChromMap[this.options.includedChroms[i].toString()]=true;
			}
		}
		else{
			for(var i=0; i<this.options.includedChroms.length; i++){
				includedChromMap[this.options.includedChroms[i].toString()]=true;
			}
			for(var i=0; i<this.options.chroms.length; i++){
				if(includedChromMap[this.options.chroms[i].name.toString()]){
					this.options.chroms[i].visible = true;
				}
				else{
					this.options.chroms[i].visible = false;
				}
			}
		}
		this.options.chromsperrow=this._chromsPerRow(nChromsIncluded);
		this.options.chromrows = nChromsIncluded / this.options.chromsperrow;
		if(this.options.chromrows != Math.floor(this.options.chromrows)){
			this.options.chromrows++;
		}
		this.options.chromrows = Math.floor(this.options.chromrows);

		// calculate chromboxheight -- total height including margins
		this.options.chromboxheight = Math.floor(this.options.ymax / this.options.chromrows);
		this.options.chromheight = this.options.chromboxheight - this.options.margin.top - this.options.margin.bottom

		// calculate # of circles across possible for a chromosome in this set
		this.options.chromboxwidth = this.options.xmax / this.options.chromsperrow;
		var avail_x = this.options.chromboxwidth - this.options.margin.left - this.options.margin.right;

		// need to calculate circle size and then chromwidth
		this.plotCircleSize = Math.round(this.options.circlesize + (this.options.maxChromsRow/2 - Math.floor(this.options.chromsperrow/2))*2);
		this.options.circleradius = Math.round(this.plotCircleSize/2);
		this.chromWidth = Math.round((1 + (12 - this.options.chromsperrow) * 0.1) * this.options.chromwidth);

		this.options.circlesPerChrom = Math.floor((avail_x-this.chromWidth-this.options.linelength)/this.plotCircleSize) + 1;

		this.maxbp = d3.max(this.options.chroms.filter(function(d){
			return d.visible;
		}), function(d){ return d.size;});

		this.yconversion = d3.scaleLinear()
			.range([0,this.options.chromheight])
			.domain([0, this.maxbp]);
		for(var i=0; i<this.options.chroms.length; i++){
			this.options.chroms[i].offset = Math.round((this.maxbp - this.options.chroms[i].size)/2);
		}

		var widget=this;
		this.chromd3 = this.canvas.selectAll("g")
			.data(function(d){
				if(widget.options.chroms.length == widget.options.includedChroms.length){
					return widget.options.chroms;
				}
				else{
					var plottedChroms = [];
					for(var i=0; i<widget.options.chroms.length; i++){
						if(widget.options.chroms[i].visible){
							plottedChroms.push(widget.options.chroms[i]);
						}
					}
					return plottedChroms;
				}
			})
			.enter()
			.append("g")
			.attr("id", function(d,i){return "chr" + (i+1).toString();})
			.attr("class", "screenchrom")
			.attr("transform", function(d,i){
				var index = i;
				if( i > widget.options.chromsperrow -1){
					index = i - widget.options.chromsperrow;
				}
				return "translate(" + (widget.options.margin.left + widget.options.chromboxwidth * index)
					+ "," + (widget.options.margin.top + widget.options.chromboxheight * Math.floor(i/widget.options.chromsperrow)).toString() + ")";
				});

	},

	/**
    * Add cytobands to chromosomes
    * @private
  */
	_addCytoBands:function(){
		for(i=0; i<this.options.cytoBands.length; i++){
			if(this.chromMap[this.options.cytoBands[i].chrom] != undefined){
				this.chromMap[this.options.cytoBands[i].chrom.toString()].cytoBands.push(this.options.cytoBands[i]);
			}
		}
	},


	/**
    * Add phenotypes to chromosomes
    * @private
  */
	_addPhenotypes:function(){
		// clear map
		this.phenoMap={};
		this.groupMap={};
		this.fillMap={};
		symbols = d3.symbols;
		groupIndex=0;

		// need to clear chrom positions
		var chromKeys = Object.keys(this.chromMap);
		for(var i=0; i<chromKeys.length; i++){
			this.chromMap[chromKeys[i]].positions = [];
			this.chromMap[chromKeys[i]].posMap = {};
			this.chromMap[chromKeys[i]].phenos = [];
		}

		for(var i=0; i<this.options.phenos.length; i++){
			var pheno = this.options.phenos[i];
			var chr = this.chromMap[pheno.chrom];
			if(chr == undefined){
				continue;
			}
			if(pheno.position > chr.size){
				continue;
			}
			if(chr.posMap[pheno.position.toString()]==undefined){
					chr.positions.push({pos: pheno.position, endPos: pheno.endPosition, count: 0,
						id:pheno.id, chrom: pheno.chrom, posColor: pheno.posColor, hasPhenos: false});
				chr.posMap[pheno.position.toString()]= chr.positions[chr.positions.length-1];
			}
			if(pheno.pheno != undefined){
				chr.phenos.push(pheno);
				this.phenoMap[pheno.pheno]=1;
				chr.posMap[pheno.position.toString()].hasPhenos=true;
			}
			
			if(this.groupMap[pheno.group] == undefined){
				this.groupMap[pheno.group] = symbols[groupIndex];
				groupIndex++;
			}
			// fillMap is a map of maps with the inner key being the phenotype name
			if(this.fillMap[pheno.shading] == undefined){
				this.fillMap[pheno.shading] = {};
			}
			
		}
		this.phenoColors=this._getPhenoColors();
		this._setupFillMap();
	},

	/**
    * Sets phenotype colors based on number of phenotypes
    * @private
    * @return {object} d3 ordinal scale
  */
	_getPhenoColors:function(){
		var phenoNames = Object.keys(this.phenoMap);
		if(phenoNames.length <= 10){
			return d3.scaleOrdinal(d3.schemeCategory10).domain(phenoNames);
		}
		else if(phenoNames.length <= 20){
			// generated the following colors from https://medialab.github.io/iwanthue/
			// using default color setting at site
			return d3.scaleOrdinal(["#d64054","#59b94c","#ad56c1","#aab647","#6c68d2","#d9a235","#697ec3","#d7602d","#45aecf","#d34595","#60c395","#9f496b","#5b8436","#ca88ca","#36835e","#e3808a","#837838","#a85036","#db9d68","#9f6f25"])
				.domain(phenoNames);
		}
		else{
			// using default color setting at site
			return d3.scaleOrdinal(["#dd572e","#4c6be0","#59bf50","#a83ba8","#7ea736","#7a53bf","#bdb834","#b775e7","#459652","#e46fce","#65c793","#c84491","#4d732d","#d94476","#42bfcf","#cd3647","#328d73","#a84c24","#5572c0","#d89333","#679dd6","#97662f","#c393d9","#79702a","#8a5b9e","#b0ac60","#9a4662","#de9b6d","#dd83a1","#d36d60"])
				.domain(phenoNames)
		}
	},


	 /**
    * Specifies which circles (phenotypes) should be filled in with color.  All others
    * will be gray.
    * @param {object} phenopos The names of the phenotypes to display in color.
    * @param {boolean} reset When true, all phenotypes are colored.
    */
	highlightPhenos: function(phenopos, ireset){
		var widget = this;
		var reset = ireset;
		if(reset){
			widget.highlightedShapes={};
			widget.highlightedFills={};
		}
		if(Object.keys(phenopos).length==0)
			reset = true;
		widget.allPhenoHighlight=false;
		widget.allShapes=widget.allFills=false;
		if(Object.keys(widget.highlightedShapes).length==0 && Object.keys(widget.highlightedFills).length==0){
			widget.allPhenoHighlight=true;
		}
		if(Object.keys(widget.highlightedShapes).length == 0){
			widget.allShapes = true;
		}
		if(Object.keys(widget.highlightedFills).length == 0){
			widget.allFills = true;
		}

		elementsize = (widget.options.circleradius*2) * (widget.options.circleradius*2);
		if(reset){
			d3.selectAll(".phenocircle")
				.transition()
				.duration(250)
				.style("fill", function(d){
					if(widget.doefunc(d.doe) && (widget.highlightedShapes[d.group] || widget.highlightedFills[d.shading] || widget.allPhenoHighlight)){
						return widget.fillMap[d.shading][d.pheno];
					}
					else{
						return widget.fillMap[d.shading]['greyfill'];
					}
				})
				.style("stroke",function(d){
					if(widget.doefunc(d.doe) && (widget.highlightedShapes[d.group] || widget.highlightedFills[d.shading] || widget.allPhenoHighlight)){
						return widget.phenoColors(d.pheno);
					}
					else{
						return widget.fillMap[d.shading]['greyfill'];
					}
				})
			 	.attr("d", d3.symbol()
			 		.type(function(d){return widget.groupMap[d.group];})
			 		.size(elementsize));
			if(ireset){		
				widget.highlightedPhenos = {};
				widget.highlightedShapes = {};
				widget.highlightedFills = {};
			}
		}
		else{
			expandedsize = elementsize*4;		
			d3.selectAll(".phenocircle")
				.transition()
				.duration(500)
				.on("start", function(d){
					if(phenopos[d.pheno]){
						d3.select(this)
							.attr("d", d3.symbol()
								.type(function(d){return widget.groupMap[d.group];})
								.size(expandedsize))
					}
				})
				.style("fill", function(d){
					if(widget.doefunc(d.doe) && phenopos[d.pheno] && ( widget.allPhenoHighlight || ((widget.allShapes || (d.group in widget.highlightedShapes) ) && (widget.allFills || (d.shading in widget.highlightedFills)) ))){
						return widget.fillMap[d.shading][d.pheno];
					}
					else{
						return widget.fillMap[d.shading]['greyfill'];
					}
				})
				.style("stroke", function(d){
					if(widget.doefunc(d.doe) && phenopos[d.pheno] && ( widget.allPhenoHighlight || ((widget.allShapes || (d.group in widget.highlightedShapes) ) && (widget.allFills || (d.shading in widget.highlightedFills)) ))){
						return widget.phenoColors(d.pheno);
					}
					else{
						return 'lightgray';
					}
				})
			 	.on("end", function(){
			 		d3.select(this)
			 		 .transition()
			 		 .duration(1000)
			 		 .attr("d", d3.symbol()
			 		 	.type(function(d){return widget.groupMap[d.group];})
			 		 	.size(elementsize))

				});
				widget.highlightedPhenos = phenopos;
			}
		},
		
		
	 /**
    * Specifies which shapes (groups) should be filled in with color.  All others
    * will be gray.
    * @param {object} groupsSelected The names of the groups to display in color.
    * @param {boolean} reset When true, all shapes are colored that are included in the highlighted phenotypes
    */
	highlightShapes: function(groupsSelected, reset){
		var widget = this;
		if(Object.keys(groupsSelected).length==0){
			widget.highlightedShapes = groupsSelected;
			reset = true;
		}
		elementsize = (widget.options.circleradius*2) * (widget.options.circleradius*2);
		// when resetting shapes return the highlighted phenotypes to color
		if(reset){
			widget.highlightPhenos(widget.highlightedPhenos, false);
		}
		else{
			expandedsize = elementsize*4;
			widget.allPhenoHighlight=false;
			widget.allPhenos=widget.allFills=false;
			if(Object.keys(widget.highlightedPhenos).length==0 && Object.keys(widget.highlightedFills).length==0){
				widget.allPhenoHighlight=true;
			}
			if(Object.keys(widget.highlightedPhenos).length == 0){
				widget.allPhenos = true;
			}
			if(Object.keys(widget.highlightedFills).length == 0){
				widget.allFills = true;
			}
				
			d3.selectAll(".phenocircle")
				.transition()
				.duration(500)
				.on("start", function(d){
					if(groupsSelected[d.group]){
						d3.select(this)
							.attr("d", d3.symbol()
								.type(function(d){return widget.groupMap[d.group];})
								.size(expandedsize))
					}
				})
				.style("fill", function(d){
					if(widget.doefunc(d.doe) && groupsSelected[d.group] && (widget.allPhenoHighlight || ( (widget.allPhenos || widget.highlightedPhenos[d.pheno]) && (widget.allFills || widget.highlightedFills[d.shading]) ) )){
						return widget.fillMap[d.shading][d.pheno];
					}
					else{
						return widget.fillMap[d.shading]['greyfill'];
					}
				})
				.style("stroke", function(d){
					if(widget.doefunc(d.doe) && groupsSelected[d.group] && (widget.allPhenoHighlight || ( (widget.allPhenos || widget.highlightedPhenos[d.pheno]) && (widget.allFills || widget.highlightedFills[d.shading]) ) )){
						return widget.phenoColors(d.pheno);
					}
					else{
						return 'lightgray';
					}
				})
			 	.on("end", function(){
			 		d3.select(this)
			 		 .transition()
			 		 .duration(1000)
			 		 .attr("d", d3.symbol()
			 		 	.type(function(d){return widget.groupMap[d.group];})
			 		 	.size(elementsize))

				});
 				widget.highlightedShapes = groupsSelected;
			}
		},
		
		
	 /**
    * Specifies which shapes (groups) should be filled in with color.  All others
    * will be gray.
    * @param {object} fillsSelected The names of the fill types to display in color.
    * @param {boolean} reset When true, all shapes are colored that are included in the highlighted phenotypes
    */
	highlightFills: function(fillsSelected, reset){
		var widget = this;
		if(Object.keys(fillsSelected).length==0){
			widget.highlightedFills = fillsSelected;
			reset = true;
		}
		elementsize = (widget.options.circleradius*2) * (widget.options.circleradius*2);
		// when resetting shapes and fills return the highlighted phenotypes to color
		if(reset){
			widget.highlightPhenos(widget.highlightedPhenos, false);
		}
		else{
			expandedsize = elementsize*4;
			widget.allPhenoHighlight=false;
			widget.allShapes=widget.allPhenos=false;
			if(Object.keys(widget.highlightedPhenos).length==0 && Object.keys(widget.highlightedShapes).length==0){
				widget.allPhenoHighlight=true;
			}
			if(Object.keys(widget.highlightedPhenos).length == 0){
				widget.allPhenos = true;
			}
			if(Object.keys(widget.highlightedShapes).length == 0){
				widget.allShapes = true;
			}
				
			d3.selectAll(".phenocircle")
				.transition()
				.duration(500)
				.on("start", function(d){
					if(fillsSelected[d.shading]){
						d3.select(this)
							.attr("d", d3.symbol()
								.type(function(d){return widget.groupMap[d.group];})
								.size(expandedsize))
					}
				})
				.style("fill", function(d){
					if(widget.doefunc(d.doe) && fillsSelected[d.shading] && (widget.allPhenoHighlight || ( (widget.allPhenos || widget.highlightedPhenos[d.pheno]) && (widget.allShapes || widget.highlightedShapes[d.group]) ) )){
						return widget.fillMap[d.shading][d.pheno];
					}
					else{
						return widget.fillMap[d.shading]['greyfill'];
					}
				})
				.style("stroke", function(d){
					if(widget.doefunc(d.doe) && fillsSelected[d.shading] && (widget.allPhenoHighlight || ( (widget.allPhenos || widget.highlightedPhenos[d.pheno]) && (widget.allShapes || widget.highlightedShapes[d.group]) ) )){
						return widget.phenoColors(d.pheno);
					}
					else{
						return 'lightgray';
					}
				})
			 	.on("end", function(){
			 		d3.select(this)
			 		 .transition()
			 		 .duration(1000)
			 		 .attr("d", d3.symbol()
			 		 	.type(function(d){return widget.groupMap[d.group];})
			 		 	.size(elementsize))

				});
 				widget.highlightedFills = fillsSelected;
			}
		},		
		
		

		/**
    	* Zoom handler
	    * @private
    	* @param {object} w csg.iphenogram widget
		* @param {object} event 
  */
		_zoom: function(w, event){
			
			var zoom_event = event.transform;
			var scale = 1/zoom_event.k;
			
			w.canvas.attr("transform", zoom_event);
// 
// 			w.canvas.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
			if(w.options.zoom_map){
				x = zoom_event.x * -scale;
				y = zoom_event.y * -scale;
				w.zoomrect.attr("transform", "translate(" + x + "," + y + ")scale(" + scale + ")");
			}

		},

	/**
    * Resets plot to original zoom and position.
    */
	resetView: function(){
// 		this.canvas.attr("transform", "translate(0,0)scale(1)");
// 		this.zm.scale(1.0);
// 		this.zm.translate([0,0]);

		this.zoomMain.call(this.zoom.transform, d3.zoomIdentity);
		if(this.options.zoom_map){
			this.zoomrect.attr("transform", "translate(0,0)scale(1)");
		}
	},

	/**
    * Resets all colored circles (phenotypes) to original color.
    */
	resetPhenoColors: function(){
		this.highlightedPhenos={};
		this.highlightPhenos(this.highlightedPhenos, true);
	},

	/**
    * Parses string in PhenoGram format to fill data within widget.
    * @param {string} inputStr PhenoGram-formatted with input data on positions and phenotypes.
  */
	readPhenoInputString: function(inputStr){
		if(inputStr){
			var widget = this;
			/* create map for converting alternate names in input string */
			var colHeaderMap={};
			colHeaderMap['id']='snp';
			colHeaderMap['pheno']='phenotype';
			colHeaderMap['chrom']='chr';
			colHeaderMap['position']='pos';
			colHeaderMap['endpos']='end';
			colHeaderMap['pval']='pvalue';
			colHeaderMap['beta']='doe'; // direction of effect
			colHeaderMap['category']='phenotype';
			colHeaderMap['data']='shading';
			colHeaderMap['fill']='shading';
			colHeaderMap['tissue']='group';
			
			inputStr=this._checkHeaders(inputStr,['chr','pos','snp'],colHeaderMap);
		
			var inputData=d3.tsvParse(inputStr, function(d){
				var ePos, pCol, ph, doeval, shade, pval;
				d.end ? ePos = +d.end : ePos = undefined;
// 				d.poscolor ? pCol = d.poscolor.toString() : pCol = undefined;
				if(d.poscolor){
					pCol = widget._mapPosColors(d.poscolor.toString());
				}
				else{
					pCol = undefined;
				}
				d.phenotype ? ph = d.phenotype.toString() : ph = undefined;
				d.group ? gr = d.group.toString() : gr = ' ';
				d.pvalue ? pval = +d.pvalue : pval=0.0;
				// if no shading specified all will be default 
				d.shading ? shade = d.shading : shade = '-';
				d.doe ? doeval = +d.doe : doeval = undefined;
				return {
					chrom: d.chr.toString(),
					position: +d.pos,
					pheno: ph,
					id: d.snp.toString(),
					endPosition: ePos,
					posColor: pCol,
					group: gr,
					pvalue: pval,
					shading: shade,
					doe: doeval
				}});
			this._setOption('phenos', inputData);
		}

	},

	/**
    * Parses string in PhenoGram format to fill cytoband information within widget.
    * @param {string} inputStr PhenoGram formatted with position and cytoband type information.
  */
	readCytoBandString: function(inputStr){
		if(inputStr){
			var colHeaderMap={};
			colHeaderMap['giestain']='type';
			colHeaderMap['chromstart'] = 'startband';
			colHeaderMap['chromend'] = 'finishband';
			colHeaderMap['#chrom'] = 'chrom';
			inputStr=this._checkHeaders(inputStr,['chrom','type','startband','finishband'],colHeaderMap);
			var inputData=d3.tsvParse(inputStr, function(d){
				return {
					chrom: d.chrom.toString().replace('chr',''),
					type: d.type.toString(),
					startBand: +d.startband,
					finishBand: +d.finishband
			}});
			this._setOption('cytoBands', inputData);
		}
	},

	/**
    * Parses string in PhenoGram format to fill genome information within widget.
    * @param {string} inputStr PhenoGram formatted with size and centromere information.
  */
	readGenomeString: function(inputStr){
		if(inputStr){
			var colHeaderMap={};
			colHeaderMap['id']='name';
			inputStr=this._checkHeaders(inputStr,['name','size','centromerestart','centromereend'], colHeaderMap);
			var inputData=d3.tsvParse(inputStr, function(d){
				return {
					name: d.name.toString(),
					size: +d.size,
					centromere_start: +d.centromerestart,
					centromere_end: +d.centromereend
			}});
			this._setOption('chroms', inputData);
		}
	},


	/**
    * Check for valid headers in input string.  Convert all to lower case.
    * @private
    * @param {string} inString Input string
    * @param {object} requiredHeaders Required headers for this input
    * @param {object} colHeaderMap Optional alternate names for column header
    * @throws Will throw an error if not all required headers are present
  */
	_checkHeaders: function(inputStr,requiredHeaders, colHeaderMap){
		// get header values
		var nline=inputStr.indexOf("\n");
		var inHeaders={};
		var error=false;
		var errorString="";
		if(nline > -1){
			var headers=inputStr.substr(0,nline).split("\t");
			for(var i=0; i<headers.length; i++){
				headers[i]=headers[i].toLocaleLowerCase();
				if(colHeaderMap[headers[i]]){
					headers[i]=colHeaderMap[headers[i]];
				}
				inHeaders[headers[i]]=true;
			};
		 }
		 requiredHeaders.forEach(function(d){
		 	if(!inHeaders.hasOwnProperty(d)){
		 		error = true;
		 		errorString += " " + d;
		 	}
		 });
		 
		 if(error){
		 	throw "Missing headers: " + errorString;
		 }
		 
		 var nHeaderLine = headers.join("\t");
		 return inputStr.replace(inputStr.substr(0,nline),nHeaderLine);
	},
	
	/**
	 * Converts numbers used as positions into colors
	 */
	_mapPosColors: function(posString){
			/* snp colors from original phenogram */
			/* ['black','blue','red','green','orange','purple','pink', 'gold'] */
/*			colHeaderMap['poscolor']=['black','blue','red','green','orange','purple','pink', 'gold']; */		
		var pcolor;
		this.options.posColorMap.hasOwnProperty(posString) ? pcolor = this.options.posColorMap[posString] : pcolor=posString;
		return pcolor;
	},


	/**
    * Clears current plot and redraws with information held in widget.
  */
	drawPlot: function(){
		var self = this;
		self._removePhenos(this);
		self._removeChroms();
		self._fillChroms();
		self._placeChroms();
		self.removeCytoBands();
		self._addCytoBands();
		self._addPhenotypes();
		self.drawCytoBands();
		if(!self.options.drawCyto){
			self.hideCytoBands();
		}
		self._drawPhenos();
		self._drawChroms();
		if(this.zoomMain){
			this.zoomMain.call(this.zoom.transform, d3.zoomIdentity);
		}
	},

	/**
    * Draw cytobands on chromosomes on plot.
  */
	drawCytoBands: function(){
			var widget = this;
			widget.chromd3.selectAll(".cytoband")
				.data(function(d){return d.cytoBands;})
				.enter()
				.append("polygon")
				.attr("class", "cytoband")
				.attr("points", function(d){
						chr = widget.chromMap[d.chrom];
					 return "1, " + Math.round(widget.yconversion(d.startBand+chr.offset))  + " "
					+ (widget.chromWidth-1).toString() + "," + Math.round(widget.yconversion(d.startBand+chr.offset)) + " "
					+ (widget.chromWidth-1).toString() + "," + Math.round(widget.yconversion(d.finishBand+chr.offset)) + " "
					+ "1," + Math.round(widget.yconversion(d.finishBand+chr.offset));})
				.attr("stroke", "none")
				.attr("pointer-events", "all")
				.attr("stroke-width", 1)
				.attr("fill", function(d){ return widget.options.cytoColor[d.type];})
				.attr("fill-opacity", 0.7)
				.attr("stroke-opacity", 0);
		},


	/**
    * Removes cytobands from plot.
  */
	removeCytoBands: function(){
		this.chromd3.selectAll(".cytoband")
			.remove();
	},

	/**
    * Hides cytobands from chromosomes.  Does not remove them.
  */
	hideCytoBands: function(){
		this.chromd3.selectAll(".cytoband")
			.attr("visibility", "hidden");
		this.options.drawCyto=false;
	},

	/**
    * Makes hidden cytobands on chromosomes visible.
  */
	showCytoBands: function(){
		this.chromd3.selectAll(".cytoband")
			.attr("visibility", "visible");
		this.options.drawCyto=true;
	},

	/**
    * Returns information on phenotypes within plot.
    * @return {array} Array with phenotype name and phenotype color objects.
  */
	getPhenoTypes: function(){
		var phenoNames = Object.keys(this.phenoMap);
		var phenoReturn = [];
		for(var i=0; i<phenoNames.length;  i++){
			phenoReturn.push({name: phenoNames[i], color: this.phenoColors(phenoNames[i])});
		}
		return phenoReturn;
	},
	
	/**
	 * 
	 * @returns string with indeligble CSS characters replaced
	 */
	cssEscape: function(text){
		const regex = /\'/;
		return text.replace(regex, '_');
	},

	
	/**
	* Set up fill map based on phenotypes
	*/
	_setupFillMap: function(){
		var fillNames = Object.keys(this.fillMap);
		if(fillNames.length == 0)
			return;
		this.canvas.select("defs").remove();
		
// 		var svg = d3.select("body").append("svg").attr("id", "d3svg")
// 		    .attr("width", 120)
// 		    .attr("height", 120);
// 		var defs = svg.append("defs");
		var defs = this.canvas.append("defs");
		var phenoNames = Object.keys(this.phenoMap);

		// first pattern -- simple fill		
		for(var j=0; j<phenoNames.length; j++){
			// escape any problematic characters for css
			// phenoNames[j] = this.cssEscape(phenoNames[j])
			this.fillMap[fillNames[0]][phenoNames[j]] = this.phenoColors(phenoNames[j]);
		}
		this.fillMap[fillNames[0]]['greyfill'] = "#f1f1f1";
			
		
		// additional fill patterns need to be created and each color 
		// added
		if(fillNames.length > 1){
			for(var j=0; j<phenoNames.length; j++){ 				
				var patternName = fillNames[1] + "_" + this.cssEscape(phenoNames[j])
				defs.append("pattern")
					.attr("id", patternName)
					.attr("width", "10") 
					.attr("height","10")
					.attr("patternUnits","userSpaceOnUse")
					.attr("patternTransform","rotate(60)")
					.append("rect")
						.attr("width","6")
						.attr("height","8")
						.attr("transform","translate(0,0)")
						.attr("fill",this.phenoColors(phenoNames[j]));

				this.fillMap[fillNames[1]][phenoNames[j]] = "url(#" + patternName + ")";
			}
			// default for key
			var patternName = fillNames[1] + "_keyDefaultPattern";
			defs.append("pattern")
				.attr("id", patternName)
				.attr("width", "10")
				.attr("height", "10")
				.attr("patternUnits", "userSpaceOnUse")
				.attr("patternTransform","rotate(60)")
				.append("rect")
					.attr("width", "6")
					.attr("height", "8")
					.attr("transform","translate(0,0)")
					.attr("fill","#000000");
			
			this.fillMap[fillNames[1]]['defaultkey'] = "url(#" + patternName + ")";
			patternName = 'greyfill';
			defs.append("pattern")
				.attr("id", patternName)
				.attr("width", "10")
				.attr("height", "10")
				.attr("patternUnits", "userSpaceOnUse")
				.attr("patternTransform","rotate(60)")
				.append("rect")
					.attr("width", "6")
					.attr("height", "8")
					.attr("transform","translate(0,0)")
					.attr("fill","#f1f1f1");							

			this.fillMap[fillNames[1]]['greyfill'] = "url(#" + patternName + ")";	
		}	
		
		
		// additional fill patterns need to be created and each color 
		// added
		if(fillNames.length > 2){
			for(var j=0; j<phenoNames.length; j++){ 				
				var patternName = fillNames[2] + "_" + this.cssEscape(phenoNames[j]);
				

				defs.append("pattern")
					.attr("id", patternName)
					.attr("width", "10") 
					.attr("height","10")
					.attr("patternUnits","userSpaceOnUse")
					.attr("patternTransform","rotate(0)")
					.append("rect")
						.attr("width","4")
						.attr("height","8")
						.attr("transform","translate(0,0)")
						.attr("fill",this.phenoColors(phenoNames[j]));				
								
// 				defs.append("pattern")
// 				.attr({id:patternName, width:"10", height:"10",patternUnits:"userSpaceOnUse", patternTransform:"rotate(0)"})
// 				.append("rect")
// 					.attr({ width:"4", height:"8", transform:"translate(0,0)", fill:this.phenoColors(phenoNames[j]) });
				this.fillMap[fillNames[2]][phenoNames[j]] = "url(#" + patternName + ")";
			}
			// default for key
			var patternName = fillNames[2] + "_keyDefaultPattern";
			defs.append("pattern")
				.attr("id", patternName)
				.attr("width", "10") 
				.attr("height","10")
				.attr("patternUnits","userSpaceOnUse")
				.attr("patternTransform","rotate(0)")
				.append("rect")
					.attr("width","4")
					.attr("height","8")
					.attr("transform","translate(0,0)")
					.attr("fill", "#000000");				
					
// 			defs.append("pattern")
// 				.attr({id:patternName, width:"10", height:"10",patternUnits:"userSpaceOnUse", patternTransform:"rotate(0)"})
// 				.append("rect")
// 					.attr({ width:"4", height:"8", transform:"translate(0,0)", fill:"#000000"});
			this.fillMap[fillNames[2]]['defaultkey'] = "url(#" + patternName + ")";
			patternName = 'greyfill';
			defs.append("pattern")
				.attr("id", patternName)
				.attr("width", "10") 
				.attr("height","10")
				.attr("patternUnits","userSpaceOnUse")
				.attr("patternTransform","rotate(0)")
				.append("rect")
					.attr("width","4")
					.attr("height","8")
					.attr("transform","translate(0,0)")
					.attr("fill", "#f1f1f1");				
			
// 			defs.append("pattern")
// 				.attr({id:patternName, width:"10", height:"10",patternUnits:"userSpaceOnUse", patternTransform:"rotate(0)"})
// 				.append("rect")
// 				.attr({ width:"4", height:"8", transform:"translate(0,0)", fill:"#f1f1f1"});
			this.fillMap[fillNames[2]]['greyfill'] = "url(#" + patternName + ")";				
		}
		
		
					
	},
	
	/**
    * Returns information on groups within plot.
    * @return {array} Array with group name and symbol type
  */
	getGroups: function(){
		var groupNames = Object.keys(this.groupMap);
		var groupReturn = [];
		for(var i=0; i<groupNames.length;  i++){
			groupReturn.push({name: groupNames[i], symbolType: this.groupMap[groupNames[i]]});
		}
		return groupReturn;
	},


	/**
    * Returns information on fills within plot.
    * @return {array} Array with fill name and fill ID
  */
	getFills: function(){
		var fillNames = Object.keys(this.fillMap);
		// first is simple black for filled in with solid color
		var fillReturn = [{name: fillNames[0], fillType: "#000000"}];
		for(var i=1; i<fillNames.length;  i++){
			fillReturn.push({name: fillNames[i], fillType: this.fillMap[fillNames[i]]['defaultkey']});
		}
		
		return fillReturn;
	},


	/**
    * Returns chromosome names within plot.
    * @return {array} Array with chromosome names
  */
	getChromNames:function(){
		var names=[];
		for(var i=0; i<this.options.chroms.length; i++){
			names.push(this.options.chroms[i].name.toString());
		}
		return names;
	},

	/**
    * Calculate and set the number of chromosomes per row on plot.
    * @private
    * @param {integer} nChr Total number of chromosomes in plot.
  */
	_chromsPerRow: function(nChr){
		var minChr=this.options.minChromsRow;
		var maxChr=this.options.maxChromsRow;
		var fewestRows = 100000;
		var val=0;


		for(var i=minChr; i<=maxChr; i++){
			val = nChr/i;
			if(Math.ceil(val) < fewestRows){
				fewestRows = Math.ceil(val);
			}
		}
		var nChrPerRow=0;
		for(var i=maxChr; i>=minChr; i--){
			val = nChr/i;
			if(Math.abs(val-fewestRows) < 0.0001){
				nChrPerRow = i;
				break;
			}
			if(Math.ceil(val) > fewestRows){
				break;
			}
			nChrPerRow=i;

		}
		return nChrPerRow;
	},

  /**
    * Position selection event.
    *
    * @event csg.iphenogram#posselected
    * @type {object}
    * @property {object} posInfo - Position selected
    */
	_triggerPosition: function(position, eventType){
 		this._trigger('posselected', {type: eventType}, {posInfo: position});
	},

  /**
    * Position selection event.
    *
    * @event csg.iphenogram#selectpheno
    * @type {object}
    * @property {string} phenoName - Phenotype name selected
    */
	_triggerPhenoSelection: function(phenoname){
		this._trigger('selectpheno',{type: 'select'}, {phenoName: phenoname});
	},


	/**
    * Generate SVG for plot.  Includes a phenotype color key at bottom of plot if
    * phenotypes are included.
    * @param {float} resize Resizes the final plot based on the value passed.
    * @return {string} HTML containing SVG representation of plot with phenotype key
    * inserted.
  */
	getSVG: function(resize){
		var s = d3.select('#'+this.canvasID);
		var origWidth = s.attr("width");
		var origHeight = s.attr("height");

		var phenos = Object.keys(this.phenoMap);
		var groups = Object.keys(this.groupMap);
		var fills = Object.keys(this.fillMap);
		var adjustment = 1.0;
		if(phenos.length > 0){
			// place block around area to "white" out the non-visible portion and
			// leave room for the phenotype color key
			var phRowInfo = this._phenosPerRow();
			
			var grRowInfo;
			if(groups.length > 1){
				grRowInfo = this._groupsPerRow();
			}
			else{
				grRowInfo = {'groupRows':0, 'groupsPerRow':0};
			}
			var fillRowInfo;
			if(fills.length > 1){
				fillRowInfo = this._fillsPerRow();
			}
			else{
				fillRowInfo = {'fillRows':0, 'fillsPerRow':0};
			}
			adjustment = 1.0 + ((phRowInfo.phenoRows + grRowInfo.groupRows + fillRowInfo.fillRows) * this.options.ymax / 40) / this.options.ymax;
			
// 			if(groups.length > 1){
// 				var grRowInfo = this._groupsPerRow();
// 				adjustment = 1.0 + ((phRowInfo.phenoRows + grRowInfo.groupRows) * this.options.ymax / 40) / this.options.ymax;
// 			}
// 			else{
// 				var grRowInfo = {'groupRows':0, 'groupsPerRow':0};
// 				adjustment = 1.0 + (phRowInfo.phenoRows * this.options.ymax / 40) / this.options.ymax;
// 			}
		}

		s.attr("width", origWidth*resize)
			.attr("height", origHeight*resize * adjustment)
			.attr("version", 1.1)
			.attr("xmlns", "http://www.w3.org/2000/svg");

		s.attr("viewBox", "0 0 " + this.options.xmax + " " + this.options.ymax*adjustment);
		
		var zoom_info = d3.zoomTransform(this.canvas.node());
		// Hb {k: 8.316220983749089, x: -13369.738100584753, y: -7413.345818957543}

		// let xStart = (0-zoom_info.x);
		let xStart = 0;
		let yStart = (this.options.ymax - zoom_info.y)/zoom_info.k;
		let scale = 1/zoom_info.k;

		var boxG = this.canvas.append("g")
			.attr("class", "phenoKeyInfo")
			.attr("transform", "translate(" + xStart + "," + yStart + ")scale(1.0)");
			// .attr("transform", "translate(" + xStart + "," + yStart + ")scale(" + 1/zoom_info.k + ")");

		boxG.append("rect")
			.attr("class", "phenoKeyInfo")
			.attr("x",0)
			.attr("y",0)
			.attr("width", this.options.xmax)
			.attr("height", this.options.ymax*adjustment - this.options.ymax)
			.attr("fill",this.options.backgroundColor);


		if(phenos.length > 0){
			this._addPhenoKey(zoom_info.x,zoom_info.y,
				zoom_info.k,phRowInfo,grRowInfo,fillRowInfo);
		}

		// add the phenotype color key here (adjusted to place in the bottom of the final plot
			var html = s.node().parentNode.innerHTML;
			s.attr("width", origWidth)
				.attr("height", origHeight)
				.attr("viewBox", "0 0 " + this.options.xmax + " " + this.options.ymax);
		d3.selectAll(".phenoKeyInfo").remove();

		// console.log(html);
		return html;
	},

	/**
    * Calculate the number of phenotypes per row for SVG color key
    * @private
    * @return {object}  Total number of rows and phenotypes per row.
  */
	_phenosPerRow:function(){
		var nameLengths = jQuery.map(Object.keys(this.phenoMap), function(name,i){
			return +name.length;
		});
		var maxName = Math.max.apply(null,nameLengths);
		var constantSpacer = this.options.circleradius * 2;

		var phPerRow = Math.ceil(this.options.xmax / ((constantSpacer + maxName) * 14));

		var rows = Math.ceil(nameLengths.length/phPerRow);
		return {'phenoRows':rows, 'phenosPerRow':phPerRow};
	},
	
	/**
    * Calculate the number of groups per row for SVG shape key
    * @private
    * @return {object}  Total number of rows and groups per row.
  */
	_groupsPerRow:function(){
		var nameLengths = jQuery.map(Object.keys(this.groupMap), function(name,i){
			return +name.length;
		});
		var maxName = Math.max.apply(null,nameLengths);
		var constantSpacer = this.options.circleradius * 2;

		var phPerRow = Math.ceil(this.options.xmax / ((constantSpacer + maxName) * 14));

		var rows = Math.ceil(nameLengths.length/phPerRow);
		return {'groupRows':rows, 'groupsPerRow':phPerRow};
	},	
	
	
	/**
    * Calculate the number of fills per row for SVG pattern key
    * @private
    * @return {object}  Total number of rows and fills per row.
  */
	_fillsPerRow:function(){
		var nameLengths = jQuery.map(Object.keys(this.fillMap), function(name,i){
			return +name.length;
		});
		var maxName = Math.max.apply(null,nameLengths);
		var constantSpacer = this.options.circleradius * 2;

		var phPerRow = Math.ceil(this.options.xmax / ((constantSpacer + maxName) * 14));

		var rows = Math.ceil(nameLengths.length/phPerRow);
		return {'fillRows':rows, 'fillsPerRow':phPerRow};
	},		
	

	/**
    * Add phenotype color key
    * @private
    * @param {integer} xoffset Offset in x
    * @param {integer} yoffset Offset in y
    * @param {float} zoom Zoom level
    * @param {object} phenoRowInfo
    * @param {object} groupRowInfo
    * @param {object} fillRowInfo
  */
	_addPhenoKey:function(xoffset, yoffset, zoom, phenoRowInfo, groupRowInfo,fillRowInfo){
		var yStart = (this.options.ymax - yoffset)/zoom;
		var xStart = (0 - xoffset)/zoom;
		var phenoKeyG = this.canvas.append("g")
			.attr("class", "phenoKeyInfo")
			.attr("transform", "translate(" + xStart + "," + yStart + ")scale(" + 1/zoom + ")");

		var currPos = 0;
		var side = this.options.circleradius*2;
		var groupNames = Object.keys(this.groupMap);
		
		var widget = this;

		var yInterval = Math.round(this.options.ymax / 40);
		var xInterval = Math.round(this.options.xmax / (groupRowInfo.groupsPerRow+1));
		
		var yPos = yInterval-side;
		var xPos = side;
		var fontSize = 44;

		
		if(groupNames.length > 1){
			var elementsize = side * side;
			for(var i=0; i<groupNames.length; i++){
					
				shapeYPos = yPos - side/1.5
				shapeXPos = xPos + side/1.5
				
				phenoKeyG.append("path")
					.attr("class", "groupKeyInfo")
					.attr("transform", "translate(" + shapeXPos + "," + shapeYPos +")")
					.attr("d", d3.symbol()
						.size(elementsize)
						.type(widget.groupMap[groupNames[i]]))
					.style("fill", "white")
					.style("stroke", "black")
					.style("stroke-width", 4);
				phenoKeyG.append("text")
					.attr("class", "phenoKeyInfo")
					.attr("x", xPos+2*side)
					.attr("y", yPos)
					.text(groupNames[i])
					.attr("font-family", "sans-serif")
					.attr("font-size", fontSize)
					.attr("fill", "black")
					.attr("text-anchor", "start");	
				xPos += xInterval;
				currPos++;
				if(currPos > groupRowInfo.groupsPerRow){
					yPos += yInterval;
					xPos = side;
					currPos=0;
				}				
			}
			yPos += yInterval;
		}
		

		var fillNames = Object.keys(this.fillMap);		
		xInterval = Math.round(this.options.xmax / (fillRowInfo.fillsPerRow+1));
		xPos = side;
		currPos = 0;
		
		if(fillNames.length > 1){
			for(var i=0; i<fillNames.length; i++){
				phenoKeyG.append("rect")
					.attr("class", "fillKeyInfo")
					.attr("x",xPos)
					.attr("y",yPos-side)
					.attr("width", side*1.2)
					.attr("height", side*1.2)
					.attr("fill",widget.fillMap[fillNames[i]]['defaultkey']);
				phenoKeyG.append("text")
					.attr("class", "phenoKeyInfo")
					.attr("x", xPos+2*side)
					.attr("y", yPos)
					.text(fillNames[i])
					.attr("font-family", "sans-serif")
					.attr("font-size", fontSize)
					.attr("fill", "black")
					.attr("text-anchor", "start");

				xPos += xInterval;
				currPos++;
				if(currPos > fillRowInfo.fillsPerRow){
					yPos += yInterval;
					xPos = side;
					currPos=0;
				}
			}
			yPos += yInterval;
		}		

		var phenoNames = Object.keys(this.phenoMap);
		xInterval = Math.round(this.options.xmax / (phenoRowInfo.phenosPerRow+1));
		currPos=0;
		xPos = side;
		
		for(var i=0; i<phenoNames.length; i++){
			phenoKeyG.append("rect")
				.attr("class", "phenoKeyInfo")
				.attr("x",xPos)
				.attr("y",yPos-side)
				.attr("width", side)
				.attr("height", side)
				.attr("fill",this.phenoColors(phenoNames[i]));
			phenoKeyG.append("text")
				.attr("class", "phenoKeyInfo")
				.attr("x", xPos+2*side)
				.attr("y", yPos)
				.text(phenoNames[i])
				.attr("font-family", "sans-serif")
				.attr("font-size", fontSize)
				.attr("fill", "black")
				.attr("text-anchor", "start");

			xPos += xInterval;
			currPos++;
			if(currPos > phenoRowInfo.phenosPerRow){
				yPos += yInterval;
				xPos = side;
				currPos=0;
			}
		}
	},
	
	/**
	 * Sort phenotypes by chromosome and base pair position
	*/
	_sortPhenos:function(){
		this.options.phenos.sort(function(a,b){
			if(a.chrom < b.chrom){
				return -1;
			}
			if(a.chrom > b.chrom){
				return 1;
			}
			return a.position - b.position;
		});
	},


	/**
	* Sets the flag for plotting based on direction of effect
	* @param {integer} 0 - for negative, 1 - for positive and 2 - for all
	*/
	setDoeFlag:function(option){
	
	// set variable to a function that sets whether the value is right
		var widget = this;
		if(option == 0){
			widget.doefunc = function(val){return val < 0.0 ? true : false;}
		}
		else if(option == 1){
			widget.doefunc = function(val){return val >= 0.0 ? true : false;}		
		}
		else{ // option == 2
			widget.doefunc = function(val){return true;}			
		}
		this.highlightPhenos(this.highlightedPhenos, false);
	},

	/**
    * Provides current zoom and all other options for widget.
    * @return {object}  Zoom and all options.
  */
	getOptions: function () {
		const zoomTransform = d3.zoomTransform(this.canvas.node());
		return {zoomState: {
			k: zoomTransform.k,
			x: zoomTransform.x,
			y: zoomTransform.y},
			options: this.options};
	},

	/**
    * Loads in options and zoom to draw plot.
    * @param {object} oldOptions Contains zoom and options as provided by getoptions method.
  */
	resetOptions: function (savedOptions) {
		// Set widget options (excluding _zoomState)
		// const { _zoomState, ...cleanOptions } = savedOptions;
		// this._setOptions(savedOptions);
		this.options = savedOptions.options;
		this.drawPlot();
		_zoomState = savedOptions.zoomState;
	
		// Restore zoom if zoom is enabled and we have saved zoom info
		if (_zoomState && this.options.zoom_enabled) {
		  const t = d3.zoomIdentity
			.translate(_zoomState.x, _zoomState.y)
			.scale(_zoomState.k);
	
		  this.zoomMain.call(this.zoom.transform, t);
	
		  if (this.options.zoom_map && this.zoomrect) {
			const scale = 1 / _zoomState.k;
			const x = _zoomState.x * -scale;
			const y = _zoomState.y * -scale;
			this.zoomrect.attr("transform", `translate(${x},${y})scale(${scale})`);
		  }
		}
	  }


});


