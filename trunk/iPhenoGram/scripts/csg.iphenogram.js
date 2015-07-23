// csg.iphenogram.js

$.widget('csg.iphenogram',{

	options: {
		height:  680, //overall height in pixels
		width:  732, //overall width in pixels
		pixelscale: 4, // scale for pixels to xmax, ymax conversion
		xmax: 2928, //scaled width
		ymax: 2720, // scaled height
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
		zoom_enabled: true,
		tip_enabled: true,
		include_map: false,
		zoom_map: false,
		lineOpacity: 0.5,
		plotname: "iphenogram",
		backgroundColor: "white",
		centromereIndentation: 4,
		chromsOnly: false,
		drawCyto: true,
		minChromsRow: 6,
		maxChromsRow: 12,
		includedChroms: [],
		background: 'rgb(255,255,255)'
	},
	
// 	_init: function(){
// 	},
	
	_create: function() {

		this._setValues();
		
		// change the ymax based on the height/width ratio
		this.options.ymax = this.options.xmax * this.options.height/this.options.width;
		
		this.options.pixelscale = this.options.xmax / this.options.width; 
		
		// set up drawing canvas
		this.widgetID = this.element.attr('id');
		this.canvasID = 'csg-iphenogram-'+this.widgetID;
		this.highlighted_pheno = "";
	  
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
	    .attr("x", -1000)
	    .attr("y", -1000)
	    .style("fill", this.options.background)

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
				.style("opacity", 0);                  // set the opacity to nil
		}

		this._addCytoBands();		
		this._addPhenotypes();

    var w = this;
    if(this.options.zoom_enabled){
			w.zm = d3.behavior.zoom().scaleExtent([1, 10]).on("zoom", function(){w._zoom(w)});
			this.canvas.call(w.zm);
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
		    .attr("fill-opacity", 0.2);
		
	},
	
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
      }
    };
    this.options[ key ] = value;

		if (key in fnMap) {
	    fnMap[key]();
	  }
  },
  
  _setIncludedChroms: function(){
  	this.includedChroms = this.options.includedChroms;
//   	this.options.includedChroms=[];
  },
  
  _removePhenos: function(self){
  	self.canvas.selectAll(".phenocircle")
  		.remove();
  	self.canvas.selectAll(".lineLinker")
  		.remove();
  	self.canvas.selectAll(".chromLinker")
  		.remove();
  },
  
  _removeChroms:function(){
		// this._removePhenos(this);
		this.chromd3.remove();
	},
	
	_destroy: function() {
		return this._superApply(arguments);
	},
	
	
	getMapRectangle: function(){
		return this.maprectangle;
	},
	
	setZoomRectangle: function(rectangle){
		this.zoomrect = rectangle;
	},
	
	// draw chromosomes on canvas
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
	
	
	_drawPhenos: function(){
		var widget = this;
		for(i=0; i<widget.options.chroms.length; i++){
		
		var links = [], nodes =[];
		// set x and y for each node based on positions on the chromosome
		for(j=0; j<widget.options.chroms[i].positions.length;j++){
			
			widget.options.chroms[i].positions[j].count=0;
			widget.options.chroms[i].positions[j].x = widget.chromWidth+widget.options.linelength;
			widget.options.chroms[i].positions[j].y = Math.round(this.yconversion(widget.options.chroms[i].positions[j].pos + widget.options.chroms[i].offset));
			widget.options.chroms[i].positions[j].orig_y = widget.options.chroms[i].positions[j].y;
			// add a node only if phenotypes associated with this position
			if(widget.options.chroms[i].positions[j].hasPhenos){
				nodes.push({x: widget.chromWidth, y:widget.options.chroms[i].positions[j].y, fixed: true})
				nodes.push(widget.options.chroms[i].positions[j]);
				links.push({source: nodes[nodes.length-2], target: nodes[nodes.length-1]});
			}
		}	
			widget.options.chroms[i].layout = d3.layout.force()
											.gravity(0.01)
											.friction(0.7)
											.linkDistance(widget.options.linelength)
											.linkStrength(1.0)
											.size([widget.options.linelength, widget.options.chromheight])
											.charge(-85)
											.nodes(nodes)
											.links(links)
											.on("tick", function(){
												for(var ni=0; ni<nodes.length; ni++){
													if(nodes[ni].x != widget.chromWidth){
														nodes[ni].x = widget.chromWidth+widget.options.linelength;
													}
													if(nodes[ni].y < widget.options.circleradius/2){
														nodes[ni].y = widget.options.circleradius/2;
													}
													else if(nodes[ni].y > widget.options.chromheight){
														nodes[ni].y = widget.options.chromheight - widget.options.circleradius/2;
													}
												}
												}
											)
											.start();
			for(var k=0; k<14; k++) widget.options.chroms[i].layout.tick();
			widget.options.chroms[i].layout.stop();
			var nodepositions=[];
		for(var z=1; z<nodes.length; z+=2){
			nodepositions.push(nodes[z].y);	
		}
		nodepositions.sort(function(a, b){return a - b;});
		// loop through nodepositions and set the positions to match in order for those that
		// have phenotypes
// 		for(z=0; z<widget.options.chroms[i].positions.length; z++){
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
				.append("circle");
	
		circles.attr("class", "phenocircle")
			.attr("cx", function(d){return d.x;})
			.attr("cy", function(d){return d.y;})
			.attr("r", widget.options.circleradius)
			.attr("fill", function(d){return widget.phenoColors(d.pheno);})
			.attr("stroke", function(d){return widget.phenoColors(d.pheno);})
			.attr("visibility", circleVis)
			.on("click", function(d){
				var ph={};
				ph[d.pheno]=1;
				widget.highlightPhenos(ph, false);
				widget._triggerPhenoSelection(d.pheno);
			});
	
		if(widget.options.tip_enabled){
			circles.on("mouseover", function(d){
				widget.tooltip.html("<ul class='no-bullet'><li>" + d.pheno + "</li></ul>")
					.style("left", (d3.event.pageX) + "px")
					.style("top", (d3.event.pageY + Math.round(widget.options.circlesize/widget.options.pixelscale)) + "px");
				widget.tooltip.transition()
					.duration(200)
					.style("opacity", 0.9);
			})
			.on("mouseout", function(d){
				widget.tooltip.transition()
					.duration(750)
					.style({opacity: 0, 'pointer-events': 'none' });
			});
			}
		
	},
	
	_changePhenoVisibility: function(){
		var circleVis;
		this.options.chromsOnly ? circleVis = "hidden" : circleVis = "visible";
		this.chromd3.selectAll(".phenocircle")
			.attr("visibility", circleVis);
		this.chromd3.selectAll(".lineLinker")
			.attr("visibility", circleVis);
	},
	
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
				dpath += " V" + Math.round(widget.yconversion(chromOffset + position.endPos)) +
					" H0 V" + position.orig_y + " Z";
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
			.attr("stroke", function(d){return d.endPos != undefined ? "none":"gray";})
			.attr("stroke-width", 4)
			.attr("fill", function(d){return d.posColor != undefined? d.posColor:"red";})
			.attr("opacity", widget.options.lineOpacity);

		if(this.options.tip_enabled){
			if(!this.options.chromsOnly){
				linkLines.on("mouseover", function(d){
					widget.tooltip.html("<ul class='no-bullet' style='text-align:left;'><li><span id='positionSpan' style='color:red;cursor:pointer;text-decoration:underline;'>ID: " + 
						d.id + "</span></li><li>pos: " + d.pos + "</li></ul>")
						.style("left", (d3.event.pageX) + "px")
						.style("top", (d3.event.pageY + Math.round(widget.options.circlesize/widget.options.pixelscale)) + "px")
						.style({'opacity': 0.8, 'pointer-events': 'all'});
					d3.select("#positionSpan").on('click', function(e){widget._triggerPosition(d, 'posID');});
				})
				.on("mouseout", function(d){
					widget.tooltip.transition()
						.duration(250)
						.delay(3000)
						.style({opacity: 0, 'pointer-events': 'none' });
				});
			}
			
			transverseLines.on("mouseover", function(d){
					widget.tooltip.html("<ul class='no-bullet' style='text-align:left;'><li><span id='positionSpan' style='color:red;cursor:pointer;text-decoration:underline;'>ID: " + 
						d.id + "</span></li><li>pos: " + d.pos + "</li></ul>")
					.style("left", (d3.event.pageX) + "px")
					.style("top", (d3.event.pageY + Math.round(widget.options.circlesize/widget.options.pixelscale)) + "px")
					.style({'opacity': 0.8, 'pointer-events': 'all'});
				d3.select("#positionSpan").on('click', function(e){widget._triggerPosition(d, 'posID');});
			})
			.on("mouseout", function(d){
				widget.tooltip.transition()
					.duration(250)
					.delay(3000)
					.style({opacity: 0, 'pointer-events': 'none' });
			});	
			
		}
		
	},
	
	
	// sets values based on options
	_setValues: function(){
// 		this.options.xmax = this._calcXmax();
// 		this.options.ymax = this._calcYmax();
	},
	
	_calcXmax:function(){
		return this.options.pixelscale * this.options.width;
	},
	
	_calcYmax:function(){
		return this.options.pixelscale * this.options.height;
	},

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
		// size is 32 for 12 chroms across and then increases from there 
		this.plotCircleSize = Math.round(this.options.circlesize + (this.options.maxChromsRow/2 - Math.floor(this.options.chromsperrow/2))*2);
		this.options.circleradius = Math.round(this.plotCircleSize/2);
		this.chromWidth = Math.round((1 + (12 - this.options.chromsperrow) * 0.1) * this.options.chromwidth);
		
		this.options.circlesPerChrom = Math.floor((avail_x-this.chromWidth-this.options.linelength)/this.plotCircleSize) + 1;
		
		this.maxbp = d3.max(this.options.chroms.filter(function(d){
			return d.visible;
		}), function(d){ return d.size;});
					
		this.yconversion = d3.scale.linear()
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
	
	_addCytoBands:function(){
		for(i=0; i<this.options.cytoBands.length; i++){ 
			if(this.chromMap[this.options.cytoBands[i].chrom] != undefined){
				this.chromMap[this.options.cytoBands[i].chrom.toString()].cytoBands.push(this.options.cytoBands[i]);
			}
		}
	},
	
	_addPhenotypes:function(){
		// clear map
		this.phenoMap={};
		
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
		}
		this.phenoColors=this._getPhenoColors()
	},
	
	_getPhenoColors:function(){
		var phenoNames = Object.keys(this.phenoMap);
		if(phenoNames.length <= 10){
			return d3.scale.category10().domain(phenoNames);
		}
		else{
			return d3.scale.category20().domain(phenoNames);
		}
	},
	
	
	highlightPhenos: function(phenopos, reset){
		var widget = this;
		if(Object.keys(phenopos).length==0)
			reset = true;
		if(reset){
	
			d3.selectAll(".phenocircle")
				.transition()
				.duration(250)
				.attr("fill", function(d){return widget.phenoColors(d.pheno);})
				.attr("stroke",function(d){return widget.phenoColors(d.pheno);})
				.attr("r",widget.options.circleradius);
			this.highlighted_pheno="";
		}
		else{
			
			d3.selectAll(".phenocircle")
				.transition()
				.duration(500)
				.each("start", function(d){
					if(phenopos[d.pheno]){
						d3.select(this)
							.attr("r", widget.options.circleradius*2);
					}
				})
				.attr("fill", function(d){
					if(phenopos[d.pheno]){
						return widget.phenoColors(d.pheno);
					}
					else{
						return "#f1f1f1";
					}
				})
				.attr("stroke", function(d){
					if(phenopos[d.pheno]){
						return widget.phenoColors(d.pheno);
					}
					else{
						return 'lightgray';
					}
				})
			 	.each("end", function(){
			 		d3.select(this)
			 		 .transition()
			 		 .duration(1000)
				 	 .attr("r", widget.options.circleradius); 

				});
				widget.highlighted_pheno=phenopos.pheno;
			}
		},
		
		_zoom: function(w){		
			var coords = d3.event.translate;
			var scale = 1/d3.event.scale;
		
			w.canvas.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
			if(w.options.zoom_map){
				coords[0] = coords[0] * -scale;
				coords[1] = coords[1] * -scale;
				w.zoomrect.attr("transform", "translate(" + coords[0] + "," + coords[1]+ ")scale(" + scale + ")");
			}

		},
		
	resetView: function(){
		this.canvas.attr("transform", "translate(0,0)scale(1)");
		this.zm.scale(1.0);
		this.zm.translate([0,0]);
		if(this.options.zoom_map){
			this.zoomrect.attr("transform", "translate(0,0)scale(1)");
		}
	},	
	
	resetPhenoColors: function(){
		this.highlightPhenos({}, true);
	},
		
	readPhenoInputString: function(inputStr){
		if(inputStr){
			var inputData=d3.tsv.parse(inputStr, function(d){
				var ePos, pCol, ph;
				d.endPosition ? ePos = +d.endPosition : ePos = undefined
				d.posColor ? pCol = d.posColor.toString() : pCol = undefined
				d.pheno ? ph = d.pheno.toString() : ph = undefined
				return {
					chrom: d.chrom.toString(),
					position: +d.position,
					pheno: ph,
					id: d.id.toString(),
					endPosition: ePos,
					posColor: pCol
				}});
			this._setOption('phenos', inputData);
		}
	},
	
	readCytoBandString: function(inputStr){
		if(inputStr){
			var inputData=d3.tsv.parse(inputStr, function(d){
				return {
					chrom: d.chrom.toString(),
					type: d.type.toString(),
					startBand: +d.startBand,
					finishBand: +d.finishBand
			}});
			this._setOption('cytoBands', inputData);
		}
	},
	
	readGenomeString: function(inputStr){
		if(inputStr){
			var inputData=d3.tsv.parse(inputStr, function(d){
				return {
					name: d.name.toString(),
					size: +d.size,
					centromere_start: +d.centromereStart,
					centromere_end: +d.centromereEnd
			}});
			this._setOption('chroms', inputData);
		}
	},
	

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
	},

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
				.attr("stroke-width", 1)
				.attr("fill", function(d){ return widget.options.cytoColor[d.type];})
				.attr("fill-opacity", 0.7)
				.attr("stroke-opacity", 0);
		},

	
	removeCytoBands: function(){
		this.chromd3.selectAll(".cytoband")
			.remove();
	},
	
	hideCytoBands: function(){
		this.chromd3.selectAll(".cytoband")
			.attr("visibility", "hidden");
		this.options.drawCyto=false;
	},
	
	showCytoBands: function(){
		this.chromd3.selectAll(".cytoband")
			.attr("visibility", "visible");
		this.options.drawCyto=true;
	},
	
	getPhenoTypes: function(){
		var phenoNames = Object.keys(this.phenoMap);
		var phenoReturn = [];
		for(var i=0; i<phenoNames.length;  i++){
			phenoReturn.push({name: phenoNames[i], color: this.phenoColors(phenoNames[i])});
		}
		return phenoReturn;
	},
	
	getChromNames:function(){
		var names=[];
		for(var i=0; i<this.options.chroms.length; i++){
			names.push(this.options.chroms[i].name.toString());
		}
		return names;
	},
	
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
	
	_triggerPosition: function(position, eventType){
 		this._trigger('simple', {type: eventType}, {posInfo: position});
	},
	
	_triggerPhenoSelection: function(phenoname){
		this._trigger('select',{type: 'select'}, {phenoName: phenoname});
		
	},
	
	
	getSVG: function(resize){
		var s = d3.select('#'+this.canvasID);
		var origWidth = s.attr("width");
		var origHeight = s.attr("height");

		var phenos = Object.keys(this.phenoMap);
		var adjustment = 1.0;
		if(phenos.length > 0){
			// place block around area to "white" out the non-visible portion and
			// leave room for the phenotype color key
			var phRowInfo = this._phenosPerRow();		
			adjustment = 1.0 + (phRowInfo.phenoRows * this.options.ymax / 40) / this.options.ymax;
		}

		s.attr("width", origWidth*resize)
			.attr("height", origHeight*resize * adjustment)
			.attr("version", 1.1)
			.attr("xmlns", "http://www.w3.org/2000/svg");
		
		s.attr("viewBox", "0 0 " + this.options.xmax + " " + this.options.ymax*adjustment);
// 		var yPos = this.options.ymax - this.zm.translate()[1] / this.zm.scale();
		var yPos = (this.options.ymax - this.zm.translate()[1])/this.zm.scale();
		// need to cut it at 2820 (taking into account scale effect)
		this.canvas.append("rect")
			.attr("class", "phenoKeyInfo")
			.attr("x",0)
			.attr("y",yPos)
			.attr("width", this.options.xmax)
			.attr("height", this.options.ymax*adjustment - this.options.ymax)
			.attr("fill",this.options.background);
		if(phenos.length > 0){
			this._addPhenoKey(this.zm.translate()[0],this.zm.translate()[1],
				this.zm.scale(),phRowInfo);
		}
		
		// add the phenotype color key here (adjusted to place in the bottom of the final plot
			var html = s.node().parentNode.innerHTML;
			s.attr("width", origWidth)
				.attr("height", origHeight)
				.attr("viewBox", "0 0 " + this.options.xmax + " " + this.options.ymax);
		d3.selectAll(".phenoKeyInfo").remove();
		return html;
	},
	
	// calculate the number of phenotypes 
	_phenosPerRow:function(){
		var nameLengths = jQuery.map(Object.keys(this.phenoMap), function(name,i){
			return +name.length;
		});
		var maxName = Math.max.apply(null,nameLengths);
		// color square worth 4 letters? 
		// add maxName and an additional space (4 letters?)
		var constantSpacer = this.options.circleradius * 2;

		var phPerRow = Math.ceil(this.options.xmax / ((constantSpacer + maxName) * 14));

		var rows = Math.ceil(nameLengths.length/phPerRow);
		return {'phenoRows':rows, 'phenosPerRow':phPerRow};
	},
	
	_addPhenoKey:function(xoffset, yoffset, zoom, phenoRowInfo){
		var yStart = (this.options.ymax - yoffset)/zoom;
		var xStart = (0 - xoffset)/zoom;
		var phenoKeyG = this.canvas.append("g")
			.attr("class", "phenoKeyInfo")
			.attr("transform", "translate(" + xStart + "," + yStart + ")scale(" + 1/zoom + ")");

		var phenoNames = Object.keys(this.phenoMap);
		var widget = this;
		var currPos = 0;
		var yInterval = Math.round(this.options.ymax / 40);
		var xInterval = Math.round(this.options.xmax / (phenoRowInfo.phenosPerRow+1));
		var side = this.options.circleradius*2;
		var yPos = yInterval-side;
		var xPos = side;
		var fontSize = 56;
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
	
	getOptions:function(){
		if(this.zm)
		 return {scale: this.zm.scale(), translate: this.zm.translate(), options: this.options};
		else
			return  {scale: 1.0, translate: [0,0], options: this.options};
	},
	
	resetOptions:function(oldOptions){
		// set options
		this.options = oldOptions.options;
		this.drawPlot();
	
		if(this.zm){
			var transformStr = "translate("+ oldOptions.translate[0] + "," + oldOptions.translate[1]+ ")scale(" + oldOptions.scale + ")";
			this.canvas.attr("transform", transformStr);
			this.zm.scale(oldOptions.scale);
			this.zm.translate(oldOptions.translate);
			if(this.options.zoom_map){
				var coords = [];
				var scale = 1/oldOptions.scale;
				coords[0] = oldOptions.translate[0] * -scale;
				coords[1] = oldOptions.translate[1] * -scale;	
// 				this.zoomrect.attr("transform", transformStr);
				this.zoomrect.attr("transform", "translate(" + coords[0] + "," + coords[1]+ ")scale(" + scale + ")");
			}
		}
	}
	
});




