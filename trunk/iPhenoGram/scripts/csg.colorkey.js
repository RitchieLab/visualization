//colorkey-widget.js

$.widget('csg.colorkey',{

	options: {
		height:  500, //overall height in pixels
		width:  150, //overall width in pixels
		barHeight: 18,
		keyInfo: []
	},
	
	_create: function() {
	
		this.widgetID = this.element.attr('id');
		this.canvasID = 'csg-colorkey-'+this.widgetID;
		
		this.canvasjquery = $('<svg id="'+this.canvasID+'"></svg>').appendTo(this.element);
		this.canvas = d3.select('#'+this.canvasID)
			.attr("width", this.options.width)
			.attr("height", this.options.height)
			.append("g");
		
		this._drawKey(this.options.keyInfo);
// 		this.selectedKeys = {};
		
		var whiteCheck = this.canvas.append("rect")
			.style("fill", "white");
		this.whiteColorString = whiteCheck.style("fill");
		whiteCheck.remove();
		this.highlightedColor = "#EDEDED";
	},
	
	// draw key with color box to left, followed by text
	// info is array of objects with name and color 
	_drawKey: function(info){
		this.selectedKeys = {};
		var n = info.length;
		d3.select('#'+this.canvasID)
			.attr("height", this.options.barHeight * n);
		var barHeight = this.options.barHeight;
		var barWidth = this.options.width * 0.2;
		var self = this;
		
		var bar = this.canvas.selectAll("g")
			.data(info)
			.enter().append("g")
			.attr("class", 'phenoColorBars')
			.attr("transform", function(d, i){ return "translate(0," + i * barHeight + ")";});
		
		self.bgBars = bar.append("rect")
			.attr("width", this.options.width)
			.attr("height", barHeight)
			.style("fill", "rgb(255, 255, 255)")
			.on("click", function(d){
				var thisBar = d3.select(this);
				var fillColor = thisBar.style("fill");
// 				if(fillColor == self.whiteColorString){// "rgb(255, 255, 255)" || fillColor == "#ffffff"){
				if(!self.selectedKeys[d.name]){
					newColor = self.highlightedColor;
					self.selectedKeys[d.name]=1;
				}
				else{
					newColor = "rgb(255, 255, 255)";
					delete self.selectedKeys[d.name];
				}
				thisBar.style("fill", newColor);
				self._triggerKeySelected();
			});
		
		bar.append("rect")
			.attr("width", barWidth)
			.attr("height", barHeight-1)
			.style("pointer-events", "none")
			.attr("fill", function(d){return d.color;});
		
		bar.append("text")
			.attr("x", barWidth + 3)
			.attr("y", barHeight / 2 - 2)
			.attr("dy", ".5em")
			.style("pointer-events", "none")
			.text(function(d){return d.name;});
	},
	
		_setOption: function( key, value ) {
		var self = this;
		fnMap = {
      'keyInfo': function () {
      	self.options[key]=value;
      	self.canvas.selectAll(".phenoColorBars")
      		.remove();
      	self.options[key] = value;
      	self._drawKey(self.options.keyInfo);
      }
    };

		if (key in fnMap) {
	    fnMap[key]();
	  }
  },
	
	_triggerKeySelected: function(keyName){
		var self=this;
		this._trigger('selected', {type: 'selected'}, {keyNames: self.selectedKeys});
	},
	
	clearSelections: function(){
		this.selectedKeys={};
// 		this.bgBars.style("fill", "rgb(255, 255, 255)");
		this.bgBars.style("fill", "white");
	},
	
	setSelection: function(keyNames){
		var self = this;
		self.selectedKeys = keyNames;
		self.bgBars.style("fill", function(d){
			if(!keyNames[d.name]){
				return self.whiteColorString;
			}
			else{
				return self.highlightedColor;
			}
		});
	}
	
});