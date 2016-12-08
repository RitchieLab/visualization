/**
* @module csg/shapekey
* @license
* Copyright (c) Marylyn Ritchie 2016
*
* This file is code for ShapeKey
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
 * @namespace csg.shapekey
 */

$.widget('csg.shapekey',{


	/**
	* @property {object} options The jQuery widget factory options
	* @property {number} options.height Overall height in pixels
	* @property {number} options.width Overall width in pixels
	* @property {number} options.barHeight Height of each row in the key
	* @property {array} options.keyInfo Array of objects containing name and symbol properties
	*/
	options: {
		height:  500, //overall height in pixels
		width:  150, //overall width in pixels
		barHeight: 18,
		keyInfo: []
	},


	/**
    * Standard jQuery widget factory _create function.  Utilizes the options object for
    * parameters and creating the widget.
    * @private
  */
	_create: function() {

		this.widgetID = this.element.attr('id');
		this.canvasID = 'csg-shapekey-'+this.widgetID;

		this.canvasjquery = $('<svg id="'+this.canvasID+'"></svg>').appendTo(this.element);
		this.canvas = d3.select('#'+this.canvasID)
			.attr("width", this.options.width)
			.attr("height", this.options.height)
			.append("g");

		this._drawKey(this.options.keyInfo);

		var whiteCheck = this.canvas.append("rect")
			.style("fill", "white");
		this.whiteColorString = whiteCheck.style("fill");
		whiteCheck.remove();
		this.highlightedColor = "#EDEDED";
	},

	/**
    * Draw key with shape to left, followed by text
    * @param {array} info Objects with name and shape
    * highlighted in the key.
    * @private
  */
	_drawKey: function(info){
		this.selectedKeys = {};
		var n = info.length;
		d3.select('#'+this.canvasID)
			.attr("height", this.options.barHeight * n);
		var barHeight = this.options.barHeight;
		var barWidth = this.options.width * 0.2;
		var self = this;

		var symbolSize = barHeight/2 * barHeight/2;

		var bar = this.canvas.selectAll("g")
			.data(info)
			.enter().append("g")
			.attr("class", 'shapeBars')
			.attr("transform", function(d, i){ return "translate(0," + i * barHeight + ")";});

		self.bgBars = bar.append("rect")
			.attr("width", this.options.width)
			.attr("height", barHeight)
			.style("fill", "rgb(255, 255, 255)")
			.on("click", function(d){
				var thisBar = d3.select(this);
				var fillColor = thisBar.style("fill");
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
		
		// add symbol 
 		bar.append("path")
 			.attr("transform", "translate(" + (barHeight/2).toString()  + "," + (barHeight/2).toString() + ")")
 			.attr("d", d3.svg.symbol()
 				.size(symbolSize)
 				.type(function(d){return d.symbolType;}))
 			.style("fill", "white")
 			.style("stroke", "black")
 			.style("stroke-width", 2);
			

// 		bar.append("rect")
// 			.attr("width", barWidth)
// 			.attr("height", barHeight-1)
// 			.style("pointer-events", "none")
// 			.attr("fill", function(d){return d.color;});

		bar.append("text")
			.attr("x", barWidth + 3)
			.attr("y", barHeight / 2 - 2)
			.attr("dy", ".5em")
			.style("pointer-events", "none")
			.text(function(d){return d.name;});
	},

	/**
    * Standard jQuery widget factory _setOption function
    * @private
  */
		_setOption: function( key, value ) {
		var self = this;
		fnMap = {
      'keyInfo': function () {
      	self.options[key]=value;
      	self.canvas.selectAll(".shapeBars")
      		.remove();
      	self.options[key] = value;
      	self._drawKey(self.options.keyInfo);
      }
    };

		if (key in fnMap) {
	    fnMap[key]();
	  }
  },

  /**
    * Selected key event
    *
    * @event csg.colorkey#selected
    * @type {object}
    * @property {object} keyNames - Object with selected key names
    */
	_triggerKeySelected: function(keyName){
		var self=this;
		this._trigger('selected', {type: 'selected'}, {keyNames: self.selectedKeys});
	},

	/**
    * Clear all selections and remove any highlighting for color.
  */
	clearSelections: function(){
		this.selectedKeys={};
		this.bgBars.style("fill", "white");
	},

	/**
    * Set selection and set highlighted colors for all selected key names.
    * @param {object} keyNames Contains all key names that will be marked as selected and
    * highlighted in the key.
  */
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