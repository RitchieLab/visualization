#!/usr/bin/env ruby

###################################################################################
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. <http://www.gnu.org/licenses/>
###################################################################################
ENV['MAGICK_CONFIGURE_PATH']='/gpfs/group/mdr23/usr/tools/etc/ImageMagick'
# Creates a plot showing LD as Haploview does
# with additional information related to the
# statistical analysis of the SNPs
# in the LD plot

# Requres rubygems
begin
  require 'rubygems'
rescue LoadError => e
  puts e
  puts "Please install rubygems -- http://docs.rubygems.org/read/chapter/3 "
  exit(1)
end

# Requires RMagick (http://rmagick.rubyforge.org/)
begin
  require 'rvg/rvg'
  rescue LoadError => e
  puts
  puts e
  puts "\nPlease install RMagick -- See documentation for ld_plus or http://rmagick.rubyforge.org/install-faq.html "
  puts
  exit(1)
end  

require 'optparse'
require 'ostruct'
include Magick


RVG::dpi=72

Version = '1.04'
# CHANGES
# 1.04 -- 3/3/15 Fixed syntax with case selection for compatibility with recent
# versions of ruby (thanks to Eric Fontanillas for identifying and providing fix).
# Added formatting for grid labels (d-prime and r-squared).

# check for windows and select alternate font based on OS
Font_family_style = RUBY_PLATFORM =~ /mswin/ ? "Verdana" : "Times"


#############################################################################
#
# Class Arg -- Parses the command-line arguments.
#
#############################################################################
class Arg

  def self.parse(args)
    options = OpenStruct.new
    options.ld_file = nil
    options.info_file = nil
    options.block_file = nil
    options.analysis_file = nil
    options.rsquared = nil
    options.dprime = nil
    options.snpinfo_file = nil
    options.out_name = 'ld_plus'
    options.highres = nil
    options.imageformat = 'png'

    help_selected = false
    version_selected = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: ld_plus.rb [options]"

      opts.on("-l [ldfile]", "Haploview LD output file") do |ldfile|
        options.ld_file = ldfile
      end
      opts.on("-i [info_file]", "Haploview information file with SNP positions") do |info_file|
        options.info_file = info_file
      end

      opts.on("-t [image_type]", "Image format for output (png default)") do |image_type|
        options.imageformat = image_type
      end

      opts.on("-r", "Draws grid for r-squared results") do |r|
        options.rsquared = true
      end

      opts.on("-d", "Draws grid for d-prime results") do |d|
        options.dprime = true
      end
      
      opts.on("-z", "Draws high resolution figure (300dpi)") do |hr|
        options.highres = true
      end


      opts.on("-b [block_file]", "Optional haploview file listing SNPs that make up blocks") do |block_file|
        options.block_file = block_file
      end
      
      opts.on("-o [output_name]", "Optional output name for the image") do |output_file|
        options.out_name = output_file
      end

      opts.on("-a [analysis_file]", "Optional file listing analyses to include in plot") do |analysis_file|
        options.analysis_file = analysis_file
      end
      
      opts.on("-s [snpinfo_file]", "Optional file listing additional snp information") do |snpinfo_file|
        options.snpinfo_file = snpinfo_file
      end
      opts.on("-f [feature_file]", "Optional file listing features such as exons") do |feature_file|
        options.feature_file = feature_file
      end
      opts.on_tail("-h", "--help", "Show this usage statement") do
        puts opts
        help_selected = true
      end
      opts.on_tail("-v", "--version", "Show version") do
        puts "\n\tVersion: #{Version}"
        version_selected = true
      end
    end

    begin
      opts.parse!(args)
    rescue StandardError => e
      puts e, "", opts
      exit(1)
    end

    if version_selected 
      puts
      exit(0)
    end

    if help_selected
      puts
      exit(0)
    end

    if !options.ld_file or !options.info_file # or !options.block_file
      help_string = opts.help
      puts "\n",help_string,"\n"
      puts "\nExample: ld_plus.rb -r -d -l haploview_output.txt -i sample.info -b block.txt\n\n"
      exit(1)
    end

    return options

  end
end


#############################################################################
#
# Class Stat -- Holds label and range information for the statistical results
#               held in the analysis file.
#
#############################################################################
class Stat
  attr_accessor :name, :min, :max, :labels

  def initialize(name, min, max, label)
    @name = name
    @min = min.to_f
    @max = max.to_f
    @labels = Array.new
    @labels << label
  end

end

#############################################################################
#
# Class Info -- Contains name for the SNP
#
#############################################################################
class Info
  attr_accessor :name

  def initialize(name)
    @name = name
  end

end

#############################################################################
#
# Class FeaturePart -- A feature specifies a region that will be displayed 
#                      above the bar showing SNP locations.
#
#############################################################################
class FeaturePart
  attr_accessor :name, :num, :startpos, :endpos

  def initialize(name, num, st, en)
    @name = name
    @num = num
    @startpos = st
    @endpos = en
  end

end


#############################################################################
#
# Class ChromosomeFeature -- contains multiple parts of a feature all the 
#                            features can be displayed on a single line above 
#                            the bar showing SNP locations
#
#############################################################################
class ChromosomeFeature
  attr_accessor :parts, :name
  
  def initialize(name)
    @parts = Array.new
    @name = name
  end
  
  def add_feature_part(num, stpos, endpos)
    @parts.each do |part|
      if num == part.num
        return
      end
    end
    @parts << FeaturePart.new(@name, num, stpos, endpos)
  end
  
  def add_part(part)
    parts << part
  end
  
end

#############################################################################
#
# Class ChromosomeFeatureList -- contains all chromosome features
#
#############################################################################
class ChromosomeFeatureList
  attr_accessor :features, :box_index
  
  def initialize 
    @features = Array.new
    @box_index = 0
  end
  
  def add_feature(name)
    @features.each do |feature|
      if name == feature.name
        return
      end
    end
    @features << ChromosomeFeature.new(name)
  end
  
  def add_feature_part(feature_part, box_set)
    @features.each_with_index do |feature, index|
      if feature.name == feature_part.name
        feature.add_part(feature_part)
        if box_set
          @box_index = index
        end
        return
      end
    end
    # if reached here need to add a new feature and then add the part
    newfeature = ChromosomeFeature.new(feature_part.name)
    newfeature.add_part(feature_part)
    @features << newfeature
    if box_set
      @box_index = @features.length - 1
    end
    
  end
  
  def move_featured_last
    if @features.length > 0
      tempfeat = @features.last
      @features[@features.length-1]= @features[@box_index]
      @features[@box_index] = tempfeat
    end
  end
  
 
end


#############################################################################
#
# Class SNP -- records information on snps such as name, location, and
#              results of statistical tests
#
#############################################################################
class SNP
  attr_accessor :name, :location, :stat_values, :snpinfo_values

  def initialize(name, location)
    @name = name
    @location = location.to_i
    @stat_values = Hash.new
    @snpinfo_values = Hash.new
  end

  def add_stat(stat, name)
    if !stat_values.has_key?(name)
      stat_values[name] = Array.new
    end
    stat_values[name] << stat
  end

  def add_snpinfo(info, name)
    snpinfo_values[name] = info
  end

end

#############################################################################
#
# Class SNPCombination -- lists combination information for two SNPs and 
#                         contains scores for that combination (r-squared,
#                         dprime, lod_score)
#
#############################################################################
class SNPCombination
  attr_accessor :snp1, :snp2, :rsquared, :dprime, :lod_score

  def initialize(snp1, snp2, dprime, lod_score, rsquared)
    @snp1 = snp1
    @snp2 = snp2
    @rsquared = rsquared.to_f
    @lod_score = lod_score.to_f
    @dprime = dprime.to_f
  end

end

#############################################################################
#
# Class SNPList -- holds SNPs along with other information such as LD 
#                  and what snps are included in the haploview output
#
#############################################################################

class SNPList
  attr_accessor :snps, :ld_scores, :stats, :snpinfo_tracks, :blocks, :snp_hash, :included_snps, :included_snps_hash,
    :index_of_included_snps, :maximum, :minimum

  def initialize()
    @snps = Array.new
    @ld_scores = Hash.new
    @min_position = 1000000000
    @max_position = 0
    @stats = Array.new
    @blocks = Hash.new
    @snpinfo_tracks = Array.new
    @snp_hash = Hash.new
    @included_snps = Array.new
    @included_snps_hash = Hash.new
    @index_of_included_snps = Hash.new
  end

  def add_snpinfo_track(name)
    @snpinfo_tracks.each do |track|
      if name == track.name
        return
      end
    end
    @snpinfo_tracks << Info.new(name)
  end

  def add_stat_info(name, min, max, label)
    @stats.each do |stat|
      if name == stat.name
        stat.labels << label
        return
      end
    end
    @stats << Stat.new(name, min, max, label)
  end

  def add_snp(snp)
    @snps << snp
    @snp_hash[snp.name] = snp
  end

  def add_block(markers)
    @blocks[snps[markers[0]-1].name] = Block.new(markers)
  end
  
  def add_new_block(block)
    @blocks[snps[block.snp_indices[0]-1].name] = block
  end

  # retrieve block
  def get_block(start_snp)
    return @blocks[start_snp]
  end

  # snps are sorted by basepair location
  def sort_snps
    @snps.sort!{|x,y| x.location <=> y.location}
  end

  # returns min position
  def get_min
    return snps[included_snps.first].location
  end

  # returns maximum number of iteractions for any SNP in list
  def get_max_interactions
    max_interaction = 0
    @ld_scores.values.each do |snp_combo_array|
      if snp_combo_array.length > max_interaction
        max_interaction = snp_combo_array.length
      end
    end
    return max_interaction
  end

  # returns max position
  def get_max
    return snps[included_snps.last].location
  end

  # return number of different statistics stored
  def get_total_stat_boxes
    return @stats.length
  end
  
  # return number of snpinfo tracks stored
  def get_total_snpinfo_boxes
    return @snpinfo_tracks.length
  end

  # adds ld score to the hash
  # key is the first snp
  def add_ld_score(snp_combination)
    if !@ld_scores.has_key?(snp_combination.snp1)
      @ld_scores[snp_combination.snp1] = Array.new
      @included_snps_hash[snp_combination.snp1] = true
      @included_snps_hash[snp_combination.snp2] = true
    end
    @ld_scores[snp_combination.snp1] << snp_combination
  end

  # sets up included snp array
  # to insure that only indicated snps are
  # included
  def process_included
    @snps.each_with_index do |snp, index|
      if @included_snps_hash[snp.name]
        @included_snps << index
        @index_of_included_snps[snp.name] = @included_snps.length-1
      end
    end
  end

  # argument is the original index of the SNP
  # returns the index of the SNP in the included list
  # these will be equal if no SNPs are excluded
  def get_snp_included_index_num(original_snp_index)
    return @index_of_included_snps[@snps[original_snp_index].name]
  end

  # returns the index number of the first included snp that is the one
  # larger than indicated position
  def get_index_first_larger_snp(pos)  
    @included_snps.each_with_index do |snp_indx, index|
      if pos < @snps[snp_indx].location.to_f
        return index
      end
    end
  end
  
  # returns the index number of the last included number smaller than
  # the position being checked
  def get_index_last_smaller_snp(pos)
    @included_snps.each_with_index do|snp_indx, index|
      if pos < @snps[snp_indx].location.to_f
        return index-1
      end
    end
  end


  def get_num_included_snps
    return @included_snps.length
  end

  # returns maximum number of haplotypes
  # for any block in the list
  def get_max_haplotype_size
    max_haplotype = 0
    @blocks.values.each do |block|
      if block.haplotypes.length > max_haplotype
        max_haplotype = block.haplotypes.length
      end
    end
    return max_haplotype
  end

end

#############################################################################
#
# Class Block -- LD block that will be highlighted on the plot.
#
#############################################################################
class Block
  attr_accessor :snp_indices, :snp_names, :haplotypes
  def initialize(snp_numbers)
    @snp_indices = snp_numbers
    @snp_names  = Array.new
    @haplotypes = Array.new
  end

  # haplotypes are sorted according to frequency
  def sort_haplotypes
    @haplotypes.sort!{|x,y| y.frequency <=> x.frequency}
  end

end

#############################################################################
#
# Class Haplotype -- Haplotype definition that stores the genotypes for the
#                    haplotypes and the frequency of the haplotypes.
#
#############################################################################
class Haplotype
  attr_accessor :bases, :frequency

  # pass hash with key as SNP name and value as base at
  # that place
  def initialize(bases, frequency)
    @bases = bases
    @frequency = frequency
  end

end


#############################################################################
#
# Class FileReader -- reads files that contain information for drawing plot.
#                     Each file type has its own method in this class.
#                     It creates objects of appropriate classes based
#                     on the file content.
# 
#############################################################################
class FileReader

  # reads file containing supplemental snp information
  def read_snpinfo_file(snp_info_file, snp_list)
    firstline = true
    header_names = Array.new

    begin
      File.open(snp_info_file, "r") do |file|
        file.each_line do |line|
          # skip blank lines
          if line !~ /\w/
            next
          end

          if !firstline
            # split line and check whether need
            # to create new holder for info
	    data=strip_and_split(line)
            data[1..data.length-1].each_with_index do |info,i|
	      snp_list.snp_hash[data[0]].add_snpinfo(info.to_i, header_names[i+1])
            end

	  else # read the headers
            firstline = false
            data=strip_and_split(line)
            header_names << data[0]
            data[1..data.length-1].each_with_index do |header,i|
              snp_list.add_snpinfo_track(header)
              header_names << data[i+1]
            end
          end

        end
      end

    rescue StandardError => e
      puts e
      exit(1)
    end
  end


  # reads file containing the definition
  # of the haplotype blocks and the freqency
  # of each haplotype
  def read_block_info_file(block_info_file, snp_list)

    haplotype_block = nil
    begin
      File.open(block_info_file) do |file|
        # skip header line
        line = file.gets
        haplotypes = Array.new
        first_line = true
        oldsnpid = "foo"

        while line = file.gets
          line.chomp!
          #Process Haplotype Definitions

            data=line.split # Splits the input line by tabs
            snpid = data[0];
            if first_line or snpid != oldsnpid
               haplotype_block = snp_list.get_block(snpid)
               first_line = false
            end

            temphap = data[1].split(//)
            frequency =  data[2]

            data_length = temphap.length # How big is the haplotype?
            haplotype_block.haplotypes << Haplotype.new(temphap, frequency.to_f)
            oldsnpid = snpid

        end
      end

    rescue StandardError => e
      puts e
      exit(1)
    end

  end

  # reads haploview format file
  # containing information on the SNPs that make
  # up the haplotype block
  def read_block_file(blockfile, snp_list)

    letters = ["N","A","C","G","T"]

    begin
      File.open(blockfile) do |file|
        while line = file.gets
          if line =~ /BLOCK/
            line.chomp!
            pieces = line.split
            pieces = pieces[3..pieces.length-1]
            integers = Array.new
            pieces.each do|piece|
              integers << piece.to_i
            end
                        
            newblock = Block.new(integers)
            line = file.gets
            
            # read haplotype frequency and block genotypes
            total_frequency = 0.0
            while line and line !~ /Multiallelic/    
              line.chomp!
              pieces = line.split
              genotypes = pieces[0].split('')
              letter_genos = []
              # convert to letters from numbers (1-4)
              genotypes.each do |g|
                g = letters[g.to_i]
                letter_genos << g
              end
              
              # frequency will be the second element
              frequency = pieces[1]
              # remove parentheses
              frequency.delete!("(")
              frequency.delete!(")")

              frequency = frequency.to_f
              
              newblock.haplotypes << Haplotype.new(letter_genos, frequency.to_f)
              line = file.gets
            end
            snp_list.add_new_block(newblock)
          end
        end

      end
    rescue StandardError => e
      puts e
      exit(1)
    end

  end

  # read in snp names and positions
  # uses same format as haploview
  def read_input_file(infofile, snp_list)

    begin
      File.open(infofile, "r") do |file|
        file.each_line do |line|
          # skip blank lines
          if line !~ /\w/
            next
          end
            # split line and check whether need
            # to create new holder for plot
            data=strip_and_split(line)
            snp_list.add_snp(SNP.new(data[0], data[1]))
        end
      end
    rescue StandardError => e
      puts e
      exit(1)
    end
  end

  # reads the optional chromosome features file
  # it has a format as follows
  # FeatureName  FeatureNumber  StartPos  EndPos
  # Exon         1              507878    509001
  def read_feature_file(featurefile, chrom_feature_list)
    begin
      box_set=false
      File.open(featurefile, "r") do |file|
        file.each_line do |line|
          # skip blank lines
          if line !~ /\w/
            next
          end
          
          if line =~ /\*/
            box_set=true
            line = line.sub(/\*/,'')
          end
          
          # split line and insert feature into feature list
          data = strip_and_split(line)
          chrom_feature_list.add_feature_part(FeaturePart.new(data[0], data[1], data[2], data[3]), box_set)
          box_set=false
        end
      end
    rescue StandardError => e
      puts e
      exit(1)
    end
    
    # move one to draw boxes into first position of array
    chrom_feature_list.move_featured_last
    
  end

  # reads file for analysis results
  # to use in creating plots
  def read_analysis_file(analysis_file, snp_list)
    firstline = true
    header_names = Array.new

    begin
      File.open(analysis_file, "r") do |file|
        file.each_line do |line|
          # skip blank lines
          if line !~ /\w/
            next
          end

          if !firstline
            # split line and check whether need
            # to create new holder for plot
            data=strip_and_split(line)
            data[1..data.length-1].each_with_index do |stat,i|         
              snp_list.snp_hash[data[0]].add_stat(stat.to_f, header_names[i+1])
            end
          else # read the headers
            firstline = false
            data=strip_and_split(line)
            header_names << data[0]
            data[1..data.length-1].each do |header|
              header_info = header.split
              header_names << header_info[0]
              snp_list.add_stat_info(header_info[0], header_info[1], header_info[2], header_info[3])
            end
          end

        end
      end
    rescue StandardError => e
      puts e
      exit(1)
    end

  end


  # reads LD file and returns array
  # of LD combinations
  def read_ld_output_file(ldfile)
    snp_lds = Array.new

    begin
      File.open(ldfile, "r") do |file|
        file.each_line  do |line|
          if line =~ /LOD/
            next
          end

          data=strip_and_split(line)
          snp_lds << SNPCombination.new(data[0], data[1], data[2], data[3], data[4])

        end
      end
    rescue StandardError => e
      puts e
      exit(1)
    end

    return snp_lds

  end

  # strips and splits the line
  # returns the data array that results
  def strip_and_split(line)
    line.strip!
    line.split(/\t/)
  end

end

##### End FileReader #####



#############################################################################
#
# Class PlotWriter -- Contains functions for drawing plots.  Most distances
#                     are adjusted based on the @box_size.  The @box_size
#                     defines the size in x and y coordinates for the 
#                     elements in the D-prime and R-squared plots.
#
#############################################################################
class PlotWriter
  attr_accessor :box_size, :canvas, :color_array, :dashed_array, :interval_array, :feature_opacity

  def initialize(box_size)
    @box_size = box_size
    @color_array = ['black','red', 'blue','gray','orange','green','purple']
    @dashed_array = [3, 10, 17, 23, 29, 34, 39]; 
    @interval_array = [3, 3, 6, 6, 9, 9, 9];
    @feature_opacity = 1.0
  end

  # took calculations from the haploview code
  # for the block colors
  def get_ld_block_color(snp_combo, use_dprime=true)

    color = 'white'
    if use_dprime then
      dprime = snp_combo.dprime
      lod_score = snp_combo.lod_score
      color = get_color_from_score(dprime, lod_score, use_dprime)
    else # use rsquared
      rsquared = snp_combo.rsquared
      color = get_color_from_score(rsquared, 0, use_dprime)
    end

    return color

  end
  
  # returns appropriate color based on 
  def get_color_from_score(score, lod_score, use_dprime)
    color = 'white'
    if use_dprime then
      if lod_score > 2
        if score < 0.5
          # high LOD, low D'
          color = "rgb(255,224,244)"
        else
          # high LOD, high D' => shades of red
          blue_green = (255-32)*2*(1-score)
          color = sprintf("rgb(255, %d, %d)", blue_green, blue_green)
        end
      elsif(score > 0.99)
        #high LD, low LOD blueish color
        color = "rgb(192,192,240)"
      else
        color = 'white'
      end
    else # use rsquared
      red_blue_green = 255 * (1.0 - score)
      color = sprintf("rgb(%d,%d,%d)", red_blue_green, red_blue_green, red_blue_green)
    end  
    return color
  end

  # draws a legend along left side for the indicated plot
  def draw_plot_legend(x_start, y_start, is_dprime=true)
    
    labels_array = Array.new
      for i in 0..10 
        labels_array << i * 0.1
      end

    current_x=0
    current_y=0
    
    # draw boxes from starting point
    @canvas.g.translate(x_start,y_start) do |legend|
      labels_array.reverse_each do |box_value|
        color_string = get_color_from_score(box_value, 2.1, is_dprime)
        legend.rect(@box_size, @box_size, current_x, current_y).styles(:stroke=>color_string,
          :stroke_width=>1, :fill=>color_string)
        if is_dprime
          color_string = get_color_from_score(box_value, 0, is_dprime)
          legend.rect(@box_size, @box_size, current_x-@box_size, current_y).styles(:stroke=>color_string,
            :stroke_width=>1, :fill=>color_string)
        end
        current_y+=@box_size
      end
    end
    
    current_y=@box_size/2
    font_size = standard_font_size
    
    if is_dprime
      current_x = current_x - @box_size
    end
    
    labels_array.reverse_each do |label|
      @canvas.g.translate(x_start, y_start).text(current_x-2, current_y) do |text|
        text.tspan("%.1f" % [label]).styles(:font_size=>font_size/1.2, :text_anchor=>'end')
      end
      current_y+=@box_size 
    end
    
    # place box around the legend
    current_y=0
    width = @box_size
    if is_dprime
      width = @box_size*2
    end
    @canvas.g.translate(x_start, y_start) do |legend|
      legend.rect(width, @box_size*11, current_x, current_y).styles(:stroke=>'black',
        :stroke_width=>1, :fill=>'none')
    end
    
    # write labels on top of columns for d-prime legend
    if is_dprime
      current_y = -2
      current_x = current_x + @box_size/2
      @canvas.g.translate(x_start, y_start).text(current_x, current_y) do | text|
        text.tspan("<=2").styles(:font_size=>font_size/1.2, :text_anchor=>'middle')
      end
      current_x = current_x + @box_size
      @canvas.g.translate(x_start, y_start).text(current_x, current_y) do | text|
        text.tspan(">2").styles(:font_size=>font_size/1.2, :text_anchor=>'middle')
      end
      current_y = current_y - @box_size/2
      @canvas.g.translate(x_start, y_start).text(current_x-@box_size/2, current_y) do | text|
        text.tspan("LOD").styles(:font_size=>font_size/1.2, :text_anchor=>'middle')
      end
    end
    
  end


  # draws transparent boxes over the stat boxes
  def draw_stat_box_blocks(x_start, ystart_stat_boxes, yend_stat_boxes, snp_list)

    x = @box_size + @box_size/2
    yheight =  yend_stat_boxes - ystart_stat_boxes

    snp_list.blocks.each do |first_snp_name, block|

      start_block_index = curr_block_member = 0
      last_snp = block.snp_indices[curr_block_member]-1

      last_index = block.snp_indices.length-1

      while curr_block_member <= last_index
        # continue with block
        if block.snp_indices[curr_block_member] == last_snp + 1

          last_snp = block.snp_indices[curr_block_member]
          curr_block_member += 1
        else
          # check to see if the missing numbers are all not in
          # the list (have been dropped by haploview)
          # if so this is still contiguous and proceed
          included_in_plot_not_block = -1
          end_block = -1

          for i in last_snp+1...block.snp_indices[curr_block_member]
            if(snp_list.included_snps_hash[snp_list.snps[i-1].name])
              # this one is included but is skipped in this block
              included_in_plot_not_block = i
              end_block = block.snp_indices[curr_block_member-1]
              break
            end
          end


          # need to draw the standard block
          if included_in_plot_not_block > 0
            if end_block > 0
              x_points = calculate_block_x(block.snp_indices[start_block_index],
                block.snp_indices[curr_block_member-1], snp_list)

              @canvas.rect(x_points[1] - x_points[0], yheight, x_start+x_points[0],
                ystart_stat_boxes).styles(:fill_opacity => 0.2, :fill => 'gray')

              last_snp = block.snp_indices[curr_block_member]-1
              start_block_index = curr_block_member
            end

            # need to draw the lighter block indicating these SNPs split the larger block
            final_included_in_plot_not_block = -1
            for i in included_in_plot_not_block...block.snp_indices[curr_block_member]
              if(snp_list.included_snps_hash[snp_list.snps[i-1].name])
                final_included_in_plot_not_block = i
              end
            end

            x_points = calculate_block_x(included_in_plot_not_block,
              final_included_in_plot_not_block, snp_list, true)
              @canvas.rect(x_points[1] - x_points[0], yheight, x_start+x_points[0],
                ystart_stat_boxes).styles(:fill_opacity => 0.2, :fill => 'skyblue')
          end
          last_snp = block.snp_indices[curr_block_member]-1
        end

      end
      # draw final block
      x_points = calculate_block_x(block.snp_indices[start_block_index],
        block.snp_indices.last, snp_list)
      @canvas.rect(x_points[1] - x_points[0], yheight, x_start+x_points[0],
        ystart_stat_boxes).styles(:fill_opacity => 0.2, :fill => 'gray')

    end


  end

  # defines start and end of the block
  # inner block should be true when the calculation is to
  # fill in the lighter ones within an outer block
  def calculate_block_x(start_snp, end_snp, snp_list, inner_block=false)
    x_points = Array.new

    # convert the start and end snps to positions in the
    # included snps

    start_snp = snp_list.index_of_included_snps[snp_list.snps[start_snp-1].name]
    end_snp = snp_list.index_of_included_snps[snp_list.snps[end_snp-1].name]

    left_inner_adjustment = 2.4
    right_inner_adjustment = 1.2
    if inner_block
      left_inner_adjustment = 1.5
      right_inner_adjustment = 1.05
    end

    box_adjustment = @box_size/left_inner_adjustment

    if start_snp == 0
      box_adjustment = @box_size/3.5
    end

    if start_snp != end_snp
      x_points << step_size * (start_snp+1) - box_adjustment
      x_points << step_size * (end_snp+1) + @box_size/right_inner_adjustment
    else
      x_points << step_size * (start_snp+1) - box_adjustment
      x_points << step_size * (start_snp+1) + @box_size/right_inner_adjustment
    end
    return x_points
  end

  # add stat box with lines
  # showing results for the current stat
  def draw_stat_box(current_stat, snp_list, x_start, y_start, x_left_boundary)

    x = @box_size + @box_size/2
    xmin = @box_size + @box_size/2
    xmax = xmin + step_size * (snp_list.included_snps.length-1)
    ymin = 0
    ymax = @box_size *3
    y_interval = ymax - ymin

    stat_max = snp_list.stats[current_stat].max
    stat_min = snp_list.stats[current_stat].min
    stat_interval = stat_max - stat_min

    # create an array to make the poly line
    stat_name = snp_list.stats[current_stat].name
    line_array = Array.new

    # need to check how many lines needed for this stat
    # use number of values in first snp
    # line_array is a 2-D array
    stat_values = snp_list.snps.first.stat_values[stat_name]
    num_lines = 0
    stat_values.each do |val|
      line_array << Array.new
      num_lines +=1
    end

    # generate a line for each
    # add labels for scale on y axis
    # add title for the box
    value_x = xmin
    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      stat_values = snp.stat_values[stat_name]
      stat_values.each_with_index do |value, i|
        value_y = ((stat_max-value) / stat_interval) * y_interval
        line_array[i] << value_x << value_y
      end
      value_x = label_step(value_x)
    end

    # draw line
    num_lines.times do |curr_line_index|
      @canvas.g.translate(x_start, y_start) do |lines|
        lines.polyline(line_array[curr_line_index]).styles(:fill=>'none', :stroke_width=>1,
          :stroke=>@color_array[curr_line_index], 
          :stroke_dasharray=>[@dashed_array[curr_line_index],@interval_array[curr_line_index]],
          :fill_opacity => 0)
      end
    end

    # draw scale (split into 5 numbers)
    stat_interval = stat_max-stat_min
    stat_break = stat_interval/4
    y_break = y_interval/4

    current_stat_value = stat_min
    y = ymax

    # need line at half-way point
    @canvas.g.translate(x_start, y_start) do |lines|
      lines.line(xmin, ymax/2, xmax, ymax/2).styles(:stroke=>'lightslategray',
        :stroke_width=>1, :stroke_dasharray=>[@box_size/5+1,@box_size/2], :fill=>'none')
    end
    
    # add tick marks to help viewing of scale on the analysis
    y_mark_start = 0
    x_mark_start = xmin - @box_size/2
    5.times do |count|
      @canvas.g.translate(x_start, y_start) do |mark|
        mark.line(x_mark_start, y_mark_start, x_mark_start+3, y_mark_start).styles(
          :stroke=>'black', :stroke_width=>1)
        y_mark_start += y_break
      end
    end

     # write labels on left side of stat box
     write_stat_box_labels(x_start, y_start, ymax, snp_list.stats[current_stat], y_interval, 0, false)
  end

  # adds stat box labels
  def write_stat_box_labels(x_start, y_start, ymax, stat, y_interval, x_left_boundary, numbers_first=false)

    x_adjustment = standard_font_size / 5 + 1

    x_start -= x_adjustment
  
    # now write them as right-aligned and horizontally 
    stat_interval = stat.max-stat.min
    stat_break = stat_interval/4
    y_break = y_interval/4
    y=ymax
    current_stat_value = stat.min
    font_size = standard_font_size   

    num_x = @box_size
    letters_x = @box_size/2
    
    4.times do |curr_break|
      @canvas.g.translate(x_start, y_start).text(num_x, y) do |text|
        label = compress_number(current_stat_value)
        text.tspan(label).styles(:font_size=>font_size/1.2, :text_anchor=>'end')
        # add tick mark
        current_stat_value += stat_break
        y -= y_break        
      end
    end
    
    @canvas.g.translate(x_start, y_start).text(num_x, y) do |text|
      text.tspan(compress_number(stat.max)).styles(:font_size=>font_size/1.2, :text_anchor=>'end')
    end

    x_start -= calculate_coordinate(0.0005 * box_size * 18)

    # add labels
    @canvas.g.translate(x_start, y_start).text(num_x, y) do |text|
      text.tspan(stat.name).styles(:font_size=>font_size*0.8, :text_anchor=>'end')
    end
    
    # need label and line
    # line needs to start right after label and continue to x_start
    if stat.labels.length > 1
      legend_x = x_left_boundary + @box_size
      y += @box_size/4
      stat.labels.each_with_index do |label, index|
        y += @box_size/3
        @canvas.g.translate(legend_x, y_start).text(0, y) do |text|
          text.tspan(label).styles(:font_size=>font_size*0.55, :text_anchor=>'end')
        end
        @canvas.g.translate(x_left_boundary, y_start) do |l|
          l.styles(:fill=>'none', :stroke_width=>1,
            :stroke=>@color_array[index],
            :stroke_dasharray=>[@dashed_array[index],@interval_array[index]],
            :fill_opacity => 0.5)
          l.line(legend_x+1, y-@box_size/6, x_start+@box_size-x_adjustment, y-@box_size/6)
        end
      end
      
    end
    
  end


  def compress_number(num)

    num_string = sprintf("%0.2f", num)
    if(num_string =~ /^00$/)
      num_string = "0  "
    end

    return num_string
  end

  # steps horizontally across image
  # lines up steps
  def label_step(current_x)
    return current_x + label_increment
  end
  
  def label_increment
    return @box_size * Math.sqrt(2)
  end

  # returns step size between snps in diagram
  def step_size
    return @box_size * Math.sqrt(2)
  end

  def standard_font_size
    return @box_size/2 + 1
  end

  # add labels for SNPs and physical distance
  # conversion
  def add_labels(snp_list, x_start, y_start, features_included=false)
    x_text = @box_size + @box_size/2
    start_x = x_text
    end_x = start_x
    y_text = @box_size*7
    font_size = standard_font_size

    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      @canvas.g.translate(x_start,y_start).text(x_text, y_text).rotate(-90) do |text|
          text.tspan(snp.name).styles(:font_size=>font_size)
      end
      end_x = x_text
      x_text = label_step(x_text)
    end

    # add the physical distance bar and
    # draw lines from that to the titles above
    min_pos = snp_list.get_min
    max_pos = snp_list.get_max
    interval_pos = max_pos - min_pos

    x_interval = end_x - start_x

    box_y_start = @box_size
    box_height = @box_size/3

    # add min and max pos to bar on either side
    @canvas.g.translate(x_start, y_start).text(start_x-@box_size/2, box_height+box_y_start) do |text|
      text.tspan(snp_list.get_min.to_s).styles(:font_size=>font_size, :text_anchor=>'end')
    end
    
    @canvas.g.translate(x_start, y_start).text(end_x+@box_size/2, box_height+box_y_start) do |text|
      text.tspan(snp_list.get_max.to_s).styles(:font_size=>font_size, :text_anchor=>'start')
    end
    

    # draw narrow rectangle with white fill
    @canvas.g.translate(x_start, y_start) do |chromosome|
      chromosome.styles(:fill=>'white', :stroke_width=>1, :stroke=>'black')
      chromosome.rect(x_interval, box_height, start_x, box_y_start)
    end

    #original multiplier was 3.5 -- Feb 2009 change
    if features_included
      y_text_line = y_start + @box_size*2
    else
      y_text_line = y_start + @box_size*3
    end

    x_text_line = @box_size + @box_size/2 + 1

    # draw lines connecting fonts to the position
    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      # determine relative position
      x_end_position = ((snp.location.to_f - min_pos) / interval_pos) * x_interval + start_x
      @canvas.g.translate(x_start, y_start) do |pos_line|
        # draw a vertical line across the box
        pos_line.styles(:stroke=>'gray', :stroke_width=>1)
        pos_line.line(x_end_position, box_y_start, x_end_position, box_y_start+box_height)
        pos_line.line(x_text_line, y_text_line, x_end_position, box_y_start+box_height)
      end

      x_text_line = label_step(x_text_line)
    end

  end


  # adds the gene features such as exons above the label track showing the gene positions
  # a feature can be made up of multiple parts
  # all parts of a feature occupy a single line on the plot
  def draw_features(current_feature, chrom_feature_list, snp_list, x_start, y_start, draw_boxes, y_start_label, y_dprime_start, y_lower_feature_start)

    x_text = @box_size + @box_size/2
    start_x = x_text
    end_x = start_x
    y_text = @box_size*7
    font_size = standard_font_size  

    snp_x_pos = Array.new

    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      end_x = x_text
      snp_x_pos << x_text
      x_text = label_step(x_text)
    end

    # get information on scale
    min_pos = snp_list.get_min
    max_pos = snp_list.get_max
    interval_pos = max_pos - min_pos
    
    x_interval = end_x - start_x
    
    box_height = @box_size/3
    box_y_start = @box_size
   
    feature = chrom_feature_list.features[current_feature]
    opacity = @feature_opacity
    
    stroke_dash_array = Array.new
    stroke_dash_array << [@box_size/2, @box_size/4]
    stroke_dash_array << [@box_size/5, @box_size/5]
    curr_dash = 0
    
    feature.parts.each do |part|
      # draw narrow rectangle with grey fill
      part_width = ((part.endpos.to_f - part.startpos.to_f) / interval_pos) * x_interval
      part_x = ((part.startpos.to_f - min_pos) / interval_pos) * x_interval + start_x

      @canvas.g.translate(x_start, y_start) do |chromosome|
        chromosome.styles(:fill=>'saddlebrown', :stroke_width=>1, :stroke=>'saddlebrown', :opacity=>opacity)
        chromosome.rect(part_width, box_height, part_x, box_y_start)
      end
      
      y_feat_text = box_y_start - box_height/2 
      text_out = "#{part.name} #{part.num}"
      @canvas.g.translate(x_start,y_start).text(part_x, y_feat_text) do |text|
          text.tspan(text_out).styles(:font_size=>font_size)
      end 
      
      # need to draw boxes that extend all the way 
      if draw_boxes
        
        # draw between labels
        left_index = snp_list.get_index_first_larger_snp(part.startpos.to_f)
        right_index = snp_list.get_index_last_smaller_snp(part.endpos.to_f)
        left_x = snp_x_pos[left_index];
        left_x = left_x - label_increment/2
        right_x = snp_x_pos[right_index]
        right_x = right_x + label_increment/2
        
        # determine position of line 
        start_x_feat = part_x
        end_x_feat = part_x + part_width

    
        @canvas.g.translate(x_start, y_start) do |pos_line|
          pos_line.styles(:fill=>'none', :stroke=>'saddlebrown', :stroke_width=>1, :stroke_dasharray=>stroke_dash_array[curr_dash],
            :stroke_opacity=>opacity)
        
          pos_line.line(start_x_feat, box_y_start+box_height ,start_x_feat, y_start_label-y_start+box_height+box_y_start)
          pos_line.line(end_x_feat, box_y_start+box_height,end_x_feat, y_start_label-y_start+box_height+box_y_start)
          
          pos_line.line(start_x_feat, y_start_label-y_start+box_height+box_y_start, left_x, y_start_label+y_text-y_start-@box_size*3)
          pos_line.line(end_x_feat, y_start_label-y_start+box_height+box_y_start, right_x, y_start_label+y_text-y_start-@box_size*3)
          
          #draw straight line from labels down to d-prime chart
          pos_line.line(left_x, y_start_label+y_text-y_start-@box_size*3,left_x, y_dprime_start-y_start);
          pos_line.line(right_x, y_start_label+y_text-y_start-@box_size*3,right_x, y_dprime_start-y_start);
                    
        end 
        
        # draw boxes just above the d-prime or r-squared plots
        @canvas.g.translate(x_start, y_start) do |lower_box|
          lower_box.styles(:fill=>'saddlebrown', :stroke_width=>1, :stroke=>'saddlebrown', :opacity=>opacity)
          lower_box.rect(right_x-left_x, box_height, left_x, y_dprime_start-y_start)
        end
        
      end
      
      if opacity == @feature_opacity
        opacity = @feature_opacity * 0.75
        curr_dash = 1
      else
        opacity = @feature_opacity
        curr_dash = 0
      end         
    end
   
  end

  # draws gene information track
  # PARAMETERIZE FOR GENE INFORMATION STRUCTURE
  def add_gene(x_start, y_start)
    box_x_start = @box_size + @box_size/2
    box_y_start = @box_size
    box_height = @box_size/3

    #Will's testing junk
    @canvas.g.translate(x_start, y_start) do |gene|
      gene.styles(:fill=>'yellow', :stroke_width=>1, :stroke=>'black')
      gene.rect(100, box_height, box_x_start, box_y_start)
    end

  end


  # draws the snp information grid
  # PARAMETERIZE FOR SNP INFORMATION STRUCTURE
  def draw_snpinfo(current_track, snp_list, x_start, y_start)

    box_y_start = 0 #@box_size * Math.sqrt(2) *0.25 #@box_size/6

    x_text = @box_size * 1.15

    color_array = ['black','red', 'blue','gray','orange','green','purple','brown','skyblue','fuchsia','lime']

    track_name = snp_list.snpinfo_tracks[current_track].name

    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      val = snp.snpinfo_values[track_name]
	if val == 1
          @canvas.g.translate(x_start,y_start) do |check|
	    check.styles(:fill=>color_array[current_track], :stroke=>'white', :stroke_width=>1)
	    check.rect(@box_size * Math.sqrt(2) *0.25, @box_size * Math.sqrt(2) *0.25, x_text, box_y_start)
          end
	end

        end_x = x_text
        #x_text += @box_size + @box_size/2 - @box_size/10
        x_text = label_step(x_text)


    end

    @canvas.g.translate(x_start, y_start) do |checks|
      #checks.line(box_x_start, box_y_start , x_text, box_y_start)
      #checks.line(box_x_start, box_y_start + @box_size * Math.sqrt(2) *0.25 , x_text, box_y_start + @box_size * Math.sqrt(2) *0.25)
    end

    font_size = standard_font_size

    @canvas.g.translate(x_start,y_start).text(@box_size, font_size * 0.8) do |text|
      text.tspan(track_name).styles(:font_size=>font_size, :text_anchor=>'end')
    end

  end



  # draws the grid
  def draw_grid(snp_list, x_start, y_start, use_dprime=true)
    y = y_start+@box_size
    x = x_start
    
    if use_dprime
      grid_label = "D-prime"
    else
      grid_label = "R-Squared"
    end
    
    label_x = x_start - standard_font_size / 5 + 1
    font_size = standard_font_size
 
    x_text = @box_size + @box_size/2 + x_start
    y_text = y_start + @box_size

    # add label for R-squared and D-prime
    @canvas.g.text(label_x, y_text) do |text|
      text.tspan(grid_label).styles(:font_size=>font_size, :text_anchor=>'end')
    end

    # draw plot legend below grid label
    draw_plot_legend(label_x, y_text+@box_size*1.5, use_dprime)

    snp_list.included_snps.each do |snp_index|
      @canvas.text(x_text, y_text) do |number|
        number.tspan((snp_index+1).to_s).styles(:text_anchor =>'middle', :font_size => standard_font_size,
          :font_family=>Font_family_style, :fill=>'black')
        x_text = label_step(x_text)
      end
    end

     y_current_start = @box_size
     x = 0

  default_stroke_color = "lightgray"
  stroke_color = default_stroke_color
  use_block_stroke = true
  if(@box_size >= 7)
    use_block_stroke = false
  end

  @canvas.g.translate(x_start,y_start).rotate(-45) do |ldbox|
    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      y_current_start += @box_size
      y = y_current_start
      if snp_list.ld_scores.has_key?(snp.name)
        snp_list.ld_scores[snp.name].each do |snp_combo|
          block_color = get_ld_block_color(snp_combo, use_dprime)
          if use_block_stroke
            stroke_color = block_color
          else
            stroke_color = default_stroke_color
          end
          ldbox.rect(@box_size, @box_size, x, y).styles(:stroke=>stroke_color, :stroke_width=>1,
            :fill=>block_color)
            y += @box_size
         end
       end
       x += @box_size
     end
   end


    snp_list.blocks.each do |first_snp_name, block|
    first_snp = snp_list.index_of_included_snps[snp_list.snps[block.snp_indices.first-1].name]
    last_snp = snp_list.index_of_included_snps[snp_list.snps[block.snp_indices.last-1].name]
      @canvas.g.translate(x_start,y_start).rotate(-45) do |draw_block|
        x_line_start = @box_size * (first_snp)
        y_line_start = @box_size * (first_snp+2) - @box_size

        x_line_end = @box_size * (last_snp) + @box_size
        y_line_end = @box_size * (last_snp+2)

      # find last point for line
      # line will span total # of snps - 1 number of blocks
        num_blocks = last_snp - first_snp
        x_intersection = x_line_start
        y_intersection = @box_size*(first_snp+2) + @box_size * num_blocks

      # now draw polyline
        draw_block.polyline(x_line_start, y_line_start, x_intersection, y_intersection,
          x_line_end, y_line_end).styles(:fill=>'none', :stroke=>'black', :stroke_width=>@box_size/5+1)
      end
    end


  end

  # adds haplotypes to the
  def draw_haplotypes(x_start, ystart_haplotypes, snp_list, include_label=true)

    label_x = x_start - standard_font_size / 5 + 1
    font_size = standard_font_size
 
    x_text = @box_size + label_x
    y_text = ystart_haplotypes + @box_size

    # add label for haplotypes
    if include_label
      @canvas.g.text(x_text, y_text) do |text|
        text.tspan("Haplotypes").styles(:font_size=>font_size, :text_anchor=>'end')
      end
    end


    y_text = 0
    snp_list.blocks.values.each do |block|
      y_text = ystart_haplotypes + @box_size/2
      block.sort_haplotypes


      #determine block length for use
      block.haplotypes.each do |haplotype|
        haplotype.bases.each_with_index do |base, snp_index|

          x_text = x_start + @box_size + @box_size/2 + step_size *
            (snp_list.get_snp_included_index_num(block.snp_indices[snp_index]-1))

          @canvas.text(x_text, y_text) do |letter|
            letter.tspan(base).styles(:text_anchor => 'middle', :font_size => standard_font_size,
              :font_family => 'arial', :fill=>'black')
	  end
        end

        start_snp = block.snp_indices.first
        end_snp = block.snp_indices.last

        x_points = calculate_block_x(start_snp, end_snp, snp_list)
        block_width = x_points[1] - x_points[0]
        frequency_x_start = x_start + x_points[0]
        frequency_x_end = haplotype.frequency * block_width
        @canvas.rect(frequency_x_end, @box_size/2, frequency_x_start, y_text-@box_size/2).styles(:fill=>'gray',
          :fill_opacity=>0.5)
        y_text += @box_size/2
      end
    end

  end


  def draw_border(width, height)
    @canvas.rect(width-2, height-2).styles(:stroke=>'black', :stroke_width=>2, :fill=>'none')
  end

  # returns plot width in inches
  def calculate_plot_width(num_snps, right_padding)
    return num_snps * @box_size.to_f/100 + right_padding
  end

  # returns value in x,y coordinate for the
  # size in inches passed
  def calculate_coordinate(size_in_inches)
    return size_in_inches * 144
  end

  def add_space_for_gene(size_in_inches, fraction)
    return size_in_inches + @box_size * fraction
  end

  def add_space_for_labels(size_in_inches, fraction)
    return size_in_inches + @box_size * fraction
  end

  def add_space_for_snpinfo(size_in_inches, fraction)
    return size_in_inches + @box_size * fraction
  end

  def add_haplotype_space(size_in_inches, max_haplotype)
    return size_in_inches + 0.00325 * @box_size * max_haplotype + @box_size.to_f/200
  end

  def add_grid_size(size_in_inches, horizontal_side_in_inches, maximum_interactions, num_snps)
    xside_fraction = Math.sqrt(maximum_interactions / Math.sqrt(2) *
      (maximum_interactions / Math.sqrt(2)) / 2) / num_snps
    return  size_in_inches + horizontal_side_in_inches * xside_fraction  + 0.05
  end

end

##### End PlotWriter #####

# Main execution begins here
options = Arg.parse(ARGV)

# increase dpi if highres is selected
if options.highres
  RVG::dpi=300
end

filereader = FileReader.new

# read snp list and ld scores from haploview
snp_list = SNPList.new
filereader.read_input_file(options.info_file, snp_list)
snp_list.sort_snps
ld_scores = filereader.read_ld_output_file(options.ld_file)
ld_scores.each do |combo|
  snp_list.add_ld_score(combo)
end
snp_list.process_included

chrom_feature_list = ChromosomeFeatureList.new

# Read block definitions file if provided
if options.block_file
  filereader.read_block_file(options.block_file, snp_list)
end

# box_size is basic unit of size in the plot
# it is adjusted for total number of snps in the set
box_size = case
           when snp_list.get_num_included_snps > 500 then 2
           when snp_list.get_num_included_snps > 350 then 4
           when snp_list.get_num_included_snps > 275 then 6
           when snp_list.get_num_included_snps > 200 then 8
           when snp_list.get_num_included_snps > 175 then 10
           when snp_list.get_num_included_snps > 150 then 12
           when snp_list.get_num_included_snps > 125 then 14
           when snp_list.get_num_included_snps > 100 then 16
           when snp_list.get_num_included_snps > 80 then 18
           when snp_list.get_num_included_snps > 60 then 20
           when snp_list.get_num_included_snps > 40 then 24
           else 28
           end

writer = PlotWriter.new(box_size)

xside_end_addition = 0.005 * box_size

# read the analysis file if one
if(options.analysis_file)
  filereader.read_analysis_file(options.analysis_file, snp_list)
  xside_end_addition = 0.012 * box_size + 0.007 * box_size
end

x_start = 0

if(options.snpinfo_file || options.analysis_file)
  max_label_length = 0
  if options.snpinfo_file 
    filereader.read_snpinfo_file(options.snpinfo_file, snp_list)
    snp_list.snpinfo_tracks.each do |track|
      if track.name.length > max_label_length
        max_label_length = track.name.length
      end
    end
  end
  
  if options.analysis_file
    snp_list.stats.each do |stat|
      if stat.name.length > max_label_length
        max_label_length = stat.name.length
      end
    end
  end
  
  max_label_length += 12 # add some space for the numbers on the stat grid
  xleftside_addition = (0.0005 * box_size * max_label_length) + 0.012 * box_size
  xside_end_addition = xside_end_addition + xleftside_addition
  x_start = writer.calculate_coordinate(xleftside_addition)
else
  # add some space for the chromosome position label
  xleftside_addition = 0.014 * box_size * 1.5
  xside_end_addition = xside_end_addition + xleftside_addition + 0.014 * box_size
  x_start = writer.calculate_coordinate(xleftside_addition)
end

# read the feature file based on physical position
if(options.feature_file)
  filereader.read_feature_file(options.feature_file, chrom_feature_list)
end

yside = 0

xside = writer.calculate_plot_width(snp_list.get_num_included_snps, xside_end_addition)
xmax = writer.calculate_coordinate(xside)

# calculate starting position for each chromosomal feature
total_chrom_features = chrom_feature_list.features.length
ystart_chrom_features = Array.new

ymax = 0

total_chrom_features.times do |i|
  ystart_chrom_features << ymax
  yside = yside + box_size/4 * 2 * 0.01667
  ymax = writer.calculate_coordinate(yside)
end

y_label_start = ymax
yside = writer.add_space_for_labels(yside, 0.05)
ymax = writer.calculate_coordinate(yside)

# add space for snpinfo boxes
total_snpinfo_boxes = snp_list.get_total_snpinfo_boxes
ystart_snpinfo_boxes = Array.new

ystart_block_stat_box = ymax

total_snpinfo_boxes.times do |i|
  ystart_snpinfo_boxes << ymax  # + box_size/8
  yside = yside + box_size/4 * 1 * 0.01666#0.0125
  ymax = writer.calculate_coordinate(yside)
end

yside = yside + box_size/5 * 1 * 0.01667
ymax = writer.calculate_coordinate(yside)

# add space for stat boxes
total_stat_boxes = snp_list.get_total_stat_boxes
ystart_stat_boxes = Array.new


ystart_haplotypes = 0
yend_haplotypes = 0


# add space for the actual haplotypes
# area used depends on the block with
# the largest number of haplotypes
if snp_list.blocks.length > 0
if snp_list.blocks.values[0].haplotypes.length > 0
  max_haplotype = snp_list.get_max_haplotype_size
  yside = writer.add_haplotype_space(yside, max_haplotype)
  ystart_haplotypes = ymax
  ymax = writer.calculate_coordinate(yside)
  yend_haplotypes = ymax
end
end

total_stat_boxes.times do |i|
  ystart_stat_boxes << ymax + box_size/2
  yside = yside + box_size * 2.2 * 0.0125
  ymax = writer.calculate_coordinate(yside)
end

yend_block_stat_box = ymax

# add space for the bars to be drawn again above the d-prime
if total_chrom_features > 0
  y_lower_feature_box_start = ymax
  yside = yside + box_size/4 * 1 * 0.01667
  ymax = writer.calculate_coordinate(yside)
end

y_dprime_start = ymax

# get largest depth of the ld interactions
# use to set the size needed for the grids
# displaying LD results
maximum_interactions = snp_list.get_max_interactions + 1

# create a new image for displaying the data
# determine how much of the xside is used
if(options.dprime)
  yside = writer.add_grid_size(yside, xside, maximum_interactions, snp_list.get_num_included_snps)
  ymax = writer.calculate_coordinate(yside)
end
y_rsquared_start = ymax

# double all vertical sizes for second grid (r-squared)
if(options.rsquared)
  yside = writer.add_grid_size(yside, xside, maximum_interactions, snp_list.get_num_included_snps)
  ymax = writer.calculate_coordinate(yside)
end

rvg = RVG.new(xside.in, yside.in).viewbox(0,0,xmax,ymax) do |canvas|
  canvas.background_fill = 'rgb(253,253,253)'
  writer.canvas = canvas
  writer.box_size = box_size

  # add labels to top of chart
  writer.add_labels(snp_list, x_start, y_label_start, total_chrom_features > 0)

  # add snpinfo to the figure
  total_snpinfo_boxes.times do |current_track|
    writer.draw_snpinfo(current_track, snp_list, x_start, ystart_snpinfo_boxes[current_track])
  end

  # add any haplotypes to the figure
  if(ystart_haplotypes > 0)
    writer.draw_haplotypes(x_start, ystart_haplotypes, snp_list, options.analysis_file != nil)
  end

  # add stat boxes
 total_stat_boxes.times do |current_stat|
    writer.draw_stat_box(current_stat, snp_list, x_start, ystart_stat_boxes[current_stat], 0)
 end

 if total_stat_boxes > 0 || ystart_haplotypes > 0
   writer.draw_stat_box_blocks(x_start, ystart_block_stat_box, yend_block_stat_box, snp_list)
 end

  # draw d' grid
  if(options.dprime)
    writer.draw_grid(snp_list, x_start, y_dprime_start, true)
  end

  # draw r-squared grid
  if(options.rsquared)
    writer.draw_grid(snp_list, x_start, y_rsquared_start, false)
  end

  # draw features at top and add boxes to the 
  # bottom for feature selected
  total_chrom_features.times do |current_feature|
      if current_feature == total_chrom_features-1
        draw_box = true
      else
        draw_box = false
      end
      writer.draw_features(current_feature, chrom_feature_list, snp_list, x_start, ystart_chrom_features[current_feature], draw_box, y_label_start, y_dprime_start, y_lower_feature_box_start)
  end

end

outfile = options.out_name + '.' + options.imageformat
print "\n\tDrawing #{outfile}..."
STDOUT.flush
img = rvg.draw
img.write(outfile)
print " Created #{outfile}\n\n"