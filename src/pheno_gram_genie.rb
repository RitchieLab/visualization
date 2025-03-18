#!/usr/bin/env ruby
# Requres rubygems

ENV['MAGICK_CONFIGURE_PATH'] = '/gpfs/group/mdr23/usr/tools/etc/ImageMagick'

Font_family = 'Verdana'

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
  puts "\nPlease install RMagick -- See documentation for pheno_gram or http://rmagick.rubyforge.org/install-faq.html "
  puts
  exit(1)
end


require 'optparse'
require 'ostruct'
include Magick

Version = '0.1.0'
Name = 'pheno_gram_genie.rb'
#GeneFile = '/gpfs/group1/m/mdr23/software/visualization/data/knownGene.symbol.txt'
#GeneFile = '/Users/dudeksm/lab/visualization/test/pheno_genie/build38/knownGene.symbol.build38.txt'
#GeneFile = '/Users/dudeksm/lab/visualization/test/pheno_genie/knownGene.symbol.txt'
GeneFile = '/Users/dudeksm/lab/visualization/test/pheno_genie/build38/knownGene.symbol.build38.txt'

# check for windows and select alternate font based on OS
Font_family_style = RUBY_PLATFORM =~ /mswin/ ? "Verdana" : "Times"
Font_plot_family = RUBY_PLATFORM =~ /darwin/ ? "Verdana" : "Helvetica"

if RUBY_PLATFORM =~ /darwin/
  Font_phenotype_names = "Geneva"
elsif RUBY_PLATFORM =~ /mswin/
  Font_phenotype_names = "Verdana"
else
  Font_phenotype_names = "Helvetica"
end

#############################################################################
#
# Class Arg -- Parses the command-line arguments.
#
#############################################################################
class Arg

  def self.parse(args)
    options = OpenStruct.new
    options.vcf = nil
		options.pheno_inputfile = nil
		options.plot_trans_separately = false
		options.imageformat = 'png'
		options.out_name = 'genie'
		options.genename = nil
		options.exon_color = 'gray'
		options.sequence_color = 'blue'
		options.include_snps = false
		options.label=nil
		options.seq_id = nil
		help_selected = false
    version_selected = false

    
    opts = OptionParser.new do |opts|
      opts.banner = "Usage: #{Name}.rb [options]"
      opts.on("-e [vcf_file]", "VCF input file"){|input_file| options.vcf = input_file}
			opts.on("-i [pheno_input]", "Phenogram input file for plotting"){|pheno_input| options.pheno_inputfile=pheno_input}
			opts.on("-f [image_type]", "Image format for output (png default)."){|image_type| options.imageformat = image_type}
			opts.on("-g [gene_name]", "Gene name for plot"){|gname| options.genename = gname}
			opts.on("-s", "--transcript-separate", "Make a separate plot for each transcript"){|trans_plot|options.plot_trans_separately=true}
			opts.on("-E [exon_color]", "Color for exons"){|exon_color|options.exon_color=exon_color}
			opts.on("-S [sequence_color]", "Color for sequence blocks"){|sequence_color|options.sequence_color=sequence_color}
			opts.on("-r [rs_sequence]", "Specify ID for SNP or region sequence to include above plot"){|id|options.seq_id=id}
			opts.on("-P", "--include-single", "Include single base positions (SNPs) specified in input"){|inc_snps| options.include_snps=true}
			opts.on("-L [label]", "Specify identifier to label in plot along top"){|label| options.label=label}
      opts.on("-o [out_name]", "Optional output name for image"){|out_name| options.out_name = out_name}
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
    rescue => e
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
    
    if !options.genename and !options.vcf
      help_string = opts.help
      puts "\n",help_string,"\n"
      puts "\nExamples: #{Name} -i pheno_input.txt -o new_plot -t \"New Image\"\n"
      puts "          #{Name} -i pheno_input.txt -f jpg\n"
      print "\n"
      exit(1)
    end
    
    return options
    
  end
  
end

		@ucscidcol = 0
		@chrcol = 1
		@transcriptcol = 2
		@startcol = 3
		@finishcol = 4
		@exonnumcol = 5
		@exonstartcol = 6
		@exonendcol = 7
		@genesymcol = 12

class Region
	attr_accessor :chr, :start, :stop, :name
	
	def initialize
		name=""
	end
end

class Gene < Region
	attr_accessor :exons, :trans_direction
	
	def initialize
		@trans_direction = 1
		@exons = Array.new
	end
	
	# ex_st is array of start positions, ex_end is matching array of end positions
	def add_exons(ex_st, ex_end)
		ex_st.each_with_index do |st, i|
			@exons << Exon.new(st.to_i, ex_end[i].to_i)
		end
	end
end


class Exon
	attr_accessor :start, :stop
	
	def initialize(st,stop)
		@start = st
		@stop = stop
	end
	
end


class Genome
	
	def self.convert_chrom(chrstring)
		return 23 if chrstring =~ /^x/i
		return 24 if chrstring =~ /^y/i
		return chrstring.to_i
	end
	
	def self.format_chrom(chr)
		return 'X' if chr==23
		return 'Y' if chr==24
		return chr.to_s
	end
	
end

class FileHandler
  
  def close()
    @file.close
  end 
     
  # strips and splits the line
  # returns the data array that results
  def strip_and_split(line)
    line.rstrip!
    line.split(/\s/)
  end 
  
  def strip_and_split_delim(line,delim)
    line.rstrip!
    line.split(/#{delim}/)
  end 
  
end




class KnownGeneParser < FileHandler
	
	def initialize
		@ucscidcol = 0
		@chrcol = 1
		@transcriptcol = 2
		@startcol = 3
		@finishcol = 4
		@exonnumcol = 7
		@exonstartcol = 8
		@exonendcol = 9
		@genesymbcol = 12
	end

	# read known gene file and find specified gene
	def parse_file(params)
		filename = params[:filename] 
		gene = params[:genename] || nil
		
		genefile = File.open(filename)
		
		# skip header
		genefile.gets
		matching = Array.new
		while line = genefile.gets
			data = strip_and_split_delim(line, "\t")
			if gene == data[@genesymbcol]
				newgene = Gene.new
				
				data[@chrcol].gsub!(/chr/,'')
				newgene.chr = Genome.convert_chrom(data[@chrcol])
				newgene.name = data[@genesymbcol]
				newgene.start = data[@startcol].to_i
				newgene.stop = data[@finishcol].to_i
				newgene.trans_direction = data[@transcriptcol]
				exon_start = strip_and_split_delim(data[@exonstartcol], ",")
				exon_end = strip_and_split_delim(data[@exonendcol], ",")
				
				newgene.add_exons(exon_start, exon_end)
				matching << newgene
			end
		end
		
		genefile.close
		
		return matching
	end
end

class VCFParser < FileHandler
	
	def initialize
		# vcf column headers
		@vcfchr = '#CHROM'
		@vcfpos = 'POS'
		@vcfid = 'ID'
		@info = 'INFO'
	end
	
	def open(filename)
    @file = File.new(filename, "r")
  end
	
	# converts file from convert_file into format matching phenogram input
	# with start and stop positions specified
	def compress_file(params)
		open(params[:filename])
		# find headers and start of data
		
		headerline = @file.gets
		outfile = params[:outname] + '.phenogram.txt'
		
		# process data (merging any that fill multiple lines
		# output to new file
		# only write ending position when reach break in contiguous locations
		start=final=1
		fout = File.open(outfile, 'w')
		fout.puts "ID\tCHR\tPOS\tENDPOS"
		oline = @file.gets
		oline = @file.gets if oline =~ /^#/
		data = strip_and_split_delim(oline, "\t")
		start = data[2].to_i
		lastpos = start
		lastChr = data[1]
		lastID = data[0]
		
		while oline=@file.gets
			next if oline =~ /^#/ 
			data = strip_and_split_delim(oline, "\t")
			if data[0] !~ /^\./
				fout.puts "#{data[0]}\t#{data[1]}\t#{data[2]}\t#{data[2]}"
			end
			if data[2].to_i - 1 != lastpos
					fout.puts "#{lastID}\t#{lastChr}\t#{start}\t#{lastpos}"
					start = data[2].to_i
					lastChr = data[1]
					lastID = data[0]
			end
			lastpos = data[2].to_i			
		end
		
		fout.close
		close
	end
	
	# converts the vcf file to a format containing information 
	# needed (chrom, pos, and ID)
	def convert_file(params)
		open(params[:filename])
		# find headers and start of data
		while oline=@file.gets
			next if oline =~ /^##/
			parse_cols(oline) if oline =~ /^#/
			break
		end
		
		outfile = params[:outname] + '.phenogenie.txt'
		
		# process data (merging any that fill multiple lines
		# output to new file
		fout = File.open(outfile, 'w') 
		fout.puts "ID\tCHR\tPOS\tINFO"
		while oline=@file.gets
			next if oline =~ /^#/
			data = strip_and_split_delim(oline, "\t")
			fout.puts "#{data[@headers[@vcfid]]}\t#{data[@headers[@vcfchr]]}\t#{data[@headers[@vcfpos]]}\t#{data[@headers[@info]]}"
		end
		
		fout.close
		close
	end
	
	def parse_cols(headerline)
		hdrs = strip_and_split_delim(headerline, "\t")
		@headers = Hash.new
		hdrs.each_with_index {|h,i| @headers[h]=i}
	end
	
end

class Sequence
	attr_accessor :start, :stop, :id, :seq5, :seq3, :seqmark
	
	def initialize(start,stop,id)
		@start = start
		@stop = stop
		@id = id
	end
end


class SequenceHolder
	attr_accessor  :seqs

	def initialize
		@seqs = Array.new
	end
	
	def addseq(start,stop,id=nil)
		@seqs << Sequence.new(start,stop,id)
		return @seqs.last
	end
	
	def getseq(seq_id)
		@seqs.each do |seq| 
			return seq if seq.id == seq_id
		end
	end
	
end
class PhenoGramFileReader < FileHandler

  def open(filename)
    @file = File.new(filename, "r")
  end
  
  def set_columns(headerline)
    @idcol = @chromcol = @bpcol  = @bpendcol = @seq5col = @seq3col= @seqmarkcol = nil
    headers = strip_and_split_delim(headerline, "\t")
    
    headers.each_with_index do |header, i|
#      if header =~ /^snp$|^snp_id$|^ID$/i
			if header =~ /ID$|^SNP/i
        @idcol = i
      elsif header =~ /^chrom|^chr$/i
        @chromcol = i
      elsif header =~ /^bp|^pos|^start/i
        @bpcol = i
      elsif header =~ /^end/i
        @bpendcol = i
			elsif header =~ /^seq5/i
				@seq5col = i
			elsif header =~ /^seq3/i
				@seq3col = i
			elsif header =~ /^seqmark/i
				@seqmarkcol = i
      end
    end
    unless @chromcol and @bpcol
			error_string = 'Input file must include chrom, pos' 
			error_string = error_string + ' columns'
      raise error_string
    end
		
		# set end to be start if no ending bp column
		@bpendcol = @bpcol unless @bpendcol
  end
  
  def chr_good?(chrnum)
    return true if chrnum == "X" || chrnum == "Y"
    return true if chrnum.to_i >=1 and chrnum.to_i <=24
    return false
  end
  
  def parse_file(filename, sequences, params)
		genes = params[:genes]
		
		genestart=1e12
		genestop=1
		genechr=1
		genes.each do |gene|
			genestart = gene.start if gene.start < genestart
			genestop = gene.stop if gene.stop > genestop
			genechr = gene.chr
		end	
		
    open(filename)
		
		line = @file.gets
		set_columns(line) 
		while(line=@file.gets)
      next unless line =~ /\w/
      data = strip_and_split_delim(line,"\t")
      # add SNP info
			seqchr = Genome.convert_chrom(data[@chromcol])
			
      next if data[@chromcol] =~ /chrM/
	    raise "Problem in #{filename} with line:\n#{line}\n#{data[@chromcol]} is not a valid chromsome number" unless chr_good?(data[@chromcol])		
			
			seqstart = data[@bpcol].to_i
			seqend = data[@bpendcol].to_i
			# skip unless specified gene region
			next unless seqchr == genechr and (seqstart >= genestart and seqstart <= genestop) or 
				(seqend >= genestart and seqend <= genestop) or (seqstart <= genestart and seqend >= genestop)
      @idcol ? name = data[@idcol] : name = "."
			new_seq = sequences.addseq(seqstart, seqend, name)
			#:seq5, :seq3, :seqmark
			new_seq.seq5 = data[@seq5col] if @seq5col
			new_seq.seq3 = data[@seq3col] if @seq3col
			new_seq.seqmark = data[@seqmarkcol].to_i if @seqmarkcol
    end
		close
  end
end


class Scale
	
	def initialize(scale, range, min=0.0)
		@scale = scale.to_f
		@range = range.to_f
		@min = min
	end
	
	def get_coordinate(val)
		return (val-@min)/@range * @scale
	end
	
end



class Plotter
	
	def standard_font_size
		return 20
	end
	
	def inches_to_coords(val)
		return val * RVG::dpi/2
	end
	
	# return width (in pixels) of a string
	def get_text_width(text, params)
		pointsize = params[:pointsize] || 20
		font_family = params[:font_family] || nil
		font_style = params[:font_style] || nil
		density_fraction = RVG::dpi
		gc = Magick::Draw.new
		img = Magick::Image.new(900,50){
			self.density="#{density_fraction}x#{density_fraction}"
			self.format='png'
		}
		gc.pointsize = pointsize
		gc.font_family = font_family if font_family
		metrics = gc.get_type_metrics(img,text)
#		gc.text(0,20,text)
		
		return metrics.width
	end
	
	
	# add gray boxes with sequences contained for display
	def draw_sequence_boxes(params)
		xstart = params[:xstart] || 0
		ystart = params[:ystart] || 0
		height = params[:height]
		width = params[:width]
		canvas = params[:canvas]
		box_seq = params[:box_seq]
		xscaler = params[:xscale]
		
		font_size = standard_font_size * 1.2
		seq_mark = box_seq.seqmark - 1 # make zero based
		number_bases_max = 44
		
		# determine which position will hold the SNP bases
		snpx = xscaler.get_coordinate(box_seq.start)
		snp_index = (snpx/width * number_bases_max).to_i
		text_step_size = width / number_bases_max
		startbase = seq_mark - snp_index
		endbase = seq_mark + (number_bases_max - snp_index)-1
		start_index=0
		
		if startbase < 0
			start_index = -startbase
			startbase = 0
		end
		if endbase > box_seq.seq5.length
			endbase = box_seq.seq5.length -1
		end
		
		seq_mark = seq_mark + start_index if seq_mark
		top_line_y = height/4 * 2
		bottom_line_y = height
		seq5 = box_seq.seq5.split("")
		seq3 = box_seq.seq3.split("")

		gray_box_width = (endbase - startbase + 1) * text_step_size
		# draw gray boxes to highlight sequences
		top_shift = 3
		
		xtext = start_index * text_step_size
		canvas.g.translate(xstart, ystart) do |draw|
			draw.rect(gray_box_width, height/2-3, xtext-text_step_size/4,top_shift).styles(:fill=>'gray')
			draw.rect(gray_box_width, height/2-3, xtext-text_step_size/4,top_line_y+3).styles(:fill=>'gray')
		end
		
		font_weight=500
		fill = 'white'
		canvas.g.translate(xstart, ystart) do |draw|
			for i in startbase..endbase
				draw.text(xtext, top_line_y, seq5[i]).styles(:font_weight=>font_weight,
					:text_anchor=>'start', :fill=>fill, :font_size=>font_size)
				draw.text(xtext, bottom_line_y, seq3[i]).styles(:font_weight=>font_weight,
					:text_anchor=>'start', :fill=>fill, :font_size=>font_size)
				xtext += text_step_size
			end
		end
		
		box_center = (snp_index-1) * text_step_size + text_step_size/1.65
		canvas.g.translate(xstart,ystart) do |draw|
			draw.rect(1.4 * text_step_size, height+4, box_center, 0).round(5).styles(:fill=>'none', :stroke=>'red',
				:stroke_width=>2)
		end
	end
	
	def draw_sequences(params)
		xstart = params[:xstart] || 0
		ystart = params[:ystart] || 0
		xscaler = params[:xscale]
		height = params[:height]
		canvas = params[:canvas]
		xwidth = params[:width]
		sequences = params[:sequences]
		sequence_color = params[:sequence_color]
		include_snps = params[:include_snps] || false	
		label = params[:label] || nil
		
		box_height = 28
		box_y_start = height - box_height
		font_weight=700
		canvas.g.translate(xstart, ystart) do |draw|
			sequences.seqs.each do |seq|
				# plot regions
				if seq.start != seq.stop
					startx = xscaler.get_coordinate(seq.start)
					startx = 0 if startx < 0
					next if startx > xwidth
					endx = xscaler.get_coordinate(seq.stop)
					endx = xwidth if endx > xwidth
					draw.rect(endx-startx, box_height, startx, box_y_start).styles(:fill=>sequence_color)
				elsif include_snps
					# draw a line at location for SNP
					linex = xscaler.get_coordinate(seq.start)
					draw.line(linex, box_y_start, linex, box_y_start+box_height).styles(:stroke=>sequence_color)
					if seq.id == label
						draw.text(linex, box_y_start-standard_font_size/5, seq.id).styles(:font_weight=>font_weight,
				:text_anchor=>'middle', :fill=>'black', :font_size=>standard_font_size)
					end
				end
			end
		end
		
	end
	
	
	def draw_gene(params)
		gene = params[:gene]
		xstart = params[:xstart]
		ystart = params[:ystart]
		xscaler = params[:xscale]
		yheight = params[:height]
		canvas = params[:canvas]
		xwidth = params[:width]
		additional_text = params[:addtext]
		exon_color = params[:exon_color]
		
		line_y = yheight.to_f*1/4
		
		line_stroke_width=2
		font_weight=700

		line_color = 'lightslategray'
		canvas.g.translate(xstart,ystart) do |draw|
			draw.line(0,line_y,xwidth,line_y).styles(:stroke=>line_color,:stroke_width=>line_stroke_width)
			formatted_n = gene.start.to_s.reverse.gsub(/...(?=.)/,'\&,').reverse
			formatted_n = Genome.format_chrom(gene.chr) + ':' + formatted_n if formatted_n =~ /\d/
			draw.text(0,line_y+line_stroke_width*20,formatted_n).styles(:font_weight=>font_weight,:font_size=>standard_font_size/1.1,
						:text_anchor=>'start')
			formatted_n = gene.stop.to_s.reverse.gsub(/...(?=.)/,'\&,').reverse
			formatted_n = Genome.format_chrom(gene.chr) + ':' + formatted_n if formatted_n =~ /\d/
			draw.text(xwidth,line_y+line_stroke_width*20,formatted_n).styles(:font_weight=>font_weight,:font_size=>standard_font_size/1.1,
						:text_anchor=>'end')
			draw.text(xwidth/2, yheight-standard_font_size, gene.name).styles(:font_weight=>font_weight,
				:text_anchor=>'middle', :fill=>'black', :font_size=>standard_font_size*1.5)
			draw.text(xwidth/2, yheight-standard_font_size, additional_text).styles(:font_weight=>font_weight,
				:text_anchor=>'middle', :fill=>'black', :font_size=>standard_font_size)
		end
		
		
		arrow_color = 'lightslategray'
		# draw 10 arrows showing direction
		arrow_interval = xwidth.to_f/10
		start_arrow_x = arrow_interval
		gene.trans_direction == '+' ? arrow_adj=-line_stroke_width*3.0 : arrow_adj=line_stroke_width*3.0
		canvas.g.translate(xstart,ystart) do |draw|
				9.times do |i|
				# draw lines to make arrow
				end_arrow_x = start_arrow_x + arrow_adj
				draw.line(end_arrow_x,line_y-line_stroke_width*3, start_arrow_x, line_y).styles(:stroke=>arrow_color,
					:stroke_width=>line_stroke_width/2.5, :stroke_opacity=>0.5)
				draw.line(end_arrow_x,line_y+line_stroke_width*3, start_arrow_x, line_y).styles(:stroke=>arrow_color,
					:stroke_width=>line_stroke_width/2.5, :stroke_opacity=>0.5)
				start_arrow_x += line_stroke_width
				end_arrow_x= start_arrow_x + arrow_adj
				draw.line(end_arrow_x,line_y-line_stroke_width*3, start_arrow_x, line_y).styles(:stroke=>arrow_color,
					:stroke_width=>line_stroke_width/2.5, :stroke_opacity=>0.5)
				draw.line(end_arrow_x,line_y+line_stroke_width*3, start_arrow_x, line_y).styles(:stroke=>arrow_color,
					:stroke_width=>line_stroke_width/2.5, :stroke_opacity=>0.5)
				start_arrow_x += arrow_interval
			end
		end
		
		y_box_offset = line_stroke_width * 10
		box_color = 'darkslategray'
		canvas.g.translate(xstart,ystart) do |draw|
			# box for each exon
			gene.exons.each do |exon|
				startx = xscaler.get_coordinate(exon.start)
				endx = xscaler.get_coordinate(exon.stop)
				draw.rect(endx-startx,y_box_offset*2,startx, line_y-y_box_offset).styles(:fill=>exon_color)
			end
		end
		
	end
	
	def plot(options, params)
		
		gene = params[:gene]
		sequences = params[:sequences]
		transcriptnum = params[:transcriptnum]
		
		RVG::dpi = 300
		# width of plot will be constant (in inches)
		width = 8
		
		# determine height - dependent on which tracks included
		# always have the gene track
		height = 1
		
		current_start_y=0
		include_seq=false
		# add height for sequence data if included
		if options.seq_id
			box_seq = sequences.getseq(options.seq_id)
			if box_seq.start >= gene.start and box_seq.stop <= gene.stop
				seq_box_ystart = current_start_y
				seq_box_height= 0.32
				sequence_boxes_yheight = inches_to_coords(seq_box_height)
				height += seq_box_height
				current_start_y += sequence_boxes_yheight
				include_seq = true
			end
		end
		
		# add height for seq track
		seq_track_height =0.45
		seq_track_ystart = current_start_y
		seq_track_yheight = inches_to_coords(seq_track_height)
		height += seq_track_height
		
		# xmax is constant with each x being 2 pixels
		xmax = inches_to_coords(width)
		ymax = inches_to_coords(height)
		margin = inches_to_coords(0.25)
		gene_y_start = seq_track_ystart + seq_track_yheight
		gene_y_height = ymax - gene_y_start
		draw_width = xmax - margin*2
		rvg=RVG.new(width.in, height.in).viewbox(0,0,xmax,ymax) do |canvas|
			canvas.background_fill = 'rgb(255,255,255)'
			x_scale = Scale.new(xmax - margin*2, gene.stop-gene.start, gene.start)
			if options.seq_id and include_seq
				box_seq = sequences.getseq(options.seq_id)
				draw_sequence_boxes(:xstart=>margin, :ystart=>seq_box_ystart, 
					:canvas=>canvas, :width=>draw_width, :height=>sequence_boxes_yheight,
					:box_seq=>box_seq, :xscale=>x_scale) 
			end			
			draw_sequences(:xstart=>margin, :ystart=>seq_track_ystart, :xscale=>x_scale,
				:height=>seq_track_yheight, :canvas=>canvas, :width=>draw_width,
				:sequences=>sequences, :sequence_color=>options.sequence_color, 
				:label=>options.label, :include_snps=>options.include_snps)
			
			draw_gene(:gene=>gene, :xstart=>margin, :ystart=>gene_y_start,
				:xscale=>x_scale, :height=>gene_y_height, :canvas=>canvas,
				:width=>draw_width,:exon_color=>options.exon_color)# :addtext=>"Transcript #{transcriptnum}")
			
		end
		

		outfile = options.out_name + "-#{transcriptnum}" + '.' + options.imageformat
		create_file(outfile, rvg)
	end
	
	def create_file(outfile, rvg)
		print "\n\tDrawing #{outfile}..."
		STDOUT.flush
		img = rvg.draw
		img.write(outfile)
		print " Created #{outfile}\n\n" 		
	end
	
	def plotall(options, params)
		
		genes = params[:genes]
		sequences = params[:sequences]
		
		RVG::dpi = 300
		# width of plot will be constant (in inches)
		width = 8
		
		# determine height - dependent on which tracks included
		# always have the gene track
		each_gene_track = 0.32
		height = each_gene_track * genes.length + 0.3 # add 0.5 for bottom gene name
		
		current_start_y=0
		# add height for sequence data if included
		if options.seq_id
			seq_box_ystart = current_start_y
			seq_box_height= 0.32
			sequence_boxes_yheight = inches_to_coords(seq_box_height)
			height += seq_box_height
			current_start_y += sequence_boxes_yheight
		end
		
		# add height for seq track
		seq_track_height =0.45
		seq_track_ystart = current_start_y
		seq_track_yheight = inches_to_coords(seq_track_height)
		height += seq_track_height
		
		buffer = 0.12
		ybuffer = inches_to_coords(buffer)
		height += buffer
		
		# xmax is constant with each x being 2 pixels
		xmax = inches_to_coords(width)
		ymax = inches_to_coords(height)
		margin = inches_to_coords(0.25)

		gene_y_start = seq_track_ystart + seq_track_yheight + ybuffer
#		gene_y_height = ymax - gene_y_start
		gene_y_height = inches_to_coords(each_gene_track)
		draw_width = xmax - margin*2
		
		genestart=1e12
		genestop=1
		genes.each do |gene|
			genestart = gene.start if gene.start < genestart
			genestop = gene.stop if gene.stop > genestop
		end	
	
		genename = genes.first.name
		
		rvg=RVG.new(width.in, height.in).viewbox(0,0,xmax,ymax) do |canvas|
			canvas.background_fill = 'rgb(255,255,255)'
			x_scale = Scale.new(xmax - margin*2, genestop-genestart, genestart)
			
			if options.seq_id
				box_seq = sequences.getseq(options.seq_id)
				draw_sequence_boxes(:xstart=>margin, :ystart=>seq_box_ystart, 
					:canvas=>canvas, :width=>draw_width, :height=>sequence_boxes_yheight,
					:box_seq=>box_seq, :xscale=>x_scale) 
			end
			
			draw_sequences(:xstart=>margin, :ystart=>seq_track_ystart, :xscale=>x_scale,
				:height=>seq_track_yheight, :canvas=>canvas, :width=>draw_width,
				:sequences=>sequences, :sequence_color=>options.sequence_color,
				:include_snps=>options.include_snps, :label=>options.label)
			genes.each_with_index do |gene,i|
				gene.name = ' '
				if i == genes.length-1 
					gene.start = genestart
					gene.stop = genestop
				else
					gene.start = gene.stop = " "
				end

				draw_gene(:gene=>gene, :xstart=>margin, :ystart=>gene_y_start,
					:xscale=>x_scale, :height=>gene_y_height, :canvas=>canvas,
					:width=>draw_width, :exon_color=>options.exon_color)
				gene_y_start += gene_y_height
			end
			
			canvas.g do |draw|
				draw.text(xmax/2, ymax-standard_font_size, genename).styles(:font_weight=>700,
				:text_anchor=>'middle', :fill=>'black', :font_size=>standard_font_size*1.5)
			end
		end
	
		outfile = options.out_name + '.' + options.imageformat
		create_file(outfile, rvg)
	end
end


options = Arg.parse(ARGV)

if options.vcf
	vcfparse=VCFParser.new
	vcfparse.convert_file(:filename=>options.vcf, :outname=>options.out_name)
	#vcfparse.compress_file(:filename=>options.vcf, :outname=>options.out_name)
else

	geneparse = KnownGeneParser.new
	genes = geneparse.parse_file(:filename => GeneFile, :genename=>options.genename)
	sequences = SequenceHolder.new
	phenoparser = PhenoGramFileReader.new
	phenoparser.parse_file(options.pheno_inputfile, sequences, :genes=>genes)

	plotter = Plotter.new

	if options.plot_trans_separately
		genes.each_with_index do |gene,i|
			plotter.plot(options, :gene=>gene, :sequences=>sequences, :transcriptnum=>i+1)
		end
	else
		plotter.plotall(options, :genes=>genes, :sequences=>sequences)
	end
end
