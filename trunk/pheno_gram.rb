#!/usr/bin/env ruby
# Requres rubygems
begin
  require 'rubygems'
rescue Exception => e
  puts e
  puts "Please install rubygems -- http://docs.rubygems.org/read/chapter/3 "
  exit(1)
end

# Requires RMagick (http://rmagick.rubyforge.org/)
begin
  require 'rvg/rvg'
  rescue Exception => e
  puts
  puts e
  puts "\nPlease install RMagick -- See documentation for pheno_gram or http://rmagick.rubyforge.org/install-faq.html "
  puts
  exit(1)
end

#require '/Users/dudeksm/Documents/lab/scripts/pleo_view/PleoView/lib/color_gen'
# adapted from http://www.codeproject.com/Tips/258405/Random-Color-Generator#alternative3
# used_colors should have key be a RGB object
module ColorGen
  class ColorGenerator
    attr_accessor :accuracy
   
    MaximumDistanceMin = 1.0
    MinimumDistanceMin = 0.01
    
    def initialize
      @used_colors = Hash.new
      # 200 is good hi diff but low speed
      # 50 is good for hi speed but less difference in colors
      @accuracy = 100
      distance_min(0.5)
      @bad_try_count = 0
    end
    
    def add_used_color(color, ratio)
      @used_colors[color] = ColorRatio.new(color, ratio)
    end
    
    def used_color_hash_reset
      distance_min(MaximumDistanceMin)
    end
    
    def used_color_hash_remove
      if @used_colors.length <= 2
        distance_min(distance_min+0.1)
      else
        distance_min(MaximumDistanceMin)
      end
    end
    
    def remove_color
      @used_colors.delete(color)
      @used_color_hash_remove
    end
    
    def get_distance_min
      return @dist_min
    end
    
    def distance_min(value)
      if(value > MaximumDistanceMin)
        @dist_min = MaximumDistanceMin
      elsif(value < MinimumDistanceMin)
        @dist_min = MinimumDistanceMin
      else
        @dist_min = value
      end
    end
    
    def get_next_color
      while(true)
        hue = rand * 360
        saturation = rand
        luminance = rand
        
        hsl = HSL.new(hue, saturation, luminance)
        col = hsl.to_color
        
        if is_far_from_existing_color(col, get_distance_min)
          @used_colors[col] = ColorRatio.new(col)
          distance_min(get_distance_min+0.02)
          @bad_try_count = 0
          return col
        end
        
        @bad_try_count+=1
        if(@bad_try_count > @accuracy)
          @bad_try_count = 0
          distance_min(get_distance_min - 0.002)
        end
      end
      
    end
    
    def is_far_from_existing_color(c, dist_min)
      @used_colors.each_pair do |color, ratio|
        distance = ColorSpaceHelper.get_color_distance_cie_lab(c, color) /100
        if (distance/ratio.keep_away_ratio < get_distance_min)
          return false
        end
      end
      return true
    end
    
  end
  
  class HSL
    attr_accessor :hue, :saturation, :luminance
    
    def initialize(hue, sat, lum)
      @hue = hue
      @saturation = sat
      @luminance = lum
    end
    
    def to_color
      return ColorSpaceHelper.hsl_to_rgb(@hue, @saturation, @luminance)
    end
    
  end
  
  
  class CIEXYZ
    attr_accessor :x,:y,:z
    
    def initialize(x,y,z)
      @x = (x>0.9505)? 0.9505 : ((x<0)? 0 : x)
			@y = (y>1.0)? 1.0 : ((y<0)? 0 : y)
			@z = (z>1.089)? 1.089 : ((z<0)? 0 : z)
    end
    
  end
  
    # white structure
  D65 = CIEXYZ.new(0.9505, 1.0, 1.0890)
  
  class CIELab
    attr_accessor :l,:a,:b
    
    def self.get_distance_between_cie2000(lab1, lab2)
      p25 = 25**7
      
      c1 = Math.sqrt(lab1.a * lab1.a + lab1.b * lab1.b)
      c2 = Math.sqrt(lab2.a * lab2.a + lab2.b * lab2.b)
      avgc = (c1+c2)/2.0
      
      powavgc = avgc ** 7
      g = (1-Math.sqrt(powavgc/(powavgc+p25)))/2.0
      
      a_1 = lab1.a * (1+g)
      a_2 = lab2.a * (1+g)
      
      c_1 = Math.sqrt(a_1 * a_1 + lab1.b * lab1.b)
      c_2 = Math.sqrt(a_2 * a_2 + lab2.b * lab2.b)
      
      avgc_ = (c_1+c_2.to_f)/2.0
      
      h1 = (Math.atan2(lab1.b, a_1) >= 0 ? Math.atan2(lab1.b, a_1) : Math.atan2(lab1.b, a_1) + 360.0)
      h2 = (Math.atan2(lab2.b, a_1) >= 0 ? Math.atan2(lab2.b, a_1) : Math.atan2(lab2.b, a_1) + 360.0)
      
      h = (h1 - h2 > 180.0 ? (h1 + h2 + 360.0) / 2.0 : (h1 + h2) / 2.0)
      
      t = 1
      t -= 0.17 * Math.cos(h - 30)
      t += 0.24 * Math.cos(2 * h)
      t += 0.32 * Math.cos(3 * h + 6)
      t -= 0.20 * Math.cos(4 * h - 63)    
      
      deltah = 0
      if (h2 - h1 <= 180)
        deltah = h2 - h1
      elsif (h2 <= h1)
        deltah = h2 - h1 + 360.0
      else
        deltah = h2 - h1 - 360.0    
      end
      
      avgl = (lab1.l + lab2.l) / 2.0
      deltal_ = lab2.l-lab1.l
      deltac_ = c_2-c_1
      deltah_ = 2 * Math.sqrt(c_1*c_2)*Math.sin(deltah/2.0)
      
      
      sl = 1 + (0.015 * ((avgl - 50) ** 2)) / Math.sqrt(20 + ((avgl - 50) ** 2))
      sc = 1 + 0.045 * avgc_
      sh = 1 + 0.015 * avgc_ * t
      
      exp = ((h - 275) / 25) ** 2;
      teta = 30 ** -exp
      
      rc = 2.0 * Math.sqrt((avgc_** 7) / ((avgc_ ** 7) + p25))
      rt = -rc * Math.sin(2*teta)
      
      deltae = 0
      deltae = (deltal_ / sl) ** 2.0
      deltae += (deltac_ / sc) ** 2.0
      deltae += (deltah_ / sh) ** 2.0
      deltae += rt * (deltac_ / sc) * (deltah_ / sh)
      deltae = Math.sqrt(deltae)
      return deltae
    end
    
  end
  
  
  class RGB
    
    attr_accessor :red, :green, :blue
    
    def initialize(r,g,b)
      @red = r
      @green = g
      @blue = b
    end
    
    def to_s
      return "rgb(#{@red},#{@green},#{@blue})"
    end
    
  end
  
  class ColorSpaceHelper
    
    def self.hsl_to_rgb(h,s,l)
      
      if(s==0)
        val = (l * 255).to_i
        return RGB.new(val,val,val)
      else
        
        q = (l<0.5)?(l * (1.0+s)):(l+s - (l*s))
        p = (2.0 * l) - q
        
        hk = h.to_f/360.0
        t = Array.new
        t[0] = hk + (1.0/3.0)
        t[1] = hk
        t[2] = hk - (1.0/3.0)
        
        3.times do |i|
					 t[i] += 1.0 if(t[i] < 0)
					 t[i] -= 1.0 if(t[i] > 1)

					if((t[i]*6) < 1)
						t[i] = p + ((q-p)*6.0*t[i])
					elsif((t[i]*2.0) < 1) #//(1.0/6.0)<=t[i] && t[i]<0.5
						t[i] = q;
					elsif((t[i]*3.0) < 2) #// 0.5<=T[i] && T[i]<(2.0/3.0)
						t[i] = p + (q-p) * ((2.0/3.0) - t[i]) * 6.0
					else 
            t[i] = p
          end
        end
      end

      return RGB.new((t[0].to_f*255).to_i,(t[1].to_f*255).to_i,(t[2].to_f*255).to_i)      
    end
    
    def self.get_color_distance_cie_lab(c1,c2)
      return CIELab.get_distance_between_cie2000(rgbtolab(c1), rgbtolab(c2))
    end
    
    def self.rgbtolab(col)
      return xyztolab(rgbtoxyz(col.red, col.green, col.blue))
    end
    
    def self.rgbtoxyz(red,green,blue)
      rLinear = red.to_f/255
      gLinear = green.to_f/255
      bLinear = blue.to_f/255
      
      r = (rLinear > 0.4045)?((rLinear+0.055)/(1+0.055))**2.4 : rLinear/12.92
      g = (gLinear > 0.4045)?((gLinear+0.055)/(1+0.055))**2.4 : gLinear/12.92
      b = (bLinear > 0.4045)?((bLinear+0.055)/(1+0.055))**2.4 : bLinear/12.92
      
      return CIEXYZ.new((r * 0.4124564 + g * 0.3575761 + b * 0.1804375),
                (r * 0.2126729 + g * 0.7151522 + b * 0.0721750),
                (r * 0.0193339 + g * 0.1191920 + b * 0.9503041))
    end
      
    #def self.xyztolab(x,y,z)
    def self.xyztolab(xyz)
      x= xyz.x
      y= xyz.y
      z= xyz.z
      
      lab = CIELab.new
      
      # uses D65 (white)
      lab.l = 116.0 * fxyz( y/D65.y ) -16;
      lab.a = 500.0 * (fxyz( x/D65.x ) - fxyz( y/D65.y) );
      lab.b = 200.0 * (fxyz( y/D65.y ) - fxyz( z/D65.z) );
      return lab      
      
    end
    
    def self.fxyz(t)
      return ((t > 0.008856)? t**(1.0/3.0) : (7.787*t + 16.0/116.0))
    end
    
    
  end
  
  class ColorRatio
    attr_accessor :color, :keep_away_ratio
    
    def initialize(color, keepAwayRatio=1.0)
      @color = color
      @keep_away_ratio = keepAwayRatio
    end
    
  end
  
end



require 'optparse'
require 'ostruct'
include Magick

RVG::dpi=600

Version = '0.3.0'
Name = 'pheno_gram.rb'

# check for windows and select alternate font based on OS
Font_family_style = RUBY_PLATFORM =~ /mswin/ ? "Verdana" : "Times"
Font_plot_family = RUBY_PLATFORM =~ /darwin/ ? "Arial" : "Helvetica"

if RUBY_PLATFORM =~ /darwin/
  Font_phenotype_names = "Geneva"
elsif RUBY_PLATFORM =~ /mswin/
  Font_phenotype_names = "Arial"
else
  Font_phenotype_names = "Helvetica"
end

ChromStyle = {:stroke=>'darkgray',:stroke_width=>2, :fill=>'none'}

#############################################################################
#
# Class Arg -- Parses the command-line arguments.
#
#############################################################################
class Arg

  def self.parse(args)
    options = OpenStruct.new
    options.input = nil
    options.out_name = 'pheno_gram'
    options.imageformat = 'png'
    options.title = "Title Here"
    options.color = 'random'
    options.pheno_spacing = 'alternative'
    options.rand_seed = 7
    help_selected = false
    version_selected = false
    
    opts = OptionParser.new do |opts|
       opts.banner = "Usage: #{Name}.rb [options]"
      opts.on("-i [input_file]", "Input file") do |input_file|
        options.input = input_file
      end
      opts.on("-o [out_name]", "Optional output name for image") do |out_name|
        options.out_name = out_name
      end
      opts.on("-t [title]", "Main title for plot (enclose in quotes)") do |title|
        options.title = title
      end
      opts.on("-c [color_range]", "Options are random (default), web, or generator") do |color_range|
        options.color = color_range
      end
      opts.on("-p [pheno_spacing]", "Options are standard or equal or alternative (default) ") do |pheno_spacing|
        options.pheno_spacing = pheno_spacing
      end
      opts.on("-r [random_seed]", "Random number generator seed (default=7") do |random_seed|
        options.rand_seed = random_seed.to_i
      end
      opts.on("-f [image_type]", "Image format for output (png default).  Other options depend on ImageMagick installation.") do |image_type|
        options.imageformat = image_type
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
    rescue Exception => e
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
    
    if !options.input
      help_string = opts.help
      puts "\n",help_string,"\n"
      puts "\nExamples: #{Name}.rb -i pheno_input.txt -o new_plot -t \"New Image\"\n"
      puts "          #{Name}.rb -i pheno_input.txt -f jpg\n"
      print "\n"
      exit(1)
    end
    
    return options
    
  end
  
end


class Phenotype
  attr_accessor :name, :color
  
  def initialize(n,c)
    @name = n
    @color = c
  end
  
end


class PhenotypeHolder
  attr_accessor :phenonames
  
  def initialize(params)
    
    if params[:color]=='web'
      @colormaker = WebColorMaker.new
    elsif params[:color]=='generator'
      @colormaker = ColorGenColorMaker.new
    else
      @colormaker = ColorMaker.new
    end
    @phenonames = Hash.new
  end
  
  def add_phenotype(name)
    unless @phenonames.has_key?(name)
      pheno = Phenotype.new(name, set_color)
      @phenonames[pheno.name] = pheno
    end
    return @phenonames[name]
  end
  
  def get_phenotype(name)
    return @phenonames[name]
  end
  
  def set_color
    return @colormaker.gen_html
  end
  
end

class Genome
  attr_accessor :chromosomes
  
  def initialize
    @chromosomes = Array.new
    (1..24).each {|i| @chromosomes[i]=Chromosome.new(i)}
  end
  
  def add_snp(params)
    @chromosomes[params[:chr].to_i].add_snp(:name=>params[:name], :pos=>params[:pos],
      :pheno=>params[:pheno], :chr=>params[:chr])
  end
  
  def snps_per_chrom
    (1..24).each do |i|
      puts "chrom #{i} has #{@chromosomes[i].snps.length} SNPs"
    end
  end
  
end

class Chromosome
  attr_accessor :number, :snps, :snpnames, :centromere, :size, :display_num
  @@centromeres = Array.new
  @@centromeres << 0
  @@centromeres << 124496354 #Array.new(121535434, 124535434) 
  @@centromeres << 92893890 #Array.new(92326171, 95326171)  
  @@centromeres << 90566355 #Array.new(90504854, 93504854)  
  @@centromeres << 50126782 #Array.new(49660117, 52660117)
  @@centromeres << 48648211 #Array.new(46405641, 49405641)
  @@centromeres << 60807186 #Array.new(58830166, 61830166)
  @@centromeres << 59726986 #Array.new(58054331, 61054331)
  @@centromeres << 45521049 #Array.new(43838887, 46838887)
  @@centromeres << 48615828 #Array.new(47367679, 50367679)
  @@centromeres << 40049365 #Array.new(39254935, 42254935)
  @@centromeres << 53363994 #Array.new(51644205, 54644205)
  @@centromeres << 35693838 #Array.new(34856694, 37856694)
  @@centromeres << 17688322 #Array.new(16000000, 19000000)
  @@centromeres << 17201315 #Array.new(16000000, 19000000)
  @@centromeres << 18930038 #Array.new(17000000, 20000000)
  @@centromeres << 36515875 #Array.new(35335801, 38335801)
  @@centromeres << 23992494 #Array.new(22263006, 25263006)
  @@centromeres << 17185073 #Array.new(15460898, 18460898)
  @@centromeres << 26387347 #Array.new(24681782, 27681782)
  @@centromeres << 27427515 #Array.new(26369569, 29369569)
  @@centromeres << 13154793 #Array.new(11288129, 14288129)
  @@centromeres << 14591276 #Array.new(13000000, 16000000)
  @@centromeres << 60340916 #Array.new(58632012, 61632012)
  @@centromeres << 12541733 #Array.new(10104553, 13104553)
  
  # sizes taken from ensembl
  @@chromsize = Array.new
  @@chromsize << 0
  @@chromsize << 249250621
  @@chromsize << 243199373
  @@chromsize << 198022430
  @@chromsize << 191154276
  @@chromsize << 180915260
  @@chromsize << 171115067
  @@chromsize << 159138663
  @@chromsize << 146364022
  @@chromsize << 141213431
  @@chromsize << 135534747
  @@chromsize << 135006516
  @@chromsize << 133851895
  @@chromsize << 115169878
  @@chromsize << 107349540
  @@chromsize << 102531392
  @@chromsize << 90354753
  @@chromsize << 81195210
  @@chromsize << 78077248
  @@chromsize << 59128983
  @@chromsize << 63025520
  @@chromsize << 48129895
  @@chromsize << 51304566
  @@chromsize << 155270560
  @@chromsize << 59373566
  
  def self.chromsize(n)
    return @@chromsize[n]
  end
  
  def initialize(n)
    @number = n
    if n < 23
      @display_num = n
    elsif n == 23
      @display_num = 'X'
    elsif n == 24
      @display_num = 'Y'
    end

    @snps = Hash.new
    @snpnames = Array.new
    @centromere = @@centromeres[n]
    @size = @@chromsize[n]
  end
  
  def add_snp(params)
    unless snp = @snps[params[:name]]
      snp = @snps[params[:name]] = SNP.new(params[:chr], params[:pos])
      @snpnames << params[:name]
    end
    snp.phenos << params[:pheno]
  end
  
  def sort_snps!
    @snpnames.sort!{|x,y| @snps[x].pos.to_i <=> @snps[y].pos.to_i}
  end
  
end

class SNP
  attr_accessor :chrom, :pos, :phenos
  
  def initialize(c,p)
    @chrom = c
    @pos = p
    @phenos = Array.new
  end
  
end


class ColorMaker
	
  def initialize
    @golden_rule_conjugate = 0.618033988749895
    @h = rand
  end
  
# HSV values in [0..1[
# returns [r, g, b] values from 0 to 255
# adapted from http://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
  def hsv_to_rgb(h, s, v)
    h_i = (h*6).to_i
    f = h*6 - h_i
    p = v * (1 - s)
    q = v * (1 - f*s)
    t = v * (1 - (1 - f) * s)
    r, g, b = v, t, p if h_i==0
    r, g, b = q, v, p if h_i==1
    r, g, b = p, v, t if h_i==2
    r, g, b = p, q, v if h_i==3
    r, g, b = t, p, v if h_i==4
    r, g, b = v, p, q if h_i==5
    color_array=[(r*256).to_i, (g*256).to_i, (b*256).to_i]
  end
  
  def gen_html 
    @h += @golden_rule_conjugate
    @h %= 1
    color_array = hsv_to_rgb(@h, 0.8, 0.95)
    return "rgb(#{color_array[0]},#{color_array[1]},#{color_array[2]})"
  end
  
end


class ColorGenColorMaker
  
  def initialize
    @generator = ColorGen::ColorGenerator.new
    # add colors to avoid (like black and white)
    black = ColorGen::RGB.new(0,0,0)
    white = ColorGen::RGB.new(255,255,255)
    @generator.add_used_color(black, 2.0)
    @generator.add_used_color(white, 1.2)
    
    # add basic colors to start
    @colors = Array.new
    blue = ColorGen::RGB.new(0,0,255)
    red = ColorGen::RGB.new(255,0,0)
    green = ColorGen::RGB.new(0,255,0)
    yellow = ColorGen::RGB.new(255,255,0)
    orange = ColorGen::RGB.new(255,165,0)
    purple = ColorGen::RGB.new(160,32,240)
    pink = ColorGen::RGB.new(255,192,203)
    brown = ColorGen::RGB.new(165,42,42)
    gray = ColorGen::RGB.new(190,190,190)
    black = ColorGen::RGB.new(0,0,0)
    
    @colors << blue
    @colors << red
    @colors << green
    @colors << yellow
    @colors << purple
    @colors << orange
    @colors << pink
    @colors << brown
    @colors << gray
    @colors << black
    
    @generator.add_used_color(blue, 1.0)
    @generator.add_used_color(red, 1.0)
    @generator.add_used_color(green, 1.0)
    @generator.add_used_color(yellow, 1.0)
    @generator.add_used_color(purple, 1.0)
    @generator.add_used_color(orange, 1.0)
    @generator.add_used_color(pink, 1.0)
    @generator.add_used_color(brown, 1.0)
    @generator.add_used_color(gray, 2.0)
    
    @colors_used = 0
    
  end
  
  def gen_html
    if @colors_used < @colors.length
      color = @colors[@colors_used]
      @colors_used += 1
    else
      color = @generator.get_next_color
    end
    return color.to_s
  end
  
end



class WebColorMaker
  
  def initialize
    @outer_loop=0
    @inner_loop=0
    @colors = Array.new
    
    @inner_loop_array  = [2,6,4,8,3,5,7,1,0]
    
    hex = ["00", "33", "66", "99", "CC", "FF"]
    hex.each do |red|
      hex.each do |green|
        hex.each do |blue|
          @colors << "##{red}#{green}#{blue}"
        end
      end
    end
    @colors.reverse!
  end
  
  
  def gen_html
    index = @outer_loop * 9 + @inner_loop_array[@inner_loop]
    @outer_loop+=1
    if @outer_loop == 24
      @inner_loop+=1
      @outer_loop=0
    end
    return @colors[index]
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
  
end



class PhenoGramFileReader < FileHandler

  def open(filename)
    @file = File.new(filename, "r")
  end
  
  def set_columns(headerline)
    @snpcol = @chromcol = @bpcol = @phenocol = nil
    headers = strip_and_split(headerline)
    
    headers.each_with_index do |header, i|
      if header =~ /^snp$/i
        @snpcol = i
      elsif header =~ /^chrom|^chr$/i
        @chromcol = i
      elsif header =~ /^bp|^pos/i
        @bpcol = i
      elsif header =~ /^pheno/i
        @phenocol = i
      end
    end
    
    unless @snpcol and @chromcol and @bpcol and @phenocol
      raise 'Input file must include SNP, chrom, pos, phenotype columns'
    end
    
  end
  
  def parse_file(filename, genome, phenoholder)
    open(filename)
    set_columns(@file.gets)
    
    while line=@file.gets
      data = strip_and_split(line)
      # add SNP info 
      pheno = phenoholder.add_phenotype(data[@phenocol])

      genome.add_snp(:name => data[@snpcol], :chr=>data[@chromcol], :pos=>data[@bpcol],
        :pheno=>pheno)
    end  
    close   
  end 
end

class Plotter
  @@circle_size=0
  @@maxchrom=0
  
  def self.set_circle(n)
    @@circle_size = n
  end
  
  def self.set_maxchrom(n)
    @@maxchrom=n
  end
  
end


class Circle
  attr_accessor :x, :y, :color
  def initialize(x,y,col)
    @x = x
    @y = y
    @color = color
  end
end

class PhenoBox < Plotter
  attr_accessor :top_y, :bottom_y, :circles, :phenocolors, :chrom_y, :height, :up
  @@circles_per_row = 0
  
  def self.set_circles_per_row(c)
    @@circles_per_row = c
  end
  
  def initialize
    @circles = Array.new
    @phenocolors = Array.new
    @up = true
  end
  
  def add_circle(x,y,color)
    @circles << Circle.new(x,y,color)
  end
  
  def add_phenocolor(p)
    @phenocolors << p
  end
  
  def estimate_height
    return (@phenocolors.length.to_f/@@circles_per_row).ceil * @@circle_size
  end
  
  def set_default_boundaries(center)
    @chrom_y = center
    
    if @up
      adjust = -@@circle_size
    else
      adjust = @@circle_size
    end
    
    @top_y = center+adjust.to_f/4
    if @phenocolors.length <= @@circles_per_row
      @bottom_y = @top_y + @@circle_size
    else # each row will overlap by 1/4 of a circle on one above
      @bottom_y = @top_y + @@circle_size * 0.75 * @phenocolors.length/@@circles_per_row
    end
    @height = @bottom_y - @top_y 
  end

  def set_even_boundaries(chrom_y, y)
    @chrom_y = chrom_y
    @top_y = y
    if @phenocolors.length <= @@circles_per_row
      @bottom_y = @top_y + @@circle_size
    else # each row will overlap by 1/4 of a circle on one above
      @bottom_y = @top_y + @@circle_size * 0.75 * @phenocolors.length/@@circles_per_row
    end   
    @height = @bottom_y - @top_y
  end
  
  # set the top and bottom location
  # top is simply halfway above first row
  def set_boundaries(center)
    @chrom_y = center
    @top_y = center#+adjust.to_f/4
    if @phenocolors.length <= @@circles_per_row
      @bottom_y = @top_y + @@circle_size
    else # each row will overlap by 1/4 of a circle on one above
      @bottom_y = @top_y + @@circle_size * 0.75 * @phenocolors.length/@@circles_per_row
    end
    @height = @bottom_y - @top_y
  end
  
end


class PhenoBin
  attr_accessor :actual_height,:height_needed, :starty, :endy, :boxes,
    :startbase, :endbase, :boxpos
  
  def initialize
    @actual_height=@height_needed=@starty=@endy=0
    @boxes = Array.new
    @boxpos = Hash.new
  end
  
  def calc_height_needed
    @height_needed=0
    @boxpos.each_value{|box| @height_needed+=box.height}
    return @height_needed
  end
  
  def estimate_height
    estimate_total=0
    @boxpos.each_value{|box| estimate_total+=box.estimate_height}
    return estimate_total  
  end
  
  def sort_boxes!
    @boxes.sort!{|x,y| x.to_i <=> y.to_i} 
  end
  
  def height_discrepancy
    return @actual_height-@height_needed
  end
  
  def add_phenotype_snp(pos, col)
    if @boxpos.has_key?(pos)
      pbox = @boxpos[pos]
    else
      pbox = PhenoBox.new
      @boxpos[pos]=pbox
      @boxes << pos
    end
    
    pbox.add_phenocolor(col)
  end
  
  def bp_from_y(y)
    return y/@actual_height.to_f * (@endbase-@startbase) + @startbase
  end
  
  def y_from_bp(bp)
    return (bp.to_f-@startbase)/(@endbase-@starbase) * @actual_height
  end
  
  def set_phenobox_y
    @height_needed=0
    @boxpos.each_pair do |pos, phenobox|
      ycenter = (pos.to_i-@startbase)/(@endbase-@startbase).to_f * @actual_height
      @height_needed += phenobox.set_boundaries(ycenter)
    end
  end

  
end

# contains the phenoboxes
class PhenoBinHolder
  attr_accessor :phenobins, :totalchromy, :totalbases, :totaly
  
  def initialize
    @phenobins = Array.new
  end
 
  def set_num_bins(n)
    n.times do |i|
     @phenobins << PhenoBin.new
     @phenobins.last.actual_height = totaly.to_f/n
    end
  end
  
  # set base intervals for the bins
  def set_bases
    height=0
    y=0
    
    y_interval = @totaly.to_f/@phenobins.length
    curr_y = 0
    base_interval = @totalbases.to_f/@phenobins.length
    currbase=0   
    @phenobins.each do |pb|
      pb.starty = curr_y
      pb.startbase = currbase
      currbase += base_interval
      curr_y += y_interval
      pb.endy = curr_y
      pb.endbase = currbase
      pb.actual_height = pb.endy-pb.starty
    end   
    
  end
  
  def get_box_array(y_from_top, chrom_size, total_chrom_y)
    phenoboxes = Array.new
    bin_y = y_from_top
    @phenobins.each do |pbin|
      pbin.boxpos.each_pair do |pos, phenobox|
        phenobox.top_y = bin_y + phenobox.top_y
        phenobox.bottom_y = bin_y + phenobox.bottom_y
        phenobox.chrom_y = pos.to_f / chrom_size * total_chrom_y
        phenoboxes << phenobox
      end
      bin_y += pbin.actual_height
    end
    return phenoboxes
  end
  
  # adds a phenotype to appropriate bin
  def add_phenotype_snp(position, color)
    @phenobins.each do |pb|
      if position.to_i >= pb.startbase and position.to_i < pb.endbase
        pb.add_phenotype_snp(position, color)
        return
      end
    end
  end
  
  def set_pheno_positions
    @phenobins.each do |pb|
      pb.set_phenobox_y
    end
  end
  
  
  def show_bin_info
    @phenobins.each do |pb|
      puts "bin endy=#{pb.endy} startbase=#{pb.startbase} endbase=#{pb.endbase} actual_height=#{pb.actual_height}"
      pb.boxpos.each do |pos, box|
        puts "box pos=#{pos} topy=#{box.top_y} bottomy=#{box.bottom_y}"
      end
      puts "---- END BIN ----"
    end 
  end
  
  # for each phenotype adjust its placement 
  # relative to its enclosing bin
  def adjust_phenos
    i=0
    @phenobins.each do |pb|
      y_span = pb.endy-pb.starty
      bp_span = pb.endbase - pb.startbase
      pb.boxes.each do |pos|
        box=pb.boxpos[pos]
        y_center = (pos.to_i-pb.startbase)/bp_span.to_f * y_span
        box.set_boundaries(y_center)
      end
      i+=1
    end
  end
 
  # change the bp represented in a bin if it the top
  # and bottom boxes have additional space and
  # there isn't enough room to fit all the phenotype boxes
  # contained in it
  def adjust_bp
    @phenobins.each_with_index do |pbin, i|
      pbin.calc_height_needed
      
      pbin.sort_boxes!
      
      # do nothing when the bin has enough space for its phenotypes
      # no need to change the bp spread    
      if pbin.actual_height <= pbin.height_needed or pbin.boxes.empty?
        return
      end
      
      topbp = pbin.startbase
      bottombp = pbin.endbase  
      alter_top = alter_bottom = false
      
      # first check top y on first box to see if it is within the bin
      if pbin.boxpos[pbin.boxes.first].top_y > 0
        topbp = pbin.bp_from_y(pbin.boxpos[pbin.boxes.first].top_y)
        alter_top=true
      end
      if pbin.boxpos[pbin.boxes.last].bottom_y < pbin.actual_height
        #bottombp = pbin.boxes.last
        bottombp = pbin.bp_from_y(pbin.boxpos[pbin.boxes.last].bottom_y)
        alter_bottom=true
      end
      
      # adjust for change in ratio
      #fract_adjust = (bottombp-topbp)/(pbin.endbase-pbin.startbase).to_f
      
      if alter_top
        pbin.startbase = topbp
      end
      if alter_bottom
        pbin.endbase = bottombp
      end
      
      # reposition phenotype boxes if any change in top or bottom
      if alter_top or alter_bottom
        pbin.set_phenobox_y
      end
    end
    
  end
  
  
  
  # place box in appropriate bin based on chromosome y value
  def add_pheno_box(box)
    # integer division works here
    index = box.chromy / @totalchromy
    @phenobins[index] << box
  end
  
  def estimate_height_needed
     total_height = 0
    @phenobins.each do |b|
      total_height+=b.estimate_height
    end
    return total_height
  end
  
  def calculate_total_height
    total_height = 0
    @phenobins.each do |b|
      total_height+=b.actual_height
    end
    return total_height
  end
  
  def calculate_height_needed
     total_height_needed = 0
    @phenobins.each do |b|
      total_height_needed+=b.calc_height_needed
    end
    return total_height_needed   
  end
  
  def total_bins   
    total_height = calculate_total_height
    small_bins = Array.new
    large_bins = Array.new
    
    total_deficit=total_excess=0
    
    # determine bins with too little height
    @phenobins.each_with_index do |b,i|
      if b.height_discrepancy < 0
        small_bins << i
        total_deficit -= b.height_discrepancy
      else
        large_bins << i
        total_excess += b.height_discrepancy
      end
    end
    
    # when all bins are large enough return or if no bins have extra 
    # space 
    if small_bins.empty? or large_bins.empty?
      return 
    end

    if total_deficit > total_excess
      fraction = 1.0
    else
      fraction = total_deficit/total_excess.to_f
    end

    released_height = total_excess * fraction

    # shrink any bins that are too large 
    large_bins.each do |i|
      # take all excess size when needed
      # only take fraction that is needed to make up deficit
      remove_amount = (@phenobins[i].height_discrepancy.to_f/total_excess) * released_height
        @phenobins[i].actual_height -= (@phenobins[i].height_discrepancy.to_f/total_excess) * released_height
    end
    
    # add height to the bins that are too small (add it proportionately 
    # based on the size needed)
    small_bins.each do |i|
      add_amount = (-@phenobins[i].height_discrepancy.to_f/total_deficit) * released_height
      @phenobins[i].actual_height += (-@phenobins[i].height_discrepancy.to_f/total_deficit) * released_height
    end
   
    curr_y=0
    phenobins.each_with_index do |pb,i|
      pb.starty = curr_y
      pb.endy = curr_y + pb.actual_height
      curr_y = pb.endy
    end
    
  end 
 
  def add_chrom(chrom)
    chrom.snps.each_value do |snp|
      snp.phenos.each do |pheno|
        add_phenotype_snp(snp.pos, pheno.color)
      end
    end
  end

end

class ChromosomePlotter < Plotter
  @@chrom_width = 0
  @@num_phenos_row = 6
  
  def self.set_chrom_width(w)
    @@chrom_width = w
  end
  
  def self.set_phenos_row(p)
    @@num_phenos_row = p
    PhenoBox.set_circles_per_row(p)
  end
  
  def self.init_phenobinholder(params)
    binholder = PhenoBinHolder.new
    binholder.totalchromy = params[:chrom_y]
    binholder.totaly = params[:available_y]
    binholder.totalbases = params[:chrom].size
    binholder.set_num_bins(5)
    binholder.set_bases
    return binholder
  end
  
    # utilizes new algorithm to place boxes
  # maintains relative location along chromosome
  def self.position_phenoboxes(params)
    
    totaly = params[:available_y]
    totalchromy = params[:chrom_y]
    chrom = params[:chrom]

    binholder=init_phenobinholder(params)
    # add phenotypes to the bins
    binholder.add_chrom(chrom)

    # check to see if can shrink the plotting area to the chromosome only
    phenobox_offset = -((totaly-totalchromy.to_f)/2)
    estimated_tot=binholder.estimate_height_needed
    
    if estimated_tot < totalchromy
      params[:available_y] = params[:chrom_y]
      binholder=init_phenobinholder(params)
      totaly = totalchromy
      binholder.add_chrom(chrom)
      phenobox_offset = 0
    elsif estimated_tot < totaly
      adjustedy = estimated_tot
      
      params[:available_y] = adjustedy
      binholder=init_phenobinholder(params)
      totaly = adjustedy
      binholder.add_chrom(chrom)
      phenobox_offset = (totalchromy-estimated_tot.to_f)/2        
    end
    
    binholder.set_pheno_positions
    
    # adjust size of bins
    binholder.total_bins 

    # remap the phenotypes to the bins
    binholder.adjust_phenos

    # change bp on bins so top of bin in bp matches first pheno/
    # and bottom of bin in bp matches last pheno
    binholder.adjust_bp

    # place phenoboxes in an array with locations relative to 
    # the absolute base position
    phenoboxes = binholder.get_box_array(phenobox_offset,chrom.size,totalchromy)
    phenoboxes.sort!{|x,y| x.top_y <=> y.top_y}
    adjust_collisions(phenoboxes, totaly)
    
    return phenoboxes
  end


  # move the boxes for collisions
  def self.adjust_collisions(phenoboxes, available_y)
    
    total_boxes = phenoboxes.length
    last_index = total_boxes-1
    notdone = true
    rounds=0

    while(notdone and rounds < 7)
      notdone = false
      for i in (0..last_index)
        for j in (i+1..last_index)
          if phenoboxes[i].bottom_y > phenoboxes[j].top_y
            move_boxes(phenoboxes, available_y, i, j)
            notdone=true
          end
        end
      end
      rounds += 1
    end 
  end
  
  
  def self.shift_boxes(phenoboxes, i, moveup)
    if moveup
      #if i == 0 or phenoboxes[i-1].bottom_y < (phenoboxes[i].top_y-phenoboxes[i].height)
      if phenoboxes[i+1].top_y > (phenoboxes[i].bottom_y-phenoboxes[i].height)
        phenoboxes[i].bottom_y = phenoboxes[i+1].top_y
        phenoboxes[i].top_y = phenoboxes[i].bottom_y - phenoboxes[i].height
      else
        phenoboxes[i].top_y -= phenoboxes[i].height
        phenoboxes[i].bottom_y -= phenoboxes[i].height
      end
    else
       if phenoboxes[i-1].bottom_y < (phenoboxes[i].top_y+phenoboxes[i].height)
         phenoboxes[i].top_y = phenoboxes[i-1].bottom_y
         phenoboxes[i].bottom_y = phenoboxes[i].top_y+phenoboxes[i].height
       else
        phenoboxes[i].bottom_y += phenoboxes[i].height
        phenoboxes[i].top_y += phenoboxes[i].height
       end
    end
    
  end
  
  
  def self.move_boxes(phenoboxes, available_y, i, j)
    last_index = phenoboxes.length-1
    # check for room on top and bottom 
    toproom=bottomroom=false
    if i==0 or (phenoboxes[i].top_y - phenoboxes[i].height > phenoboxes[i-1].bottom_y)
      toproom=true
    end
    if j==last_index or (phenoboxes[j].bottom_y+phenoboxes[j].height < phenoboxes[j].top_y)
      bottomroom=true
    end
    # move in direction with more space available
    # use xor to move in direction when only one available
    if (toproom and bottomroom) or (!toproom and !bottomroom)
      topdist = phenoboxes[i].top_y.to_f / i
      bottomdist = (available_y.to_f - phenoboxes[j].bottom_y)/(last_index-j)
      if topdist > bottomdist
        shift_boxes(phenoboxes,i,true)
      else
        shift_boxes(phenoboxes,j,false)
      end
      
    # when both are ok or neither is, move in direction of more available space
    elsif toproom
      # shift top box up
      shift_boxes(phenoboxes, i, true)
    else
      # shift lower box down
      shift_boxes(phenoboxes, j, false)
    end   
  end
  
  
  def self.plot_chrom(params)
   
    padding = @@circle_size*2
    # leave some space at top and bottom
    available_y = params[:height]-padding
    
    chrom = params[:chrom]
       
    # determine location of the chromosome -- leave a circle at top and bottom
    total_chrom_y = chrom.size.to_f/ @@maxchrom * available_y
    
    start_chrom_y = padding.to_f/2 + (available_y-total_chrom_y).to_f/2
    end_chrom_y = start_chrom_y + total_chrom_y
    
    canvas = params[:canvas]
    xbase = params[:xstart]
    ybase = params[:ystart]
    
    centromere = chrom.centromere
    centromere_y = total_chrom_y * (centromere/chrom.size.to_f) + start_chrom_y
    
    draw_chr(:canvas=>canvas, :centromere_y=>centromere_y, :start_chrom_y=>start_chrom_y, 
      :end_chrom_y=>end_chrom_y, :xbase=>xbase, :ybase=>ybase, :chromnum=>chrom.display_num)

    circle_start_x = @@circle_size*3

    if params[:alt_spacing]==:alternative
      phenoboxes = position_phenoboxes(:chrom=>chrom, :available_y=>available_y, :chrom_y=>total_chrom_y)

      phenoboxes.sort!{|x,y| x.top_y <=> y.top_y}
      phenoboxes.each do |box|
        draw_phenos(canvas, box, circle_start_x, xbase, ybase + start_chrom_y)
      end

    elsif params[:alt_spacing]==:equal
      orig_circle=@@circle_size
      phenoboxes = get_pheno_boxes_equal_spacing(available_y, total_chrom_y, chrom)
      # sort by y location
      phenoboxes.sort!{|x,y| x.top_y <=> y.top_y}
      phenoboxes.each do |box|
        draw_phenos_equal(canvas, box, circle_start_x, xbase, ybase + padding.to_f/2)
      end
      @@circle_size=orig_circle
    else
      phenoboxes = get_pheno_boxes(total_chrom_y, chrom)
      phenoboxes.sort!{|x,y| x.top_y <=> y.top_y}
      phenoboxes.each do |box|
        draw_phenos(canvas, box, circle_start_x, xbase, ybase + start_chrom_y)
      end
    end
  end

 
  # ybase is from start of chromosome drawing
  # need start of chromosome and start of bins
  def self.draw_phenos(canvas, phenobox, start_x, xbase, ybase)
    
    y = phenobox.top_y - @@circle_size * 0.75
    x = start_x
    
    canvas.g.translate(xbase,ybase) do |draw|
      draw.line(0, phenobox.chrom_y, @@chrom_width, phenobox.chrom_y).styles(:stroke=>'black', :stroke_width=>1)
      draw.line(@@chrom_width,phenobox.chrom_y,start_x,phenobox.top_y).styles(:stroke=>'black',:stroke_width=>1)
    end
    
    phenobox.phenocolors.each_with_index do |color, i|
      if i % @@num_phenos_row == 0
        y += @@circle_size * 0.75
        x = start_x
      end
      canvas.g.translate(xbase, ybase) do |draw|
        draw.circle(@@circle_size.to_f/2, x, y).styles(:fill=>color, :stroke=>'black')
      end
      
      x += @@circle_size
    end
    
  end
  
  
  # ybase is from start of chromosome drawing
  def self.draw_phenos_equal(canvas, phenobox, start_x, xbase, ybase)
    
    y = phenobox.top_y - @@circle_size * 0.75
    x = start_x
    
    canvas.g.translate(xbase,ybase) do |draw|
      draw.line(0, phenobox.chrom_y, @@chrom_width, phenobox.chrom_y).styles(:stroke=>'black', :stroke_width=>1)
      draw.line(@@chrom_width,phenobox.chrom_y,start_x,phenobox.top_y).styles(:stroke=>'black',:stroke_width=>1)
    end
    
    phenobox.phenocolors.each_with_index do |color, i|
      if i % @@num_phenos_row == 0
        y += @@circle_size * 0.75
        x = start_x
      end
      canvas.g.translate(xbase, ybase) do |draw|
        draw.circle(@@circle_size.to_f/2, x, y).styles(:fill=>color, :stroke=>'black')
      end
      
      x += @@circle_size
    end
    
  end
  
  # adds phenotypes to sets matching same SNP
  def self.get_pheno_boxes(total_chrom_y, chrom)
    pheno_boxes = Array.new
    
    chrom.snps.each_value do |snp|
      # center of first circle that will be pointed at 
      pos_fraction = snp.pos.to_f/chrom.size
      y_offset = pos_fraction * total_chrom_y
      phenobox = PhenoBox.new
      if pos_fraction <= 0.5
        phenobox.up = true
      else
        phenobox.up = false
      end
      snp.phenos.each do |pheno|
        phenobox.add_phenocolor(pheno.color)
      end
      phenobox.set_default_boundaries(y_offset)
      pheno_boxes << phenobox
    end
    
    return pheno_boxes
  end


  # for this case can spread the snps out along the chromosome and the entire
  # box dedicated to the chromosomes
  def self.get_pheno_boxes_equal_spacing(available_y, total_chrom_y, chrom)
    
    pheno_boxes = Array.new
    # determine amount of y dedicated to each snp
    y_per_snp = available_y/chrom.snpnames.length.to_f
    y = 0
    chrom_offset = (available_y - total_chrom_y.to_f)/2   

    # change circle size based on number of snps
    if y_per_snp < 20
      if y_per_snp > 10
        @@circle_size = y_per_snp
      else
        @@circle_size = 10
      end
    end

    # sort by position
    chrom.sort_snps!
    chrom.snpnames.each do |snpname|
      snp = chrom.snps[snpname]
      pos_fraction = snp.pos.to_f/chrom.size
      y_offset = pos_fraction * total_chrom_y + chrom_offset
      phenobox = PhenoBox.new
      snp.phenos.each do |pheno|
        phenobox.add_phenocolor(pheno.color)
      end
      phenobox.set_even_boundaries(y_offset, y)
      pheno_boxes << phenobox
      y += y_per_snp
    end
    
    return pheno_boxes
  end
  
  
  def self.draw_chr(params)
    canvas = params[:canvas]
    centromere_y = params[:centromere_y]
    start_chrom_y = params[:start_chrom_y]
    end_chrom_y = params[:end_chrom_y]
    xbase = params[:xbase]
    ybase = params[:ybase]
    number = params[:chromnum]
    centromere_offset = @@circle_size.to_f/2
    tpath = "M0,#{start_chrom_y} C0,#{start_chrom_y-@@circle_size/2} #{@@chrom_width},#{start_chrom_y-@@circle_size/2} #{@@chrom_width},#{start_chrom_y}"
    bpath = "M0,#{end_chrom_y} C0,#{end_chrom_y+@@circle_size/2} #{@@chrom_width},#{end_chrom_y+@@circle_size/2} #{@@chrom_width},#{end_chrom_y}"

    canvas.g.translate(xbase,ybase) do |draw|
      draw.line(0,start_chrom_y,0,centromere_y-centromere_offset).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.line(0,centromere_y-centromere_offset,centromere_offset.to_f/2,centromere_y).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.line(centromere_offset.to_f/2,centromere_y,0,centromere_y+centromere_offset).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.line(0,centromere_y+centromere_offset,0,end_chrom_y).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.path(tpath).styles(ChromStyle)
      
      draw.line(@@chrom_width,start_chrom_y,@@chrom_width,centromere_y-centromere_offset).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.line(@@chrom_width,centromere_y-centromere_offset,@@chrom_width-centromere_offset.to_f/2,centromere_y).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.line(@@chrom_width-centromere_offset.to_f/2,centromere_y,@@chrom_width,centromere_y+centromere_offset).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.line(@@chrom_width,centromere_y+centromere_offset,@@chrom_width,end_chrom_y).styles(:stroke=>'darkgray',:stroke_width=>2)
      draw.path(bpath).styles(ChromStyle)
      
    end
    
    font_size = @@circle_size

    canvas.g.translate(xbase,ybase).text(@@chrom_width.to_f/2,end_chrom_y+2*@@circle_size) do |write|
      write.tspan(number.to_s).styles(:font_size=>font_size, :text_anchor=>'middle')
    end
    
  end
end

class Title < Plotter
  
  def self.draw(params)
    xbase = params[:ybase]
    ybase = params[:xbase]
    anchor = params[:anchor] || 'middle'
    
  end
  
  def self.draw_center(params)
    xtotal = params[:xtotal]
    xcenter = xtotal.to_f/2
    title = params[:title]
    ypos = params[:ypos]
    canvas = params[:canvas]
    
    font_size = @@circle_size * 3
    
    canvas.g.translate(xcenter,ypos).text(0,0) do |write|
      write.tspan(title).styles(:font_size=>font_size, :text_anchor=>'middle')
    end
    
  end
  
end


class PhenotypeLabels < Plotter
  
  def self.draw(params)
    canvas=params[:canvas]
    ystart=params[:ystart]
    xstart=params[:xstart]
    phenoholder=params[:phenoholder]
    phenos_per_row=params[:pheno_row]
    xtotal=params[:xtotal]
    
    pheno_space = xtotal.to_f/phenos_per_row
    
    radius = @@circle_size.to_f/2
    y = 0
    x = 0
    i=0
    
    phenokeys = phenoholder.phenonames.keys.sort{|a,b| a.capitalize <=> b.capitalize}
    
    phenokeys.each do |key|
      pheno = phenoholder.phenonames[key]
      if i % phenos_per_row == 0
        x = 0
        y += @@circle_size * 2
      end
      canvas.g.translate(xstart,ystart) do |draw|
        draw.circle(radius, x, y).styles(:fill=>pheno.color, :stroke=>'black')
      end
      canvas.g.translate(xstart,ystart).text(x+@@circle_size*1.5,y+@@circle_size.to_f/4) do |text|
        text.tspan(pheno.name).styles(:font_size=>22)
      end
      x += pheno_space
      i += 1
    end
  end
end

options = Arg.parse(ARGV)

srand(options.rand_seed)

genome = Genome.new
phenoholder = PhenotypeHolder.new(:color=>options.color)
filereader = PhenoGramFileReader.new

filereader.parse_file(options.input, genome, phenoholder)

if options.pheno_spacing == 'standard'
  alternative_pheno_spacing = :standard
elsif options.pheno_spacing == 'equal'
  alternative_pheno_spacing = :equal
else
  alternative_pheno_spacing = :alternative
  
end

# 6 phenotypes in a row (when more wrap around with some overlap to show they are
# linked) -- main problem will be the verical spacing of the dots
# will be done in matched pairs (1 & 13, 2 & 14 etc.)
# 
# 1.0 inches at top
# 4.5 inches for top row of chromosomes
# 2.5 inches for bottom row of chromosomes
# 0.5 inches between two rows
# 1.0 inches to list below
# 8.0 inches to display phenotypes (probably make this dynamic)

# so 19.5 inches vertical X 10 inches horizontal

# a circle will be 5Y in height
# max chrom size (1) will be 40 circles (or 200Y) in height
# chromosome will be 2 circles wide
circle_size = 20
num_circles_in_row=7
Plotter.set_circle(circle_size)
Plotter.set_maxchrom(Chromosome.chromsize(1))
chrom_width = circle_size * 1.5
chrom_circles_width = circle_size * num_circles_in_row
chrom_box_width = chrom_circles_width + chrom_width
ChromosomePlotter.set_chrom_width(chrom_width)
ChromosomePlotter.set_phenos_row(num_circles_in_row-1)
title_margin = circle_size * 7
first_row_start = title_margin

max_chrom_height = 40 * circle_size
max_chrom_box = max_chrom_height + circle_size * 4
# X chromosome will be largest chromosome in second row
second_row_start = max_chrom_box + first_row_start

second_row_start -= circle_size*5 if alternative_pheno_spacing == :standard

second_row_box_max = max_chrom_box * Chromosome.chromsize(23)/Chromosome.chromsize(1)

phenotypes_per_row = 5
phenotype_rows = phenoholder.phenonames.length/phenotypes_per_row
phenotype_rows += 1 unless phenoholder.phenonames.length % phenotypes_per_row == 0

total_y = second_row_start + second_row_box_max

# each row should be 2 circles high + 2 circle buffer on top
phenotype_labels_total = circle_size * 2 * (phenotype_rows+1)
phenotype_labels_y_start = total_y
total_y += phenotype_labels_total
# total y for now
width_in = 8

padded_width = circle_size
# total_y is 2 for chromsome width 6 for circles * number of chroms + space on sides
total_x = (chrom_width + num_circles_in_row * circle_size )* 12 + padded_width * 2
# height can now be determined based on a width of 10 and the ratios
# of the total x and total y
height_in = width_in * total_y / total_x.to_f

xmax = total_x
ymax = total_y

inches_ratio = height_in/width_in.to_f
coord_ratio = ymax/xmax.to_f

rvg=RVG.new(width_in.in, height_in.in).viewbox(0,0,xmax,ymax) do |canvas|
  canvas.background_fill = 'rgb(255,255,255)'
  xstart = padded_width
 
  Title.draw_center(:canvas=>canvas, :title=>options.title, :ypos=>title_margin-circle_size*3,
    :xtotal=>xmax)
  
  # draw each chromosome (first 12 on top row and then second 12 below)
  (1..12).each do |chr|
    ChromosomePlotter.plot_chrom(:canvas=>canvas, :chrom=>genome.chromosomes[chr], 
      :xstart=>xstart, :ystart=>first_row_start, :height=>max_chrom_height, :alt_spacing=>alternative_pheno_spacing)
    xstart += chrom_box_width 
  end
  
  xstart = padded_width
  (13..24).each do |chr|
    ChromosomePlotter.plot_chrom(:canvas=>canvas, :chrom=>genome.chromosomes[chr], 
      :xstart=>xstart, :ystart=>second_row_start, :height=>second_row_box_max, :alt_spacing=>alternative_pheno_spacing)    
    xstart += chrom_box_width
  end
  
  PhenotypeLabels.draw(:canvas=>canvas, :xstart=>padded_width, :ystart=>phenotype_labels_y_start,
    :phenoholder=>phenoholder, :pheno_row=>phenotypes_per_row, :xtotal=>xmax-padded_width)
  
end


# produce output file
outfile = options.out_name + '.' + options.imageformat
print "\n\tDrawing #{outfile}..."
STDOUT.flush
img = rvg.draw

img.write(outfile)

print " Created #{outfile}\n\n"