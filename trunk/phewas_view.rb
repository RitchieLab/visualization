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

# Requires rubygems
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
  puts "\nPlease install RMagick -- See documentation for PheWAS-View or http://rmagick.rubyforge.org/install-faq.html "
  puts
  exit(1)
end

require 'optparse'
require 'ostruct'
include Magick

RVG::dpi=600

Version = '0.3.2'

# check for windows and select alternate font based on OS
Font_family_style = RUBY_PLATFORM =~ /mswin/ ? "Verdana" : "Times"
Font_plot_family = RUBY_PLATFORM =~ /darwin/ ? "Helvetica" : "Helvetica"

#Font_plot_family = 'Helvetica'

#Font_phenotype_names = RUBY_PLATFORM =~ /darwin/ ? "Geneva" : "Helvetica"

Name = 'phewas_view.rb'

if RUBY_PLATFORM =~ /darwin/
  Font_phenotype_names = "Geneva"
elsif RUBY_PLATFORM =~ /mswin/
  Font_phenotype_names = "Arial"
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
    options.phewasfile = nil
    options.out_name = 'phewas_view'
    options.highres = nil
    options.imageformat = 'png'
    options.groupfile = nil
    options.genefile = nil
    options.phenotype_listfile = nil
    options.phenotype_correlations_file = nil
    options.title = " "
    options.largetext = false
		options.boldfont = false
    options.p_thresh = 0.0
    options.phenofile = nil
		options.group_colors = false
    options.ethlist = Array["AA","EA","MA", "H", "AP", "NA"]
    options.snpid = nil
    options.maxp_to_plot = 1.0
    options.beta = false
    options.rotate = false
    options.ethmapfile = nil
    options.lowres = false
    options.plot_sampsizes =false
    options.nolines = false
    options.p_cut = 1.0
    options.phenoname = nil
    options.sun_file = false
    options.genename = nil
    options.labelgene = false
    options.include_eth = false
    options.redline = nil
    options.classname = nil
    options.showbest = false
		options.display_pval = false
    options.include_all_eths = true
		options.group_symbols = false
    options.label_cols = Array['phenotype', 'phenotype_long', 'substudy']
		options.include_genename = false
		options.show_direction = false
		options.narrow = false
    options.tilted_text = false
    help_selected = false
    version_selected = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: phewas_view.rb [options]"
      
      # original phewas_view options
      opts.on("-e [phewas_file]", "PheWAS format file for input") do |phewas_file|
        options.phewasfile = phewas_file
      end
      opts.on("-o [output_name]", "Optional output name for the image") do |output_file|
        options.out_name = output_file
      end
      opts.on("-t [title_str]", "Main title for plot (enclose in quotes)") do |title_str|
        if title_str and title_str =~ /\w/
          options.title = title_str
        else
          options.title = "PheWAS"
        end
      end
      opts.on("-f [image_type]", "Image format for output (png default).  Other options depend on ImageMagick installation.") do |image_type|
        options.imageformat = image_type
      end
      opts.on("-w","--lowres", "Low resolution image (72 dpi)") do |lowres|
        options.lowres = true
      end
      opts.on("-a","--rotate", "Rotates final image 90 degrees") do |rot|
        options.rotate=true
      end
			opts.on("-n","--incl-gene-name","Include italicized gene name in label"){|gname|options.include_genename=true}
			opts.on("-W","--show-direction","Plot p-values as triangles displaying direction of effect"){|triangle|options.show_direction=true}
			opts.on("-j","--narrow","Decrease distance between points on track"){|narrow|options.narrow=true}
      opts.on("-Q","--tilted","Tilt labels"){|narrow|options.tilted_text=true}
      opts.on("-p [pthresh]", "p value threshold, values less significant will be plotted in gray") do |pthresh|
        options.p_thresh = pthresh.to_f
      end
      opts.on("-m, --maxp [max_p_value]", "Maximum p value to plot.  Values greater are not plotted") do |max_p_value|
        options.maxp_to_plot = max_p_value.to_f
      end
      opts.on("-R, --redline [redline]", "Draw red line at designated p value") do |redline|
        options.redline = redline
      end
      opts.on("-b", "Include beta on plot") do |beta|
        options.beta = true
      end    
      opts.on("-A", "--samp-size", "Include sample size plot") do |samp|
        options.plot_sampsizes = true
      end
      opts.on("-l [group_map]", "Optional group map file") do |ancestry_map|
        options.ethmapfile = ancestry_map
      end
			opts.on("-y", "--group-symbols", "Groups displayed with different shapes") do |sym|
				options.group_symbols=true
			end
			opts.on("-D", "--group-color", "Use colors from group file for point colors"){|gcolor| options.group_colors=true}
      opts.on("-c [class_name]", "Only results matching this phenotype class name are plotted") do |class_name|
        options.classname = class_name
      end   
       opts.on("-B", "--showbest", "Display detailed information for best score at each phenotype") do |best|
        options.showbest = true
      end     
      opts.on("-x [pheno_file]", "PheWAS expected phenotypes file") do |pheno_file|
        options.phenofile = pheno_file
      end
      opts.on("-r [group]", Array, "List of groups to include (AA, EA, MA)") do |ancestry|
        options.ethlist = ancestry
        options.include_all_eths = false
      end
      opts.on("-s [snp_id]", "SNP to display from input file") do |snp_id|
        options.snpid = snp_id
      end
      opts.on("-L [phenotype_list_file]", "Optional phenotype list for inclusion") do |phenotype_list_file|
        options.phenotype_listfile = phenotype_list_file
      end
      opts.on("-N","--no-lines", "No lines drawn on plot") do |nolines|
        options.nolines = true
      end
      opts.on("-C [phenotype_correlation_file]", "Optional file with phenotype correlations") do |phenotype_correlation_file|
        options.phenotype_correlations_file = phenotype_correlation_file
      end
      opts.on("-S", "Produce sun plot") do |splot|
        options.sun_file = true
      end
      opts.on("-P [pheno_name]", "Phenotype to display in center of sun plot from input file") do |pheno_name|
        options.phenoname = pheno_name
      end
      opts.on("-g [gene_name]", "Gene to display in center of sun plot from input file") do |gene_name|
        options.genename = gene_name
      end
			opts.on("-d", "--display-pval","Display p-value on sun plot for each 'ray'"){|display_val| options.display_pval=true}
      opts.on("-F", "--large-text", "Use larger text") do |large_text|
        options.largetext = true
      end
			opts.on("-H", "--bod-font") {|boldfont| options.boldfont = true}
      opts.on("-G", "Include gene name along with SNP when SNP is selected for sun plot") do |gene|
        options.labelgene = true
      end
      opts.on("-E", "Include ancestry as description of result for sun plot") do |eth|
        options.include_eth = true
      end
      opts.on("-T [text_columns]",Array, "List columns for inclusion in text label") do |text_columns|
        options.label_cols = text_columns
      end
      opts.on_tail("-h", "--help", "Show this usage statement") do
        #puts opts
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
      #puts
      option_list
      exit(0)
    end

    if options.sun_file
      if !options.phewasfile or (!options.phenoname and !options.snpid and !options.genename)
        option_list
        puts "\nExamples: #{Name}.rb -e NHANES3_PheWAS_Whites_Merge_sorted.txt -o NHANES3_PheWas\n\n"
        puts "          #{Name}.rb -e NHANES3_PheWAS_combined.txt -p 0.1 -t \"rs2197089 results\" -s rs2197089 -o rs2197089_EA_AA\n\n"
        exit(1)
      end      
    else
      if !options.phewasfile
        option_list
        puts "\nExamples: #{Name} -e NHANES3_PheWAS_Whites_Merge_sorted.txt -o NHANES3_PheWas\n\n"
        puts "          #{Name} -e NHANES3_PheWAS_combined.txt -p 0.1 -b -t \"EA AA rs2197089\" -s rs2197089 -r EA,AA -o rs2197089_EA_AA\n\n"
        exit(1)
      end
    end

    return options
  end
  
  # output option list with some options repeated for both Standard and Sun plots
  # need this function as otherwise options will only appear once in the list
  def self.option_list
  print <<EOF
  
Usage: phewas_view.rb [options]

    -h, --help                       Show this usage statement
    -v, --version                    Show version
    -e [phewas_file]                 PheWAS format file for input
    -o [output_name]                 Optional output name for the image
    -t [title_str]                   Main title for plot (enclose in quotes)
    -f [image_type]                  Image format for output (png default).  Other options 
                                     depend on ImageMagick installation.
    -w, --lowres                     Low resolution image (72 dpi)
    -a, --rotate                     Rotates final image 90 degrees
    -p [pthresh]                     p value threshold, values less significant will be 
                                     plotted in gray
    -m, --maxp [max_p_value]         Maximum p value to plot.  Values greater are not 
                                     plotted
    -R, --redline [redline]          Draw red line at designated p value
    -b                               Include beta on plot
    -A, --samp-size                  Include sample size plot
    -l [group_map]                   Optional group file
    -c [class_name]                  Only results matching this phenotype class name 
                                     are plotted
    -B, --showbest                   Display detailed information for best score at each 
                                     phenotype
    -n, --incl-gene-name             Include italicized gene name in label
    -W, --show-direction             Plot p-values as triangles displaying direction of effect
    -j, --narrow                     Decrease distance between points on track
		-y, --group-symbols              Groups displayed with different shapes
    -D, --group-color                Use colors from group file for point colors
    -x [pheno_file]                  PheWAS expected phenotypes file
    -r [groups]                      List of groups to include (AA, EA, MA)
    -s [snp_id]                      SNP to display from input file
    -L [phenotype_list_file]         Optional phenotype list for inclusion
    -N, --no-lines                   No lines drawn on plot
    -C [phenotype_correlation_file]  Optional file with phenotype correlations
    -S                               Produce sun plot
    -s [snp_id]                      SNP to display in center of sun plot from input file
    -P [pheno_name]                  Phenotype to display in center of sun plot from input file
    -g [gene_name]                   Gene to display in center of sun plot from input file
    -G                               Include gene name along with SNP when SNP is selected for 
                                     sun plot
    -F                               Plot with larger text
    -E                               Include ancestry as description of result for sun plot
    -m, --maxp [max_p_value]         Choose a p value threshold, p values less significant will 
                                     not be plotted
    -p [pthresh]                     For plotted results, any results more significant will be 
                                     plotted in red
    -b                               To apply direction of effect for phenotypes in sun plot, 
                                     - is negative direction, + is positive direction

EOF
  end
  
end


############################################################################
#
# Class Ethmap -- Holds labels for race/ethnic group along with keywords
# identifying them and colors to use for display
#
############################################################################
class EthMap
  attr_accessor :ethincluded, :ethindata, :eths, :include_all

  def initialize
    @eths = Hash.new
    @eths["AA"] = Ethnicity.new('blue', 'AA')
    @eths["EA"] = Ethnicity.new('red', 'EA')
    @eths["MA"] = Ethnicity.new('green', 'MA')
    @eths["H"] = Ethnicity.new('purple', 'H')
    @eths["AP"] = Ethnicity.new('orange', 'AP')
    @eths["NA"] = Ethnicity.new('turquoise', 'NA')

    @names = Hash.new
    @names["WHITE"] = @eths["EA"]
    @names["NON-HISPANIC WHITE"] = @eths["EA"]
    @names["BLACK"] = @eths["AA"]
    @names["NON-HISPANIC BLACK"] = @eths["AA"]
    @names["AA"] = @eths["AA"]
    @names["MEXICAN AMERICAN"] = @eths["MA"]
    @names["HISPANIC"] = @eths["H"]
    @names["ASIAN/PACIFIC ISLANDER"] = @eths["AP"]
    @names["AMERICAN INDIAN"] = @eths["NA"]
    @names["EA"] = @eths["EA"]
    @names["MA"] = @eths["MA"]
    @names["H"] = @eths["H"]
    @names["AP"] = @eths["AP"]
    @names["NA"] = @eths["NA"]

    @ethincluded = Hash.new
    @ethincluded["AA"] = true
    @ethincluded["EA"] = true
    @ethincluded["MA"] = true
    @ethincluded["H"] = true
    @ethincluded["AP"] = true
    @ethincluded["NA"] = true

    @ethindata = Hash.new
    @ethindata["AA"] = false
    @ethindata["EA"] = false
    @ethindata["MA"] = false
    @ethindata["H"] = false
    @ethindata["AP"] = false
    @ethindata["NA"] = false
    
    @include_all = true
		@shapes = ['circle', 'square', 'diamond']
  end

  def clear
    @eths.clear
    @names.clear
  end

  # adds major group to map
  def add_group(lab, col)
    @eths[lab] = Ethnicity.new(col, lab)
    @ethincluded[lab] = true  if @include_all
  end

  # adds keyword for identifying group
  def add_keyword(group, key)
    @names[key] = @eths[group]
  end

	def set_shapes
		i=0
		@eths.each do |key,value|
			@eths[key].shape = @shapes[i]
			i+=1
		end
	end
	
  # returns group based on the keyword passed
  def get_group(key)
    @names.each do |keyword, ethnicity|
      if key =~ /^#{keyword}$/i
        return ethnicity
      end
    end
    print "No match for ethnicity #{key} in input file -- ethnicity must match exactly with ethnicity map file\n"
    exit(1)
  end

  def get_group_color(gname)
    return @eths[gname].colorstr
  end
	
	def get_group_shape(gname)
		return @eths[gname].shape
	end

  # sets list of restricted ethnicities for this plot
  def set_restricted_eths(etharray)
    @ethincluded.clear
    etharray.each do |groupname|
      @ethincluded[groupname] = true
    end
  end

  # marks an ethnicity as present in dataset
  def mark_group_present(gname)
    @ethindata[gname] = true
  end

  def included_eth?(eth)
    return @ethincluded[eth]
  end

end #EthMap

############################################################################
#
# Class Ethnicity -- Holds label and color string for a given ethnicity/race
#
############################################################################
class Ethnicity
  attr_accessor :colorstr, :label, :shape

  def initialize(col, lab)
    @colorstr = col
    @label = lab
		@shape = 'circle'
  end

end #Ethnicity


############################################################################
#
# Class CorrelationMatrix -- Holds phenotype correlation matrix
#
############################################################################
class CorrelationMatrix
#  attr_accessor

  def initialize
    @correlationmap = Hash.new
  end

  # stores score in hash
  # try storing each score
  def add_correlation(first_pheno, second_pheno, score)

    @correlationmap[first_pheno] = Hash.new unless @correlationmap.has_key?(first_pheno)
    @correlationmap[first_pheno][second_pheno] = score

    @correlationmap[second_pheno] = Hash.new unless @correlationmap.has_key?(second_pheno)
    @correlationmap[second_pheno][first_pheno] = score
    
  end

  # returns the score from hash
  def get_correlation(first_pheno, second_pheno)
    a = [first_pheno, second_pheno]
    a.sort!
    if @correlationmap.has_key?(a[0])
      return @correlationmap[a[0]][a[1]]
    else
      return nil
    end
  end

  def get_max_correlations
    max=0
    @correlationmap.keys.each do |key|
      if @correlationmap[key].values.length > max
        max = @correlationmap[key].values.length
      end
    end
    return max
  end

end


############################################################################
#
# Class ResultHolder -- Holds results that can be retrieved either by
#
#
############################################################################
class ResultHolder
  attr_accessor :snp_list, :pheno_list, :results, :minpval, :maxpval, :unexpectedcolor,
    :nonsigcolor, :expectedcolor, :single_snp, :maxbeta, :minbeta, :ethcolors, :ethhash,
    :ethmap, :included_phenos, :include_all_phenos, :correlations, :minsampsize, :maxsampsize,
    :max_best_title, :best_results

  def initialize
    @snp_list = SnpList.new
    @pheno_list = PhenoList.new
    @results = Array.new
    @results_by_snp_pheno = Hash.new
    @minpval = @maxpval = 0
    @nonsigcolor = 'gray'
    @expectedcolor = 'purple'
    @unexpectedcolor = 'blue'
    @single_snp = nil
    @correlations = CorrelationMatrix.new
    @maxbeta = 0
    @minbeta = 0
    @maxsampsize = 0
    @minsampsize = 100000
    @ethmap = nil
    @include_all_phenos = true
    @included_phenos = Hash.new
    @best_results = Hash.new
    @max_best_title = 1
    @snpgenes = Hash.new
    @gene_name=nil
  end

  # return true when phenotype is present
  def contains_pheno?(phenoname)
    return @pheno_list.valid_pheno?(phenoname)
  end

  # return a standard ethnicity code for string passed
  def get_matching_eth(eth)
    ethn = @ethmap.get_group(eth)
    return ethn.label
  end

  def set_included_phenos(phenohash)
    @included_phenos = phenohash
    @include_all_phenos = false
  end

  # returns maximum number of interactions for any phenotype
  def get_max_correlations
    return @correlations.get_max_correlations
  end
  
  def add_snp_gene(snp, gene)
    @snpgenes[snp]=gene
  end

  # sets best values for each phenotype
  def set_best_values
    
    @best_scores=Hash.new
    
    # key will be phenotype name and value will be name 
    @pheno_list.pheno_order.each do |phenoname|
      @best_scores[phenoname]=0
      @snp_list.snp_order.each do |snpname|
        if !@results_by_snp_pheno[snpname].has_key?(phenoname)
          next
        end
        result = @results[@results_by_snp_pheno[snpname][phenoname]]
        result.pvalues.each do |ethnicity, resultScore|
          next unless resultScore
          if resultScore.pval > @best_scores[phenoname]
            betamark = ''
            if resultScore.beta
              if resultScore.beta.to_f > 0
                betamark = '+'
              else
                betamark = '-'
              end
            end
            genename = @snpgenes[snpname] || ''
            @best_results[phenoname] = "#{snpname} #{@snpgenes[snpname]} #{betamark}"
            @best_scores[phenoname]=resultScore.pval
          end
        end
      end
    end
    
    @best_results.each_value do |title|
      @max_best_title = title.length if title.length > @max_best_title
    end
    
  end
  
  
  # creates array of arrays of ResultValues to plot
  # each sub array contains the values for one phenotype
  def generate_result_values(pval_thresh, params)
		params[:group_colors] ? group_colors = true : group_colors = false
    @pval_threshold = pval_thresh

    resultvalues = Array.new

    @pheno_list.pheno_order.each do |phenoname|
      pheno_resval = Array.new
      @snp_list.snp_order.each do |snpname|
        if !@results_by_snp_pheno[snpname].has_key?(phenoname)
          next
        end
        result = @results[@results_by_snp_pheno[snpname][phenoname]]

        # result exists so create ResultValues
        result.pvalues.each do |ethnicity, resultScore|
          if !resultScore
            next
          end
          resval = ResultValue.new
          resval.pval = resultScore.pval
          resval.betaval = resultScore.beta
          resval.sampsize = resultScore.sampsize
          resval.snpname = resultScore.snpname
					if @single_snp or group_colors
            resval.colorstr = get_single_pval_color(resultScore.pval, result, ethnicity)
						resval.shape = get_single_pval_shape(ethnicity)
          else
            resval.colorstr = get_pval_color(resultScore.pval, result)
          end

          if resval.pval < @pval_threshold
            resval.belowthresh = true
          end

          pheno_resval << resval
        end

      end
      resultvalues << pheno_resval
    end

    return resultvalues

  end

  # return colors based on whether score is above or below pvalue threshold
  def get_pval_color(pval, result)
    if pval < @pval_threshold.to_f
      return @nonsigcolor
    elsif @pheno_list.expected_phenotype?(result.pheno, result.SNP)
      return @expectedcolor
    else
      return @unexpectedcolor
    end
  end

  # resturns color for plotting based on ethnicity
  def get_single_pval_color(pval, result, ethnicity)
    return @ethmap.get_group_color(ethnicity)
  end

	def get_single_pval_shape(ethnicity)
		return @ethmap.get_group_shape(ethnicity)
	end

  # adds the pvalue result for a SNP-phenotype combination
  def add_result(snpname, phenoname, pval, ethnicity, phenogroup, betaval, 
			sampsize,genename,notes)
    # don't insert blanks
    if pval !~ /\d/
      return
    end

    if pval.to_f == 0.0
      print "P value of 0.0 for #{phenoname} #{snpname}\n"
      @pheno_list.add_pheno(phenoname, genename,notes)
      return
    end

    pval = get_log10(pval)

    if pval > @maxpval
      @maxpval = pval
    end

    if betaval
      if betaval > @maxbeta
        @maxbeta = betaval
      end

      if betaval < @minbeta
        @minbeta = betaval
      end
    end

    if sampsize
      if sampsize < @minsampsize
        @minsampsize = sampsize
      end
      if sampsize > @maxsampsize
        @maxsampsize = sampsize
      end
    end

    if !@results_by_snp_pheno.has_key?(snpname)
      @results_by_snp_pheno[snpname] = Hash.new
    end

    if !@results_by_snp_pheno[snpname][phenoname]
      @results_by_snp_pheno[snpname][phenoname] = @results.length

      phenogroup = pheno_list.get_group_name(phenogroup)

      @results << Result.new(phenoname, snpname, phenogroup)
    end

    if @results[@results_by_snp_pheno[snpname][phenoname]].pvalues.has_key?(ethnicity)
      if pval.to_f == @results[@results_by_snp_pheno[snpname][phenoname]].pvalues[ethnicity].pval
        print "\nSNP #{snpname} has duplicate entries for\nphenotype = #{phenoname} ethnicity = #{ethnicity} -- Plotted only one\n\n"
      else
        print "\nSNP #{snpname} has duplicate entries for\nphenotype = #{phenoname}, ethnicity = #{ethnicity} with different p values -- Not plotted\n\n"
        @results[@results_by_snp_pheno[snpname][phenoname]].pvalues[ethnicity] = nil
        return
      end
    end

    @results[@results_by_snp_pheno[snpname][phenoname]].pvalues[ethnicity] =  ResultScore.new(pval, betaval, sampsize, snpname)

    @snp_list.add_snp(snpname)
    @pheno_list.add_pheno(phenoname, genename, notes)
  end

  # return result based on snp name and phenotype
  def get_result(snpname, phenoname)
    if @results_by_snp_pheno[snpname].has_key?(phenoname)
      return @results[@results_by_snp_pheno[snpname][phenoname]]
    else
      return nil
    end
  end

  # returns longest label to help adjust size of bottom buffer
  def get_longest_label
    return pheno_list.get_longest_name
  end

  # returns number of phenotypes in holder
  def get_num_phenos
    return @pheno_list.pheno_order.length
  end

  # converts value to -log10
  def get_log10(val)
    if val.to_f <= 0
      print "p value out of range #{val}\n"
    end
    return -Math.log10(val.to_f)
  end

end

############################################################################
#
# Class Result -- Holds pvalue for the phenotype and SNP combination
#
############################################################################
class Result
  attr_accessor :pheno, :SNP, :pvalues, :pheno_run_for

  def initialize(phenoname, snpname, pheno_group)
    @pvalues = Hash.new
    @pheno = phenoname
    @SNP = snpname
    @pheno_run_for = pheno_group
  end

end


############################################################################
#
# Class ResultScore -- Holds pvalue for the phenotype and SNP combination and
# beta value if present
#
############################################################################
class ResultScore
  attr_accessor :pval, :beta, :sampsize, :snpname

  def initialize(pv, b, s, snp)
    @pval = pv
    @beta = b
    @sampsize = s
    @snpname = snp
  end

end


############################################################################
#
# Class SnpList -- Holds SNP list with an array that keeps the order of the
# SNPs and a hash that can be used to extract that order based on the
# SNP name.
#
############################################################################
class SnpList
  attr_accessor :snp_order, :snp_hash

  def initialize
    @snp_order = Array.new
    @snp_hash = Hash.new
  end

  def add_snp(snpname)
    if !@snp_hash.has_key?(snpname)
      @snp_hash[snpname] = @snp_order.length
      @snp_order << snpname
    end
  end

  def get_snp_order(snpname)
    return @snp_hash[snpname]
  end

end


############################################################################
#
# Class PhenoList -- Holds Phenotype list
#
############################################################################
class PhenoList
  attr_accessor :pheno_order, :pheno_hash, :gene_name_order, :note_hash

  def initialize
    @pheno_hash = Hash.new
    @pheno_order = Array.new
    @pheno_group = Hash.new
    @expected = Hash.new
		@gene_name_order = Array.new
		@note_hash = Hash.new
  end

  def add_pheno(phenoname,genename,notes)
    if !@pheno_hash.has_key?(phenoname)
      @pheno_hash[phenoname] = @pheno_order.length
      @pheno_order << phenoname
			@gene_name_order << genename
			@note_hash[phenoname]=Hash.new
    end
		if !notes.nil?
			notes.each do |note|
				@note_hash[phenoname][note]=1
			end
		end
		
  end

  def valid_pheno?(phenoname)
    if @pheno_hash.has_key?(phenoname)
      return true
    else
      return false
    end
  end

  def add_expected_pheno(phenotype,snps)
    @expected[phenotype]=Hash.new unless @expected[phenotype]
    snps.each {|snpname|@expected[phenotype][snpname]=true}
  end

  # adds keyword to hash that includes groupname
  def add_pheno_keyword(group, keyword)
    if !@pheno_group.has_key?(group)
      @pheno_group[group] = Hash.new
    end
    @pheno_group[group][keyword]=1
  end


  # returns phenotype group name based on phenotype passed (which should be one from candidate
  # study) 
  def get_group_name(phenostr)
    val = "none"
    @pheno_group.each do |phenogroup, keywordhash|
      keywordhash.each_key do |keyword|
        if phenostr =~ /#{keyword}/i
          val = phenogroup
          break
        end
      end
    end
    return val
  end

  
  def expected_phenotype?(phenotype, snpname)
    return @expected[phenotype][snpname] if @expected[phenotype]
    return false 
  end
  

  # returns true when result is expected by the phenotypes
  # otherwise returns false
  def phenotype_expected?(group, keywords)
    val = false   
    if @pheno_group.has_key?(group)
      @pheno_group[group].each_key do |phenokey|
        # check to see if the phenotype keywords appear in the phenotype used
        if keywords =~ /#{phenokey}/i
          val = true
          break
        end
      end
    end
    return val
  end


  # return the longest name of the phenotypes (in characters)
  def get_longest_name
    max = 0
    @pheno_order.each do |pheno|
      if pheno.length > max
        max = pheno.length
      end
    end
    return max
  end

end


############################################################################
#
# Class FileReader -- Base class
#
############################################################################
class FileReader

  # strips and splits the line
  # returns the data array that results
  def strip_and_split(line)
    line.rstrip!
    line.split(/\t/)
  end

  # strips and then splits on comma
  def strip_and_split_comma(line)
    line.rstrip!
    line.split(/,/)
  end

end


############################################################################
#
# Class ResultFileReader -- Read result file and store in ResultHolder
#
############################################################################
class ResultFileReader < FileReader

  # main function for reading file
  #def read_file(resultholder, filename, grouprequired, samprequired)
  def read_file(params)
    resultholder = params[:resultholder]
    filename= params[:filename]
    grouprequired =params[:grouprequired]
    samprequired = params[:samprequired]
		include_genename = params[:include_genename] || nil
    classname = params[:classname] || nil

    firstline = true

    restrict_snp = resultholder.single_snp

    @snpcol=@phenocol=@pvalcol=@ethcol=@groupcol=@betacol=
      @phenolongcol=@sampsizecol=@phenoclasscol=@genecol= nil    
    lineno = ethcounter = 0
    resultsfound=0
    File.open(filename, "r") do |file|

      file.each_line do |oline|
        oline.each_line("\r") do |line|
        lineno += 1

        # skip blank lines
        if line !~ /\w/
          next
        end

        # for first line read the column headers
        if firstline
          get_columns(line)
          firstline = false
          next
        end

        # standard line of input - add result to ResultHolder
        data = strip_and_split(line)
        if !@snpcol
          print "\nNo SNP column header in file #{filename}\n\n"
          exit
        elsif !@phenocol
          print "\nNo Phenotype column header in file #{filename}\n\n"
          exit
        elsif !@pvalcol
          print "\nNo P_value column header in file #{filename}\n\n"
          exit
        elsif !@sampsizecol and samprequired
          print "\nNo sample_size column in file #{filename}\n\n"
          exit
        elsif !@phenoclasscol and classname
          print "\nNo phenotype_class column in #{filename}\n\n"
          exit          
        elsif !@ethcol
          resultholder.single_snp = nil
        end

        # skip results when only a single SNP is needed
        if (restrict_snp and restrict_snp != data[@snpcol]) or 
            (classname and data[@phenoclasscol] != classname)
          next
        end
        
        resultsfound += 1  
          
        resultholder.add_snp_gene(data[@snpcol], data[@genecol]) if @genecol

        betaval = nil
        betaval = data[@betacol].to_f if @betacol
        assoc_pheno = ""
        assoc_pheno = data[@groupcol] if @groupcol
        sampsize = nil
        sampsize = data[@sampsizecol].to_i if @sampsizecol
        
        if @ethcol
          ethnicity = resultholder.get_matching_eth(data[@ethcol])
          # check that ethnicity should be included before continuing
          if resultholder.ethmap.included_eth?(ethnicity)
            resultholder.ethmap.mark_group_present(ethnicity)
          else
            next
          end
        else
          ethnicity = ethcounter.to_s
          ethcounter +=1
        end
        
        if @phenolongcol
          phenotypename = data[@phenocol] + ": " + data[@phenolongcol]
        else
          phenotypename = data[@phenocol]
        end
				if include_genename and @genecol
					phenotypename = data[@genecol] + ' ' + phenotypename
				end
        phenotypename.gsub!('"','')
        
				if @notecol 
					notes = data[@notecol].split(/,/)
				end	
					
				(@genecol.nil? or !include_genename) ? genename = nil : genename = data[@genecol]
        if (resultholder.include_all_phenos or resultholder.included_phenos.has_key?(phenotypename))
          begin
            resultholder.add_result(data[@snpcol], phenotypename, data[@pvalcol], ethnicity,
              assoc_pheno, betaval, sampsize, genename,notes)
          rescue StandardError => problem
            print "Problem adding result for line# #{lineno}:\n#{line}\n\n"
            exit
          end
        else
        end
        end #line
      end #oline
    end #file
    
    if resultsfound==0 and classname
      print "\ntNo results found for class name #{classname}\n\n"
      exit
    end
    
  end


  # gets columns from tab delimited strings in line and sets variables for identifying columns
  def get_columns(line)
    data = strip_and_split(line)

    data.each_with_index do |header, i|
      if header =~ /^snpID|^snp/i
        @snpcol = i
      elsif header =~ /^\s*phenotype\s*$/i
        @phenocol = i
      elsif header =~ /^p_value|^pval|^p-value/i
        @pvalcol = i
      elsif header =~ /Race_ethnicity|ancestry|group/i
        @ethcol = i
      elsif header =~ /Associated_Phenotype/i
        @groupcol = i
      elsif header =~ /^es$|^beta$/i #es for effect size
        @betacol = i
      elsif header =~ /phenotype_long/i
        @phenolongcol = i
      elsif header =~ /sample_size|^N$/i
        @sampsizecol = i
      elsif header =~ /phenotype_class/i
        @phenoclasscol = i
      elsif header =~ /^gene/i
        @genecol = i
			elsif header =~/^notes_beta_feature/i
				@notecol = i
      end
    end #data
  end


end #ResultFileReader



############################################################################
#
# Class ExpectedPhenoReader -- Reads expected phenotype information
#
############################################################################
class ExpectedPhenoReader < FileReader

  # file consists of a phenotype followed by the SNPs expected to be associated
  # each line is a phenotype
  def read_file(resultholder, filename)
    
    # no header in file
    File.open(filename, "r") do |file|
      file.each_line do |oline|
        oline.each_line("\r") do |line|
          #skip blank lines
          next if line !~ /\w/
          
          # split the line with first being the phenotype
          data = strip_and_split(line)
          resultholder.pheno_list.add_expected_pheno(data[0], data[1..data.length-1])
        end
      end
    end
    
  end


  # store column headers for use in setting keywords for membership in the group
  def get_columns(line)
    data = strip_and_split(line)

    data.each_with_index do |header, i|
      @column_headers << header
    end #data
  end

end #ResultFileReader


############################################################################
#
# Class EthMapReader -- Reads race/ethnicity information file
#
############################################################################
class EthMapReader < FileReader

  # main function for reading file
  def read_file(ethmap, filename)

    ethmap.clear

    firstline = true
    @column_headers = Array.new

    File.open(filename, "r") do |file|

      file.each_line do |oline|
        oline.each_line("\r") do |line|
        # skip blank lines
        if line !~ /\w/
          next
        end

        # skip column headers
        if firstline
          firstline = false
          next
        end

        # standard line of input - add result to ResultHolder
        data = strip_and_split(line)

        ethmap.add_group(data[0], data[1])
          
        last_index = data.size-1

        data.values_at(2..last_index).each do |keyword|
          ethmap.add_keyword(data[0], keyword)
        end
        end
      end #line
    end #file

  end


end #EthMapReader


############################################################################
#
# Class CorrelationReader -- Reads phenotype correlations and stores in the
# resultholder.
#
############################################################################
class CorrelationReader < FileReader

  def read_file(filename, resultholder)

    pheno_long_present=false
    
    # determine if have a phenotype long column 
    # along with a phenotype one
    linenum=0
    firstdata = Array.new
    seconddata = Array.new
    startline=2
    tfile = File.new(filename, "r")
    while(oline=tfile.gets)
      oline.each_line("\r") do |line|
        next if line !~ /\w/
        linenum+=1
        if linenum==1
          firstdata=strip_and_split(line)
          firstdata.shift
        elsif linenum==2
          seconddata=strip_and_split(line)
          seconddata.shift
        elsif linenum==3
          data = strip_and_split(line)
          if data[1] =~ /[a-z][A-Z]/
            pheno_long_present=true
            firstdata.shift
            seconddata.shift
            startline=3
          end
          break
        end
      end
    end
    tfile.close
    
    columns = Array.new
    unless pheno_long_present
      firstdata.each do |first|
        name = first
        name.gsub!('"', '')
        columns << name
      end
    else
      firstdata.each_with_index do |first, i|
        name = first + ": " + seconddata[i]
        name.gsub!('"', '')
        columns << name
      end     
    end
    
    # first two lines contain the phenotype names for the columns
    # first two columns are for identifying the horizontal phenotype names
    f = File.new(filename, "r")    
 
    # process all remaining lines
    linenum=0
    while(oline=f.gets)
      oline.each_line("\r") do |line|

      next if line !~ /\w/
      linenum+=1
      next if linenum < startline
           
      # get name first
      data = strip_and_split(line)
      if pheno_long_present
        name = data.shift + ": " + data.shift
      else
        name = data.shift
      end
      name.gsub!('"','')

      # iterate through the columns and only store those scores where phenotypes are requested
      data.each_with_index do |score, i|
        if resultholder.contains_pheno?(name) and resultholder.contains_pheno?(columns[i]) and
            score=~ /\d/
          resultholder.correlations.add_correlation(name, columns[i], score)
        end
      end
      end
    end
  end


end #CorrelationReader

############################################################################
#
# Class PhenoListReader -- Reads phenotypes that are to be included in plot.
# Format should be phenotype and phenotype long name separated by a tab.
#
############################################################################
class PhenoListReader < FileReader

  def read_file(filename)
    phenolist = Hash.new

    firstline = true

    File.open(filename, "r") do |file|
      file.each_line do |oline|
        oline.each_line("\r") do |line|
        if line !~ /\w/ or firstline
          firstline = false
          next
        end
        data = strip_and_split(line)
        #concatenate names
        phenolist[data[0] + ": " + data[1]] = true
      end
    end
    end
    return phenolist
  end

end


############################################################################
#
# Class CodeMapReader -- Reads phenotype code map and uses that to
# assign the correct phenotypes to codes
#
############################################################################
class CodeMapReader < FileReader

  # reads phenotype map file
  def read_file(phenohash, filename)

    phenohash.clear

    firstline = true
    @column_headers = Array.new

    File.open(filename, "r") do |file|

      file.each_line do |oline|
        oline.each_line("\r") do |line|
        # skip blank lines
        if line !~ /\w/
          next
        end

        # skip column headers
        if firstline
          firstline = false
          next
        end

        # standard line of input - add result to ResultHolder
        data = strip_and_split_comma(line)

        if data[1]
          phenohash[data[0]] = data[0] + ": " + data[1]
        else
          phenohash[data[0]] = data[0]
        end
        end
      end #line
    end #file
  end

end



############################################################################
#
# Class Point -
#
############################################################################
class Point
  attr_accessor :x, :y, :uptriangle, :shape

  def initialize(xpt, ypt, sh='circle', uptriang=nil)
    @x = xpt
    @y = ypt
    @uptriangle = uptriang
		@shape = sh
  end

end


############################################################################
#
# Class ResultValue - Contains information about plotting such as color
# and whether there is a beta value for plotting triangles
#
############################################################################
class ResultValue
  attr_accessor :pval, :colorstr, :betaval, :belowthresh, :sampsize, :snpname,
		:shape

  def initialize
    @pval = 0
    @sampsize = 0
    @colorstr = 'gray'
    @betaval = nil
    @belowthresh = false
    @snpname = ""
		@shape = 'circle'
  end

end


############################################################################
#
# Class SunResultHolder -- Hold results
#
############################################################################
class SunResultHolder
  attr_accessor :results, :center_name, :max_title_length, :phenorequired,
    :appendgene, :appendeth, :gene_name

  def initialize
    @results = SunResultList.new
    @center_name = ''
    @max_title_length = 0
    @phenorequired = false
    @appendgene = false
    @appendeth = false
    @gene_name = nil
  end

  def add_score(n,val,betaval)
    if(n.length > @max_title_length)
      @max_title_length = n.length
    end
    @results.add_score(n,val,betaval)
  end

  def total_values
    return @results.scores.length
  end


end  #ResultHolder


############################################################################
#
# Class SunResultList -- Hold results - can be either SNPs or phenotypes
#
############################################################################
class SunResultList
  attr_accessor :scores, :ordered, :max_result

  def initialize
    @scores = Hash.new
    @ordered = Array.new
    @max_result = 0
  end

  def add_score(n, value, betaval)
    @scores[n] = SunScore.new(get_neg_log(value), betaval)
    if @scores[n].pval > @max_result
      @max_result = @scores[n].pval
    end
  end

  # sorts in descending order
  def sort_scores_by_value
    @ordered = @scores.keys
    @ordered.sort! {|x,y| @scores[y].pval <=> @scores[x].pval}
  end

  # returns -log10 of parameter
  def get_neg_log(pval)
    if pval == 0
      return 50
    elsif pval == 1
      return 0.0
    else
      return -Math.log10(pval.to_f)
    end
  end
  
  def get_pval(log10val)
    return 10**-log10val
  end

end #ResultList


############################################################################
#
# Class SunResultList -- Hold results - can be either SNPs or phenotypes
#
############################################################################
class SunScore
  attr_accessor :pval, :beta
  
  def initialize(p,b)
    @pval = p
    @beta = b
  end

end


############################################################################
#
# Class ResultFileReader -- Read result file and store in ResultHolder
#
############################################################################
class SunResultFileReader < FileReader

  # main function for reading file
  def read_file(resultholder, filename, p_cut, label_headers)

    firstline = true

    restrict_name = resultholder.center_name.downcase

    @snpcol = nil
    @phenocol = nil
    @pvalcol = nil
    @ethcol = nil
    @assoc_phenocol = nil
    @betacol = nil
    @phenolongcol = nil
    @sampsizecol = nil
    @studycol = nil
    @substudycol = nil
    @genecol = nil
    @phenoclasscol = nil

    @label_columns = Array.new
    
    genename = nil
    lineno = 0

    File.open(filename, "r") do |file|
      file.each_line do |oline|
        oline.each_line("\r") do |line|
        lineno += 1

        # skip blank lines
        if line !~ /\w/
          next
        end

        # for first line read the column headers
        if firstline
          get_columns(line, label_headers)
          firstline = false
          next
        end

        # standard line of input - add result to ResultHolder
        data = strip_and_split(line)

        if !@snpcol
          print "\nNo SNP column header in file #{filename}\n\n"
          exit
        elsif !@phenocol
          print "\nNo Phenotype column header in file #{filename}\n\n"
          exit
        elsif !@pvalcol
          print "\nNo pval column header in file #{filename}\n\n"
          exit
        end

       resultname = ''

        # skip if it isn't the right SNP or phenotype for center
        snpname = data[@snpcol].downcase
        if data[@pvalcol].to_f <= p_cut
          if snpname.downcase == restrict_name# and data[@pvalcol].to_f <= p_cut
              resultname = ''
              @label_columns.each {|col| resultname += ':' + data[col] if data[col]}
              resultname.slice!(0)        
            genename = data[@genecol] if @genecol && data[@genecol]=~/[a-zA-Z0-9]/
          elsif data[@phenocol].downcase == restrict_name
            resultname = snpname
            if @genecol
              resultname += ':' + data[@genecol]
            end
            if @assoc_phenocol
              resultname += ':' + data[@assoc_phenocol]
            end
          elsif @genecol and restrict_name == data[@genecol].downcase
            resultname = data[@snpcol] +": "+data[@phenocol] 
            resultname += " " + data[@phenolongcol] if @phenolongcol
          else
            # no matching information
            next
          end

          betaval = nil
          if @betacol
            betaval = data[@betacol].to_f
          end
          
          if resultholder.appendeth && @ethcol
            resultname += ' (' + data[@ethcol] + ')'
          end
          resultname.gsub!('"','')          
          resultholder.add_score(resultname, data[@pvalcol].to_f, betaval);
        end # p-value <= p_cut
        end
      end #line
    end #file
    # sort based on value
    resultholder.results.sort_scores_by_value
    
#    resultholder.center_name += "\n#{genename}" if resultholder.appendgene
    resultholder.gene_name = genename if resultholder.appendgene
  end


  # gets columns from tab delimited strings in line and sets variables for identifying columns
  def get_columns(line, label_headers)
    data = strip_and_split(line)

    data.each_with_index do |header, i|
      if header =~ /^\s*snpID\s*$|^snp$/i
        @snpcol = i
      elsif header =~ /^\s*phenotype\s*$/i
        @phenocol = i
      elsif header =~ /p_value|^pval/i
        @pvalcol = i
      elsif header =~ /Race_ethnicity|ancestry|group/i
        @ethcol = i
      elsif header =~ /Associated_Phenotype/i
        @assoc_phenocol = i
      elsif header =~ /beta|es/i
        @betacol = i
      elsif header =~ /phenotype_long/i
        @phenolongcol = i
      elsif header =~ /sample_size/i
        @sampsizecol = i
      elsif header =~ /^study$/i
        @studycol = i
      elsif header =~ /^substudy$/i
        @substudycol = i
      elsif header =~ /^gene$/i
        @genecol = i
      elsif header =~ /phenotype_class/i
        @phenoclasscol = i
      end
    end #data
    
    # ID the label headers
    label_headers.each do |label|
      data.each_with_index do |header,i|
        if header =~ /^\s*#{label}\s*$/i
          @label_columns << i
          break
        end
      end
    end
    
  end

end #ResultFileReader

############################################################################
#
# Class RadialPlotter -- Draws plot
#
############################################################################
class RadialPlotter
  attr_accessor :radius, :font_size_multiple, :canvas, :line_size, :center_size_multiple,
		:title_text_adj

  def initialize
    @radius = 2
    @font_size_multiple = 1
		@center_size_multiple = 1
    @canvas = nil
    @line_size = 20
		@title_text_adj = 1
  end

  # returns plot width in inches
  def calculate_plot_width(radii, padding)
    return radii * @line_size.to_f/97 + padding
  end

  # sets standard font size on the plot
  def standard_font_size
    return @font_size_multiple * @line_size/2 + 1
  end

  # adds space for title
  def add_space_for_title(size_in_inches, fraction)
    return size_in_inches + @line_size * fraction
  end

  # returns value in x,y coordinate for the
  # size in inches passed
  def calculate_coordinate(size_in_inches)
    return size_in_inches * 144
  end

  # add title to the plot
  def write_main_title(title_str, x_start, y_start, x_end, y_end, split_title=false)

    font_size = standard_font_size * 1.4

    if split_title
      line_length = (title_str.length/2).to_i
      # check for space at or beyond line_length char
      start_char = line_length
      while title_str[start_char,1] != " " and start_char < title_str.length-1
        start_char += 1
      end
      if start_char == title_str.length-1
        start_char = line_length
        while title_str[start_char,1] != ' ' and start_char > 0
          start_char -= 1
        end
      end

      # start_char will now contain the split position and
      # a "\n" should be inserted there
      title_str.insert(start_char, "\n")
    end

    @canvas.g.translate(x_start, y_start).text((x_end-x_start)/2, (y_end-y_start)/1.5) do |text|
      text.tspan(title_str).styles(:font_size=>font_size, :text_anchor=>'middle')
    end
  end #write_main_title


	def convert_pval(pval)
		if pval < 0.0001
			label = (sprintf "%.2g", pval)
    elsif pval < 0.001
      label = (sprintf "%.4f", pval)
    elsif pval < 0.01
      label = (sprintf "%.3f", pval)
    else
      label = (sprintf "%.2f", pval)
     end		
		 return label
	end
	
  # draw main plot results with center circle and lines radiating out from that
 # def draw_main_plot(center_name, reslist, midx, ystart, yend, p_thresh, total_radii)
  def draw_main_plot(params)
    center_name = params[:center_name]
    reslist = params[:reslist]
    midx = params[:midx]
    ystart = params[:ystart]
    yend = params[:yend]
    p_thresh = params[:p_thresh]
    total_radii = params[:total_radii]
    use_beta = params[:use_beta] || false
		use_pval = params[:use_pval] || false
		params[:bold_font] ? font_weight = 'bold' : font_weight = 'normal'
#		font_weight = params[:bold_font] || false
    gene_name = params[:gene_name]
    largetext = params[:largetext]
    
    midy = (ystart+yend)/2

    # find maximum value
    maxp = reslist.max_result

    # determine amount to move down and then back up based on number of results
    num_results = reslist.scores.length
    # figure out top to bottom of the range and then use that as the amount to change y
    yinterval =  @line_size * total_radii
    # divide by half as there are two sides to the result plotting
    yincrement = yinterval/(num_results.to_f/2)

    half_total_radii = total_radii.to_f/2

    radius_interval = @line_size * half_total_radii * 0.85 - @line_size * @radius
    standard_interval = @line_size * @radius + @line_size * half_total_radii * 0.15

    max_length = radius_interval + standard_interval

    y_result = midy - radius_interval - standard_interval
    y_result_start = y_result

    midindex = num_results/2-1

    rotate_increment = 360 / num_results.to_f
    rotation = 0

    first_below_rotate = 360

    nonsigcolor = '#858585'
    sigcolor = 'red'

    y_correction = 2.0 + (0.02 * (num_results-50))
    num_results.times do |i|
      stroke_color = 'black'
      if p_thresh > 0
        if reslist.scores[reslist.ordered[i]].pval > p_thresh
          stroke_color = sigcolor
        else
          stroke_color = nonsigcolor
          if first_below_rotate == 360
            first_below_rotate = rotation
          end
        end
      end

      canvas.g.translate(midx, midy).rotate(rotation) do |rotated|
        length = (reslist.scores[reslist.ordered[i]].pval/maxp.to_f) * radius_interval + standard_interval
        rotated.line(0,0,0,-(length)).styles(:stroke=>stroke_color, :stroke_width=>1)
				
      if rotation < 180 and rotation > 0
        txt_anchor = 'start'
      else
        txt_anchor = 'end'
      end

      y_shift = -length-line_size/5
   
      if rotation >= 140 and rotation < 180
        y_shift -= line_size * (rotation-140)/20 * (@font_size_multiple*y_correction*@title_text_adj-1.0)
      elsif rotation >= 180 and rotation < 220
        y_shift -= line_size * (220-rotation)/20 * (@font_size_multiple*y_correction*@title_text_adj-1.0)
      elsif rotation == 0
				largetext ? y_shift -= line_size * (120-rotation)/120 * (@font_size_multiple*2.0*@title_text_adj-1.0) : y_shift -= line_size * 0.18 * @font_size_multiple 
      elsif rotation <= 40
				y_shift -= line_size * (40-rotation)/120 * (@font_size_multiple*@title_text_adj-1.0)
        # adjust it based on length of line size at that position
        # make farther from original if line is shorter
      elsif rotation >= 320
				y_shift -= line_size * (rotation-320)/20 * (@font_size_multiple*2*@title_text_adj-1.0)
      end
				label = reslist.ordered[i]
        if use_beta
          if reslist.scores[reslist.ordered[i]].beta < 0
            label += " -"
          else
            label += " +"
          end
        end
        if i==0 
          pval = reslist.get_pval(reslist.scores[reslist.ordered[i]].pval)
          if pval < 0.0001
            label += "\n" + (sprintf "%.2g", reslist.get_pval(reslist.scores[reslist.ordered[i]].pval))
          elsif pval < 0.001
            label += "\n" + (sprintf "%.4f", reslist.get_pval(reslist.scores[reslist.ordered[i]].pval))
          elsif pval < 0.01
            label += "\n" + (sprintf "%.3f", reslist.get_pval(reslist.scores[reslist.ordered[i]].pval))
          else
            label += "\n" + (sprintf "%.2f", reslist.get_pval(reslist.scores[reslist.ordered[i]].pval))
					end
				elsif use_pval
					formatted_pval = convert_pval(reslist.get_pval(reslist.scores[reslist.ordered[i]].pval))
					rotation < 180 ? label = formatted_pval + " - " + label : label = label +  " - " + formatted_pval
        end
				
				midlabel = label.length/2
				newlabel = label.dup
				
					rotated.text(0,y_shift,newlabel).rotate(-rotation).styles(:font_weight=>font_weight,:font_size=>standard_font_size*0.65,
						:text_anchor=>txt_anchor) #:font_family=>Font_plot_family
      end
      rotation += rotate_increment
    end

    fill_color = '#FF7F24'

    circle_size = @radius * @line_size
    # draw circle of appropriate size in center of the plot
    canvas.circle(circle_size, midx, midy).styles(:fill=>fill_color,
      :stroke=>'red', :stroke_width=>1)

    # if the rotation that fell below the threshold is less than 360
    # calculate the start of that location on the arc and draw over with
    # an arc that is the non-significant color
    # formula is
    # x = x_center + d cos(angle)
    # y = y_center + d sin (angle)
    if first_below_rotate < 360
      if first_below_rotate > 180
        arc_flag = 0
      else
        arc_flag = 1
      end

      if first_below_rotate >= 90 and first_below_rotate < 360
        first_below_rotate -= 90
      else
        first_below_rotate += 270
      end
      # convert degrees to radians
      radians = first_below_rotate * Math::PI / 180

      end_arc_x = midx + @radius * @line_size * Math.cos(radians)
      end_arc_y = midy + @radius * @line_size * Math.sin(radians)

      arc_radius = @radius * @line_size
      final_arc_x = end_arc_x - midx
      final_arc_y = end_arc_y - midy+arc_radius

      canvas.path("M#{midx},#{midy-arc_radius} a#{arc_radius},#{arc_radius} 0 #{arc_flag},0 #{final_arc_x},#{final_arc_y}").
        styles(:fill=>'none', :stroke=>nonsigcolor, :stroke_width=>1)
    end

    # add text for the center circle
    y_circle_text = midy + @line_size/5
    @canvas.text(midx, y_circle_text) do |text|
      text.tspan(center_name).styles(:font_family=>Font_plot_family, :font_size=>standard_font_size*0.8*@center_size_multiple, :text_anchor=>'middle')
    end
    if gene_name
      y_circle_text += circle_size/4.to_f
      @canvas.text(midx, y_circle_text) do |text|
        text.tspan(gene_name).styles(:font_family=>Font_plot_family, :font_size=>standard_font_size*0.8*@center_size_multple, 
          :text_anchor=>'middle', :font_style=>'italic')
      end
    end

  end

end


############################################################################
#
# Class DotPlotter -- Draws actual dot plot for PheWAS results
#
############################################################################
class DotPlotter
  attr_accessor :x_offset, :y_pval_zero, :y_offset, :diameter, :pval_threshold, :x_interior_offset,
    :y_plot_height, :offset_mult, :first_color, :font_adjuster, :font_size_multiple,
		:pheno_dist_mult, :note_font_multiple

  def initialize
    @diameter = 8
    @x_offset = 10
    @y_pval_zero = 120
    @y_offset = 10
    # basic circle size in plot
    @pval_threshold = 0
    @offset_mult = 4
    @first_color = 'gray'
    @font_adjuster = 1
    @font_size_multiple =1.0
		@pheno_dist_mult = 2.0
    @note_font_multiple = 2.5
  end

  # returns plot width in inches
  def calculate_plot_width(num_phenos, padding)
    #return num_phenos * @diameter.to_f/97 + padding
    if num_phenos < 10
      divisor = 140
    elsif num_phenos < 25
      divisor = 142
    elsif num_phenos < 50
      divisor = 144
    elsif num_phenos < 150
      divisor = 144
    elsif num_phenos < 250
      divisor = 146
    elsif num_phenos < 800
      divisor = 148
    else
      divisor = 149
    end
		divisor *= 2.0/@pheno_dist_mult
    return num_phenos * @diameter.to_f/divisor + padding
  end


  def add_grid_size(size_in_inches, horizontal_side_in_inches, max_correlations, num_phenos)

    max_correlations = num_phenos-1

    correlation_fraction = Math.sqrt(max_correlations / Math.sqrt(2) *
      (max_correlations / Math.sqrt(2)) / 2) / num_phenos
   return size_in_inches + horizontal_side_in_inches * correlation_fraction -0.4
  end

  # returns value in x,y coordinate for the
  # size in inches passed
  def calculate_coordinate(size_in_inches)
    return size_in_inches * 300
  end

  def add_space_for_title(size_in_inches, fraction)
    return size_in_inches + @diameter * fraction
  end

  def draw_diameter
    return @diameter / 2
  end

  # adds title to top of plot
  def draw_title(canvas, xstart, ystart, width, title_height, title, rotate=false)
    font_size = standard_font_size

    if rotate
      canvas.g.translate(xstart,ystart).text(title_height, width/2).rotate(-90) do |text|
        text.tspan(title).styles(:font_size=>font_size*1.5, :text_anchor=>'middle')
      end
    else
      canvas.g.translate(xstart,ystart).text(width/2, title_height) do |text|
        text.tspan(title).styles(:font_size=>font_size*1.5, :text_anchor=>'middle')
      end
    end
  end

  def standard_font_size
    return @font_size_multiple * @diameter/0.8 + 1
  end
	
	# return space taken by text in coordinate system
	def space_for_text(coordinates_per_pixel,label)	
		font_size = standard_font_size
		font_size *= 2.14
		pixels = get_text_width(label, :pointsize=>font_size, :font_family=>'arial')
		return pixels * coordinates_per_pixel
  end

	# return text width
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
		return metrics.width
	end
	
  # draws legend identifying groups and colors associated with those groups
  def draw_legend(canvas, ethcolors, ethshapes,xstart, ystart, params)
		rotate = params[:rotate]
		coordinates_per_pixel = params[:coordinates_per_pixel]
		max_label = params[:max_label]
    legend_start = xstart + @offset_mult * @x_offset + @x_interior_offset

    length_x = @diameter * 5
    if rotate
      ystart -= @diameter/2.5
    else
      ystart += -@diameter/3
    end

    font_size = standard_font_size * 1
    text_x = @offset_mult * @x_offset + @x_interior_offset
    y_text_adjust=0

    if rotate
      rotation = -90
      anchor = 'start'
      length_x = @diameter * 2
      y_text_adjust=@diameter*2
    else
      rotation = 0
      anchor = 'end'
			text_x = @offset_mult * @x_offset
    end

    ethcolors.each do |eth, colorstr|
			length_x = space_for_text(coordinates_per_pixel,max_label) unless rotate
			if rotate
				canvas.g.translate(xstart, ystart).text(text_x, @y_offset-y_text_adjust).rotate(rotation) do |text|
					text.tspan(eth).styles(:font_size => font_size, :text_anchor=>anchor)
				end
			else
				canvas.g.translate(xstart, ystart).text(text_x+length_x-@diameter/5, @y_offset-y_text_adjust).rotate(rotation) do |text|
					text.tspan(eth).styles(:font_size => font_size, :text_anchor=>anchor)
				end
			end
			
      canvas.g.translate(xstart, ystart) do |box|
        box.styles(:fill=>colorstr, :stroke=>'none', :stroke_width=>1, :fill_opacity=>0.8)
        if rotate
					if ethshapes[eth]=='square' or ethshapes[eth]=='circle'
						box.rect(@diameter, @diameter, text_x-@diameter, @y_offset-@diameter)
					elsif ethshapes[eth]=='diamond'
#						offset = @diameter/2.to_f 
						offset = @diameter*0.6
						centerx=text_x-@diameter/2.to_f
						centery=@y_offset-@diameter/2.to_f
						xpts =[centerx,centerx+offset,centerx,centerx-offset]
						ypts =[centery-offset,centery,centery+offset,centery]
						box.polygon(xpts,ypts)
					else
						box.circle(@diameter/2.to_f,text_x-@diameter/2.to_f,@y_offset-@diameter/2.to_f)
					end
        else
					if ethshapes[eth]=='square' or ethshapes[eth]=='circle'
						box.rect(@diameter, @diameter, text_x+length_x, @y_offset-@diameter)
					end
					length_x += @diameter + space_for_text(coordinates_per_pixel, " ")
        end

      end
      text_x += length_x
    end
  end


  # based on absolute value of score
  # ranges from blue (for zero) to yellow (for one)
  def get_correlation_color(score)
    # blue is 0-0-255
    # yellow is 255-255-0
    if score == nil
      return 'white'
    else

      # stretch this by compressing original bottom half to
      # 0-0.25 and upper half from 0.25 to 1.0

      score = score.to_f
      score_r = score.abs / 1 * 255
      score_g = score.abs / 1 * 255
      score_b = 255 - score.abs/1 *255

      color = sprintf("rgb(%d, %d, %d)", score_r, score_g, score_b)
      return color
    end
  end

  # draws grid for correlation values of phenotypes
  def draw_grid(canvas, resultholder, x_start, y_start, rotate)
    block_size = Math.sqrt(2 * @diameter * @diameter)
    x_grid_start = x_start + @offset_mult * @x_offset + @x_interior_offset 
    x_grid_start += block_size*0.5 if rotate
    y = y_start #+ @diameter/4
    
    if rotate
      rotation = 45
      #y_new_start = y_start + 10*block_size + block_size * resultholder.pheno_list.pheno_order.length * @correlation_fraction
      y_vertical = ((resultholder.pheno_list.pheno_order.length-1) ** 2) - ((resultholder.pheno_list.pheno_order.length.to_f/2*Math.sqrt(2)) ** 2)
      y_new_start = y_start + Math.sqrt(y_vertical)*block_size
      
      increment = -block_size
      legend_y_start = y_new_start - @diameter*4
    else
      rotation = -45
      y_new_start = y
      increment = block_size
      legend_y_start = y_start+@diameter*2
    end
    
    draw_plot_legend(canvas, x_grid_start, legend_y_start, rotate)
    # iterate through the phenotypes and set blocks for each correlation
    canvas.g.translate(x_grid_start, y_new_start).rotate(rotation) do |ldbox|
      # draw scores for each combination present
      x=0
      y_current_start = -increment
      resultholder.pheno_list.pheno_order.each_with_index do |pheno, i|
        y_current_start += increment
        y=y_current_start
        (i+1..resultholder.pheno_list.pheno_order.length-1).each do |j|
          # get the score for this combination
          score = resultholder.correlations.get_correlation(pheno, resultholder.pheno_list.pheno_order[j])
          color = get_correlation_color(score)

          stroke_color = 'lightgray'

          ldbox.rect(block_size, block_size, x, y).styles(:stroke=>stroke_color, :stroke_width=>1,
            :fill=>color)

          y += increment
        end
        x += block_size
      end
    end
    
  end

   def draw_plot_legend(canvas, x_start, y_start, rotate)

    labels_array = Array.new
      for i in 0..10
        labels_array << i * 0.1
      end

    current_x=0
    current_y=0

    #legend_box_size = @diameter *0.75
    legend_box_size = @diameter * 1.0
    
    if rotate
      increment = -legend_box_size
      text_rotate = 180
      x_text_addition = -legend_box_size*1.5
      y_text_addition = -legend_box_size
    else
      increment = legend_box_size
      text_rotate = 0
      y_text_addition = x_text_addition = 0
    end
     
    # draw boxes from starting point
    canvas.g.translate(x_start,y_start) do |legend|
      labels_array.reverse_each do |box_value|
        color_string = get_correlation_color(box_value)
        legend.rect(legend_box_size, legend_box_size, current_x, current_y).styles(:stroke=>color_string,
          :stroke_width=>1, :fill=>color_string)
        current_y+=increment
      end
    end

    current_y=legend_box_size
    font_size = standard_font_size

    labels_array.reverse_each do |label|
      canvas.g.translate(x_start+x_text_addition, y_start+y_text_addition).text(current_x-2, current_y) do |text|
        text.tspan(label).rotate(text_rotate).styles(:font_size=>font_size/1.3, :text_anchor=>'end')
      end
      current_y+=increment
    end

    # place box around the legend
    current_y = 0 unless rotate
    width = legend_box_size

    canvas.g.translate(x_start, y_start) do |legend|
      legend.rect(width, legend_box_size*11, current_x, current_y).styles(:stroke=>'black',
        :stroke_width=>1, :fill=>'none')
    end

  end


  # draws standard dot with phenotypes along horizontal axis
	def draw_standard_dot(params)
		
		canvas = params[:canvas]
		pheno_list = params[:pheno_list]
		plot_values = params[:plot_values]
		xstart = params[:xstart]
		ystart = params[:ystart]
		xmax = params[:xmax]
		ymax = params[:ymax]
		minpval = params[:minpval]
		maxpval = params[:maxpval]
		rotate = params[:rotate] || false
		draw_triangle = params[:draw_triangle] || false
		plot_beta = params[:plot_beta] || false
		minbeta = params[:minbeta] || 0
		maxbeta = params[:maxbeta] || 1
		labels_on_top = params[:labels_on_top] || false
		plot_sampsizes = params[:plot_sampsizes] || false
		minsize = params[:minsize] || 0
		maxsize = params[:maxsize] || 1000
		draw_lines = params[:draw_lines] || true
		redlinep = params[:redlinep] || nil
		y_best_offset = params[:y_best_offset] || 0
		best_results = params[:best_results] || nil
		coordinates_per_pixel = params[:coordinates_per_pixel]
		narrow_plot = params[:narrow] || false
    tilted_labels = params[:slanted_names] || false
		add_notes = params[:add_notes] || false

		pheno_order = pheno_list.pheno_order
		gene_name_order = pheno_list.gene_name_order
		
		# find p values that are farthest apart
		all_pvals = Array.new
    pheno_order.each_with_index do |phenoname, i|
      plot_values[i].each do |resval|
        if resval.pval < minpval
          next
        end
				all_pvals << resval.pval
			end
		end
		all_pvals.sort!
		pval1=0
		pval2=0
		maxdist = 0
		(1..all_pvals.length-1).each do |i|
			if maxdist < (all_pvals[i]-all_pvals[i-1]).abs
				maxdist = (all_pvals[i]-all_pvals[i-1]).abs
				pval1 = all_pvals[i-1]
				pval2 = all_pvals[i]
			end
		end
		minpscales = Array.new
		maxpscales = Array.new
		pvalintervals = Array.new
		if(maxdist > (maxpval-minpval).to_f/3)
			if pval1 < pval2
				splitmin = pval1
				splitmax = pval2
			else
				splitmin = pval2
				splitmax = pval1
			end
			# split in middle
			minpscales << minpval
			maxpscales << splitmin.ceil
			pvalintervals << maxpscales[0] - minpscales[0]
			minpscales << splitmax.floor
			maxpscales << maxpval
			pvalintervals << maxpscales[1] - minpscales[1]
		else
			minpscales << minpval
			maxpscales << maxpval
			pvalintervals << maxpval - minpval
		end
		
    # determine vertical placement based on distance from top to bottom of plot
    y_interval = @y_plot_height
		y_half_interval = y_interval/2

    pval_interval = maxpval - minpval
    x_circle_start = x_circle = xstart + @offset_mult * @x_offset + @x_interior_offset

    if plot_beta and plot_sampsizes
      y_label_mult = 3
      y_offset_mult = 2
    elsif plot_beta or plot_sampsizes
      y_label_mult = 2
      y_offset_mult = 1.5
    else
      y_label_mult = 1
      y_offset_mult = 1
    end
    
    ystart+=y_best_offset

    if rotate and plot_beta and plot_sampsizes
      y_plots_start = ystart+@y_offset*2+2*y_interval
    elsif rotate and (plot_beta or plot_sampsizes)
      y_plots_start = ystart+@y_offset*1.5+y_interval
    elsif !labels_on_top
      y_plots_start = ystart + @y_offset
    else
        if rotate
          y_plots_start= ystart + @y_offset
        else
          if plot_beta and plot_sampsizes
            y_plots_start = ymax - @y_offset * 2 - @y_plot_height * 3
          else
            y_plots_start = ymax - y_interval * y_label_mult - ystart - y_offset * y_offset_mult + 7 * @diameter/3
          end
        end
    end
    
    if plot_beta or draw_triangle or draw_lines
      x = x_circle_start
      canvas.g.translate(xstart, ystart) do |draw|
        draw.styles(:stroke=>'#CCCCCC', :stroke_width=>1, :opacity=>0.05)
        pheno_order.each do |phenoname|
          draw.line(x, y_plots_start-ystart, x, y_plots_start+y_interval-ystart)
          x=increment_x_phenotype(x)
        end
      end
    end

    # key is color and value is array of y locations
    color_hash = Hash.new
    
    pheno_order.each_with_index do |phenoname, i|
      plot_values[i].each do |resval|
        if resval.pval < minpval
          next
        end

				if(minpscales.length > 1)
					if(resval.pval < maxpscales[0])
						y_point = (1-((resval.pval- minpscales[0]) / pvalintervals[0])) * y_half_interval + y_half_interval
					else
						y_point = (1-((resval.pval- minpscales[1]) / pvalintervals[1])) * y_half_interval * 0.8
					end
				else
					y_point = (1-((resval.pval- minpscales[0]) / pvalintervals[0])) * y_interval
				end
        colorstr = resval.colorstr
        if !color_hash.has_key?(colorstr)
          color_hash[colorstr] = Array.new
        end

        uptriangle = nil

        if draw_triangle and !resval.betaval.nil? and !resval.belowthresh
          if resval.betaval > 0
            uptriangle = false
          else
            uptriangle = true
          end
        end

        color_hash[colorstr] << Point.new(x_circle, y_point, resval.shape, uptriangle)
      end 
      x_circle = increment_x_phenotype(x_circle)
    end

    color_hash.each do |colorstr, points|
      canvas.g.translate(xstart, y_plots_start) do |pen|
        if !draw_triangle
          insert_points(points, colorstr, pen)
        else
          insert_triangles(points, colorstr, pen)
        end
      end
    end

    # start phenotypes at beginning
    x_text = x_circle_start + @diameter/2.5

    if !labels_on_top or rotate
      y_text = y_interval * y_label_mult + y_offset * y_offset_mult + @diameter/3 * y_label_mult
    elsif best_results
      y_text = ymax - y_interval * y_label_mult - ystart - y_offset * y_offset_mult - @diameter/3 - y_best_offset
    else
      y_text = ymax - y_interval * y_label_mult - ystart - y_offset * y_offset_mult - @diameter/3
    end

    # determine font size (needs to be in points and scaled to match the diameter
    # of the circles
    # increase font_size by a multiple
    font_size = standard_font_size

    txt_anchor = 'end'
    txt_anchor = 'start' if labels_on_top and !rotate
    if rotate
      y_txt_best = 0
      best_anchor = 'start'
			txt_rotate=-90
    else
      y_txt_best = ymax-ystart-y_best_offset
      best_anchor = 'end'
			tilted_labels ? txt_rotate=-75 : txt_rotate=-90
    end
    
    x_note = x_circle_start
    note_start=0
    
    # write phenotype labels across bottom of plot
    pheno_order.each_with_index do |phenoname,i|
			if gene_name_order[i].nil?
				canvas.g.translate(xstart, ystart).text(x_text, y_text).rotate(txt_rotate) do |text|
					text.tspan(phenoname).styles(:font_size=>font_size, :font=>Font_phenotype_names,
						:text_anchor =>txt_anchor)
				end
			else
				regexresult = (phenoname =~ /gene_name_order[i]/)
				substrIndex = phenoname.index(gene_name_order[i])
				phenoname_plot = phenoname[gene_name_order[i].length+1..(phenoname.length-1)]
				canvas.g.translate(xstart, ystart).text(x_text, y_text).rotate(txt_rotate) do |text|
					text.tspan(phenoname_plot).styles(:font_size=>font_size, :font=>Font_phenotype_names,
						:text_anchor =>txt_anchor)
				end

				pheno_space = space_for_text(coordinates_per_pixel,phenoname_plot) * 0.98
				# calculate final x,y for line (x2,y2)=(x1+l⋅cos(a),y1+l⋅sin(a))
				 x2 = x_text - pheno_space * Math.cos((360 + txt_rotate) * Math::PI/180)
				 y2 = y_text + pheno_space * Math.sin(-(txt_rotate) * Math::PI/180)
				 
				canvas.g.translate(xstart,ystart).text(x2, y2).rotate(txt_rotate) do |text|
					text.tspan(gene_name_order[i]).styles(:font_size=>font_size, :font=>Font_phenotype_names,
						:text_anchor =>txt_anchor, :font_style=>'italic', :font_weight=>'bold')
				end
      end
      if best_results
        canvas.g.translate(xstart, ystart).text(x_text, y_txt_best).rotate(-90) do |text|
          text.tspan(best_results[phenoname]).styles(:font_size=>font_size, :font=>Font_phenotype_names,
          :text_anchor =>best_anchor)
        end
      end
			
			# add notes if selected
			if add_notes
				# add blocks of notes from bottom up
				# y distance is 1/20 of total height
				y_block_interval = y_interval/20
				y_block_draw = y_interval - y_block_interval*10 
        #y_block_draw = y_block_draw + (note_start*y_block_interval*0.5)
        y_block_start = y_block_draw
				x_block = x_text-@diameter/2.5
				x_block_adj = (x_block - decrement_x_phenotype(x_block)).to_f/2
        alternating = -1
        block_counter = 0
#f_m = 3.9
f_size = font_size/@note_font_multiple
#puts "font size is #{f_size} mult is #{f_m}"
				pheno_list.note_hash[phenoname].each_key do |note|
					x_block_start = x_block - x_block_adj
#					canvas.g.translate(xstart,y_plots_start).rotate(-30).line(x_block_start,y_block_draw,
#						x_block_start+x_block_adj*2, y_block_draw).
#						styles(:stroke=>'darkgray', :opacity=>0.2)
#					canvas.g.translate(xstart, y_plots_start).text(x_text, y_block_draw) do |text|
#						text.tspan(note).styles(:font_size=>f_size, :font=>Font_phenotype_names,
#							:text_anchor=>'start', :font_weight=>'lighter')
#          canvas.g.translate(xstart, y_plots_start).text(x_note, y_block_draw) do |text|
#						text.tspan(note).styles(:font_size=>f_size, :font=>Font_phenotype_names,
#							:text_anchor=>'middle', :font_weight=>'lighter')
          canvas.g.translate(xstart, y_plots_start).text(x_note, y_block_draw).rotate(-30) do |text|
						text.tspan(note).styles(:font_size=>f_size, :font=>Font_phenotype_names,
							:text_anchor=>'middle', :font_weight=>'lighter')
					end
          alternating *= -1
          block_counter += 1
          y_block_draw = y_block_start - (block_counter.to_f/2).ceil * y_block_interval * alternating
            
#					y_block_draw -= y_block_interval
					break if y_block_draw <= 0
				end
			end
      x_note = increment_x_phenotype(x_note)
      x_text = increment_x_phenotype(x_text)
      note_start = 1 - note_start
    end
#exit
    # calculate end positon of box around plot
    x_finish = decrement_x_phenotype(x_text)
    x_end_box = x_finish + @x_interior_offset - xstart
     # if selected draw red line at indicated p value
		if redlinep and (1-((redlinep.to_f) / pval_interval) >= 0.0)
#    if redlinep and (1-((redlinep.to_f) / pval_interval) >= 0.0)
				if(minpscales.length > 1 and (redlinep < maxpscales[0] or redlinep > minpscales[1] ))
					if(redlinep < maxpscales[0])
						y_point = (1-((redlinep- minpscales[0]) / pvalintervals[0])) * y_half_interval + y_half_interval
					else
						y_point = (1-((redlinep- minpscales[1]) / pvalintervals[1])) * y_half_interval * 0.8
					end
				else
					y_point = (1-((redlinep- minpscales[0]) / pvalintervals[0])) * y_interval
				end
		 canvas.g.translate(xstart, y_plots_start) do |draw|
				draw.styles(:stroke=>'red')
				draw.line(@offset_mult*@x_offset,y_point,x_end_box,y_point)
			end
    end

    canvas.g.translate(xstart, y_plots_start) do |draw|
      draw.styles(:stroke_width=>1, :stroke=>'black')
      draw.line(@offset_mult*@x_offset,0,x_end_box,0)
      draw.line(@offset_mult*@x_offset,y_interval,x_end_box,y_interval)

      draw.line(@offset_mult*@x_offset, 0, @offset_mult*@x_offset, y_interval)
      draw.line(x_end_box, 0, x_end_box, y_interval)
    end

		if maxpscales.length == 1
			pval_labels(canvas, minpval, maxpval, 10, xstart,  y_plots_start, offset_mult*@x_offset, y_plots_start+y_interval,
				font_size, "-log10(p value)", rotate)
		else
			pval_labels(canvas, minpscales[0], maxpscales[0], 5, xstart,  y_plots_start+y_half_interval, offset_mult*@x_offset, y_plots_start+y_interval,
				font_size, "-log10(p value)", rotate)
			# draw break lines
      canvas.g.translate(xstart, y_plots_start) do |draw|
        draw.styles(:stroke=>'black', :stroke_width=>2, :opacity=>1.0)
        draw.line(xstart+@offset_mult * @x_offset-@diameter, y_half_interval - y_half_interval*0.125, xstart+@offset_mult * @x_offset+@diameter, y_half_interval - y_half_interval*0.1)
				draw.line(xstart+@offset_mult * @x_offset-@diameter, y_half_interval- y_half_interval*0.1, xstart+@offset_mult * @x_offset+@diameter, y_half_interval - y_half_interval*0.075)
      end
			pval_labels(canvas, minpscales[1], maxpscales[1], 4, xstart,  y_plots_start, offset_mult*@x_offset, y_plots_start+y_half_interval*0.8,
				font_size, " ", rotate)	
		end

    beta_ystart = ystart
    if plot_beta
      if rotate and !labels_on_top
        if !plot_sampsizes
          beta_ystart = ystart+@y_offset
        else
          beta_ystart = ystart+@y_offset*1.5+y_interval
        end
      elsif rotate
        beta_ystart = ystart+@y_offset
        #beta_ystart = ystart+@y_offset*1.5+y_interval
      else
        beta_ystart = y_plots_start + y_interval + @y_offset * 0.75
      end

      draw_beta_plot(canvas, pheno_order, plot_values, xstart, beta_ystart, xmax, y_interval,
        minbeta, maxbeta, minpval, x_end_box, rotate)
    end

    if plot_sampsizes
      if rotate and !labels_on_top
        samp_ystart = ystart+@y_offset
      elsif rotate
        if plot_beta
          samp_ystart = beta_ystart + @y_offset * 0.5 +y_interval
        else
          samp_ystart = ystart+@y_offset
        end
      else
        if plot_beta
          samp_ystart = beta_ystart + y_interval + @y_offset * 0.75
        else
          samp_ystart = y_plots_start + y_interval + @y_offset * 0.75
        end
      end
      draw_sampsize_plot(canvas, pheno_order, plot_values, xstart, samp_ystart, xmax, y_interval, minsize, maxsize,
        minpval, x_end_box, rotate)
    end

  end

  # draws beta plot values
  def draw_beta_plot(canvas, pheno_order, plot_values, xstart, ystart, xmax, y_plot_size, minbeta, maxbeta, minpval,
    x_end_box, rotate=false)
    y_interval = @y_plot_height

    beta_interval = maxbeta - minbeta

    x_circle_start = x_circle = xstart + @offset_mult * @x_offset + @x_interior_offset
    y_plots_start = 0

    x = x_circle_start
    canvas.g.translate(xstart, ystart) do |draw|
      draw.styles(:stroke=>'#CCCCCC', :stroke_width=>1, :opacity=>0.05)
      pheno_order.each do |phenoname|
        draw.line(x, 0, x, y_plot_size)
        x=increment_x_phenotype(x)
      end
    end

    # key is color and value is array of y locations
    color_hash = Hash.new

    pheno_order.each_with_index do |phenoname, i|
      plot_values[i].each do |resval|
        if resval.pval < minpval
          next
        end

        y_point = (1-((resval.betaval-minbeta) / beta_interval)) * y_plot_size

        uptriangle = nil
        colorstr = resval.colorstr
        if !color_hash.has_key?(colorstr)
          color_hash[colorstr] = Array.new
        end
        color_hash[colorstr] << Point.new(x_circle, y_point, resval.shape, uptriangle)
      end
      x_circle = increment_x_phenotype(x_circle)
    end

    # first do any of the first_plot color
    if color_hash.has_key?(first_color)
      canvas.g.translate(xstart, ystart) do |pen|
        insert_points(color_hash[first_color], first_color, pen)
      end
    end


    color_hash.each do |colorstr, points|
      if colorstr != first_color
        canvas.g.translate(xstart, ystart) do |pen|
            insert_points(points, colorstr, pen)
        end
      end
    end

    canvas.g.translate(xstart, ystart) do |draw|
      draw.styles(:stroke_width=>1, :stroke=>'black')
      draw.line(@offset_mult*@x_offset,0,x_end_box,0)
      draw.line(@offset_mult*@x_offset,y_plot_size,x_end_box,y_plot_size)

      draw.line(@offset_mult*@x_offset, 0, @offset_mult*@x_offset, y_plot_size)
      draw.line(x_end_box, 0, x_end_box, y_plot_size)
    end

    # draw line where zero level is on beta plot
    y_point = (1-((0-minbeta) / beta_interval)) * y_plot_size
    dash1 = 2
    dash2 = 4
    dash_array = Array[dash1, dash2]

    canvas.g.translate(xstart, ystart) do |draw|
      draw.styles(:fill=>'none', :stroke_width=>1, :stroke=>'gray', :fill_opacity=>0.0,
        :stroke_dasharray=>dash_array)
      draw.line(@offset_mult*@x_offset, y_point, xmax-xstart-@x_offset, y_point)
    end

    # add labels to side
    font_size = standard_font_size

    pval_labels(canvas, minbeta, maxbeta, 10, xstart,  ystart, 2*@x_offset,
      ystart+y_plot_size, font_size, "beta", rotate)

  end

  # generic version of plot to draw a track for sample size, etc.
  def draw_sampsize_plot(canvas, pheno_order, plot_values, xstart, ystart, xmax, y_plot_size, minsize, maxsize,
    minpval, x_end_box, rotate=false)

    y_interval = @y_plot_height

    sampsize_interval = maxsize - minsize
    x_circle_start = x_circle = xstart + @offset_mult * @x_offset + @x_interior_offset
    y_plots_start = 0

    x = x_circle_start
    canvas.g.translate(xstart, ystart) do |draw|
      draw.styles(:stroke=>'#CCCCCC', :stroke_width=>1, :opacity=>0.05)
      pheno_order.each do |phenoname|
        draw.line(x, 0, x, y_plot_size)
        x=increment_x_phenotype(x)
      end
    end

    # key is color and value is array of y locations
    color_hash = Hash.new

    pheno_order.each_with_index do |phenoname, i|
      plot_values[i].each do |resval|
        if resval.pval < minpval
          next
        end
        y_point = (1-((resval.sampsize-minsize) / sampsize_interval.to_f)) * y_plot_size

        uptriangle = nil
        colorstr = resval.colorstr
        if !color_hash.has_key?(colorstr)
          color_hash[colorstr] = Array.new
        end
        color_hash[colorstr] << Point.new(x_circle, y_point, resval.shape, uptriangle)
      end
      x_circle = increment_x_phenotype(x_circle)
    end

    # first do any of the first_plot color
    if color_hash.has_key?(first_color)
      canvas.g.translate(xstart, ystart) do |pen|
        insert_points(color_hash[first_color], first_color, pen)
      end
    end

    color_hash.each do |colorstr, points|
      if colorstr != first_color
        canvas.g.translate(xstart, ystart) do |pen|
            insert_points(points, colorstr, pen)
        end
      end
    end

    canvas.g.translate(xstart, ystart) do |draw|
      draw.styles(:stroke_width=>1, :stroke=>'black')
      draw.line(@offset_mult*@x_offset,0,x_end_box,0)
      draw.line(@offset_mult*@x_offset,y_plot_size,x_end_box,y_plot_size)

      draw.line(@offset_mult*@x_offset, 0, @offset_mult*@x_offset, y_plot_size)
      draw.line(x_end_box, 0, x_end_box, y_plot_size)
    end

    font_size = standard_font_size

    pval_labels(canvas, minsize, maxsize, 10, xstart,  ystart, 2*@x_offset,
      ystart+y_plot_size, font_size, "sample size", rotate)

  end


  # pass array of points, color and canvas
  def insert_points(points, colorstr, pen)
    pen.styles(:fill=>colorstr, :stroke=>'none', :stroke_width=>1, :fill_opacity=>0.7)
		half_diameter=draw_diameter/2.to_f
    points.each do |point|
			if point.shape == 'square'
				pen.rect(draw_diameter, draw_diameter, point.x-half_diameter, point.y-half_diameter)
			elsif point.shape == 'diamond'
				offset = half_diameter * 1.2
				xpts = [point.x,point.x+offset,point.x,point.x-offset]
				ypts = [point.y-offset,point.y,point.y+offset,point.y]
				pen.polygon(xpts,ypts)
#				pen.polygon(point.x,point.y-half_diameter,point.x+half_diameter,point.y,
#					point.x,point.y-half_diameter,point.x-half_diameter,point.y)
			else  #point.shape == 'circle
				pen.circle(draw_diameter/2, point.x, point.y)
			end
    end
  end

  # draw triangles for values with beta values and circles for those without
  def insert_triangles(points, colorstr, pen)
    circles = Array.new
    pen.styles(:fill=>colorstr, :stroke=>'none', :stroke_width=>1, :fill_opacity=>0.8)

    points.each do |point|
      if point.uptriangle.nil?
        circles << point
      else

        xpoints = Array.new
        ypoints = Array.new

        xpoints << point.x-@diameter/2
        xpoints << point.x+@diameter/2
        xpoints << point.x

        y_adjust = @diameter
        if !point.uptriangle
          y_adjust = -y_adjust
        end
        ypoints << point.y
        ypoints << point.y
        ypoints << point.y+y_adjust

        pen.polygon(xpoints,ypoints)

      end
    end

    if !circles.empty?
      insert_points(circles, colorstr, pen)
    end

  end


  # draw and label p values along the vertical axis
  def pval_labels(canvas, min, max, num_intervals, xstart, ystart, xend, yend, font_size, plot_title,rotate)
    pval_interval = max-min
    stat_break = pval_interval.to_f/num_intervals
    y_interval = yend - ystart
    y_break = y_interval.to_f/num_intervals
    num_x = @x_offset * @offset_mult - @x_offset * 0.5
    current_stat_value = min
    precision=1

    num_intervals+=1

    if rotate
      rotation = -180
      anchor = 'start'
      y=y_interval - @diameter/4
    else
      rotation = 0
      anchor = 'end'
      y=y_interval+@diameter/4
    end

    max_label_length=0
    
    num_intervals.times do |curr_break|
      canvas.g.translate(xstart,ystart).text(num_x, y).rotate(rotation) do |text|
        label = compress_number(current_stat_value, precision)
        max_label_length = label.length if label.length > max_label_length
        text.tspan(label).styles(:font_size=>font_size/1.4, :text_anchor=>anchor)
        current_stat_value += stat_break
        y -= y_break
      end
    end

    y=y_interval

    num_intervals.times do |curr_break|
      canvas.g.translate(xstart,ystart) do |draw|
        draw.styles(:stroke_width=>1, :stroke=>'black')
        draw.line(@offset_mult*@x_offset-@x_offset/4,y,@offset_mult*@x_offset,y)
        y -= y_break
      end
    end

#    dist_mult = max_label_length.to_i/2 + 1.25
		dist_mult = 3.25
    # add title for this portion of plot
    canvas.g.translate(xstart,ystart).text(num_x-@diameter*dist_mult,y_interval/2).rotate(-90) do |text|
      text.tspan(plot_title).styles(:font_size=>font_size, :text_anchor=>'middle')
    end

  end


  def compress_number(num, precision=2)

    if precision > 0
      num_string = sprintf("%0.#{precision}f", num)
    else
      num_string = sprintf("%d", num)
    end
    if(num_string =~ /^00$/)
      num_string = "0  "
    end

    return num_string
  end

  # moves x to next location along plot
  def increment_x_phenotype(x)
#    return x + @diameter * 2
		return x + @diameter * @pheno_dist_mult
  end

  def decrement_x_phenotype(x)
#    return x - @diameter * 2
		return x - @diameter * @pheno_dist_mult
  end

end


def draw_phewas(options)
  # need to read in phenotype file to determine
  # number of phenotypes for proper sizing of the plot
  resultholder = ResultHolder.new

  if options.phenotype_listfile
    phenolistreader = PhenoListReader.new
    phenolisthash = phenolistreader.read_file(options.phenotype_listfile)
    resultholder.set_included_phenos(phenolisthash)
  else
    resultholder.include_all_phenos = true
  end

  ethmap = EthMap.new
  ethmap.set_restricted_eths(options.ethlist)
  ethmap.include_all = options.include_all_eths

  if options.ethmapfile
    ethreader = EthMapReader.new
    ethreader.read_file(ethmap, options.ethmapfile)
		if options.group_symbols
			ethmap.set_shapes
		end
  end

  resultholder.ethmap = ethmap

  if options.phenofile
    phenoreader = ExpectedPhenoReader.new
    phenoreader.read_file(resultholder, options.phenofile)
    resultholder.unexpectedcolor = 'purple'
    resultholder.expectedcolor = 'blue'
  end

  resultholder.single_snp = options.snpid
  phenoconverthash = Hash.new

  resreader = ResultFileReader.new
  if options.phenofile
    grouprequired = true
  else
    grouprequired = false
  end
  resreader.read_file(:resultholder=>resultholder, :filename=>options.phewasfile, 
    :grouprequired=>grouprequired, :samprequired=>options.plot_sampsizes, 
    :classname=>options.classname, :include_genename=>options.include_genename)

  # set up titles for best results
  resultholder.set_best_values if options.showbest
    
  # read the correlation matrix if included
  if options.phenotype_correlations_file
    cor_reader = CorrelationReader.new
    cor_reader.read_file(options.phenotype_correlations_file, resultholder)
  end

  # use synthesis_view like approach to sizing the plot
  # based on number of phenotypes in the list
  total_phenotypes = resultholder.pheno_list.pheno_order.length

  diameter_size = case
                when total_phenotypes > 200 then 16
                when total_phenotypes > 160 then 18
                when total_phenotypes > 120 then 20
                when total_phenotypes > 80 then 22
                else 24
                end

  dotter = DotPlotter.new
  dotter.diameter = diameter_size
  dotter.font_size_multiple = 1.25 if options.largetext
  pheno_narrow_multiplier = 1.0
	if options.narrow
    # no notes to include
    if resultholder.pheno_list.note_hash.empty?
      dotter.font_size_multiple = 0.85
      dotter.pheno_dist_mult = 1.2
    else
      # notes to include
      dotter.pheno_dist_mult = 2.0
      dotter.font_size_multiple = 1.0
      dotter.note_font_multiple = 3.9
      pheno_narrow_multiplier = 0.97
      options.largetext=false
#      exit
    end
	elsif !resultholder.pheno_list.note_hash.empty?
		dotter.pheno_dist_mult = 5.0
		dotter.font_size_multiple = 1.2
    pheno_narrow_multiplier = 1.2
    options.largetext=false
	end

  xside_end_addition = 0.005 * dotter.diameter * 2
  xside_end_addition *= 2.5 if dotter.diameter > 20

  # determine the left side of the plot
  xleft_addition = 0.001 * dotter.diameter * 10 + 0.012 * dotter.diameter
  xside_end_addition += xleft_addition
  xside = dotter.calculate_plot_width(total_phenotypes, xside_end_addition)
  xmax = dotter.calculate_coordinate(xside)

  sidebuffer = dotter.diameter * 3
  dotter.x_interior_offset = dotter.diameter * 2
  dotter.x_offset = dotter.diameter*2
  dotter.y_pval_zero = dotter.diameter * 50

  # calculate vertical size (need to account for length of phenotype label)
  dotter.y_offset =  dotter.diameter * 2
  dotter.y_plot_height = dotter.diameter * 18

  ymax = 0
  yside = 0

  # add offset for title
  y_title_start = ymax
  title_fract = 0.008
  title_fract = 0.008 if options.rotate
  yside = dotter.add_space_for_title(yside, title_fract)
  ymax = dotter.calculate_coordinate(yside)
	y_title_start = ymax

	 if (!resultholder.single_snp.nil? or options.group_colors) and options.rotate
    yside = yside + dotter.diameter * 0.00355
  end

	dotter.x_offset = sidebuffer
	legend_y_space=0
	# add space for legend when rotated (add to Y)
	if (!resultholder.single_snp.nil? or options.group_colors) 
			# need max label for ethnicity
			max_label=""
			edata = Array.new
			resultholder.ethmap.ethindata.each do |ethname, tf|
				if tf
					max_label = ethname if ethname.length > max_label.length
					edata.push(ethname)
				end
			end
		if options.rotate
			legend_y_space = dotter.space_for_text(ymax/yside.in, max_label)
#			y_pval_start += legend_y_space
			yside += legend_y_space/(ymax/yside.in) / RVG::dpi
			ymax += legend_y_space
		else
			legend_x_space = dotter.space_for_text(xmax/xside.in, max_label)
			legend_x_start = dotter.offset_mult * dotter.x_offset + dotter.x_interior_offset
			# determine end
			distance_increment=dotter.increment_x_phenotype(0)
			legend_x_end = resultholder.pheno_list.pheno_order.length * distance_increment
			legend_x_space = legend_x_space + dotter.space_for_text(xmax/xside.in, " ") + dotter.diameter;
			
			# per labels per row 
			labels_per_row = (legend_x_end / legend_x_space).to_i
			# if labels_per_row < 1 extend the right side of plot
			if labels_per_row < 1
				labels_per_row = 1
				additional_x = legend_x_space - legend_x_end 
				xside += additional_x/(xmax/xside.in) / RVG::dpi
				xmax += additional_x
			end
			nlabel_rows = edata.length / labels_per_row.to_f
			if((nlabel_rows - nlabel_rows.to_i) > 0.000000000001)
				nlabel_rows +=1
			end
			nlabel_rows = nlabel_rows.to_i
			legend_y_row =  dotter.diameter * 1.25
			legend_y_space = (nlabel_rows-1) * legend_y_row
#			y_pval_start += legend_y_space
			yside += legend_y_space/(ymax/yside.in) / RVG::dpi
			ymax += legend_y_space		
		end
	end
	
  # add plot height for beta
  plot_height_mult=1
  plot_offset_mult=1
  if options.beta
    plot_height_mult += 1
    plot_offset_mult += 0.5
  end
  # add plot height for sample sizes
  if options.plot_sampsizes
    plot_height_mult +=1
    plot_offset_mult += 0.5
  end

  # exit
  y_pval_start = ymax

  yside = yside + dotter.diameter * 2 * 0.0355 * plot_height_mult
  ymax = dotter.calculate_coordinate(yside)
  

  # add room for the labels across bottom of plot
  max_pheno_length = resultholder.get_longest_label
  options.largetext ? text_multiplier = 0.0020 * dotter.font_size_multiple : text_multiplier = 0.0020
  if total_phenotypes < 10
    text_multiplier = 0.0025
  elsif total_phenotypes < 20
    text_multiplier = 0.0024
  elsif total_phenotypes < 30
    text_multiplier = 0.0023
  elsif total_phenotypes < 40
    text_multiplier = 0.0022
  end
  
  yside += (text_multiplier* dotter.diameter * max_pheno_length* pheno_narrow_multiplier) + 0.012 * dotter.diameter
  ymax = dotter.calculate_coordinate(yside)

  ybest_offset=0
  # add space for titles for best
  if options.showbest
    yside += (text_multiplier*dotter.diameter*resultholder.max_best_title)+0.012 * dotter.diameter
    ynewmax = dotter.calculate_coordinate(yside)
    ybest_offset = ynewmax-ymax
    ymax = ynewmax
  end
  
 
  # add space for the correlation matrix boxes if correlation information
  # was provided
  ygrid_start = ymax
  if options.phenotype_correlations_file
    max_correlations = resultholder.get_max_correlations
    yside = dotter.add_grid_size(yside, xside, max_correlations-1, resultholder.pheno_list.pheno_order.length)
    ymax = dotter.calculate_coordinate(yside)
  end
	
  # set first color to be plotted
  dotter.first_color = resultholder.nonsigcolor

  rvg = RVG.new(xside.in, yside.in).viewbox(0,0,xmax,ymax) do |canvas|
    canvas.background_fill = 'rgb(255,255,255)'

#    dotter.x_offset = sidebuffer
    pval_threshold = 0
    if options.p_thresh > 0
      pval_threshold = resultholder.get_log10(options.p_thresh)
    end

    plot_values = resultholder.generate_result_values(pval_threshold, :group_colors=>options.group_colors)

    if options.rotate
      dotter.draw_title(canvas, 0, 0, ymax, y_pval_start+dotter.diameter/4-legend_y_space, options.title, options.rotate)
    else
      dotter.draw_title(canvas, 0, 0, xmax, y_title_start+dotter.diameter/4, options.title)
    end

    if options.maxp_to_plot.to_f < 1.0
      resultholder.minpval = resultholder.get_log10(options.maxp_to_plot.to_f)
    end

    rotate_grid_offset = dotter.diameter*3.0
    # add legend when only single SNP plotted
    if !resultholder.single_snp.nil? or options.group_colors
      # only draw legend including ethnicities and colors selected
      ethcolorhash = Hash.new
			ethshapehash = Hash.new
			ethorder = Array.new
      resultholder.ethmap.ethindata.each do |ethname, tf|
        if tf
          ethcolorhash[ethname] = resultholder.ethmap.eths[ethname].colorstr
					ethshapehash[ethname] = resultholder.ethmap.eths[ethname].shape
					ethorder.push(ethname)
        end
      end
      rotate_grid_offset = 0
			#labels_per_row legend_y_row
			if options.rotate
				dotter.draw_legend(canvas, ethcolorhash, ethshapehash,0, y_pval_start, :rotate=>options.rotate,
					:coordinates_per_pixel=>xmax/xside.in, :max_label=>max_label)
			else
				current_row=nlabel_rows - 1
				current_eth=0
				while current_row >= 0
					tempEths = Hash.new
					tempShapes = Hash.new
					ethIndx=0
					while(ethIndx < labels_per_row and ethIndx < ethorder.length)
						tempEths[ethorder[current_eth]]=ethcolorhash[ethorder[current_eth]]
						tempShapes[ethorder[current_eth]]=ethshapehash[ethorder[current_eth]]
						ethIndx+=1
						current_eth+=1
					end
					dotter.draw_legend(canvas, tempEths, ethshapehash,0, y_pval_start-(legend_y_row*current_row), :rotate=>options.rotate,
						:coordinates_per_pixel=>xmax/xside.in, :max_label=>max_label)
					current_row -=1
				end
			end
    end

    options.redline = resultholder.get_log10(options.redline.to_f) if options.redline
  
    if options.rotate and options.phenotype_correlations_file
      ygrid_orig = ygrid_start
      ygrid_start = y_pval_start + rotate_grid_offset
      y_pval_start = ymax- ygrid_orig + y_pval_start
    end

	
		dotter.draw_standard_dot(:canvas=>canvas, :pheno_list=>resultholder.pheno_list,
			:plot_values=>plot_values, :xstart=>0, :ystart=>y_pval_start, :xmax=>xmax,
			:ymax=>ygrid_start, :minpval=>resultholder.minpval, :maxpval=>resultholder.maxpval.ceil,
			:rotate=>options.rotate, :draw_triangle=>(!resultholder.single_snp.nil? or options.show_direction), 
			:plot_beta=>options.beta, :minbeta=>resultholder.minbeta, :maxbeta=>resultholder.maxbeta, 
			:labels_ton_top=>options.phenotype_correlations_file, :plot_sampsizes=>options.plot_sampsizes,
			:minsize=>resultholder.minsampsize, :maxsize=>resultholder.maxsampsize, 
			:draw_lines=>!options.nolines, :redlinep=>options.redline, :y_best_offset=>ybest_offset,
			:best_results=>resultholder.best_results, :coordinates_per_pixel=>ymax/yside.in,
			:narrow=>options.narrow, :add_notes=>!resultholder.pheno_list.note_hash.empty?, 
      :slanted_names=>(options.narrow || options.tilted_text))
    
    if options.phenotype_correlations_file #and !options.rotate
      dotter.draw_grid(canvas, resultholder, 0, ygrid_start, options.rotate)
    end
  end

  # produce output file
  outfile = options.out_name + '.' + options.imageformat
  print "\n\tDrawing #{outfile}..."

  STDOUT.flush
  begin 
    img = rvg.draw
  rescue => e
    puts "ERROR: #{e.to_s}"
    exit!
  end
  if options.rotate
    img.rotate!(90)
  end
  img.write(outfile){self.quality = 80}
  print " Created #{outfile}\n\n"
end



def draw_sun(options)
  resultholder = SunResultHolder.new
  
  if options.snpid =~ /\w/
    resultholder.center_name = options.snpid
    resultholder.phenorequired = false
    resultholder.appendgene = options.labelgene
  elsif options.phenoname =~ /\w/
    resultholder.center_name = options.phenoname
    resultholder.phenorequired = true
  else options.genename =~ /\w/
    resultholder.center_name = options.genename
    resultholder.phenorequired = false
    resultholder.appendgene = options.labelgene
  end

  resultholder.appendeth = options.include_eth
  
  resreader = SunResultFileReader.new

  resreader.read_file(resultholder, options.phewasfile, options.maxp_to_plot.to_f, 
    options.label_cols)

  # need to figure out height and width for plot
  # width is more straightforward
  # 2r for the circle plus 10r for the lines
  total_values = resultholder.total_values

  # need to calculate a radius size based on the number of values to plot
  radius_size = case
              when total_values > 100 then 20
              when total_values > 80 then 28
              when total_values > 60 then 28
              when total_values > 40 then 28
              when total_values > 20 then 28
              when total_values > 10 then 28
              else 28
              end

  plotter = RadialPlotter.new
  plotter.line_size = radius_size 
  yside_mult = 1.45
  title_height_mult = 1

  if total_values > 70
    total_radii = 44
    plotter.radius = 8
    if options.largetext
      yside_mult = 1.2
      plotter.font_size_multiple *= 1.8
      title_size_fit = 63
      title_adjustment = 0.009 * plotter.font_size_multiple
    else
      plotter.font_size_multiple = 1.2
      title_size_fit = 88
      title_adjustment = 0.015 * plotter.font_size_multiple
    end
  elsif total_values > 50
    total_radii = 34
    plotter.radius = 6
    title_size_fit = 75
    if options.largetext
      yside_mult = 1.22
      plotter.font_size_multiple = 1.65
      title_size_fit = 53
      title_adjustment = 0.014 * plotter.font_size_multiple
    else
      plotter.font_size_multiple = 1.1
      title_size_fit = 80
      title_adjustment = 0.017 * plotter.font_size_multiple
    end
  elsif total_values > 40
    total_radii = 20
    plotter.radius = 4
    if options.largetext
      yside_mult = 1.25
      plotter.font_size_multiple = 1.5
      title_size_fit = 34
      title_adjustment = 0.016 * plotter.font_size_multiple
    else
      plotter.font_size_multiple = 1.05
      title_size_fit = 48
      title_adjustment = 0.019 * plotter.font_size_multiple
    end
  elsif total_values > 30
    total_radii = 20
    plotter.radius = 4
    if options.largetext
			title_height_mult = 1.05
      yside_mult = 1.25
      plotter.font_size_multiple = 1.9
			plotter.center_size_multiple = 0.9
			plotter.title_text_adj = 0.8
      title_size_fit = 34
      title_adjustment = 0.016 * plotter.font_size_multiple
    else
      plotter.font_size_multiple = 1.05
      title_size_fit = 48
      title_adjustment = 0.019 * plotter.font_size_multiple
    end
  elsif total_values > 20
    total_radii = 18
    plotter.radius = 4
    title_size_fit = 50
    if options.largetext
      title_height_mult = 1.15
      yside_mult = 1.3
      plotter.font_size_multiple = 2.1
			plotter.center_size_multiple = 0.9
			plotter.title_text_adj = 0.75
      title_size_fit = 27
      title_adjustment = 0.017 * plotter.font_size_multiple
    else
      plotter.font_size_multiple = 1.0
      title_size_fit = 43
      title_adjustment = 0.020 * plotter.font_size_multiple
    end
  elsif total_values > 12
    total_radii = 12
    plotter.radius = 2
    if options.largetext
      title_height_mult = 1.25
      yside_mult = 1.35
      plotter.font_size_multiple = 2.1
			plotter.center_size_multiple = 0.7
			plotter.title_text_adj = 0.6
      title_size_fit = 20
      title_adjustment = 0.027 * plotter.font_size_multiple
    else
      plotter.font_size_multiple = 1.0
      title_size_fit = 30
      title_adjustment = 0.029 * plotter.font_size_multiple
    end
    title_size_fit = 25
  else
#    total_radii = 28 # 12 is previous value
    plotter.radius = 1.5 # 2 is original
    if options.largetext
			total_radii = 20
#      title_height_mult = 1.6
      yside_mult = 1.4
			options.title =~ /\w/ ? title_height_mult = 1.6 : title_height_mult = 0.35
      plotter.font_size_multiple = 1.8
			plotter.center_size_multiple = 0.59 # 0.5 is original
			plotter.title_text_adj = 0.5
      title_size_fit = 35
      title_adjustment = 0.028 * plotter.font_size_multiple
    else
			total_radii=14
			options.title =~ /\w/ ? title_height_mult = 0.8 : title_height_mult = 0.35
			plotter.radius=1.5
      plotter.font_size_multiple = 1.1
			plotter.center_size_multiple = 0.8
			plotter.title_text_adj = 0.8
      title_size_fit = 35
      title_adjustment = 0.032 * plotter.font_size_multiple
    end
  end

title_size_fit = 28

if total_values > 30
	plotter.font_size_multiple = 2.0
elsif total_values > 24
	plotter.font_size_multiple = 2.0
elsif total_values > 18
	plotter.font_size_multiple = 2.0
	title_size_fit = 28
elsif total_values > 12
	plotter.font_size_multiple = 2.0
else
	plotter.font_size_multiple = 2.8
	title_size_fit = 28
end	
	
plotter.radius = 2
total_radii = 20
#      title_height_mult = 1.6
yside_mult = 1.4
options.title =~ /\w/ ? title_height_mult = 1.6 : title_height_mult = 0.4
#plotter.font_size_multiple = 1.8
plotter.center_size_multiple = 0.59 # 0.5 is original
plotter.title_text_adj = 0.5

title_adjustment = 0.028 * plotter.font_size_multiple	
	
  xside_end_addition = 0
  resultholder.max_title_length += 2 if options.beta
  resultholder.max_title_length += 3 if options.include_eth
	resultholder.max_title_length += 7 # 0.123 or 1e281

  # add some space on both sides outside of plot
  if resultholder.max_title_length > title_size_fit
    xside_end_addition = title_adjustment * plotter.radius * (resultholder.max_title_length - title_size_fit)
  end

  xside = plotter.calculate_plot_width(total_radii, xside_end_addition)
  xmax = plotter.calculate_coordinate(xside)
  ymax = 0
  yside = 0

  # set space for the main title
  y_title_start = ymax
  yside = plotter.add_space_for_title(yside, 0.016*title_height_mult)
  ymax = plotter.calculate_coordinate(yside);

  y_plot_start = ymax
  # add vertical space for the main circle plot (which should be the same as plot width)
  yside += plotter.calculate_plot_width(total_radii, 0)/yside_mult
  ymax += plotter.calculate_coordinate(yside)
	
  rvg = RVG.new(xside.in, yside.in).viewbox(0,0,xmax,ymax) do |canvas|
    canvas.background_fill = 'rgb(255,255,255)'
    plotter.canvas = canvas

    # for now just write title without splitting line
    plotter.write_main_title(options.title, 0, y_title_start, xmax, y_plot_start, false)

    p_val_thresh = resultholder.results.get_neg_log(options.p_thresh)

    plotter.draw_main_plot(:center_name=>resultholder.center_name, :reslist=>resultholder.results,
      :midx=>xmax/2, :ystart=>y_plot_start, :yend=>ymax, :p_thresh=>p_val_thresh, 
      :total_radii=>total_radii, :use_beta=>options.beta, :gene_name=>resultholder.gene_name,
      :largetext=>options.largetext, :use_pval=>options.display_pval, :bold_font=>options.boldfont)
  end #end RVG.new

  # produce output file
  outfile = options.out_name + '.' + options.imageformat
  print "\n\tDrawing #{outfile}..."
  STDOUT.flush
  begin 
    img = rvg.draw
  rescue => e
    puts "ERROR: #{e.to_s}"
    exit!
  end
  img.write(outfile)
  print " Created #{outfile}\n\n"
#  
#  smallfile = options.out_name + '.small.' + options.imageformat
#  smallimg = img.scale(0.5)
#  smallimg.write(smallfile)
#  print " Created #{smallfile}\n\n"
  
end


#######################################
#######################################
#######################################
#######################################
# Main execution begins here
options = Arg.parse(ARGV)

if options.lowres
  RVG::dpi = 72
end
 
if options.sun_file
  draw_sun(options)
else
  draw_phewas(options)
end

