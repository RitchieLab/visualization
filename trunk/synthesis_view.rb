#!/usr/bin/env ruby

###################################################################################
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. <http://www.gnu.org/licenses/>
###################################################################################

###################################################################################
# Creates plots using RMagick that integrate multiple pieces of information such
# as p values, beta (effect sizes) and mean allele frequencies.  Displays SNPs by
# chromosome and position and allows for display of gene information to show
# features near SNPs.  Also can display LD information in a Haploview-like format.
###################################################################################

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
  puts "\nPlease install RMagick -- See documentation for synthesis_view or http://rmagick.rubyforge.org/install-faq.html "
  puts
  exit(1)
end

require 'optparse'
require 'ostruct'
include Magick

RVG::dpi=72
Version = '1.05'
Maximum_html_image_x = 1000

# check for windows and select alternate font based on OS
Font_family_style = RUBY_PLATFORM =~ /mswin/ ? "Verdana" : "Times"


#############################################################################
#
# Class Arg -- Parses the command-line arguments.
#
#############################################################################
class Arg

  def self.init_options
    options = OpenStruct.new

    options.rsquared = nil
    options.dprime = nil
    options.eaglefile = nil
    options.out_name = 'synthesis_view'
    options.highres = nil
    options.abbrevfile = nil
    options.imageformat = 'png'
    options.groupfile = nil
    options.ld_file = nil
    options.genefile = nil
    options.title = "-log10(p-value) plot"
    options.largetext = false
    options.p_thresh = 0.0
    options.phenotitle = "Avg phenotype"
    options.pmin = nil
    options.effect_name = "beta"
    options.no_triangles = false
    options.rotate = false
    options.ormin_zero = false
    options.casecontrol = false
    options.circlecasecon = false
    options.cafcasecontrol = false
    options.circlecafcasecontrol = false
    options.clean_axes = false
    options.plot_beta = false
    options.jitter = false
    options.plot_pval = false
    options.plot_studynum = false
    options.large_or = false
    options.option_logfile = nil
    options.forest_legend =false
    options.highlighted_group = ""
    options.htmlfile = false
    options.add_snp_locs = false
    options.manhattanfile = nil
    options.manhattanthresh = nil
    options.include_dist_track = false
    options.plot_power = false
    options.grayscale = false
    options.additional_columns = Array.new

    return options
  end

  def self.parse(args, options=nil)

    if options == nil
      options=init_options
    end

    options.eaglefile = nil
    options.groupfile = nil
    options.ld_file = nil
    options.genefile = nil
    options.option_logfile = nil

    help_selected = false
    version_selected = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: synthesis_view.rb [options]"

      opts.on_tail("-h", "--help", "Show this usage statement") do
        puts opts
        help_selected = true
      end

      opts.on_tail("-v", "--version", "Show version") do
        puts "\n\tVersion: #{Version}"
        version_selected = true
      end

      opts.on("-e [synth_file]", "Synthesis-View format file for input") do |synth_file|
        options.eaglefile = synth_file
      end

      opts.on("-G [log_file_name]", "Options log file") do |log_file_name|
        options.option_logfile = log_file_name
      end

      opts.on("-s [phenotypic_summary_file]", "Group summary file for input") do |phenotypic_summary_file|
        options.groupfile = phenotypic_summary_file
      end
      opts.on("-g [gene_file]", "Gene summary file") do |gene_file|
        options.genefile = gene_file
      end

      opts.on("-l [ldfile]", "Linkage disequilibrium file") do |ldfile|
        options.ld_file = ldfile
      end

      opts.on("-b [abbrev_file]", "Abbreviation definition file") do |abbrev_file|
        options.abbrevfile = abbrev_file
      end

      opts.on("-t [title_str]", "Main title for plot (enclose in quotes)") do |title_str|
        if title_str 
          options.title = title_str
        else
          options.title = "-log10(p-value) plot"
        end
      end

      opts.on("-X [manhattan_file]", "Manhattan plot file") do |manhattan_file|
        options.manhattanfile = manhattan_file
      end
      opts.on("-Q [manhattan_thresh]", "Manhattan threshold line value (p-value)") do |manhattan_thresh|
        options.manhattanthresh = manhattan_thresh
      end

      opts.on("-M", "--more-locs", "Label SNP locations at additional spots along chromosome") do |more_locs|
        options.add_snp_locs = true
      end

      opts.on("-E", "--dist-track", "Includes the distance track in the plot") do |dist|
        options.include_dist_track = true
      end

      opts.on("-T", "Larger plot font") do |tsize|
        options.largetext = true
      end
      opts.on("-S", "--cleaner-axis", "Axes on plots \"cleaner\".  Default uses min and max values as axis boundaries") do |clean|
        options.clean_axes = true
      end

      opts.on("-J", "--jitter", "Offsets results from different groups horizontally") do |jitter|
        options.jitter = true
      end
      opts.on("-c", "--no-triangle", "P values will be drawn as circles instead of triangles that show direction of effect") do |c|
        options.no_triangles = true
      end
      opts.on("-U [highlighted]", "Group will be drawn as diamonds") do |highlighted|
        options.highlighted_group = highlighted.upcase
      end
      opts.on("-A", "--pval", "Plot p value track") do |pval|
        options.plot_pval = true
      end

      opts.on("-B", "--beta", "Plot effect size track") do |beta|
        options.plot_beta = true
      end

      opts.on("-u [columns]", Array, "Columns to plot that are not standard") do |columns|
        options.additional_columns = columns
      end
      
      opts.on("-n [eff_name]", "Effect size name (beta or es)") do |eff_name|
        options.effect_name = eff_name
      end
      opts.on("-p [pthresh]", "p value threshold") do |pthresh|
        options.p_thresh = pthresh.to_f
      end
      opts.on("-P [pheno_title]", "Phenotype title for phenotype summary plot (enclose in quotes)") do |pheno_title|
        if pheno_title and pheno_title =~ /\w/
          options.phenotitle = pheno_title
        else
          options.phenotitle = "Avg phenotype"
        end
      end
      opts.on("-m", "--pmin [pmin_val]", "Maximum y-axis setting for p-value track") do |pmin_val|
        options.pmin = pmin_val
      end
      opts.on("-N", "--plot-study-total", "Plot number of studies in meta analyis track") do |plot|
        options.plot_studynum = true
      end
      opts.on("-d", "Draws d-prime result plot") do |d|
        options.dprime = true
      end
      opts.on("-r", "Draws r-squared result plot") do |r|
        options.rsquared = true
      end
      opts.on("-a", "--rotate", "Draw Forest plot") do |rot|
        options.rotate = true
      end
      opts.on("-D", "--forest-legend", "Draw legend on Forest plot") do |for_legend|
        options.forest_legend = true
      end
      opts.on("-W", "--power-plot", "Draw power plot") do |pow_plot|
        options.plot_power = true
      end
      opts.on("-L", "--large-or", "Significant odds ratios larger in size on plot") do |large_or|
        options.large_or = true
      end
      opts.on("-R", "--or-zero", "Set minimum value on Forest plot x-axis to zero") do |orzero|
        options.ormin_zero = true
      end
      opts.on("-C", "--case-control", "Plot case and control totals on separate plots") do |cascon|
        options.casecontrol= true
      end
      opts.on("-O", "--circle-case-con", "Plots cases and controls as closed and open circles on same track") do |ccon|
        options.circlecasecon = true
      end
      opts.on("-H", "--caf-case-control", "Plots case and control caf on separate plots") do |cafcascon|
        options.cafcasecontrol = true
      end
      opts.on("-K", "--caf-circle-case-control", "Plots case and control caf as closed and open circles on same track") do |cafcascon|
        options.circlecafcasecontrol = true
      end
      opts.on("-z", "Draws high resolution figure (300dpi)") do |hr|
        options.highres = true
      end
      opts.on("-Y", "Create grayscale image") do |gray|
        options.grayscale = true
      end
      opts.on("-w", "Produce html file with links to dbSNP") do |www|
        options.htmlfile = true
      end
      opts.on("-f [image_type]", "Image format for output (png default)") do |image_type|
        options.imageformat = image_type
      end
      opts.on("-o [output_name]", "Optional output name for the image") do |output_file|
        options.out_name = output_file
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

    if !options.eaglefile
      help_string = opts.help
      puts "\n",help_string,"\n"
      puts "\nExample: synthesis_view.rb -e OVERLAPPING_SNPS_LipidStudy.txt -o lipid_study\n\n"
      exit(1)
    end

    return options

  end

  # write options out to a file
  def self.write_options(options, outfilename)

    options_dump = options.marshal_dump

    File.open(outfilename, 'w') do |fwrite|
      options_dump.each do |key, value|
        fwrite.puts "#{key}\t#{value}\n"
      end
    end
  end

  # read options from file
  def self.read_options(infilename)
    options = OpenStruct.new
    options_dump = Hash.new
    begin
      File.open(infilename, "r") do |file|

      file.each_line do |line|
        line.strip!
        data = line.split(/\t/)
        if data[1] =~ /^true$/
          data[1] = true
        elsif data[1] =~ /^false$/
          data[1] = false
        elsif data[1] =~ /^\d+\.*\d*$/
          data[1] = data[1].to_f
        end

        options_dump[data[0]] = data[1]
      end

     end
    rescue Exception => e
      puts e
      exit(1)
    end

    options.marshal_load(options_dump.inject({ }) { |h, (k,v)| h[k.to_sym] = v; h })

  return options

  end

end


#############################################################################
#
# Class ChromosomeList -- Holds all chromosomes.
#
#############################################################################
class ChromosomeList
  attr_accessor :chromhash, :chromarray, #:maxpscore, :minpscore, :minbeta, :maxbeta,
    :maxscore, :minscore, :genes

  def initialize
    @chromhash = Hash.new
    @chromarray = Array.new
    @maxscore = Hash.new
    @minscore = Hash.new
    @maxscore['pvalue'] = 0.0
    @minscore['pvalue'] = 1.0
    @maxscore['beta'] = -10000
    @minscore['beta'] = 10000
    @maxscore['maf'] = 0.5
    @minscore['maf'] = 0.0
    @maxscore['N'] = 0
    @minscore['N'] = 1000000
    @minscore['or'] = 10.0
    @maxscore['or'] = 0.0
    @maxscore['uci'] = 0.0
    @minscore['uci'] = 10.0
    @maxscore['lci'] = 0.0
    @minscore['lci'] = 10.0
    @maxscore['cases'] = 0
    @minscore['cases'] = 10000000
    @maxscore['controls'] = 0
    @minscore['controls'] = 10000000
    @maxscore['cafcases'] = 0.0
    @minscore['cafcases'] = 1.0
    @maxscore['cafcontrols'] = 0.0
    @minscore['cafcontrols'] = 1.0
    @maxscore['study'] = 0
    @minscore['study'] = 10000000
    @maxscore['betauci'] = -100000
    @minscore['betauci'] = 1000000
    @minscore['betalci'] = 1000000
    @maxscore['betalci'] = -100000
  end

  
  def set_columns(add_columns)
    add_columns.each do |column|
      @minscore[column]=1e20
      @maxscore[column]=-1e20
    end
  end
  
  # returns length of longest SNP name in chromosome list
  def get_max_snp_name
    max_name = 0
    chromhash.each_value do |chrom|
      if chrom.get_max_snp_name > max_name
        max_name = chrom.get_max_snp_name
      end
    end
    return max_name
  end

  def add_chrom(chrom)
    @chromhash[chrom.chrom_num] = chrom
    @chromarray << chrom.chrom_num
  end

  def get_chrom(chromnum)
    return @chromhash[chromnum]
  end

  def get_nsnps
    total_snps=0
    @chromarray.each do |num|
      total_snps = total_snps + @chromhash[num].snp_list.snps.length
    end
    return total_snps
  end

  # return maximum number of snps in any chromosome
  def get_max_included_snps
    max_snps = 0
    @chromarray.each do |num|
      if @chromhash[num].snp_list.get_num_included_snps > max_snps
        max_snps = @chromhash[num].snp_list.get_num_included_snps
      end
    end
    return max_snps
  end

  def sort_chroms!
    @chromarray.sort!{|x,y| x <=> y}
  end

  def reverse!
    @chromarray.reverse!
  end

  def get_max_interactions
    max_interactions = 0
    @chromarray.each do |num|
      if @chromhash[num].snp_list.get_max_interactions > max_interactions
        max_interactions = @chromhash[num].snp_list.get_max_interactions
      end
    end
    return max_interactions
  end

  # set max and min p val scores
  def set_max_min(grouplisthash)
    @maxscore.each do |scorekey, val|
      @chromhash.each do |chrnum, chrom|
        chrom.snp_list.snps.each do |snp|
          snp.results.each do |groupname, res|
            if groupname !~ /:/
              currgroup = grouplisthash[GroupList.get_default_name].grouphash[groupname]
            else
              namepcs = groupname.split /:/
              currgroup = grouplisthash[namepcs[1]].grouphash[groupname]
            end

            if (scorekey =~ /beta/i and currgroup.has_beta?) or scorekey !~ /beta/i
              if res.values[scorekey] !~ /\d/
                next
              end
              if res.values[scorekey].to_f > @maxscore[scorekey].to_f
                @maxscore[scorekey] = res.values[scorekey]
              end
              if res.values[scorekey].to_f < @minscore[scorekey].to_f
                @minscore[scorekey] = res.values[scorekey]
              end
            end
          end
        end
      end
    end
  end

  # add ld combination to appropriate chromosome
  def add_ld_score(snp_combination)
    @chromhash.each do |key, chrom|
      if chrom.snp_list.snp_hash.has_key?(snp_combination.snp1) and
        chrom.snp_list.snp_hash.has_key?(snp_combination.snp2)
        chrom.snp_list.add_ld_score(snp_combination)
      end
    end
  end

  # adds gene to appropriate chromosome
  def add_gene(gene)
    if @chromhash.has_key?(gene.chr.to_i)
      @chromhash[gene.chr.to_i].add_gene(gene)
    end
  end

  # sorts genes in each chromosome by start position
  def sort_genes
    @chromhash.each_value do |chrom|
      chrom.sort_genes
    end
  end

  # calculate number of rows for the gene information
  def num_gene_rows

    # need width of 1.1 snps to display a gene
    original_gene_space = 1.1
    max_rows = 0

    # remove any genes that are outside the bounds of the
    @chromhash.each_value do |chrom|
      minpos = chrom.snp_list.get_min
      maxpos = chrom.snp_list.get_max
      deletions = Array.new
      chrom.genes.each do |gene|
        if gene.start < minpos or gene.end > maxpos
          deletions << gene
        end
      end

      deletions.each do |d_gene|
        chrom.genes.delete(d_gene)
      end
    end

    # for each chromosome need to calculate how many rows needed
    # and return max across all chromosomes
    @chromhash.each_value do |chrom|
      minpos = chrom.snp_list.get_min
      maxpos = chrom.snp_list.get_max

      num_snps = chrom.snp_list.get_num_included_snps
      snp_size = (maxpos-minpos)/num_snps

      # need approximately space for 3 snps for each label
      # but need to make sure no overlap between adjacent

      # sweep through and take out as many as will fit in a row
      # then make second and additional passes for output
      num_genes = chrom.genes.length
      total_indexes = num_genes
      last_snp = 0

      # mark ones that are used here
      used_gene = Hash.new
      num_rows = 1

      while used_gene.length != num_genes
        total_indexes.times do |i|
          # when a gene is placed mark it as done

          if used_gene[i] == 1
            next
          end
          # locate ending point (in terms of snps for this gene)
          end_snp = (chrom.genes[i].end.to_f - minpos)/snp_size
          start_snp = (chrom.genes[i].start.to_f - minpos) /snp_size

          # line_start and line_end mark location of actual bar on plot
          chrom.genes[i].line_start = start_snp / num_snps
          chrom.genes[i].line_end = end_snp / num_snps
          snp_length = end_snp - start_snp
          chrom.genes[i].alignment = 'middle'
          chrom.genes[i].text_anchor = (end_snp-start_snp)/num_snps + chrom.genes[i].line_start

          gene_space = original_gene_space
          # increase gene space for long gene names
          if chrom.genes[i].name.length >=7
            # gene_space += (chrom.genes[i].name.length-6) * 0.15
            gene_space = 2
          end


          # pad for label around actual location when needed
          if snp_length < gene_space
            start_snp = start_snp - gene_space/2 + snp_length/2
            end_snp = end_snp + gene_space/2 - snp_length/2
          end

          # shift start_snp and end_snp if this would be off edge
          if start_snp < 0
            end_snp = end_snp + start_snp.abs
            start_snp = 0
            chrom.genes[i].alignment = 'start'
            chrom.genes[i].text_anchor = start_snp / num_snps
          end

          if end_snp > num_snps
            start_snp = start_snp - (end_snp - num_snps)
            end_snp = num_snps
            chrom.genes[i].alignment = 'end'

            chrom.genes[i].text_anchor = 1.0
          end

          chrom.genes[i].row = num_rows
          # overlap here so need to skip and check again in next row
          if start_snp < last_snp
            next
          end

          # this gene will fit in row
          used_gene[i] = 1
          # reset last_snp
          last_snp = end_snp

        end # end loop through genes
        num_rows = num_rows + 1
        last_snp = 0
      end # end while loop
      if num_rows > max_rows
        max_rows = num_rows
      end
    end #end chromosome
    return max_rows-1

  end

  # returns true if all snps have results for every SNP in set
  def results_complete?(max_results)
    complete = true
    @chromhash.each_value do |chrom|
      if !chrom.results_complete?(max_results)
        complete = false
        break
      end
    end
    return complete
  end

end

#############################################################################
#
# Class Chromosome -- Holds all information associated with a chromosome.
#                     Includes SNPs and genes.
#
#############################################################################
class Chromosome
	attr_accessor :snp_list, :chrom_num, :genes

	# requires chromosome number
	def initialize(chr_num)
		@chrom_num = chr_num
		@snp_list = SNPList.new
		@genes = Array.new
	end

	def add_gene(gene)
	  @genes << gene
	end

	# sort genes by starting position
	def sort_genes
  	@genes.sort!{|x,y| x.start <=> y.start}
  end

  # checks that all SNPs have results
  def results_complete?(max_results)
    return @snp_list.results_complete?(max_results)
  end

  def get_max_snp_name
    return @snp_list.max_snp_name
  end

end


#############################################################################
#
# Class Gene -- Holds all information associated with a gene such as
#               start position, end position, and name.
#
#############################################################################
class Gene
  attr_accessor :start, :end, :pos_strand, :name, :chr, :alignment, :text_anchor,
    :line_start, :line_end, :row

  def initialize()
    @chr = 0
    @start = @end = 0
    @pos_strand = true
    @name = ""
  end

  def fill_info(st, en, pos, nam, chromosome)
    @start = st
    @end = en
    @pos_strand = pos
    @name = nam
    @chr = chromosome
  end

end


#############################################################################
#
# Class RegressResult -- Holds p value and beta value for result
#
#############################################################################
class RegressResult
	attr_accessor :values

	def initialize
	  @values = Hash.new
	end

end


#############################################################################
#
# Class SNP -- Holds information and results for SNP
#
#############################################################################
class SNP
	attr_accessor :name, :location, :results, :groups

	def initialize(name, location)
		@name = name
		@location = location.to_i
		@results = Hash.new
		@groups = Hash.new
	end

	def group_calculated?(groupname)
		return @results[groupname]
	end

	def add_result(groupname, pval, betaval, maf, sampsizeval, orval, lci, uci, ca, con, cafca, cafcon,
	  studynum, pownum, betauci, betalci, add_columns)

		@results[groupname] =  RegressResult.new
		@results[groupname].values['pvalue'] = pval
	  @results[groupname].values['beta'] = betaval
	  @results[groupname].values['maf'] = maf
	  @results[groupname].values['N'] = sampsizeval
	  @results[groupname].values['or'] = orval
	  @results[groupname].values['uci'] = uci
	  @results[groupname].values['lci'] = lci
	  @results[groupname].values['cases'] = ca
	  @results[groupname].values['controls']=con
	  @results[groupname].values['cafcases'] = cafca
	  @results[groupname].values['cafcontrols'] = cafcon
	  @results[groupname].values['study'] = studynum
    @results[groupname].values['power'] = pownum
    @results[groupname].values['betauci'] = betauci
    @results[groupname].values['betalci'] = betalci
    add_columns.each_pair{|key,value| @results[groupname].values[key]=value}
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
  attr_accessor :snps, :ld_scores, :stats, :snpinfo_tracks, :snp_hash, :included_snps, :included_snps_hash,
    :index_of_included_snps, :maximum, :minimum, :blocks, :max_snp_name

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
    @max_snp_name = 0
  end

  # returns nil when not found
  def get_snp(name)
    return snp_hash[name]
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
    if snp.name.length > @max_snp_name
      @max_snp_name = snp.name.length
    end
  end

  # snps are sorted by basepair location
  def sort_snps
    @snps.sort!{|x,y| x.location <=> y.location}
  end

  def reverse_snps!
    @snps.reverse!
  end

  # returns min position
  def get_min
    return snps[included_snps.first].location
  end

  # returns maximum number of iteractions for any SNP in list
  def get_max_interactions
    max_interaction = 0
#    @ld_scores.values.each do |snp_combo_array|
#      if snp_combo_array.length > max_interaction
#        max_interaction = snp_combo_array.length
#      end
#    end
    
    @ld_scores.each_key do |snpname|
      if @ld_scores[snpname].length > max_interaction
        max_interaction = @ld_scores[snpname].length
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
      @ld_scores[snp_combination.snp1] = Hash.new
      @included_snps_hash[snp_combination.snp1] = true
      @included_snps_hash[snp_combination.snp2] = true
    end
    @ld_scores[snp_combination.snp2]=Hash.new unless @ld_scores.has_key?(snp_combination.snp2)
    #@ld_scores[snp_combination.snp1] << snp_combination
    @ld_scores[snp_combination.snp1][snp_combination.snp2]=snp_combination
    @ld_scores[snp_combination.snp2][snp_combination.snp1]=snp_combination
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


  # sets all snps to be included for this plot
  def make_all_included
    @snps.each_with_index do |snp, index|
      @included_snps_hash[snp.name] = snp
      @included_snps << index
      @index_of_included_snps[snp.name] = @included_snps.length-1
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

  # checks all SNPs to see if every snp is complete
  def results_complete?(max_results)
    complete = true
    @included_snps.each do |snp|
      if max_results != @snps[snp].results.length
        complete = false
        break
      end
    end
    return complete
  end


end


#############################################################################
#
# Class Group -- Group information and location in Synthesis-View input file
#
#############################################################################
class Group
  attr_accessor :name, :pcol, :betacol, :colorstr, :mafcafcol, :values, :Ncol, :subgroups,
    :orcol, :ucicol, :lcicol, :casescol, :controlscol, :cafcasescol, :cafcontrolscol, :studycol,
    :fullname, :highlighted, :powercol, :betaucicol, :betalcicol, :additional_cols

  #def initialize(n, pc, bc, freqcol, col, sampsizecol, oddsratio, uci, lci, cacol, concol,
  #  cafcacol, cafconcol, studycol, powcol, beta_ucicol, beta_lcicol, hilite=false)
  def initialize(n, hilite, col=nil)
    @name = n
    @pcol = -1
    @betacol = -1
    @colorstr = col
    @mafcafcol = -1
    @values = Hash.new
    @values['num'] = 0
    @values['pheno_avg'] = 0
    @Ncol = -1
    @subgroups = Hash.new
    @orcol = -1
    @ucicol = -1
    @lcicol = -1
    @casescol = -1
    @controlscol = -1
    @cafcasescol = -1
    @cafcontrolscol = -1
    @powercol = -1
    @studycol = -1
    @highlighted = hilite
    @betaucicol = -1
    @betalcicol = -1
    @additional_cols = Hash.new
  end

  def has_beta?
    retval = true
    if @betacol < 0
      retval = false
    end
    return retval
  end

  def has_median?
    if @values[:median].to_f > 0
      return true
    else
      return false
    end
  end

  def has_overall_sample_size?
    if @values['num'].to_i > 0
      return true
    else
      return false
    end
  end

  def has_pheno_avg?
    if @values['pheno_avg'].to_i > 0
      return true
    else
      return false
    end
  end

  # adds subgroup to this group such as AA or EA
  def add_subgroup(g)
    @subgroups[g.name] = g
    # needs same color as the main q
    @subgroups[g.name].colorstr = @colorstr
  end

  # returns all the subgroup keys for this group
  def get_subgroup_names
    names = Array.new
    @subgroup.keys.each do |key|
      keypcs = key.split/:/
      names << keypcs[1]
    end
    return names
  end

end


#############################################################################
#
# Class GroupList -- Lists groups and locations for each group in file
#
#############################################################################
class GroupList
  @@color_array = ['rgb(55,126,184)', 'rgb(228,26,28)', 'rgb(255,127,0)', 'rgb(152,78,163)', 'rgb(77,175,74)', #'rgb(255,255,51)',
    'rgb(166,86,40)', 'rgb(247,129,191)', 'yellow']
  @@color_index = 0
  @@defaultname = ''
  attr_accessor :groups, :grouphash, :mafcoltitle

  def initialize()
    @groups = Array.new
    @grouphash = Hash.new
    @color_index = -1
    @mafcoltitle = 'MAF'
  end

  def self.get_default_name
    return @@defaultname
  end

  def self.get_next_color
    colorstr = @@color_array[@@color_index]
    @@color_index += 1
    if @@color_index > @@color_array.length - 1
      @@color_index = 0
    end
    return colorstr
  end
  
  def self.set_grayscale
    @@color_array = ['rgb(82,82,82)', 'rgb(204,204,204)', 'rgb(150,150,150)', 'rgb(247,247,247)']
  end

  def add_group(g)
    @groups << g
    @grouphash[g.name] = @groups.last
  end

  # returns group by name
  def get_group_by_name(groupname)
    return @grouphash[groupname]
  end

  # returns true when the groups have values for overall sample size
  def plot_summary_size?
    make_plot = false
    @groups.each do |g|
      if g.has_overall_sample_size?
        make_plot = true
        break
      end
    end
    return make_plot
  end

  # return true when the groups have values for the phenotype averages
  def plot_pheno_avg?
    make_plot = false
    @groups.each do |g|
    if g.has_pheno_avg?
        make_plot = true
        break
      end
    end
    return make_plot
  end

  def plot_box_plot?
    make_plot = false
    @groups.each do |g|
    if g.has_median?
        make_plot = true
        break
      end
    end
    return make_plot
  end

  def plot_oddsratio?
    plot = false
    @groups.each do |g|
      if g.orcol > -1
        plot = true
        break
      end
    end
    return plot
  end

  def plot_betas?
    plot = false
    @groups.each do |g|
      if g.betacol > -1
        plot = true
        break
      end
    end
    return plot
  end

  def plot_pvals?
    plot = false
    @groups.each do |g|
      if g.pcol > -1
        plot = true
        break
      end
    end
    return plot
  end

  def plot_sample_sizes?
    plot = false
    @groups.each do |g|
      if g.Ncol > -1
        plot = true
        break
      end
    end
    return plot
  end

  def plot_maf?
    plot = false
    @groups.each do |g|
      if g.mafcafcol > -1
        plot = true
        break
      end
    end
    return plot
  end

  def plot_study?
    plot = false
    @groups.each do |g|
      if g.studycol > -1
        plot = true
        break
      end
    end
    return plot
  end

end

#############################################################################
#
# Class DBSNPhtml -- Create html file with image map for linking to dbSNP information for each SNP
#
#############################################################################
class DBSNPhtml

# create html file with image map for linking to dbSNP information for each SNP
def self.write_html(imgfile, xsize, ysize, linkPositions, htmlfile, rotate)

  File.open(htmlfile, 'w') do |fwrite|
    fwrite.puts "<html><head><title>Synthesis-View</title></head>\n<body>"
    fwrite.puts "<img src=\"#{imgfile}\" width=\"#{xsize.to_i}\" ysize=\"#{ysize.to_i}\" alt=\"Synthesis-View\" usemap=\"#synthmap\" />"

    fwrite.puts "<map name=\"synthmap\">"
    # loop through and write out position for each SNP with link to dbSNP
    # http://www.ncbi.nlm.nih.gov/sites/entrez?db=snp&cmd=search&term=rs802019
    linkPositions.each do |snp|
        fwrite.puts "<area shape=\"rect\" coords=\"#{snp.point1.x.to_i},#{snp.point1.y.to_i},#{snp.point2.x.to_i},#{snp.point2.y.to_i}\"" +
         " href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=snp&cmd=search&term=#{snp.linkid}\" alt=\"#{snp.linkid}\" target=\"_blank\" />"
    end
    fwrite.puts "</map>"
    fwrite.puts "</body></html>"
  end

end

end

#############################################################################
#
# Class Point -- Records x and y locations
#
#############################################################################
class Point
  attr_accessor :x, :y
  def initialize(xtemp,ytemp)
    @x=xtemp
    @y=ytemp
  end
end

#############################################################################
#
# Class ImagePositions -- Holds information for link on image map
#
#############################################################################
# records image positions for use in creating links to dbSNP
class ImagePositions
  attr_accessor :point1, :point2, :linkid

  def initialize(p1, p2, id)
    @point1 = p1
    @point2 = p2
    @linkid = id
  end

end


#############################################################################
#
# Class Shape -- Records basic information for plotting
#
#############################################################################
class Shape
  attr_accessor :color, :ymidpoint, :xmidpoint, :stroke, :interval
  @@grayscale = false
  
  def self.set_grayscale
    @@grayscale=true
    @interval = nil
  end
  
  def self.get_grayscale
    @@grayscale
  end
  
  def draw(pen)
  end

  def adjust_x(xval)
  end

end


#############################################################################
#
# Class Circle -- Records basic information for plotting a circle
#
#############################################################################
class Circle < Shape
  attr_accessor :x,:y,:radius

  def initialize(xpt, ypt, rad, col)
    @x = xpt
    @y = ypt
    @radius = rad
    @color = col
    @ymidpoint = @y
    @xmidpoint = @x
  end

  def draw(pen)
    @stroke = 'none'
    if @@grayscale
      @stroke = 'black'
    end
    
    pen.styles(:fill=>@color, :stroke=>@stroke)
    pen.circle( @radius, @x, @y)
  end

  def adjust_x(xval)
    @x += xval
    @xmidpoint = @x
  end

end

#############################################################################
#
# Class Line -- Records basic information for plotting a line
#
#############################################################################
class Line < Shape
  attr_accessor :xpoints
  
  def initialize(xpts, ypts, col)
    @xpoints = xpts
    @ypoints = ypts
    @color = col
    @ymidpoint = (ypts[0] + ypts[1])/2
    @xmidpoint = (@xpoints[0]+@xpoints[1])/2
  end  
  
  def draw(pen)
    @stroke = @color
    if @@grayscale 
      @stroke = 'black'
    end
    pen.styles(:stroke=>@stroke)
    
    pen.line(@xpoints[0], @ypoints[0], @xpoints[1], @ypoints[1])
  end
 
  def adjust_x(xval)
    @xpoints.each_with_index do |x,i|
      @xpoints[i] += xval
    end
    @xmidpoint = (@xpoints[0]+@xpoint[1])/2    
  end
  
end


#############################################################################
#
# Class Triangle -- Records basic information for plotting a triangle
#
#############################################################################
class Diamond < Shape
  attr_accessor :xpoints, :ypoints

  def initialize(xpts, ypts, col)
    @xpoints = xpts
    @ypoints = ypts
    @color = col
    @ymidpoint = (ypts[1] + ypts[2])/2
    @xmidpoint = @xpoints[0]
  end

  def draw(pen)
    @stroke = @color
    if @@grayscale
      @stroke = 'black'
    end
    pen.styles(:fill=>@color, :stroke=>@stroke)
    pen.polygon(@xpoints, @ypoints)
  end

  def draw_open(pen)
    pen.styles(:fill=>'none', :stroke=>@color)
    pen.polygon(@xpoints, @ypoints)
  end

  def adjust_x(xval)
    @xpoints.each_with_index do |x,i|
      @xpoints[i] += xval
    end
    @xmidpoint = @xpoints[0]
  end

end

#############################################################################
#
# Class Triangle -- Records basic information for plotting a triangle
#
#############################################################################
class Triangle < Shape
  attr_accessor :xpoints,:ypoints

  def initialize(xpts, ypts, col)
    @xpoints = xpts
    @ypoints = ypts
    @color = col
    @ymidpoint = (ypts[1] + ypts[2])/2
    @xmidpoint = @xpoints[2]
  end

  def draw(pen)
    @stroke = @color
    if @@grayscale
      @stroke = 'black'
    end
    pen.styles(:fill=>@color, :stroke=>@stroke)
    pen.polygon(@xpoints, @ypoints)
  end

  def adjust_x(xval)
    @xpoints.each_with_index do |x,i|
      @xpoints[i] += xval
    end
    @xmidpoint = @xpoints[2]
  end

end


#############################################################################
#
# Class PlotWriter -- Contains functions for drawing plots.  Most distances
#                     are adjusted based on @box_size.  @box_size
#                     defines the size in x and y coordinates for the
#                     elements in the D-prime and R-squared plots.
#
#############################################################################
class PlotWriter
  attr_accessor :box_size, :canvas, :color_array, :dashed_array, :interval_array, :feature_opacity,
    :font_size_multiple, :first_grid, :map_pos, :pixels_per_coordinate, :maxy, :maxx

  def initialize(box_size)
    @box_size = box_size
    @color_array = ['black','red', 'blue','gray','orange','green','purple']
    @dashed_array = [3, 10, 17, 23, 29, 34, 39];
    @interval_array = [3, 3, 6, 6, 9, 9, 9];
    @feature_opacity = 1.0
    @font_size_multiple = 1.0
    @first_grid = true
    @map_pos = Array.new
    @pixels_per_coordinate = 0
    @maxy =0
    @maxx=0
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


  # adjusts points to offset points
  def adjust_horizontal(points, y_interval)
    # determine minimum distance for offset (1/25 of overall height)
    jitter_interval = y_interval / 40
    offset_hash = Hash.new
    offset_group = Array.new
    curr_offset_group = -1
    points.length.times do |i|
      j = i+1
      while j < points.length
        if (points[i].ymidpoint - points[j].ymidpoint).abs < jitter_interval
          if !offset_hash[i]
            curr_offset_group += 1
            offset_group << Array.new
            offset_group[curr_offset_group] << i
          end
          if !offset_hash[j]
            offset_group[curr_offset_group] << j
          end
          offset_hash[i]=true
          offset_hash[j]=true
        end
        j += 1
      end
    end
    # adjust all points in the offset groups
    # first adjusted left, then to the right
    offset_group.each do |offset_array|
      adjust_x = -@box_size.to_f/5 * (offset_array.length.to_f-1)/2
      offset_array.each do |offset_point|
        points[offset_point].adjust_x(adjust_x)
        adjust_x+=@box_size.to_f/5
      end
    end
  end

  # returns appropriate color based on score
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

  # draws gene information from chromosome genes
  def add_gene_rows(chrom, x_start, y_start, x_end, y_end, nrows)

    if chrom.snp_list.get_num_included_snps < 2 or chrom.genes.length < 1
      return
    end

    start_plot_x = @box_size
    end_plot_x = x_end - x_start -@box_size * 0.7
    start_plot_y = 0
    end_plot_y = y_end - y_start - @box_size * 0.5

    # needs a box all the way around to denote area
    @canvas.g.translate(x_start, y_start) do |plot|
      plot.styles(:stroke_width=>1, :stroke=>'black')
      plot.line(start_plot_x, start_plot_y, end_plot_x, start_plot_y)
      plot.line(start_plot_x, end_plot_y, end_plot_x, end_plot_y)
      plot.line(start_plot_x, start_plot_y, start_plot_x, end_plot_y)
      plot.line(end_plot_x, start_plot_y, end_plot_x, end_plot_y)
    end

    # make row be the size of one box
    row_size = @box_size
    x_total_width = end_plot_x - start_plot_x
    # for each gene place in appropriate row with label above
    chrom.genes.each do |gene|
      y_line = gene.row * row_size
      x_line_start = start_plot_x + gene.line_start * x_total_width
      x_line_end = start_plot_x + gene.line_end * x_total_width

      if x_line_end - x_line_start < 3
        x_line_start = x_line_start-1
        x_line_end = x_line_end+1
      end

      # draw line at appropriate row and position
      @canvas.g.translate(x_start, y_start) do |plot|
        plot.styles(:stroke_width=>2, :stroke=>'navy')
        plot.line(x_line_start, y_line, x_line_end, y_line)
      end

      # add label above
      text_anchor_x = start_plot_x + gene.text_anchor * x_total_width
      font_size = standard_font_size * 0.68
      @canvas.g.translate(x_start, y_start).text(text_anchor_x, y_line-@box_size/4) do |text|
        text.tspan(gene.name).styles(:font_size=>font_size, :text_anchor=>gene.alignment)
      end

    end
  end

  # returns new triangle shape
  def create_triangle(box_y_start, value_x, up_triangle, box_mult, colorstr)
    x_points = Array.new
    x_points << value_x - @box_size*box_mult * Math.sqrt(2)*0.2
    x_points << value_x + @box_size*box_mult * Math.sqrt(2)*0.2
    # third point is midway between other 2
    x_points << (x_points[0]+x_points[1])/2

     y_points = Array.new
	   y_points << box_y_start
	   y_points << box_y_start
	   if up_triangle
       y_points << box_y_start - @box_size*box_mult * Math.sqrt(2) *0.35
     else
       y_points << box_y_start + @box_size*box_mult * Math.sqrt(2) *0.35
	   end

	   return Triangle.new(x_points, y_points, colorstr)
  end

  # returns a new diamond shape
  def create_diamond(box_y_start, value_x, box_mult, colorstr)

    # establish x coordinates
    x_points = Array.new
    left_x = value_x - @box_size*box_mult * Math.sqrt(2)*0.2
    right_x = value_x + @box_size*box_mult * Math.sqrt(2)*0.2
    x_points << (left_x + right_x)/2

    x_points << right_x
    x_points << x_points[0]
    x_points << left_x

    # establish y coordinates
    y_points = Array.new
    y_points << box_y_start - @box_size*box_mult * Math.sqrt(2) *0.35
	  y_points << box_y_start
    y_points << box_y_start + @box_size*box_mult * Math.sqrt(2) *0.35
    y_points << box_y_start

    return Diamond.new(x_points, y_points, colorstr)
  end

  # adds plot with p values where betas determine direction of the plot
#  def draw_pvalue_plot(jitter, no_triangles, grouplist, snp_list, x_start, y_start, stat_max,
#    stat_min, original_min, first_plot, prefix_title='', rotate=false)
  def draw_pvalue_plot(params)
    jitter = params[:jitter]
    no_triangles = params[:no_triangles]
    grouplist = params[:grouplist]
    snp_list = params[:snp_list]
    x_start = params[:x_start]
    y_start = params[:y_start]
    stat_max = params[:stat_max]
    stat_min = params[:stat_min]
    original_min = params[:original_min]
    first_plot = params[:first_plot]
    prefix_title = params[:prefix_title] || ''
    rotate = params[:rotate] || false
    clean_axis = params[:clean_axis] || false
    
    x = @box_size
    xmin = x

    ymin = 0
    ymax = @box_size*9
    y_interval = ymax-ymin

    check_max = false

    stat_max = get_neg_log(stat_min)
    top_max = get_neg_log(stat_min)

    #############
    #######  CHECK ALL Draw_pvalue calls ###################

    # check to see if there is a cut-off imposed for the pvalue plot
    if stat_min.to_f > original_min.to_f
      check_max = true
      box_mult=1.8
      # calculate size of triangle (or circle) and adjust the p value plot max size
      y_vertical = @box_size*box_mult * Math.sqrt(2) *0.35
      y_fraction = 1- y_vertical / y_interval

      # adjust stat_max to be the value needed for a triangle to just touch top
      stat_max = top_max / y_fraction
    # if no cut-off imposed check for the use of a cleaner axis
    elsif clean_axis
      increment, stat_min, stat_max = calculate_increments(0, stat_max.to_f)
    end

    stat_min = 0

    # set stat_min to be zero
    if !check_max and stat_max - stat_max.floor > 0
      stat_max = stat_max.floor + 1
    end

    # store for later use based when adding features
    @functional_max = stat_max

    stat_interval = stat_max - stat_min
    value_x = xmin * 0.7 + (@box_size * Math.sqrt(2)*0.25)/2

    # draw box (change to triangle) with correct color
    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      draw_separation(x_start, y_start, 0, ymax, value_x)

      # store points here for outputting
      points = Array.new

      grouplist.groups.each do |group|

        if snp.results.has_key?(group.name)
          over_max = false
          box_mult=1
          result = snp.results[group.name]
          pvalue = get_neg_log(result.values['pvalue'].to_f)
          if pvalue > top_max
            over_max = true
            pvalue = top_max
            box_mult=1.8
          end

          box_y_start = ((stat_max-pvalue) / stat_interval) * y_interval
	        if (!no_triangles and group.has_beta?)
	          if over_max and result.values['beta'].to_f <= 0
  	            box_y_start = 0
	          end

	          if rotate
              box_y_start = y_interval - box_y_start
            end

            if (result.values['beta'].to_f > 0 and !rotate) or (result.values['beta'].to_f <= 0 and rotate)
              up_triangle = true
            else
              up_triangle = false
            end

            points << create_triangle(box_y_start, value_x, up_triangle, box_mult, group.colorstr)

  	      elsif !group.highlighted

  	        if over_max
  	          box_y_start = @box_size * box_mult * Math.sqrt(2)*0.125/2
  	        end

  	        if rotate
              box_y_start = y_interval - box_y_start
            end
    	      points << Circle.new(value_x, box_y_start, @box_size * box_mult * Math.sqrt(2)*0.125, group.colorstr)
    	    else
    	      #draw diamond for highlighted group
    	      if over_max
  	            box_y_start = 0
	          end
	          if rotate
              box_y_start = y_interval - box_y_start
            end
            points << create_diamond(box_y_start, value_x, box_mult, group.colorstr)
	        end
        end

      end

      # sort by y-midpoint
      points.sort! {|x,y| y.ymidpoint <=> x.ymidpoint}

      # when jitter selected check for overlap and offset horizontally
      # those points that overlap significantly
      if jitter
        adjust_horizontal(points, y_interval)
      end

	    # iterate through points and draw all of them
	    points.each do |point|
  	    @canvas.g.translate(x_start,y_start) do |check|
         check.styles(:stroke_width=>1, :fill_opacity=>0.9)
  	      point.draw(check)
  	    end
  	  end
      value_x = label_step(value_x)
    end

      # write labels on left side of stat box
    if first_plot
      if rotate
        x_label_start = value_x + x_start+ @box_size * 0.4
      else
        x_label_start = x_start
      end

        label_precision = 1

      write_plot_labels(x_label_start, y_start, ymax, stat_max, stat_min, y_interval, 0,
        prefix_title + '-log10(p value)', false, label_precision, rotate)
    end
  end

  # for drawing plots for beta values and mean allele frequencies
  # draws a dot in place of the triangle produced for the p value plot    
  def draw_basic_plot(params)  
    jitter = params[:jitter]
    grouplist = params[:grouplist]
    snp_list = params[:snp_list]
    x_start = params[:x_start]
    y_start = params[:y_start]
    stat_max = params[:stat_max]
    stat_min = params[:stat_min]
    data_key = params[:data_key]
    plot_labels = params[:plot_labels]
    title = params[:title]
    precision = params[:precision]
    rotate = params[:rotate] || false
    size_mult = params[:size_mult] || 1
    prefix_title = params[:prefix_title] || ""
    lci_key = params[:lci_key] || nil
    uci_key = params[:uci_key] || nil
    
    x = @box_size
    xmin = x
    ymin = 0
    ymax = @box_size*9 * size_mult
    y_interval = ymax-ymin
    # set stat_min to be zero
    stat_interval = stat_max - stat_min

    value_x = xmin * 0.7 + (@box_size * Math.sqrt(2)*0.25)/2
    
    fill_opacity = 0.9
    if Shape.get_grayscale
      fill_opacity = 1.0
    end
    
    # draw box (change to triangle) with correct color
    snp_list.included_snps.each do |snp_index|
      draw_separator = true
      snp = snp_list.snps[snp_index]
      
      #draw_separation(x_start, y_start, 0, ymax, value_x) unless lci_key

      # store points here for output
      points = Array.new
      lines = Array.new

      grouplist.groups.each do |group|
        if snp.results.has_key?(group.name)
          result = snp.results[group.name]
          if !result.values[data_key] or result.values[data_key] !~ /\d/
            next
          end
          box_y_start = ((stat_max-result.values[data_key].to_f) / stat_interval) * y_interval

          if rotate
            box_y_start = y_interval - box_y_start
          end

          if !group.highlighted
            points << Circle.new(value_x, box_y_start, @box_size * Math.sqrt(2)*0.125, group.colorstr)
          else
            # create triangle here
            points << create_diamond(box_y_start, value_x, 1.0, group.colorstr)
          end
          
          # create a line if there is a confidence interval
          if result.values[uci_key]
            draw_separator = false
            uci_y = ((stat_max-result.values[uci_key].to_f) / stat_interval) * y_interval
            lci_y = ((stat_max-result.values[lci_key].to_f) / stat_interval) * y_interval
            if rotate
              uci_y = y_interval - uci_y
              lci_y = y_interval - lci_y
            end           

            points.last.interval = Line.new([value_x,value_x], [lci_y,uci_y], group.colorstr)
          end
          
        end

      end
      draw_separation(x_start, y_start, 0, ymax, value_x)# if draw_separator
      # sort by y-midpoint
      points.sort! {|x,y| y.ymidpoint <=> x.ymidpoint}
      #lines.sort! {|x,y| y.ymidpoint <=> x.ymidpoint} unless lines.empty?

      # those points that overlap significantly
      if jitter
        adjust_horizontal(points, y_interval)
      end

	    # iterate through points and draw all of them
	    points.each do |point|
        if point.interval            
          point.interval.xpoints[0] = point.interval.xpoints[1] = point.xmidpoint
          @canvas.g.translate(x_start,y_start) do |check|
            check.styles(:stroke_width=>2)
            point.interval.draw(check)
          end
        end 
  	  end

	    points.each do |point|
  	    @canvas.g.translate(x_start,y_start) do |check|
    	    check.styles(:stroke_width=>1, :fill_opacity=>fill_opacity)
  	      point.draw(check)
  	    end 
      end
    
      value_x = label_step(value_x)
    end

    if plot_labels
      if rotate
        x_label_start = value_x + x_start+ @box_size * 0.4 
      else
        x_label_start = x_start 
      end

      if size_mult >= 1
        write_plot_labels(x_label_start, y_start, ymax, stat_max, stat_min, y_interval, 0, prefix_title+title, false, precision, rotate)
      else
        write_plot_labels(x_label_start, y_start, ymax, stat_max, stat_min, y_interval, 0, prefix_title+title,
          false, precision, rotate, 5)
      end
    end
  end

  # draw plot where one set of numbers are presented as filled circles and
  # the second as closed circles
  def draw_circles_plot(grouplist, snp_list, x_start, y_start, stat_max, stat_min, closed_key, open_key,
    plot_labels, title, precision, rotate=false)
    x = @box_size
    xmin = x
    ymin = 0
    ymax = @box_size*9
    y_interval = ymax-ymin
    # set stat_min to be zero
    stat_interval = stat_max - stat_min
    value_x = xmin * 0.7 + (@box_size * Math.sqrt(2)*0.25)/2

    # draw box (change to triangle) with correct color
    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      draw_separation(x_start, y_start, 0, ymax, value_x)
      grouplist.groups.each do |group|
        if snp.results.has_key?(group.name)
          result = snp.results[group.name]
          if result.values[closed_key] and result.values[closed_key] =~ /\d/

            box_y_start = ((stat_max-result.values[closed_key].to_f) / stat_interval) * y_interval

            if rotate
              box_y_start = y_interval - box_y_start
            end

            @canvas.g.translate(x_start,y_start) do |check|
  	          check.styles(:fill=>group.colorstr, :stroke=>group.colorstr, :stroke_width=>1, :fill_opacity=>0.9)

              if !group.highlighted
 	            # need to draw a circle with point above or below line
               check.circle(@box_size * Math.sqrt(2)*0.125, value_x, box_y_start)
             else
               diamond = create_diamond(box_y_start, value_x, 1.0, group.colorstr)
               diamond.draw(check)
             end

	          end
	        end

	        # plot the open circle if present
          if !result.values[open_key] or result.values[open_key] !~ /\d/
            next
          end

          box_y_start = ((stat_max-result.values[open_key].to_f) / stat_interval) * y_interval

          if rotate
            box_y_start = y_interval - box_y_start
          end

          @canvas.g.translate(x_start,y_start) do |check|
  	        check.styles(:fill=>'none', :stroke=>group.colorstr, :stroke_width=>1, :fill_opacity=>0.0)
            if !group.highlighted
 	            # need to draw a circle with point above or below line
               check.circle(@box_size * Math.sqrt(2)*0.125, value_x, box_y_start)
             else
               diamond = create_diamond(box_y_start, value_x, 1.0, group.colorstr)
               diamond.draw_open(check)
             end
	        end

        end
      end
      value_x = label_step(value_x)
    end

    if plot_labels
      if rotate
        x_label_start = value_x + x_start+ @box_size * 0.4
      else
        x_label_start = x_start
      end
      write_plot_labels(x_label_start, y_start, ymax, stat_max, stat_min, y_interval, 0, title, false, precision, rotate)
    end
  end

  # adds forest plot legend
  def add_forest_legend(x_legend_rows, y_row_start, grouplist, xstart, ystart)
    font_size = standard_font_size * 1.3
    mult = 1

    text_rotate = 90

    max_abbrev_length=0
    grouplist.groups.each do |group|
      if group.name.length > max_abbrev_length
         max_abbrev_length = group.name.length
      end
    end

    i=0
    grouplist.groups.reverse_each do |group|
      # draw box
      # want box to have bottom match the bottom of font
      @canvas.g.translate(xstart,ystart) do |check|
        check.styles(:fill=>group.colorstr, :stroke=>group.colorstr, :stroke_width=>2, :fill_opacity=>0.9)
        # adjust start of box so that mid-point marks actual odds ratio
     	  check.rect(box_size*mult,  box_size*mult,x_legend_rows[i], y_row_start)
     	end

 	    text_y_start = y_row_start + box_size * mult + box_size/2

      # draw abbreviation
	    @canvas.g.translate(xstart, ystart).text(x_legend_rows[i], text_y_start).rotate(text_rotate) do |text|
	      text.tspan(group.name).styles(:font_size=>font_size, :text_anchor=>'start')
	    end

	    text_y_start += (2 * box_size + max_abbrev_length * box_size/4) *  @font_size_multiple

      # draw full name
      if group.fullname
        @canvas.g.translate(xstart, ystart).text(x_legend_rows[i], text_y_start).rotate(text_rotate) do |text|
	        text.tspan("#{group.fullname}").styles(:font_size=>font_size, :text_anchor=>'start')
	      end
	    end
      i+=1
    end
  end


  # adds legend for the plots where some values are closed circles and some are open circles
  def add_circle_plot_legend(closed_label, open_label, rotated, x_start, y_start, plot_below=false)

     circle_center_y = @box_size
    circle_center_x = @box_size
    font_size = standard_font_size

    # Place legend after or before end of the x_start
    if rotated
      x_start -= @box_size/1.1
      text_rotate = 90
      anchor = 'start'
      text_adjust = circle_center_x-@box_size * Math.sqrt(2)*0.125
      y_adjust = @box_size * Math.sqrt(2)*0.125+@box_size/2
    elsif plot_below
      if @font_size_multiple == 1.0
        x_start -= @box_size * 7
      else
        x_start -= @box_size * 10
      end
      y_start += @box_size * 9
      text_adjust = circle_center_x+@box_size * Math.sqrt(2)*0.125 + @box_size/4
      y_adjust = @box_size * Math.sqrt(2)*0.125
      text_rotate = 0
      anchor = 'start'
    else
      x_start += -@box_size
      text_rotate = 0
      anchor = 'start'
      text_adjust = circle_center_x+@box_size * Math.sqrt(2)*0.125 + @box_size/4
      y_adjust = @box_size * Math.sqrt(2)*0.125
    end

    # need circles followed by text
    @canvas.g.translate(x_start,y_start) do |check|
      check.styles(:fill=>'black', :stroke=>'black', :stroke_width=>1, :fill_opacity=>0.9)

	    # need to draw a circle with point above or below line
	    check.circle(@box_size * Math.sqrt(2)*0.125, circle_center_x, circle_center_y)
	  end

	  @canvas.g.translate(x_start, y_start).text(text_adjust, circle_center_y + y_adjust).rotate(text_rotate) do |text|
	    text.tspan(closed_label).styles(:font_size=>font_size, :text_anchor=>anchor)
	  end

	  if rotated
	    circle_center_x -= @box_size/1.5
	    text_adjust = circle_center_x-@box_size * Math.sqrt(2)*0.125
	  elsif plot_below
	    if @font_size_multiple == 1.0
  	    circle_center_x += @box_size * 3
  	  else
  	    circle_center_x += @box_size * 4
  	  end
	    text_adjust = circle_center_x+@box_size * Math.sqrt(2)*0.125 + @box_size/4
	  else
      circle_center_y += @box_size
	  end

    # need circles followed by text
    @canvas.g.translate(x_start,y_start) do |check|
      check.styles(:fill=>'none', :stroke=>'black', :stroke_width=>1, :fill_opacity=>0.0)
	    # need to draw a circle with point above or below line
	    check.circle(@box_size * Math.sqrt(2)*0.125, circle_center_x, circle_center_y)
	  end

	  @canvas.g.translate(x_start, y_start).text(text_adjust, circle_center_y + y_adjust).rotate(text_rotate) do |text|
	    text.tspan(open_label).styles(:font_size=>font_size, :text_anchor=>anchor)
	  end

  end

  # for drawing plots for beta values and mean allele frequencies
  # draws a dot in place of the triangle produced for the p value plots
  def draw_basic_or_plot(jitter, grouplist, snp_list, x_start, y_start, stat_max, stat_min, data_key, plot_labels, title,
    precision, large_sig=false, rotate=false, size_mult=1)
    x = @box_size
    xmin = x
    ymin = 0
    ymax = @box_size*9 * size_mult
    y_interval = ymax-ymin
    # set stat_min to be zero
    stat_interval = stat_max - stat_min

    value_x = xmin * 0.7 + (@box_size * Math.sqrt(2)*0.25)/2

    # draw box (change to triangle) with correct color
    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      draw_separation(x_start, y_start, 0, ymax, value_x)

      # store points here for output
      points = Array.new

      grouplist.groups.each do |group|
        if snp.results.has_key?(group.name)
          result = snp.results[group.name]
          if !result.values[data_key] or result.values[data_key] !~ /\d/
            next
          end
          box_y_start = ((stat_max-result.values[data_key].to_f) / stat_interval) * y_interval

          if rotate
            box_y_start = y_interval - box_y_start
          end

          diameter = @box_size * Math.sqrt(2)*0.125
          if large_sig and ((result.values[data_key].to_f > 1.0 and result.values['lci'] == nil) or
            result.values['lci'].to_f > 1.0) or ((result.values[data_key].to_f < 1.0 and result.values['uci'] == nil) or
            result.values['uci'].to_f < 1.0)
            diameter *= 2
          end

          points << Circle.new(value_x, box_y_start, diameter, group.colorstr)
        end
      end

      # sort by y-midpoint
      points.sort! {|x,y| y.ymidpoint <=> x.ymidpoint}

      # those points that overlap significantly
      if jitter
        adjust_horizontal(points, y_interval)
      end

	    # iterate through points and draw all of them
	    points.each do |point|
  	    @canvas.g.translate(x_start,y_start) do |check|
    	    check.styles(:stroke_width=>1, :fill_opacity=>0.9)
  	      point.draw(check)
  	    end
  	  end

      value_x = label_step(value_x)
    end

    if plot_labels
      if rotate
        x_label_start = value_x + x_start+ @box_size * 0.4
      else
        x_label_start = x_start
      end
      if size_mult >= 1
        write_plot_labels(x_label_start, y_start, ymax, stat_max, stat_min, y_interval, 0, title, false, precision, rotate)
      else
        write_plot_labels(x_label_start, y_start, ymax, stat_max, stat_min, y_interval, 0, title,
          false, precision, rotate, 5)
      end
    end
  end


  # draw forest plot for odds ratio results
  # will always be rotated
  def draw_or_plot(group, snp_list, x_start, y_start, ormin, ormax, plot_labels, sig_larger)
    x = @box_size
    xmin = x

    ymin = 0
    ymax = @box_size*9
    y_interval = ymax-ymin

    stat_interval = ormax.to_f - ormin.to_f
    value_x = xmin * 0.7 + (@box_size * Math.sqrt(2)*0.25)/2

    temp_x = label_step(value_x)
    block_x_size = (value_x-temp_x).abs


    data_key = 'or'

    # determine location of odds ratio 1.0
    y_line = y_interval - (ormax-1.0) / stat_interval * y_interval

    # draw line across plot area at 1.0
    @canvas.g.translate(x_start,y_start) do |plot|
      plot.line(value_x, y_line, block_x_size * (snp_list.included_snps.length+1),y_line).styles(:stroke_width=>1, :stroke_opacity=>'lightgray', :opacity=>0.75)
    end

    value_x = value_x - (value_x - temp_x).abs

    # draw circle with line for the upper and lower bounds of the confidence interval
    # for every SNP that has an odds ratio for the indicated group at that SNP
    snp_list.included_snps.each do |snp_index|
      value_x = label_step(value_x)
      snp = snp_list.snps[snp_index]
      draw_separation(x_start, y_start, 0, ymax, value_x)
      if snp.results.has_key?(group.name)
          result = snp.results[group.name]
          if !result.values[data_key] or result.values[data_key] !~ /\d/
            next
          end

          # check for larger box if lower confidence interval is greater than 1.0
          this_box = @box_size
          if sig_larger and (result.values['lci'].to_f > 1.0 or result.values['uci'].to_f < 1.0)
            this_box *= 2
          end

          box_y_start = y_interval - ((ormax.to_f-result.values[data_key].to_f) / stat_interval) * y_interval
          box_uci =  y_interval - ((ormax.to_f-result.values['uci'].to_f) / stat_interval) * y_interval
          box_lci =  y_interval - ((ormax.to_f-result.values['lci'].to_f) / stat_interval) * y_interval
          @canvas.g.translate(x_start,y_start) do |check|
  	        check.styles(:fill=>group.colorstr, :stroke=>group.colorstr, :stroke_width=>2, :fill_opacity=>0.9)

            or_box_x = value_x - Math.sqrt(2) *0.25/2 * this_box

            # adjust start of box so that mid-point marks actual odds ratio
 	          check.rect(this_box * Math.sqrt(2) *0.25, this_box * Math.sqrt(2) *0.25, or_box_x, box_y_start-this_box * Math.sqrt(2) *0.25/2)

	          # draw line for confidence interval
	            check.line(value_x, box_lci, value_x, box_uci)
	        end
        end
    end

    value_x = label_step(value_x)

    if plot_labels
      x_label_start = value_x + x_start + @box_size * 0.4
      write_plot_labels(x_label_start, y_start, ymax, ormax, ormin, y_interval, 0, group.name + " Odds Ratio", false, 2, true)
    end

  end

  # draws faint line to mark space between SNPs in the plot
  def draw_separation(x_start, y_start, y_min, y_max, x)
    @canvas.g.translate(x_start, y_start) do |lines|
      lines.styles(:stroke=>'#CCCCCC', :stroke_width=>1, :opacity=>0.05)
      lines.line(x, y_min, x, y_max)
    end
  end

  # adds stat box labels
  def write_plot_labels(x_start, y_start, ymax, stat_max, stat_min, y_interval, x_left_boundary, title,
    numbers_first=false, precision=2, rotate=false, number_intervals = 10)
    x_adjustment = standard_font_size / 5 + 1

    x_start -= x_adjustment
    # now write them as right-aligned and horizontally
    stat_interval = stat_max-stat_min
    stat_break = stat_interval.to_f/number_intervals
    y_break = y_interval.to_f/number_intervals
    y=ymax
    current_stat_value = stat_min

    if rotate
      current_stat_value = stat_max
      stat_break = -stat_break
    end


    font_size = standard_font_size
    if font_size_multiple > 1
      font_size *= 0.85
    end

    num_x = @box_size/2.8
    letters_x = @box_size/2
    number_intervals.times do |curr_break|
      @canvas.g.translate(x_start, y_start).text(num_x, y) do |text|
        label = compress_number(current_stat_value, precision)
        text.tspan(label).styles(:font_size=>font_size/1.2, :text_anchor=>'end')
        # add tick mark
        current_stat_value += stat_break
        y -= y_break
      end
    end

    final_stat_value = stat_max
    if rotate
      final_stat_value = stat_min
    end

    @canvas.g.translate(x_start, y_start).text(num_x, y) do |text|
      text.tspan(compress_number(final_stat_value, precision)).styles(:font_size=>font_size/1.2, :text_anchor=>'end')
    end

    x_start -= calculate_coordinate(0.0005 * box_size * 18)

    # add labels
    # center around middle of the vertical axis
    if rotate
      title_adjust = 2 * @box_size + (font_size_multiple-1.0) * 1.2 * @box_size
      rotation = 90
    else
      title_adjust = -@box_size - (font_size_multiple-1.0) * 1.2 * @box_size
      rotation = -90
    end

    @canvas.g.translate(x_start, y_start).text(num_x+title_adjust, y_interval/2).rotate(rotation) do |text|
      text.tspan(title).styles(:font_size=>font_size*1.2, :text_anchor=>'middle')
    end

  end

  # adds title to plot
  def write_main_title(title_str, x_start, y_start, x_end, y_end, split_title,
    max_title_length)
    font_size = standard_font_size

    # divide into two lines of roughly equal size
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

    @canvas.g.translate(x_start, y_start).text((x_end-x_start)/2, (y_end-y_start)/2) do |text|
      text.tspan(title_str).styles(:font_size=>font_size*2, :text_anchor=>'middle')
    end
  end

  # draws lines across the plot to mark top and bottom
  def draw_plot_boundaries(x_start, x_end, y_start, size_mult=1.0)
    y_interval = @box_size*9*size_mult
    x_begin = @box_size * 0.7
    @canvas.g.translate(x_start, y_start) do |l|
      l.styles(:fill=>'none', :stroke_width=>1, :stroke=>'gray', :fill_opacity=>0.0)
      l.line(x_begin,y_interval,x_end-x_start,y_interval)
      l.line(x_begin,0,x_end-x_start,0)
    end
  end

  # draw red line across plot at indicated position
  def draw_red_line(x_start, y_start, x_end, stat_max, stat_min, p_thresh, rotate)
    x = @box_size
    xmin = x

    ymin = 0
    ymax = @box_size*9
    y_interval = ymax-ymin

    # set stat_min to be zero
    stat_min = 0
    stat_max = @functional_max

    stat_interval = stat_max - stat_min

    p_thresh = get_neg_log(p_thresh)
    y=((stat_max-p_thresh) / stat_interval) * y_interval
    if rotate
      y = y_interval - y
    end

    x_begin = @box_size * 0.7

    @canvas.g.translate(x_start, y_start) do |l|
      l.styles(:fill=>'none', :stroke_width=>1, :stroke=>'red', :fill_opacity=>0.1, :stroke_opacity=>0.5)
      l.line(x_begin, y, x_end-x_start, y)
    end

  end

  # places faint dashed line across plot at indicated value
  def draw_dashed(x_start, y_start, x_end, value, stat_max, stat_min)
    y_interval = @box_size*9
    x_begin = @box_size * 0.7
    stat_interval = stat_max - stat_min
    dash1 = @box_size/5+1
    dash2 = @box_size/2
    dash_array = Array.new
    dash_array << dash1
    dash_array << dash2
    y = ((stat_max-value) / stat_interval) * y_interval
    @canvas.g.translate(x_start, y_start) do |l|
      l.styles(:fill=>'none', :stroke_width=>1, :stroke=>'gray', :fill_opacity=>0.0, :stroke_dasharray=>dash_array)
      l.line(x_begin, y, x_end-x_start, y)
    end
  end

  # reduces number representation to appropriate precision
  def compress_number(num, precision=2)
    if precision > 0
      num_string = sprintf("%0.#{precision}f", num)
    else
      num_string = sprintf("%d", num.round)
    end
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
    return @font_size_multiple * @box_size/2 + 1
  end

  # creates an entry for the label with x coordinates marking the edges
  # of the map
  # fraction is the overall size relative to the original for a reduced file
  def map_label(loc1, loc2, name, rotate, fraction=1.0)

    if rotate #and loc1 < 0
      # convert coordinates to
      loc1 = (@maxx - pixels_from_coordinate(loc1)) * fraction
      point1 = Point.new(0, loc1)
      loc2 = (@maxx - pixels_from_coordinate(loc2)) * fraction
      point2 = Point.new(@maxy*fraction, loc2)
    else
      point1 = Point.new(pixels_from_coordinate(loc1)*fraction,0)
      point2 = Point.new(pixels_from_coordinate(loc2)*fraction, @maxy*fraction)
    end

    ipos = ImagePositions.new(point1, point2, name)
    @map_pos << ipos
  end

  # add labels for SNPs and physical distance
  # conversion -- store y coordinates for use in creating
  # an html map if needed
  def add_labels(chrom_num, snp_list, x_start, y_start, y_label_end, features_included=false,
    rotate=false, add_snp_locs=true, draw_dist=true)

    fraction = 1.0

    if @maxx > @maxy
      max_dim = @maxx
    else
      max_dim = @maxy
    end

    if max_dim > Maximum_html_image_x
      fraction = Maximum_html_image_x / max_dim
    else
      fraction = 1.0
    end
    x_text = @box_size

    start_x = x_text
    end_x = start_x
    # change for Sarah's longer titles
    y_text = y_label_end - y_start - @box_size/6

    font_size = standard_font_size

    rotate_angle = -90
    txt_anchor = 'start'
    box_adjust=0
    if rotate
      rotate_angle = 90
      txt_anchor = 'end'
      box_adjust = @box_size * Math.sqrt(2) *0.25 - (@box_size * Math.sqrt(2) *0.25)*0.25
    end

    x_map_pixel = pixels_from_coordinate(x_start-@box_size.to_f/2)

    snp_list.included_snps.each do |snp_index|
      snp = snp_list.snps[snp_index]
      @canvas.g.translate(x_start,y_start).text(x_text-box_adjust, y_text).rotate(rotate_angle) do |text|
          text.tspan(snp.name).styles(:font_size=>font_size, :text_anchor=>txt_anchor)
      end
      end_x = x_text

      x1 = x_start + x_text - @box_size.to_f/1.75
      x2 = x1 + @box_size
      map_label(x1, x2, snp.name, rotate, fraction)

      x_text = label_step(x_text)
    end

    # add the physical distance bar and
    # draw lines from that to the titles above
    if draw_dist
      min_pos = snp_list.get_min
      max_pos = snp_list.get_max
      interval_pos = max_pos - min_pos

      x_interval = end_x - start_x

      box_height = @box_size/3

      # start lower when larger text
      if @font_size_multiple <= 1
        box_y_start = @box_size*3
        y_text_line = @box_size*5.9
      else
        box_y_start = @box_size*4.5
        y_text_line = @box_size*7.4
      end

      font_style = 'italic'
      font_family = 'Arial'
      if snp_list.included_snps.length > 1
        @canvas.g.translate(x_start, y_start).text(end_x+@box_size/4-box_adjust,box_y_start-box_height).rotate(rotate_angle) do |text|
          text.tspan(snp_list.get_max.to_s).styles(:font_size=>font_size/1.3, :text_anchor=>txt_anchor, :font_style=>font_style, :font_family=>font_family)
        end
      else
      end

      @canvas.g.translate(x_start, y_start).text(start_x+@box_size/4-box_adjust, box_y_start-box_height).rotate(rotate_angle) do |text|
        text.tspan(snp_list.get_min.to_s).styles(:font_size=>font_size/1.3, :text_anchor=>txt_anchor, :font_style=>font_style, :font_family=>font_family)
      end

      @canvas.g.translate(x_start, y_start).text((end_x-start_x)/2+@box_size*1.22-box_adjust-@box_size.to_f/5,0) do |text|
        text.tspan("chr#{chrom_num.to_s}").styles(:font_size=>font_size/1.5, :text_anchor=>'middle', :font_style=>font_style, :font_family=>font_family)
      end

      # draw narrow rectangle with white fill
      @canvas.g.translate(x_start, y_start) do |chromosome|
        chromosome.styles(:fill=>'white', :stroke_width=>1, :stroke=>'black')
        chromosome.rect(x_interval, box_height, start_x, box_y_start)
      end

      x_text_line = @box_size + 1

      last_x_position = ((snp_list.snps[snp_list.included_snps.first].location.to_f - min_pos) / interval_pos) * x_interval + start_x

      # draw lines connecting fonts to the position
      if snp_list.included_snps.length > 1
        final_x_position = ((snp_list.snps[snp_list.included_snps.last].location.to_f - min_pos) / interval_pos) * x_interval + start_x

        snp_list.included_snps.each do |snp_index|
          snp = snp_list.snps[snp_index]
          # determine relative position
          x_end_position = ((snp.location.to_f - min_pos) / interval_pos) * x_interval + start_x

          # when enough room to label add the current basepair location
          # distance should be 1/2 of a box size
          if add_snp_locs
            if x_end_position - last_x_position > @box_size.to_f/1.8 and final_x_position - x_end_position > @box_size.to_f/1.8
              @canvas.g.translate(x_start, y_start).text(x_end_position+@box_size/4-box_adjust, box_y_start-box_height).rotate(rotate_angle) do |text|
                text.tspan(snp.location.to_s).styles(:font_size=>font_size/1.3, :text_anchor=>txt_anchor, :font_style=>font_style, :font_family=>font_family)
              end
              last_x_position = x_end_position
            end
          end

          @canvas.g.translate(x_start, y_start) do |pos_line|
            # draw a vertical line across the box
            pos_line.styles(:stroke=>'gray', :stroke_width=>1)
            pos_line.line(x_end_position, box_y_start, x_end_position, box_y_start+box_height)
            pos_line.line(x_text_line, y_text_line, x_end_position, box_y_start+box_height)
          end
          x_text_line = label_step(x_text_line)
        end
      else
        # only on SNP to draw line
        @canvas.g.translate(x_start, y_start) do |pos_line|
          pos_line.line(start_x, y_text_line, start_x, box_y_start+box_height)
        end
      end
    end

    # add space to end and return for marking next chromosome
    return label_step(x_start+end_x) - @box_size * 0.7

  end


  # draws gene information track
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

  # draws boundaries around plot and add title
  def group_plot_boundaries(x_start, y_start, num_rows, title)
    start_plot_x = @box_size * 0.7
    end_plot_x = @box_size * 14

    start_plot_y = @box_size
    end_y = @box_size * (num_rows+1)

    x_interval = end_plot_x - start_plot_x

    # draw box around plot area
    @canvas.g.translate(x_start, y_start) do |plot|
      plot.styles(:stroke_width=>1, :stroke=>'black')
      plot.line(start_plot_x, start_plot_y, end_plot_x, start_plot_y)
      plot.line(start_plot_x, end_y, end_plot_x, end_y)
      plot.line(start_plot_x, start_plot_y, start_plot_x, end_y)
      plot.line(end_plot_x, start_plot_y, end_plot_x, end_y)
    end

    # add title
    font_size = standard_font_size*1.2
    if @font_size_multiple > 1
      font_size = font_size * 0.8
    end

    @canvas.g.translate(x_start,y_start).text((start_plot_x + x_interval)/2, start_plot_y-@box_size/2) do |text|
      text.tspan(title).styles(:font_size=>font_size, :text_anchor=>'middle')
    end

    return start_plot_x, end_plot_x, start_plot_y, end_y

  end

  # calculate the values interval on a group summary plot
  def group_value_interval(max_value, min_value, offset)
    #determine min and max values based on the max to min differences
    value_difference = max_value - min_value

    # make increments scale to values portrayed
    if value_difference < 2
      max_value = max_value.ceil.to_f
      min_value = min_value.floor.to_f
    else
      if min_value.floor - offset > 0
        min_value = min_value.floor - offset
      else
        min_value = 0
      end
      max_value = max_value.floor + offset
    end

    return max_value, min_value, max_value - min_value
  end

  # labels bottom of plot with markers for determining values
  def group_label_points(max_value, min_value, start_plot_x, x_interval, x_start, y_start, end_y)
    increment = (max_value - min_value.to_f) / 4
    x_increment = x_interval/4
    current_value = min_value
    font_size = standard_font_size
    number_x = start_plot_x
    font_size = standard_font_size
    if @font_size_multiple > 1
      font_size = font_size * 0.8
    end
    5.times do |i|
      @canvas.g.translate(x_start, y_start).text(number_x,end_y+@box_size/2) do |text|
        text.tspan(current_value).styles(:font_size=>font_size, :text_anchor=>'middle')
      end
      number_x = number_x + x_increment
      current_value = current_value + increment
    end
  end

  # displays information related to a specific group
  # in a box plot
  def group_box_plot(grouplist, x_start, y_start, offset)
    start_plot_x, end_plot_x, start_plot_y, end_y = group_plot_boundaries(x_start, y_start,
      grouplist.groups.length, 'Box Plot')
    x_interval = end_plot_x - start_plot_x

    # plot names and values
    font_size = standard_font_size * 0.9

    current_y = start_plot_y + @box_size
    text_x = 0

    # determine min and max values
    min_value =1e50
    max_value =-1e50

    max_title_length = 0

    grouplist.groups.each do |group|
      if group.name.length > max_title_length
        max_title_length = group.name.length
      end
      @canvas.g.translate(x_start,y_start).text(@box_size/2.8, current_y-@box_size/3) do |text|
        text.tspan(group.name).styles(:font_size=>font_size, :text_anchor=>'end')
      end
      current_y = current_y + @box_size
      if group.values[:max].to_f > max_value
        max_value = group.values[:max].to_f
      end
      if group.values[:min].to_f < min_value
        min_value = group.values[:min].to_f
      end
    end

    max_value, min_value, value_interval = group_value_interval(max_value, min_value, offset)
    y_label_points = end_y + (@font_size_multiple-1.0) *@box_size.to_f/2
    group_label_points(max_value, min_value, start_plot_x, x_interval, x_start, y_start, y_label_points)

    current_y = start_plot_y + @box_size/2

    y_box_offset = @box_size.to_f/4
    y_box_height = @box_size.to_f/2

    # draw box plot values
    grouplist.groups.each do |group|
      # calculate x position
      x_median = (group.values[:median].to_f-min_value)/value_interval * x_interval + start_plot_x
      x_box_end = (group.values[:percent75].to_f-min_value)/value_interval * x_interval + start_plot_x
      x_box_start = (group.values[:percent25].to_f-min_value)/value_interval * x_interval + start_plot_x
      x_max = (group.values[:max].to_f-min_value)/value_interval * x_interval + start_plot_x
      x_min = (group.values[:min].to_f-min_value)/value_interval * x_interval + start_plot_x
      y_box_top = current_y - y_box_offset
      y_box_bottom = current_y + y_box_offset

      @canvas.g.translate(x_start, y_start) do |box|
        box.styles(:stroke_width=>1, :stroke=>group.colorstr, :fill=>'none')
        # draw box
        box.rect(x_box_end-x_box_start, y_box_height, x_box_start, y_box_top)

        # draw median line
        box.line(x_median, y_box_top, x_median, y_box_bottom)

        # draw max and min lines
        box.line(x_min, y_box_top, x_min, y_box_bottom)
        box.line(x_max, y_box_top, x_max, y_box_bottom)

        # draw whispers to max and min
        stroke_array = [2,2]
        box.line(x_min, current_y, x_box_start, current_y).styles(:stroke_dasharray=>stroke_array)
        box.line(x_box_end, current_y, x_max, current_y).styles(:stroke_dasharray=>stroke_array)

      end

      current_y = current_y+@box_size
    end

    x_font_adjustment = 0
    if font_size_multiple > 1
      x_font_adjustment = (end_plot_x - x_start) * 0.5
    end

    return x_start + end_plot_x + @box_size/3.2 * max_title_length *font_size_multiple

  end


  # displays information related to a specific group
  # such as Number of individuals or average phenotype
  # includes standard deviation lines
  def group_summary_plot(grouplist, x_start, y_start, title, datakey, offset)
    start_plot_x, end_plot_x, start_plot_y, end_y = group_plot_boundaries(x_start, y_start,
      grouplist.groups.length, title)
    x_interval = end_plot_x - start_plot_x

    # plot names and values
    font_size = standard_font_size * 0.9

    current_y = start_plot_y + @box_size
    text_x = 0

    # determine min and max values
    min_value =1e50
    max_value =-1e50

    max_title_length = 0

    grouplist.groups.each do |group|
      if group.name.length > max_title_length
        max_title_length = group.name.length
      end
      @canvas.g.translate(x_start,y_start).text(@box_size/2.8, current_y-@box_size/3) do |text|
        text.tspan(group.name).styles(:font_size=>font_size, :text_anchor=>'end')
      end
      current_y = current_y + @box_size
      if group.values[datakey].to_f + group.values[datakey + "_stddev"].to_f > max_value
        max_value = group.values[datakey].to_f + group.values[datakey + "_stddev"].to_f
      end
      if group.values[datakey].to_f + group.values[datakey + "_stddev"].to_f < min_value
        min_value = group.values[datakey].to_f - group.values[datakey + "_stddev"].to_f
      end
    end

    max_value, min_value, value_interval = group_value_interval(max_value, min_value, offset)

    group_label_points(max_value, min_value, start_plot_x, x_interval, x_start, y_start, end_y)

    current_y = start_plot_y + @box_size/2

    # draw points with standard deviation bars
    grouplist.groups.each do |group|
      # calculate x position
      x = (group.values[datakey].to_f-min_value)/value_interval * x_interval + start_plot_x
      @canvas.g.translate(x_start, y_start) do |point|
  	    point.styles(:fill=>group.colorstr, :stroke=>group.colorstr, :stroke_width=>1, :fill_opacity=>0.9)
	      point.circle(@box_size * Math.sqrt(2)*0.125, x, current_y)
	      # draw standard dev line
	      line_x1 = (group.values[datakey].to_f-group.values[datakey+"_stddev"].to_f-min_value)/value_interval * x_interval + start_plot_x
	      line_x2 = (group.values[datakey].to_f+group.values[datakey+"_stddev"].to_f-min_value)/value_interval * x_interval + start_plot_x
	      point.line(line_x1, current_y, line_x2, current_y)
      end
      current_y = current_y+@box_size
    end

    x_font_adjustment = 0
    if font_size_multiple > 1
      x_font_adjustment = (end_plot_x - x_start) * 0.5
    end

    return x_start + end_plot_x + @box_size/3.2 * max_title_length *font_size_multiple

  end

  # draws the group information grid -- shows whether the group has
  # information about the SNPs in the list
  def group_membership_plot(groupname, colorstr, snp_list, x_start, y_start, first_plot=false, show_one=false,
    grayscale=false)
    box_y_start = 0 #@box_size * Math.sqrt(2) *0.25 #@box_size/6

    x_text = @box_size * 0.7

    if grayscale
      stroke = 'black'
      box_y_start += @box_size/7
    else
      stroke = 'white'
    end
    
    snp_list.included_snps.each do |snp_index|
      if show_one
        break
      end
      snp = snp_list.snps[snp_index]
	    if snp.results[groupname] != nil
          @canvas.g.translate(x_start,y_start) do |check|
	          check.styles(:fill=>colorstr, :stroke=>stroke, :stroke_width=>1)
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

    font_size = standard_font_size * 0.9

    if first_plot
      namepcs = groupname.split /:/
      @canvas.g.translate(x_start,y_start).text(@box_size/2.8, font_size * 0.8) do |text|
        text.tspan(namepcs[0]).styles(:font_size=>font_size, :text_anchor=>'end')
      end
      if show_one
        box_y_start = 0
        box_x_start = @box_size
        @canvas.g.translate(x_start,y_start) do |check|
          check.styles(:fill=>colorstr, :stroke=>stroke, :stroke_width=>1)
          check.rect(@box_size * Math.sqrt(2) *0.25 * font_size_multiple, @box_size * Math.sqrt(2) *0.25 * font_size_multiple, box_x_start, box_y_start)
        end
      end
      
    end

  end


  # draws the grid showing LD values
  def draw_grid(snp_list, x_start, y_start, use_dprime=true)
    x_start = x_start - @box_size*0.7
    y = y_start+@box_size

    x = x_start-@box_size
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
    if @first_grid
      @canvas.g.text(label_x, y_text) do |text|
        text.tspan(grid_label).styles(:font_size=>font_size, :text_anchor=>'end')
      end
    end

    # draw plot legend below grid label
    if @first_grid
      draw_plot_legend(label_x, y_text+@box_size*1.5, use_dprime)
      @first_grid = false
    end

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
      start_index = snp_index +1
      end_index = snp_list.included_snps.length-1
      if snp_list.ld_scores.has_key?(snp.name)
        (start_index..end_index).each do |snp2_index|
          snp2 = snp_list.snps[snp2_index]
          snp_combo = snp_list.ld_scores[snp.name][snp2.name]
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

  def draw_plot_legend(x_start, y_start, is_dprime=true)

    labels_array = Array.new
      for i in 0..10
        labels_array << i * 0.1
      end

    current_x=0
    current_y=0

    legend_box_size = @box_size *0.75

    # draw boxes from starting point
    @canvas.g.translate(x_start,y_start) do |legend|
      labels_array.reverse_each do |box_value|
        color_string = get_color_from_score(box_value, 2.1, is_dprime)
        legend.rect(legend_box_size, legend_box_size, current_x, current_y).styles(:stroke=>color_string,
          :stroke_width=>1, :fill=>color_string)
        if is_dprime
          color_string = get_color_from_score(box_value, 0, is_dprime)
          legend.rect(legend_box_size, legend_box_size, current_x-legend_box_size, current_y).styles(:stroke=>color_string,
            :stroke_width=>1, :fill=>color_string)
        end
        current_y+=legend_box_size
      end
    end

    current_y=legend_box_size
    font_size = standard_font_size

    if is_dprime
      current_x = current_x - legend_box_size
    end

    labels_array.reverse_each do |label|
      @canvas.g.translate(x_start, y_start).text(current_x-2, current_y) do |text|
        text.tspan(label).styles(:font_size=>font_size/1.3, :text_anchor=>'end')
      end
      current_y+=legend_box_size
    end

    # place box around the legend
    current_y=0
    width = legend_box_size
    if is_dprime
      width = legend_box_size*2
    end
    @canvas.g.translate(x_start, y_start) do |legend|
      legend.rect(width, legend_box_size*11, current_x, current_y).styles(:stroke=>'black',
        :stroke_width=>1, :fill=>'none')
    end

    # write labels on top of columns for d-prime legend
    if is_dprime
      current_y = -2
      current_x = current_x + legend_box_size/2
      @canvas.g.translate(x_start, y_start).text(current_x, current_y) do | text|
        text.tspan("<=2").styles(:font_size=>font_size/1.7, :text_anchor=>'middle')
      end
      current_x = current_x + legend_box_size + legend_box_size/2
      @canvas.g.translate(x_start, y_start).text(current_x, current_y) do | text|
        text.tspan(">2").styles(:font_size=>font_size/1.7, :text_anchor=>'end')
      end
      current_x -= legend_box_size/2
      current_y = current_y - legend_box_size/1.5
      @canvas.g.translate(x_start, y_start).text(current_x-legend_box_size/2, current_y) do | text|
        text.tspan("LOD").styles(:font_size=>font_size/1.5, :text_anchor=>'middle')
      end
    end

  end


  # set up max and minimum so that 0 appears as indicated increment when
  # plotted with labels at every 0.1 of the interval
  def calculate_increments_include_zero(min, max)
    # increase max and min by 10% to allow for adjusting the location of zero
    max *= 1.1
    min = min - (min*0.1).abs

    interval = max-min

    new_max = max
    new_min = min
    # check for zero being outside range (pretty unlikely for beta values)
    if max < 0
      new_max = 0
    elsif min > 0
      new_min = 0
    else interval > 5
      zero_offset = max
      fraction = zero_offset/interval
      # convert to an evenly divided value from 0-1.0
      # will be 0,0.1,0.2,0.3...1.0

      increment_fraction = ((fraction*10).round.to_f)/10

      if interval > 5
          increment = (new_max / (increment_fraction*10)).ceil
          new_max = increment * (increment_fraction*10)
          new_min = -(increment * ((1-increment_fraction)*10))
      elsif interval > 4
          increment = 0.5
          new_max = increment * (increment_fraction*10)
          new_min = -(increment * ((1-increment_fraction)*10))
      elsif interval > 3
          increment = 0.4
          new_max = increment * (increment_fraction*10)
          new_min = -(increment * ((1-increment_fraction)*10))
      elsif interval > 2 #(here can have interval of 0.5 at least or smaller)
          increment = 0.3
          new_max = increment * (increment_fraction*10)
          new_min = -(increment * ((1-increment_fraction)*10))
      elsif interval > 1
          increment = 0.2
          new_max = increment * (increment_fraction*10)
          new_min = -(increment * ((1-increment_fraction)*10))
      else
        increment = ((new_max / (increment_fraction*10)) * 10).ceil/10.to_f
        new_max = increment * (increment_fraction*10)
        new_min = -(increment * ((1-increment_fraction)*10))
      end
    end

    increment = (new_max - new_min)/10
    return (new_max-new_min.to_f)/10, new_min, new_max
  end


  # calculates and returns the minimum and increment to use in plotting 10 increments
  # based on the minimum and max passed
  def calculate_increments(min, max)
    interval = max - min

    # use log10 to determine the proper interval
    power_of_ten = Math.log10(interval)
    if power_of_ten < 0
      power_of_ten -=1
    end

    #need check for when the power is exactly equal
    if power_of_ten.to_i == power_of_ten
      power_of_ten -= 1
    end

    increment = 10 ** power_of_ten.to_i

    # reduce increment if range too large
    if power_of_ten > 0
      if power_of_ten - power_of_ten.to_i < Math.log10(2.5)
        increment = increment.to_f/4
      elsif power_of_ten - power_of_ten.to_i < Math.log10(5)
        increment = increment.to_f / 2
      end
    else
      if power_of_ten - power_of_ten.to_i < Math.log10(0.25)
        increment = increment.to_f/4
      elsif power_of_ten - power_of_ten.to_i < Math.log10(0.5)
        increment = increment.to_f / 2
      end
    end

    # calculate the new minimum
    # start at zero and then increase increments until greater than min
    start = 0
    if min > 0
      10.times do |i|
        if start + increment > min
          break
        end
        start += increment
      end
    else # for min < 0
      while start > min
        start -= increment
      end
    end
    # if can't fit both in then need to adjust new_min until they fit
    while start+10*increment < max
      start += increment.to_f/10
    end

    if start + 10*increment < max
      print "PROBLEM\n"
    end


    # return increment, new_min, and new_max
    return increment, start, start+10*increment
  end

  # returns plot width in inches
  def calculate_plot_width(num_snps, right_padding)
    return num_snps * @box_size.to_f/97 + right_padding
  end

  #returns number of pixels based on the coordinate
  def pixels_from_coordinate(x)
    return pixels_per_coordinate * x
  end

  def set_pixels_per_coordinate(pixels, coordinates)
    @pixels_per_coordinate = pixels / coordinates.to_f
  end

  # returns value in x,y coordinate for the
  # size in inches passed
  def calculate_coordinate(size_in_inches)
    return size_in_inches * 144
  end

  def add_space_for_gene(size_in_inches, fraction)
    return size_in_inches + @box_size * fraction
  end

  def add_space_for_title(size_in_inches, fraction)
    return size_in_inches + @box_size * fraction
  end

  def add_space_for_labels(size_in_inches, fraction, distance_track_included)
    dist = @box_size * fraction + (@font_size_multiple - 1.0) * 0.3 * @box_size * fraction
    unless distance_track_included
      dist /= 3
    end
    return size_in_inches + dist
  end

  def add_space_for_snpinfo(size_in_inches, fraction)
    return size_in_inches + @box_size * fraction
  end

  def add_haplotype_space(size_in_inches, max_haplotype)
    return size_in_inches + 0.00325 * @box_size * max_haplotype + @box_size.to_f/200
  end

  def add_grid_size(size_in_inches, horizontal_side_in_inches, maximum_interactions, num_snps)
    if maximum_interactions < 10
      maximum_interactions = 15
    end
    xside_fraction = Math.sqrt(maximum_interactions / Math.sqrt(2) *
      (maximum_interactions / Math.sqrt(2)) / 2) / num_snps
    return  size_in_inches + horizontal_side_in_inches * xside_fraction  + 0.05
  end

end

##### End PlotWriter #####



#############################################################################
#
# Class FileReader -- Base class for file reading.  Contains functions useful
#                     to all other file readers.
#
#############################################################################
class FileReader

  # strips and splits the line
  # returns the data array that results
  def strip_and_split(line)
    line.strip!
    line.split(/\t/)
  end

end

#############################################################################
#
# Class SynthesisViewReader -- reads file for Synthesis-View format.
#                     The data are divided into chromosomes for purposes of display.
#
#############################################################################
class SynthesisViewReader<FileReader
  @@mysql_support = false

  # set column constants
  def initialize
    @snpid = nil
    @snpname = -1
    @chromnum = nil
    @location = nil
  end

  def self.mysql_support(tf)
    @@mysql_support = tf
  end

  # fills chromosome list with data from synthesis_view format file
  def read_synthesisview_file(synthfile, chromlist, glisthash, highlighted_group="")
    defaultkey = GroupList.get_default_name
    firstline = true
   begin
      File.open(synthfile, "r") do |file|

      # read in all lines
      lines = Array.new
      #file.each_line("\r") do |line|
      file.each_line("\n") do |line|
        line.each("\r") do |splitline|
          lines << splitline unless splitline =~ /^\n$/
        end
      end
      group_column, subgroup_column = check_subgroup(lines[0])  
      unless subgroup_column || group_column
        set_groups(glisthash, lines[0], defaultkey, highlighted_group)
      else
        set_groups_subgroup(glisthash, lines, defaultkey, group_column, subgroup_column, highlighted_group)
      end

      # if no SNP name in file set SNP name to be same as SNPID
      if @snpname == -1
        @snpname = @snpid
      end
        
      # if there is built-in MySQL support for looking up position
      # information
      snp_positions = get_locations(lines, synthfile)

      lines.each do |line|
        # skip blank lines
        if line !~ /\w/
          next
        end

        # read column headers to create groups
        if firstline
          firstline = false
          next
        end

        # split line and check whether need to create new chromosome
        data = strip_and_split(line)

        # if no location information present skip
        if !snp_positions.has_key?(data[@snpid])
          next
        end

        chrnum = snp_positions[data[@snpid]]['chr']

        if (chromosome = chromlist.get_chrom(chrnum.to_i)) == nil
          # create a new SNP
          newchrom = Chromosome.new(chrnum.to_i)
          chromlist.add_chrom(newchrom)
          chromosome = chromlist.get_chrom(chrnum.to_i)
        end

        # location will be first num listed
        location = snp_positions[data[@snpid]]['pos']

        location = location.to_i
        snp = chromosome.snp_list.get_snp(data[@snpname])
        
        unless snp
          snp = SNP.new(data[@snpname], location)
          chromosome.snp_list.add_snp(snp)
          snp = chromosome.snp_list.get_snp(data[@snpname])
        end

        # add results to SNP
        if !subgroup_column and !group_column
          glisthash.each_value do |glist|
            glist.groups.each do |group|
              read_group_data(snp, group, data)
            end
          end
        elsif subgroup_column and group_column
          glisthash.each do |name, glist|
            if name == data[subgroup_column]
              glist.groups.each do |group|
                if group.name == data[group_column] + ":" + data[subgroup_column]
                  read_group_data(snp, group, data)
                end
              end
            end
          end          
        elsif group_column
          glisthash.each_value do |glist|
            glist.groups.each do |group|
              if group.name == data[group_column]
                read_group_data(snp, group, data)
              end
            end
          end
         end
      end
     end
    rescue Exception => e
      puts e
      exit(1)
    end
  end


# searches a MySQL database for locations of SNPs that do not
# have chromosome and position information in input file
def get_locations(lines, filename)
  positions = Hash.new
  missing_pos = Array.new
  first = true
  lines.each do |line|
    if first
      first = false
      next
    end
    data = strip_and_split(line)
    snpid = data[@snpid]
    if data[@location] =~ /\d/
      positions[snpid] = Hash.new
      positions[snpid]['chr'] = data[@chromnum]
      positions[snpid]['chr'] =~ /(\d+)/
      positions[snpid]['chr'] = $1
      positions[snpid]['pos'] = data[@location]
    elsif data[@snpid] =~ /\w/
      # missing so will try to get information from database
      missing_pos << snpid
    end
  end

  # if any missing try database for location if have mysql support
  if @@mysql_support and !missing_pos.empty?
    finder = SNPpos.new
    snpsinfo = finder.get_pos(missing_pos)
    still_missing = Array.new
    missing_pos.each do |rsid|
      if snpsinfo.has_key?(rsid) and snpsinfo[rsid]['chr'] != nil
        positions[rsid] = snpsinfo[rsid]
      else
        still_missing << rsid
      end
    end

    missing_pos = still_missing
  end

  if !missing_pos.empty?
    write_missing_log(filename, missing_pos)
  end

  return positions
end


# writes a log file of all SNPs that are skipped because no location
# information is available
def write_missing_log(filename, missing_pos)
  logname = filename + ".miss.log"
  File.open(logname, 'w') do |fwrite|
    missing_pos.each do |missing|
    fwrite.puts "#{missing}"
    end
  end
  print "Produced file listing SNPs with no position information:  #{logname}\n"
end



# reads in group info using the data array passed
# and stores in SNP
def read_group_data(snp, group, data)

  if data[group.pcol] =~ /\d/ and data[group.pcol].to_f > 0.0
    betaval = nil
    mafval = nil
    sampsizeval = nil
    orval = lowerci = upperci = cases = controls = studynum = nil
    powernum = cafcases = cafcontrols = betauci = betalci = nil
    add_columns = Hash.new
    
    if group.betacol > -1
      betaval = data[group.betacol]
    end
    if group.mafcafcol > -1
      mafval = data[group.mafcafcol]
    end
    if group.Ncol > -1
      sampsizeval = data[group.Ncol]
    end
    if group.orcol > -1
      orval = data[group.orcol]
    end
    if group.lcicol > -1
      lowerci = data[group.lcicol]
    end
    if group.ucicol > -1
      upperci = data[group.ucicol]
    end
    if group.casescol > -1
      cases = data[group.casescol]
    end
    if group.controlscol > -1
      controls = data[group.controlscol]
    end
    if group.casescol > -1
      cafcases = data[group.cafcasescol]
    end
    if group.controlscol > -1
      cafcontrols = data[group.cafcontrolscol]
    end
    if group.studycol > -1
      studynum = data[group.studycol]
    end
    if group.powercol > -1
      powernum = data[group.powercol]
    end
    if group.betaucicol > -1
      betauci = data[group.betaucicol]
    end
    if group.betalcicol > -1
      betalci = data[group.betalcicol]
    end

    group.additional_cols.each_pair {|key,index| add_columns[key]=data[index] if data[index]}
    
    snp.add_result(group.name, data[group.pcol], betaval, mafval, sampsizeval, orval, lowerci, upperci,
      cases, controls, cafcases, cafcontrols, studynum, powernum, betauci, betalci, add_columns)

  end

end


# returns the indexes for the group and subgroup columns
def check_subgroup(headerline)
  data = strip_and_split(headerline)
  subgroup_column = group_column = nil
  data.each_with_index do |header, i|
    if header =~ /^subgroup$/i
      subgroup_column = i
    elsif header =~ /^group$/i
      group_column = i
    end
  end
  return group_column, subgroup_column
end

# creates groups and sets the columns in all groups based on header line
# this version uses a subgroup column
def set_groups_subgroup(glisthash, lines, defaultkey, groupcol, subgroupcol, highlighted_group="")
  # set up hash for holding columns for main group names
  subnames = Hash.new
  subgrouporder = Array.new
  groupnames = Hash.new
  groupnameorder = Array.new
  
  for i in (1..lines.length-1)
    data = strip_and_split(lines[i])
    if subgroupcol
      unless subnames.has_key?(data[subgroupcol])
        subnames[data[subgroupcol]]=1
        subgrouporder << data[subgroupcol]
      end
    end
    if groupcol
      unless groupnames.has_key?(data[groupcol])
        groupnames[data[groupcol]]=1
        groupnameorder << data[groupcol]
      end
    end
  end
  
  groups = Hash.new
  grouporder = Array.new
  groupkeys = Array.new 
  
  # construct the groups
  groupnameorder.each do |gname|
    gname == highlighted_group ? highlighted = true : highlighted = false
    if subnames.empty?
      key = gname
      groups.store(key, Group.new(key, highlighted,GroupList.get_next_color))
      grouporder << key
      groupkeys << key
    else
      subnames.each_key do |sname|  
        key = gname + ':' + sname
        groups.store(key, Group.new(key, highlighted))
        grouporder << key
        groupkeys << key
      end
    end
  end
  
  mafcoltitle = 'MAF'
  headers = strip_and_split(lines[0])

  # create the groups using the headers 
  headers.each_with_index do |header, i|
    
    if header =~ /snp\s*id/i || header=~ /snp_id/i || header =~ /^snp$/i
      @snpid = i
    elsif header =~ /snp\s*name/i
      @snpname = i
    elsif header =~ /chromosome|CHR/i
      @chromnum = i
    elsif header =~ /location/i
      @location = i
    elsif header =~ /^subgroup$/i || header =~ /^group$/i  # skip if no _ to mark name
      next
    else
      header.strip!
      column_type = header
      if column_type =~ /pval/i
        groupkeys.each {|key| groups[key].pcol = i}
      elsif column_type =~ /beta_uci|betauci/
        groupkeys.each {|key| groups[key].betaucicol = i}
      elsif column_type =~ /beta_lci|betalci/
        groupkeys.each {|key| groups[key].betalcicol = i}
      elsif column_type =~ /beta/i or column_type =~ /^es$/
        groupkeys.each {|key| groups[key].betacol = i}
      elsif column_type =~ /^n$/i
        groupkeys.each {|key| groups[key].Ncol = i}
      elsif column_type =~ /cafcases/i
        groupkeys.each {|key| groups[key].cafcasescol = i}
      elsif column_type =~ /cafcontrols/i
        groupkeys.each {|key| groups[key].cafcontrolscol = i}
      elsif column_type =~ /^maf|caf$/i
        if column_type =~ /caf/i
          mafcoltitle = 'CAF'
        end
        groupkeys.each {|key| groups[key].mafcafcol = i}
      elsif column_type =~ /^or$/i
        groupkeys.each {|key| groups[key].orcol = i}
      elsif column_type =~ /^upper_ci|uci$/i
        groupkeys.each {|key| groups[key].ucicol = i}
      elsif column_type =~ /lower_ci|lci/i
        groupkeys.each {|key| groups[key].lcicol = i}
      elsif column_type =~ /cases/i
        groupkeys.each {|key| groups[key].casescol = i }
      elsif column_type =~ /controls/i
        groupkeys.each {|key| groups[key].controlscol = i}
      elsif column_type =~ /study/i
        groupkeys.each {|key| groups[key].studycol = i}
      elsif column_type =~ /^power$/i
        groupkeys.each {|key| groups[key].powercol = i}
      else
        groupkeys.each {|key| groups[key].additional_cols[column_type] = i}       
      end
    end
  end

  unless @snpid and @location and @chromnum
    puts "ERROR:  Need SNPID, CHROMOSOME, and LOCATION columns in input file"
    exit
  end
  # add groups to the grouplist
  grouporder.each do |g|
    namepcs = g.split /:/

    # add to default grouplist
    if namepcs.length == 1
      if !glisthash.has_key?(defaultkey)
        glisthash[defaultkey] = GroupList.new
        glisthash[defaultkey].mafcoltitle = mafcoltitle
      end
      glisthash[defaultkey].add_group(groups[g])
    else
      if !glisthash.has_key?(namepcs[1])
        glisthash[namepcs[1]] = GroupList.new
        glisthash[namepcs[1]].mafcoltitle = mafcoltitle
      end
      glisthash[namepcs[1]].add_group(groups[g])
    end
  end

  # need to match all colors when multiple grouplists
  if glisthash.length > 1
    # determine number of unique groups
    unique_names = Hash.new
    glisthash.each_value do |glist|
      glist.grouphash.each_key do |name|
        namepcs = name.split /:/
        unique_names[namepcs[0]] = 1
      end
    end
    colorhash = Hash.new
    unique_names.each_key do |name|
      colorhash[name] = GroupList.get_next_color
    end

    glisthash.each_value do |glist|
      glist.grouphash.each do |name, group|
        namepcs = name.split /:/
        group.colorstr = colorhash[namepcs[0]]
      end
    end

  end
  
end

# sets groups in the grouplist based on header columns
# from Synthesis-View file
def set_groups(glisthash, line, defaultkey, highlighted_group="")
  data = strip_and_split(line)
  groups = Hash.new
  grouporder = Array.new

  mafcoltitle = 'MAF'

  data.each_with_index do |header, i|
    if header =~ /snp\s*id/i || header =~ /snp_id/i || header =~ /^snp$/i
      @snpid = i
    elsif header =~ /snp\s*name/i
      @snpname = i
    elsif header =~ /chromosome|CHR/i
      @chromnum = i
    elsif header =~ /location/i
      @location = i
    elsif header !~ /:/  # skip if no _ to mark name
      next
    else # information for a group
      # split name on ':'
      header.strip!
      pcs = header.split(/:/)
      pcs[0].upcase!
      column_type = pcs[1]
      # get or create this group
      if pcs.length == 2
        if !groups.has_key?(pcs[0])
          pcs[0] == highlighted_group ? highlighted = true : highlighted = false
          #groups.store(pcs[0], Group.new(pcs[0],-1,-1,-1, GroupList.get_next_color,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,highlighted))
          groups.store(pcs[0], Group.new(pcs[0], highlighted, GroupList.get_next_color))
          grouporder << pcs[0]
        end
        currgroup = groups.fetch(pcs[0])
        column_type = pcs[1].strip
      else # create subgroup if necessary and set current group to be this subgroup
        pcs[1].upcase!
        key = pcs[0] + ':' + pcs[1]
        key == highlighted_group ? highlighted = true : highlighted = false
        if !groups.has_key?(key)
          #groups.store(key, Group.new(key, -1, -1,-1, nil, -1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,highlighted))
          groups.store(key, Group.new(key, highlighted))
          grouporder << key
        end
        currgroup = groups.fetch(key)
        column_type = pcs[2].strip
      end
      
      if column_type =~ /pval/i
        currgroup.pcol = i
      elsif column_type =~ /beta_uci|betauci/
        currgroup.betaucicol = i
      elsif column_type =~ /beta_lci|betalci/
        currgroup.betalcicol = i
      elsif column_type =~ /beta/i or column_type =~ /^es$/
        currgroup.betacol = i
      elsif column_type =~ /^n$/i
        currgroup.Ncol = i
      elsif column_type =~ /cafcases/i
        currgroup.cafcasescol = i
      elsif column_type =~ /cafcontrols/i
        currgroup.cafcontrolscol = i
      elsif column_type =~ /^maf|caf$/i
        currgroup.mafcafcol = i
        if column_type =~ /caf/i
          mafcoltitle = 'CAF'
        end
      elsif column_type =~ /or/i
        currgroup.orcol = i
      elsif column_type =~ /upper_ci|uci/i
        currgroup.ucicol = i
      elsif column_type =~ /lower_ci|lci/i
        currgroup.lcicol = i
      elsif column_type =~ /cases/i
        currgroup.casescol = i
      elsif column_type =~ /controls/i
        currgroup.controlscol = i
      elsif column_type =~ /study/i
        currgroup.studycol = i
      elsif column_type =~ /^power$/i
        currgroup.powercol = i
      else
        currgroup.additional_cols[column_type]=i
      end
    end
  end

  unless @snpid and @location and @chromnum
    puts "ERROR:  Need SNPID, CHROMOSOME, and LOCATION columns in input file"
    exit
  end
  
  # add groups to the grouplist
  grouporder.each do |g|
    namepcs = g.split /:/

    # add to default grouplist
    if namepcs.length == 1
      if !glisthash.has_key?(defaultkey)
        glisthash[defaultkey] = GroupList.new
        glisthash[defaultkey].mafcoltitle = mafcoltitle
      end
      glisthash[defaultkey].add_group(groups[g])
    else
      if !glisthash.has_key?(namepcs[1])
        glisthash[namepcs[1]] = GroupList.new
        glisthash[namepcs[1]].mafcoltitle = mafcoltitle
      end
      glisthash[namepcs[1]].add_group(groups[g])
    end
  end

  # need to match all colors when multiple grouplists
  if glisthash.length > 1
    # determine number of unique groups
    unique_names = Hash.new
    glisthash.each_value do |glist|
      glist.grouphash.each_key do |name|
        namepcs = name.split /:/
        unique_names[namepcs[0]] = 1
      end
    end
    colorhash = Hash.new
    unique_names.each_key do |name|
      colorhash[name] = GroupList.get_next_color
    end

    glisthash.each_value do |glist|
      glist.grouphash.each do |name, group|
        namepcs = name.split /:/
        group.colorstr = colorhash[namepcs[0]]
      end
    end

  end

end

# reads the group total file headers and sets the indexes in a hash that is returned
def read_group_file_headers(line)
  data = strip_and_split(line)

  header_cols = Hash.new

  data.each_with_index do |header, i|
    if header =~ /group/i
      header_cols[:group] = i
    elsif header =~ /color/i
      header_cols[:color] = i
    elsif header =~ /avg pheno/i
      header_cols[:phenoavg] = i
    elsif header =~ /std dev/i
      header_cols[:stddev] = i
    elsif header =~ /samples/i
      header_cols[:samples] = i
    elsif header =~ /Median/i
      header_cols[:median] = i
    elsif header =~ /25/i
      header_cols[:percent25] = i
    elsif header =~ /75/i
      header_cols[:percent75] = i
    elsif header =~ /Min/i
      header_cols[:min] = i
    elsif header =~ /Max/i
      header_cols[:max] = i
    end
  end

  return header_cols
end

  # fills group info about totals and the phenotype averages and standard deviations
  def read_group_total_file(filename, grouplisthash)
    firstline = true

    header_cols = nil

    begin
      File.open(filename, "r") do |file|
      file.each_line do |oline|

        oline.each("\r") do |line|
        next if line =~ /^\n$/

        # skip blank lines
        if firstline
          header_cols = read_group_file_headers(line)
          firstline = false
          next
        end

        if line !~ /\w/
          next
        end

        # split line and strip whitespace
        data = strip_and_split(line)
        # group name is first column
        data[header_cols[:group]].upcase!
        # find appropriate group
        namepcs = data[header_cols[:group]].split /:/
        group = nil
        if namepcs.length == 1
          group = grouplisthash[GroupList.get_default_name].grouphash[data[header_cols[:group]]]
        else
          group = grouplisthash[namepcs[1]].grouphash[data[header_cols[:group]]]
        end
        if header_cols.has_key?(:samples)
          group.values['num'] = data[header_cols[:samples]]
          group.values['num_stddev'] = 0
        end

        if header_cols.has_key?(:phenoavg)
          group.values['pheno_avg'] = data[header_cols[:phenoavg]]
        end
        if header_cols.has_key?(:stddev)
          group.values['pheno_avg_stddev'] = data[header_cols[:stddev]]
        end
        if header_cols.has_key?(:color)
          group.colorstr = data[header_cols[:color]]
        end

        if header_cols.has_key?(:median)
          if !(header_cols.has_key?(:percent25) and header_cols.has_key?(:percent75) and
            header_cols.has_key?(:min) and header_cols.has_key?(:max))
            raise "Group file must have all columns for box plot: Median 25% 75% Min Max"
          end
          group.values[:median] = data[header_cols[:median]]
          group.values[:percent25] = data[header_cols[:percent25]]
          group.values[:percent75] = data[header_cols[:percent75]]
          group.values[:min] = data[header_cols[:min]]
          group.values[:max] = data[header_cols[:max]]
        end
      end
      end
    end
    rescue Exception => e
      puts e
      exit(1)
    end
  end

  #reads abbreviation file
  def read_abbrev_file(abbrevfile)
    definitions = Hash.new
    begin
      File.open(abbrevfile) do |file|
        file.each_line do |oline|
          oline.each("\r") do |line|
          next if line =~ /^\n$/
          if line !~ /\w/
            next
          end
          data = strip_and_split(line)
          definitions[data[0]] = data[1]
          end
        end
      end

    rescue Exception => e
      puts e
      exit(1)
    end
    return definitions
  end

  # reads LD file and returns array
  # of LD combinations
  def read_ld_output_file(ldfile)
    snp_lds = Array.new

    begin
      File.open(ldfile, "r") do |file|
        file.each_line  do |oline|
          oline.each("\r") do |line|
          next if line =~ /^\n$/
          if line =~ /LOD/
            next
          end

          data=strip_and_split(line)         
          snp_lds << SNPCombination.new(data[0], data[1], data[2], data[3], data[4])
          end
        end
      end
    rescue Exception => e
      puts e
      exit(1)
    end
    return snp_lds

  end

  # reads gene file information and returns array of genes
  def read_gene_file(genefile)

    first_line = true
    newgenes = Array.new

    begin
      File.open(genefile, "r") do |file|
        file.each_line do |oline|
          oline.each("\r") do |line|
          next if line =~ /^\n$/
          # skip header line
          if first_line
            first_line = false
            next
          elsif line !~ /\w/
            next
          end
          data=strip_and_split(line)
          ngene = Gene.new
          # check for presence of + in last column
          positive = false
          if data[4] =~ /\+/
            positive = true
          end
          ngene.fill_info(data[2].to_f, data[3].to_f, positive, data[0], data[1])
          newgenes << ngene
          end
        end
      end
    rescue Exception => e
      puts e
      exit(1)
    end
    return newgenes
  end

end

#############################################################
##
## Main execution begins here
##
#############################################################

SynthesisViewReader.mysql_support(true)

begin
  require '/home/dudek/synthesisview/SNPpos'
  rescue Exception => e
  puts
  puts e
  puts "\nNo MySQL support"
  puts
  SynthesisViewReader.mysql_support(false)
end


begin
  require './synthManhattan.rb'
  manhattan_support = true
rescue Exception => e
  puts e
  puts "\nNo Manhattan plot support"
  manhattan_support = false
end

options = Arg.init_options
command_args = Array.new
command_args.replace(ARGV)

Arg.parse(ARGV, options)

if options.option_logfile
  newoptions = Arg.read_options(options.option_logfile)
  Arg.parse(command_args, newoptions)
  options = newoptions
end
Arg.write_options(options, options.out_name + ".opt.log")

if options.grayscale
  GroupList.set_grayscale
  Shape.set_grayscale
end

# increase dpi if highres is selected
if options.highres
  RVG::dpi=300
end


chromlist = ChromosomeList.new
chromlist.set_columns(options.additional_columns)

synth_reader = SynthesisViewReader.new

grouplisthash = Hash.new

synth_reader.read_synthesisview_file(options.eaglefile, chromlist, grouplisthash, 
  options.highlighted_group)
chromlist.sort_chroms!
if options.rotate
  chromlist.reverse!
end

if options.manhattanfile
  manhattan_reader = Manhattan::ManhattanReader.new
  manhattan_results = manhattan_reader.read_manhattan(options.manhattanfile)
end

# make all snps included -- for use in plot -- later may wish
# to alter this to exclude some snps
# in ld_plus only snps that are listed with LD are part of the plot
chromlist.chromarray.each do |chrnum|
  chromlist.chromhash[chrnum].snp_list.sort_snps
 if options.rotate
   chromlist.chromhash[chrnum].snp_list.reverse_snps!
 end
  chromlist.chromhash[chrnum].snp_list.make_all_included
end

# read in ld scores
if(options.ld_file)
  ld_scores = synth_reader.read_ld_output_file(options.ld_file)
  ld_scores.each do |combo|
    chromlist.add_ld_score(combo)
  end
end

chromlist.set_max_min(grouplisthash)
total_snps = chromlist.get_nsnps

# box_size is basic unit of size in the plot
# it is adjusted for total number of snps in the set
box_size = case
           when total_snps > 500 then 2
           when total_snps > 350 then 4
           when total_snps > 275 then 6
           when total_snps > 200 then 8
           when total_snps > 175 then 10
           when total_snps > 150 then 12
           when total_snps > 125 then 14
           when total_snps > 100 then 16
           when total_snps > 80 then 18
           when total_snps > 60 then 20
           when total_snps > 40 then 22
           else 24
           end

if total_snps < 12
  total_snps = 12
  vertical_totals = true
end

# writer actually draws the plots
writer = PlotWriter.new(box_size)

if options.largetext
  writer.font_size_multiple = 1.8
  max_title_length = (1.75 * total_snps + 8).to_i
else
  writer.font_size_multiple = 1.0
  max_title_length = (3 * total_snps + 12).to_i
end

if options.title.length > max_title_length
  split_title = true
else
  split_title = false
end

xside_end_addition = 0.005 * box_size * 2
if box_size > 20
  xside_end_addition *=2.5
end

x_start = 0

if options.groupfile
  synth_reader.read_group_total_file(options.groupfile,grouplisthash)
end

if options.abbrevfile
  definitions = synth_reader.read_abbrev_file(options.abbrevfile)
  definitions.each_pair do |key, value|
    if grouplisthash[GroupList.get_default_name].get_group_by_name(key)
      grouplisthash[GroupList.get_default_name].get_group_by_name(key).fullname = value
    end
  end
end

# store x locations for forest legend lines of information
x_forest_legend_lines = Array.new

if grouplisthash.length > 0
  max_label_length=0
  last_forest_pos=0
  grouplisthash.each do |key, grouplist|
    x_total_forest_legend_addition = 0
    grouplist.groups.each do |group|
      if group.name.length > max_label_length
        max_label_length = group.name.length
      end
      if options.forest_legend
        xside_inc = box_size/2 * 1 * 0.01666 + (writer.font_size_multiple-1.0)*box_size/4*1*0.01667
        x_forest_legend_lines << writer.calculate_coordinate(xside_inc) + last_forest_pos
        last_forest_pos = x_forest_legend_lines.last
        x_total_forest_legend_addition += x_forest_legend_lines.last
        xside_end_addition += xside_inc
      end
    end

  end
  max_label_length += 10
  xleftside_addition = (0.001 * box_size * max_label_length) + 0.012 * box_size
  xleftside_addition += (writer.font_size_multiple-1.0) * 0.5 *  xleftside_addition
  x_start = writer.calculate_coordinate(xside_end_addition)
  xside_end_addition = xside_end_addition + xleftside_addition
  if !options.rotate
    x_start = writer.calculate_coordinate(xleftside_addition)
  end
end


gene_rows = 0
if options.genefile
  genes = synth_reader.read_gene_file(options.genefile)
  genes.each do |g|
    chromlist.add_gene(g)
  end
  chromlist.sort_genes
  gene_rows = chromlist.num_gene_rows
end

yside = 0

addition_total_snps = 0
circle_legend_below = false
if options.circlecasecon and !options.rotate and (options.dprime or
  options.rsquared or options.plot_study)
  addition_total_snps = 2
elsif !options.rotate
  circle_legend_below = true
end

xside = writer.calculate_plot_width(total_snps+addition_total_snps, xside_end_addition)
xmax = writer.calculate_coordinate(xside)

ymax = 0

y_title_start = ymax
title_height_mult=1
# add space for the title
if split_title
  title_height_mult = 2
end
yside = writer.add_space_for_title(yside, 0.032*title_height_mult)
ymax = writer.calculate_coordinate(yside)

# add height for rows of information on genes
# will appear above the location track
y_gene_start = ymax
if gene_rows > 0
  if vertical_totals
    offset = 0.5 * y_title_start
  else
    offset = 0
  end

  multiplier = 1
  if gene_rows == 1
    multiplier = 1.5
  end

  yside = yside + box_size/4 * 2 * 0.01667 * gene_rows * multiplier + offset
  ymax = writer.calculate_coordinate(yside)
end

if options.manhattanfile
 #add space for manhattan plot
  y_manhattan_start = ymax
  yside = writer.add_space_for_labels(yside, 0.065)
  ymax = writer.calculate_coordinate(yside)
  y_manhattan_end = ymax - ((ymax-y_manhattan_start)*0.1)
end

y_label_start = ymax

max_snp_name = chromlist.get_max_snp_name

# add more space in general if it is large text option
modification = 0

if options.largetext
  modification = 0.005
  #mod_mult = 0.0080
  mod_mult = 0.0075
else
  modification = 0
  #mod_mult = 0.0038
  mod_mult = 0.0034
end

if max_snp_name > 7
 modification += (max_snp_name - 7) * mod_mult
 if options.include_dist_track 
   if options.largetext
   modification /= 2.5
    else
   modification /= 1.5
   end
 end
end
yside = writer.add_space_for_labels(yside, 0.065 + modification, options.include_dist_track)

ymax = writer.calculate_coordinate(yside)
y_label_end = ymax

# add space for snpinfo boxes
ygroup_boxes_start = Array.new

ystart_block_stat_box = ymax

# add if need to show which groups have data for each SNP
if grouplisthash.length == 1 and !options.rotate
  grouplist = grouplisthash[GroupList.get_default_name]
  if !chromlist.results_complete?(grouplist.groups.length)
    grouplist.groups.length.times do |i|
      ygroup_boxes_start << ymax  # + box_size/8
      yside = yside + box_size/4 * 1 * 0.01667 + (writer.font_size_multiple-1.0)*box_size/4*0.5*0.01667#0.0125
      ymax = writer.calculate_coordinate(yside)
    end
  end
end

yside = yside + box_size/5 * 1 * 0.01667
ymax = writer.calculate_coordinate(yside)

yforest_legend_start = ymax

# when only one group list do standard plot
if grouplisthash.length == 1 and !options.rotate
  total_stat_boxes = 0
  if options.plot_pval
    total_stat_boxes +=1
  end

  grouplist = grouplisthash[GroupList.get_default_name]
  if grouplist.plot_betas? and options.plot_beta
    total_stat_boxes +=1
  end

  if options.plot_power
    total_stat_boxes +=1
  end

  if grouplist.plot_maf?
    total_stat_boxes +=1
  end

  if grouplist.plot_sample_sizes?
    total_stat_boxes +=1
  end

  if options.casecontrol
    total_stat_boxes += 2
  end

  if options.circlecasecon
    total_stat_boxes += 1
  end

  if options.cafcasecontrol
    total_stat_boxes +=2
  end

  if options.circlecafcasecontrol
    total_stat_boxes +=1
  end

  if grouplist.plot_oddsratio?
    total_stat_boxes +=1
  end
  
  # add any additional columns to plot
  total_stat_boxes += options.additional_columns.length

elsif options.rotate
  total_stat_boxes = grouplisthash[GroupList.get_default_name].groups.length
  if grouplisthash[GroupList.get_default_name].plot_pvals? and options.plot_pval
    total_stat_boxes +=1
  end

  if grouplisthash[GroupList.get_default_name].plot_betas?  and options.plot_beta
    total_stat_boxes +=1
  end

  if options.plot_power
    total_stat_boxes +=1
  end

  if grouplisthash[GroupList.get_default_name].plot_maf?
    total_stat_boxes +=1
  end
  if options.casecontrol
    total_stat_boxes += 2
  end
  if options.circlecasecon
    total_stat_boxes += 1
  end
  if options.cafcasecontrol
    total_stat_boxes +=2
  end

  if options.circlecafcasecontrol
    total_stat_boxes +=1
  end
  
  total_stat_boxes += options.additional_columns.length
  
else # multiple group lists usually for ethnicity so one plot per ethnicity
  total_stat_boxes = grouplisthash.length

  # iterate through groups and add a box for each beta that needs to be plotted
  grouplisthash.each_value do |group|
    if group.plot_betas? and options.plot_beta
      total_stat_boxes += 1
    end
  end

end

ystart_stat_boxes = Array.new

ystart_haplotypes = 0
yend_haplotypes = 0

total_stat_boxes.times do |i|
  ystart_stat_boxes << ymax + box_size/2
  yside = yside + box_size * 2 * 0.0355
  ymax = writer.calculate_coordinate(yside)
end

if circle_legend_below
  yside = yside + box_size * 2 * 0.0355/8
  ymax = writer.calculate_coordinate(yside)
end

# check for 1/2 sized boxes (i.e. study)
ystart_half_boxes = Array.new
if grouplisthash.length == 1
  if grouplisthash[GroupList.get_default_name].plot_study? and options.plot_studynum
    ystart_half_boxes << ymax + box_size/2
    yside = yside + box_size * 0.0360
    ymax = writer.calculate_coordinate(yside)
  end
end

yend_block_stat_box = ymax

y_dprime_start = ymax

# get largest depth of the ld interactions
# use to set the size needed for the grids
# displaying LD results
maximum_interactions = chromlist.get_max_interactions+1
max_num_included_snps = chromlist.get_max_included_snps

# create a new image for displaying the data
# determine how much of the xside is used
if(options.dprime)
  yside = writer.add_grid_size(yside, xside, maximum_interactions, total_snps)
  ymax = writer.calculate_coordinate(yside)
end
y_rsquared_start = ymax

# double all vertical sizes for second grid (r-squared)
if(options.rsquared)
  yside = writer.add_grid_size(yside, xside, maximum_interactions, total_snps)
  ymax = writer.calculate_coordinate(yside)
end

y_group_legend = ymax
y_group_legend_rows = Array.new
# add if need to show which groups have data for each SNP
if grouplisthash.length == 1 and !options.rotate
  grouplist = grouplisthash[GroupList.get_default_name]
  if chromlist.results_complete?(grouplist.groups.length)
    grouplist.groups.length.times do |i|
      y_group_legend_rows << ymax  # + box_size/8
      yside = yside + box_size/4 * 1 * 0.01666 + (writer.font_size_multiple-1.0)*box_size/4*1*0.01667#0.0125
      ymax = writer.calculate_coordinate(yside)
    end
  end
else #when working with multiple ethnicities include legend anyway
  keys = grouplisthash.keys
  glist = grouplisthash[keys[0]]
  glist.groups.length.times do |i|
      y_group_legend_rows << ymax  # + box_size/8
      yside = yside + box_size/4 * 1 * 0.01666 + (writer.font_size_multiple-1.0)*box_size/4*1*0.01667#0.0125
      ymax = writer.calculate_coordinate(yside)
  end
end

if options.groupfile
  thisgrouplist = nil
  if grouplisthash.length == 1
    thisgrouplist = grouplisthash[GroupList.get_default_name]
  else
    thisgrouplist = GroupList.new
    grouplisthash.each_value do |glist|
      glist.grouphash.each_value do |group|
        thisgrouplist.add_group(group)
      end
    end
  end

  # add space for group total plots
  if thisgrouplist.plot_pheno_avg? || thisgrouplist.plot_summary_size? || thisgrouplist.plot_box_plot?
    y_group_total_start = Array.new
    y_group_total_start << ymax + box_size/2
    yside = yside + box_size * thisgrouplist.groups.length * 0.01 + box_size * 0.01
    ymax = writer.calculate_coordinate(yside)
  end

  if vertical_totals and thisgrouplist.plot_summary_size?
    y_group_total_start << ymax + box_size/2
    yside = yside + box_size * thisgrouplist.groups.length * 0.01 + box_size * 0.01
    ymax = writer.calculate_coordinate(yside)
  end
end

if !options.rotate
  xmax += 10
else
  x_start -= 14
end

x_original_start = x_start
x_line_end = 0

or_index = -1
or_min = 0
or_max = 0

# set pixels per inch in the plotwriter
writer.set_pixels_per_coordinate(xside.in, xmax)
writer.maxy = writer.pixels_from_coordinate(ymax)
writer.maxx = writer.pixels_from_coordinate(xmax)

chrom_x_starts = Array.new
chrom_x_ends = Array.new

# Draw features of plot here
rvg = RVG.new(xside.in, yside.in).viewbox(0,0,xmax,ymax) do |canvas|
  canvas.background_fill = 'rgb(253,253,253)'
  writer.canvas = canvas
  writer.box_size = box_size
  # draw for each chromosome
  first_chrom = true

  if vertical_totals
    writer.write_main_title(options.title, 0, y_title_start, xmax,y_gene_start, split_title, max_title_length)
  else
    writer.write_main_title(options.title, 0, y_title_start, xmax,y_gene_start, split_title, max_title_length)
  end

  chromlist.chromarray.each_with_index do |chrnum, chromindex|
    current_chrom = chromlist.chromhash[chrnum]
    # add labels to top of chart
    chrom_x_starts << x_start
    x_new_start= writer.add_labels(current_chrom.chrom_num, current_chrom.snp_list, x_start,
      y_label_start, y_label_end, true, options.rotate,  options.add_snp_locs, options.include_dist_track)
    chrom_x_ends << x_new_start
  # add gene information if necessary
  if gene_rows > 0
    writer.add_gene_rows(current_chrom, x_start, y_gene_start, x_new_start, y_label_start, gene_rows)
  end

  # display which groups were genotyped -- only for single grouplist plot
  if grouplisthash.length == 1 and !options.rotate
    grouplist = grouplisthash[GroupList.get_default_name]
    if !chromlist.results_complete?(grouplist.groups.length)
      contains_meta = false
      meta_group = nil
      meta_index = 0
      grouplist.groups.each_with_index do |group, index|
        if group.name !~ /^meta$/i
          thisindex = index
          if contains_meta
            thisindex = thisindex -1
          end
          writer.group_membership_plot(group.name, group.colorstr, current_chrom.snp_list, x_start, ygroup_boxes_start[thisindex], first_chrom, false,
            options.grayscale)
        else
          meta_group = group
          contains_meta = true
          meta_index = index
        end
      end
      if contains_meta
        writer.group_membership_plot(meta_group.name, meta_group.colorstr, current_chrom.snp_list, x_start, ygroup_boxes_start.last, first_chrom, false,options.grayscale)
      end
    end
  end

  curr_stat_box = 0


  # when single group list from input file
  if grouplisthash.length == 1 and !options.rotate
    grouplist = grouplisthash[GroupList.get_default_name]
  # display p value results

    pvalmin = chromlist.minscore['pvalue']
    if options.pmin
      pvalmin=options.pmin
    end

    if options.plot_pval
      pmaxscore = chromlist.maxscore['pvalue']
      pminscore = chromlist.minscore['pvalue']
      
      writer.draw_pvalue_plot(:jitter=>options.jitter, :no_triangles=>options.no_triangles, :grouplist=>grouplist,
        :snp_list=>current_chrom.snp_list, :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box],
        :stat_max=>pmaxscore, :original_min=>pminscore, :first_plot=>first_chrom,
        :rotate=>options.rotate, :stat_min=>pvalmin, :clean_axis=>options.clean_axes)
      curr_stat_box += 1
    end

    if grouplist.plot_betas? and options.plot_beta
      minbeta = chromlist.minscore['beta'].to_f
      maxbeta = chromlist.maxscore['beta'].to_f
      
      if grouplist.groups.first.betaucicol > -1
        minbeta = chromlist.minscore['betalci'].to_f
        maxbeta = chromlist.maxscore['betauci'].to_f
      end  
        
        
      if options.clean_axes
        increment, minbeta, maxbeta = writer.calculate_increments_include_zero(minbeta, maxbeta)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>maxbeta, :stat_min=>minbeta,
        :data_key=>'beta', :plot_labels=>first_chrom, :title=>options.effect_name, :precision=>2,
        :rotate=>options.rotate, :lci_key=>'betalci', :uci_key=>'betauci')
      curr_stat_box+=1
    end

    # plot the odds ratio results
    if grouplist.plot_oddsratio?

      if !options.ormin_zero
        ormin = chromlist.minscore['or']
      else
        ormin = 0
      end
      ormax = chromlist.maxscore['or']

      if options.clean_axes
        increment, or_min, or_max = writer.calculate_increments(ormin.to_f, ormax.to_f)
      else
        or_min = ormin.to_f
        or_max = ormax.to_f
      end
      writer.draw_basic_or_plot(options.jitter, grouplist, current_chrom.snp_list, x_start, ystart_stat_boxes[curr_stat_box], or_max, or_min, 'or', first_chrom, 'Odds Ratio',2, options.large_or)
      or_index = curr_stat_box
      curr_stat_box += 1
    end

    if grouplist.plot_maf?
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>1.0, :stat_min=>0.0,
        :data_key=>'maf', :plot_labels=>first_chrom, :title=>grouplist.mafcoltitle, :precision=>2)
        
      curr_stat_box+=1
    end

    if options.cafcasecontrol
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>1.0, :stat_min=>0.0,
        :data_key=>'cafcases', :plot_labels=>first_chrom, :title=>'CAF Cases', :precision=>1)
      curr_stat_box+=1
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>1.0, :stat_min=>0.0,
        :data_key=>'cafcontrols', :plot_labels=>first_chrom, :title=>'CAF Controls', :precision=>1)
      curr_stat_box+=1
    end

    if options.circlecafcasecontrol
      writer.draw_circles_plot(grouplist, current_chrom.snp_list, x_start, ystart_stat_boxes[curr_stat_box],
        1, 0, 'cafcases', 'cafcontrols', first_chrom, 'CAF Cases/Controls',1)
      curr_stat_box+=1
    end

    if grouplist.plot_sample_sizes?
      nmin = chromlist.minscore['N'].to_f
      nmax = chromlist.maxscore['N'].to_f
      if options.clean_axes
        increment, nmin, nmax = writer.calculate_increments(nmin, nmax)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>nmax, :stat_min=>nmin,
        :data_key=>'N', :plot_labels=>first_chrom, :title=>'Sample Size', :precision=>0)
      curr_stat_box+=1
    end

    if options.plot_power
      powmin = 0
      powmax = 100
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>powmax,
        :stat_min=>powmin, :data_key=>'power', :plot_labels=>first_chrom, :title=>'Power', :precision=>0)     
      curr_stat_box+=1
    end

    if options.casecontrol
      caseconmax = chromlist.maxscore['cases'].to_f
      caseconmin = chromlist.minscore['cases'].to_f
      if options.clean_axes
        increment, caseconmin, caseconmax = writer.calculate_increments(caseconmin, caseconmax)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>caseconmax,
        :stat_min=>caseconmin, :data_key=>'cases', :plot_labels=>first_chrom, :title=>'Cases', :precision=>0)
      curr_stat_box+=1
      caseconmax = chromlist.maxscore['controls'].to_f
      caseconmin = chromlist.minscore['controls'].to_f
      if options.clean_axes
        increment, caseconmin, caseconmax = writer.calculate_increments(caseconmin, caseconmax)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>caseconmax,
        :stat_min=>caseconmin, :data_key=>'controls', :plot_labels=>first_chrom, :title=>'Controls', :precision=>0)
      curr_stat_box+=1
    end

    if options.circlecasecon
      if chromlist.maxscore['cases'].to_f > chromlist.maxscore['controls'].to_f
        caseconmax = chromlist.maxscore['cases'].to_f
      else
        caseconmax = chromlist.maxscore['controls'].to_f
      end
      if chromlist.minscore['cases'].to_f < chromlist.minscore['controls'].to_f
        caseconmin = chromlist.minscore['cases'].to_f
      else
        caseconmin = chromlist.minscore['controls'].to_f
      end
      if options.clean_axes
        increment, caseconmin, caseconmax = writer.calculate_increments(caseconmin, caseconmax)
      end
      writer.draw_circles_plot(grouplist, current_chrom.snp_list, x_start, ystart_stat_boxes[curr_stat_box],
        caseconmax, caseconmin, 'cases', 'controls', first_chrom, 'Cases/Controls',0)
      plot_labels = (chromindex == chromlist.chromarray.length-1)
      if plot_labels
        writer.add_circle_plot_legend('Cases', 'Controls', options.rotate, x_new_start,
          ystart_stat_boxes[curr_stat_box], circle_legend_below)
      end
      curr_stat_box+=1
    end

    options.additional_columns.each do |column|
      cmin = chromlist.minscore[column].to_f
      cmax = chromlist.maxscore[column].to_f
      if options.clean_axes
        increment, cmin, cmax = writer.calculate_increments(cmin, cmax)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>cmax, :stat_min=>cmin,
        :data_key=>column, :plot_labels=>first_chrom, :title=>column, :precision=>2)
      curr_stat_box+=1
    end  
      
      
    curr_half_box = 0
    if grouplist.plot_study? and options.plot_studynum
      studymax = chromlist.maxscore['study'].to_f
      studymin = chromlist.minscore['study'].to_f
      if options.clean_axes
        increment, studymax, studymin = writer.calculate_increments(studymax, studymin)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_half_boxes[curr_half_box], :stat_max=>10,
        :stat_min=>0, :data_key=>'study', :plot_labels=>first_chrom, :title=>'Study #',
        :precision=>0, :rotate=>options.rotate, :size_mult=>0.5)
    end

  elsif options.rotate # for rotated plots for odds ratio
    plot_labels = (chromindex == chromlist.chromarray.length-1)

    grouplist = grouplisthash[GroupList.get_default_name]

    pvalmin = chromlist.minscore['pvalue']
    if options.pmin
      pvalmin=options.pmin
    end

    curr_stat_box = 0

    # plot p values if present
    if grouplist.plot_pvals? and options.plot_pval
      writer.draw_pvalue_plot(:jitter=>options.jitter, :no_triangles=>options.no_triangles, :grouplist=>grouplist,
        :snp_list=>current_chrom.snp_list, :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box],
        :stat_max=>chromlist.maxscore['pvalue'], :stat_min=>pvalmin, :original_min=>chromlist.minscore['pvalue'],
        :first_plot=>plot_labels, :rotate=>options.rotate, :clean_axis=>options.clean_axes)
      curr_stat_box +=1
    end

    # plot betas if present and requested
    if grouplist.plot_betas? and options.plot_beta
      minbeta = chromlist.minscore['beta'].to_f
      maxbeta = chromlist.maxscore['beta'].to_f
      if options.clean_axes
        increment, minbeta, maxbeta = writer.calculate_increments_include_zero(minbeta, maxbeta)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>maxbeta,
        :stat_min=>minbeta, :data_key=>'beta', :plot_labels=>plot_labels, :title=>options.effect_name,
        :precision=>2, :rotate=>options.rotate, :lci_key=>'betalci', :uci_key=>'betauci')      
      curr_stat_box+=1
    end

    if !options.ormin_zero
      if chromlist.minscore['lci'].to_f != 10
        ormin = chromlist.minscore['lci']
      else
        ormin = chromlist.minscore['or']
      end
    else
      ormin = 0
    end

    if chromlist.maxscore['uci'].to_f > 0
      ormax = chromlist.maxscore['uci']
    else
      ormax = chromlist.maxscore['or']
    end

    # plot each groups odds ratio values
    grouplist.groups.each do |group|
      if options.clean_axes
        increment, newormin, newormax = writer.calculate_increments(ormin.to_f, ormax.to_f)
      else
        newormin = ormin.to_f
        newormax = ormax.to_f
      end
      writer.draw_or_plot(group, current_chrom.snp_list, x_start, ystart_stat_boxes[curr_stat_box],newormin, newormax, plot_labels, options.large_or)
      curr_stat_box += 1
    end

    if grouplist.plot_maf?
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>1.0,
        :stat_min=>0, :data_key=>'maf', :plot_labels=>plot_labels, :title=>grouplist.mafcoltitle,
        :precision=>2, :rotate=>options.rotate)      
      curr_stat_box+=1
    end

    if options.plot_power
      powmin = 0
      powmax = 100
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>powmax,
        :stat_min=>powmin, :data_key=>'power', :plot_labels=>plot_labels, :title=>'Power',
        :precision=>0, :rotate=>options.rotate)        
      curr_stat_box+=1
    end

    if options.cafcasecontrol
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>1,
        :stat_min=>0, :data_key=>'cafcases', :plot_labels=>plot_labels, :title=>'CAF Cases',
        :precision=>1, :rotate=>options.rotate)  
      curr_stat_box+=1
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>1,
        :stat_min=>0, :data_key=>'cafcontrols', :plot_labels=>plot_labels, :title=>'CAF Controls',
        :precision=>1, :rotate=>options.rotate) 
      curr_stat_box+=1
    end

    if options.circlecafcasecontrol
      if first_chrom
        writer.add_circle_plot_legend('Cases', 'Controls', options.rotate, x_start, ystart_stat_boxes[curr_stat_box])
      end
      writer.draw_circles_plot(grouplist, current_chrom.snp_list, x_start, ystart_stat_boxes[curr_stat_box],
        1.0, 0.0, 'cafcases', 'cafcontrols', plot_labels, 'CAF Cases/Controls',1, options.rotate)
      curr_stat_box+=1
    end

    if options.casecontrol
      caseconmax = chromlist.maxscore['cases'].to_f
      caseconmin = chromlist.minscore['cases'].to_f
      if options.clean_axes
        increment, caseconmin, caseconmax = writer.calculate_increments(caseconmin, caseconmax)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>caseconmax,
        :stat_min=>caseconmin, :data_key=>'cases', :plot_labels=>plot_labels, :title=>'Cases',
        :precision=>0, :rotate=>options.rotate) 
      curr_stat_box+=1
      caseconmax = chromlist.maxscore['controls'].to_f
      caseconmin = chromlist.minscore['controls'].to_f
      if options.clean_axes
        increment, caseconmin, caseconmax = writer.calculate_increments(caseconmin, caseconmax)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_stat_boxes[curr_stat_box], :stat_max=>caseconmax,
        :stat_min=>caseconmin, :data_key=>'controls', :plot_labels=>plot_labels, :title=>'Controls',
        :precision=>0, :rotate=>options.rotate) 
      curr_stat_box+=1
    end

    if options.circlecasecon
      if first_chrom
        writer.add_circle_plot_legend('Cases', 'Controls', options.rotate, x_start, ystart_stat_boxes[curr_stat_box])
      end
      if chromlist.maxscore['cases'].to_f > chromlist.maxscore['controls'].to_f
        caseconmax = chromlist.maxscore['cases'].to_f
      else
        caseconmax = chromlist.maxscore['controls'].to_f
      end
      if chromlist.minscore['cases'].to_f < chromlist.minscore['controls'].to_f
        caseconmin = chromlist.minscore['cases'].to_f
      else
        caseconmin = chromlist.minscore['controls'].to_f
      end
      if options.clean_axes
        increment, caseconmin, caseconmax = writer.calculate_increments(caseconmin, caseconmax)
      end
      writer.draw_circles_plot(grouplist, current_chrom.snp_list, x_start, ystart_stat_boxes[curr_stat_box],
        caseconmax, caseconmin, 'cases', 'controls', plot_labels, 'Cases/Controls',0, options.rotate)
      curr_stat_box+=1
    end

    curr_half_box = 0
    if grouplist.plot_study? and options.plot_studynum
      studymax = chromlist.maxscore['study'].to_f
      studymin = chromlist.minscore['study'].to_f
      if options.clean_axes
        increment, studymax, studymin = writer.calculate_increments(studymax, studymin)
      end
      writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
        :x_start=>x_start, :y_start=>ystart_half_boxes[curr_half_box], :stat_max=>10,
        :stat_min=>0, :data_key=>'study', :plot_labels=>first_chrom, :title=>'Study #',
        :precision=>0, :rotate=>options.rotate, :size_mult=>0.5)     
    end

  else # multiple grouplists with each one being a different ethnicity

    pvalmin = chromlist.minscore['pvalue']
    if options.pmin
      pvalmin=options.pmin
    end
    grouplistkeys = grouplisthash.keys
    grouplistkeys.length.times do |i|
      writer.draw_pvalue_plot(:jitter=>options.jitter, :no_triangles=>options.no_triangles, :grouplist=>grouplisthash[grouplistkeys[i]],
        :snp_list=>current_chrom.snp_list, :x_start=>x_start, :y_start=>ystart_stat_boxes[i],
        :stat_max=>chromlist.maxscore['pvalue'], :stat_min=>pvalmin, :original_min=>chromlist.minscore['pvalue'],
        :first_plot=>first_chrom, :prefix_title=>grouplistkeys[i] + "\n", :rotate=>options.rotate,
        :clean_axis=>options.clean_axes)
    end

    # draw beta values on plots
    grouplistkeys.each_with_index do |grouplistname, i|
      i = grouplistkeys.length + i
      grouplist = grouplisthash[grouplistname]
      if grouplist.plot_betas? and options.plot_beta
        minbeta = chromlist.minscore['beta'].to_f
        maxbeta = chromlist.maxscore['beta'].to_f
        if options.clean_axes
          increment, minbeta, maxbeta = writer.calculate_increments_include_zero(minbeta, maxbeta)
        end
        writer.draw_basic_plot(:jitter=>options.jitter, :grouplist=>grouplist, :snp_list=>current_chrom.snp_list,
          :x_start=>x_start, :y_start=>ystart_stat_boxes[i], :stat_max=>maxbeta,
          :stat_min=>minbeta, :data_key=>'beta', :plot_labels=>first_chrom, :title=>options.effect_name,
          :precision=>2, :rotate=>options.rotate, :prefix_title=>grouplistname + "\n", :lci_key=>'betalci', :uci_key=>'betauci') 
        curr_stat_box+=1
      end
    end

  end

  # draw d' grid when more than one snp to draw
  if current_chrom.snp_list.get_num_included_snps > 0
    if(options.dprime)
     writer.first_grid = first_chrom
     writer.draw_grid(current_chrom.snp_list, x_start, y_dprime_start, true)
    end

    # draw r-squared grid
    if(options.rsquared)
      writer.first_grid = first_chrom
      writer.draw_grid(current_chrom.snp_list, x_start, y_rsquared_start, false)
    end
  end

  first_chrom = false
  x_start=x_new_start
  x_line_end = x_new_start
  end

end

# add manhattan plot if selected
if options.manhattanfile
  man_plot = Manhattan::ManhattanPlotter.new
  man_plot.plot(writer.canvas, manhattan_results, x_original_start, y_manhattan_start,
      x_start, y_manhattan_end-box_size, box_size/10, options.manhattanthresh)

  chromlist.chromarray.each do |chrnum|
    current_chrom = chromlist.chromhash[chrnum]
    locs = Array.new
    locs << current_chrom.snp_list.get_min
    locs << current_chrom.snp_list.get_max

    manhattan_x_pos = man_plot.manhattan_x_locations(manhattan_results, chrnum, locs)
    writer.canvas.g do |man_line|
      man_line.styles(:fill=>'black', :stroke_width=>'2')
      man_line.line(x_original_start+manhattan_x_pos[0], y_manhattan_end, x_original_start+manhattan_x_pos[0], y_manhattan_end-box_size/2)
      man_line.line(x_original_start+manhattan_x_pos[1], y_manhattan_end, x_original_start+manhattan_x_pos[1], y_manhattan_end-box_size/2)
      man_line.line(x_original_start+manhattan_x_pos[0], y_manhattan_end, chrom_x_starts[0]+box_size, y_label_start)
      man_line.line(x_original_start+manhattan_x_pos[0], y_manhattan_end, chrom_x_ends[0]-box_size, y_label_start)
    end

  end

end

if grouplisthash.length == 1 and !options.rotate
  grouplist = grouplisthash[GroupList.get_default_name]
  # draw p value threshold across pvalue plot
  if options.p_thresh > 0
    writer.draw_red_line(x_original_start, ystart_stat_boxes[0], x_line_end, chromlist.maxscore['pvalue'],
      chromlist.minscore['pvalue'], options.p_thresh, options.rotate)
  end

  # draw dotted line across beta
  if grouplist.plot_betas? and options.plot_beta
    minbeta = chromlist.minscore['beta'].to_f
    maxbeta = chromlist.maxscore['beta'].to_f
    if options.clean_axes
      increment, minbeta, maxbeta = writer.calculate_increments_include_zero(minbeta, maxbeta)
    end
    if 0 > minbeta and 0 < maxbeta
      writer.draw_dashed(x_original_start, ystart_stat_boxes[1], x_line_end, 0, maxbeta, minbeta)
    end
  end

  # draw dotted line across odds ratio at 1.0
  if grouplist.plot_oddsratio?
    writer.draw_dashed(x_original_start, ystart_stat_boxes[or_index], x_line_end, 1.0, or_max, or_min)
  end

  # draw lines across value plots
  total_stat_boxes.times do |i|
    writer.draw_plot_boundaries(x_original_start, x_line_end, ystart_stat_boxes[i])
  end

  # draw lines around half-size plots
  ystart_half_boxes.length.times do |i|
    writer.draw_plot_boundaries(x_original_start, x_line_end, ystart_half_boxes[i], 0.5)
  end

  # draw group legend when needed
  if y_group_legend_rows.length > 0
    grouplist.groups.each_with_index do |group, index|
      writer.group_membership_plot(group.name, group.colorstr,  chromlist.chromhash[chromlist.chromarray.first].snp_list, x_original_start, y_group_legend_rows[index], true, true, options.grayscale)
    end
  end

  # draw group summary plots
  if options.groupfile
    x_number_plot = x_original_start
    if grouplist.plot_pheno_avg?
      x_number_plot = writer.group_summary_plot(grouplist, x_original_start, y_group_total_start[0], options.phenotitle, 'pheno_avg',1)
    end
    if grouplist.plot_box_plot?
      x_number_plot = writer.group_box_plot(grouplist, x_number_plot, y_group_total_start[0], 1)
    end
    if grouplist.plot_summary_size?
      if vertical_totals
        writer.group_summary_plot(grouplist, x_original_start, y_group_total_start[1], 'N', 'num',50)
      else
        writer.group_summary_plot(grouplist, x_number_plot, y_group_total_start[0], 'N', 'num',50)
      end
    end
  end
# for rotation and display of forest plots
elsif options.rotate
  #draw lines across value plots
  total_stat_boxes.times do |i|
    writer.draw_plot_boundaries(x_original_start, x_line_end, ystart_stat_boxes[i])
  end

  # draw lines around half-size plots
  ystart_half_boxes.length.times do |i|
    writer.draw_plot_boundaries(x_original_start, x_line_end, ystart_half_boxes[i], 0.5)
  end

  if options.p_thresh > 0
    writer.draw_red_line(x_original_start, ystart_stat_boxes[0], x_line_end, chromlist.maxscore['pvalue'],
      chromlist.minscore['pvalue'], options.p_thresh, options.rotate)
  end

  # draw dotted line across beta
  if grouplist.plot_betas? and options.plot_beta
    minbeta = chromlist.minscore['beta'].to_f
    maxbeta = chromlist.maxscore['beta'].to_f
    if options.clean_axes
       increment, minbeta, maxbeta = writer.calculate_increments_include_zero(minbeta, maxbeta)
    end
    if 0 > minbeta and 0 < maxbeta
      writer.draw_dashed(x_original_start, ystart_stat_boxes[1], x_line_end, 0, maxbeta, minbeta)
    end
  end

  # add legend at bottom when requested
  if options.forest_legend
    writer.add_forest_legend(x_forest_legend_lines, yforest_legend_start, grouplisthash[GroupList.get_default_name], 0, 0)
  end

else # multiple grouplists for ethnicity so different lines need to be drawn
  grouplistkeys = grouplisthash.keys
  grouplistkeys.length.times do |i|
    if options.p_thresh > 0
      writer.draw_red_line(x_original_start, ystart_stat_boxes[i], x_line_end, chromlist.maxscore['pvalue'],
        chromlist.minscore['pvalue'], options.p_thresh, options.rotate)
    end
    writer.draw_plot_boundaries(x_original_start, x_line_end, ystart_stat_boxes[i])
  end


  grouplistkeys.each_with_index do |grouplistname, i|
    i = grouplistkeys.length + i
    grouplist = grouplisthash[grouplistname]
    # draw dotted line across beta
    if grouplist.plot_betas? and options.plot_beta
      minbeta = chromlist.minscore['beta'].to_f
      maxbeta = chromlist.maxscore['beta'].to_f
      if options.clean_axes
        increment, minbeta, maxbeta = writer.calculate_increments_include_zero(minbeta, maxbeta)
      end
      if 0 > minbeta and 0 < maxbeta
        writer.draw_dashed(x_original_start, ystart_stat_boxes[i], x_line_end, 0, maxbeta, minbeta)
      end
      writer.draw_plot_boundaries(x_original_start, x_line_end, ystart_stat_boxes[i])
    end
  end

  if options.groupfile
    combinedgrouplist = GroupList.new
    grouplisthash.each_value do |glist|
      glist.grouphash.each_value do |group|
        combinedgrouplist.add_group(group)
      end
    end
    x_number_plot = x_original_start
    if combinedgrouplist.plot_pheno_avg?
      x_number_plot = writer.group_summary_plot(combinedgrouplist, x_original_start, y_group_total_start[0], options.phenotitle, 'pheno_avg',1)
    end
    if grouplist.plot_box_plot?
      x_number_plot = writer.group_box_plot(grouplist, x_number_plot, y_group_total_start[0], 1)
    end

    if combinedgrouplist.plot_summary_size?
      if vertical_totals
        writer.group_summary_plot(combinedgrouplist, x_original_start, y_group_total_start[1], 'N', 'num',50)
      else
        writer.group_summary_plot(combinedgrouplist, x_number_plot, y_group_total_start[0], 'N', 'num',50)
      end
    end
  end

  # draw group legend when needed
  if y_group_legend_rows.length > 0
    keys = grouplisthash.keys
    glist = grouplisthash[keys[0]]

    glist.groups.each_with_index do |group, index|
      writer.group_membership_plot(group.name, group.colorstr,  chromlist.chromhash[chromlist.chromarray.first].snp_list, x_original_start, y_group_legend_rows[index], true, true, options.grayscale)
    end
  end

end


# produce output file
outfile = options.out_name + '.' + options.imageformat
print "\n\tDrawing #{outfile}..."
STDOUT.flush
img = rvg.draw
#if options.rotate
img.rotate!(-90) if options.rotate
#end
img.write(outfile)

if options.htmlfile
  if options.rotate
    width = yside.in
    height = xside.in
  else
    width = xside.in
    height = yside.in
  end

  if width > height
    max_dim = width
  else
    max_dim = height
  end

  imgname = outfile
  if max_dim > Maximum_html_image_x
    fraction =  Maximum_html_image_x.to_f / max_dim
  else
    fraction = 1.0
  end

  img.scale!(fraction)
  imgname = options.out_name + '.' + 'small' + '.' + options.imageformat
  img.write(imgname)
  DBSNPhtml.write_html(imgname, width*fraction, height*fraction, writer.map_pos, options.out_name + '.html', options.rotate)
end

print " Created #{outfile}\n\n"
