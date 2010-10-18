####
# Analysis of the Mikkelsen '07 Genome wide histone modification data
# in ES, NPC, MEF cells
####

require 'rubygems'
require 'fileutils'
require 'rake/clean'

directory "publication"
CLOBBER.include('publication')

directory "scripts"
CLOBBER.include('scripts')
CLEAN.include('scripts')

directory "lib"
CLOBBER.include('lib')
CLEAN.include('lib')

directory 'results'
CLOBBER.include('results')

directory 'results/WindowIntervals'
directory 'results/HMMIntervals'


####
# Fetch the publication and supplemental data
###
window_fl = FileList['publication/WindowIntervals/*.txt']
hmm_fl = FileList['publication/HMMIntervals/*.txt']

task :fetch_publication => ['publication'] do 
  #TODO - actually get the data :)

  #tidy up a bit
  sh "mv publication/WindowIntervals/readme.txt publication/WindowIntervals/README"
  sh "mv publication/HMMIntervals/readme.txt publication/HMMIntervals/README"

  #and reset notes of what we have
  window_fl = FileList['publication/WindowIntervals/*.txt']
  hmm_fl = FileList['publication/HMMIntervals/*.txt']
end


##scripts

file 'scripts/window2bed.R' => ['scripts'] do
  sh "wget -O scripts/window2bed.R  'http://github.com/cassj/my_bioinfo_scripts/raw/master/Mikkelsen07/window2bed.R'"
end

file 'scripts/hmm2bed.R' => ['scripts'] do
  sh "wget -O scripts/hmm2bed.R  'http://github.com/cassj/my_bioinfo_scripts/raw/master/Mikkelsen07/hmm2bed.R'"
end

file 'scripts/liftOver.R' => ['scripts'] do
    sh "wget -O scripts/liftOver.R  'http://github.com/cassj/my_bioinfo_scripts/raw/master/liftOver.R'"
end

file 'scripts/qw.R' => ['scripts'] do
    sh "wget -O scripts/qw.R  'http://github.com/cassj/my_bioinfo_scripts/raw/master/qw.R'"
end


file 'scripts/mm8tomm9.R' => ['scripts', 'scripts/liftOver.R', 'scripts/qw.R'] do
  sh "wget -O scripts/mm8tomm9.R  'http://github.com/cassj/my_bioinfo_scripts/raw/master/Mikkelsen07/mm8tomm9.R'"
end


file 'scripts/mm9toiranges.R' => ['scripts', 'scripts/qw.R'] do
  sh "wget -O scripts/mm9toiranges.R  'http://github.com/cassj/my_bioinfo_scripts/raw/master/Mikkelsen07/mm9toiranges.R'"
end

desc "Grab all the scripts needed"
task :fetch_scripts => ['scripts/window2bed.R', 'scripts/hmm2bed.R', 'scripts/liftOver.R', 'scripts/qw.R', 'scripts/mm8tomm9.R', 'scripts/mm9toiranges.R'] do end

#file rules for the raw data files
rule(/^publication\/WindowIntervals\/.*\.txt$/ => ['publication']) do end
rule(/^publication\/HMMIntervals\/.*\.txt$/ => ['publication']) do end




# supporting  libs etc

file 'lib/mm8ToMm9.over.chain' => ['lib'] do
  sh "wget -O lib/mm8ToMm9.over.chain.gz  'http://hgdownload.cse.ucsc.edu/goldenPath/mm8/liftOver/mm8ToMm9.over.chain.gz'"
  sh "gunzip -c lib/mm8ToMm9.over.chain.gz > lib/mm8ToMm9.over.chain"
end


desc 'fetch all supporting libs' 
task :fetch_libs => ['lib/mm8ToMm9.over.chain'] do end

###
# Generate mm8 BED files from Window Data
###



#eg rake results/WindowIntervals/ES.K27.mm8.bed
rule(/^results\/WindowIntervals\/.*\.mm8\.bed$/ => [proc {|tn| tn.sub(/results/,'publication').sub(/\.mm8\.bed/,'.txt')}, 'results/WindowIntervals', 'scripts/window2bed.R'
]) do |t|
  sh %Q(R --vanilla --args filename=\\"#{t.source}\\" < scripts/window2bed.R)
end

#make a list of the bed files we can generate from the windowinterval raw data we have
window_mm8_bed = window_fl.gsub(/publication/,'results').gsub(/\.txt/,'.mm8.bed')

#and generate them all
desc 'make window-interval mm8 bed files'
task :mm8_bed_windowdata => ['scripts/window2bed.R', 'results/WindowIntervals'] | window_mm8_bed






###
# Generate mm8 BED files from HMM Data
###

#eg rake results/HMMIntervals/HMM_ES_K27.mm8.bed
rule(/^results\/HMMIntervals\/.*\.mm8\.bed$/ => [proc {|tn| tn.sub(/results/,'publication').sub(/\.mm8\.bed/,'.txt')}, 'results/HMMIntervals', 'scripts/hmm2bed.R'
]) do |t|
  sh %Q(R --vanilla --args filename=\\"#{t.source}\\" < scripts/hmm2bed.R)
end

#make a list of the bed files we can generate from the windowinterval raw data we have
hmm_mm8_bed = hmm_fl.gsub(/publication/,'results').gsub(/\.txt/,'.mm8.bed')

#and generate them all
desc 'make hmm-interval mm8 bed files'
task :mm8_bed_hmmdata => ['scripts/hmm2bed.R', 'results/HMMIntervals'] | hmm_mm8_bed




###
# Generate mm9 BED files from mm8 (both Window and HMM)
###

#eg rake results/WindowIntervals/ES.K27.mm9.bed
rule(/^results\/.*\.mm9\.bed$/ => [proc {|tn| tn.sub(/mm9/,'mm8')}, 'scripts/mm8tomm9.R'
]) do |t|
  sh %Q(R --vanilla --args filename=\\"#{t.source}\\" < scripts/mm8tomm9.R)

end

#make a list of the bed files we can generate from the windowinterval raw data we have
window_mm9_bed = window_mm8_bed.gsub(/mm8/,'mm9')

#make a list of the bed files we can generate from the hmminterval raw data we have
hmm_mm9_bed = hmm_mm8_bed.gsub(/mm8/,'mm9')

#and generate them all
desc 'make all mm9 bed files'
task :mm9_bed => ['scripts/mm8tomm9.R'] | window_mm9_bed | hmm_mm9_bed



###
# Generate IRanges RangedData from mm9 data
###



#### NOTE - if we're assuming the original data are BED format, then they should be 0-based, half open,
#in which case this script should Start++ to make the IRanges 1-based fully closed.

#eg rake results/WindowIntervals/ES.K27.RangedData.R
rule(/^results\/.*\.RangedData\.R$/ => [proc {|tn| tn.sub(/\.RangedData\.R/,'.mm9.bed')}, 'scripts/mm9toiranges.R'
]) do |t|
  sh %Q(R --vanilla --args filename=\\"#{t.source}\\" < scripts/mm9toiranges.R)
end

#make a list of the bed files we can generate from the windowinterval raw data we have
iranges = window_mm9_bed.gsub(/\.mm9\.bed/,'.RangedData.R') |  hmm_mm9_bed.gsub(/\.mm9\.bed/,'.RangedData.R')

#and generate them all
desc 'make all iranges files'
task :iranges => iranges


