#!/usr/bin/ruby
#
require 'bio'
include Bio

class GetPosHash # GetPosHash

  def initialize( speAgid, speA_prot_f, blastDBspeA, blastDBspeB )
    @speApeph         = {}
    @cutting_region_h = {}
    @speAgid          = speAgid
    @blastDBspeA = blastDBspeA
    @blastDBspeB = blastDBspeB
    check_makeblastdb
    create_protein_h( speA_prot_f ) # speAh, speBh
  end

  attr_reader :cutting_region_h

  def check_makeblastdb
 
    unless FileTest.exist?( "#{@blastDBspeA}.nhr" )
      ` /home/wadaken/bin/makeblastdb -in #{@blastDBspeA} -dbtype nucl `
#      system(" makeblastdb -db #{@blastDBspeA} -dbtype nucl ")
    end
    unless FileTest.exist?( "#{@blastDBspeB}.nhr" )
      ` /home/wadaken/bin/makeblastdb -in #{@blastDBspeB} -dbtype nucl `
#      system(" makeblastdb -db #{@blastDBspeB} -dbtype nucl ")
    end
    STDERR.puts "makeblasdb check DONE"
 
  end

  def create_protein_h( speA_protein_f )

    speAf = FlatFile.new( FastaFormat, open( speA_protein_f ) ) # protein seqs
    speAf.each{ |e| @speApeph[ e.definition ] = e.aaseq }
    STDERR.puts "Method create_protein_h Done. (step 2) "

  end # create_protein_h in GetPosHash 

  def pos_for_target
    get_pos_by_blast( @blastDBspeB )
  end

  def pos_for_query_self
    get_pos_by_blast( @blastDBspeA )
  end
    
  def get_pos_by_blast( blastDB )
    # example of qry_hash is @speAh or @speBh

    qry_tmp_filename = "tmp.#{Time.now.strftime("%Y%m%d%R%N")}"
    qry_tmp = File.new( qry_tmp_filename, "w+" )
    qry_tmp.print @speApeph[ @speAgid ].to_fasta( @speAgid, 60 )
    qry_lengt = @speApeph[ @speAgid ].length
    qry_tmp.close
    qry_tbn_out = ` cat #{qry_tmp_filename} | tblastn -query - -db #{blastDB} -max_target_seqs 1 -outfmt 6 `
    qry_array = qry_tbn_out.split("\t")
    qry_stt   = 0
    qry_end   = 0
    tgt_str   = ""
    if qry_array != []
      tgt_sca = qry_array[1]
      qry_stt = qry_array[6].to_i
      qry_end = qry_array[7].to_i
      tgt_stt = qry_array[8].to_i
      tgt_end = qry_array[9].to_i
      tgt_str = "+" if tgt_end - tgt_stt > 0
      tgt_str = "-" if tgt_stt - tgt_end > 0
      tgt_stt = qry_array[9].to_i if tgt_str == "-"
      tgt_end = qry_array[8].to_i if tgt_str == "-"
    end
    File.delete( qry_tmp )
    return [ @speAgid, tgt_sca, tgt_stt, tgt_end, tgt_str ]
    
  end # get_pos_by_blast

  def sampletest
    "good result"
  end

end # GetPosHash

if __FILE__ == $0

  speAgid        = "SPU_025302"
#  speA_protein_f = ARGV.shift
  speA_protein_f = "./sample_data/spur_protein.fas"
#  blastDBspeA    = ARGV.shift # scaffold nucleotide sequences
  blastDBspeA    = "./sample_data/spur_scaffold.fas"
#  blastDBspeB    = ARGV.shift # scaffold nucleotide sequences
  blastDBspeB    = "./sample_data/hpul_scaffold.fas"

  gph = GetPosHash.new( speAgid, speA_protein_f, blastDBspeA, blastDBspeB )
  p gph.pos_for_query_self
  p gph.pos_for_target

end


