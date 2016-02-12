#!/usr/bin/ruby
#
require 'bio'
include Bio

class GetPosHash # GetPosHash

  def initialize( listf, speA_prot_f, speB_prot_f, blstDBspeA, blstDBspeB )
    @h = {}
    @cutting_region_h = {}
    @speAh = {}; 
    @speBh = {};
    create_h( listf )                                  # h
    create_protein_h( speA_prot_f, speB_prot_f ) # speAh, speBh
    @BlastDBspeA = blstDBspeA
    @BlastDBspeB = blstDBspeB
    get_pos_for_each_species
  end

  attr_reader :cutting_region_h, :h

  def create_h( listf )

    @h = {}
    ortho_id = ""
    speA_gid = ""
    speB_gid = ""
    a = []
    listf.each do |x|
      next if x =~ /^\#/
      a = x.chomp.split("\t")
      ortho_id = "#{a[0]}_#{a[1]}"
      speA_gid = a[2]
      speB_gid = a[3]
      speB_gid = speB_gid.slice(/(\S+)gn/, 1) if speB_gid =~ /gn/
      @h[ ortho_id ] = [ speA_gid, speB_gid ]
      @cutting_region_h[ ortho_id ] = [] if @cutting_region_h[ ortho_id ] == nil
    end
    STDERR.puts "Method create_h Done. (step 1) "

  end

  def create_protein_h( speA_protein_f, speB_protein_f )

    speAf = FlatFile.new( FastaFormat, open( speA_protein_f ) ) # protein seqs
    speBf = FlatFile.new( FastaFormat, open( speB_protein_f ) ) # protein seqs
    speAf.each{ |e| @speAh[ e.definition ] = e.aaseq }
    speBf.each{ |e| @speBh[ e.definition ] = e.aaseq }
    STDERR.puts "Method create_protein_h Done. (step 2) "
  end # create_protein_h in GetPosHash 
    
  def get_pos_by_blast( qry_hash, qry_gid, qry_length, blastdb ) 
    # example of qry_hash is @speAh or @speBh

    qry_tmp_filename = "tmp.#{Time.now.strftime("%Y%m%d%R%N")}"
    qry_tmp = File.new( qry_tmp_filename, "w+" )
    qry_tmp.print qry_hash[ qry_gid ].to_fasta( qry_gid, 60 )
    qry_length  = qry_hash[ qry_gid ].length
    qry_tmp.close
    qry_tbn_out = ` cat #{qry_tmp_filename} | tblastn -query - -db #{blastdb} -max_target_seqs 1 -outfmt 6 `
    qry_array = qry_tbn_out.split("\t")
    qry_stt   = 0
    qry_end   = 0
    qry_str   = ""
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
      if    tgt_str == "+"
        reg_stt = tgt_stt - ( qry_stt * 3 )
        reg_stt = 1 if reg_stt < 1
        reg_end = tgt_end + ( ( qry_length - qry_end ) * 3 )
      elsif tgt_str == "-"
        reg_stt = tgt_stt - ( qry_stt * 3 )
        reg_stt = 1 if reg_stt < 1
        reg_end = tgt_end + ( ( qry_length - qry_end ) * 3 )
      end
    end
    File.delete( qry_tmp )
    return [ qry_gid, tgt_sca, reg_stt, reg_end, tgt_str ]
    
  end # get_pos_by_blast

  def get_pos_for_each_species
    @h.each_key do |ortho_id|
      speA_gid, speB_gid = @h[ ortho_id ]
      speA_gid_length = @speAh[ speA_gid ].length
      speB_gid_length = @speBh[ speB_gid ].length
      @cutting_region_h[ ortho_id ][0] = get_pos_by_blast( @speAh, speA_gid, speA_gid_length, @BlastDBspeB ) 
      @cutting_region_h[ ortho_id ][1] = get_pos_by_blast( @speBh, speB_gid, speB_gid_length, @BlastDBspeA ) 
    end
  end

end # GetPosHash

if __FILE__ == $0

  listf = open( ARGV.shift )
  speA_protein_f = ARGV.shift
  speB_protein_f = ARGV.shift
  BlastDBspeA = ARGV.shift # scaffold nucleotide sequences
  BlastDBspeB = ARGV.shift # scaffold nucleotide sequences

  gph = GetPosHash.new( listf, 
                        speA_protein_f, speB_protein_f, 
                        BlastDBspeA, BlastDBspeB )
  gph.cutting_region_h.each_key do |ortho_id|
    speA_gid, speA_tgt_sca, speA_reg_stt, speA_reg_end, speA_tgt_str = gph.cutting_region_h[ ortho_id ][0]
    speB_gid, speB_tgt_sca, speB_reg_stt, speB_reg_end, speB_tgt_str = gph.cutting_region_h[ ortho_id ][1]
    p [ ortho_id, speA_gid, speA_tgt_sca, speA_reg_stt, speA_reg_end, speA_tgt_str, speB_gid, speB_tgt_sca, speB_reg_stt, speB_reg_end, speB_tgt_str ]
  end

end


