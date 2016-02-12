#!/usr/bin/ruby

# require 'minitest/unit'
require 'minitest/autorun'
require './GetStartPos.rb'

class TC_GetStartPos < MiniTest::Unit::TestCase

  def setup
    @speAgid        = "SPU_025302"
    @path2speApepfa = "spur_protein.fas"
    @path2speAscafa = "spur_scaffold.fas"
    @path2speBscafa = "hpul_scaffold.fas"
    @gph = GetPosHash.new( @speAgid, @path2speApepfa, 
                           @path2speAscafa, @path2speBscafa )
  end

  def test_sampletest
    assert_equal( "good result", @gph.sampletest)
  end

  def test_check_makeblastdb
    @gph.check_makeblastdb
    assert( FileTest.exist?( "#{@path2speAscafa}.nhr" ), message = "OK" )
  end

end

#  def create_protein_h( speA_protein_f )
#  def pos_for_target
#  def pos_for_query_self
#  def get_pos_by_blast( blastDB )


