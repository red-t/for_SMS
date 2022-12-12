import unittest
import os
import inspect
import sys

 # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)
     
     
from TESequenceBuilder import SequenceContainer


class Test_BasicStartSeqList(unittest.TestCase):
    def setUp(self):
        seqlist=["AAAAAAAAAA","TTT","CCC","GGG"]
        self.tt=SequenceContainer(seqlist)
        
    def test_getseq(self):
        tt=self.tt
        self.assertEqual(tt.getTESequence("$1").sequence,"AAAAAAAAAA")
        self.assertEqual(tt.getTESequence("$2").sequence,"TTT")
        self.assertEqual(tt.getTESequence("$3").sequence,"CCC")
        self.assertEqual(tt.getTESequence("$4").sequence,"GGG")
        
        
    def test_definition(self):
        tt=self.tt
        tt.addDefinition("d=$1")
        self.assertEqual(tt.getTESequence("d").sequence,"AAAAAAAAAA")
        tt.addDefinition("e=$1-")
        self.assertEqual(tt.getTESequence("e").sequence,"TTTTTTTTTT")
        tt.addDefinition("f=$1+")
        self.assertEqual(tt.getTESequence("f").sequence,"AAAAAAAAAA")
        

    def test_unusualdefinition(self):
        tt=self.tt
        tt.addDefinition("d =$1")
        self.assertEqual(tt.getTESequence("d").sequence,"AAAAAAAAAA")
        tt.addDefinition("e= $1-")
        self.assertEqual(tt.getTESequence("e").sequence,"TTTTTTTTTT")
        tt.addDefinition("f = $1+")
        self.assertEqual(tt.getTESequence("f").sequence,"AAAAAAAAAA")
        tt.addDefinition(" g = $1+ ")
        self.assertEqual(tt.getTESequence("g").sequence,"AAAAAAAAAA")
        
    def test_variousdefinitions(self):
        tt=self.tt
        tt.addDefinition("a=$1")
        self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAAAA")
        tt.addDefinition("b=a-")
        self.assertEqual(tt.getTESequence("b").sequence,"TTTTTTTTTT")
        tt.addDefinition("c=\"ATCG\"")
        self.assertEqual(tt.getTESequence("c").sequence,"ATCG")
        
    def test_tsdreading(self):
        tt=self.tt
        tt.addDefinition("a=$1+3bp")
        self.assertEqual(tt.getTESequence("a").tsd,3)
        tt.addDefinition("b=a-100bp")
        self.assertEqual(tt.getTESequence("b").tsd,100)
        tt.addDefinition("c=\"ATCG\"+12bp")
        self.assertEqual(tt.getTESequence("c").tsd,12)

class Test_Deletions(unittest.TestCase):
     def setUp(self):
          seqlist=["AAATTTAAA","CCCCAAAACCCC","AAACCCTTTGGGAAA","AAATTTCCC"]
          self.tt=SequenceContainer(seqlist)
          self.tt.addDefinition("a=$1")
          self.tt.addDefinition("b=$2")
          self.tt.addDefinition("c=$3")
          self.tt.addDefinition("x=$4")
        
     def test_singledeletion(self):
          tt=self.tt
          tt.addDefinition("d=a[4..6]")
          self.assertEqual(tt.getTESequence("a").sequence,"AAATTTAAA")
          self.assertEqual(tt.getTESequence("d").sequence,"AAAAAA")
          tt.addDefinition("e=b[5..8]")
          self.assertEqual(tt.getTESequence("b").sequence,"CCCCAAAACCCC")
          self.assertEqual(tt.getTESequence("e").sequence,"CCCCCCCC")
          
     def test_multideletion(self):
          tt=self.tt
          tt.addDefinition("d=c[4..6,10..12]")
          self.assertEqual(tt.getTESequence("c").sequence,"AAACCCTTTGGGAAA")
          self.assertEqual(tt.getTESequence("d").sequence,"AAATTTAAA")
     
     def test_overlappingdeletion(self):
          tt=self.tt
          tt.addDefinition("d=c[4..6,4..12]")
          self.assertEqual(tt.getTESequence("c").sequence,"AAACCCTTTGGGAAA")
          self.assertEqual(tt.getTESequence("d").sequence,"AAAAAA")
     
     def test_multideletion_strand1(self):
          tt=self.tt
          tt.addDefinition("d=c[4..6,10..12]+")
          tt.addDefinition("e=c[4..6,10..12]-")
          self.assertEqual(tt.getTESequence("c").sequence,"AAACCCTTTGGGAAA")
          self.assertEqual(tt.getTESequence("d").sequence,"AAATTTAAA")
          self.assertEqual(tt.getTESequence("e").sequence,"TTTAAATTT")
     
     def test_multideletion_strand2(self):
          tt=self.tt
          tt.addDefinition("y=x[4..6]+")
          tt.addDefinition("z=x[4..6]-")
          self.assertEqual(tt.getTESequence("x").sequence,"AAATTTCCC")
          self.assertEqual(tt.getTESequence("y").sequence,"AAACCC")
          self.assertEqual(tt.getTESequence("z").sequence,"GGGTTT")
          
     def test_carryoverofTSD(self):
          tt=self.tt
          tt.addDefinition("d=a+4bp")
          tt.addDefinition("e=d[4..6]")
          tt.addDefinition("f=d[4..6]+")
          tt.addDefinition("g=d[4..6]-")
          self.assertEqual(tt.getTESequence("d").tsd,4)
          self.assertEqual(tt.getTESequence("e").tsd,4)
          self.assertEqual(tt.getTESequence("f").tsd,4)
          self.assertEqual(tt.getTESequence("g").tsd,4)


     def test_overwriteDefaultTSD(self):
          tt=self.tt
          tt.addDefinition("d=a+4bp")
          tt.addDefinition("e=d[4..6]+2bp")
          tt.addDefinition("f=d[4..6]-6bp")
          self.assertEqual(tt.getTESequence("d").tsd,4)
          self.assertEqual(tt.getTESequence("e").tsd,2)
          self.assertEqual(tt.getTESequence("f").tsd,6)

          
     def test_weirdoposition(self):
          tt=self.tt
          tt.addDefinition("d=c[^..3,13..$]")
          self.assertEqual(tt.getTESequence("c").sequence,"AAACCCTTTGGGAAA")
          self.assertEqual(tt.getTESequence("d").sequence,"CCCTTTGGG")
          tt.addDefinition("e=b[|..$]")
          tt.addDefinition("f=b[^..|]")
          self.assertEqual(tt.getTESequence("b").sequence,"CCCCAAAACCCC")
          self.assertEqual(tt.getTESequence("e").sequence,"CCCCA")
          self.assertEqual(tt.getTESequence("f").sequence,"AACCCC")
          
          
     
     def test_alternativeDefinitions(self):
          tt=self.tt
          tt.addDefinition("d=$1[4..6]")
          self.assertEqual(tt.getTESequence("$1").sequence,"AAATTTAAA")
          self.assertEqual(tt.getTESequence("d").sequence,"AAAAAA")
          
          tt.addDefinition("e=\"CCCGGGAAA\"[4..6]")
          self.assertEqual(tt.getTESequence("e").sequence,"CCCAAA")


class Test_SimpleNested(unittest.TestCase):
     def setUp(self):
          seqlist=["AAAGGG","AAATTTGGG","CCC","TTT","GGG"]
          self.tt=SequenceContainer(seqlist)
          self.tt.addDefinition("a=$1")
          self.tt.addDefinition("b=$2")
          self.tt.addDefinition("c=$3")
          self.tt.addDefinition("t=$4")
          self.tt.addDefinition("g=$5")

     def test_singleinsertion(self):
          tt=self.tt
          tt.addDefinition("e=a+{3:c}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAGGG")
          self.assertEqual(tt.getTESequence("e").sequence,"AAACCCGGG")
          
          tt.addDefinition("f=a-{3:c}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAGGG")
          self.assertEqual(tt.getTESequence("f").sequence,"CCCCCCTTT")
          
     def test_singleinsertion_instanddefinition(self):
          tt=self.tt
          tt.addDefinition("e=a+{3:\"TAAT\"}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAGGG")
          self.assertEqual(tt.getTESequence("e").sequence,"AAATAATGGG")
          
          tt.addDefinition("f=a-{3:\"TAAT\"}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAGGG")
          self.assertEqual(tt.getTESequence("f").sequence,"CCCTAATTTT")
     
     def test_singleinsertion_tsd(self):
          tt=self.tt
          tt.addDefinition("e=\"AATCC\"")
          tt.addDefinition("f=\"GG\"+0bp")
          tt.addDefinition("f1=\"GG\"+1bp")
          tt.addDefinition("f2=\"GG\"+2bp")
          tt.addDefinition("f3=\"GG\"+3bp")
          tt.addDefinition("f4=\"GG\"+4bp")
          tt.addDefinition("f5=\"GG\"+5bp")
          tt.addDefinition("x=e+{3:f}")
          tt.addDefinition("x1=e+{3:f1}")
          tt.addDefinition("x2=e+{3:f2}")
          tt.addDefinition("x3=e+{3:f3}")

          
          self.assertEqual(tt.getTESequence("e").sequence,"AATCC")
          self.assertEqual(tt.getTESequence("x").sequence, "AATGGCC")
          self.assertEqual(tt.getTESequence("x1").sequence,"AATGGTCC")
          self.assertEqual(tt.getTESequence("x2").sequence,"AATGGATCC")
          self.assertEqual(tt.getTESequence("x3").sequence,"AATGGAATCC")
          
     def test_singleinsertion_instanddefinition_minusstrand(self):
          tt=self.tt
          tt.addDefinition("e=a+{3:\"TAAT\"-}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAGGG")
          self.assertEqual(tt.getTESequence("e").sequence,"AAAATTAGGG")
          
          tt.addDefinition("f=a-{3:\"TAAT\"-}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAGGG")
          self.assertEqual(tt.getTESequence("f").sequence,"CCCATTATTT")
          

class Test_MultiNested(unittest.TestCase):
     def setUp(self):
          seqlist=["AAAAGGGGAAAA"]
          self.tt=SequenceContainer(seqlist)
          self.tt.addDefinition("a=$1")


     def test_twoinsertions(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("res=a+{2:i1,6:i2}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAGGGGAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AATTTTAAGGCCCCGGAAAA")
          
     def test_twoinsertions_strand(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("res=a-{2:i1,6:i2}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAGGGGAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"TTTTTTTTCCCCCCCCTTTT")
          

     def test_twoinsertions_tsd(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+1bp")
          tt.addDefinition("i2=\"CCCC\"+2bp")
          tt.addDefinition("res=a+{2:i1,6:i2}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAGGGGAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AATTTTAAAGGCCCCGGGGAAAA")

class Test_DeepNested(unittest.TestCase):
     def setUp(self):
          seqlist=["AAAAAAAA"]
          self.tt=SequenceContainer(seqlist)
          self.tt.addDefinition("a=$1")


     def test_twoinsertions(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("res=a+{4:i1+{2:i2+}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AAAATTCCCCTTAAAA")

     def test_twoinsertions_strand(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("res=a+{4:i1+{2:i2-}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AAAATTGGGGTTAAAA")
     
     def test_threeinsertions_strand(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("i3=\"AAAA\"+0bp")
          tt.addDefinition("res=a+{4:i1+{2:i2-{2:i3+}}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AAAATTGGAAAAGGTTAAAA")

     def test_fourinsertions_strand(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("i3=\"AAAA\"+0bp")
          tt.addDefinition("res=a+{4:i1+{2:i2-{2:i3+{2:i2}}}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AAAATTGGAACCCCAAGGTTAAAA")
     
     def test_fourinsertions_multi(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"CCCC\"+0bp")
          tt.addDefinition("i3=\"AAAA\"+0bp")
          tt.addDefinition("res=a+{4:i1+{2:i2-{2:i3+{2:i2}}},6:i1}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AAAATTGGAACCCCAAGGTTAATTTTAA")


class Test_CombineNestedDeletions(unittest.TestCase):
     def setUp(self):
          seqlist=["AAAAAAAA"]
          self.tt=SequenceContainer(seqlist)
          self.tt.addDefinition("a=$1")

     def test_twoinsertions_onedeletion(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"ACGT\"+0bp")
          tt.addDefinition("res1=a+{4:i1+{2:i2+}}")
          tt.addDefinition("res2=a+{4:i1+{2:i2[1..2]+}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res1").sequence,"AAAATTACGTTTAAAA")
          self.assertEqual(tt.getTESequence("res2").sequence,"AAAATTGTTTAAAA")

     def test_twoinsertions_twodel(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"ACGT\"+0bp")
          tt.addDefinition("res1=a+{4:i1+{2:i2+}}")
          tt.addDefinition("res2=a+{4:i1+{2:i2[1..1,4..4]+}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res1").sequence,"AAAATTACGTTTAAAA")
          self.assertEqual(tt.getTESequence("res2").sequence,"AAAATTCGTTAAAA")


     def test_twoinsertions_fourdel(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"AAACCCGGGTTT\"+0bp")
          tt.addDefinition("res1=a+{4:i1+{2:i2+}}")
          tt.addDefinition("res2=a+{4:i1+{2:i2[1..1,4..4,7..7,9..9]+}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res1").sequence,"AAAATTAAACCCGGGTTTTTAAAA")
          self.assertEqual(tt.getTESequence("res2").sequence,"AAAATTAACCGTTTTTAAAA")
     
     def test_twoinsertions_fourdel_twoparellel(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"AAACCCGGGTTT\"+0bp")
          tt.addDefinition("res1=a+{4:i1+{2:i2+}}")
          tt.addDefinition("res2=a+{4:i1+{2:i2[1..1,4..4,7..7,9..9]+},6:i1,8:i1}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res1").sequence,"AAAATTAAACCCGGGTTTTTAAAA")
          self.assertEqual(tt.getTESequence("res2").sequence,"AAAATTAACCGTTTTTAATTTTAATTTT")

     def test_twoinsertions_fourdel_delmain(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"AAACCCGGGTTT\"+0bp")
          tt.addDefinition("res1=a[1..1,8..8]+{4:i1+{2:i2[1..1,4..4,7..7,9..9]+},6:i1}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res1").sequence,"AAAATTAACCGTTTTTAATTTT")
          

     def test_twoinsertions_fourdel_supernest1(self):
          tt=self.tt
          tt.addDefinition("i1=\"TTTT\"+0bp")
          tt.addDefinition("i2=\"AAACCCGGGTTT\"+0bp")
          tt.addDefinition("res1=a[1..1,8..8]+{4:i1+{2:i2[1..1,4..4,7..7,9..9]+},6:i1-{2:i1}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AAAAAAAA")
          self.assertEqual(tt.getTESequence("res1").sequence,"AAAATTAACCGTTTTTAAAATTTTAA")


class Test_CombineNestedTSD(unittest.TestCase):
     def setUp(self):
          seqlist=["AATTCCGGCCTTAA"]
          self.tt=SequenceContainer(seqlist)
          self.tt.addDefinition("a=$1")

     def test_twoinsertions(self):
          tt=self.tt
          tt.addDefinition("i1=\"ATTA\"+2bp")
          tt.addDefinition("i2=\"CCCC\"+2bp")
          tt.addDefinition("res=a+{7:i1+{2:i2+}}")
          self.assertEqual(tt.getTESequence("a").sequence,"AATTCCGGCCTTAA")
          self.assertEqual(tt.getTESequence("res").sequence,"AATTCCGATCCCCATTACGGCCTTAA")

                          
if __name__ == '__main__':
    unittest.main()
    
    
    
    